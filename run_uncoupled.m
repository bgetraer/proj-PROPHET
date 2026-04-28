steps=[1];

% location of proj-PROPHET
prophdir      = fileparts(mfilename('fullpath')); % the directory where this file lives
expdir        = fullfile(prophdir,'exp');
experimentdir = fullfile(prophdir,'experiments','uncoupled_experiments');

% directories where experiments live
modeldir     = fullfile(experimentdir,'Models');
outputdir    = fullfile(experimentdir,'output');
processeddir = fullfile(experimentdir,'processed');

experiment_options = { ...
	'control', ...           % 1        
	'collapse_400mpy', ...   % 2
	'budd_zero_melt', ...    % 3
	'budd_ABUM', ...         % 4
	'weertman_zero_melt', ...% 5
	'weertman_ABUM', ...     % 6
	'n4_zero_melt', ...      % 7
	'n4_ABUM', ...           % 8
	'};
experiment = experiment_options{7};
prefix=sprintf('PROPHET_%s_',experiment);

org=organizer('repository',modeldir,'prefix',prefix,'steps',steps);

if perform(org,'Param') % {{{
	% load md structure
	md=loadmodel(fullfile(prophdir,'experiments/ISSM_initialization/Models/PROPHET_issm_init_TransientPrep.mat')); % md structure

	% Basal forcings {{{
	disp('   -- Set basalforcings');
	switch experiment
		case {'control','budd_zero_melt','weertman_zero_melt','n4_zero_melt'};
			md.basalforcings = basalforcings();
			md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);
			md.basalforcings.geothermalflux           = zeros(md.mesh.numberofvertices,1);
			md.basalforcings.floatingice_melting_rate = zeros(md.mesh.numberofvertices,1); % m/year
		case {'collapse_400mpy','budd_ABUM','weertman_ABUM','n4_ABUM'}
			md.basalforcings.floatingice_melting_rate = 400*ones(md.mesh.numberofvertices,1); % m/year
		otherwise
			error('experiment undefined in Basal forcings');
	end % }}}
	% Rheology n {{{
	disp('   -- Set rheology_n')
	switch experiment
		case {'control','collapse_400mpy','budd_zero_melt','budd_ABUM','weertman_zero_melt','weertman_ABUM'}
			disp('  nothing to do');
		case {'n4_zero_melt','n4_ABUM'}
			rheology_n = 4;
			% extract B from the existing model (calculated for n=3)
			B3 = md.materials.rheology_B;
			% calculate effective strain rate from the MODELED velocities, do NOT average over the elements
			[strainrate] = strainrate_SSA(md,md.results.StressbalanceSolution.Vx, md.results.StressbalanceSolution.Vy,0);

			% calculate approximation for B, where mu(3) = mu(n), ie B3/eps_e^(2/3) = Bn/eps_e^((n-1)/n)
			disp(['Converting B for n=' num2str(rheology_n)]);
			Bn = B3.*(strainrate.eff/md.constants.yts).^((rheology_n-1)/rheology_n - 2.0/3.0);

			% update with new parameters for n for all ice vertices
			ind=strainrate.eff>0; % index the ice, not the ocean, ignore zero strainrate
			md.materials.rheology_B(ind)=Bn(ind); % update B
			md.materials.rheology_n = rheology_n * ones(size(md.materials.rheology_n)); % set n

			% test if the  model immediately converges to the same modeled velocity solution
			testconvergence = 1;
			if testconvergence
				fprintf(['\n' ...
					'*********************************************************\n' ...
					'TESTING CONVERGENCE: solution should converge in one step\n' ...
					'*********************************************************\n']);
				md.inversion.iscontrol=0;
				md.verbose.convergence = 1;
				mdtest = solve(md,'sb');
				clear mdtest;
			end
		otherwise
			error('experiment undefined in Rheology n');
	end % }}}
	savemodel(org,md);
end % }}}
if perform(org,'TransientRun') % {{{
	% load model
	md=loadmodel(org,'BasalForcings');

	% set options
	disp('setting transient options');
	executionpath = '/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/ISSM_control_run/run';
	md.cluster=generic('name',oshostname(),'np',75,'executionpath',executionpath,'interactive',0);
	md.verbose.solution=true;
	md.timestepping.start_time=2013;
	md.timestepping.final_time=2300;
	md.settings.output_frequency=20; % every 20 timesteps, about 1 per year



	% solve
	md.miscellaneous.name=sprintf('PROPH-%s-TransientRun',experiment);
	md=solve(md,'tr');

	% save
	savemodel(org,md);
end % }}}
if perform(org,'ProcessResults') % {{{
	disp('Loading results...');
	md = loadmodel(org,'BasalForcings');	
	md = loadresultsfromdisk(md,fullfile(outputdir,sprintf('PROPH-%s-TransientRun.outbin',experiment)));

	% constants
	rho_ice = 917; % kg/m^3
	kg2gt = 1e-12; % Gt/kg
	gt2mmsle = 1/361.8; % mmSLE/Gt


	% basin flags
	flagfile = fullfile(expdir,'basin_flags.mat');
	if exist(flagfile)
		disp('Loading sub-basin element flags');
		load(flagfile);
	else
		disp('Flagging elements in sub-basins');
		flag_basin21=FlagElements(md,fullfile(expdir,'reg21_thwaites.exp')); % flags for basin 21
		flag_basin22=FlagElements(md,fullfile(expdir,'reg22_pineisland.exp')); % flags for basin 22
		fprintf('Saving flags to %s\n',flagfile);
		save(flagfile,'flag_basin21','flag_basin22');
	end

	t = [md.results.TransientSolution.time]; % time (y)
	vaf = [md.results.TransientSolution.IceVolumeAboveFloatation]; % ice volume above flotation (m^3)

	ind = 1:numel(t);% which time-steps to take

	vaf_basin21 = zeros(size(vaf));
	vaf_basin22 = zeros(size(vaf));
	fprintf('Calculating VAF in sub-basins: %i timesteps\n',numel(ind));
	for j = 1:numel(ind)
		vaf_basin21(j) = VolumeAboveFloatation(md,j,flag_basin21); % basin volume above flotation at TransientSolution(j) (m^3)
		vaf_basin22(j) = VolumeAboveFloatation(md,j,flag_basin22); % basin volume above flotation at TransientSolution(j) (m^3)
	end

	% interpolating onto timesteps
	tq   = 2013:2300;
	vafq = interp1(t,vaf,tq);
	maf = vaf .* rho_ice .* kg2gt; % mass above flotation (Gt)
	sle = maf .* gt2mmsle; % sea level equivalence (mm)
	vaf_basin21q = interp1(t,vaf_basin21,tq);
	vaf_basin22q = interp1(t,vaf_basin22,tq);

	vafq(1)=vaf(1);
	vaf_basin21q(1)=vaf_basin21(1);
	vaf_basin22q(1)=vaf_basin22(1);


	vaf=vafq;
	vaf_basin21=vaf_basin21q;
	vaf_basin22=vaf_basin22q;

	filename = sprintf('issm_%s_interpolated_results',experiment);
	fprintf('   saving results to %s\n',fullfile(processeddir,filename));
	save(fullfile(processeddir,filename), 'tq','vaf','maf','sle','vaf_basin21','vaf_basin22');
	disp('   done');
end % }}}
