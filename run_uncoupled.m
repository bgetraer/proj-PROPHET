steps=[4];

% location of proj-PROPHET
prophdir      = fileparts(mfilename('fullpath')); % the directory where this file lives
expdir        = fullfile(prophdir,'exp');
experimentdir = fullfile(prophdir,'experiments','uncoupled_experiments');

% directories where experiments live
modeldir     = fullfile(experimentdir,'Models');
rundir       = fullfile(experimentdir,'run');
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
	};
experiment = experiment_options{6};
prefix=sprintf('PROPHET_%s_',experiment);

org=organizer('repository',modeldir,'prefix',prefix,'steps',steps);

if perform(org,'ConvergeSB') % {{{
	% load md structure
	md=loadmodel(fullfile(prophdir,'experiments/ISSM_initialization/Models/PROPHET_issm_init_TransientPrep.mat')); % md structure
	md.cluster = generic('np',80);
	disp('  resolving stressbalance');
	md.inversion.iscontrol=0;
	md.verbose.convergence = 1;
	md = solve(md,'sb');
	% Initialize modeled velocities from stress balance
	md.initialization.vx = md.results.StressbalanceSolution.Vx;
	md.initialization.vy = md.results.StressbalanceSolution.Vy;
	md.initialization.vel = md.results.StressbalanceSolution.Vel;

	% save file
	filename = 'PROPHET_ConvergeSB';
	save(fullfile(modeldir,filename),'md');
end % }}}
if perform(org,'Param') % {{{
	% load md structure
	filename = 'PROPHET_ConvergeSB';
	md=loadmodel(fullfile(modeldir,filename));
	md.cluster = generic('np',25);

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

			% NOTE: the initial model had md.materials.rheology_B defined on the vertices, and so 
			% rescaling for n=4 will not give an exact match. 
			% I deel with this by first recalculating md.materials.rheology_B on the elements, then
			% rescaling. While this approach is not a perfect initial condition match, it minimizes 
			% the distortion compared to keeping is on the vertices or re-averaging back.
			% Here is code for exploring the mismatch {{{
			%% test if the  model immediately converges to the same modeled velocity solution
			% md=loadmodel(fullfile(prophdir,'experiments/ISSM_initialization/Models/PROPHET_issm_init_InversionC.mat')); % md structure
			%testconvergence(md);

			%md1 = md;

			%% test with element-wise defintion
			%md.materials.rheology_B = mean(md.materials.rheology_B(md.mesh.elements),2);
			%md.inversion.iscontrol=0;
			%md.verbose.convergence = 1;
			%md = solve(md,'sb');
			%% Initialize modeled velocities from stress balance
			%md.initialization.vx = md.results.StressbalanceSolution.Vx;
			%md.initialization.vy = md.results.StressbalanceSolution.Vy;
			%md.initialization.vel = md.results.StressbalanceSolution.Vel;
			%md.initialization.vx = md.results.StressbalanceSolution.Vx;


			%rheology_n = 4;
			%% extract B from the existing model (calculated for n=3)
			%B3 = md.materials.rheology_B;
			%% calculate effective strain rate from the MODELED velocities, do NOT average over the elements
			%[strainrate] = strainrate_SSA(md,md.initialization.vx, md.initialization.vy,0);

			%% calculate approximation for B, where mu(3) = mu(n), ie B3/eps_e^(2/3) = Bn/eps_e^((n-1)/n)
			%disp(['Converting B for n=' num2str(rheology_n)]);
			%Bn = B3.*(strainrate.eff/md.constants.yts).^((rheology_n-1)/rheology_n - 2.0/3.0);

			%% update with new parameters for n for all ice vertices
			%ind=strainrate.eff>0; % index the ice, not the ocean, ignore zero strainrate
			%md.materials.rheology_B(ind)=Bn(ind); % update B
			%md.materials.rheology_n = rheology_n * ones(size(md.materials.rheology_n)); % set n

			%% back to vertices 
			%md.materials.rheology_B = averaging(md,md.materials.rheology_B,0);
			%% test if the  model immediately converges to the same modeled velocity solution
			%md = testconvergence(md);
			% }}}

			% set rheology_B per element
			disp('  setting rheology_B per element')
			md.materials.rheology_B = mean(md.materials.rheology_B(md.mesh.elements),2);
			md.inversion.iscontrol=0;
			md.verbose.convergence = 1;
			disp('  resolving stressbalance');
			md = solve(md,'sb');
			% Initialize modeled velocities from stress balance
			md.initialization.vx = md.results.StressbalanceSolution.Vx;
			md.initialization.vy = md.results.StressbalanceSolution.Vy;
			md.initialization.vel = md.results.StressbalanceSolution.Vel;
			md.initialization.vx = md.results.StressbalanceSolution.Vx;

			rheology_n = 4;
			% extract B from the existing model (calculated for n=3)
			B3 = md.materials.rheology_B;
			% calculate effective strain rate from the MODELED velocities, on the elements!
			[strainrate] = strainrate_SSA(md,md.initialization.vx, md.initialization.vy,0);

			% calculate approximation for B, where mu(3) = mu(n), ie B3/eps_e^(2/3) = Bn/eps_e^((n-1)/n)
			disp(['Converting B for n=' num2str(rheology_n)]);
			Bn = B3.*(strainrate.eff/md.constants.yts).^((rheology_n-1)/rheology_n - 2.0/3.0);

			% update with new parameters for n for all ice vertices
			ind=strainrate.eff>0; % index the ice, not the ocean, ignore zero strainrate
			md.materials.rheology_B(ind)=Bn(ind); % update B
			md.materials.rheology_n = rheology_n * ones(size(md.materials.rheology_n)); % set n

			% check convergence
			md = testconvergence(md);
		otherwise
			error('experiment undefined in Rheology n');
	end % }}}
	% Friction C {{{
	% NOTE: Trying to rescale frictionschoof() to other laws has been hard to get perfect convergence.
	% Despite the calculated basal stress agreeing within ~1E-9, the resulting stress balance is not
	% converged, and differences in the modeled velocities of up to 400 m/year appear.
	disp('   -- Set friction C')
	switch experiment
		case {'control','collapse_400mpy','n4_zero_melt','n4_ABUM'}
			disp('  nothing to do');
		case {'budd_zero_melt','budd_ABUM'}	
			%Change sliding law from Schoof to Budd, rescaling changing the friction coefficient	
			tau_b = basalstress(md); % basal stress for schoof (Pa)

			% Budd law: tau_b = C_budd.^2 .* N .* ub
			N = max(0.1, effectivepressure(md)); % (Pa)
			ub = sqrt(md.initialization.vx.^2 + md.initialization.vy.^2) /md.constants.yts; % m/s
			C_budd = sqrt(tau_b ./ (N .* ub)); 

			% set new friction
			md2 = md;
			md2.friction = friction();
			md2.friction.coefficient = C_budd;
			md2.friction.p = ones(md.mesh.numberofelements,1); % s = 1/p
			md2.friction.q = ones(md.mesh.numberofelements,1); % r = q/p
			md2.friction.coupling = 2; % keep the pressure calculation the same
			md2.friction.linearize = 1; % make linear over elements

			% compare the new basal stress to the old one is the same
			tau_b2 = basalstress(md2);
         check_nan = isnan(tau_b2) == isnan(tau_b);
         assert(all(check_nan));
			tau_b_diff = max(abs(tau_b2 - tau_b));
         fprintf('RESCALING: max tau_b difference of: %0.3e\n',tau_b_diff);

			% mask singular values
			ind = isnan(md2.friction.coefficient);
			md2.friction.coefficient(ind) = md.friction.C(ind);

			% check convergence
			md2.verbose.convergence = 1;
			md2 = solve(md2,'sb');

			% how off are the velocities
			vel_diff = max(abs(md2.results.StressbalanceSolution.Vel - md2.initialization.vel));
			fprintf('RECONVERGED: max vel difference of %0.3f\n',vel_diff);

			% Initialize modeled velocities from stress balance
         md2.initialization.vx  = md2.results.StressbalanceSolution.Vx;
         md2.initialization.vy  = md2.results.StressbalanceSolution.Vy;
         md2.initialization.vel = md2.results.StressbalanceSolution.Vel;

			% reset md
			md = md2;

		case {'weertman_zero_melt','weertman_ABUM'}
			%Change sliding law from Schoof to Weertman, rescaling changing the friction coefficient
			tau_b = basalstress(md); % basal stress for schoof (Pa)

			% Weertman law: tau_b = C_weertman.^2 .* ub.^(1./m-1) .* ub
			m = 3; % weertman exponent
			ub = sqrt(md.initialization.vx.^2 + md.initialization.vy.^2) /md.constants.yts; % m/s
			C_weertman = sqrt(tau_b ./ ub.^(1./m));

			% set new friction
         md2 = md;
			md2.friction = frictionweertman();
			md2.friction.C = C_weertman;
			md2.friction.m = m .* ones(md.mesh.numberofelements,1);
			md2.friction.linearize = 0; % CHANGED to 0 because 1 would not converge...

			% compare the new basal stress to the old one is the same
         tau_b2 = basalstress(md2);
         check_nan = isnan(tau_b2) == isnan(tau_b);
         assert(all(check_nan));
         tau_b_diff = max(abs(tau_b2 - tau_b));
         fprintf('RESCALING: max tau_b difference of: %0.3e\n',tau_b_diff);

			% mask singular values
			ind = isnan(md2.friction.C);
         md2.friction.C(ind) = md.friction.C(ind);

			% check convergence
         md2.verbose.convergence = 1;
         md2 = solve(md2,'sb');

			% how off are the vel
			vel_diff = max(abs(md2.results.StressbalanceSolution.Vel - md2.initialization.vel));
         fprintf('RECONVERGED: max vel difference of %0.3f\n',vel_diff);

         % Initialize modeled velocities from stress balance
         md2.initialization.vx  = md2.results.StressbalanceSolution.Vx;
         md2.initialization.vy  = md2.results.StressbalanceSolution.Vy;
         md2.initialization.vel = md2.results.StressbalanceSolution.Vel;

         % reset md
         md = md2;

		otherwise 
			error('experiment undefined in Friction C');
	end % }}}
	savemodel(org,md);
end % }}}
if perform(org,'TransientRun') % {{{
	% load model
	md=loadmodel(org,'Param');

	% set options
	disp('setting transient options');
	md.cluster=generic('name',oshostname(),'np',25,'executionpath',rundir,'interactive',0);
	md.verbose.convergence = false;
	md.verbose.solution=true;
	md.timestepping.start_time=2013;
	md.timestepping.final_time=2100;
	md.settings.output_frequency=20; % every 20 timesteps, about 1 per year

	% save
	savemodel(org,md);

	% solve
	md.miscellaneous.name=sprintf('PROPH-%s-TransientRun',experiment);
	md=solve(md,'tr');
end % }}}
if perform(org,'ProcessResults') % {{{
	disp('Loading results...');
	md = loadmodel(org,'Param');	
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

	plot(t,sle-sle(1),'-x')
end % }}}

% functions
function md = testconvergence(md) %{{{
	% testconvergence - test if the  model immediately converges to the same modeled velocity solution
	% test_converged(md)

	fprintf(['\n' ...
		'*********************************************************\n' ...
		'TESTING CONVERGENCE: solution should converge in one step\n' ...
		'*********************************************************\n']);
	md.inversion.iscontrol=0;
	md.verbose.convergence = 1;
	md = solve(md,'sb');
end % }}}
function md = get_frictionbudd(md) % {{{
	% initial guess from driving stress
	[sx,sy,s]=slope(md,md.geometry.surface); % slope 's' comes on elements
	sslope=averaging(md,s,1); % average the slope once on the vertices, because 's' comes on elements, we need this data on vertices

	vel=md.initialization.vel; % m/yr
	vel=vel/md.constants.yts;  % m/s

	N = (md.materials.rho_ice*md.geometry.thickness+md.materials.rho_water*md.geometry.base)*md.constants.g;
	N = max(0,N); % setting minimum positive pressure

	driving_stress=md.materials.rho_ice*md.constants.g*md.geometry.thickness.*(sslope);
	C=sqrt(driving_stress./(N.*vel));
end
% }}}
