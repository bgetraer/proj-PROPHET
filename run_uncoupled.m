steps=[3];

% directories {{{
% location of proj-PROPHET
prophdir      = fileparts(mfilename('fullpath')); % the directory where this file lives
expdir        = fullfile(prophdir,'exp');
experimentdir = fullfile(prophdir,'experiments','uncoupled_experiments');

% directories where experiments live
modeldir     = fullfile(experimentdir,'Models');
rundir       = fullfile(experimentdir,'run');
outputdir    = fullfile(experimentdir,'output');
processeddir = fullfile(experimentdir,'processed');
inputdir     = fullfile(experimentdir,'inputdata');

% data directories
bedmachinepath='/totten_1/ModelData/Antarctica/BedMachine/BedMachineAntarctica-v4.0.nc'; % path to dataset
% }}}
% experiments {{{
experiment_options = { ...
	'control_zero_melt', ...   % 1        
	'control_ABUM', ...        % 2
	'budd_zero_melt', ...      % 3
	'budd_ABUM', ...           % 4
	'weertman_zero_melt', ...  % 5
	'weertman_ABUM', ...       % 6
	'n4_zero_melt', ...        % 7
	'n4_ABUM', ...             % 8
	'control_ISMIP_PW600', ... % 9 
	'control_ISMIP_PW800', ... % 10 
	'budd_ISMIP_PW600', ...    % 11 
	'budd_ISMIP_PW800', ...    % 12 
	'weertman_ISMIP_PW600', ...% 13
	'weertman_ISMIP_PW800', ...% 14
	'n4_ISMIP_PW600',...       % 15
	'n4_ISMIP_PW800'           % 16
};
experiment = experiment_options{10};
prefix=sprintf('PROPHET_%s_',experiment);
% }}}

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
		case {'control_zero_melt','budd_zero_melt','weertman_zero_melt','n4_zero_melt'};
			md.basalforcings = basalforcings();
			md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);
			md.basalforcings.geothermalflux           = zeros(md.mesh.numberofvertices,1);
			md.basalforcings.floatingice_melting_rate = zeros(md.mesh.numberofvertices,1); % m/year
		case {'control_ABUM','budd_ABUM','weertman_ABUM','n4_ABUM'}
			md.basalforcings.floatingice_melting_rate = 400*ones(md.mesh.numberofvertices,1); % m/year
		case {'control_ISMIP_PW600','n4_ISMIP_PW600','budd_ISMIP_PW600','weertman_ISMIP_PW600'}
			% calculate ISMIP6 style forcings from piecewise thermocline with set structure
			% see get_ismip_piecewise_basalforcing(), ISMIP_style_melt/runme.m, and run_sensitivity_experiments.m
			mit=loadmodel(fullfile(prophdir,'experiments/MITgcm_initialization/Models/PROPHET_mitgcm_init_CompileMITgcm.mat')); % mit structure
			depth_thermocline = 600; % bottom of thermocline (depth in m)
			md.basalforcings = get_ismip_piecewise_basalforcing(mit,md,bedmachinepath,depth_thermocline);

			% check visually that things make sense
			plotavgtf = 0;
			if plotavgtf
				figure(100); clf; hold on;
				tf = [md.basalforcings.tf{1,1,:}];
				avgtf = mean(tf(1:end-1,:)); % dont include date
				tf_depths = md.basalforcings.tf_depths;
				scatter(avgtf,tf_depths);
				xlabel('thermal forcing');
				ylabel('depth')
				yline(-depth_thermocline,'--k')
				yline(-depth_thermocline+400,'--k')
			end

		case {'control_ISMIP_PW800','n4_ISMIP_PW800','budd_ISMIP_PW800','weertman_ISMIP_PW800'}
			% calculate ISMIP6 style forcings from piecewise thermocline with set structure
			% see get_ismip_piecewise_basalforcing(), ISMIP_style_melt/runme.m, and run_sensitivity_experiments.m
			mit=loadmodel(fullfile(prophdir,'experiments/MITgcm_initialization/Models/PROPHET_mitgcm_init_CompileMITgcm.mat')); % mit structure
			depth_thermocline = 800; % bottom of thermocline (depth in m)
			md.basalforcings = get_ismip_piecewise_basalforcing(mit,md,bedmachinepath,depth_thermocline);

			% check visually that things make sense
			plotavgtf = 0;
			if plotavgtf
				figure(100); clf; hold on;
				avgtf = mean([md.basalforcings.tf{1,1,:}],1);
				tf_depths = md.basalforcings.tf_depths;
				scatter(avgtf,tf_depths);
				xlabel('thermal forcing');
				ylabel('depth')
				yline(-depth_thermocline,'--k')
				yline(-depth_thermocline+400,'--k')
			end

		otherwise
			error('experiment undefined in Basal forcings');
	end
	% }}}
	% Rheology n {{{
	disp('   -- Set rheology_n')
	switch experiment
		case {'control','collapse_400mpy','budd_zero_melt','budd_ABUM','weertman_zero_melt','weertman_ABUM', ...
				'control_ISMIP_PW600','control_ISMIP_PW800', 'budd_ISMIP_PW600','budd_ISMIP_PW800', ...
				'weertman_ISMIP_PW600','weertman_ISMIP_PW800'}
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
		case {'control','collapse_400mpy','n4_zero_melt','n4_ABUM',...
				'control_ISMIP_PW600','control_ISMIP_PW800','n4_ISMIP_PW600','n4_ISMIP_PW800'}
			disp('  nothing to do');
		case {'budd_zero_melt','budd_ABUM','budd_ISMIP_PW600','budd_ISMIP_PW800'}	
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

		case {'weertman_zero_melt','weertman_ABUM','weertman_ISMIP_PW600','weertman_ISMIP_PW800'}
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
	md.cluster=generic('name',oshostname(),'np',50,'executionpath',rundir,'interactive',0);
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
function [basalforcings] = get_ismip_piecewise_basalforcing(mit,md,bedmachinepath,depth_thermocline) % {{{
	% make ismip forcing 
	disp('Building ISMIP6-style basal forcings');

	disp('   -- Compute Boundary Conditions'); 
	[Tobw,Sobw,Tobs,Sobs] = get_mitgcm_boundary_conditions(mit,depth_thermocline);
	% T and S in linear index order 
	T = [Tobw(:);Tobs(:)]; % vector of T
	S = [Sobw(:);Sobs(:)]; % vector of S

	disp('   -- Compute Boundary Grid'); 
	% domain geometry in linear index order (matches MITgcm because we use 3D meshgrid)
	Xobw = squeeze(mit.mesh.XC(:,1,:)); % X OBW points
	Yobw = squeeze(mit.mesh.YC(:,1,:)); % Y OBW points
	Zobw = squeeze(mit.mesh.ZC(:,1,:)); % Z OBW points
	Xobs = squeeze(mit.mesh.XC(1,:,:)); % X OBS points
	Yobs = squeeze(mit.mesh.YC(1,:,:)); % Y OBS points
	Zobs = squeeze(mit.mesh.ZC(1,:,:)); % Z OBS points
	X=[Xobw(:);Xobs(:)]; % vector of X points
	Y=[Yobw(:);Yobs(:)]; % vector of Y points
	Z=[Zobw(:);Zobs(:)]; % vector of Z points

	% only look at cells with data
	ind = find((S~=0));
	T=T(ind);
	S=S(ind);
	X=X(ind);
	Y=Y(ind);
	Z=Z(ind);

	% define the real levels of the data
	levels=flip(unique(Z));

	disp('   -- Making Large Interpolation Grid');
	xq=min(mit.mesh.xp):1e3:max(mit.mesh.xp);
	yq=min(mit.mesh.yp):1e3:max(mit.mesh.yp);
	[Xq Yq]=meshgrid(xq,yq);

	% interpolate BedMachine bed and mask
	BED = interpBedmachineAntarctica(Xq,Yq,'bed','linear',bedmachinepath);
	M   = interpBedmachineAntarctica(Xq,Yq,'mask','nearest',bedmachinepath);

	disp('   -- Define open ocean pixels');
	[~,indx0] =min(abs(xq-mit.mesh.xc(1))); % find the nearest x in the ISMIP6 Grid to the OBW
	[~,indy0] =min(abs(yq-mit.mesh.yc(1))); % find the nearest y in the ISMIP6 Grid to the OBS
	[~,indxend] = min(abs(xq-max(X))); % find the nearest x in the ISMIP6 Grid to boundary end
	[~,indyend] = min(abs(yq-max(Y))); % find the nearest y in the ISMIP6 Grid to boundary end
	FAR = (Xq>=xq(indx0) & Xq<=xq(indxend) & Yq==yq(indy0)) | (Yq>=yq(indy0) & Yq<=yq(indyend) & Xq==xq(indx0));

	% optional plotting
	plot_domain=0;
	if plot_domain
		figure(100);
		subplot(1,3,1);
		imagesc(BED); title('BED'); set(gca,'ydir','normal');axis equal tight off;
		subplot(1,3,2);
		imagesc(M); title('MASK'); set(gca,'ydir','normal');axis equal tight off;
		subplot(1,3,3);
		imagesc(FAR); title('FAR'); set(gca,'ydir','normal');axis equal tight off;
		colormap(flip(gray));
	end

	% define tf_depths, the levels that we will interpolate onto
	tf_depths = -(25:50:950)';

	%Initialize output
	disp('   -- Build connectivity and ID matrices');
	connectivity = false(size(Xq,1),size(Xq,2),numel(tf_depths));
	ID=zeros(size(Xq)); % deepest index
	T_level=zeros(numel(levels),1); % level averaged temp, before interpolation
	S_level=zeros(numel(levels),1); % level averaged salinity, before interpolation

	% label depths
	for i=1:numel(tf_depths),
		%disp(['   -- Depth = ' num2str(tf_depths(i)) ' ' num2str(i) '/' num2str(numel(tf_depths))]);
		%disp('    ... Labelling');
		CC=bwlabel(BED<=tf_depths(i));
		%disp('    ... Finding unique labels');
		pos=find(FAR & BED<tf_depths(i));
		list = unique(CC(pos));
		connectivity(:,:,i)=ismember(CC,list);
		ID(connectivity(:,:,i))=i;
	end
	% average over levels
	for i=1:numel(levels),
		%disp(['   -- Depth = ' num2str(levels(i)) ' ' num2str(i) '/' num2str(numel(levels))]);
		%disp('    ... Level averaged temp and salt');
		ind = find(Z==levels(i));
		T_level(i,:)=mean(T(ind),1);
		S_level(i,:)=mean(S(ind),1);
	end

	% interpolate onto tf_depths
	Tq = interp1(levels,T_level,tf_depths);
	Sq = interp1(levels,S_level,tf_depths);

	% constants for thermal forcing
	rho_w=1029.35; % approximate value for estimating pressure (kg/m^3)
	% Calculate in situ freezing point (from Holland, Jenkins, and Holland 2008)
	g = 9.81;    % m/s^2
	a = -0.0573; % deg C
	b = 0.0832;  % deg C
	c = 7.53E-3*1E-5*rho_w.*g; % deg C/Pa

	% build thermal forcings structure
	disp('   -- Build thermal forcings structure');
	tf = cell(1,1,numel(tf_depths));

	for i=1:numel(tf_depths)
		%disp(['   -- Depth = ' num2str(tf_depths(i)) ' ' num2str(i) '/' num2str(numel(tf_depths))]);

		% min depth indexing
		levelID=min(ID,i);
		% initialize the state fields for this level
		temperature  =NaN([numel(yq),numel(xq)]);
		salinity     =NaN([numel(yq),numel(xq)]);
		% initialize the md.basalforcings.tf field time-series for this level
		tf{1,1,i}=single(zeros(md.mesh.numberofvertices+1,1));

		[posi posj]=find(levelID~=0);
		ind = sub2ind(size(levelID),posi,posj);
		temperature(ind) = Tq(levelID(ind)); % in-situ temp at min depth (deg C)
		salinity(ind)    = Sq(levelID(ind)); % in-situ salt at min depth
		freezingpoint = a.*salinity + b + c.*tf_depths(i); % deg C
		theta = temperature - freezingpoint; % thermal forcing (deg C)
		theta(isnan(theta)) = 0; % replace NaN values for min level zero
		tf{1,1,i}(1:end-1) = single(interp2(xq,yq,theta,md.mesh.x,md.mesh.y,'linear',0)); % tf at md vert. for this time (deg C)
		tf{1,1,i}(end) = 2013;
	end

	disp('   -- Define basalforcings');
	% build basalforcings field
	basalforcings = basalforcingsismip6();
	basalforcings.num_basins = 1;
	basalforcings.basin_id   = 1*ones(md.mesh.numberofelements,1);
	basalforcings.gamma_0    = 14477.3368;
	basalforcings.tf_depths  = tf_depths';
	basalforcings.tf         = tf;
	basalforcings.delta_t    = 1.066526770591736; % taken from basin 10, ISMIP6 Antarctica
	basalforcings.islocal    = 0;
	basalforcings.geothermalflux           = zeros(md.mesh.numberofvertices,1);
	basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);
end % }}}
function [OBWtheta,OBWsalt,OBStheta,OBSsalt] = get_mitgcm_boundary_conditions(mit,depth_thermocline) % {{{
	% SEE run_sensitivity_experiments.m
	disp('   -- Calculate piecewise boundary conditions for the ocean domain')

	% Define baseline ocean conditions for synthetic boundary conditions (after De Rydt et al., 2014 and Bett et al., 2024)
	z_top       = -300; % baseline z-coordinate of the top of the thermocline (winter water layer ends) (m)
	z_bot       = -700; % baseline z-coordinate of the bottom of the thermocline (circumpolar deep water layer begins) (m)
	theta_ww    = -1.0; % potential temperature of winter water (deg C)
	theta_cdw   = +1.2; % potential temperature of circumpolar deep water (deg C)
	salt_ww     = 34;   % salinity of winter water (g/kg)
	salt_cdw    = 34.7; % salinity of circumpolar deep water (g/kg)
	dtheta_cdw_OBS = -0.6; % potential temperature anomaly at "southern" open boundary (deg C)
	dsalt_cdw_OBS  = -0.1; % salinity anomaly at "southern" open boundary (g/kg)

	% define dummy z-coordinate to reflect thermocline depth anomaly
	dz_thermocline = (-depth_thermocline) - z_bot; % difference between experiment thermocline bottom depth and baseline depth (m)
	zc_prime = mit.mesh.zc - dz_thermocline; % dummy z-coordinates for calculating isocline (m)

	% reference profiles of potential temperature and salinity
	tRef = get_isocline(zc_prime,z_top,z_bot,theta_ww,theta_cdw);
	sRef = get_isocline(zc_prime,z_top,z_bot,salt_ww,salt_cdw);

	% boundary profiles of potential temperature and salinity (Nz x Nrecords)
	OBWtheta_prof = tRef;
	OBWsalt_prof  = sRef;

	OBStheta_prof = get_isocline(zc_prime,z_top,z_bot,theta_ww,theta_cdw + dtheta_cdw_OBS);
	OBSsalt_prof  = get_isocline(zc_prime,z_top,z_bot,salt_ww, salt_cdw  + dsalt_cdw_OBS);

	% Western Open Boundary
	OBWtheta = repmat(OBWtheta_prof,mit.mesh.Ny,1);
	OBWsalt  = repmat(OBWsalt_prof, mit.mesh.Ny,1);
	% apply mask
	OBWtheta = OBWtheta.*squeeze(mit.geometry.open_mask(:,1,:));
	OBWsalt  = OBWsalt .*squeeze(mit.geometry.open_mask(:,1,:));

	% Southern Open Boundary
	OBStheta = repmat(OBStheta_prof,mit.mesh.Nx,1);
	OBSsalt  = repmat(OBSsalt_prof, mit.mesh.Nx,1);
	% apply mask
	OBStheta = OBStheta.*squeeze(mit.geometry.open_mask(1,:,:));
	OBSsalt  = OBSsalt .*squeeze(mit.geometry.open_mask(1,:,:));

	% check sizes to confirm
	assert(all(size(OBWtheta)==[mit.mesh.Ny,mit.mesh.Nz]));
	assert(all(size(OBWsalt) ==[mit.mesh.Ny,mit.mesh.Nz]));
	assert(all(size(OBStheta)==[mit.mesh.Nx,mit.mesh.Nz]));
	assert(all(size(OBSsalt) ==[mit.mesh.Nx,mit.mesh.Nz]));
end	% }}}
function A = get_isocline(z,z0,z1,A0,A1) % {{{
	m = [z0,1;z1,1]\[A0;A1];
	A = m(1).*z + m(2);
	A(z>z0)=A0;
	A(z<z1)=A1;
end % }}}
