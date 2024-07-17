% RUN_ISSM_INIT is a script to initialize the ISSM ice sheet model for the 
% PROPHET Amundsen Sea Coupling experiment.
% 
% Overview of experiment scripts:
%  RUN_ISSM_INIT.m      Initialize ice sheet model
%  RUN_MITGCM_INIT.m    Initialize ocean model
%  RUN_ISSMxMITGCM.m    Run coupled model
%
% Overview of this script:
%  1) initialize the domain mesh and parameters of the ice sheet model
%  2) invert for the flow parameter B and the friction parameter C
% 
% Output:
%  md is a structure of the ISSM model class that contains the initial state
%      of the model.
%
% https://github.com/bgetraer/proj-PROPHET.git

steps=[3];

experiment.name='ISSM_initialization';
% directory structure {{{
proph_dir =pwd; % base directory for this project
% make Exp directory if needed
% this will hold the .exp files for ISSM domain
if ~exist(fullfile(proph_dir,'Exp'))
	mkdir(fullfile(proph_dir,'Exp'));
end
% }}}
% define experiment directories {{{
expdir=fullfile(proph_dir,'experiments',experiment.name);
% make this experiment directory if needed
if ~exist(expdir)
   mkdir(expdir);
end
% make model directory if needed
% this will hold md structures from ISSM
modeldir=fullfile(expdir,'Models');
if ~exist(modeldir)
   mkdir(modeldir);
end
% }}}

% Ice model files
expfile='./Exp/domain.exp';

org=organizer('repository',modeldir,'prefix','PROPHET_test','steps',steps);
if perform(org,'MeshParam'),   % Build ISSM mesh and parameterize{{{ 
	md=model(); % initialize ISSM model structure
	% Exp/amundsenicedomain.exp {{{
	EXP=struct; % initialize domain exp structure
	% outline of ice domain in x (m)
	EXP.x=[-1669315.4719493026,-1669315.4719493026,-1193987.0047960179,...
		-1026979.7055259449,-1026979.7055259449,-1556906.7128252152,...
		-1772089.1945770399,-1772089.1945770399,-1669315.4719493026]; 
	% outline of ice domain in y (m)
	EXP.y=[-420940.0927398402,-829553.2314259715,-829553.2314259715,...
		-530867.1000391102,-58750.3117179424,170008.8123696489,...
		70446.7685740285,-420940.0927398402,-420940.0927398402]; 
	% write to file
	write_expfile(expfile,EXP);
	% }}}
	% Mesh {{{
	hmin=1000; % min edge length (m)
	hmax=15e3; % max edge length (m)
	md=triangle(md,expfile,hmin); % initialize unstructured triangular mesh

	nsteps = 2; % number of mesh adaptation steps
	for i=1:nsteps,
		disp(['--- Performing static mesh adaptation. Step ' num2str(i) '/' num2str(nsteps)]);
		% using a priori analysis (observed velocity)
		disp('   -- Interpolating some data');
		[velx vely] = interpMouginotAnt2017(md.mesh.x,md.mesh.y); % interpolate observed velocities (nominal year 2013) (m/yr)
		ocean_levelset=-ones(size(md.mesh.x));% all floating
		ocean_levelset(find(interpBedmachineAntarctica(md.mesh.x,md.mesh.y,'mask')==2))=1; % grounded from BedMachine
		ice_levelset=-ones(size(md.mesh.x));
		ice_levelset(find(interpBedmachineAntarctica(md.mesh.x,md.mesh.y,'mask')==0))=1; % grounded from BedMachine

		pos=find(isnan(velx) | isnan(vely) | ice_levelset>0);% | ocean_levelset<0);
		velx(pos)=0; vely(pos)=0; vel=sqrt(velx.^2+vely.^2);

		hVertices = NaN(md.mesh.numberofvertices,1);
		hVertices(find(vel>200)) = hmin;
		md=bamg(md,'gradation',1.6180,'anisomax',1.e6,'KeepVertices',0,'Hessiantype',0,'Metrictype',0,...
			'hmax',hmax,'hmin',hmin,'hVertices',hVertices,'field',vel,'err',3);

		md.private.bamg=struct();
	end
	% }}}
	% Param {{{
	%  Set parameters and options for the ISSM model of the initial ice state
   
	% Model name
   md.miscellaneous.name='PROPHET_ISSM';

	% Material propoerties
	% Englacial temperatures from PIG: Truffer & Stanton (2015); Mulvaney (2017). Median temp ~ -22 deg C
	temp0 = 273.15-22; % approximate englacial temperature of the ice (for A and initialization) (deg K) 
   md.materials.rheology_B=cuffey(temp0)*ones(md.mesh.numberofvertices,1); % Cuffey & Paterson, for T=-10degC (Pa s^1/n)
   md.materials.rho_water=1028; % ocean water density (kg/m^3)

	% Geometry from BedMachine Antarctica
	md.geometry.surface = interpBedmachineAntarctica(md.mesh.x,md.mesh.y,'surface','linear'); % surface elevation data (m)
	md.geometry.bed     = interpBedmachineAntarctica(md.mesh.x,md.mesh.y,'bed','linear'); % bed elevation data (m)
	% Grounded ice base
	md.geometry.base    = md.geometry.bed; % assume initially all grounded (m)
	% Open ocean minimum ice surface
	min_surface = md.masstransport.min_thickness*(1-md.materials.rho_ice/md.materials.rho_water); % minimum surface allowed on open ocean (m)
	md.geometry.surface(md.geometry.surface<min_surface) = min_surface; % set minimum ice surface on the ocean part (m)
	% Floating ice base
	flotation_base = md.geometry.surface*md.materials.rho_ice/(md.materials.rho_ice-md.materials.rho_water); % base if all ice was floating (m)
	md.geometry.base(flotation_base>md.geometry.bed) = flotation_base(flotation_base>md.geometry.bed); % floating ice base (m)
	% Ice thickness
	md.geometry.thickness            = md.geometry.surface-md.geometry.base; % ice thickness (m)
	% Grounded ice minimum thickness
	ind = md.geometry.thickness<md.masstransport.min_thickness;
	md.geometry.thickness(ind) = md.masstransport.min_thickness; % assert minimum thickness (m)
	md.geometry.surface(ind) = md.geometry.thickness(ind) + md.geometry.base(ind); % recalculate surfaces

	% Ocean mask
	md.mask.ocean_levelset = ones(md.mesh.numberofvertices,1);  % set 'grounded ice' everywhere
	ind_floating = md.geometry.surface*md.materials.rho_ice/(md.materials.rho_ice-md.materials.rho_water) > md.geometry.bed; % updated flotation criterion
	md.mask.ocean_levelset(ind_floating) = -1; % set 'floating' ice mask based on interpolated data
	md.mask.ocean_levelset = reinitializelevelset(md, md.mask.ocean_levelset);
	% Ice mask
	md.mask.ice_levelset = -1*ones(md.mesh.numberofvertices,1); % set 'presence of ice' everywhere
	ind_no_ice = (md.mask.ocean_levelset<0 & abs(md.geometry.thickness-md.masstransport.min_thickness)<1E-10); % floating and minimum thickness
	md.mask.ice_levelset(ind_no_ice) = +1; % set 'no ice' on the only ocean part 
   % Some checks on the geometry setup {{{
   if any(isnan(md.geometry.bed)) | any(isnan(md.geometry.surface))
      error('NaN was found in the data set!')
   end
   if any(md.geometry.surface<0)
      error('surface < 0')
   end
   if any(md.geometry.thickness<md.masstransport.min_thickness)
      error('thickness < min_thickness)!')
   end
   if any(md.geometry.base~=md.geometry.bed & md.mask.ocean_levelset>0)
      error('base is not equal to bed on grounded ice!')
   end
   if any(abs(md.geometry.thickness-(md.geometry.surface-md.geometry.base))>1E-10)
      error('thickness is not equal to surface-base!');
   end
   if any(md.geometry.base<md.geometry.bed & md.mask.ocean_levelset<0)
      error('base < bed on floating ice')
   end
	if any(md.mask.ocean_levelset<0 & abs(md.geometry.surface - md.geometry.thickness*(1-md.materials.rho_ice/md.materials.rho_water)) > 1E-10)
		error('floating ice is not floating at hydrostatic equilibrium')
	end
	if any(md.mask.ice_levelset>0 & abs(md.geometry.thickness-md.masstransport.min_thickness)>1E-10)
		error('ice is too thick over the ocean only domain')
	end
   % }}}
	% Adjust ice mask
	pos = max(md.mask.ice_levelset(md.mesh.elements),[],2)>0; % elements with at least 1 vertex that has no ice
	md.mask.ice_levelset(md.mesh.elements(pos,:))= +1; % set 'no ice' on the partial ice elements to force ice front boundary upstream of the cliff
	md.mask.ice_levelset   = reinitializelevelset(md, md.mask.ice_levelset);

	% Velocities from MEaSUREs
	[md.inversion.vx_obs, md.inversion.vy_obs] = interpMouginotAnt2017(md.mesh.x,md.mesh.y); % interpolated surface velocity component data (m yr^-1)
	pos_vel_nan = find(isnan(md.inversion.vx_obs) | isnan(md.inversion.vy_obs)); % find location of NaN values
	save(fullfile(modeldir, 'pos_vel_nan.mat'),'pos_vel_nan');
	md.inversion.vx_obs(pos_vel_nan) = 0; % fill NaN values with zero (m yr^-1)
	md.inversion.vy_obs(pos_vel_nan) = 0; % fill NaN values with zero (m yr^-1)
   md.inversion.vel_obs  = sqrt(md.inversion.vx_obs.^2+md.inversion.vy_obs.^2); % velocity magnitude (m yr^-1)
	% initialization velocities
	md.initialization.vx = md.inversion.vx_obs; % x component of initialization velocity (m yr^-1)
	md.initialization.vy = md.inversion.vy_obs; % y component of initialization velocity (m yr^-1)
   md.initialization.vz  = zeros(md.mesh.numberofvertices,1); % z component of initialization velocity (m yr^-1)
   md.initialization.vel  = md.inversion.vel_obs; % initialization velocity magnitude (m yr^-1)	

	% Initialization pressure
	md.initialization.pressure=md.constants.g*md.materials.rho_ice*md.geometry.thickness; % initial pressure field (Pa)
	md.initialization.temperature=(temp0)*ones(md.mesh.numberofvertices,1); % initial temperature field (K)

	% Basal Friction
   md.friction = frictionschoof(); % Coulomb-limited sliding law used in MISMIP+ (Schoof 2005, 2007, Cornford et. al., 2020)
   md.friction.C = sqrt(3.16E6).*ones(md.mesh.numberofvertices,1); % C^2 = 3.16E6 (Pa m^−(friction.m) s^(friction.m) ) see (Cornford et. al., 2020)
	md.friction.C(md.mask.ocean_levelset<0) = 0; % C^2 = 0 (Pa m^−(friction.m) s^(friction.m) )
	md.friction.Cmax = 0.5.*ones(md.mesh.numberofvertices,1); % Iken's bound 
   md.friction.m = 1/3.*ones(md.mesh.numberofelements,1); % assume to follow Weertman (non-dimensional)
	md.friction.coupling=2;
	
	% Surface Mass Balance from RACMO
   md.smb.mass_balance = interpRACMOant(md.mesh.x,md.mesh.y); % interpolate accumulation rate data from RACMO (SMB_RACMO2.3_1979_2011.nc) (m/yr ice equivalence)

	% Geothermal heat flux
   md.basalforcings.geothermalflux  = interpSeaRISE(md.mesh.x,md.mesh.y,'bheatflx_shapiro',-1); % intperolate geothermal flux from Shapiro et al. (W/m^2)

	% Boundary Conditions
	% reset b.c. on the vertices
   md.stressbalance.spcvx        = NaN(md.mesh.numberofvertices,1);
   md.stressbalance.spcvy        = NaN(md.mesh.numberofvertices,1);
   md.stressbalance.spcvz        = NaN(md.mesh.numberofvertices,1);
   md.stressbalance.referential  = NaN(md.mesh.numberofvertices,6);
   md.stressbalance.loadingforce = zeros(md.mesh.numberofvertices,3);
   md.masstransport.spcthickness = NaN(md.mesh.numberofvertices,1);
	% set dirichlet velocity b.c. on along the edge of the domain
   pos = find((md.mask.ice_levelset<0).*(md.mesh.vertexonboundary)); % the outer domain contour part of the ice
   md.stressbalance.spcvx(pos)   = md.initialization.vx(pos); % set x component b.c. (m yr^-1)
   md.stressbalance.spcvy(pos)   = md.initialization.vy(pos); % set y component b.c. (m yr^-1)
	md.stressbalance.spcvx(isnan(md.stressbalance.spcvx(pos)))=0;
	md.stressbalance.spcvy(isnan(md.stressbalance.spcvy(pos)))=0;

   % Use SSA equations for ice flow
   md=setflowequation(md,'SSA','all');
	% }}}
	savemodel(org,md);
end%}}}
if perform(org,'InversionB'),  % Invert for flow law parameter B{{{
   md=loadmodel(org,'MeshParam');

   % new inversion with M1QN3 package
   md.inversion=m1qn3inversion(md.inversion);

   % Control general
   md.inversion.iscontrol=1; % flag to turn on inversion
   md.inversion.control_parameters={'MaterialsRheologyBbar'}; % invert for B
   md.inversion.maxsteps=50; % maximum number of iterations (gradient computation)
   md.inversion.maxiter=50;  % maximum number of function evaluations (forward run)
   md.inversion.dxmin=0.1;   % convergence criterion: two points less than dxmin from each other (sup-norm) are considered identical
   md.inversion.gttol=1E-6;  % convergence criterion: ||g(X)||/||g(X0)|| (g(X0): gradient at initial guess X0)
   md.inversion.incomplete_adjoint=0; % 0 non linear viscosity; 1 linear viscosity 

   % Cost functions
	% 101: SurfaceAbsVelMisfit (fit in linear space)
	% 103: SurfaceLogVelMisfit (fit in log space)
	% 502: RheologyBbarAbsGradient (regularization)
   md.inversion.cost_functions=[101 103 502]; % set the cost functions
   md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,length(md.inversion.cost_functions)); % cost function coefficient at every node
	
	% Cost function coefficients: cost functions 101 and 103 should have about the same contribution at the end of the inversion
   md.inversion.cost_functions_coefficients(:,1) = 4.5; % coefficient for linear fit
   md.inversion.cost_functions_coefficients(:,2) = 1; % coefficient for log space fit (always 1)
   md.inversion.cost_functions_coefficients(:,3) = 0.5E-20; % coefficient for regularization
	load(fullfile(modeldir, 'pos_vel_nan.mat'),'pos_vel_nan'); % load locations of NaN data
   md.inversion.cost_functions_coefficients(pos_vel_nan,:)=0; % do not try to fit positions with NaN in the velocity data set

	% Controls on max/min B allowed: see Truffer & Stanton, 2015 for englacial ice temperature bounds, Cuffey & Paterson, 2010 for A and E*
   % min: A=2.4E-24 for   0degC ice, E*=10  for max enhancement factor of Antarctic ice in strong shear
   % max: A=6.8E-26 for -25degC ice, E*=0.6 for min enhancement factor of Antarctic ice shelves
   lim_A = [6.8E-26, 2.4E-24];                               % limits on A (Pa^-n s^-1)
   lim_Estar = [0.6, 10];                                    % limits on enhancement factor E* (non-dim.)
   lim_B = (lim_A.*lim_Estar).^(-1/md.materials.rheology_n); % limits on B (Pa s^1/n)
   md.inversion.min_parameters= min(lim_B)*ones(size(md.materials.rheology_B)); % lower bound on B (Pa s^1/n)
   md.inversion.max_parameters= max(lim_B)*ones(size(md.materials.rheology_B)); % upper bound on B (Pa s^1/n)
   md.inversion.min_parameters(md.mask.ocean_levelset>0) = md.materials.rheology_B(md.mask.ocean_levelset>0); % do not allow B to change on grounded ice
   md.inversion.max_parameters(md.mask.ocean_levelset>0) = md.materials.rheology_B(md.mask.ocean_levelset>0); % do not allow B to change on grounded ice

   % Stress balance parameters
   md.stressbalance.maxiter=50;   %
   md.stressbalance.reltol=NaN;   %
   md.stressbalance.abstol=NaN;   %

	% Extract the floating model subdomain and prepare to solve
	% This sets Dirichlet velocity b.c. at the grounding line
   mds=extract(md,md.mask.ocean_levelset<0);
   mds.cluster=generic('name',oshostname(),'np',45); %for totten 45 ideal
   mds.verbose=verbose('solution',false,'control',true);
   mds.miscellaneous.name='inversion_B';
   mds.friction.C(:)=0; % make sure there is no basal friction

   % Solver parameters
   mds.toolkits.DefaultAnalysis=bcgslbjacobioptions();% biconjugate gradient with block Jacobi preconditioner
   mds.toolkits.DefaultAnalysis.ksp_max_it=500;
   mds.settings.solver_residue_threshold=1e-6;

	% Solve
   mds=solve(mds,'Stressbalance'); % only extracted model

   % Update the full model rheology_B accordingly
   md.materials.rheology_B(mds.mesh.extractedvertices) = mds.results.StressbalanceSolution.MaterialsRheologyBbar;
   savemodel(org,md);
end % }}}
if perform(org,'InversionC'),  % Invert for friction coefficient C {{{
   md=loadmodel(org,'InversionB');

	% new inversion with M1QN3 package
   md.inversion=m1qn3inversion(md.inversion);

	% Control general
   md.inversion.iscontrol=1; % flag to turn on inversion
   md.inversion.control_parameters={'FrictionC'}; % invert for C
   md.inversion.maxsteps=50; % maximum number of iterations (gradient computation)
   md.inversion.maxiter=50;  % maximum number of function evaluations (forward run)
   md.inversion.dxmin=0.1;   % convergence criterion: two points less than dxmin from each other (sup-norm) are considered identical
   md.inversion.gttol=1E-6;  % convergence criterion: ||g(X)||/||g(X0)|| (g(X0): gradient at initial guess X0)
   md.inversion.incomplete_adjoint=0; % 0 non linear viscosity; 1 linear viscosity

   % Cost functions
   % 101: SurfaceAbsVelMisfit (fit in linear space)
   % 103: SurfaceLogVelMisfit (fit in log space)
	% 502: FrictionCAbsGradient (regularization)
   md.inversion.cost_functions=[101 103 501]; % set the cost functions
   md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,length(md.inversion.cost_functions)); % cost function coefficient at every node

   % Cost function coefficients: cost functions 101 and 103 should have about the same contribution at the end of the inversion
   md.inversion.cost_functions_coefficients(:,1) = 600; % coefficient for linear fit
   md.inversion.cost_functions_coefficients(:,2) = 1; % coefficient for log space fit (always 1)
   md.inversion.cost_functions_coefficients(:,3) = 5E-10; % coefficient for regularization
   load(fullfile(modeldir, 'pos_vel_nan.mat'),'pos_vel_nan'); % load locations of NaN data
   md.inversion.cost_functions_coefficients(pos_vel_nan,:)=0; % do not try to fit positions with NaN in the velocity data set

   % Controls on max/min C allowed
   md.inversion.min_parameters = 0*md.friction.C;
   md.inversion.max_parameters =  4*md.friction.C;
   % Keep basal friction constant and 0 in the floating part
   md.inversion.min_parameters(md.mask.ocean_levelset<0)=0;
   md.inversion.max_parameters(md.mask.ocean_levelset<0)=0;
   
	% Stress balance parameters
   md.stressbalance.maxiter=50;   %
   md.stressbalance.reltol=NaN;   %
   md.stressbalance.abstol=NaN;   %

	% Prepare to solve
   md.cluster=generic('name',oshostname(),'np',55);
   md.verbose=verbose('solution',false,'control',true);
	md.miscellaneous.name='inversion_friction_C';

   % Solver parameters
   md.toolkits.DefaultAnalysis=bcgslbjacobioptions();% biconjugate gradient with block Jacobi preconditioner
   md.toolkits.DefaultAnalysis.ksp_max_it=500;
   md.settings.solver_residue_threshold=1e-5;

   % Solve inversion for C
   md=solve(md,'Stressbalance'); 

   % Update friction coefficient and velocity field
   md.friction.C=md.results.StressbalanceSolution.FrictionC;
   md.initialization.vx= md.results.StressbalanceSolution.Vx;
   md.initialization.vy=md.results.StressbalanceSolution.Vy;
   md.initialization.vel=md.results.StressbalanceSolution.Vel;

   % Solve for the converged modeled velocity field 
   md.inversion.iscontrol=0;
   md.verbose.convergence = 1;
   md = solve(md,'sb');
   % Update velocity field
   md.initialization.vx= md.results.StressbalanceSolution.Vx;
   md.initialization.vy=md.results.StressbalanceSolution.Vy;
   md.initialization.vel=md.results.StressbalanceSolution.Vel;

   savemodel(org,md);
end%}}}
if perform(org,'InversionB2'),  % Re-invert for flow law parameter B{{{
   md=loadmodel(org,'InversionC');

   % new inversion with M1QN3 package
   md.inversion=m1qn3inversion(md.inversion);

   % Control general
   md.inversion.iscontrol=1; % flag to turn on inversion
   md.inversion.control_parameters={'MaterialsRheologyBbar'}; % invert for B
   md.inversion.maxsteps=50; % maximum number of iterations (gradient computation)
   md.inversion.maxiter=50;  % maximum number of function evaluations (forward run)
   md.inversion.dxmin=0.1;   % convergence criterion: two points less than dxmin from each other (sup-norm) are considered identical
   md.inversion.gttol=1E-6;  % convergence criterion: ||g(X)||/||g(X0)|| (g(X0): gradient at initial guess X0)
   md.inversion.incomplete_adjoint=0; % 0 non linear viscosity; 1 linear viscosity 

   % Cost functions
	% 101: SurfaceAbsVelMisfit (fit in linear space)
	% 103: SurfaceLogVelMisfit (fit in log space)
	% 502: RheologyBbarAbsGradient (regularization)
   md.inversion.cost_functions=[101 103 502]; % set the cost functions
   md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,length(md.inversion.cost_functions)); % cost function coefficient at every node
	
	% Cost function coefficients: cost functions 101 and 103 should have about the same contribution at the end of the inversion
   md.inversion.cost_functions_coefficients(:,1) = 200; % coefficient for linear fit
   md.inversion.cost_functions_coefficients(:,2) = 1; % coefficient for log space fit (always 1)
   md.inversion.cost_functions_coefficients(:,3) = 1E-16; % coefficient for regularization
	load(fullfile(modeldir, 'pos_vel_nan.mat'),'pos_vel_nan'); % load locations of NaN data
   md.inversion.cost_functions_coefficients(pos_vel_nan,:)=0; % do not try to fit positions with NaN in the velocity data set

   % Controls on max/min B allowed: see Truffer & Stanton, 2015 for englacial ice temperature bounds, Cuffey & Paterson, 2010 for A and E*
	% min: A=2.4E-24 for   0degC ice, E*=10  for max enhancement factor of Antarctic ice in strong shear
	% max: A=6.8E-26 for -25degC ice, E*=0.6 for min enhancement factor of Antarctic ice shelves
	lim_A = [6.8E-26, 2.4E-24];                               % limits on A (Pa^-n s^-1)
	lim_Estar = [0.6, 10];                                    % limits on enhancement factor E* (non-dim.)
	lim_B = (lim_A.*lim_Estar).^(-1/md.materials.rheology_n); % limits on B (Pa s^1/n)
	md.inversion.min_parameters= min(lim_B)*ones(size(md.materials.rheology_B)); % lower bound on B (Pa s^1/n)
   md.inversion.max_parameters= max(lim_B)*ones(size(md.materials.rheology_B)); % upper bound on B (Pa s^1/n)
	md.inversion.min_parameters(md.mask.ocean_levelset>0) = md.materials.rheology_B(md.mask.ocean_levelset>0); % do not allow B to change on grounded ice
	md.inversion.max_parameters(md.mask.ocean_levelset>0) = md.materials.rheology_B(md.mask.ocean_levelset>0); % do not allow B to change on grounded ice

   % Stress balance parameters
   md.stressbalance.maxiter=50;   %
   md.stressbalance.reltol=NaN;   %
   md.stressbalance.abstol=NaN;   %

   md.cluster=generic('name',oshostname(),'np',45); %for totten 45 ideal
   md.verbose=verbose('solution',false,'control',true);
   md.miscellaneous.name='inversion_B2';

   % Solver parameters
   md.toolkits.DefaultAnalysis=bcgslbjacobioptions();% biconjugate gradient with block Jacobi preconditioner
   md.toolkits.DefaultAnalysis.ksp_max_it=500;
   md.settings.solver_residue_threshold=1e-6;

	% Solve
   md=solve(md,'Stressbalance'); % only extracted model

   % Update the full model rheology_B accordingly
   md.materials.rheology_B = md.results.StressbalanceSolution.MaterialsRheologyBbar;
   savemodel(org,md);
end % }}}

% local functions {{{
% file writing
function write_expfile(fname,EXP) % {{{
	% WRITE_EXPFILE writes an exp domain outline from EXP to fname
	% INPUT: fname    .exp file to be written to
	%        EXP      structure with fields (x,y) 
	%          x      field with x outline values
	%          y      field with y outline values

	disp(['writing ISSM exp file to ' fname])
	fileID = fopen(fname,'w');
	fprintf(fileID,'## Name:domainoutline\n');
	fprintf(fileID,'## Icon:0\n');
	fprintf(fileID,'# Points Count  Value\n');
	fprintf(fileID,'%i 1.0\n',length(EXP.x));
	fprintf(fileID,'# X pos Y pos\n');
	fprintf(fileID,'%f %f\n',[EXP.x; EXP.y]);
	fclose(fileID);
end 
% }}}
% }}}
