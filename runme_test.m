% PROPHET Amundsen Sea Coupling
% Outline:
%	Initialize and spin up ocean model
%  Initialize ice sheet model
%  Run coupled 
%
% Experiment:
%  Run with monthly climatology for ocean boundary conditions
%    Each month is different but gets repeated every year
%    The calendar is a "model" calendar:
%      12 months of 30 days = 360 days per year
%      Either 

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
% }}}

% Ice Model
expfile='./Exp/domain.exp';
% ISSM domain setup

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
	temp0 = 273.15-10; % approzimate temperature of the ice (for A and initialization) (deg K)
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
	pos = find(max(md.mask.ice_levelset(md.mesh.elements),[],2)>0); % elements with at least 1 vertex that has no ice
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
	md.friction = frictioncoulomb2(); % Coulomb-limited sliding law used in MISMIP+ (Cornford et. al., 2020)
	md.friction.C = sqrt(3.16E6).*ones(md.mesh.numberofelements,1); % C^2 = 3.16E6 (Pa m^−(friction.m) s^(friction.m) ) see (Cornford et. al., 2020)
	md.friction.m = 1/3.*ones(md.mesh.numberofelements,1); % assume to follow Weertman (non-dimensional)
	
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
   md.inversion.cost_functions_coefficients(:,3) = 1E-18; % coefficient for regularization
	load(fullfile(modeldir, 'pos_vel_nan.mat'),'pos_vel_nan'); % load locations of NaN data
   md.inversion.cost_functions_coefficients(pos_vel_nan,:)=0; % do not try to fit positions with NaN in the velocity data set

   % Controls on max/min B allowed
   md.inversion.min_parameters=0.8*cuffey(273.15)*ones(size(md.materials.rheology_B)); % from Seroussi et al, 2014
   md.inversion.max_parameters=1*cuffey(273.15-30)*ones(size(md.materials.rheology_B)); % from Seroussi et al, 2014

   % Stress balance parameters
   md.stressbalance.maxiter=50;   %
   md.stressbalance.reltol=NaN;   %
   md.stressbalance.abstol=NaN;   %

	% Extract the floating model subdomain and prepare to solve
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

	% Basal Friction
   md.friction = frictionschoof(); % Coulomb-limited sliding law used in MISMIP+ (Cornford et. al., 2020)
   md.friction.C = sqrt(3.16E6).*ones(md.mesh.numberofvertices,1); % C^2 = 3.16E6 (Pa m^−(friction.m) s^(friction.m) ) see (Cornford et. al., 2020)
	md.friction.Cmax = 0.5.*ones(md.mesh.numberofvertices,1);
   md.friction.m = 1/3.*ones(md.mesh.numberofelements,1); % assume to follow Weertman (non-dimensional)
	md.friction.coupling=2;

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
   md.inversion.cost_functions_coefficients(:,1) = 300; % coefficient for linear fit
   md.inversion.cost_functions_coefficients(:,2) = 1; % coefficient for log space fit (always 1)
   md.inversion.cost_functions_coefficients(:,3) = 1E-16; % coefficient for regularization
   load(fullfile(modeldir, 'pos_vel_nan.mat'),'pos_vel_nan'); % load locations of NaN data
   md.inversion.cost_functions_coefficients(pos_vel_nan,:)=0; % do not try to fit positions with NaN in the velocity data set

   % Controls on max/min C allowed
   md.inversion.min_parameters = 0*md.friction.C;
   md.inversion.max_parameters =  4*md.friction.C;
   % Keep basal friction constant in the floating part
   md.inversion.min_parameters(md.mask.ocean_levelset<0)=md.friction.C(md.mask.ocean_levelset<0);
   md.inversion.max_parameters(md.mask.ocean_levelset<0)=md.friction.C(md.mask.ocean_levelset<0);
   
	% Stress balance parameters
   md.stressbalance.maxiter=50;   %
   md.stressbalance.reltol=NaN;   %
   md.stressbalance.abstol=NaN;   %

	% Prepare to solve

   md.cluster=generic('name',oshostname(),'np',55); %for totten 45 ideal
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


% local functions {{{
% file writing
function write_sizefile(fname,SZ) % {{{
	% Reads from SZ.reffile and writes to fname
	% INPUT
	%    fname   file to write to 
	%    SZ      struct with fields: sNx,sNy,OLx,OLy,nSx,nSy,nPx,nPy,Nx,Ny,Nr, and reffile
	%     SZ.reffile   MITgcm reference file to read from

	if SZ.Nx~=(SZ.sNx*SZ.nSx*SZ.nPx) | SZ.Ny~=(SZ.sNy*SZ.nSy*SZ.nPy)
		error('MITgcm domain discretization inconsistent');
	end

	disp(['writing SIZE     file to ' fname]);
	writeID=fopen(fname,'w');
	readID=fopen(SZ.reffile,'r');
	% read through the template file, write to the new file
	formatSpec='%s\n'; % new line after each string is written
	values=[SZ.sNx, SZ.sNy, SZ.OLx, SZ.OLy, SZ.nSx, SZ.nSy, SZ.nPx, SZ.nPy, SZ.Nx, SZ.Ny, SZ.Nr]; % ensure correct ordering of values
	% read through any uncommented header, do NOT write to new file {{{
	tline = fgetl(readID);
	while ~strcmp(tline,'CBOP')
		tline = fgetl(readID);
	end % }}}
	% read through the commented header, write to new file {{{
	while tline(1)=='C'
		fprintf(writeID,formatSpec,tline);
		tline = fgetl(readID);
	end %}}}
	% read through the variable declarations, write to new file {{{
	while contains(tline,'INTEGER') | contains(tline,'PARAMETER')
		fprintf(writeID,formatSpec,tline);
		tline = fgetl(readID);
	end % }}}
	% read through the variable values, write to new file {{{
	i=1;
	while contains(tline,'&')
		% assumes Nx  = sNx*nSx*nPx and Ny  = sNy*nSy*nPy
		if ~any(i==[9,10])
			% extract the string of char before and after the template variable value
			[tempvalue,sline] = regexp(tline,'\d*','Match','split');
			% insert the variable value
			tline=[sline{1} num2str(values(i)) sline{2}];
		end
		% write to new file
		fprintf(writeID,formatSpec,tline);
		% advance to the next line
		i=i+1;
		tline = fgetl(readID);
	end % }}}
	% read the rest of the template file, write to new file (assumes MAX_OLX = OLx, and MAX_OLY = OLy) {{{
	while isstr(tline)
		fprintf(writeID,formatSpec,tline);
		tline = fgetl(readID);
	end % }}}
	fclose(writeID);
	fclose(readID);
end % }}}
function write_pkgconffile(fname,PKGCONF) % {{{
	% WRITE_PKGCONFFILE writes the requested package names to the MITgcm compile-time configuration file
	disp(['writing config.  file to ' fname]);
	fileID = fopen(fname,'w');
	fprintf(fileID,'# %s\n',PKGCONF.description); % write description
	fprintf(fileID,'%s\n',PKGCONF.pkg{:}); % write packages
	fclose(fileID);
end % }}}
function write_datafile(fname,C,head) % {{{
	% WRITE_DATAFILE writes structures in C to fname {{{
	% INPUT: fname   string file to write
	%        C       cell array of structures (P) to write
	%        head    string of file header
	%
	% OUTPUT: writes to file with the following form:
	%
	% # head
	% # C{1}.description
	%  &C{1}.header
	%  C{1}.field1=value1,
	%  C{1}.field2=value2,
	%  ...
	%  & 
	%
	% # C{2}.description
	%  &C{2}.header
	%  C{2}.field1=value1,
	%  C{2}.field2=value2,
	%  ...
	%  &
	%  ...
	% }}}
	disp(['writing namelist file to ' fname])
	fileID = fopen(fname,'w');
	fprintf(fileID,'# %s\n',head); % write descriptive file header 

	% loop through structures in C
	for i=1:length(C)
		% write a descriptive comment if it exists {{{
		if isfield(C{i},'description')
			if iscell(C{i}.description)
				fprintf(fileID,'# %s\n',C{i}.description{:}); % write multi-line description
			else
				fprintf(fileID,'# %s\n',C{i}.description); % write description
			end
			C{i}=rmfield(C{i},'description'); % do not write as field
		end % }}}
		% write PARM header {{{
		fprintf(fileID,' &%s\n',C{i}.header); % write header
		C{i}=rmfield(C{i},'header'); % do not write as field
		% }}}
		% write parameter fields and end section {{{
		if isfield(C{i},'N') % if diagnostic fields
			writediagfields(fileID,C{i}.N);
			C{i}=rmfield(C{i},'N'); % done with diag fields
		end
		writefields(fileID,C{i}); % write non-diagnostic fields
		fprintf(fileID,' &\n\n'); % end section
		% }}}
	end 
	fclose(fileID);
end % }}}
function writefields(fileID,P) % {{{
	% WRITEFIELDS writes each field in P to fileID as a new line
	%  Each line takes the form ' fieldname=value,'
	fields=fieldnames(P);
	for i=1:length(fields)
		val=getfield(P,fields{i});
		fprintf(fileID,' %s=%s,\n',fields{i},num2str(val));
	end
end  % }}}
function writediagfields(fileID,N) % {{{
	% WRITEDIAGFIELDS writes diagnostic fields for each output stream
	% INPUT: N is a struct. array with an element for each output stream n
	%  Each line takes the form fieldname(n)=value except for 
	%  'fields' with take the form fields(1:length,n)='fields{1}', 'fields{2}', ...

	% loop over each output stream
	for n=1:length(N)
		subfields=fieldnames(N(n));
		% loop over each subfield
		for i=1:length(subfields)
			switch subfields{i}
				case 'fields'
					LHS=[subfields{i} '(1:' num2str(length(N(n).fields)) ',' num2str(n) ')'];
					dfields=getfield(N(n),subfields{i});
					dfields=strcat('''',dfields,''', ');
					RHS=strcat(dfields{:});
					fprintf(fileID,' %s=%s\n',LHS,RHS); % write to file
				case 'levels'
					error('diag. levels not supported');
				otherwise
					LHS=[subfields{i} '(' num2str(n) ')'];
					RHS=num2str(getfield(N(n),subfields{i}));
					fprintf(fileID,' %s=%s,\n',LHS,RHS); % write to file
			end
		end
		fprintf(fileID,'\n'); % line break
	end
end % }}}
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
