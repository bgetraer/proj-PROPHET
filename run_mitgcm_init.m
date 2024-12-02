steps=[10];

% PROPHET Amundsen Sea Coupling
% Outline:
%	Initialize MITgcm
%  Compile MITgcm
%  Run MITgcm
%
% Experiments:
%   'RCP85'   % ensemble average forcings from ISMIP-6
%   'Paris2C' % ensemble average forcings from Paris 2
%   'CLIM'    % monthly climatology 2001-2012 from Paris 2

experiment.name='Paris2C';
experiment.init='MITgcm_initialization';
mit_dT      = 100; % s
coupling_dT = 15*24*60*60; % s
% directory structure {{{
mitgcm_dir='/nobackup/bgetraer/MITgcm'; % MITgcm directory (pleaides)
%mitgcm_dir='/totten_1/bgetraer/MITgcm'; % MITgcm directory (totten)
proph_dir = pwd; % base directory for this project
% define experiment directories {{{
% make experiments directory if needed
% this will hold subdirectories for each model run
if ~exist(fullfile(proph_dir,'experiments'))
	mkdir(fullfile(proph_dir,'experiments'));
end
% make the initialization experiment directory if needed
initdir=fullfile(proph_dir,'experiments',experiment.init);
if ~exist(initdir)
	mkdir(initdir);
end
% make code directory if needed
% this will hold the current compilation of MITgcm
if ~exist(fullfile(initdir,'code'))
	mkdir(fullfile(initdir,'code'));
end
% make input directory if needed
% this will hold the runtime options for MITgcm
if ~exist(fullfile(initdir,'input'))
	mkdir(fullfile(initdir,'input'));
end
% make model directory if needed
% this will hold md structures from ISSM
modeldir=fullfile(initdir,'Models');
if ~exist(modeldir)
	mkdir(modeldir);
end
% make experiment directory if needed
% this is where the different runs will occur
expdir=fullfile(proph_dir,'experiments',experiment.name);
if ~exist(expdir)
   mkdir(expdir);
end
% }}}
% check for climate forcing directory {{{
% make climateforcings directory if needed
% this will hold subdirectories for the climate forcing data
if ~exist(fullfile(proph_dir,'climateforcings'))
	warning('no climateforcings directory, making one now');
	mkdir(fullfile(proph_dir,'climateforcings'));
end
% }}}
% }}}
% suppress warnings {{{
wid = {'MATLAB:hg:AutoSoftwareOpenGL','MATLAB:polyshape:repairedBySimplify'};
for i=1:length(wid)
	warning('off',wid{i});
end
% }}}
org=organizer('repository',modeldir,'prefix','PROPHET_mitgcm_init_','steps',steps);
% Move into the experiment directory
disp(['Moving to experiment directory: ', initdir]);
cd(initdir); 

% These steps are the same for each experiment
if perform(org,'MeshInit'), % {{{	
	mit = struct(); % initialize model parameter structure
	% Directory and filename information
	% climate forcing filenames {{{
	% There are two experiments from Naughten et al., 2023 (Paris2C and RCP85), each with 10 members.
	% There are two types of files: 
	%    bdry files which give the ocean state interpolated in x and y at the bottom and left
	%       boundaries, for each month between 2010--2100
	%    init files which give the ocean state interpolated in x and y across the entire domain for 
	%       the initial state of Jan 2010.
	%  Files have been interpolated from the original output data by Dan Goldberg.
	mit.forcing = struct();
	mit.forcing.Ndir = fullfile(proph_dir,'climateforcings/naughten2023_data'); % directory path for unprocessed Naughten data
	mit.forcing.bdry_prefix = 'bdryDatapaholPy'; % prefix for boundary forcing files
	mit.forcing.init_prefix = 'initpaholPy'; % prefix for initial state files
	mit.forcing.exp = {'Paris2C', 'RCP85'}; % string for climate experiment 
	mit.forcing.member = sprintfc('%03.0f',[1:10]); % string for member
	mit.forcing.suffix = 'benjy.mat'; % suffix for all files
	mit.forcing.bdry_example = [mit.forcing.bdry_prefix mit.forcing.exp{1} mit.forcing.member{1} mit.forcing.suffix]; % one of the bdry files
	mit.forcing.init_example = [mit.forcing.init_prefix mit.forcing.exp{1} mit.forcing.member{1} mit.forcing.suffix]; % one of the init files
	mit.forcing.field_bdry = {'Tleftbdry','Sleftbdry','Uleftbdry','Vleftbdry','Tbotbdry','Sbotbdry','Ubotbdry','Vbotbdry'};  % field names for bdry files
	mit.forcing.field_init = {'Tinit','Sinit'}; % field names for init files
	mit.forcing.Ddir = fullfile(proph_dir,'climateforcings/data'); % directory path for processed boundary forcing data
	mit.forcing.field_bdry_out = {'theta.obw','salt.obw','uvel.obw','vvel.obw','theta.obs','salt.obs','uvel.obs','vvel.obs'}; % data fields
	mit.forcing.field_init_out = {'theta.init','salt.init'}; % data fields
	% }}}
	% MITgcm filenames {{{
	mit.fname = struct();
	% compile-time files
	mit.fname.sizefile='code/SIZE.h';
	mit.fname.pkgconffile='code/packages.conf';
	mit.fname.diagsizefile='code/DIAGNOSTICS_SIZE.h';
	mit.fname.obcsoptionsfile='code/OBCS_OPTIONS.h';
	mit.fname.cppoptionsfile='code/CPP_OPTIONS.h';
	% run-time namelist files
	mit.fname.datafile='input/data';
	mit.fname.dataobcsfile='input/data.obcs';
	mit.fname.datashelficefile='input/data.shelfice';
	mit.fname.datacalfile='input/data.cal';
	mit.fname.dataexffile='input/data.exf';
	mit.fname.datadiagfile='input/data.diagnostics';
	mit.fname.datapkgfile='input/data.pkg';
	mit.fname.eedatafile='input/eedata';
	% binary input files
	mit.fname.treffile='tref.bin';
	mit.fname.sreffile='sref.bin';
	mit.fname.delrfile = 'delr.bin';
	mit.fname.bathyfile='bathy.bin';
	% binary OBCS files  WITHOUT THE EXPERIMENT SUFFIX (i.e. Paris2, RCP85)
	mit.fname.thetaOBWfile= ['theta.obw.' ];
	mit.fname.saltOBWfile = ['salt.obw.'  ];
	mit.fname.uvelOBWfile = ['uvel.obw.'  ];
	mit.fname.vvelOBWfile = ['vvel.obw.'  ];
	mit.fname.thetaOBSfile= ['theta.obs.' ];
	mit.fname.saltOBSfile = ['salt.obs.'  ];
	mit.fname.uvelOBSfile = ['uvel.obs.'  ];
	mit.fname.vvelOBSfile = ['vvel.obs.'  ];
	% binary SHELFICE files
	mit.fname.draftfile='draft.bin';
	% }}}
	% Generate mesh
	% Vertical discretization {{{
	% The vertical grid is defined by the "cell centered" approach (see MITgcm readthedocs 2.11.5) and is defined by the vector delzF
	% which provides the thickness of the cells. Here, the vertical grid is non-constant, varying from 10m at the surface to a maximum
	% defined by delzF_max at depth. Between these cell thicknesses, the grid spacing is lifted directly from Naughten et al., 2023.
	fname = fullfile(mit.forcing.Ndir, mit.forcing.init_example);
	zc_N = getfield(load(fname,'z'),'z');     % Naughten2023 cell center locations in z (m)
	zp_N = zeros(size(zc_N));                 % Initialize vertical cell edge locations in z (m)
	for i=1:length(zp_N)
		zp_N(i+1) = 2*zc_N(i) - zp_N(i);       % Naughten2023 vertical cell edge locations in z (m)
	end
	delzF_N = abs(diff(zp_N));                     % Naughten2023 vertical cell thickness (m)
	delzC_N = abs(diff([zp_N(1) zc_N zp_N(end)])); % Naughten2023 vertical cell center spacing (m)

	% Here, we cap the cell thickness at delzF_max, and use the same cell spacing as Naughten up
	% until that limit, followed by cells of uniform delzF_max thickness until we reach the lowest
	% bathymetry we need to capture.
	zmin = -2.067E3; % minimum bathymetry, from B_dagger, which actually is calculated later.(m)
	% Calculate the vertical cell thickness
	delzF_max = 32.5; % max vertical spacing (m)
	delzF_upper = delzF_N(delzF_N<=delzF_max); % use the same vertical spacing as Naughten in upper column (m)
	delzF_lower = repmat(delzF_max,1,abs(round((zmin/delzF_max)))); % use constant spacing of delzF_max in the upper column (m)
	delzF = [delzF_upper delzF_lower]; % vertical cell thickness (m)
	% Calculate the vertical cell edge locations and center spacing
	zp = [0 cumsum(-delzF)];
	z_end = find(zp<zmin,1); % index to terminate the vertical grid
	zp = zp(1:z_end);                 % vertical cell edge locations in z (m)
	delzF = delzF(1:z_end-1);         % vertical cell thickness in z (m)
	zc = zp(1:end-1) - 0.5*delzF;     % cell center locations in z (m)
	% Calculate the vertical cell center spacing
	delzC = abs(diff([zp(1) zc zp(end)])); % vertical cell center spacing (m)
	% }}}
	% Horizontal discretization {{{
	% A cartesian coordinate system is used with x and y following the coordinates of the Polar Stereographic projection with a standard
	% parallel of 71 deg S.
	% At the Amundsen sea, the projection places North roughly in the -x direcion, East in +y, South in +x, West in -y.
	% Where MITgcm convention denotes a "north" or "south" (for instance in the OBCS package) these actually refer to Cartesisan 
	% +y and -y respectively. Because this model implements an f-plane approximation for the Coriolis force, no geographic sense of north
	% or south is nessecary, and the domain generally follows a "North up" "y up" convention, despite the fact that in our implementation 
	% that direction is not geographic north.

	% ocean domain
	x0=-1700000; % x origin of the active ocean domain (m)
	y0=-710000;  % y origin of the active ocean domain (m)
	Lx=300E3; % length of total ocean domain in x (m) 
	Ly=504E3; % length of total ocean domain in y (m)
	%Lz=-2075;  % height of total ocean domain in z (m)

	dx=1.2E3; % hor. resolution in x (m)
	dy=dx;    % hor. resolution in y (m)
	%dz=25;    % ver. resolution in z (m)

	% SETUP:
%	%	5 nodes, 140 processes, 10x14 tiles, 25x30 cells, 300E3x505E3 m domain
%	sNx=25;  % Number of X points in tile
%	sNy=30;  % Number of Y points in tile

%	%  5 nodes, 70 processes, 5x14 tiles, 50x30 cells
%	sNx=50;  % Number of X points in tile
%	sNy=30;  % Number of Y points in tile

%	%	10 nodes, 280 processes, 10x28 tiles, 25x15 cells, 300E3x505E3 m domain
%	sNx=25;  % Number of X points in tile
%	sNy=15;  % Number of Y points in tile

%	%	4 nodes, 100 processes, 10x10 tiles, 25x42 cells, 300E3x505E3 m domain
%	sNx=25;  % Number of X points in tile
%	sNy=42;  % Number of Y points in tile

%	%	6 nodes, 150 processes, 10x15 tiles, 25x28 cells, 300E3x505E3 m domain
%	sNx=25;  % Number of X points in tile
%	sNy=28;  % Number of Y points in tile

	%	7 nodes, 175 processes, 5x35 tiles, 50x12 cells, 300E3x505E3 m domain
	sNx=50;  % Number of X points in tile
	sNy=12;  % Number of Y points in tile

	Nx=Lx/dx;      % number of cells in x
	Ny=Ly/dy;      % number of cells in y
	%Nz=Lz/dz;      % number of total cells in z
	Nz=length(zc); % number of cells in z

	xp=((0:Nx)*dx)+x0; % location of all cell edges in x (m)
	yp=((0:Ny)*dy)+y0; % location of all cell edges in y (m)
	%zp=(0:Nz)*dz;      % location of all cell edges in z (m)
	zg=zp(1:end-1);   % location of z edge points, not including the lowest (redundant) boundary (m)

	xc=xp(1:end-1)+0.5*dx; % location of cell centers in x (m)
	yc=yp(1:end-1)+0.5*dy; % location of cell centers in y (m)
	%zc=zp(1:end-1)+0.5*dz; % location of cell centers in z (m)
	% }}}
	% Save mesh to mit structure {{{
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PROPHET mesh
	mit.mesh = struct();
	mit.mesh.delxF = dx; % horizontal spacing between cell faces in x (m)
	mit.mesh.delyF = dy; % horizontal spacing between cell faces in y (m)
	mit.mesh.delzF = delzF; % vertical spacing between cell faces (m)
	mit.mesh.delzC = delzC; % vertical spacing between cell centers (m)
	mit.mesh.xp = xp; % cell face locations in x (m)
	mit.mesh.yp = yp; % cell face locations in y (m)
	mit.mesh.zp = zp; % cell face locations in z (m)
	mit.mesh.xc = xc; % cell center locations in x (m)
	mit.mesh.yc = yc; % cell center locations in y (m)
	mit.mesh.zc = zc; % cell center locations in z (m)
	[mit.mesh.XC, mit.mesh.YC, mit.mesh.ZC] = meshgrid(xc,yc,zc); % all cell center locations in x, y, and z (m)
	[mit.mesh.hXC, mit.mesh.hYC] = meshgrid(xc,yc);	% horizontal cell center locations in x and y (m)
	mit.mesh.Nx = numel(xc); % number of cells in x dim
	mit.mesh.Ny = numel(yc); % number of cells in y dim
	mit.mesh.Nz = numel(zc); % number of cells in z dim
	% PROPHET tiling
	mit.mesh.sNx = sNx; % Number of X points in tile
	mit.mesh.sNy = sNy; % Number of Y points in tile
	% Naughten mesh
	mit.mesh.zc_N = zc_N; % Naughten2023 cell center locations in z (m)
	mit.mesh.zp_N = zp_N; % Naughten2023 cell face locations in z (m)
	mit.mesh.delzF_N = delzF_N; % Naughten2023 vertical spacing between cell faces (m)
	mit.mesh.delzC_N = delzC_N; % Naughten2023 vertical spacing between cell centers (m)
	[mit.mesh.XC_N,mit.mesh.YC_N,mit.mesh.ZC_N]=meshgrid(xc,yc,zc_N);  % Naughten2023 cell center locations in x, y, and z (m)
	% write delRfile
	fname = ['input/' mit.fname.delrfile];
	disp(['writing mit.mesh.delzF to ' fname]);
	write_binfile(fname,mit.mesh.delzF);
	% }}}
	% Run-time options
	% input/eedata {{{
	EEP=struct; % initialize EEDATA structure
	% EEPARMS: Execution Environment Parameters {{{
	EEP.header='EEPARMS';
	% }}}
	mit.inputdata.EEP={EEP};
	% }}}
	% input/data {{{
	P1=struct;P2=struct;P3=struct;P4=struct;P5=struct; % initialize PARM structures
	% PARM01: Continuous equation parameters {{{
	% structure information
	P1.header='PARM01';
	P1.description='Continuous equation parameters';

	% physical constants
	P1.rhoConst = 1030; % vertically constant reference density (Boussinesq) (kg/m3)
	P1.gravity=9.81; % gravitational acceleration (m/s2)

	% equation of state
	P1.eostype='''JMD95Z''';
	P1.tRefFile = ['''' mit.fname.treffile '''']; % filename for vertical profile of init. and ref. salin. (g/kg)
	P1.sRefFile = ['''' mit.fname.sreffile '''']; % filename for vertical profile of init. and ref. potential temp. (deg C)
	P1.rhoNil=1000; % ref. density for linear EOS (kg/m^3)
	P1.tAlpha = 3.733E-5; % thermal expansion coefficient (deg C)^-1
	P1.sBeta  = 7.8434E-4; % haline expansion coefficient (g/kg)^-1
	P1.HeatCapacity_cp = 3974.0; % specific heat capacity Cp (ocean) (J/kg/K)
	% coriolis parameters
	P1.selectCoriMap=0; % 0=f-plane; 1=b-plane
	lat0=-75.0; % reference latitude (deg)
	phi0=lat0/360*2*pi; % reference latitude (rad)
	rotationPeriod=8.6164E+04; % MITgcm default (s)
	omega=2*pi/rotationPeriod; % angular velocity of Earth (rad/s)
	P1.f0=round(2*omega*sin(phi0),2,'significant'); % reference Coriolis parameter (1/s)
	if P1.selectCoriMap==1 % beta-plane for MITgcm coriolis
		radiusA=6.378E6; % equatorial radius (m)
		radiusB=6.357E6; % polar radius (m)
		radius0=1/sqrt((cos(phi0)/radiusA)^2 + (sin(phi0)/radiusB)^2); % ref. radius (m)
		P1.beta= 2*omega*cos(phi0)/radius0; % (partial f)/(partial y) (m*s)^-1
	end

	% free surface
	P1.rigidLid='.FALSE.'; % rigid lid off
	P1.implicitFreeSurface='.TRUE.'; % implicit free surface on
	P1.exactConserv='.TRUE.'; % exact total volume conservation (recompute divergence after pressure solver) on
	%P1.nonlinFreeSurf=0; % use the linear free surface
	P1.useRealFreshWaterFlux = '.FALSE.'; % Conserve volume with freshwater flux (changes free surface/sea level)

	% full and partial cells
	P1.hFacMin = 0.2; % minimum fractional height of partial gridcell
	%P1.hFacSup = 2.0; % maximum fractional height of surface cell
	%P1.hFacInf = 0.2; % minimum fractional height of surface cell

	% Momentum Equations
	%P1.vectorInvariantMomentum = '.TRUE.'; % use vector-invariant form of momentum equations
	P1.implicitViscosity = '.TRUE.'; % compute vertical diffusive fluxes of momentum implicitly
	P1.viscAr=1.E-4;    % Vertical Eddy Viscosity m^2/s
	viscAhscheme='modifiedLeith'; % choose between constant horizontal viscosity and gridscale/flow aware viscosity
	switch viscAhscheme
		case 'constant'
			P1.viscAh=10.0; % Horizontal Eddy Viscosity m^2/s
		case 'modifiedLeith'
			P1.viscAhGrid=.01;  % non-dimensional Laplacian grid viscosity parameter (background, constant)
			P1.viscAhGridMax=1; % maximum non-dimensional gridscale viscosity (CFL constraint)
			P1.viscC2leith=2; % Leith harmonic visc. factor on grad(vort), non-dim. (grid/flow aware viscosity)
			P1.viscC2leithD=2; % Leith harmonic visc. damping factor on grad(div), non-dim. (grid/flow aware viscosity)
		otherwise
			error('choose valid viscAhscheme')
	end
	P1.bottomDragQuadratic=2.5E-3; % no-slip bottom, quadratic bottom drag
	P1.bottomVisc_pCell='.TRUE.';  % account for partial-cell in bottom viscosity
	P1.selectImplicitDrag=2;       % fully implicit top/bottom drag
	P1.selectBotDragQuadr=2;       % quadr. bottom drag discretization: use local velocity norm @u,v location, w wet-point averaging of other velocity component

	% Tracer equations
	P1.implicitDiffusion = '.TRUE.'; % compute vertical diffusive fluxes of tracers implicitly
	P1.diffKhT= 3.; % horizontal Laplacian diffusivity of heat m^2/s
	P1.diffKrT=5.E-5; % vertical Laplacian diffusivity of heat m^2/s
	P1.diffKhS= 3.; % horizontal Laplacian diffusivity of salt m^2/s
	P1.diffKrS=5.E-5; % vertical Laplacian diffusivity of salt m^2/s
	P1.tempAdvScheme=77; % 2nd order flux limiters
	P1.saltAdvScheme=77; % 2nd order flux limiters
	P1.staggerTimestep = '.TRUE.'; % use staggered time stepping (thermodynamic vs. flow variables)

	% input/output
	P1.useSingleCpuIO='.TRUE.'; % only root MPI process does I/O
	P1.globalFiles='.TRUE.'; % write global files, not per tile
	P1.readBinaryPrec=64; % precision used for reading binary files (32 or 64)
	P1.debuglevel = 1; % level of printing of MITgcm activity messages/statistics (1-5)
	% }}}
	% PARM02: Elliptic solver parameters {{{
	% structure information
	P2.header='PARM02';
	P2.description='Elliptic solver parameters';

	% pressure solver
	P2.cg2dMaxIters=1000; % upper limit on 2D conjugate gradient solver iterations
	P2.cg2dTargetResidual=1.E-11; % 2D conjugate gradient target residual (non-dim. due to RHS normalization )
	% }}}
	% PARM03: Time stepping parameters {{{
	% structure information
	P3.header='PARM03';
	P3.description='Time stepping parameters';

	% Run Start and Duration
	P3.nIter0=[];     % starting timestep iteration number
	P3.deltaT=[];     % model time step (s)
	P3.nTimeSteps=[]; % number of model clock timesteps to execute

	% Ocean Convection
	P3.cAdjFreq = -1; % <0 sets frequency of convective adj. scheme to deltaT (s)

	% Restart/Pickup Files
	P3.pChkptFreq=[]; % permanent restart/pickup checkpoint file write interval (s)
	P3.ChkptFreq=[]; % rolling restart/pickup checkpoint file write interval (s)

	% Frequency/Amount of Output
	P3.monitorFreq=[];  % interval to write monitor output (s)
	P3.monitorSelect=1; % 1: dynamic variables only, 2: add vorticity variables, 3: add surface variables
	P3.dumpInitAndLast='.FALSE.'; % write out initial and last iteration model state (OFF)
	% }}}
	% PARM04: Gridding parameters {{{
	% structure information
	P4.header='PARM04';
	P4.description='Gridding parameters';

	% coordinate system
	P4.usingCartesianGrid='.TRUE.'; % use Cartesian coordinates
	%P4.delr=[num2str(Nz) '*' num2str(dz)]; % vertical grid spacing 1D array (m)
	P4.delRFile=['''' mit.fname.delrfile '''']; % vertical grid spacing 1D array (m)
	P4.dxSpacing=dx; % x-axis uniform grid spacing (m)
	P4.dySpacing=dy; % y-axis uniform grid spacing (m)
	P4.xgOrigin=x0; % west edge x-axis origin (m)
	P4.ygOrigin=y0; % south edge y-axis origin (m)
	% }}}
	% PARM05: Input datasets {{{
	% structure information
	P5.header='PARM05';
	P5.description='Input datasets';
	% input files
	P5.bathyfile    =['''' mit.fname.bathyfile '''']; % filename for 2D ocean bathymetry (m)
	% }}}
	mit.inputdata.PARM={P1,P2,P3,P4,P5};
	% }}}
	% input/data.pkg {{{
	PKG=struct; % initialize PACKAGES structure
	% PACKAGES: run-time flags for packages to use {{{
	PKG.header='PACKAGES';

	PKG.useDiagnostics='.TRUE.';
	PKG.useOBCS='.TRUE.';
	PKG.useShelfIce='.TRUE.';
	PKG.useCAL='.TRUE.';
	PKG.useEXF='.TRUE.';
	% }}}
	mit.inputdata.PKG={PKG};
	% }}}
	% Run-time package options
	% input/data.obcs {{{
	OBCS_P1=struct;OBCS_P2=struct;OBCS_P3=struct; % initialize OBCS_PARM structures
	% OBCS_PARM01: Open boundaries {{{
	OBCS_P1.header='OBCS_PARM01';
	OBCS_P1.description='Open boundaries';

	% Southern Boundary (bottom y boundary, not geographic south)
	OB_J=1; % j index wrt Ny where the southern OBCS is set
	OBCS_P1.OB_Jsouth=[num2str(Nx) '*' num2str(OB_J)]; % Nx-vector of J-indices (w.r.t. Ny) of Southern OB at each I-position (w.r.t. Nx)
	OBCS_P1.OBSvFile=[]; % Nx by Nz matrix of v velocity at Southern OB
	OBCS_P1.OBStFile=[]; % Nx by Nz matrix of pot. temp. at Southern OB
	OBCS_P1.OBSsFile=[]; % Nx by Nz matrix of salin. at Southern OB

	% Western Boundary (lefthand x boundary, not geographic west)
	OB_I=1; % i index wrt Nx where the western OBCS is set
	OBCS_P1.OB_Iwest=[num2str(Ny) '*' num2str(OB_I)];  % Ny-vector of I-indices (w.r.t. Nx) of Western OB at each J-position (w.r.t. Ny)
	OBCS_P1.OBWuFile=[]; % Nx by Nz matrix of u velocity at Western OB
	OBCS_P1.OBWtFile=[]; % Nx by Nz matrix of pot. temp. at Western OB
	OBCS_P1.OBWsFile=[]; % Nx by Nz matrix of salin. at Western OB

	OBCS_P1.useOBCSprescribe='.TRUE.'; % prescribe OB conditions
	OBCS_P1.useOBCSbalance='.FALSE.'; % use OB balance
	%OBCS_P1.OBCS_balanceFacN=1.0; % balance factor for N OB
	%OBCS_P1.OBCS_balanceFacW=1.0; % balance factor for W OB

	% Sponge layer
	useOBCSsponge=1; % on or off
	if useOBCSsponge==1
		OBCS_P1.useOBCSsponge='.TRUE.';    % use sponge layers
		OBCS_P1.useLinearSponge='.TRUE.';  % use linear sponge layer
	end
	% }}}
	% OBCS_PARM02: Orlanski parameters {{{
	OBCS_P2.header='OBCS_PARM02';
	OBCS_P2.description='Orlanski parameters';
	% }}}
	% OBCS_PARM03: Sponge layer parameters {{{
	OBCS_P3.header='OBCS_PARM03';
	OBCS_P3.description='Sponge layer parameters';

	if useOBCSsponge==1
		OBCS_P3.spongeThickness=5; % sponge layer thickness (in grid points)
		OBCS_P3.Vrelaxobcsbound=[]; % relaxation time scale at the outermost sponge layer point of a zonal OB (s)
		OBCS_P3.Urelaxobcsbound=[]; % relaxation time scale at the outermost sponge layer point of a meridional OB (s) 
	end
	% }}}
	mit.inputdata.OBCS = {OBCS_P1,OBCS_P2, OBCS_P3};
	% }}}
	% input/data.shelfice {{{
	SHELFICE_PARM01=struct; % initialize SHELFICE_PARM structure
	% SHELFICE_PARM01: Parameters for SHELFICE package {{{
	SHELFICE_P1.header='SHELFICE_PARM01';

	% general options
	SHELFICE_P1.SHELFICEwriteState='.TRUE.'; % write ice shelf state to file
	SHELFICE_P1.SHELFICEconserve='.TRUE.'; % use conservative form of temperature boundary conditions
	%SHELFICE_P1.SHELFICEMassStepping = '.TRUE.'; % recalculate ice shelf mass at every time step

	% input files
	SHELFICE_P1.SHELFICEtopoFile=['''' mit.fname.draftfile '''']; % filename for under-ice topography of ice shelves

	% boundary layer options
	SHELFICE_P1.SHELFICEboundaryLayer='.TRUE.'; % use simple boundary layer mixing parameterization
	%SHELFICE_P1.SHI_withBL_realFWflux='.TRUE.'; % use real-FW flux from boundary layer
	%SHELFICE_P1.SHI_withBL_uStarTopDz='.TRUE.'; % compute uStar from uVel, vVel averaged over top Dz thickness

	% thermodynamic exchange options (SHELFICEsaltTransCoef set automatically through SHELFICEsaltToHeatRatio)
	SHELFICE_P1.SHELFICEuseGammaFrict='.TRUE.'; % use velocity dependent exchange coefficients (Holland and Jenkins 1999) 
	SHELFICE_P1.SHELFICEheatTransCoeff=1.2E-4; % transfer coefficient for temperature (m/s)
	SHELFICE_P1.SHELFICE_transition_gamma='.TRUE.'; % use Dan's transition from constant to vel dependent gamma
	SHELFICE_P1.SHELFICETransGammaThickness=200.; % water column thickness at which to transition

	% drag options
	SHELFICE_P1.SHELFICEDragQuadratic=.006; % quadratic drag coefficient at bottom ice shelf (non-dim.)
	SHELFICE_P1.shiCdrag=SHELFICE_P1.SHELFICEDragQuadratic; % set to be the same
	% }}}
	mit.inputdata.SHELFICE = {SHELFICE_P1};
	% }}}
	% input/data.cal {{{
	CAL=struct; % intialize CAL structure
	% CAL_NML: The calendar package namelist {{{
	CAL.header='CAL_NML';

	CAL.TheCalendar='''model'''; % choose 'model' calendar
	CAL.startdate_1=[]; % yyyymmdd of start date
	%CAL.startDate_2=[]; % hhMMss of start date
	% Benjy: I am not sure exactly what calendarDumps does.
	CAL.calendarDumps='.FALSE.'; % align the output with calendar months?
	% }}}
	mit.inputdata.CAL = {CAL};
	% }}}
	% input/data.exf {{{
	EXF1=struct; EXF2=struct; EXF3=struct; EXF4=struct; EXF_OBCS=struct; % intialize EXF structures
	% EXF_NML_01: External Forcings namelist 1 {{{
	EXF1.header='EXF_NML_01';
	EXF1.useExfCheckRange  = '.FALSE.'; % check range of input fields and stop if out of range
	EXF1.exf_iprec = 64; % precision of input fields (32-bit or 64-bit)
	% }}}
	% EXF_NML_02: External Forcings namelist 2 {{{
	EXF2.header='EXF_NML_02';
	% }}}
	% EXF_NML_03: External Forcings namelist 3 {{{
	EXF3.header='EXF_NML_03';
	% }}}
	% EXF_NML_04: External Forcings namelist 4 {{{
	EXF4.header='EXF_NML_04';
	% }}}
	% EXF_OBCS: External Forcings OBCS options {{{
	EXF_OBCS.header='EXF_NML_OBCS';
	EXF_OBCS.useOBCSYearlyFields = '.TRUE.'; % append current year postfix of form _YYYY on filename
	obcsWstartdate1 = []; % W boundary - YYYYMMDD; start year (YYYY), month (MM), day (DD) of field to determine record number
	obcsWperiod     = -1;     % interval between two records: the special value -1 means non-repeating (calendar) monthly records
	obcsSstartdate1 = []; % S boundary - YYYYMMDD; start year (YYYY), month (MM), day (DD) of field to determine record number
	obcsSperiod     = -1;     % interval between two records: the special value -1 means non-repeating (calendar) monthly records
	% }}}
	mit.inputdata.EXF = {EXF1,EXF2,EXF3,EXF4,EXF_OBCS};
	% }}}
	% input/data.diagnostics {{{
	DIAG_LIST=struct; DIAG_STATIS=struct; % intialize DIAG structures
	% DIAGNOSTICS_LIST: Diagnostic Package Choices {{{
	DIAG_LIST.header='DIAGNOSTICS_LIST';
	DIAG_LIST.description={'Diagnostic Package Choices', ...
		'--------------------', ...
		'  dumpAtLast (logical): always write output at the end of simulation (default=F)', ...
		'  diag_mnc   (logical): write to NetCDF files (default=useMNC)', ...
		'--for each output-stream:', ...
		'  fileName(n) : prefix of the output file name (max 80c long) for outp.stream n', ...
		'  frequency(n):< 0 : write snap-shot output every |frequency| seconds', ...
		'               > 0 : write time-average output every frequency seconds', ...
		'  timePhase(n)     : write at time = timePhase + multiple of |frequency|', ...
		'    averagingFreq  : frequency (in s) for periodic averaging interval', ...
		'    averagingPhase : phase     (in s) for periodic averaging interval', ...
		'    repeatCycle    : number of averaging intervals in 1 cycle', ...
		'  levels(:,n) : list of levels to write to file (Notes: declared as REAL)', ...
		'                when this entry is missing, select all common levels of this list', ...
		'  fields(:,n) : list of selected diagnostics fields (8.c) in outp.stream n', ...
		'                (see "available_diagnostics.log" file for the full list of diags)', ...
		'  missing_value(n) : missing value for real-type fields in output file "n"', ...
		'  fileFlags(n)     : specific code (8c string) for output file "n"', ...
		'--------------------'};
	% THESE ALL GET SET IN THE TIMESTEPPING STEP
	% }}}
	% DIAG_STATIS_PARMS: Diagnostic Per Level Statistics {{{
	DIAG_STATIS.header='DIAG_STATIS_PARMS';
	DIAG_STATIS.description={'Parameter for Diagnostics of per level statistics:', ...
		'--------------------', ...
		'  diagSt_mnc (logical): write stat-diags to NetCDF files (default=diag_mnc)', ...
		'  diagSt_regMaskFile : file containing the region-mask to read-in', ...
		'  nSetRegMskFile   : number of region-mask sets within the region-mask file', ...
		'  set_regMask(i)   : region-mask set-index that identifies the region "i"', ...
		'  val_regMask(i)   : region "i" identifier value in the region mask', ...
		'--for each output-stream:', ...
		'  stat_fName(n) : prefix of the output file name (max 80c long) for outp.stream n', ...
		'  stat_freq(n):< 0 : write snap-shot output every |stat_freq| seconds', ...
		'               > 0 : write time-average output every stat_freq seconds', ...
		'  stat_phase(n)    : write at time = stat_phase + multiple of |stat_freq|', ...
		'  stat_region(:,n) : list of "regions" (default: 1 region only=global)', ...
		'  stat_fields(:,n) : list of selected diagnostics fields (8.c) in outp.stream n', ...
		'                (see "available_diagnostics.log" file for the full list of diags)', ...
		'--------------------'};
	% }}}
	mit.inputdata.DIAG = {DIAG_LIST,DIAG_STATIS};
	% }}}
	savedata(org,mit); % save the model
end % }}}
if perform(org,'Bathymetry'), % {{{
	mit=load(org,'MeshInit');
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Bathymetry
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% The bathymetry is interpolated from Bedmachine Antarctica, but adjusted along the bottom and left boundaries of the domain to 
	% match the boundary condition forcing files from Naughten et al., 2023 exactly. This transition is accomplished by a piece-wise 
	% linear transition over a length Lb. Along the boundary, we define the bed as being the bottom of the first real valued cell 
	% when viewed from below. This definition does not address nan data where an ice shelf exists along the boundary, and maintaining 
	% constant ice thickness along the boundary is accomplished elsewhere. 
	% Due to the way OBCS is executed, the bathymetry must match the OBCS files at the INSIDE edge of the cell, not the center. This
	% means that we set the first TWO cells along the boundary to match the Naughten et al., 2023 bathymetry.
	% Bear Ridge, which extends as a seafloor prominence north from Bear Island/Peninsula, is treated in the Naughten 2023 model with
	% a solid wall extending along the eastern 300m depth contour from the sea floor to the surface. This wall is intended to represent
	% the tendency of large ice bergs to ground against Bear Ridge, diverting the current. Here we follow a similar implementation.
	% The wall extends from the seafloor to sea level, and is implemented along Bear Ridge between the eastern 300m depth contour and 
	% the top of the ridge, as defined by the change in ridge aspect. A major difference between the two implementations is that Naughten
	% 2023 have a continuous wall from Bear Peninsula, while in our model a gap with a depth of around 400m separates the wall from 
	% the grounded ice on Bear Peninsula, allowing water flow over the ridge at the 400m depth saddle. This gap is not allowed in the
	% the Naughten model despite the depth. The wall outline is also transitioned over Lb to fit the location of the wall at the 
	% boundary.

	% Naughten 2023 Bathymetry
	disp('Processing Naughten 2023 Bathymetry');
	fname = fullfile(mit.forcing.Ndir,mit.forcing.init_example); % one of the init files
	fieldname = mit.forcing.field_init{1}; % one of the init state variable fields (Tinit)
	A = getfield(load(fname,fieldname),fieldname); % load state variable
	A = pagetranspose(reshape(A,size(A,2),size(A,1),size(A,3))); % correct the indexing to column-major
	A(end:mit.mesh.Ny,:) = NaN; % extend to actual boundary
	nanbed = flip(~isnan(A),3); % binary mask of real-valued cells
	% define bed from bottom up (ignore the ice shelf)
	flipzg = flip(circshift(mit.mesh.zp_N(1:end-1),-1)); % these are the bottoms of the cells, when viewed from below
	[~,~,flipZG] = meshgrid(mit.mesh.xc,mit.mesh.yc,flipzg);
	% pull out the first cell center in each column that is given an nan value
	[~,i]=max(nanbed,[],3,'linear'); 
	B_prime = flipZG(i); % B_prime is the bed from Naughten et al. without using the partial cells at the base
	% find the Naughten wall contour around Bear Ridge
	x_bear = -16.80E5; % x location on Bear Ridge (m)
	y_bear = -6.33E5;  % y location on Bear Ridge (m)
	[cc] = contourc(mit.mesh.xc,mit.mesh.yc,B_prime,[0 0]); % 0m depth contours for Naughten
	i=1; % index of this contour
	while i<length(cc)
		Ni = cc(2,i); % number of vertices in contour
		xcc = cc(1,(i+1):(i+Ni)); % vertex locations in x (m)
		ycc = cc(2,(i+1):(i+Ni)); % vertex locations in y (m)
		if(inpolygon(x_bear,y_bear,xcc,ycc))
			break
		end
		i = i+Ni+1; % go to the next contour
	end
	polyB_prime = polyshape(xcc,ycc); % polyshape of 0m depth contour for Naughten
	wallmaskB_prime = reshape(isinterior(polyB_prime,mit.mesh.hXC(:),mit.mesh.hYC(:)),size(B_prime)); % the wall mask for Naughten

	% Bedmachine Bathymetry
	disp('Processing Bedmachine Bathymetry');
	bedmachinepath='/nobackup/bgetraer/ModelData/BedMachine/BedMachineAntarctica-v4.0.nc'; % path to dataset
	B = interpBedmachineAntarctica(mit.mesh.hXC,mit.mesh.hYC,'bed','nearest',bedmachinepath); % B is the bed from Bedmachine Antarctica

	% The bathymetry interpolation scheme
	B_prime(wallmaskB_prime)=min(B(wallmaskB_prime),0); % do not interpolate over Naughten's wall, just use Bedmachine
	% shifting the bed to ensure that the first two boundary cells will match the bathymetry
	B_prime(2,:) = B_prime(1,:); 
	B_prime(:,2) = B_prime(:,1); 
	Lb = 10E3; % distance over which to transition the beds (m)
	a1 = max(min(min(1/Lb.*(mit.mesh.hYC-mit.mesh.yc(2)),1/Lb.*(mit.mesh.hXC-mit.mesh.xc(2))),1),0); % weighting between the two beds (0--1)
	% The adjusted bathymetry
	B_dagger = min(a1.*B + (1-a1).*B_prime, 0); % Adjusted bathymetry, capped at sea level (m)

	% Our Bear Ridge wall
	disp('Processing Bear Ridge Wall');
	% find the Bedmachine contour around Bear Ridge
	[cc] = contourc(mit.mesh.xc,mit.mesh.yc,B,[-300 -300]); % 300m depth contours for Bedmachine
	i=1; % index of this contour
	while i<length(cc)
		Ni = cc(2,i); % number of vertices in contour
		xcc = cc(1,(i+1):(i+Ni)); % vertex locations in x (m)
		ycc = cc(2,(i+1):(i+Ni)); % vertex locations in y (m)
		if(inpolygon(x_bear,y_bear,xcc,ycc))
			break
		end
		i = i+Ni+1; % go to the next contour
	end
	polyB = polyshape(xcc,ycc); % polyshape of 300m depth contour for Bedmachine

	% find the ridge line along Bear Ridge (the ridge runs roughly SW to NE so we use the gradient direction to extract the ridge)
	Bflat = B; % Bflat is the bathymetry with everything flattened besides Bear Ridge
	Bflat(~inpolygon(mit.mesh.hXC,mit.mesh.hYC,xcc,ycc))=min(B(inpolygon(mit.mesh.hXC,mit.mesh.hYC,xcc,ycc)),[],[1,2]);
	Bfilt = imgaussfilt(Bflat,1); % Gaussian smoothing before gradient extraction
	[Gmag,Gdir] = imgradient(Bfilt); % extract gradient and gradient direction
	[cc]=contourc(mit.mesh.xc,mit.mesh.yc,Gdir,[0,0]); % find the zero contour of the gradient direction
	% extract the largest contour
	i=1; Ni=0;
	while i<length(cc)
		thisNi = cc(2,i);
		if thisNi > Ni
			Ni = thisNi;
			xcdir = cc(1,(i+1):(i+Ni));
			ycdir = cc(2,(i+1):(i+Ni));
		end
		i = i+thisNi+1;
	end
	polyGdir = polyshape(xcdir,ycdir); % polyshape of Gdir 0 contour
	polywall = subtract(polyB,polyGdir); % difference (where the wall goes in B)
	polywall = polywall.regions; % make sure we only take the largest region
	polywall = polywall(1); % the largest region
	wallmaskB = reshape(isinterior(polywall,mit.mesh.hXC(:),mit.mesh.hYC(:)),size(B)); % the wall mask in Bedmachine

	% interpolate the wall close to the boundary to fit Naughten's wall at the boundary (like the Bathymetry interpolation)
	bdryVertices_N = polyB_prime.Vertices([1,end],:); % the wall polygon vertices at the boundary for Naughten
	% fit upper wall to Naughtens wall at boundary 
	ind = polywall.Vertices(:,2)<bdryVertices_N(1,2) & polywall.Vertices(:,2)>bdryVertices_N(2,2);
	ind = 1:find(~ind,1)-1;
	a1 = max(min(1/Lb.*(polywall.Vertices(ind,1)-mit.mesh.xc(2)),1),0);
	polywall.Vertices(ind,2) = a1.*polywall.Vertices(ind,2) + (1-a1).*bdryVertices_N(1,2);
	% fit lower wall to Naughtens wall at boundary
	ind = polywall.Vertices(:,2)<bdryVertices_N(2,2); % fit wall to Naughten's lower wall at boundary 
	a1 = max(min(1/Lb.*(polywall.Vertices(ind,1)-mit.mesh.xc(2)),1),0);
	polywall.Vertices(ind,2) = a1.*polywall.Vertices(ind,2) + (1-a1).*bdryVertices_N(2,2);

	% enforce wall in bathymetry
	wallmask = reshape(isinterior(polywall,mit.mesh.hXC(:),mit.mesh.hYC(:)),size(B_prime)); % the Bear Ridge wall mask
	B_dagger(wallmask)=0; % set Bear Ridge wall
	B_dagger(:,end)=0;    % set E boundary wall
	B_dagger(end,:)=0;    % set N boundary wall

	% find the partial topography cells and make a mask of the open gridcells
	disp('Processing Partial Topography Cells');
	dbathy = reshape(mit.mesh.zp,1,1,length(mit.mesh.zp)) - B_dagger; % distance between cell edge and bathymetry (m)
	ddbathy = abs(diff(dbathy>0,[],3)); % find where dbathy flips from negative to positive (boolean)
	[~,linear_ind]=max(ddbathy,[],3); % find which cell the bathymetry falls in (linear index)
	[k_ind] = ind2sub(size(ddbathy),linear_ind); % find upper boundary of cell bathymetry falls in (or on) (z-dim index)
	hFacDim = mit.mesh.zp(k_ind)-B_dagger; % the dimensional thickness of the final cell (m)
	hFacMin = 0.20; % the non-dimensional minimum fractional thickness of a gridcell
	hFacC = hFacDim./mit.mesh.delzF(k_ind); % the non-dimensional fractional thickness of the final cell
	ind = hFacC<hFacMin; % find small hFac
	roundValues = [0, hFacMin]; % the hFac values that we will round to 
	hFacC(ind) = interp1(roundValues,roundValues,hFacC(ind),'nearest'); % round small hFac to the nearest of 0 or hFacMin
	hFacDim = hFacC.*mit.mesh.delzF(k_ind); % the dimensional thickness of the final cell (m)
	opencell_mask = (reshape(mit.mesh.zp,1,1,length(mit.mesh.zp)) -  mit.mesh.zp(k_ind) + hFacDim) > 0; % mask of open cells in MITgcm

	% Save Bathymetry data
	mit.geometry=struct();
	mit.geometry.bathy   = B_dagger; % the bathymetry (m)
	mit.geometry.polywall= polywall; % the polyshape defining the bear ridge wall
	mit.geometry.hFacMin = hFacMin; % the non-dimensional minimum fractional thickness of a gridcell
	mit.geometry.hFacDim = hFacDim; % the dimensional thickness of the final cell (m)
	mit.geometry.hFacC   = hFacC;   % the non-dimensional fractional thickness of the final cell
	mit.geometry.k_ind   = k_ind;   % index of the upper cell edge above hFac (z-dim index)
	mit.geometry.open_mask = opencell_mask(:,:,1:end-1); % mask of open cells in MITgcm
	mit.forcing.obs_mask_N = isnan(squeeze(A(1,:,:)))'; % Naughten nanmask for bottom boundary
	mit.forcing.obw_mask_N = isnan(squeeze(A(:,1,:)))'; % Naughten nanmask for left boundary

	% Write to file
	fname = ['input/' mit.fname.bathyfile];
	write_binfile(fname,permute(mit.geometry.bathy,[2,1])); % write to binary input file with size Nx Ny
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% find where the boundary cells are open and closed and save to mit.inputdata.OBCS
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% S boundary
	nopenOBS = sum(mit.geometry.open_mask(1,:,:),3); % number of open cells in column 
	pind = find(diff(nopenOBS~=0)==1); % index from closed to open cell
	nind = find(diff(nopenOBS~=0)==-1); % index from open to closed cell
	ref = [zeros(size(pind)) ones(size(nind))]; 
	[ind, si] = sort([pind nind]); % where the transitions happen
	div = diff([0,ind,length(nopenOBS)]); % length of each division of open cells and closed cells
	if min(pind)<min(nind)
		ref=mod(0:length(div)-1,2); % start with closed cell, ie 0101010...
	else
		ref=mod(1:length(div),2); % start with open cell, ie 1010101...
	end
	ob_jsouth = '';
	for i=1:length(div)
		ob_jsouth = [ob_jsouth sprintf('%i*%i,',div(i),ref(i))];
	end
	mit.inputdata.OBCS{1}.OB_Jsouth=ob_jsouth(1:end-1);
	% W boundary
	nopenOBW = sum(mit.geometry.open_mask(:,1,:),3)'; % number of open cells in column 
	pind = find(diff(nopenOBW~=0)==1); % index from closed to open cell
	nind = find(diff(nopenOBW~=0)==-1); % index from open to closed cell
	ref = [zeros(size(pind)) ones(size(nind))]; 
	[ind, si] = sort([pind nind]); % where the transitions happen
	div = diff([0,ind,length(nopenOBW)]); % length of each division of open cells and closed cells
	if min(pind)<min(nind)
		ref=mod(0:length(div)-1,2); % start with closed cell, ie 0101010...
	else
		ref=mod(1:length(div),2); % start with open cell, ie 1010101...
	end
	ob_iwest = '';
   for i=1:length(div)
      ob_iwest = [ob_iwest sprintf('%i*%i,',div(i),ref(i))];
   end
   mit.inputdata.OBCS{1}.OB_Iwest=ob_iwest(1:end-1);

	% Save mit structure
	savedata(org,mit);
end % }}}
if perform(org,'BoundaryForcings'), % {{{
	mit=load(org,'Bathymetry');

	processBF = 0; % flag for actually reprocessing the boundary condition files
	if processBF
		D={}; % Initialize cell array for climate forcing data structures
		for i = 1:length(mit.forcing.exp)
			% Load all members of the experiment
			for j = 1:length(mit.forcing.member)
				fname = [mit.forcing.bdry_prefix mit.forcing.exp{i} mit.forcing.member{j} mit.forcing.suffix]; % file to be loaded
				fprintf('Loading '); ls(fullfile(mit.forcing.Ndir,fname)); % print out the file being read
				D{j} = load(fullfile(mit.forcing.Ndir,fname)); % load data into array of structures
			end
			% Loop through all fields
			for k=1:length(mit.forcing.field_bdry)
				disp(['Processing ' mit.forcing.field_bdry{k}]);
				% collect field from all members
				A = []; % matrix to collect field. dim1 is z, dim2 is x or y, dim3 is time (month), dim4 is member (1-10)
				for j = 1:length(D)
					A(:,:,:,j) = getfield(D{j},mit.forcing.field_bdry{k});
				end
				% take mean of field across members
				B = nanmean(A,4);

				% enforce nanmask on boundary condition (avoid interpolation error in velocity)
				switch mit.forcing.field_bdry_out{k}(end-2:end)
					case 'obw' % left boundary
						% the domain extent is actually further into the grounded ice: pad this with nan to reach the correct size
						B(:,end+1:mit.mesh.Ny,:)=NaN; % pad the domain with nan 
						% mask nan values for every page of the matrix
						for m = 1:size(B,3)
							ind0 = sub2ind(size(B),1,1,m); % linear index of (1,1,m) for page m
							B(find(mit.forcing.obw_mask_N) + ind0 - 1) = NaN; % set the NaN values for page m
						end
						% partial cell indexing
						k_ind = mit.geometry.k_ind(:,1); % vertical index where we need to correct for partial cell
						zc_hFac = mit.mesh.zp(mit.geometry.k_ind(:,1)) - mit.geometry.hFacDim(:,1)./2; % center for partial cell (m)
						open_ind = mit.geometry.hFacDim(:,1)>0;
					case 'obs'
						% mask nan values for every page of the matrix
						for m = 1:size(B,3)
							ind0 = sub2ind(size(B),1,1,m); % linear index of (1,1,m) for page m
							B(find(mit.forcing.obs_mask_N) + ind0 - 1) = NaN; % set the NaN values for page m
						end
						% partial cell indexing
						k_ind = mit.geometry.k_ind(1,:); % vertical index where we need to correct for partial cell
						zc_hFac = mit.mesh.zp(mit.geometry.k_ind(1,:)) - mit.geometry.hFacDim(1,:)./2; % center for partial cell (m)
						open_ind = mit.geometry.hFacDim(1,:)>0;
					otherwise 
						error('parsing failure on boundary identifier');
				end

				% interpolate from Naughten z grid onto the refined z grid
				Bq = interp1(mit.mesh.zc_N,B,mit.mesh.zc,'nearest'); % assign nearest real value to all cell centers
				Bq(isnan(Bq)) = 0; % set NaN values to zero (resulting in NaN in solution)

				% assign the correct value for the partial topography cells at the bed
				for m = 1:size(B,3)
					for ii = 1:size(B,2)
						if isnan(Bq(k_ind(ii),ii,m)) & open_ind(ii)
							Bq(k_ind(ii),ii,m) =  Bq(k_ind(ii)-1,ii,m); % replace NaN values with the cell value above 
						end
					end
				end

				% divide B into one file per year and write to binary file for MITgcm
				years = double(D{j}.ystart:D{j}.yend);
				for n = 1:length(years)
					% the filename to write to (this will be read by OBCS)
					fname_out = [mit.forcing.field_bdry_out{k} '.' mit.forcing.exp{i} '_' num2str(years(n))];
					% print out the file being written
					if n==1 | n==length(years)
						disp(['    Saving ' fname_out]);
					elseif n==2
						disp('     ...');
					end
					% write the file
					ind = (1:12) + (n-1)*12; % the index for the months of this year
					write_binfile(fullfile(mit.forcing.Ddir,fname_out),permute(Bq(:,:,ind),[2,1,3]));
				end
			end
		end
	else % do not process the files, just set the start date and move on
		disp('processBF = FALSE, skipping file processing');
		D = load(fullfile(mit.forcing.Ndir,mit.forcing.bdry_example));
		years = double(D.ystart:D.yend);
	end
	% save mit
	mit.forcing.startdate = datetime(years(1),1,1);
	savedata(org,mit);
end % }}}
if perform(org,'InitialConditions'), % {{{
	mit=load(org,'BoundaryForcings');

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Initial Conditions
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% MITgcm needs U,V,T,S, and EtaN to be defined in order to run. Additionally, SHELFICE needs
	% the draft to be defined.
	%   U, V, and EtaN are not explicitly defined and default to zero.
	%   T and S are given as tRef and sRef, vertical profiles which are then applied over the 
	% entire domain. 
	%   draft is defined by interpolation from ISSM, using the ice and ocean masks to ensure
	% that grounded ice has a draft equal to the MITgcm bathymetry, and open ocean has a draft 
	% of zero.	

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Define T and S profiles
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Dinit = []; % inital data for T and S 
	for i = 1:length(mit.forcing.exp)
		% Load all members of the experiment
		D={}; % Initialize cell array for climate forcing data structures
		for j = 1:length(mit.forcing.member)
			fname = [mit.forcing.init_prefix mit.forcing.exp{i} mit.forcing.member{j} mit.forcing.suffix]; % file to be loaded
			fprintf('Loading '); ls(fullfile(mit.forcing.Ndir,fname)); % print out the file being read
			D{j} = load(fullfile(mit.forcing.Ndir,fname)); % load data into array of structures
		end
		% Loop through all fields
		for k=1:length(mit.forcing.field_init)
			disp(['Processing ' mit.forcing.field_init{k}]);
			% collect field from all members
			A = []; % matrix to collect field. dim1 is x, dim2 is y, dim3 is z, dim4 is member (1-10)
			for j = 1:length(D)
				A(:,:,:,j) = getfield(D{j},mit.forcing.field_init{k});
			end
			% take mean of field across members and horizontal dimensions
			B = squeeze(nanmean(A,[1,2,4]));
			z = mit.mesh.zc_N(~isnan(B));
			b = B(~isnan(B));
			zq = mit.mesh.zc;
			bq = interp1(z,b,zq,'linear','extrap');
			Dinit(k,:,i) = bq;
		end
	end
	% average to get Tref and Sref
	tRef = mean(Dinit(1,:,:),[3]);
	sRef = mean(Dinit(2,:,:),[3]);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Find the draft along the boundaries that must be enforced throughout the model
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% find the ice shelf thickness along the bottom boundary of the domain from Naughten
	dMask = diff(mit.forcing.obs_mask_N); % where the mask goes from real value to NaN (0 no change; 1 real value to NaN; -1 NaN to real value)
	[K,J] = ind2sub(size(dMask),find(dMask==-1)); % location where NaN values overhang real values (K is z dim index, J is x dim index)
	draftOBS = zeros(1,mit.mesh.Nx); % default to zero draft (m)
	draftOBS(J) = mit.mesh.zp_N(K+1); % the bottom of the ice draft of the OBS ice shelf in z (m)
	% find the ice shelf thickness along the left boundary of the domain from Naughten
	dMask = diff(mit.forcing.obw_mask_N); % where the mask goes from real value to NaN (0 no change; 1 real value to NaN; -1 NaN to real value)
	[K,I] = ind2sub(size(dMask),find(dMask==-1)); % location where NaN values overhang real values (K is z dim index, J is x dim index)
	draftOBW = zeros(1,mit.mesh.Ny); % default to zero draft (m)
	draftOBW(I) = mit.mesh.zp_N(K+1); % the bottom of the ice draft of the OBW ice shelf in z (m)
	% save to mit structure
	mit.geometry.draftOBS = draftOBS; % draft at bottom boundary (m)
	mit.geometry.draftOBW = draftOBW; % draft at left boundary (m)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Define initial ice shelf draft from ISSM
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('reading ice shelf draft from ISSM');
	fname = fullfile(proph_dir,'experiments/ISSM_initialization/Models/PROPHET_issm_init_TransientPrep.mat');
	md = loadmodel(fname);
	% interpolate initial ice draft and masks from ISSM
	draft = zeros(mit.mesh.Ny,mit.mesh.Nx);    % initialize matrix for ice draft depth (m)
	draft(:)=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.geometry.base,mit.mesh.hXC(:),mit.mesh.hYC(:),'default',0); % m
	mask_oce=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ocean_levelset,mit.mesh.hXC(:),mit.mesh.hYC(:),'default',-1); % -1 ocean, 1 grounded
	mask_ice=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ice_levelset,mit.mesh.hXC(:),mit.mesh.hYC(:),'default',1); % -1 ice, 1 no ice	
	draft(mask_oce>0)=mit.geometry.bathy(mask_oce>0); % set all grounded ice to have a draft equal to the bathymetry (m)	
	draft(mask_ice>0)=0; % set all open ocean to have zero draft (m)
	draft(1,:) = draftOBS; % set draft at bottom boundary (m)
	draft(:,1) = draftOBW; % set draft at left boundary (m)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% write input data to files
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	fname = ['input/' mit.fname.treffile]; 
	disp(['writing tRef to ' fname]);
	write_binfile(fname,tRef);
	fname = ['input/' mit.fname.sreffile]; 
	disp(['writing sRef to ' fname]);
	write_binfile(fname,sRef);
   fname = ['input/' mit.fname.draftfile];
	disp(['writing draft to ' fname]);
   write_binfile(fname,permute(draft,[2,1]));

	savedata(org,mit);
end % }}}
if perform(org,'CompileMITgcm'), % {{{
	mit=load(org,'InitialConditions');
	% Compile-time options {{{
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% code/SIZE.h
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	SZ=struct; % initialize SIZE.h structure
	SZ.reffile=fullfile(mitgcm_dir,'model/inc/SIZE.h'); % MITgcm example file
	% define SZ structure fields
	% note domain decomposition must follow: Nx= sNx*nSx*nPx, Ny = sNy*nSy*nPy
	SZ.sNx=mit.mesh.sNx;   % Number of X points in tile. 
	SZ.sNy=mit.mesh.sNy;   % Number of Y points in tile.
	SZ.OLx=3;   % Tile overlap extent in X.                
	SZ.OLy=3;   % Tile overlap extent in Y.                
	SZ.nSx=1;   % Number of tiles per process in X.        
	SZ.nSy=1;   % Number of tiles per process in Y.        
	SZ.nPx=mit.mesh.Nx/SZ.sNx/SZ.nSx; % Number of processes to use in X.         
	SZ.nPy=mit.mesh.Ny/SZ.sNy/SZ.nSy; % Number of processes to use in Y.         
	SZ.Nx =mit.mesh.Nx;  % Number of points in X for the full domain
	SZ.Ny =mit.mesh.Ny;  % Number of points in Y for the full domain
	SZ.Nr =mit.mesh.Nz;  % Number of points in vertical direction.
	% write to ./code/SIZE.h
	write_sizefile(mit.fname.sizefile,SZ);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% code/packages.conf
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	PKGCONF=struct; % initialize Package Configuration structure
	% define Package Configuration structure fields
	PKGCONF.description='Configuration File for Package Inclusion';
	PKGCONF.pkg={'gfd','obcs','shelfice','cal','exf','diagnostics'};
	% write to ./code/packages.conf
	write_pkgconffile(mit.fname.pkgconffile,PKGCONF);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% code/DIAGNOSTICS_SIZE.h
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	DIAG_SZ=struct; % initialize DIAGNOSTICS_SIZE.h structure
	DIAG_SZ.reffile=fullfile(mitgcm_dir,'pkg/diagnostics/DIAGNOSTICS_SIZE.h'); % MITgcm example file
	DIAG_SZ.numDiags = 10*numel(mit.mesh.zc)+5; % maximum size of the storage array for active 2D/3D diagnostics
	% write to ./code/DIAGNOSTICS_SIZE.h
   write_diagsizefile(mit.fname.diagsizefile,DIAG_SZ);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% code/OBCS_OPTIONS.h
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	OPT=struct; % initialize OBCS_OPTIONS.h structure
   OPT.reffile=fullfile(mitgcm_dir,'pkg/obcs/OBCS_OPTIONS.h'); % MITgcm example file
	% OBCS options to allow
   OPT.define = {'ALLOW_OBCS_SOUTH','ALLOW_OBCS_WEST','ALLOW_OBCS_PRESCRIBE',...
		'ALLOW_OBCS_SPONGE'}; 
	% OBCS options to block
   OPT.undef = {'ALLOW_OBCS_NORTH','ALLOW_OBCS_EAST','ALLOW_ORLANSKI','ALLOW_OBCS_BALANCE'};
   % write to ./code/OBCS_OPTIONS.h
	fname = 'code/OBCS_OPTIONS.h';
   write_optionsfile(fname,OPT);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% code/SHELFICE_OPTIONS.h
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	OPT=struct; % initialize SHELFICE_OPTIONS.h structure
   OPT.reffile=fullfile(mitgcm_dir,'pkg/shelfice/SHELFICE_OPTIONS.h'); % MITgcm example file
	% options to allow
   OPT.define = {};
	% options to block
   OPT.undef = {'ALLOW_ISOMIP_TD','ALLOW_SHELFICE_REMESHING'};
   % write to ./code/SHELFICE_OPTIONS.h
	fname = 'code/SHELFICE_OPTIONS.h';
   write_optionsfile(fname,OPT);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % code/CPP_OPTIONS.h
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   OPT=struct; % initialize OPTIONS.h structure
   OPT.reffile=fullfile(mitgcm_dir,'model/inc/CPP_OPTIONS.h'); % MITgcm example file
   % options to allow
   OPT.define = {'ALLOW_SOLVE4_PS_AND_DRAG','SOLVE_DIAGONAL_LOWMEMORY'};
   % options to block
   OPT.undef = {'NONLIN_FRSURF'};
   % write to ./code/CPP_OPTIONS.h
	fname = 'code/CPP_OPTIONS.h';
   write_optionsfile(fname,OPT);


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% code files with non-trivial alterations
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	shelfice_dir = fullfile(initdir, 'shelfice'); % directory where modified shelfice code is kept
	fnames = {'SHELFICE.h','shelfice_readparms.F', ...
		'shelfice_thermodynamics.F','shelfice_init_fixed.F'}; % the files we need
	nerr = 0;
	for i=1:length(fnames)
		fprintf(' - checking for file code/% -30s     ',[fnames{i} '...']);
		if exist(fullfile(shelfice_dir, fnames{i}))
			fprintf('YES\n');
			% link the found file
			file_path = fullfile(shelfice_dir, fnames{i}); % file location
			dest_path = fullfile(initdir,'code',fnames{i}); % file destination
			command = ['cp ' file_path ' ' dest_path];
			system(command);
		else
			fprintf('NO\n');
			nerr = nerr+1;
		end
	end
	if nerr>0
		error([num2str(nerr) ' compile-time files are missing from code/ directory'])
	end

	mit.build=struct();
	mit.build.SZ=SZ;
	mit.build.PKGCONF=PKGCONF;

	savedata(org,mit);
	% }}}
	% Compile {{{
	prompt = 'Compile MITgcm now? (''y'' or ''Y'' to proceed, ''n'' or ''N'' to skip)\n';
	txt=0;
	while txt==0;
		txt = input(prompt,'s');
		switch txt
			case {'y','Y'}
				cont = 1;
			case {'n','N'}
				cont = 0;
			otherwise
				txt=0;
		end
	end
	
	if cont
		disp('Compile!');

		% locate files and scripts
		genmake2=fullfile(mitgcm_dir,'/tools/genmake2');
		optfile=fullfile(mitgcm_dir,'tools/build_options/linux_amd64_ifort+mpi_ice_nas'); % (pleaides)
		%optfile=fullfile(mitgcm_dir, 'tools/build_options/linux_amd64_gfortran'); % (totten)
		% clear the build directory
		!rm -r ./build
		mkdir('build');
		cd('./build');
		% make the MITgcm executable
		command=[genmake2 ' -mpi -mo ../code -optfile ' optfile ' -rd ' mitgcm_dir];
		system(command); % generate Makefile
		system('make CLEAN');  % prepare for new compilation
		system('make -j10 depend'); % create symbolic links from the local directory to the source file locations
		system('make -j10');        % compile code and create executable file mitgcmuv
	end % }}}	
	cd(initdir);
end % }}}

% These steps diverge between experiments
if perform(org,'RuntimeOptionsOcean') % {{{
	mit=load(org,'CompileMITgcm');

	% TIME STEPPING
	disp(' - Setting timestepping options');
	% coupling time step parameters
	mit.timestepping=struct();
	mit.timestepping.y2d = 360; % use 12 months of 30 days each ('model' calendar, see input/data.cal) (d/yr) 
	mit.timestepping.y2s = mit.timestepping.y2d*24*60*60; % y2s using 'model' calendar (s/yr)
	mit.timestepping.spinupduration = 3*mit.timestepping.y2s; % spinup duration: 3 Model years (2010,2011,2012) (s)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% input/data Time stepping parameters
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Run Start and Duration
	mit.inputdata.PARM{3}.nIter0=0;        % starting timestep iteration number
	mit.inputdata.PARM{3}.deltaT=100.;     % model time step (s)
	mit.inputdata.PARM{3}.nTimeSteps=(mit.timestepping.spinupduration/mit.inputdata.PARM{3}.deltaT); % number of model clock timesteps to execute
	% Restart/Pickup Files
	mit.inputdata.PARM{3}.pChkptFreq=0;								% permanent pickup checkpoint file write interval (s)
	mit.inputdata.PARM{3}.ChkptFreq=mit.timestepping.y2s/24; % temporary pickup checkpoint file write interval - twice per model month (s)
	% Frequency/Amount of Output
	mit.inputdata.PARM{3}.monitorFreq=mit.timestepping.y2s/24; % interval to write monitor output - twice per model month (s)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% input/data.obcs Sponge layer parameters
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	mit.inputdata.OBCS{3}.Vrelaxobcsbound=1*(24*60*60); % relaxation time scale at the outermost sponge layer point of a zonal OB (s)
	mit.inputdata.OBCS{3}.Urelaxobcsbound=1*(24*60*60); % relaxation time scale at the outermost sponge layer point of a meridional OB (s)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% input/data.cal Calendar parameters
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	cal_startdate=mit.forcing.startdate;
	%cal_startdate.Year=cal_startdate.Year+1;
	mit.inputdata.CAL{1}.startdate_1=string(cal_startdate,'yyyyMMdd'); % yyyyMMdd of start date
	mit.inputdata.CAL{1}.startDate_2=string(cal_startdate,'HHmmss');   % HHmmss of start date

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% input/data.exf External forcing parameters
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	mit.inputdata.EXF{end}.obcsWstartdate1 = string(mit.forcing.startdate,'yyyyMMdd'); % W boundary start year (YYYY), month (MM), day (DD) to determine record number
	mit.inputdata.EXF{end}.obcsWperiod     = -1.0;     % interval between two records: the special value -1 means non-repeating (calendar) monthly records
	mit.inputdata.EXF{end}.obcsSstartdate1 = string(mit.forcing.startdate,'yyyyMMdd'); % S boundary start year (YYYY), month (MM), day (DD) to determine record number
	mit.inputdata.EXF{end}.obcsSperiod     = -1.0;     % interval between two records: the special value -1 means non-repeating (calendar) monthly records

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% input/data.diagnostics Diagnostic output parameters
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Output Stream 1: surfDiag (snapshot every month)
	mit.inputdata.DIAG{1}.N(1).filename  = '''surfDiag''';
	mit.inputdata.DIAG{1}.N(1).frequency = -mit.timestepping.y2s/12; % (s)
	mit.inputdata.DIAG{1}.N(1).fields    = {'SHIfwFlx','ETAN    ','SHIuStar','SHIForcT','SHItrans'};

	% Output Stream 2: dynDiag (time-average every month)
	mit.inputdata.DIAG{1}.N(2).filename  = '''dynDiag''';
	mit.inputdata.DIAG{1}.N(2).frequency = mit.timestepping.y2s/12; % (s)
	mit.inputdata.DIAG{1}.N(2).fields    = {'UVEL    ','VVEL    ','WVEL    ','THETA   ','SALT    '};

	% Output Stream 3: SHICE_fwFluxtave (time average twice per month)
	mit.inputdata.DIAG{1}.N(3).filename  = '''SHICE_fwFluxtave''';
	mit.inputdata.DIAG{1}.N(3).frequency = mit.timestepping.y2s/24; % (s)
	mit.inputdata.DIAG{1}.N(3).fields    = {'SHIfwFlx'};

	% print settings
	parm={'spinupduration',mit.timestepping.spinupduration./mit.timestepping.y2s,' y';...
		'deltaT',mit.inputdata.PARM{3}.deltaT,' s';...
		'ChkptFreq',mit.inputdata.PARM{3}.ChkptFreq./24/60/60,' d';...
		'relaxobcsbound',mit.inputdata.OBCS{3}.Vrelaxobcsbound./24/60/60,' d';...
		'startdate',[],string(cal_startdate);...
		'surfDiagfreq',mit.inputdata.DIAG{1}.N(1).frequency./24/60/60,' d';...
		'dynDiagfreq',mit.inputdata.DIAG{1}.N(2).frequency./24/60/60,' d';...
		'SHICE_fwFluxtavefrq',mit.inputdata.DIAG{1}.N(3).frequency./24/60/60,' d'};
	formatstr='% 30s = %0.1f%s\n';
	for i=1:size(parm,1)
		fprintf(formatstr,parm{i,:});
	end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% write all of the input data files
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp(' - Writing runtime options to input/');
	write_datafile(mit.fname.eedatafile,       mit.inputdata.EEP,      'EXECUTION ENVIRONMENT PARAMETERS');
	write_datafile(mit.fname.datafile,         mit.inputdata.PARM,     'MODEL PARAMETERS');
	write_datafile(mit.fname.datapkgfile,      mit.inputdata.PKG,      'PACKAGES');
	write_datafile(mit.fname.datashelficefile, mit.inputdata.SHELFICE, 'SHELFICE RUNTIME PARAMETERS');
	write_datafile(mit.fname.datacalfile,      mit.inputdata.CAL,      'CALENDAR PARAMETERS');
	write_datafile(mit.fname.dataexffile,      mit.inputdata.EXF,      'EXTERNAL FORCINGS PARAMETERS');
	write_datafile(mit.fname.datadiagfile,     mit.inputdata.DIAG,     'DIAGNOSTICS RUNTIME PARAMETERS');

	% Make one data.obcs file for each experiment
	for i = 1:length(mit.forcing.exp)
		% Bottom boundary
		mit.inputdata.OBCS{1}.OBSuFile=['''' mit.fname.uvelOBSfile  mit.forcing.exp{i} '''']; % Nx by Nz matrix of u velocity at Southern OB
		mit.inputdata.OBCS{1}.OBSvFile=['''' mit.fname.vvelOBSfile  mit.forcing.exp{i} '''']; % Nx by Nz matrix of v velocity at Southern OB
		mit.inputdata.OBCS{1}.OBStFile=['''' mit.fname.thetaOBSfile mit.forcing.exp{i} '''']; % Nx by Nz matrix of pot. temp. at Southern OB
		mit.inputdata.OBCS{1}.OBSsFile=['''' mit.fname.saltOBSfile  mit.forcing.exp{i} '''']; % Nx by Nz matrix of salin. at Southern OB
		% Left boundary
		mit.inputdata.OBCS{1}.OBWuFile=['''' mit.fname.uvelOBWfile  mit.forcing.exp{i} '''']; % Nx by Nz matrix of u velocity at Western OB
		mit.inputdata.OBCS{1}.OBWvFile=['''' mit.fname.vvelOBWfile  mit.forcing.exp{i} '''']; % Nx by Nz matrix of v velocity at Western OB
		mit.inputdata.OBCS{1}.OBWtFile=['''' mit.fname.thetaOBWfile mit.forcing.exp{i} '''']; % Nx by Nz matrix of pot. temp. at Western OB
		mit.inputdata.OBCS{1}.OBWsFile=['''' mit.fname.saltOBWfile  mit.forcing.exp{i} '''']; % Nx by Nz matrix of salin. at Western OB

		fname = [mit.fname.dataobcsfile mit.forcing.exp{i}];
		write_datafile(fname, mit.inputdata.OBCS, 'OBCS RUNTIME PARAMETERS');
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Diverge run directories for each experiment
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	prompt = 'Reset runocean directories now? (''y'' or ''Y'' to proceed, ''n'' or ''N'' to skip)\n';
	txt=0;
	while txt==0;
		txt = input(prompt,'s');
		switch txt
			case {'y','Y'}
				cont = 1;
			case {'n','N'}
				cont = 0;
			otherwise
				txt=0;
		end
	end

	if cont
		disp(' - Diverging models for climate experiments');
		for i = 1:length(mit.forcing.exp)
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% Directory management
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% make this experiment subdirectory if needed
			subdir=fullfile(proph_dir,'experiments',mit.forcing.exp{i});
			if ~exist(subdir)
				mkdir(subdir);
			end
			% rename previous run directory and create new one
			dirname='runocean';
			rundir=fullfile(subdir,dirname);
			oldrundir=fullfile(subdir,[dirname '.old']);
			if exist(oldrundir)
				system(['\rm -rf ' oldrundir]);
			end
			if exist(rundir)
				system(['\mv ' rundir ' ' oldrundir]);
			end
			% make the run directory in subdir
			mkdir(rundir);

			disp(['    linking files to ' rundir])

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% Link to run directory
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% make links to input files
			S = dir(fullfile(initdir,'input/*'));
			for j=1:length(S)
				if ~S(j).isdir
					if contains(S(j).name,'data.obcs')
						if contains(S(j).name,mit.forcing.exp{i})
							% only link data.obcs for the experiment we are running
							file_path = fullfile(S(j).folder, S(j).name); % file location
							link_path = fullfile(rundir,'data.obcs'); % link location
							command = ['ln -s ' file_path ' ' link_path];
							system(command);
						end
					else
						file_path = fullfile(S(j).folder, S(j).name); % file location
						link_path = rundir; % link location
						command = ['ln -s ' file_path ' ' link_path];
						system(command);
					end
				end
			end
			% make links to forcing files
			S = dir(mit.forcing.Ddir);
			for j=1:length(S)
				if ~S(j).isdir & contains(S(j).name,mit.forcing.exp{i})
					file_path = fullfile(S(j).folder, S(j).name); % file location
					link_path = rundir; % link location
					command = ['ln -s ' file_path ' ' link_path];
					system(command);
					if contains(S(j).name,num2str(mit.forcing.startdate.Year))
						file_path = fullfile(S(j).folder, S(j).name); % file location
						splstr = strsplit(S(j).name,'_'); % split the filename base and the year
						yr = num2str(mit.forcing.startdate.Year-1); % link to previous year
						link_path = fullfile(rundir,[splstr{1} '_' yr]); % link location from yr0-1 to yr0
						command = ['ln -s ' file_path ' ' link_path];
						system(command);
					end
				end
			end
			% make link to mitgcmuv executable
			file_path = fullfile(initdir,'build/mitgcmuv'); % file location
			link_path = rundir; % link location
			command = ['ln -s ' file_path ' ' link_path];
			system(command);
		end
	else
		disp('Skipping reset of runocean directories!');
	end

	savedata(org,mit);
end % }}}
if perform(org,'RunOcean') % {{{
	mit=load(org,'RuntimeOptionsOcean');

	% set run parameters for PBS queue file
	rundir = fullfile(expdir,'runocean'); % which experiment directory to run
	%rundir = '/nobackupp18/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/test/run';
	grouplist = 's2950'; % account on Pleiades
	npMIT=mit.build.SZ.nPx*mit.build.SZ.nPy; % number of processors for MITgcm
	queuename = 'long'; % which queue to submit to (long or devel)
	walltime = duration(120,0,0); % walltime to request
	% write the .queue file
	fname = write_queuefile(rundir,grouplist,npMIT,'queuename',queuename,'walltime',walltime,'iscoupled',0); % returns the name of the .queue file
	%fname = write_queuefile(rundir,grouplist,1,'HelloWorld'); % returns the name of the .queue file
	fprintf('Submitting queue file:   ')
	command=['qsub ' fullfile(rundir,fname)];
	system(command);
end % }}}
if perform(org,'RuntimeOptionsCoupled') % {{{
	mit=load(org,'RuntimeOptionsOcean');

	% TIME STEPPING
	% During the coupled phase our ocean model does two runs per coupled step:
	%		1) "Relaxation run"
	%			- immediately after updating the ice shelf geometry
	%			- uses a much smaller deltaT
	%			- runs for duration defined by relaxT
	%		2) "Continuation run"
	%			- uses pickup from the relaxation run
	%			- uses a larger deltaT
	%			- runs until the end of the coupledTimeStep
	% The two runs are defined by different sets of data parameters that are set 
	%	during the run
	disp(' - Setting timestepping options');
	mit.timestepping.ispickup         = 0;                                                         % are we running a pickup or an initial coupling [0,1]
	mit.timestepping.coupled_basetime = mit.timestepping.spinupduration;                           % the model time that we are starting coupling from (s)
	mit.timestepping.coupled_endtime  = 90*mit.timestepping.y2s;                                   % the end model time that we are aiming for (s)
	mit.timestepping.startTime        = mit.timestepping.coupled_basetime;                         % modeltime of the start of the simulation (s)
	mit.timestepping.deltaT_coupled   = coupling_dT;                                               % coupling time step (s)
	mit.timestepping.coupledT         = mit.timestepping.coupled_endtime;                          % duration of coupling simulation (s)
	mit.timestepping.nsteps           = mit.timestepping.coupledT/mit.timestepping.deltaT_coupled; % number of coupled time steps to take
	mit.timestepping.relaxT           = 1*60*60;                                                   % duration of relaxation time after coupling (s)
	mit.timestepping.deltaT_relax     = 10;                                                        % small deltaT to use for relaxation run (s)
	mit.timestepping.contT            = mit.timestepping.deltaT_coupled-mit.timestepping.relaxT;  % duration of continuation run after relaxation (s)
	mit.timestepping.deltaT_cont      = 100;                                                       % large deltaT to use for continuation run (s)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% ./data
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% These fields get set at each run within the coupling function. here are just placeholders.
	% Run Start and Duration 
	mit.inputdata.PARM{3}.deltaT      = 0; % mitgcm deltaT (s) 
	mit.inputdata.PARM{3}.startTime   = 0; % run start time for this integration (s)
   mit.inputdata.PARM{3}.nIter0      = 0; % starting timestep iteration number
   mit.inputdata.PARM{3}.nTimeSteps  = 0; % number of timesteps to execute
   % Restart/Pickup Files
   mit.inputdata.PARM{3}.pChkptFreq  = 0; % permanent pickup checkpoint file write interval (s)
   mit.inputdata.PARM{3}.ChkptFreq   = 0; % temporary pickup checkpoint file write interval (s)
   % Frequency/Amount of Output
   mit.inputdata.PARM{3}.monitorFreq = 0; % interval to write monitor output - every coupled time step (s)
	% Initialization files
	mit.fname.uvelfile  = 'uvel.bin';
	mit.fname.vvelfile  = 'vvel.bin';
	mit.fname.thetafile = 'theta.bin';
	mit.fname.saltfile  = 'salt.bin';
	mit.fname.etanfile  = 'etan.bin';
	mit.inputdata.PARM{5}.uVelInitFile    = ['''' mit.fname.uvelfile ''''];
	mit.inputdata.PARM{5}.vVelInitFile    = ['''' mit.fname.vvelfile ''''];
	mit.inputdata.PARM{5}.hydrogThetaFile = ['''' mit.fname.thetafile ''''];
	mit.inputdata.PARM{5}.hydrogSaltFile  = ['''' mit.fname.saltfile ''''];
	mit.inputdata.PARM{5}.pSurfInitFile   = ['''' mit.fname.etanfile ''''];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % ./data.diagnostics
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	mit.inputdata.DIAG{1}.N(1).frequency=0;
	mit.inputdata.DIAG{1}.N(2).frequency=0;
	mit.inputdata.DIAG{1}.N(3).frequency=0;

	% print settings
	disp('mit.timestepping:');
	disp(mit.timestepping);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Directory management
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	builddir = fullfile(initdir,'build'); % initalization directory where model was compiled
	inputdir=fullfile(initdir,'input'); % initialization directory for runtime input options
	runoceandir = fullfile(expdir,'runocean'); % run directory for the ocean model spinup
	runcoupleddir = fullfile(expdir,sprintf('runcoupled_dt%03i_ct%07i',mit_dT,coupling_dT)); % run directory
	% rename previous run directory and create new one
	oldruncoupleddir=[runcoupleddir '.old'];
	if exist(oldruncoupleddir), rmdir(oldruncoupleddir,'s'); end
	if exist(runcoupleddir), movefile(runcoupleddir,oldruncoupleddir); end
	% make the run directory in subdir
	mkdir(runcoupleddir);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Build run directory
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp(['  - Initializing coupledrun directory in ' runcoupleddir])
	% copy files from runoceandir: get the files matching the right start time, 
	% and copy to the runcoupleddir, using the modeltime as the new reference suffix
	modeltime_pickup = mit.timestepping.startTime;                                                       % the modeltime we want to start from
	pickup_fname     = searchMITgcmFile(runoceandir,'pickup','timeInterval',modeltime_pickup);           % find the matching pickup file
	melt_fname       = searchMITgcmFile(runoceandir,'SHICE_fwFluxtave','timeInterval',modeltime_pickup); % find the matching melt file
	% rename files to modeltime
	source_fnames={[pickup_fname '.meta'],[pickup_fname '.data'],[melt_fname '.meta'],[melt_fname '.data'],...
		'hFacC.meta','hFacC.data'};
	suffix=sprintf('save.%010i',modeltime_pickup);
	disp(['   copying ' num2str(numel(source_fnames)) ' files from']);
	disp(['       ' runoceandir ' to']);
	disp(['       ' runcoupleddir]);
	for i=1:numel(source_fnames)
		splstr=strsplit(source_fnames{i},'.');
		destination_fname=[splstr{1} '.'  suffix '.' splstr{end}];
		disp(sprintf('         - %-17s  ->  %s', source_fnames{i},destination_fname));
		copyfile(fullfile(runoceandir,source_fnames{i}),fullfile(runcoupleddir,destination_fname));
	end

	% COPY the input files
	dataobcsfile=['data.obcs' experiment.name];
	filelist={'eedata','data','data.cal','data.diagnostics','data.exf',dataobcsfile,'data.pkg','data.shelfice',...
		'bathy.bin','delr.bin','sref.bin','tref.bin'};
	cp_filelist(inputdir,filelist,runcoupleddir); % copy all files from inputdir to runcoupleddir
	% rename data.obcs
	source_fnames={'draft.bin',dataobcsfile};
	destination_fname={['draft.' suffix '.bin'], 'data.obcs'};
   disp(['   copying ' num2str(numel(source_fnames)) ' files from']);
   disp(['       ' inputdir ' to']);
   disp(['       ' runcoupleddir]);
   for i=1:numel(source_fnames)
      disp(sprintf('         - %-17s  ->  %s', source_fnames{i},destination_fname{i}));
      copyfile(fullfile(inputdir,source_fnames{i}),fullfile(runcoupleddir,destination_fname{i}));
   end

	% link mitgcmuv executable
	filelist={'mitgcmuv'};
	ln_filelist(builddir,filelist,runcoupleddir); % link file from builddir to runcoupleddir

	% link the boundary forcing files
	filelist={dir(fullfile(mit.forcing.Ddir,['*' experiment.name '*'])).name};
	ln_filelist(mit.forcing.Ddir,filelist,runcoupleddir);

	save(fullfile(runcoupleddir,'RuntimeOptionsCoupled'),'mit');
end % }}}
if perform(org,'RunCoupled') % {{{
	rundir = fullfile(expdir,sprintf('runcoupled_dt%03i_ct%07i',mit_dT,coupling_dT)); % run directory
	mitfile=fullfile(rundir,'RuntimeOptionsCoupled.mat');
	%rundir = fullfile(expdir,'runcoupled'); % which experimen directory to run
	mit=loadmodel(mitfile);
	
	% set run parameters for PBS queue file
	mccdir = fullfile(proph_dir,'runcouple/mccfiles/');
	mdfile = fullfile(proph_dir,'experiments/ISSM_initialization/Models/PROPHET_issm_init_TransientPrep.mat');
	grouplist = 's2950'; % account on Pleiades
	npMIT=mit.build.SZ.nPx*mit.build.SZ.nPy; % number of processors for MITgcm
	queuename = 'long'; % which queue to submit to (long or devel)
	%queuename = 'devel'; % which queue to submit to (long or devel)
	walltime = duration(4*24+11,0,0); % walltime to request
	%% write the .queue file
	fname = write_queuefile(rundir,grouplist,npMIT,...
		'queuename',queuename,'walltime',walltime,'iscoupled',1,'mccdir',mccdir,'mccargin',{mdfile,mitfile}); % returns the name of the .queue file
	fprintf('Submitting queue file:   ')
	command=['qsub ' fullfile(rundir,fname)];
	system(command);	
	% if interactive!!
	%cd(rundir);
	%%runalone(mdfile,mitfile);
	%addpath(fullfile(proph_dir,'runcouple'));
	%%runcouple(mdfile,mitfile);
	%runcouple_beta(mdfile,mitfile);
end % }}}
if perform(org,'RunPickup') % {{{
	rundir = fullfile(expdir,sprintf('runcoupled_dt%03i_ct%07i',mit_dT,coupling_dT)); % run directory
	mit=loadmodel(fullfile(rundir,'RuntimeOptionsCoupled'));

	modeltime_pickup = 839808000; % (s)

	% TIME STEPPING
   % During the coupled phase our ocean model does two runs per coupled step:
   %     1) "Relaxation run"
   %        - immediately after updating the ice shelf geometry
   %        - uses a much smaller deltaT
   %        - runs for duration defined by relaxT
   %     2) "Continuation run"
   %        - uses pickup from the relaxation run
   %        - uses a larger deltaT
   %        - runs until the end of the coupledTimeStep
   % The two runs are defined by different sets of data parameters that are set
   %  during the run
	disp(' - Setting timestepping options');
   mit.timestepping.ispickup         = 1;                                                         % are we running a pickup or an initial coupling [0,1]
   mit.timestepping.startTime        = modeltime_pickup;                                          % modeltime of the start of the simulation (s)
   mit.timestepping.nsteps           = mit.timestepping.coupledT/mit.timestepping.deltaT_coupled; % number of coupled time steps to take
   mit.timestepping.relaxT           = 12*60*60;                                                   % duration of relaxation time after coupling (s)
   mit.timestepping.deltaT_relax     = 5;                                                        % small deltaT to use for relaxation run (s)
   mit.timestepping.contT            = mit.timestepping.deltaT_coupled-mit.timestepping.relaxT;  % duration of continuation run after relaxation (s)
   mit.timestepping.deltaT_cont      = 100;                                                       % large deltaT to use for continuation run (s)

	% print settings
   disp('mit.timestepping:');
   disp(mit.timestepping);

	% save mit structure
	mitfile = fullfile(rundir,'RunPickup.mat');
	save(mitfile,'mit');

	% set run parameters for PBS queue file
	mccdir = fullfile(proph_dir,'runcouple/mccfiles/');
	mdfile = fullfile(proph_dir,'experiments/ISSM_initialization/Models/PROPHET_issm_init_TransientPrep.mat');
	grouplist = 's2950'; % account on Pleiades
	npMIT=mit.build.SZ.nPx*mit.build.SZ.nPy; % number of processors for MITgcm
	queuename = 'long'; % which queue to submit to (long or devel)
	walltime = duration(5*24,0,0); % walltime to request
	%queuename = 'devel'; % which queue to submit to (long or devel)
	%walltime = duration(1,0,0); % walltime to request
	% write the .queue file
	fname = write_queuefile(rundir,grouplist,npMIT,...
		'queuename',queuename,'walltime',walltime,'iscoupled',1,'mccdir',mccdir,'mccargin',{mdfile,mitfile}); % returns the name of the .queue file
	fprintf('Submitting queue file:   ')
	command=['qsub ' fullfile(rundir,fname)];
	system(command);	
	% if interactive!!
	%cd(rundir);
	%addpath(fullfile(proph_dir,'runcouple'));
	%runcouple_beta(mdfile,mitfile);
end % }}}

% Move back to root directory
disp(['Moving to root directory: ', proph_dir]);
cd(proph_dir);

return
% local functions 
function [fname] = searchMITgcmFile(parentdir,prefix,fieldname,value); % {{{
	% searchMITgcmFile finds all files in parentdir that match prefix*.meta and returns the 
	% filename (without extension) of the file which contains fieldname = value
	% Example: 
	%    searchMITgcmFile(rundir,'pickup','timeStepNumber',niter0);
	%    searchMITgcmFile(rundir,'SHICE_fwFluxtave','timeInterval',startTime);
	fnames=flip({dir(fullfile(parentdir, [prefix '*.meta'])).name}); % match prefix to filenames, search in reverse order 
	i=1;
	while i<=numel(fnames)
		fid=fopen(fullfile(parentdir,fnames{i}));
		tline=fgetl(fid); % read the next line
		while ischar(tline)
			if contains(tline, fieldname)
				break;
			end
			tline = fgetl(fid); % read the next line
		end
		fclose(fid);
		thisvalue=str2num(extractBefore(extractAfter(tline,'['),']'));
		if thisvalue(end)==value
			break;
		else
			i=i+1;
		end
	end
	if i>numel(fnames)
		error('No pickup file is found for modeltime_pickup!');
	else
		fname=extractBefore(fnames{i},'.meta');
	end
end % }}}
function D=readmeta(parentdir,fname,varargin) % {{{
%READMETA looks for a file of the form fname*.meta in parentdir
% if multiple files are matched, it finds all of them.
% The output is a cell array of structures containing the metadata
% of each matched file
   S=dir(fullfile(parentdir,[fname '*.meta']));
   D=[];
   command=''; % initialize blank command

   if numel(varargin)>0
      fields=varargin;
      for i=1:numel(S)
         fid=fopen(S(i).name); % open file
         D(i).fname=S(i).name; % save filename
         tline=fgetl(fid); % read line
         while ischar(tline)
            command = [command tline]; % build commmand
            if strcmp(tline(end),';')
               thisfield=strip(extractBefore(tline,'='));
               if any(strcmp(thisfield,fields))
                  command=strip(command,' '); % remove whitespace
                  command=['D(' num2str(i) ').' command]; % save in the structure
                  disp(command); % print
                  eval(command); % evaluate command
               end
               command=''; % reset command
            end
            tline=fgetl(fid); % read next line
         end
         fclose(fid); % close file
      end
   else
		fields=[];
      for i=1:numel(S)
         fid=fopen(S(i).name); % open file
         D(i).fname=S(i).name; % save filename
         tline=fgetl(fid); % read line
         while ischar(tline)
            command = [command tline]; % build commmand
            if strcmp(tline(end),';')
               command=strip(command,' '); % remove whitespace
               command=['D(i).' command]; % save in the structure
               disp(command); % print
               eval(command); % evaluate command
               command=''; % reset command
            end
            tline=fgetl(fid); % read next line
         end
         fclose(fid); % close file
      end
   end
end % }}}
function ln_filelist(parentdir,filelist,targetdir) % {{{
% LN_FILELIST soft-links a list of files located in parentdir to targetdir
	if ~isdir(parentdir)
		error('parentdir must be a directory');
	elseif any(~isfile(fullfile(parentdir,filelist)))
		error('filelist contains files which do not exist in parentdir');
	elseif  ~isdir(targetdir)
		error('targetdir must be a directory');
	end
	% link the files
	disp(['   linking ' num2str(numel(filelist)) ' files from ']);
   disp(['       ' parentdir ' to']);
   disp(['       ' targetdir]);
   for i=1:numel(filelist)
      file_path=fullfile(parentdir,filelist{i}); % file location
      command = ['ln -s ' file_path ' ' targetdir];
		if numel(filelist)<20
			disp(['         - ' filelist{i}]);
		elseif i==1
			disp(['         ...']);
		end
			system(command);
   end
end	% }}}
function cp_filelist(parentdir,filelist,targetdir) % {{{
% CP_FILELIST copies a list of files located in parentdir to targetdir
	if ~isdir(parentdir)
		error('parentdir must be a directory');
	elseif any(~isfile(fullfile(parentdir,filelist)))
		error('filelist contains files which do not exist in parentdir');
	elseif  ~isdir(targetdir)
		error('targetdir must be a directory');
	end
	% copy the files
	disp(['   copying ' num2str(numel(filelist)) ' files from ']);
	disp(['       ' parentdir ' to']);
	disp(['       ' targetdir]);
   for i=1:numel(filelist)
      file_path=fullfile(parentdir,filelist{i}); % file location
      command = ['cp ' file_path ' ' targetdir];
      disp(['         - ' filelist{i}]);
      system(command);
   end
end	% }}}
function fname=write_queuefile(rundir,grouplist,ncpus,varargin) % {{{
%WRITE_QUEUEFILE generates a .queue file to launch an MITgcm or coupled MITgcmXISSM model
% run on Pleiades using PBS
% 
% EXAMPLES:
%    fname = write_queuefile(rundir,grouplist,ncpus); % configures using defaults (uncoupled, devel queue, etc)
%    fname = write_queuefile(rundir,grouplist,ncpus,varargin); % configures using specified options
%    fname = write_queuefile(rundir,grouplist,1,'HelloWorld'); % configures a minimal working example of a queue file for testing
%    fname = write_queuefile(rundir,grouplist,ncpus,'queuename','long','walltime',duration(1,0,0),'iscoupled',0);
%
% INPUT:
%    rundir     string  - the full path of the MITgcm directory to run in
%    grouplist  string  - the acct grouplist for Pleiades to charge to
%    ncpus      numeric - number of cpus to request 
%    varargin:
%        queuename    string   - the name of the queue ('low','normal','long','debug','devel')
%        walltime     duration - the walltime requested (HH:MM:SS)
%        iscoupled    [0,1]    - 0: only MITgcm, 1: MITgcm and ISSM
%        mccdir       dir      - directory path where mcc files for coupled run are compiled
% OUTPUT:
%    fname      string  - filename of the .queue file which is written to rundir
% SUBFUNCTIONS:
%    resourcestring=buildresourcestring(ncpus,nodemodel)    returns a PBS resource string based on requested number of cpus
%    INPUT:
%       ncpus       numeric  - input from WRITE_QUEUEFILE
%       nodemodel   string   - defines # cpus per node ('bro' is only nodemodel supported currently)
%    OUTPUT:
%       resourcestring  string - formatted for PBS defining number of nodes and how many cpus per node
	
	% parse input
	% create inputParser object
	p = inputParser;
	% add inputs to the scheme
	defaultQueue='devel';
	validQueue={'low','normal','long','debug','devel'}; % see https://www.nas.nasa.gov/hecc/support/kb/pbs-job-queue-structure_187.html
	defaultWalltime=[duration(0,30,0),duration(1,0,0),duration(1,0,0),duration(0,30,0),duration(0,20,0)]; % based on queue chosen
	maxWalltime=[duration(4,0,0),duration(8,0,0),duration(120,0,0),duration(2,0,0),duration(2,0,0)]; % based on queue chosen
	checkQueue=@(x) any(validatestring(x,validQueue));

	checkIsbinary=@(x) any(x==[0,1]);
	checkIsnatural=@(x) (x>0 & mod(x,1)==0);
	checkMccdir=@(x) (isempty(x) | isdir(x));
	checkIsHelloWorld=@(x) (isempty(x) | strcmp(x,'HelloWorld'));

	checkMccargin=@(x) (isempty(x) | (iscell(x) & numel(x)==2 & isfile(x{1}) & isfile(x{2})));

	addRequired(p,'rundir',@isdir);
	addRequired(p,'grouplist',@ischar);
	addRequired(p,'ncpus',checkIsnatural);
	addOptional(p,'HelloWorld',0,checkIsHelloWorld);
	addParameter(p,'queuename',defaultQueue,checkQueue);
	addParameter(p,'walltime',[],@isduration);
	addParameter(p,'iscoupled',0,checkIsbinary);
	addParameter(p,'mccdir',[],checkMccdir);
	addParameter(p,'mccargin',[],checkMccargin);

	% parse the inputs and save locally
	parse(p,rundir,grouplist,ncpus,varargin{:})
	rundir    =p.Results.rundir;
	grouplist =p.Results.grouplist;
	ncpus     =p.Results.ncpus;
	queuename =p.Results.queuename;
	walltime  =p.Results.walltime;
	iscoupled =p.Results.iscoupled;
	mccdir    =p.Results.mccdir;
	mccargin  =p.Results.mccargin;
	if strcmp(p.Results.HelloWorld,'HelloWorld')
		ishelloworld=1;
		rundir    =p.Results.rundir;
		grouplist =p.Results.grouplist;
		ncpus     =p.Results.ncpus;
		queuename =p.Results.queuename;
		walltime  =[];
		iscoupled =0;
		mccdir    =[];
	else
		ishelloworld=0;
	end
	clear p;

	% deal with walltime
	if isempty(walltime)
		walltime=defaultWalltime(strcmp(queuename,validQueue));
	end
	if walltime<=0 | walltime>maxWalltime(strcmp(queuename,validQueue))
		error(['Walltime ' char(walltime) ' is not valid for ' queuename ' queue!']);
	end

	% deal with mccdir
	if iscoupled & isempty(mccdir)
		error('Coupled scheme requires mccdir to be defined!');
	elseif ~iscoupled & isdir(mccdir)
		error('mccdir is defined for a non-coupled scheme!');
	end

	%set .queue filename
	pathparts=strsplit(rundir,'/'); % split directory name
	prefix=pathparts{end-1}; % experiment prefix
	if ishelloworld
      prefix='HelloWorld';
      fname=[prefix '.queue']; % filename for hellow world file
	elseif iscoupled
		fname=[prefix '_runcoupled.queue']; % filename for coupled queue file
	else
		fname=[prefix '_runocean.queue']; % filename for uncoupled queue file
	end
	% build string for resource allocation
	resourcestring=buildresourcestring(ncpus,'bro');
	% set ouput and err files
	outlogfname = ['run' prefix '.outlog'];
   errlogfname = ['run' prefix '.errlog'];
	% set the modules we need
	modules = {'mpi-hpe/mpt','comp-intel','hdf5/1.8.18_mpt hdf4/4.2.12 netcdf/4.4.1.1_mpt'};
	if iscoupled
		modules = {modules{:},'matlab/2022b','petsc/3.17.3_intel_mpt_py'};
   end
	modulelines = strcat({'module load '},modules);

	% print the inputs
	disp(['Preparing queue file:']);
	disp(['  rundir:    ' rundir]);
	disp(['  grouplist: ' grouplist]);
	disp(['  ncpus:     ' num2str(ncpus)]);
	disp(['  queuename: ' queuename]);
	disp(['  walltime:  ' char(walltime)]);
	if ishelloworld
      disp(['  ishelloworld: ' num2str(ishelloworld)]);
   else
		disp(['  iscoupled: ' num2str(iscoupled)]);
	end
	if iscoupled
		disp(['  mccdir:    ' mccdir]);
	end
	disp(['  fname:     ' fname]);
	disp(['  resources: ' resourcestring]);
	
	%write the .queue file 
	disp(['Writing .queue file ' fname]);
	lines =	{...
		'#PBS -S /bin/bash', ...
		['#PBS -l ' resourcestring], ...
		['#PBS -q ' queuename], ...
		['#PBS -l walltime=' char(walltime)], ...
		'#PBS -m e', ...
		['#PBS -W group_list=' grouplist], ...
		['#PBS -o ' fullfile(rundir,outlogfname)], ...
		['#PBS -e ' fullfile(rundir,errlogfname)], ...
		'', ...
		'. /usr/share/modules/init/bash', ...
		'', ...
		'#load modules', ...
		modulelines{:}, ...
		'',...
		'#Export some variables', ...
		['export PATH=''' getenv("PATH") ':.'''], ...
		'export MPI_LAUNCH_TIMEOUT=800', ...
		'export MPI_GROUP_MAX=800'};
	if ishelloworld
      lines = [lines {...
         '',...
         ['cd ' rundir], ...
         '',...
         'echo "Hello, World"'}];
   elseif iscoupled
		mcc_command='./run_MCCexecutable.sh';
		lib_command=['/nasa/netcdf/4.4.1.1_mpt/lib:',...
			getenv("ISSM_DIR") '/lib:',...
			getenv("PETSC_DIR") '/lib:',...
			getenv("MPI_ROOT") '/lib:',...
			getenv("MKLROOT") '/lib/intel64_lin:',...
			getenv("MKLROOT") '/../compiler/lib/intel64_lin:',...
			getenv("ISSM_DIR") '/externalpackages/triangle/install/lib:',...
			'/nasa/matlab/2022b'];
		mcc_argin=[mccargin{1} ' ' mccargin{2}];
		lines = [lines {...
			'',...
			'#ISSM stuff', ...
			'export ISSM_DIR="/nobackup/bgetraer/trunk-jpl"', ...
			['source /nobackup/bgetraer/trunk-jpl/etc/environment.sh'], ...
			'#move to the run directory, link the MCC files', ...
			['cd ' rundir], ...
			['ln -s ' fullfile(mccdir,'run_MCCexecutable.sh') ' ./'], ...
			['ln -s ' fullfile(mccdir,'MCCexecutable') ' ./'], ...
			'', ...
			'#run the runcouple executable with the envfile input',...
			[mcc_command ' ' lib_command ' ' mcc_argin]}];
	else
		lines = [lines {...
         '', ...
			['cd ' rundir], ...
			'', ...
			'#run the MITgcm executable with MPI', ...
			['mpirun -np ' num2str(ncpus) ' ./mitgcmuv > out 2> err']}];
	end
	fid=fopen(fullfile(rundir,fname),'w+');
	fprintf(fid,'%s\n',lines{:});
	fclose(fid);
	function resourcestring=buildresourcestring(ncpus,nodemodel) % {{{
		% determine how many cpus to use per node
		switch nodemodel
			case 'bro' % broadwell node
				cpupernode = 25;
			otherwise 
				error(['nodemodel ' nodemodel ' is not supported by this queue script yet.']);
		end
		% divide number of processes into whole nodes and partial node
		wholenodes=floor(ncpus/cpupernode);
		partialnodecpus=rem(ncpus,cpupernode);
		% build the string for PBS script
		% make the string for the partial nodes
		partialnode_string='';
		if partialnodecpus>0
			partialnode_string=['1:ncpus=' num2str(partialnodecpus) ':model=' nodemodel];
			% add plus sign if needed
			if wholenodes>0
				partialnode_string=[partialnode_string '+'];
			end
		end
		% make the string for the whole nodes
		wholenode_string='';
		if wholenodes>0
			wholenode_string=[num2str(wholenodes) ':ncpus=' num2str(cpupernode) ':model=' nodemodel];
		end
		% assemble the resource allocation string
		resourcestring=['select=' partialnode_string wholenode_string];
	end
	% }}}
end % }}}
function write_sizefile(fname,SZ) % {{{
	% Reads from SZ.reffile and writes to fname
	% INPUT
	%    fname   file to write to 
	%    SZ      struct with fields: sNx,sNy,OLx,OLy,nSx,nSy,nPx,nPy,Nx,Ny,Nr, and reffile
	%     SZ.reffile   MITgcm reference file to read from

	if SZ.Nx~=(SZ.sNx*SZ.nSx*SZ.nPx) | SZ.Ny~=(SZ.sNy*SZ.nSy*SZ.nPy)
		error('MITgcm domain discretization inconsistent');
	end

	disp([' - writing SIZE     file to ' fname]);
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
	disp([' - writing config.  file to ' fname]);
	fileID = fopen(fname,'w');
	fprintf(fileID,'# %s\n',PKGCONF.description); % write description
	fprintf(fileID,'%s\n',PKGCONF.pkg{:}); % write packages
	fclose(fileID);
end % }}}
function write_diagsizefile(fname,DIAG_SZ) % {{{
	% WRITE_DIAGSIZEFILE copies the example DIAGNOSTICS_SIZE.h file from the MITgcm directory and 
	% changes the numDiags parameter to DIAG_SZ.numDiags 
	disp([' - writing DIAG_SZ  file to ' fname]);
	writeID=fopen(fname,'w');
   readID=fopen(DIAG_SZ.reffile,'r');
	 % read through the template file, write to the new file
   formatSpec='%s\n'; % new line after each string is written
	% read through the commented header, write to new file {{{
	tline = fgetl(readID);
   while tline(1)=='C'
      fprintf(writeID,formatSpec,tline);
      tline = fgetl(readID);
   end %}}}
	% read through the variable declarations, write to new file {{{
	while contains(tline,'INTEGER')
		fprintf(writeID,formatSpec,tline);
		tline = fgetl(readID);
	end % }}}
 % read through the variable values, write to new file {{{
   while contains(tline,'PARAMETER')
		if contains(tline,'numDiags')
			tline = ['      PARAMETER( numDiags = ' num2str(DIAG_SZ.numDiags) ' )'];
		end
      fprintf(writeID,formatSpec,tline);
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
function write_optionsfile(fname,OPT) % {{{
	% WRITE_OPTIONSFILE copies the example OPT.reffile from the MITgcm directory and 
	% changes the defined/undefined options as requested
	disp([' - writing OPTTIONS file to ' fname]);
	writeID=fopen(fname,'w');
   readID=fopen(OPT.reffile,'r');
	 % read through the template file, write to the new file
   formatSpec='%s\n'; % new line after each string is written
	% read through the entire file, write to new file
	tline = fgetl(readID);
   while isstr(tline)
		if startsWith(tline,'#')
			tlinesplit = strsplit(tline,' '); % split the setting and the parameter 
			if contains(tlinesplit{2},OPT.define)
				tline = ['#define ' tlinesplit{2}];
			elseif contains(tlinesplit{2},OPT.undef)
				tline = ['#undef ' tlinesplit{2}];
			end
		end
		fprintf(writeID,formatSpec,tline);
		tline = fgetl(readID);
	end		
	fclose(writeID);
	fclose(readID);
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
	disp(['    writing namelist file to ' fname])
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
		fprintf(fileID,'  %s=%s,\n',fields{i},num2str(val));
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
					fprintf(fileID,'  %s=%s\n',LHS,RHS); % write to file
				case 'levels'
					error('diag. levels not supported');
				otherwise
					LHS=[subfields{i} '(' num2str(n) ')'];
					RHS=num2str(getfield(N(n),subfields{i}));
					fprintf(fileID,'  %s=%s,\n',LHS,RHS); % write to file
			end
		end
		fprintf(fileID,'\n'); % line break
	end
end % }}}
function write_binfile(fname,field) % {{{
	% write to binary input file
	fid = fopen(fname,'w','b'); 
	fwrite(fid,field,'real*8'); 
	fclose(fid); 
end % }}}
function [nc_data, nc_x, nc_y, nc_z]= load_ncdata(fname,varname) % {{{
	nc_x = ncread(fname,'x');   % load the x coordinates (m)
	nc_y = ncread(fname,'y');   % load the y coordinates (m)
	nc_z = ncread(fname,'z');   % load the z coordinates (m)
	nc_data = ncread(fname,varname); % load the data field
end % }}}
