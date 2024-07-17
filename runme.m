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

steps=[];

% experiment control{{{
% 'CLIM'   monthly climatology 2001-2012 from Paris 2
% 'RCP85'  ensemble average forcings from ISMIP-6
% 'PARIS2' ensemble average forcings from Paris 2
expname='CLIM'; % CLIM, RCP85, PARIS2
expprefix=sprintf('PROPHET_%s_',expname); %set model name prefix
% }}}
% directory structure {{{
mitgcm_dir='/nobackup/bgetraer/MITgcm'; % MITgcm directory
proph_dir ='/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET'; % base directory for this project
% define base directories {{{
% make code directory if needed
% this will hold the current compilation of MITgcm
if ~exist(fullfile(proph_dir,'code'))
   mkdir(fullfile(proph_dir,'code'));
end
% make experiments directory if needed
% this will hold subdirectories for each model run
if ~exist(fullfile(proph_dir,'experiments'))
   mkdir(fullfile(proph_dir,'experiments'));
end
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
% Ocean Model
% MITgcm filenames {{{
% compile-time files
sizefile='./code/SIZE.h';
pkgconffile='./code/packages.conf';
% run-time namelist files
datafile='./input/data';
dataobcsfile='./input/data.obcs';
datashelficefile='./input/data.shelfice';
datacalfile='./input/data.cal';
dataexffile='./input/data.exf';
datadiagfile='./input/data.diagnostics';
datapkgfile='./input/data.pkg';
eedatafile='./input/eedata';
% binary input files
bathyfile='bathy.bin';
thetainitfile='theta.init';
saltinitfile='salt.init';
etainitfile='eta.init';
% binary OBCS files
uvelOBSfile='uvel.obs';
vvelOBSfile='vvel.obs';
thetaOBSfile='theta.obs';
saltOBSfile='salt.obs';
uvelOBWfile='uvel.obw';
vvelOBWfile='vvel.obw';
thetaOBWfile='theta.obw';
saltOBWfile='salt.obw';
% binary SHELFICE files
shelficetopofile='shelficetopo.bin';
shelficemassfile='shelficemass.bin';
shelficedmdtfile='shelficedmdt.bin';
shelficesubglfluxfile='shelficesubglflux.bin';
% }}}
% coordinate system {{{
% A cartesian coordinate system is used with x and y following the coordinates of the Polar Stereographic projection.
% At the Amundsen sea, the projection places North roughly in the -x direcion, East in +y, South in +x, West in -y.
% Where MITgcm convention denotes a "north" or "south" (for instance in the OBCS package) these actually refer to Cartesisan 
% +y and -y respectively. Because this model implements an f-plane approximation for the Coriolis force, no geographic sense of north
% or south is nessecary, and the domain generally follows a "North up" "y up" convention, despite the fact that in our implementation 
% that direction is not geographic north.

% ocean domain
x0_OC=-1700000; % x origin of the active ocean domain (m)
y0_OC=-712000;  % y origin of the active ocean domain (m)
Lx=350E3; % length of total ocean domain in x (m)
Ly=518E3; % length of total ocean domain in y (m)
Lz=1500;  % height of active ocean domain in z (m)

dx=1E3; % hor. resolution in x (m)
dy=dx;  % hor. resolution in y (m)
dz=50;  % ver. resolution in z (m)

% SETUP:
%	5 nodes, 140 processes, 10x14 tiles, 35x37 cells, 350E3x518E3 m domain
nWx=1; % number of wall cells in x
nWy=1; % number of wall cells in y
sNx=13;  % Number of X points in tile
sNy=15;  % Number of Y points in tile
%
4LxOC=Lx-dx*nWx; % length of active ocean domain in x (m)
%LyOC=Ly-dy*nWy; % length of active ocean domain in y (m)

nWw=ceil(nWx/2);  % number of wall cells before x=x0_OC (excess placed on east)
nWs=ceil(nWy/2);  % number of wall cells before y=y0_OC (excess placed on north)
x0=x0_OC-nWw*dx; % origin of entire domain in x (m)
y0=y0_OC-nWs*dy; % origin of entire domain in y (m)

NxOC=LxOC/dx;  % number of ocean cells in x
NyOC=LyOC/dy;  % number of ocean cells in y
Nx=NxOC+nWx;   % number of total cells in x
Ny=NyOC+nWy;   % number of total cells in y
Nz=Lz/dz;      % number of total cells in z

xp=((0:Nx)*dx)+x0; % location of all cell edges in x (m)
yp=((0:Ny)*dy)+y0; % location of all cell edges in y (m)
zp=(0:Nz)*dz;      % location of all cell edges in z (m)
zg=zp(1:end-1);   % location of z edge points, not including the lowest (redundant) boundary (m)

xc=xp(1:end-1)+0.5*dx; % location of cell centers in x (m)
yc=yp(1:end-1)+0.5*dy; % location of cell centers in y (m)
zc=zp(1:end-1)+0.5*dz; % location of cell centers in z (m)

xcOC=xc(xc>=0 & xc<=LxOC); % location of filled ocean cell centers in x (m)
ycOC=yc(yc>=0 & yc<=LyOC); % location of filled ocean cell centers in y (m)

% horizontal grid center points
[XC YC]=meshgrid(xc,yc);
XC=XC(:);
YC=YC(:);
% }}}
% time stepping {{{
nsteps=1;               % number of coupled time steps to take
coupledTimeStep=24*60*60; % coupling time step (s)
y2d=360; % use 12 months of 30 days each ('model' calendar, see input/data.cal) (d/yr)
y2s=y2d*24*60*60; % y2s using 'model' calendar
cal_t0=datetime(2013,1,1,0,0,0); % starting calendar datetime
% }}}
% Run-time options
% input/eedata {{{
EEP=struct; % initialize EEDATA structure
% EEPARMS: Execution Environment Parameters {{{
EEP.header='EEPARMS';
% }}}
write_datafile(eedatafile,{EEP},'EXECUTION ENVIRONMENT PARAMETERS');
% }}}
% input/data {{{
P1=struct;P2=struct;P3=struct;P4=struct;P5=struct; % initialize PARM structures
% PARM01: Continuous equation parameters {{{
% structure information
P1.header='PARM01';
P1.description='Continuous equation parameters';

% physical constants
P1.rhoconst = 1030; % kg/m^3
P1.gravity=9.81; % m/s^2

% equation of state
P1.eostype='''JMD95Z''';
P1.tRef = [num2str(Nz) '*' num2str(-1.9)]; % vertical array reference potential temp. (deg C)
P1.sRef = [num2str(Nz) '*' num2str(34.4)]; % vertical array reference salin. (g/kg)
P1.rhoNil=1000; % ref. density (kg/m^3)
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
P1.nonlinFreeSurf=4; % What option is this?
useRealFreshWaterFlux = '.TRUE.'; % Conserve volume with freshwater flux (changes free surface/sea level) 

% full and partial cells
P1.hFacMin = 0.2; % minimum fractional height of partial gridcell
P1.hFacSup=2.0;   % maximum fractional height of surface cell
P1.hFacInf = 0.2; % minimum fractional height of surface cell

% Momentum Equations
P1.vectorInvariantMomentum = '.TRUE.'; % use vector-invariant form of momentum equations
P1.implicitViscosity = '.TRUE.'; % compute vertical diffusive fluxes of momentum implicitly
P1.viscAr=1.E-4;    % Vertical Eddy Viscosity m^2/s
viscAhscheme='constant'; % choose between constant horizontal viscosity and gridscale/flow aware viscosity
switch viscAhscheme
	case 'constant'
		P1.viscAh=10.0; % Horizontal Eddy Viscosity m^2/s
	case 'modifiedLeith'
		P1.viscAhGrid=.01;  % non-dimensional Laplacian grid viscosity parameter (background, constant)
		%P1.viscA4Grid=.001; % non-dimensional bi-harmonic grid viscosity parameter (background, constant)
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

% options set by Dan. I do not know what they do
%P2.cg2dmincolumneps = 5;
%P2.pReleaseVisc=2;
%P2.thincolDamp=100.;
% }}}
% PARM03: Time stepping parameters {{{
% structure information
P3.header='PARM03';
P3.description='Time stepping parameters';

% Run Start and Duration
P3.nIter0=0;        % starting timestep iteration number
P3.deltaT=100.;     % model time step (s) 
P3.nTimeSteps=(coupledTimeStep/P3.deltaT); % number of model clock timesteps to execute

% Ocean Convection
P3.cAdjFreq = -1;   % <0 sets frequency of convective adj. scheme to deltaT (s)

% Restart/Pickup Files
P3.pChkPtFreq=coupledTimeStep; % permanent restart/pickup checkpoint file write interval (s)

% Frequency/Amount of Output
P3.monitorFreq=coupledTimeStep; % interval to write monitor output (s)
P3.monitorSelect=1; % 1: dynamic variables only, 2: add vorticity variables, 3: add surface variables
P3.dumpInitAndLast='.FALSE.'; % write out initial and last iteration model state (OFF)

% Unused settings
%% IN CHANNEL SCRIPT THIS HAD BEEN SET TO coupledTimeStep/deltaT + 1 BUT I THINK UNNECESSARY EXTRA STEP
%P3.chkPtFreq=coupledTimeStep*10; % rolling restart/pickup checkpoint file write interval (s)
%P3.dumpFreq=coupledTimeStep*10; % interval to write model state/snapshot data (s)
% }}}
% PARM04: Gridding parameters {{{
% structure information
P4.header='PARM04';
P4.description='Gridding parameters';

% coordinate system
P4.usingCartesianGrid='.TRUE.'; % use Cartesian coordinates 
P4.delr=[num2str(Nz) '*' num2str(dz)]; % vertical grid spacing 1D array (m) 
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
P5.bathyfile=['''' bathyfile '''']; % filename for 2D ocean bathymetry (m)
P5.hydrogthetafile=['''' thetainitfile '''']; % filename for 3D specification of initial potential temperature (deg C)
P5.hydrogsaltfile=['''' saltinitfile '''']; % filename for 3D specification of initial salinity (g/kg)
P5.psurfinitfile=['''' etainitfile '''']; % filename for 2D specification of initial free surface position (m)
% }}}
write_datafile(datafile,{P1,P2,P3,P4,P5},'MODEL PARAMETERS');
% }}}
% input/data.pkg {{{
PKG=struct; % initialize PACKAGES structure
% PACKAGES: run-time flags for packages to use {{{
PKG.header='PACKAGES';

PKG.useDiagnostics='.TRUE.';
PKG.useOBCS='.TRUE.';
PKG.useShelfIce='.TRUE.';
PKG.useCAL='.FALSE.';
PKG.useEXF='.FALSE.';
% }}}
write_datafile(datapkgfile,{PKG},'Packages');
% }}}
% Run-time package options
% input/data.obcs {{{
OBCS_P1=struct;OBCS_P2=struct;OBCS_P3=struct; % initialize OBCS_PARM structures
% OBCS_PARM01: Open boundaries {{{
OBCS_P1.header='OBCS_PARM01';
OBCS_P1.description='Open boundaries';

% Southern Boundary (bottom y boundary, not geographic south)
OB_J=nWs+1; % j index wrt Ny where the southern OBCS is set
OBCS_P1.OB_Jsouth=[num2str(Nx) '*' num2str(OB_J)]; % Nx-vector of J-indices (w.r.t. Ny) of Northern OB at each I-position (w.r.t. Nx)
OBCS_P1.OBSuFile=['''' uvelOBSfile '''']; % Nx by Nz matrix of u velocity at Northern OB
OBCS_P1.OBSvFile=['''' vvelOBSfile '''']; % Nx by Nz matrix of v velocity at Northern OB
OBCS_P1.OBStFile=['''' thetaOBSfile '''']; % Nx by Nz matrix of pot. temp. at Northern OB
OBCS_P1.OBSsFile=['''' saltOBSfile '''']; % Nx by Nz matrix of salin. at Northern OB

% Western Boundary (lefthand x boundary, not geographic west)
OB_I=nWw+1; % i index wrt Nx where the western OBCS is set
OBCS_P1.OB_Iwest=[num2str(Ny) '*' num2str(OB_I)];  % Ny-vector of I-indices (w.r.t. Nx) of Western OB at each J-position (w.r.t. Ny)
OBCS_P1.OBWuFile=['''' uvelOBWfile '''']; % Nx by Nz matrix of u velocity at Western OB
OBCS_P1.OBWvFile=['''' vvelOBWfile '''']; % Nx by Nz matrix of v velocity at Western OB
OBCS_P1.OBWtFile=['''' thetaOBWfile '''']; % Nx by Nz matrix of pot. temp. at Western OB
OBCS_P1.OBWsFile=['''' saltOBWfile '''']; % Nx by Nz matrix of salin. at Western OB

OBCS_P1.useOBCSprescribe='.TRUE.'; % prescribe OB conditions
OBCS_P1.useOBCSbalance='.TRUE.'; % use OB balance
OBCS_P1.OBCS_balanceFacN=1.0; % balance factor for N OB
OBCS_P1.OBCS_balanceFacW=1.0; % balance factor for W OB

% Sponge layer
useOBCSsponge=0; % on or off
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
	OBCS_P3.Vrelaxobcsbound=coupledTimeStep*10; % relaxation time scale at the outermost sponge layer point of a zonal OB (s)
	OBCS_P3.Urelaxobcsbound=coupledTimeStep*10; % relaxation time scale at the outermost sponge layer point of a meridional OB (s) 
end
% }}}
write_datafile(dataobcsfile,{OBCS_P1,OBCS_P2, OBCS_P3},'OBCS RUNTIME PARAMETERS');
% }}}
% input/data.shelfice {{{
SHELFICE_PARM01=struct; % initialize SHELFICE_PARM structure
% SHELFICE_PARM01: Parameters for SHELFICE package {{{
SHELFICE_P1.header='SHELFICE_PARM01';

% general options
SHELFICE_P1.SHELFICEwriteState='.TRUE.'; % write ice shelf state to file
SHELFICE_P1.SHELFICEconserve='.TRUE.'; % use conservative form of temperature boundary conditions
SHELFICE_P1.SHELFICEMassStepping = '.TRUE.'; % recalculate ice shelf mass at every time step

% input files
SHELFICE_P1.SHELFICEtopoFile=['''' shelficetopofile '''']; % filename for under-ice topography of ice shelves
SHELFICE_P1.SHELFICEmassFile=['''' shelficemassfile '''']; % filename for mass of ice shelves
SHELFICE_P1.SHELFICEMassDynTendFile=['''' shelficedmdtfile '''']; % filename for mass tendency of ice shelves

% boundary layer options
SHELFICE_P1.SHELFICEboundaryLayer='.TRUE.'; % use simple boundary layer mixing parameterization
SHELFICE_P1.SHI_withBL_realFWflux='.TRUE.'; % use real-FW flux from boundary layer
%SHIELFICE_P1.SHI_withBL_uStarTopDz='.TRUE.'; % compute uStar from uVel, vVel averaged over top Dz thickness

% thermodynamic exchange options
SHELFICE_P1.SHELFICEuseGammaFrict='.TRUE.'; % use velocity dependent exchange coefficients (Holland and Jenkins 1999) 
SHELFICE_P1.SHELFICEheatTransCoeff=0.0135; % transfer coefficient for temperature (m/s)
SHELFICE_P1.SHELFICEsaltTransCoeff=0.000265; % transfer coefficient for salinity (m/s)

% drag options
SHELFICE_P1.SHELFICEDragQuadratic=.006; % quadratic drag coefficient at bottom ice shelf (non-dim.)
SHELFICE_P1.shiCdrag=SHELFICE_P1.SHELFICEDragQuadratic; % set to be the same

% ice draft free-surface remeshing options
SHELFICE_P1.SHELFICERemeshFrequency=3600.0; % how often to remesh the surface cells (s)
SHELFICE_P1.SHELFICESplitThreshold=1.3; % maximum fractional height of ice shelf cavity surface cell
SHELFICE_P1.SHELFICEMergeThreshold=.29; % minimum fractional height of ice shelf cavity surface cell 

% subglacial discharge options
addrunoff=0;
if addrunoff==1
	SHELFICE_P1.SHELFICEaddrunoff='.TRUE.'; % use runoff
	SHELFICE_P1.SHELFICESubglFluxFile=['''' shelficesubglfluxfile '''']; % filename for subglacial runoff specification
end

% Dan's options, i think related to the modified melt rate at depth
%SHELFICE_P1.SHELFICE_massmin_trueDens='.TRUE.';
%SHELFICE_P1.SHELFICE_conserve_ssh='.TRUE.';
%SHELFICE_P1.SHELFICEdepthMinMelt=10.;
%SHELFICE_P1.SHELFICE_transition_gamma='.FALSE.'
%SHELFICE_P1.SHELFICETransGammaThickness=0.
% }}}
write_datafile(datashelficefile,{SHELFICE_P1},'SHELFICE RUNTIME PARAMETERS');
% }}}
% input/data.cal {{{
CAL=struct; % intialize CAL structure
% CAL_NML: The calendar package namelist {{{
CAL.header='CAL_NML';

CAL.TheCalendar='model'; % choose 'model' calendar
CAL.startdate_1=string(cal_t0,'yyyyMMdd'); % YYYYMMDD of start date
CAL.startDate_2=string(cal_t0,'hhmmss');   % HHMMSS of start date
% Benjy: I am not sure exactly what calendarDumps does.
CAL.calendarDumps='.FALSE.'; % align the output with calendar months?
% }}}
write_datafile(datacalfile,{CAL},'CALENDAR PARAMETERS');
% }}}
% input/data.exf {{{
switch experiment.name
	case 'CLIM'   % monthly climatology 2001-2012 from Paris 2
	case 'RCP85'  % ensemble average forcings from ISMIP-6
	case 'PARIS2' % ensemble average forcings from Paris 2
end
EXF1=struct; EXF2=struct; EXF3=struct; EXF4=struct; EXF_OBCS=struct; % intialize EXF structures
% EXF_NML_01: External Forcings namelist 1 {{{
EXF1.header='EXF_NML_01';
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
EXF_OBCS.header='EXF_NML_EXF_OBCS';
% }}}
write_datafile(dataexffile,{EXF1,EXF2,EXF3,EXF4,EXF_OBCS},'EXTERNAL FORCINGS PARAMETERS');
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

% Output Stream 1: surfDiag
DIAG_LIST.N(1).filename='''surfDiag''';
DIAG_LIST.N(1).frequency=86400;
DIAG_LIST.N(1).fields={'SHIfwFlx','SHI_mass','SHIRshel','ETAN    ','SHIuStar','SHIForcT'};

% Output Stream 2: dynDiag
DIAG_LIST.N(2).filename='''dynDiag''';
DIAG_LIST.N(2).frequency=2592000.;
DIAG_LIST.N(2).fields={'UVEL    ','VVEL    ','WVEL    ','THETA   ','SALT    '};

% Output Stream 3: meltDiag
DIAG_LIST.N(3).filename='''meltDiag''';
DIAG_LIST.N(3).frequency=864000.;
DIAG_LIST.N(3).fields={'SHIfwFlx'};
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
write_datafile(datadiagfile,{DIAG_LIST,DIAG_STATIS},'DIAGNOSTICS RUNTIME PARAMETERS');
% }}}

% Ice Model
expfile='./Exp/domain.exp';
% ISSM domain setup

org=organizer('repository',modeldir,'prefix',experiment.prefix,'steps',steps);
if perform(org,'BuildMITgcm'), % Compile MITgcm{{{
	% Compile-time options {{{
	% code/SIZE.h {{{
	SZ=struct; % initialize SIZE.h structure
	SZ.reffile=fullfile(mitgcm_dir,'model/inc/SIZE.h'); % MITgcm example file
	% define SZ structure fields
	% note domain decomposition must follow: Nx= sNx*nSx*nPx, Ny = sNy*nSy*nPy
	C     sNx :: Number of X points in tile.
	C     sNy :: Number of Y points in tile.
	C     OLx :: Tile overlap extent in X.
	C     OLy :: Tile overlap extent in Y.
	C     nSx :: Number of tiles per process in X.
	C     nSy :: Number of tiles per process in Y.
	C     nPx :: Number of processes to use in X.
	C     nPy :: Number of processes to use in Y.
	C     Nx  :: Number of points in X for the full domain.
	C     Ny  :: Number of points in Y for the full domain.
	C     Nr  :: Number of points in vertical direction.
	SZ.sNx=1;   % Number of X points in tile.
	SZ.sNy=1;   % Number of Y points in tile.              
	SZ.OLx=3;   % Tile overlap extent in X.                
	SZ.OLy=3;   % Tile overlap extent in Y.                
	SZ.nSx=1;   % Number of tiles per process in X.        
	SZ.nSy=1;   % Number of tiles per process in Y.        
	SZ.nPx=100; % Number of processes to use in X.         
	SZ.nPy=100; % Number of processes to use in Y.         
	SZ.Nx =Nx;  % Number of points in X for the full domain
	SZ.Ny =Ny;  % Number of points in Y for the full domain
	SZ.Nr =Nz;  % Number of points in vertical direction.
	% write to ./code/SIZE.h
	write_sizefile(sizefile,SZ);
	% }}}
	% code/packages.conf {{{
	PKGCONF=struct; % initialize Package Configuration structure
	% define Package Configuration structure fields
	PKGCONF.description='Configuration File for Package Inclusion';
	PKGCONF.pkg={'gfd','obcs','shelfice','cal','exf','diagnostics'};
	% write to ./code/packages.conf
	write_pkgconffile(pkgconffile,PKGCONF);
	% }}}
	% }}}
	% Compile {{{
	% locate files and scripts
	genmake2='${MITGCM_ROOTDIR}/tools/genmake2';
	optfile='${MITGCM_ROOTDIR}/tools/build_options/linux_amd64_ifort+mpi_ice_nas';
	% clear the build directory
	cd(proph_dir);
	!rm -r ./build
	mkdir('build');
	cd('./build');
	% make the MITgcm executable
	command=[genmake2 ' -mpi -mo ../code -optfile ' optfile ' -rd ${MITGCM_ROOTDIR}'];
	system(command); % generate Makefile
	system('make CLEAN');  % prepare for new compilation
	system('make depend'); % create symbolic links from the local directory to the source file locations
	system('make');        % compile code and create executable file mitgcmuv
	cd(proph_dir);
	% }}}
end % }}}
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
		hVertices(find(vel>200)) = fine_vel;
		md=bamg(md,'gradation',1.6180,'anisomax',1.e6,'KeepVertices',0,'Hessiantype',0,'Metrictype',0,...
			'hmax',coarse,'hmin',fine_vel,'hVertices',hVertices,'field',vel,'err',3);

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
	md.friction.m = 1/3; % assume to follow Weertman (non-dimensional)
	
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

   % set M1QN3 package
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
   mds.cluster=cluster;
   mds.verbose=verbose('solution',false,'control',true);
   mds.miscellaneous.name='inversion_B';
   mds.friction.coefficient(:)=0; % make sure there is no basal friction

   % Solver parameters
   mds.toolkits.DefaultAnalysis=bcgslbjacobioptions();% biconjugate gradient with block Jacobi preconditioner
   mds.toolkits.DefaultAnalysis.ksp_max_it=500;
   mds.settings.solver_residue_threshold=1e-6;

	% Solve
   mds=solve(mds,'Stressbalance'); % only extracted model

   % Update the full model rheology_B accordingly
   md.materials.rheology_B(mds.mesh.extractedelements) = mds.results.StressbalanceSolution.MaterialsRheologyBbar;
   savemodel(org,md);
end % }}}
if perform(org,'InversionC'),  % Invert for friction coefficient C {{{
   md=loadmodel(org,'InversionB');

   
	% Friction Law
   disp('   -- Setting up a Budd''s sliding law (with m=1, q=1)');
   md.friction.p=ones(md.mesh.numberofelements,1);
   md.friction.q=ones(md.mesh.numberofelements,1);
   md.inversion.cost_functions=[101 103 501];
   md.inversion.cost_functions_coefficients = ones(md.mesh.numberofvertices,3); % IMPORTANT: set as 0 again the vertices with no data (NaN)
   md.inversion.cost_functions_coefficients(:,1)=2.2452e+03;% for budd change this such that the cost functions 101 and 103 have the same contribution
   md.inversion.cost_functions_coefficients(:,2)=1; % always 1
   md.inversion.cost_functions_coefficients(:,end)=1e-8;  % for budd
	load(fullfile(modeldir, 'pos_vel_nan.mat'),'pos_vel_nan'); % load locations of NaN data
   md.inversion.cost_functions_coefficients(pos_vel_nan,:)=0; % do not try to fit positions with NaN in the velocity data set

   %Controls
   md.friction.coefficient=EstimateFric_Budd(md); % initial guess from Driving Stress (using r=1 and s=1)
   md.inversion.control_parameters={'FrictionCoefficient'};
   md.inversion.min_parameters=0*ones(md.mesh.numberofvertices,1);
   md.inversion.max_parameters=3000*ones(md.mesh.numberofvertices,1);
   md.miscellaneous.name='inversion_drag_budd';

   % Control
   md.inversion.iscontrol=1;
   md.inversion.maxsteps=80;
   md.inversion.maxiter=80;
   md.inversion.dxmin=0.01;
   md.inversion.gttol=1.0e-6;
   md.inversion.incomplete_adjoint=0; % 0: non linear viscosity, 1: linear viscosity 04/29/2019 changed to non linear

   % Keep friction 0 in the floating part
   pos=find(md.mask.ocean_levelset<0);
   md.friction.coefficient(pos)=0;
	 % Additional parameters for stress balance
   md.stressbalance.restol=0.002; % 04/29/2019
   %md.stressbalance.reltol=0.01; % 04/29/2019
   %md.stressbalance.abstol=10; % 04/29/2019
   %md.stressbalance.maxiter=30; % 04/29/2019
   % updated?
   md.stressbalance.maxiter=50; % 10/24/2019
   md.stressbalance.reltol=NaN; % 11/05/2019
   md.stressbalance.abstol=NaN; % 11/05/2019

   % Set cluster
   md.cluster=cluster;
   md.verbose=verbose('solution',false,'control',true);

   % Define solver
   md.toolkits.DefaultAnalysis=bcgslbjacobioptions();% biconjugate gradient with block Jacobi preconditioner
   md.toolkits.DefaultAnalysis.ksp_max_it = 500;

   md=solve(md,'Stressbalance');

   % Update friction coefficient and velocity field
   md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;
   md.initialization.vx= md.results.StressbalanceSolution.Vx;
   md.initialization.vy=md.results.StressbalanceSolution.Vy;
   md.initialization.vel=md.results.StressbalanceSolution.Vel;

   % Get the modeled velocity field that N3 and N4 models have in common
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
