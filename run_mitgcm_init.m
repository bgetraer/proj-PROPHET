% PROPHET Amundsen Sea Coupling
% Outline:
%  Compile MITgcm
%	Initialize MITgcm
%

steps=[4];

experiment.name='MITgcm_initialization';
experiment.forcing = 'CLIM'; % 'CLIM', 'RCP85', or 'Paris2C'
% directory structure {{{
%mitgcm_dir='/nobackup/bgetraer/MITgcm'; % MITgcm directory
proph_dir = pwd; % base directory for this project
% define experiment directories {{{
% make experiments directory if needed
% this will hold subdirectories for each model run
if ~exist(fullfile(proph_dir,'experiments'))
	mkdir(fullfile(proph_dir,'experiments'));
end
% make this experiment directory if needed
expdir=fullfile(proph_dir,'experiments',experiment.name);
if ~exist(expdir)
	mkdir(expdir);
end
% make code directory if needed
% this will hold the current compilation of MITgcm
if ~exist(fullfile(expdir,'code'))
	mkdir(fullfile(expdir,'code'));
end
% make input directory if needed
% this will hold the runtime options for MITgcm
if ~exist(fullfile(expdir,'input'))
	mkdir(fullfile(expdir,'input'));
end
% make model directory if needed
% this will hold md structures from ISSM
modeldir=fullfile(expdir,'Models');
if ~exist(modeldir)
	mkdir(modeldir);
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
% Move into the experiment directory
disp(['Moving to experiment directory: ', expdir]);
cd(expdir); 
if 0
	% Ocean Model
	% climate forcing filenames {{{
	switch experiment.forcing
		case 'CLIM'   % monthly climatology 2001-2012 from Paris 2
		case 'RCP85'  % ensemble average forcings from ISMIP-6
		case 'PARIS2' % ensemble average forcings from Paris 2
		otherwise
			error('');
	end
	% }}}
	% MITgcm filenames {{{
	% compile-time files
	sizefile='code/SIZE.h';
	pkgconffile='code/packages.conf';
	% run-time namelist files
	datafile='input/data';
	dataobcsfile='input/data.obcs';
	datashelficefile='input/data.shelfice';
	datacalfile='input/data.cal';
	dataexffile='input/data.exf';
	datadiagfile='input/data.diagnostics';
	datapkgfile='input/data.pkg';
	eedatafile='input/eedata';
	% binary input files
	delrfile = 'delr.bin';
	bathyfile='bathy.bin';
	thetainitfile='theta.init';
	saltinitfile='salt.init';
	etainitfile='eta.init';
	% binary OBCS files
	uvelOBSfile = ['uvel.obs.'  experiment.forcing];
	vvelOBSfile = ['vvel.obs.'  experiment.forcing];
	thetaOBSfile= ['theta.obs.' experiment.forcing];
	saltOBSfile = ['salt.obs.'  experiment.forcing];
	uvelOBWfile = ['uvel.obw.'  experiment.forcing];
	vvelOBWfile = ['vvel.obw.'  experiment.forcing];
	thetaOBWfile= ['theta.obw.' experiment.forcing];
	saltOBWfile = ['salt.obw.'  experiment.forcing];
	% binary SHELFICE files
	shelficetopofile='shelficetopo.bin';
	shelficemassfile='shelficemass.bin';
	shelficedmdtfile='shelficedmdt.bin';
	shelficesubglfluxfile='shelficesubglflux.bin';
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
	%P4.delr=[num2str(Nz) '*' num2str(dz)]; % vertical grid spacing 1D array (m) 
	P4.delRFile=['''' delrfile '''']; % vertical grid spacing 1D array (m) 
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
	OBCS_P1.OB_Jsouth=[num2str(Nx) '*' num2str(OB_J)]; % Nx-vector of J-indices (w.r.t. Ny) of Southern OB at each I-position (w.r.t. Nx)
	OBCS_P1.OBSuFile=['''' uvelOBSfile '''']; % Nx by Nz matrix of u velocity at Southern OB
	OBCS_P1.OBSvFile=['''' vvelOBSfile '''']; % Nx by Nz matrix of v velocity at Southern OB
	OBCS_P1.OBStFile=['''' thetaOBSfile '''']; % Nx by Nz matrix of pot. temp. at Southern OB
	OBCS_P1.OBSsFile=['''' saltOBSfile '''']; % Nx by Nz matrix of salin. at Southern OB

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
	EXF_OBCS.useOBCSYearlyFields = '.TRUE.'; 
	obcsWstartdate2 = 20051216; % W boundary - YYYYMMDD of start date
	obcsWperiod     = -1.0;     % period
	obcsSstartdate1 = 20051216; % S boundary - YYYYMMDD of start date
	obcsSperiod     = -1.0;     % period
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
end
% suppress warnings {{{
wid = {'MATLAB:hg:AutoSoftwareOpenGL','MATLAB:polyshape:repairedBySimplify'};
for i=1:length(wid)
	warning('off',wid{i});
end
% }}}
org=organizer('repository',modeldir,'prefix','PROPHET_mitgcm_init_','steps',steps);
if perform(org,'Mesh'), % {{{
	mit=mitgcmstruct(); % initialize model parameter structure that will store information for MITgcm

	% Vertical discretization
	% The vertical grid is defined by the "cell centered" approach (see MITgcm readthedocs 2.11.5) and is defined by the vector delzF
	% which provides the thickness of the cells. Here, the vertical grid is non-constant, varying from 10m at the surface to a maximum
	% defined by delzF_max at depth. Between these cell thicknesses, the grid spacing is lifted directly from Naughten et al., 2023.
	fname = fullfile(proph_dir,'climateforcings/naughten2023_data/initpaholPyParis2C001benjy.mat');
	zc_N = getfield(load(fname,'z'),'z');     % Naughten2023 cell center locations in z (m)
	zp_N = zeros(size(zc_N));                 % Initialize vertical cell edge locations in z (m)
	for i=1:length(zp_N)
		zp_N(i+1) = 2*zc_N(i) - zp_N(i);       % Naughten2023 vertical cell edge locations in z (m)
	end
	delzF_N = diff(zp_N);                     % Naughten2023 vertical cell thickness (m)
	delzC_N = diff([zp_N(1) zc_N zp_N(end)]); % Naughten2023 vertical cell center spacing (m)

	% Here, we cap the cell thickness at delzF_max, and use the same cell spacing as Naughten up
	% until that limit, followed by cells of uniform delzF_max thickness until we reach the lowest
	% bathymetry we need to capture.
	zmin = -2.067E3; % minimum bathymetry, from B_dagger, which actually is calculated later.(m)
	% Calculate the vertical cell thickness
	delzF_max = 32.5; % max vertical spacing (m)
	delzF_upper = delzF_N(abs(delzF_N)<=delzF_max); % use the same vertical spacing as Naughten in upper column (m)
	delzF_lower = repmat(-delzF_max,1,abs(round((zmin/delzF_max)))); % use constant spacing of delzF_max in the upper column (m)
	delzF = [delzF_upper delzF_lower]; % vertical cell thickness (m)
	% Calculate the vertical cell edge locations and center spacing
	zp = [0 cumsum(delzF)];
	z_end = find(zp<zmin,1); % index to terminate the vertical grid
	zp = zp(1:z_end);                 % vertical cell edge locations in z (m)
	delzF = delzF(1:z_end-1);         % vertical cell thickness in z (m)
	zc = zp(1:end-1) + 0.5*delzF;     % cell center locations in z (m)
	% Calculate the vertical cell center spacing
	delzC = diff([zp(1) zc zp(end)]); % vertical cell center spacing (m)

	% Horizontal discretization
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
	%	5 nodes, 140 processes, 10x14 tiles, 25x30 cells, 300E3x505E3 m domain
	sNx=25;  % Number of X points in tile
	sNy=30;  % Number of Y points in tile

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

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Save mesh to mit structure
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PROPHET mesh
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
	% Naughten mesh
	mit.mesh.zc_N = zc_N; % Naughten2023 cell center locations in z (m)
	mit.mesh.zp_N = zp_N; % Naughten2023 cell face locations in z (m)
	mit.mesh.delzF_N = delzF_N; % Naughten2023 vertical spacing between cell faces (m)
	mit.mesh.delzC_N = delzC_N; % Naughten2023 vertical spacing between cell centers (m)
	[mit.mesh.XC_N,mit.mesh.YC_N,mit.mesh.ZC_N]=meshgrid(xc,yc,zc_N);  % Naughten2023 cell center locations in x, y, and z (m)

	savedata(org,mit); % save the model
end % }}}
if perform(org,'Param'), % {{{
	mit=load(org,'Mesh');
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% climate forcing filenames
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
	field_prefix = {'T','S','U','V'}; % data fields
	mit.forcing.field_bdry = [strcat(field_prefix,'leftbdry') strcat(field_prefix,'botbdry')]; % field names for bdry files
	mit.forcing.field_init = strcat(field_prefix(1:2),'init'); % field names for init files
	mit.forcing.Ddir = fullfile(proph_dir,'climateforcings/data'); % directory path for processed boundary forcing data
	field_prefix_out = {'theta','salt','uvel','vvel'}; % data fields
	mit.forcing.field_bdry_out = [strcat(field_prefix_out,'.obw'),strcat(field_prefix_out,'.obs')]; % field names for processed forcing files

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% MITgcm filenames
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	mit.fname = struct();

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Bathymetry
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% The bathymetry is interpolated from Bedmachine Antarctica, but adjusted along the bottom and left boundaries of the domain to 
	% match the boundary condition forcing files from Naughten et al., 2023 exactly. This transition is accomplished by a piece-wise 
	% linear transition over a length Lb. Along the boundary, we define the bed as being the bottom of the first real valued cell 
	% when viewed from below. This definition does not address nan data where an ice shelf exists along the boundary, and maintaining 
	% constant ice thickness along the boundary is accomplished elsewhere.
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
	fname = fullfile(mit.forcing.Ndir,mit.forcing.init_example); % one of the init files
	fieldname = mit.forcing.field_init{1}; % one of the init state variable fields (Tinit)
	A = getfield(load(fname,fieldname),fieldname); % load state variable
	A = pagetranspose(reshape(A,size(A,2),size(A,1),size(A,3))); % correct the indexing to column-major
	A(end:length(mit.mesh.yc),:) = NaN; % extend to actual boundary
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
	B = interpBedmachineAntarctica(mit.mesh.hXC,mit.mesh.hYC,'bed','nearest'); % B is the bed from Bedmachine Antarctica

	% The bathymetry interpolation scheme
	B_prime(wallmaskB_prime)=min(B(wallmaskB_prime),0); % do not interpolate over Naughten's wall, just use Bedmachine
	Lb = 10E3; % distance over which to transition the beds (m)
	a1 = min(min(1/Lb.*(mit.mesh.hYC-mit.mesh.yc(1)),1/Lb.*(mit.mesh.hXC-mit.mesh.xc(1))),1); % weighting between the two beds (0--1)
	% The adjusted bathymetry
	B_dagger = min(a1.*B + (1-a1).*B_prime, 0); % Adjusted bathymetry, capped at sea level (m)

	% Our Bear Ridge wall
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
	a1 = min(1/Lb.*(polywall.Vertices(ind,1)-mit.mesh.xc(1)),1);
	polywall.Vertices(ind,2) = a1.*polywall.Vertices(ind,2) + (1-a1).*bdryVertices_N(1,2);
	% fit lower wall to Naughtens wall at boundary
	ind = polywall.Vertices(:,2)<bdryVertices_N(2,2); % fit wall to Naughten's lower wall at boundary 
	a1 = min(1/Lb.*(polywall.Vertices(ind,1)-mit.mesh.xc(1)),1);
	polywall.Vertices(ind,2) = a1.*polywall.Vertices(ind,2) + (1-a1).*bdryVertices_N(2,2);

	% enforce wall in bathymetry
	wallmask = reshape(isinterior(polywall,mit.mesh.hXC(:),mit.mesh.hYC(:)),size(B_prime)); % the Bear Ridge wall mask
	B_dagger(wallmask)=0; 

	% Check if our B_dagger bathymetry accurately masks the correct cells
	mask = (mit.mesh.ZC_N)<B_dagger; % nanmask from B_dagger
	mask_N = isnan(A); % Naughten nanmask
	mismatch = (mask ~= isnan(A)); % the difference
	mask_zc = (mit.mesh.ZC)<B_dagger; % the mask on our z grid

	% Save Bathymetry data
	mit.geometry=struct();
	mit.geometry.bathy=B_dagger; % the bathymetry (m)
	mit.geometry.polywall=polywall; % the polyshape defining the bear ridge wall
	mit.forcing.obs_mask=squeeze(mask_zc(1,:,:))';
	mit.forcing.obw_mask=squeeze(mask_zc(:,1,:))';
	mit.forcing.obs_mask_N=squeeze(mask_N(1,:,:))';
	mit.forcing.obw_mask_N=squeeze(mask_N(:,1,:))';

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Plot
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	figure(1);clf;
	subplot(2,1,1);hold on;
	j = 1;
	h=pcolor(squeeze(mit.mesh.XC_N(j,:,:)),squeeze(mit.mesh.ZC_N(j,:,:)),squeeze(1*(mask_N(j,:,:))));
	set(h,'edgecolor','flat');
	plot(mit.mesh.xc,B_dagger(j,:),'r')
	axis tight; ylabel('z');xlabel('x');
	title('isnan mask of Naughten');
	subtitle('our bed follows beneath floating ice');
	subplot(2,1,2);hold on;
	h=pcolor(squeeze(mit.mesh.XC_N(j,:,:)),squeeze(mit.mesh.ZC_N(j,:,:)),squeeze(1*(mismatch(j,:,:))));
	set(h,'edgecolor','flat');
	plot(mit.mesh.xc,B_dagger(j,:),'r')
	axis tight; ylabel('z');xlabel('x');
	title('nanmask of Naughten != our nanmask');
	subtitle('our nanmask does not include floating ice');

	figure(2);clf;
	subplot(2,1,1);hold on;
	i = 1;
	h=pcolor(squeeze(mit.mesh.YC_N(:,i,:)),squeeze(mit.mesh.ZC_N(:,i,:)),squeeze(1*(mask_N(:,i,:))));
	set(h,'edgecolor','flat');
	plot(mit.mesh.yc,B_dagger(:,i),'r')
	axis tight; ylabel('z');xlabel('x');
	title('isnan mask of Naughten');
	subtitle('our bed follows beneath floating ice');
	subplot(2,1,2);hold on;
	h=pcolor(squeeze(mit.mesh.YC_N(:,i,:)),squeeze(mit.mesh.ZC_N(:,i,:)),squeeze(1*(mismatch(:,i,:))));
	set(h,'edgecolor','flat');
	plot(mit.mesh.yc,B_dagger(:,i),'r')
	axis tight; ylabel('z');xlabel('x');
	title('nanmask of Naughten != our nanmask');
	subtitle('our nanmask does not include floating ice');
	caxis([0,1]);

	% Plot Bear Ridge Wall
	figure(3);clf;
	hold on;
	imagesc(mit.mesh.xc,mit.mesh.yc,B)
	set(gca,'ydir','normal');axis equal tight;
	xlim(mit.mesh.xc([1,end]));
	caxis([min(B,[],[1,2]),0]);
	h0N=plot(polyB_prime,'facecolor','none','edgecolor','r','linewidth',2);
	hwall=plot(polywall,'facecolor','none','edgecolor','k','linewidth',2);
	h300=plot(polyB,'facecolor','none','edgecolor','w','linestyle','--');
	legend([h0N,hwall,h300],'Naughten''s wall','My wall','300m depth contour')
	cb=colorbar;
	cb.Label.String = 'depth (m)';
	title('Bear Ridge wall')
	ylabel('y (m)');xlabel('x (m)');

	% Save mit structure
	savedata(org,mit);
end % }}}
if perform(org,'BoundaryForcings'), % {{{
	mit=load(org,'Param');

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
			if strcmp( mit.forcing.field_bdry_out{k}(end-2:end),'obw')
				B(mit.forcing.obw_mask_N)=NaN;
				% the domain extent is actually further into the grounded ice: pad this with nan to reach the correct size
				B(:,end+1:length(mit.mesh.yc),:)=NaN; % pad the domain with nan 
			else
				B(mit.forcing.obs_mask_N)=NaN;
			end

			% interpolate from Naughten z grid onto the refined z grid
			Bq = interp1(mit.mesh.zc_N,B,mit.mesh.zc,'nearest');

			% plot the vertical grid and bathymetry
			if k==1 & i==1
				figure(1);clf;
				subplot(2,1,1); hold on;
				h=pcolor(mit.mesh.yp(1:end-1),mit.mesh.zp_N(1:end-1),B(:,:,1));
				set(h,'edgecolor','flat');
				xlim(mit.mesh.yc([1,end]));
				ylim([-1500,0]);
				plot(mit.mesh.yc,mit.geometry.bathy(:,1),'r','linewidth',2);
				yline(mit.mesh.zp_N);xline(mit.mesh.yp);
				xlim([-5.6 -4.6]*1E5);ylim([-1000 -500]);
				title('Naughten vertical grid')
				subtitle([mit.forcing.field_bdry{k} ', red is bathymetry'])

				subplot(2,1,2); hold on;
				h=pcolor(mit.mesh.yp(1:end-1),mit.mesh.zp(1:end-1),Bq(:,:,1));
				%h=pcolor(mit.mesh.yc,mit.mesh.zc,1*mit.forcing.obw_mask);
				%h=pcolor(mit.mesh.yc,mit.mesh.zc,1*(A-mit.forcing.obw_mask));
				set(h,'edgecolor','flat');
				xlim(mit.mesh.yc([1,end]));
				ylim([-1500,0]);
				plot(mit.mesh.yc,mit.geometry.bathy(:,1),'r','linewidth',2);
				yline(mit.mesh.zp);xline(mit.mesh.yp);
				xlim([-5.6 -4.6]*1E5);ylim([-1000 -500]);
				title('My vertical grid')
				subtitle([mit.forcing.field_bdry{k} ', red is bathymetry'])
			end

			% divide B into one file per year and write to binary file for MITgcm
			years = D{j}.ystart:D{j}.yend;
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
				write_binfile(fullfile(mit.forcing.Ddir,fname_out),Bq(:,:,ind));
			end
		end
	end
end % }}}
if perform(org,'InitialConditions'), % {{{
	mit=load(org,'Param');
	mit.initialization=struct();
	
	D={}; % Initialize cell array for climate forcing data structures
	for i = 1:length(mit.forcing.exp)
		% Load all members of the experiment
		for j = 1:length(mit.forcing.member)
			fname = [mit.forcing.init_prefix mit.forcing.exp{i} mit.forcing.member{j} mit.forcing.suffix]; % file to be loaded
			fprintf('Loading '); ls(fullfile(mit.forcing.Ndir,fname)); % print out the file being read
			D{j} = load(fullfile(mit.forcing.Ndir,fname)); % load data into array of structures
		end
		% Loop through all fields
      for k=1:length(mit.forcing.field_init)
			 disp(['Processing ' mit.forcing.field_init{k}]);
         % collect field from all members
         A = []; % matrix to collect field. dim1 is z, dim2 is x or y, dim3 is time (month), dim4 is member (1-10)
         for j = 1:length(D)
				A(:,:,:,j) = getfield(D{j},mit.forcing.field_init{k});
         end
			% take mean of field across members and horizontal dimensions
         B = squeeze(nanmean(A,[1,2,4]));
			%sigma = squeeze(nanstd(A,[],[1,2,4]));

			x = mit.mesh.zc_N(~isnan(B));
			b = B(~isnan(B));
			%sigma = sigma(~isnan(B));
			xq = mit.mesh.zc;
			bq = interp1(x,b,xq,'linear','extrap');
			
			% save the average state 
			fieldname = [mit.forcing.field_init{k}(1),'_',mit.forcing.exp{i}];
			mit.initialization.(fieldname) = bq;
		end
	end
	figure(1);clf;
	ax1=subplot(1,2,1); hold on;
	plot(mit.mesh.zc,mit.initialization.T_Paris2C,'linewidth',2);
	plot(mit.mesh.zc,mit.initialization.T_RCP85,'linewidth',2);
	xlabel('z (m)');,ylabel('T (deg C)');
	ax2=subplot(1,2,2); hold on;
	plot(mit.mesh.zc,mit.initialization.S_Paris2C,'linewidth',2);
	plot(mit.mesh.zc,mit.initialization.S_RCP85,'linewidth',2);
	xlabel('z (m)');,ylabel('S (ppt)');
	legend('Paris2C','RCP85','location','sw');

	savedata(org,mit);
end % }}}
if perform(org,'CompileMITgcm'), % compile MITgcm{{{
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
	SZ.sNx=sNx;   % Number of X points in tile.
	SZ.sNy=sNy;   % Number of Y points in tile.              
	SZ.OLx=3;   % Tile overlap extent in X.                
	SZ.OLy=3;   % Tile overlap extent in Y.                
	SZ.nSx=1;   % Number of tiles per process in X.        
	SZ.nSy=1;   % Number of tiles per process in Y.        
	SZ.nPx=XX; % Number of processes to use in X.         
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
if perform(org,'RunPrepMITgcm'), % prepare MITgcm for run {{{
	% input/data binary files {{{
	% define bathymetry {{{ 
	bathymetry = interpBedmachineAntarctica(XC,YC,'bed','linear'); % interpolate bed elevation data from BedMachine (m)
	walls=(XC<x0_OC | XC>(LxOC+x0_OC) | YC<y0_OC | YC>(LyOC+y0_OC)); % index of wall locations
	bathymetry(walls) = 0; % set wall bathymetry (m)
	write_binfile(bathyfile,permute(bathymetry,[2,1])); % write to binary input file with size Nx Ny
	% }}}
	% define initial T, S conditions {{{
	% thetainitfile
	[nc_data, nc_x, nc_y, nc_z] = load_ncdata(thetaforcingfile,'temperature'); % read data from nc file
	zq=-zc'; % the z query points (m)
	zq(zq>max(nc_z)) = max(nc_z); % nearest neighbor extrapolation for the upper layers
	thetainit = interp3(nc_x,nc_y,nc_z',nc_data,xc,yc,zq); % interpolate theta onto cell centers (deg C)
	write_binfile(thetainitfile,permute(thetainit,[2,1,3])); % write to binary input file with size Nx, Ny, Nz 
	% saltinitfile
	[nc_data, nc_x, nc_y, nc_z] = load_ncdata(saltforcingfile,'salinity'); % read data from nc file
	zq=-zc; % the z query points (m)
	zq(zq>max(nc_z)) = max(nc_z); % nearest neighbor extrapolation for the upper layers
	saltinit = interp3(nc_x,nc_y,nc_z',nc_data,xc,yc,zq'); % interpolate theta onto cell centers (deg C)
	write_binfile(saltinitfile,permute(saltinit,[2,1,3])); % write to binary input file
	% some checks {{{
	if any(isnan(thetainit(:)))
		error('nan value detected in thetainit, potential problem with interpolation');
	end
	if any(isnan(saltinit(:)))
		error('nan value detected in saltinit, potential problem with interpolation');
	end
	% }}}
	% }}}
	% etainitfile
	% }}}
	% input/data.OBCS binary files {{{

	% Coordinates on the Western (left side) boundary: these are all size (Nz,Ny)
	% U points (on the Arakawa C grid) are on xg, yc
	% V points (on the Arakawa C grid) are on xc, yg
	% Center points are on xc, yc
	OBW_xc = repmat(xc(OB_I),    Nz, Ny); % x cell center for OBW
	OBW_xg = repmat(xp(OB_I),    Nz, Ny); % x cell edge for OBW
	OBW_yc = repmat(yc,          Nz, 1 ); % y cell center for OBW
	OBW_yg = repmat(yp(1:end-1), Nz, 1 ); % y cell edge for OBW
	OBW_zc = repmat(-zc',         1,  Ny); % z cell center for OBW

	% Coordinates on the Southern (bottom side) boundary: these are all size (Nz,Nx)
	% U points (on the Arakawa C grid) are on xg, yc
	% V points (on the Arakawa C grid) are on xc, yg
	% Center points are on xc, yc
	OBS_xc = repmat(xc,          Nz, 1 ); % x cell center for OBS
	OBS_xg = repmat(xp(1:end-1), Nz, 1 ); % x cell edge for OBS
	OBS_yc = repmat(yc(OB_J),    Nz, Nx); % y cell center for OBS
	OBS_yg = repmat(yp(OB_J),    Nz, Nx); % y cell edge for OBS
	OBC_zc = repmat(-zc',         1,  Nx); % z cell center for OBS

	warning('Using placeholder zero values for OBCS vel data');
	uvelOBW = zeros(Nz,Ny); % x component of left boundary velocity
	vvelOBW = zeros(Nz,Ny); % y component of left boundary velocity

	uvelOBS = zeros(Nz,Nx); % % x component of bottom boundary velocity
	vvelOBS = zeros(Nz,Nx); % % y component of bottom boundary velocity

	% uvelOBWfile vvelOBWfile thetaOBWfile saltOBWfile
	% uvelOBSfile vvelOBSfile thetaOBSfile saltOBSfile 
	% }}}







end % }}}

if perform(org,'PlotMITgcmDomain'), % Plot the MITgcm domain {{{
	x = ncread('/totten_1/ModelData/ISMIP6/Projections/AIS/Ocean_Forcing/climatology_from_obs_1995-2017/obs_salinity_1995-2017_8km_x_60m.nc','x');
	y = ncread('/totten_1/ModelData/ISMIP6/Projections/AIS/Ocean_Forcing/climatology_from_obs_1995-2017/obs_salinity_1995-2017_8km_x_60m.nc','y');
	S = ncread('/totten_1/ModelData/ISMIP6/Projections/AIS/Ocean_Forcing/climatology_from_obs_1995-2017/obs_salinity_1995-2017_8km_x_60m.nc','salinity');

	EXP.x=[-1669315.4719493026,-1669315.4719493026,-1193987.0047960179,...
		-1026979.7055259449,-1026979.7055259449,-1556906.7128252152,...
		-1772089.1945770399,-1772089.1945770399,-1669315.4719493026];
	% outline of ice domain in y (m)
	EXP.y=[-420940.0927398402,-829553.2314259715,-829553.2314259715,...
		-530867.1000391102,-58750.3117179424,170008.8123696489,...
		70446.7685740285,-420940.0927398402,-420940.0927398402];

	% outline of Naughten et al., 2023
	lon1 = -140;
	lon2 = -80;
	lon = [lon1 lon1:lon2 lon2];
	lat1 = 62;
	lat = [90, repmat(lat1,1,length(lon)-2), 90];
	[x_kn,y_kn] = ll2xy(lat,lon,-1);

	% all of antarctica 
	AntEXP = expread('/totten_1/ModelData/Antarctica/Exp/Antarctica.exp');

	figure(1); clf; hold on;
	%imagesc(x,y,S(:,:,1))
	%bathymetry = interpBedmachineAntarctica(XC,YC,'bed','linear'); % interpolate bed elevation data from BedMachine (m)
	%imagesc(XC,YC,reshape(bathymetry,Ny,Nx));
	plot(AntEXP.x,AntEXP.y)
	set(gca,'YDir','normal')
	axis equal;

	[xLim yLim] = ll2xy([60 60 60 60],[270 90 180 0],-1);
	xlim(xLim(1:2));
	ylim(yLim(3:4));
	%plot(EXP.x,EXP.y);

	plot([min(xp) min(xp) max(xp) max(xp) min(xp)], [min(yp) max(yp) max(yp) min(yp) min(yp)],'k');

	plot(x_kn,y_kn);

	U_coordOBW = [xp(OB_I).*ones(1,Ny); yc         ]; % (1,:) x coordinates for U on the C Arakawa grid, (2,:) y coordinates for U on the C Arakawa grid
	V_coordOBW = [xc(OB_I).*ones(1,Ny); yp(1:end-1)]; % (1,:) x coordinates for V on the C Arakawa grid, (2,:) y coordinates for V on the C Arakawa grid

	U_coordOBS = [xp(1:end-1); yc(OB_J).*ones(1,Nx)]; % (1,:) x coordinates for U on the C Arakawa grid, (2,:) y coordinates for U on the C Arakawa grid
	V_coordOBS = [xc;          yp(OB_J).*ones(1,Nx)]; % (1,:) x coordinates for V on the C Arakawa grid, (2,:) y coordinates for V on the C Arakawa grid
end % }}}

% Move back to root directory
disp(['Moving to root directory: ',proph_dir]);
cd(proph_dir);

% local functions 
function mit = mitgcmstruct() % {{{
	% initialize model parameter structure
	mit = struct();
	mit.mesh = struct();
	% initialize parameters for climate forcing model from Naughten et al.,2023
	%mit.climateforcings.mesh = struct();
end	% }}}
function A = loadinitfield(fname,fieldname) % {{{
	A = getfield(load(fname,fieldname),fieldname); % load state variable
	A = pagetranspose(reshape(S,size(S,2),size(S,1),size(S,3))); % correct the indexing to column-major
end	% }}}
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
