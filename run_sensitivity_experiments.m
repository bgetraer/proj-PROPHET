steps=[6];
experiment_numbers=[5,6];

for exp_no = experiment_numbers
% PROPHET Amundsen Sea Coupling sensitivity testing
% Each sensitivity experiment lives in its own directory with its own input directory.
% Some files are taken from experiments/MITgcm_initialization/ and were created by the runme.m file controlling
% the initializaiton of the PROPHET model. These experiments assume the existence of the ocean model, ice model,
% and corresponding files setup in the other runme.m files.
% Outline:
%	Set boundary forcings
%  Initialize ocean state
%  Run MITgcm
% 
% Assumes the existence of the ocean model, ice model, and corresponding files setup in the the other runme.m files

% directory structure {{{
mitgcm_dir='/nobackup/bgetraer/MITgcm'; % MITgcm directory (pleaides)
proph_dir ='/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET'; % base directory for this project
exp_root_dir=fullfile(proph_dir,'experiments/sensitivity_experiments/'); % base directory for all of the sensitivity results
init_dir = fullfile(proph_dir,'experiments/MITgcm_initialization/');
% }}}
% experiments {{{
% table index for all experiments
varNames = ["Name", "is_coupled", "forcing_period", "Description"];
varTypes = ["cellstr", "logical", "cellstr", "cellstr"];
experiments = table('Size',[1,numel(varNames)], 'VariableTypes', varTypes, 'VariableNames', varNames);

experiments(1, varNames) = {'No_melt',					false,	'',			'Ice only model, no coupling, no melt.'};
experiments(2, varNames) = {'KN_constant_clim',		true,		'constant',	'Constant forcing fields, 2010--2020 average from KN Paris2C model results.'};
experiments(3, varNames) = {'KN_monthly_clim',		true,		'monthly',	'Monthly forcing fields, 2010--2020 monthly averages from KN Paris2C model results.'};
experiments(4, varNames) = {'PW700_constant',		true,		'constant',	'Constant forcing fields, piecewise, thermocline bot. 700m.'};
experiments(5, varNames) = {'PW600_constant',		true,		'constant',	'Constant forcing fields, piecewise, thermocline bot. 600m.'};
experiments(6, varNames) = {'PW800_constant',		true,		'constant',	'Constant forcing fields, piecewise, thermocline bot. 800m.'};
experiments(7, varNames) = {'PW700_amp50_per2',		true,		'periodic',	'Periodic forcing fields, piecewise, thermocline bot. 700m, 50m amplitude, 2 yr period.'};
experiments(8, varNames) = {'PW700_amp50_per5',		true,		'periodic',	'Periodic forcing fields, piecewise, thermocline bot. 700m, 50m amplitude, 5 yr period.'};
experiments(9, varNames) = {'PW700_amp50_per10',	true,		'periodic',	'Periodic forcing fields, piecewise, thermocline bot. 700m, 50m amplitude, 10 yr period.'};
experiments(10, varNames) = {'PW700_amp100_per2',	true,		'periodic',	'Periodic forcing fields, piecewise, thermocline bot. 700m, 100m amplitude, 2 yr period.'};
experiments(11, varNames) = {'PW700_amp100_per5',	true,		'periodic',	'Periodic forcing fields, piecewise, thermocline bot. 700m, 100m amplitude, 5 yr period.'};
experiments(12, varNames) = {'PW700_amp100_per10',	true,		'periodic',	'Periodic forcing fields, piecewise, thermocline bot. 700m, 100m amplitude, 10 yr period.'};

% setup experiment directories
expdirname = cell(size(experiments,1),1);
for i=1:numel(expdirname)
	expdirname{i} = sprintf('se%03i_%s',i,experiments.Name{i});
	% directory for this experiments
	exp_dir=fullfile(exp_root_dir,expdirname{i});
	if ~exist(exp_dir)
		mkdir(exp_dir);
	end
	input_dir = fullfile(exp_dir,'input');
	if ~exist(input_dir)
		mkdir(input_dir);
	end
	model_dir = fullfile(exp_dir,'Models');
	if ~exist(model_dir)
		mkdir(model_dir);
	end
end
clear exp_dir input_dir model_dir;
experiments.dirname = expdirname;
% }}}

% choose the experiment to run
experiment = table2struct(experiments(exp_no,:));
exp_dir = fullfile(exp_root_dir,experiment.dirname);

% suppress warnings {{{
wid = {'MATLAB:hg:AutoSoftwareOpenGL','MATLAB:polyshape:repairedBySimplify'};
for i=1:length(wid)
	warning('off',wid{i});
end
% }}}
org=organizer('repository',fullfile(exp_dir,'Models'),'prefix',experiment.dirname,'steps',steps);
disp_experiment(experiment);

if perform(org,'BoundaryAndInitialConditions'), % {{{
	mit=loadmodel(fullfile(init_dir,'Models/PROPHET_mitgcm_init_CompileMITgcm.mat'));

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

	% default to turning piecewise off, must be turned on for each case

	% Define the OBCS forcing fields and theta/salt initialization profiles
	switch char(experiment.dirname)
		case 'se001_No_melt' % {{{
			disp('No ocean forcing to initialize.');
			% }}}
		case 'se002_KN_constant_clim' % {{{
			% Boundary Forcings
			% re-process bounary forcings from Kaitlyn Naughten (these files have already been averaged across
			% the 10 member ensemble, and interpolated onto the correct grid for my boundary).
			years = 2010:2019; % years to average over
			n_records = 12; % number of records per file (1 per month)
			D = cell(numel(mit.forcing.field_bdry_out),1);
			% loop over each boundary field
			for i = 1:numel(mit.forcing.field_bdry_out)
				fprintf('Processing %s\n',mit.forcing.field_bdry_out{i})
				% files names
				fname_base = sprintf('%s.%s_',mit.forcing.field_bdry_out{i},'Paris2C');
				filenames = append(fname_base,string(years));

				% initialize matrix
				if contains(mit.forcing.field_bdry_out{i},'obw'); Nh = mit.mesh.Ny; else; Nh = mit.mesh.Nx; end
				A = zeros(Nh,mit.mesh.Nz,n_records,numel(years));

				for j = 1:numel(years)
					A(:,:,:,j) = binread(fullfile(mit.forcing.Ddir,filenames(1)),8,Nh,mit.mesh.Nz,n_records);
				end

				% extract desired climatology
				D{i} = mean(A,[3,4]);	% average over all months and years
			end

			% Initial conditions:
			% extract all theta and salt data
			theta_ref = cell2mat(D(contains(mit.forcing.field_bdry_out,'theta')));
			salt_ref = cell2mat(D(contains(mit.forcing.field_bdry_out,'salt')));
			% take average over horizontal dimension, ignoring missing data (0)
			dims = [1]; 
			tRef = sum(theta_ref, dims)./sum(theta_ref~=0,dims);
			sRef = sum(salt_ref, dims)./sum(salt_ref~=0,dims);
			% replace missing data with linear interpolation or nearest extrapolation
			tRef = fillmissing(tRef,'linear','EndValues','nearest');
			sRef = fillmissing(sRef,'linear','EndValues','nearest');
			% }}}
		case 'se003_KN_monthly_clim' % {{{
			% Boundary Forcings
			% re-process bounary forcings from Kaitlyn Naughten (these files have already been averaged across
			% the 10 member ensemble, and interpolated onto the correct grid for my boundary).
			years = 2010:2019; % years to average over
			n_records = 12; % number of records per file (1 per month)
			D = cell(numel(mit.forcing.field_bdry_out),1);
			% loop over each boundary field
			for i = 1:numel(mit.forcing.field_bdry_out)
				fprintf('Processing %s\n',mit.forcing.field_bdry_out{i})
				% files names
				fname_base = sprintf('%s.%s_',mit.forcing.field_bdry_out{i},'Paris2C');
				filenames = append(fname_base,string(years));

				% initialize matrix
				if contains(mit.forcing.field_bdry_out{i},'obw'); Nh = mit.mesh.Ny; else; Nh = mit.mesh.Nx; end
				A = zeros(Nh,mit.mesh.Nz,n_records,numel(years));

				for j = 1:numel(years)
					A(:,:,:,j) = binread(fullfile(mit.forcing.Ddir,filenames(1)),8,Nh,mit.mesh.Nz,n_records);
				end

				% extract desired climatology
				D{i} = nanmean(A,4);
			end

			% Initial conditions:
			% extract all theta and salt data
			theta_ref = cell2mat(D(contains(mit.forcing.field_bdry_out,'theta')));
			salt_ref = cell2mat(D(contains(mit.forcing.field_bdry_out,'salt')));
			% take average over horizontal dimension and months, ignoring missing data (0)
			dims = [1,3]; 
			tRef = sum(theta_ref, dims)./sum(theta_ref~=0,dims);
			sRef = sum(salt_ref, dims)./sum(salt_ref~=0,dims);
			% replace missing data with linear interpolation or nearest extrapolation
			tRef = fillmissing(tRef,'linear','EndValues','nearest');
			sRef = fillmissing(sRef,'linear','EndValues','nearest');
			% }}}
		otherwise % assume piecewise forcing, extract parameters from name {{{
			use_piecewise_forcing = contains(experiment.Name,'PW');
			assert(use_piecewise_forcing,'Non-piecewise experiments must be defined separately!');

			% extract parameters from name
			depth_thermocline = str2double(regexp(experiment.Name,'(?<=PW)\d+','match','once'));  % 'PWX_'  depth of bottom of thermocline given by X (m)
			amplitude			= str2double(regexp(experiment.Name,'(?<=amp)\d+','match','once')); % 'ampX_' amplitude of periodic thermocline depth given by X (m)
			period				= str2double(regexp(experiment.Name,'(?<=per)\d+','match','once')); % 'perX'  period of periodic forcing given by X (years)

			% additional flags
			use_periodic_forcing = strcmp(experiment.forcing_period,'periodic');

			% Define baseline ocean conditions for synthetic boundary conditions (after De Rydt et al., 2014 and Bett et al., 2024)
			z_top       = -300; % baseline z-coordinate of the top of the thermocline (winter water layer ends) (m)
			z_bot       = -700; % baseline z-coordinate of the bottom of the thermocline (circumpolar deep water layer begins) (m)
			theta_ww    = -1.0; % potential temperature of winter water (deg C)
			theta_cdw   = +1.2; % potential temperature of circumpolar deep water (deg C)
			salt_ww		= 34;	  % salinity of winter water (g/kg)
			salt_cdw		= 34.7; % salinity of circumpolar deep water (g/kg)
			dtheta_cdw_OBS = -0.6; % potential temperature anomaly at "southern" open boundary (deg C) 
			dsalt_cdw_OBS  = -0.1; % salinity anomaly at "southern" open boundary (g/kg)

			% define dummy z-coordinate to reflect thermocline depth anomaly
			dz_thermocline = (-depth_thermocline) - z_bot; % difference between experiment thermocline bottom depth and baseline depth (m)
			zc_prime = mit.mesh.zc - dz_thermocline; % dummy z-coordinates for calculating isocline (m)

			% reference profiles of potential temperature and salinity
			tRef = get_isocline(zc_prime,z_top,z_bot,theta_ww,theta_cdw);
			sRef = get_isocline(zc_prime,z_top,z_bot,salt_ww,salt_cdw);
			
			if use_periodic_forcing
				% Time is non-dimensionalized by the forcing period T. Frequencies are non-dimensionalized as cycles per period.
				% The actual period is set in MITgcm EXF as the "field period" which corresponds to the actual sampling interval.
				% All forcing files will have the number of records corresponding to fs_norm, repeating over them.
				fs_norm = 24; % normalized sampling frequency (samples per period)
				dt_norm = 1/fs_norm; % normalized sampling interval (fraction of period per sample)
				t_norm  = (dt_norm/2):dt_norm:(1-dt_norm/2); % normalized time grid spanning one period (0,1), with forcing at the middle of the time step
				dz_norm = sin(2*pi*t_norm); % normalized periodic forcing (unit forcing, one period) 

				% periodically changing dummy coordinates for calculating isocline (Nz x Nrecords) 
				zc_sampling = reshape(zc_prime,mit.mesh.Nz,1) - amplitude.*reshape(dz_norm,1,fs_norm); % (m)

				% boundary profiles of potential temperature and salinity (Nz x Nrecords)
				OBWtheta_prof = get_isocline(zc_sampling,z_top,z_bot,theta_ww,theta_cdw);
            OBWsalt_prof  = get_isocline(zc_sampling,z_top,z_bot,salt_ww, salt_cdw);

				OBStheta_prof = get_isocline(zc_sampling,z_top,z_bot,theta_ww,theta_cdw + dtheta_cdw_OBS);
            OBSsalt_prof  = get_isocline(zc_sampling,z_top,z_bot,salt_ww, salt_cdw  + dsalt_cdw_OBS);

				% Western Open Boundary
				OBWtheta = repmat(reshape(OBWtheta_prof,1,mit.mesh.Nz,fs_norm),mit.mesh.Ny,1,1);
				OBWsalt  = repmat(reshape(OBWsalt_prof, 1,mit.mesh.Nz,fs_norm),mit.mesh.Ny,1,1);
				% apply mask
				OBWtheta = OBWtheta.*squeeze(mit.geometry.open_mask(:,1,:));
				OBWsalt  = OBWsalt .*squeeze(mit.geometry.open_mask(:,1,:));

				% Southern Open Boundary
				OBStheta = repmat(reshape(OBStheta_prof,1,mit.mesh.Nz,fs_norm),mit.mesh.Nx,1,1);
				OBSsalt  = repmat(reshape(OBSsalt_prof, 1,mit.mesh.Nz,fs_norm),mit.mesh.Nx,1,1);
				% apply mask
				OBStheta = OBStheta.*squeeze(mit.geometry.open_mask(1,:,:));
				OBSsalt  = OBSsalt .*squeeze(mit.geometry.open_mask(1,:,:));

				% check sizes to confirm
				assert(all(size(OBWtheta)==[mit.mesh.Ny,mit.mesh.Nz,fs_norm]));
				assert(all(size(OBWsalt) ==[mit.mesh.Ny,mit.mesh.Nz,fs_norm]));
				assert(all(size(OBStheta)==[mit.mesh.Nx,mit.mesh.Nz,fs_norm]));
				assert(all(size(OBSsalt) ==[mit.mesh.Nx,mit.mesh.Nz,fs_norm]));
			else
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
				OBWsalt	= OBWsalt .*squeeze(mit.geometry.open_mask(:,1,:));

				% Southern Open Boundary
				OBStheta = repmat(OBStheta_prof,mit.mesh.Nx,1);
				OBSsalt	= repmat(OBSsalt_prof, mit.mesh.Nx,1);
				% apply mask
				OBStheta = OBStheta.*squeeze(mit.geometry.open_mask(1,:,:));
				OBSsalt  = OBSsalt .*squeeze(mit.geometry.open_mask(1,:,:));

				% check sizes to confirm
				assert(all(size(OBWtheta)==[mit.mesh.Ny,mit.mesh.Nz]));
				assert(all(size(OBWsalt) ==[mit.mesh.Ny,mit.mesh.Nz]));
				assert(all(size(OBStheta)==[mit.mesh.Nx,mit.mesh.Nz]));
				assert(all(size(OBSsalt) ==[mit.mesh.Nx,mit.mesh.Nz]));
			end

			% velocity from constant climate forcings
			OBvel_dirname = 'se002_KN_constant_clim';
			OBvel_dir = fullfile(exp_root_dir,OBvel_dirname,'input');
			OBvel_fmtstr = '%s.%s';

			Nrecords = size(OBWtheta,3); % duplicate the field to match the number of records needed

			% W boundary velocity
			fname = fullfile(OBvel_dir,sprintf(OBvel_fmtstr,'uvel.obw',OBvel_dirname));
			OBWuvel = binread(fname,8,mit.mesh.Ny,mit.mesh.Nz);
			OBWuvel = repmat(OBWuvel,1,1,Nrecords);
			fname = fullfile(OBvel_dir,sprintf(OBvel_fmtstr,'vvel.obw',OBvel_dirname));
			OBWvvel = binread(fname,8,mit.mesh.Ny,mit.mesh.Nz);
			OBWvvel = repmat(OBWvvel,1,1,Nrecords);

			% S boundary velocity
			fname = fullfile(OBvel_dir,sprintf(OBvel_fmtstr,'uvel.obs',OBvel_dirname));
			OBSuvel = binread(fname,8,mit.mesh.Nx,mit.mesh.Nz);
			OBSuvel = repmat(OBSuvel,1,1,Nrecords);
			fname = fullfile(OBvel_dir,sprintf(OBvel_fmtstr,'vvel.obs',OBvel_dirname));
			OBSvvel = binread(fname,8,mit.mesh.Nx,mit.mesh.Nz);
			OBSvvel = repmat(OBSvvel,1,1,Nrecords);

			% structure the OB fields
			D = cell(numel(mit.forcing.field_bdry_out),1);
			D{strcmp(mit.forcing.field_bdry_out,'theta.obw')}	= OBWtheta;
			D{strcmp(mit.forcing.field_bdry_out,'salt.obw')}	= OBWsalt;
			D{strcmp(mit.forcing.field_bdry_out,'uvel.obw')}	= OBWuvel;
			D{strcmp(mit.forcing.field_bdry_out,'vvel.obw')}	= OBWvvel;
			D{strcmp(mit.forcing.field_bdry_out,'theta.obs')}	= OBStheta;
			D{strcmp(mit.forcing.field_bdry_out,'salt.obs')}	= OBSsalt;
			D{strcmp(mit.forcing.field_bdry_out,'uvel.obs')}	= OBSuvel;
			D{strcmp(mit.forcing.field_bdry_out,'vvel.obs')}	= OBSvvel;
			% }}}
	end

	% start date
	mit.forcing.startdate = datetime(2010,1,1);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% write input data to files
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp(['writing files to ' exp_dir]);
	fname = ['input/' mit.fname.treffile]; 
	disp(['writing tRef to ' fname]);
	write_binfile(fullfile(exp_dir,fname),tRef);
	fname = ['input/' mit.fname.sreffile]; 
	disp(['writing sRef to ' fname]);
	write_binfile(fullfile(exp_dir,fname),sRef);
	fname = ['input/' mit.fname.draftfile];

	% get number of records for each field
	mit.forcing.Nrecords = nan(size(mit.forcing.field_bdry_out));

	% write forcing files or copy them from provided file path
	for i = 1:numel(mit.forcing.field_bdry_out)
		fname = sprintf('input/%s.%s',mit.forcing.field_bdry_out{i},experiment.dirname);
		if isnumeric(D{i})
			% write boundary forcing files
			mit.forcing.Nrecords(i) = size(D{i},3); % number of records for this field
			disp(['writing ' fname ': nrecords = ' num2str(mit.forcing.Nrecords(i))])
			write_binfile(fullfile(exp_dir,fname),D{i});
		elseif isfile(D{i})
			disp(['copying ' dir(D{i}).name ' to ' fname])
			copyfile(D{i},fullfile(exp_dir,fname));
		else
			error('Forcing field is not numeric or a valid file path. File path is incorrect or has not been created yet!')
		end
	end

	% copy files
	filenames = {'bathy.bin','draft.bin','delr.bin'};
	for i=1:numel(filenames)
		fname = ['input/' filenames{i}];
		disp(['copying ' filenames{i} ' to ' fname ]);
		source_path = fullfile(init_dir,fname); % file location
		destination_path = fullfile(exp_dir,fname); % link location
		copyfile(source_path,destination_path)
	end	

	% save mit
	savedata(org,mit);
end % }}}
if perform(org,'RuntimeOptionsOcean') % {{{
	mit=load(org,'BoundaryAndInitialConditions');

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% set runtime options
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp(' - Setting runtime options');
	% runtime options for all experiments {{{
	% TIME STEPPING
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
	mit.inputdata.CAL{1}.startdate_1=string(cal_startdate,'yyyyMMdd'); % yyyyMMdd of start date
	mit.inputdata.CAL{1}.startDate_2=string(cal_startdate,'HHmmss');   % HHmmss of start date

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
	% }}}
	% input/data.exf External forcing parameters {{{
	disp(['  - setting EXF options for experiment.forcing_period = ' experiment.forcing_period]);
	switch experiment.forcing_period
		case 'yearly'
			mit.inputdata.EXF{5}.useOBCSYearlyFields = '.TRUE.'
			mit.inputdata.EXF{end}.obcsWstartdate1 = string(mit.forcing.startdate,'yyyyMMdd'); % W boundary start year (YYYY), month (MM), day (DD) to determine record number
			mit.inputdata.EXF{end}.obcsWperiod     = -1.0;     % interval between two records: the special value -1 means non-repeating (calendar) monthly records
			mit.inputdata.EXF{end}.obcsSstartdate1 = string(mit.forcing.startdate,'yyyyMMdd'); % S boundary start year (YYYY), month (MM), day (DD) to determine record number
			mit.inputdata.EXF{end}.obcsSperiod     = -1.0;     % interval between two records: the special value -1 means non-repeating (calendar) monthly records
		case 'monthly'
			mit.inputdata.EXF{5}.useOBCSYearlyFields = '.FALSE.'
			mit.inputdata.EXF{end}.obcsWstartdate1 = string(mit.forcing.startdate,'yyyyMMdd'); % W boundary start year (YYYY), month (MM), day (DD) to determine record number
			mit.inputdata.EXF{end}.obcsWperiod     = -12.0;     % interval between two records: the special value -12 means 12 repeating (calendar) monthly records
			mit.inputdata.EXF{end}.obcsSstartdate1 = string(mit.forcing.startdate,'yyyyMMdd'); % S boundary start year (YYYY), month (MM), day (DD) to determine record number
			mit.inputdata.EXF{end}.obcsSperiod     = -12.0;     % interval between two records: the special value -12 means 12 repeating (calendar) monthly records
		case 'constant'
			mit.inputdata.EXF{5}.useOBCSYearlyFields = '.FALSE.'
			mit.inputdata.EXF{end}.obcsWstartdate1 = string(mit.forcing.startdate,'yyyyMMdd'); % W boundary start year (YYYY), month (MM), day (DD) to determine record number
			mit.inputdata.EXF{end}.obcsWperiod     = 0.0;     % one file, one record provided
			mit.inputdata.EXF{end}.obcsSstartdate1 = string(mit.forcing.startdate,'yyyyMMdd'); % S boundary start year (YYYY), month (MM), day (DD) to determine record number
			mit.inputdata.EXF{end}.obcsSperiod     = 0.0;     % one file, one record provided
		case 'periodic'
			period = str2double(regexp(experiment.Name,'(?<=per)\d+','match','once')); % 'perX'  period of periodic forcing given by X (years)
			fs_norm = unique(mit.forcing.Nrecords(~isnan(mit.forcing.Nrecords))); % normalized sampling frequency (samples per period)
			assert(numel(fs_norm)==1, 'non-unique sampling frequency!');
			dt_sampling = (1/fs_norm) .* period .* mit.timestepping.y2s;  % actual sampling interval/interval between records (s)
			rep_cycle   = period .* mit.timestepping.y2s;                 % actual time period over which to repeat records (s)

			mit.inputdata.EXF{5}.useOBCSYearlyFields = '.FALSE.'
			mit.inputdata.EXF{end}.obcsWstartdate1 = string(mit.forcing.startdate,'yyyyMMdd'); % W boundary start year (YYYY), month (MM), day (DD) to determine record number
			mit.inputdata.EXF{end}.obcsWperiod     = dt_sampling;     % interval between two records: >0 means cycle through repeating records
			mit.inputdata.EXF{end}.obcsWrepCycle   = rep_cycle;
			mit.inputdata.EXF{end}.obcsSstartdate1 = string(mit.forcing.startdate,'yyyyMMdd'); % S boundary start year (YYYY), month (MM), day (DD) to determine record number
			mit.inputdata.EXF{end}.obcsSperiod     = dt_sampling;     % interval between two records: >0 means cycle through repeating records
			mit.inputdata.EXF{end}.obcsSrepCycle   = rep_cycle;
		otherwise
			error('not set up');
	end % }}}
	% input/data.obcs Open boundary conditions parameters {{{
	% Bottom boundary
	mit.inputdata.OBCS{1}.OBSuFile=['''' mit.fname.uvelOBSfile  experiment.dirname '''']; % Nx by Nz matrix of u velocity at Southern OB
	mit.inputdata.OBCS{1}.OBSvFile=['''' mit.fname.vvelOBSfile  experiment.dirname '''']; % Nx by Nz matrix of v velocity at Southern OB
	mit.inputdata.OBCS{1}.OBStFile=['''' mit.fname.thetaOBSfile experiment.dirname '''']; % Nx by Nz matrix of pot. temp. at Southern OB
	mit.inputdata.OBCS{1}.OBSsFile=['''' mit.fname.saltOBSfile  experiment.dirname '''']; % Nx by Nz matrix of salin. at Southern OB
	% Left boundary
	mit.inputdata.OBCS{1}.OBWuFile=['''' mit.fname.uvelOBWfile  experiment.dirname '''']; % Ny by Nz matrix of u velocity at Western OB
	mit.inputdata.OBCS{1}.OBWvFile=['''' mit.fname.vvelOBWfile  experiment.dirname '''']; % Ny by Nz matrix of v velocity at Western OB
	mit.inputdata.OBCS{1}.OBWtFile=['''' mit.fname.thetaOBWfile experiment.dirname '''']; % Ny by Nz matrix of pot. temp. at Western OB
	mit.inputdata.OBCS{1}.OBWsFile=['''' mit.fname.saltOBWfile  experiment.dirname '''']; % Ny by Nz matrix of salin. at Western OB
	% }}}

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% write all of the input data files
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp(' - Writing runtime options to input/');
	write_datafile(fullfile(exp_dir,mit.fname.eedatafile),       mit.inputdata.EEP,      'EXECUTION ENVIRONMENT PARAMETERS');
	write_datafile(fullfile(exp_dir,mit.fname.datafile),         mit.inputdata.PARM,     'MODEL PARAMETERS');
	write_datafile(fullfile(exp_dir,mit.fname.datapkgfile),      mit.inputdata.PKG,      'PACKAGES');
	write_datafile(fullfile(exp_dir,mit.fname.datashelficefile), mit.inputdata.SHELFICE, 'SHELFICE RUNTIME PARAMETERS');
	write_datafile(fullfile(exp_dir,mit.fname.datacalfile),      mit.inputdata.CAL,      'CALENDAR PARAMETERS');
	write_datafile(fullfile(exp_dir,mit.fname.dataexffile),      mit.inputdata.EXF,      'EXTERNAL FORCINGS PARAMETERS');
	write_datafile(fullfile(exp_dir,mit.fname.datadiagfile),		 mit.inputdata.DIAG,		 'DIAGNOSTICS RUNTIME PARAMETERS');
	write_datafile(fullfile(exp_dir,mit.fname.dataobcsfile),     mit.inputdata.OBCS,		 'OBCS RUNTIME PARAMETERS');

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Diverge run directories for each experiment
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	prompt = 'Reset runocean directory now? (''y'' or ''Y'' to proceed, ''n'' or ''N'' to skip)\n';
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
		disp(' - Preparing runocean directory');
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Directory management
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% rename previous run directory and create new one
		dirname='runocean';
		rundir=fullfile(exp_dir,dirname);
		oldrundir=fullfile(exp_dir,[dirname '.old']);
		if exist(oldrundir)
			system(['\rm -rf ' oldrundir]);
		end
		if exist(rundir)
			system(['\mv ' rundir ' ' oldrundir]);
		end
		% make the run directory in exp_dir
		mkdir(rundir);

		disp(['    linking files to ' rundir])

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Link to run directory
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% make links to input files
		S = dir(fullfile(exp_dir,'input/*'));
		for j=1:length(S)
			if ~S(j).isdir
				file_path = fullfile(S(j).folder, S(j).name); % file location
				command = ['ln -s ' file_path ' ' rundir];
				system(command);
			end
		end
		% make link to mitgcmuv executable
		file_path = fullfile(init_dir,'build/mitgcmuv'); % file location
		link_path = rundir; % link location
		command = ['ln -s ' file_path ' ' rundir];
		system(command);
	else
		disp('Skipping reset of runocean directories!');
	end

	savedata(org,mit);
end % }}}
if perform(org,'RunOcean') % {{{
	mit=load(org,'RuntimeOptionsOcean');
	% check if this experiment needs its own spinup
	switch experiment.forcing_period
		case 'periodic'
			disp('Periodic forcing experiments use ocean spinup from constant piecewise forcing. Skipping RunOcean spinup.');
		otherwise
			% set run parameters for PBS queue file
			rundir = fullfile(exp_dir,'runocean'); % which experiment directory to run
			%rundir = '/nobackupp18/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/test/run';
			grouplist = 's2013'; % account on Pleiades
			npMIT=mit.build.SZ.nPx*mit.build.SZ.nPy; % number of processors for MITgcm
			queuename = 'long'; % which queue to submit to (long or devel)
			walltime = duration(120,0,0); % walltime to request
			% write the .queue file
			fname = write_queuefile(rundir,grouplist,npMIT,'queuename',queuename,'walltime',walltime,'iscoupled',0); % returns the name of the .queue file
			%fname = write_queuefile(rundir,grouplist,1,'HelloWorld'); % returns the name of the .queue file
			fprintf('Submitting queue file:   ')
			command=['qsub ' fullfile(rundir,fname)];
			system(command);
	end
end % }}}
if perform(org,'RuntimeOptionsCoupled') % {{{
	mit=load(org,'RuntimeOptionsOcean');

	% set directories
	builddir      = fullfile(init_dir,'build'); % initalization directory where model was compiled
	inputdir      = fullfile(exp_dir,'input'); % initialization directory for runtime input options
	runcoupleddir = fullfile(exp_dir,'runcoupled'); % run directory
	switch experiment.forcing_period
		case 'periodic'
			 disp('Period forcing experiment: taking initial state from constant spinup.')
			 depth_thermocline = str2double(regexp(experiment.Name,'(?<=PW)\d+','match','once'));  % 'PWX_'  depth of bottom of thermocline given by X (m)
			 spinup_name = sprintf('PW%i_constant',depth_thermocline);
			 spinup_dir = experiments.dirname{contains(experiments.dirname,spinup_name)};
			 runoceandir   = fullfile(exp_root_dir,spinup_dir,'runocean');
		 otherwise
			 runoceandir   = fullfile(exp_dir,'runocean'); % run directory for the ocean model spinup
	 end

	% Shared runtime parameters {{{
	% TIME STEPPING
	% During the coupled phase the ocean model does one run per coupled step:
	% The data parameters are set during the run
	disp(' - Setting timestepping options');
	mit.timestepping.ispickup         = 0;                                 % are we running a pickup or an initial coupling [0,1]
	mit.timestepping.coupled_basetime = mit.timestepping.spinupduration;   % the model time that we are starting coupling from (s)
	mit.timestepping.deltaT_coupled   = 15*24*60*60;                       % coupling time step (s)
	mit.timestepping.nsteps           = 600;                               % number of coupled time steps to take
	mit.timestepping.coupled_endtime  = mit.timestepping.nsteps*mit.timestepping.deltaT_coupled; % the end model time that we are aiming for (s)
	mit.timestepping.startTime        = mit.timestepping.coupled_basetime; % run start time for this integration (s)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% ./data
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% everything gets set here, and during the loop all that is changed is the startTime
	mit.inputdata.PARM{3} = struct();
	% structure information
	mit.inputdata.PARM{3}.header='PARM03';
	mit.inputdata.PARM{3}.description='Time stepping parameters';
	% Run Start and Duration
	mit.inputdata.PARM{3}.nIter0      = 0;                                                            % starting timestep iteration number
	mit.inputdata.PARM{3}.deltaT      = 100;                                                          % mitgcm deltaT (s)
	mit.inputdata.PARM{3}.nEndIter    = mit.timestepping.deltaT_coupled/mit.inputdata.PARM{3}.deltaT; % end timestep iteration number
	mit.inputdata.PARM{3}.startTime   = mit.timestepping.startTime;                                   % run start time for this integration (s)
	% Restart/Pickup Files
	mit.inputdata.PARM{3}.pChkptFreq  = mit.timestepping.deltaT_coupled;                              % permanent pickup checkpoint file write interval (s)
	mit.inputdata.PARM{3}.ChkptFreq   = 0;                                                            % temporary pickup checkpoint file write interval (s)
	% Frequency/Amount of Output
	mit.inputdata.PARM{3}.monitorFreq = mit.timestepping.deltaT_coupled;                              % interval to write monitor output - every coupled time step (s)
	mit.inputdata.PARM{3}.cAdjFreq        = -1;                                                       % frequency of convective adj. scheme
	mit.inputdata.PARM{3}.monitorSelect   = 1;                                                        % group of monitor variables to output
	mit.inputdata.PARM{3}.dumpInitAndLast = '.FALSE.';                                                % write out initial and last iteration model state

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
	filelist={'eedata','data','data.cal','data.diagnostics','data.exf','data.obcs','data.pkg','data.shelfice',...
		'bathy.bin','delr.bin','sref.bin','tref.bin'};
	cp_filelist(inputdir,filelist,runcoupleddir); % copy all files from inputdir to runcoupleddir
	% rename data.obcs, draft file, and make a bathy_ref.bin file
	source_fnames={'draft.bin','bathy.bin'};
	destination_fname={['draft.' suffix '.bin'],'bathy_ref.bin'};
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
	% }}}

	% link the boundary forcing files
	filelist={dir(fullfile(inputdir,['*' experiment.dirname])).name};
	ln_filelist(inputdir,filelist,runcoupleddir);

	save(fullfile(runcoupleddir,'RuntimeOptionsCoupled'),'mit');
end % }}}
if perform(org,'RunCoupled') % {{{
	rundir = fullfile(exp_dir,'runcoupled'); % run directory
	mitfile=fullfile(rundir,'RuntimeOptionsCoupled.mat');
	mit=loadmodel(mitfile);

	% set run parameters for PBS queue file
	mccdir = fullfile(proph_dir,'runcouple/mccfiles/');
	mdfile = fullfile(proph_dir,'experiments/ISSM_initialization/Models/PROPHET_issm_init_TransientPrep.mat');
	grouplist = 's2013'; % account on Pleiades
	npMIT=mit.build.SZ.nPx*mit.build.SZ.nPy; % number of processors for MITgcm
	queuename = 'long'; % which queue to submit to (long or devel)
	%queuename = 'devel'; % which queue to submit to (long or devel)
	walltime = duration(5*24,0,0); % walltime to request

	interactive = 0; % run interactive?
	if interactive
		% INTERACTIVE RUN
		% see pbs scripts folder
		cd(rundir);
		addpath(fullfile(proph_dir,'runcouple'));
		runcouple(mdfile,mitfile);
	else
		% NON-INTERACTIVE RUN
		% write the .queue file
		fname = write_queuefile(rundir,grouplist,npMIT,...
			'queuename',queuename,'walltime',walltime,'iscoupled',1,'mccdir',mccdir,'mccargin',{mdfile,mitfile}); % returns the name of the .queue file
		fprintf('Submitting queue file:   ')
		command=['qsub ' fullfile(rundir,fname)];
		system(command);	
	end

end % }}}

if perform(org,'RunPickup') % {{{
	rundir = fullfile(exp_dir,'runcoupled'); % run directory
	mit=loadmodel(fullfile(rundir,'RuntimeOptionsCoupled'));

	% get the pickup time from the last issmDiag file
	filenames = {dir(fullfile(rundir,'issmDiag*.mat')).name};
   modeltimes = sort(str2double(extractBetween(filenames, 'issmDiag.', '.mat')));
	modeltime_pickup = modeltimes(end);
	fprintf('Experiment: %s Year: %i\n',exp_no,round(modeltime_pickup/mit.timestepping.y2s+2010))	

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
	mit.timestepping.ispickup         = 1;                % are we running a pickup or an initial coupling [0,1]
	mit.timestepping.startTime        = modeltime_pickup; % modeltime of the start of the simulation (s)
	mit.timestepping.nsteps           = (290*360 - modeltime_pickup/3600/24)/15;   % number of coupled time steps to take

	% print settings
	disp('mit.timestepping:');
	disp(mit.timestepping);

	% save mit structure
	mitfile = fullfile(rundir,'RunPickup.mat');
	save(mitfile,'mit');

	% set run parameters for PBS queue file
	mccdir = fullfile(proph_dir,'runcouple/mccfiles/');
	mdfile = fullfile(proph_dir,'experiments/ISSM_initialization/Models/PROPHET_issm_init_TransientPrep.mat');
	grouplist = 's2013'; % account on Pleiades
	npMIT=mit.build.SZ.nPx*mit.build.SZ.nPy; % number of processors for MITgcm
	queuename = 'long'; % which queue to submit to (long or devel)
	%queuename = 'devel'; % which queue to submit to (long or devel)
	walltime = duration(5*24,0,0); % walltime to request


	interactive = 0; % run interactive?
	if interactive
		% INTERACTIVE RUN
		% see pbs scripts folder
		cd(rundir);
		addpath(fullfile(proph_dir,'runcouple'));
		runcouple(mdfile,mitfile);
	else
		% NON-INTERACTIVE RUN
		% write the .queue file
		fname = write_queuefile(rundir,grouplist,npMIT,...
			'queuename',queuename,'walltime',walltime,'iscoupled',1,'mccdir',mccdir,'mccargin',{mdfile,mitfile}); % returns the name of the .queue file
		fprintf('Submitting queue file:   ')
		command=['qsub ' fullfile(rundir,fname)];
		system(command);	
	end
end % }}}

% Move back to root directory
disp(['Moving to root directory: ', proph_dir]);
cd(proph_dir);

end
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
			['ln -sf ' fullfile(mccdir,'run_MCCexecutable.sh') ' ./'], ...
			['ln -sf ' fullfile(mccdir,'MCCexecutable') ' ./'], ...
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
	disp([' - writing OPTIONS file to ' fname]);
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
function A = get_isocline(z,z0,z1,A0,A1) % {{{
	m = [z0,1;z1,1]\[A0;A1];
	A = m(1).*z + m(2);
   A(z>z0)=A0;
   A(z<z1)=A1;
end % }}}
function disp_experiment(experiment) % {{{ 
	msg = ['   experiment: ' experiment.dirname '   '];
	esc = char(27); % ANSI escape
	fprintf('%c[37;41m%s%c[0m\n', esc, msg, esc);  % 37=white fg, 41=red bg, 0=reset
end % }}}
