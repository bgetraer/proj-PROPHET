steps=[2];
%exp_name='Paris2C';
exp_name='RCP85';

prophdir='/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET';
prefix=['PROPHET_ISMIPstyle_' exp_name '_'];
modeldir='./Models';
if ~isdir(modeldir)
	mkdir(modeldir);
end
	
org=organizer('repository',modeldir,'prefix',prefix,'steps',steps);

% Initialize the time vector 
years=2013:2100; % all available years (years)
duration=88; % duration to load (years)
nfiles = ceil(duration); % the number of files we need to load
time=years(1) + ([1:duration*12]-1)./12; % the time vector with (years)

if perform(org,'BasalForcings') % {{{
	% load mit structure
	mit=loadmodel(fullfile(prophdir,'experiments/MITgcm_initialization/Models/PROPHET_mitgcm_init_CompileMITgcm.mat')); % mit structure
	md=loadmodel(fullfile(prophdir,'experiments/ISSM_initialization/Models/PROPHET_issm_init_TransientPrep.mat')); % md structure

	processforcing=0;
	forcingfile=['./Models/' exp_name 'forcingdata.mat'];
	if (~isfile(forcingfile) | processforcing)
		disp('   -- Loading Boundary Conditions');
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

		% in-situ temperature and salinity
		temp_matrix=[]; % initialize temperature matrix
		salt_matrix=[]; % initialize salinity matrix
		for i=1:nfiles % Iterate over years
			% TEMP
			% get obw data in linear index order 
			fname=fullfile(mit.forcing.Ddir,['theta.obw.' exp_name '_' num2str(years(i))]); % filepath
			Tobw = binread(fname,8,mit.mesh.Ny*mit.mesh.Nz,12); % load the temperature records for this year
			% get obs data in linear index order
			fname=fullfile(mit.forcing.Ddir,['theta.obs.' exp_name '_' num2str(years(i))]); % filepath
			Tobs = binread(fname,8,mit.mesh.Nx*mit.mesh.Nz,12); % load the temperature records for this year
			T = [Tobw;Tobs]; % put them together
			rowind=find(sum(T,2)~=0); % find rows with data
			% put in matrix 
			temp_matrix = [temp_matrix T(rowind,:)]; % save to temp matrix

			% SALT
			% get obw data in linear index order
			fname=fullfile(mit.forcing.Ddir,['salt.obw.' exp_name '_' num2str(years(i))]); % filepath
			Sobw = binread(fname,8,mit.mesh.Ny*mit.mesh.Nz,12); % load the temperature records for this year
			% get obs data in linear index order
			fname=fullfile(mit.forcing.Ddir,['salt.obs.' exp_name '_' num2str(years(i))]); % filepath
			Sobs = binread(fname,8,mit.mesh.Nx*mit.mesh.Nz,12); % load the temperature records for this year
			S = [Sobw;Sobs]; % put them together
			rowind=find(sum(S,2)~=0); % find rows with data
			% put in matrix
			salt_matrix = [salt_matrix S(rowind,:)]; % save to salt matrix
		end

		% mask geometry points with no values
		X=X(rowind);
		Y=Y(rowind);
		Z=Z(rowind);
		levels=flip(unique(Z));

		% save
		disp(['Saving forcing data to ' forcingfile]);
		save(forcingfile,'X','Y','Z','levels','temp_matrix','salt_matrix','time');
	else
		disp(['Loading forcing data from ' forcingfile]);
		load(forcingfile);
	end

	disp('   -- Making Large Interpolation Grid');
	xq=min(mit.mesh.xp):1e3:max(mit.mesh.xp);
	yq=min(mit.mesh.yp):1e3:max(mit.mesh.yp);
	[Xq Yq]=meshgrid(xq,yq);

	disp('   -- Interpolating BedMachine bed and mask');
	bedmachinepath='/nobackup/bgetraer/ModelData/BedMachine/BedMachineAntarctica-v4.0.nc'; % path to dataset
	BED = interpBedmachineGreenland(Xq,Yq,'bed','linear',bedmachinepath);
	M   = interpBedmachineGreenland(Xq,Yq,'mask','nearest',bedmachinepath);

	disp('   -- Define open ocean pixels');
	[~,indx0] =min(abs(xq-mit.mesh.xc(1))); % find the nearest x in the ISMIP6 Grid to the OBW
	[~,indy0] =min(abs(yq-mit.mesh.yc(1))); % find the nearest y in the ISMIP6 Grid to the OBS
	[~,indxend] = min(abs(xq-max(X))); % find the nearest x in the ISMIP6 Grid to boundary end
	[~,indyend] = min(abs(yq-max(Y))); % find the nearest y in the ISMIP6 Grid to boundary end
	FAR = (Xq>=xq(indx0) & Xq<=xq(indxend) & Yq==yq(indy0)) | (Yq>=yq(indy0) & Yq<=yq(indyend) & Xq==xq(indx0));

	% define tf_depths, the levels that we will interpolate onto
	tf_depths = -(25:50:950)';

	%Initialize output
	disp('   -- Build connectivity and ID matrices');
	connectivity = false(size(Xq,1),size(Xq,2),numel(tf_depths));
	ID=zeros(size(Xq)); % deepest index
	T=zeros(numel(levels),numel(time)); % level averaged temp, before interpolation
	S=zeros(numel(levels),numel(time)); % level averaged salinity, before interpolation

	for i=1:numel(tf_depths),
		disp(['   -- Depth = ' num2str(tf_depths(i)) ' ' num2str(i) '/' num2str(numel(tf_depths))]);
		disp('    ... Labelling');
		CC=bwlabel(BED<=tf_depths(i));
		disp('    ... Finding unique labels');
		pos=find(FAR & BED<tf_depths(i));
		list = unique(CC(pos));
		connectivity(:,:,i)=ismember(CC,list);
		ID(connectivity(:,:,i))=i;
	end
	for i=1:numel(levels),
		disp(['   -- Depth = ' num2str(levels(i)) ' ' num2str(i) '/' num2str(numel(levels))]);
		disp('    ... Level averaged temp and salt');
		T(i,:)=mean(temp_matrix(find(Z==levels(i)),:),1);
		S(i,:)=mean(salt_matrix(find(Z==levels(i)),:),1);
	end

	% interpolate onto tf_depths
	Tq = interp1(levels,T,tf_depths);
	Sq = interp1(levels,S,tf_depths);

	% constants
   rho_w=1029.35; % approximate value for estimating pressure (kg/m^3)
   % Calculate in situ freezing point (from Holland, Jenkins, and Holland 2008)
   g = 9.81;    % m/s^2
   a = -0.0573; % deg C
   b = 0.0832;  % deg C
   c = 7.53E-3*1E-5*rho_w.*g; % deg C/Pa

	disp('   -- Get temperature and salinity maps for each tf_depths for each decade');
	% divide into decades
	decade_endyear = years(find(mod(years,10)==0));
	for d=1:numel(decade_endyear)
		d_time = time(time>=(decade_endyear(d)-10) &  time<=decade_endyear(d));
		disp(sprintf('   -- DOING %4.0i--%4.0i',floor(d_time(1)),floor(d_time(end))));
		tf = cell(1,1,numel(tf_depths));
		for i=1:numel(tf_depths)
			disp(['   -- Depth = ' num2str(tf_depths(i)) ' ' num2str(i) '/' num2str(numel(tf_depths))]);

			% min depth indexing
			levelID=min(ID,i);
			% initialize the state fields for this level
			temperature  =NaN([numel(yq),numel(xq)]);
			salinity     =NaN([numel(yq),numel(xq)]);
			% initialize the md.basalforcings.tf field time-series for this level
			tf{1,1,i}=single(zeros((md.mesh.numberofvertices+1),numel(d_time)));

			[posi posj]=find(levelID~=0);
			ind = sub2ind(size(levelID),posi,posj);
			for t=1:numel(d_time)
				temperature(ind) = Tq(levelID(ind),t); % in-situ temp at min depth (deg C)
				salinity(ind)    = Sq(levelID(ind),t); % in-situ salt at min depth
				freezingpoint = a.*salinity + b + c.*tf_depths(i); % deg C
				theta = temperature - freezingpoint; % thermal forcing (deg C)
				theta(isnan(theta)) = 0; % replace NaN values for min level zero
				tf{1,1,i}(1:end-1,t) = single(interp2(xq,yq,theta,md.mesh.x,md.mesh.y,'linear',0)); % tf at md vert. for this time (deg C)
				tf{1,1,i}(end,t)     = single(d_time(t)); % the time stamp for this column (y)
			end
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
	
		basalforcings_fname = sprintf('basalforcings_%s_%4.0i-%4.0i.mat',exp_name,floor(d_time(1)),floor(d_time(end)));
		disp(['   -- saving ' basalforcings_fname]);
		save(fullfile(org.repository,basalforcings_fname),'basalforcings');
	end
end % }}}
if perform(org,'TransientRun') % {{{
	% load model
	mdname=fullfile(prophdir,'experiments/ISSM_initialization/Models/PROPHET_issm_init_TransientPrep.mat'); % md structure
	disp(['loading ' dir(mdname).name]);
   md=loadmodel(mdname);

	% set options
	disp('setting transient options');
	md.cluster=generic('name',oshostname(),'np',75);
	md.verbose.solution=true;

	% save
	savemodel(org,md);

	% run each decade
	decade_endyear = years(find(mod(years,10)==0));
   for d=1:length(decade_endyear)
		% time 
		d_time = time(time>=(decade_endyear(d)-10) &  time<=decade_endyear(d));
		disp(sprintf('Begin run: %4.0i--%4.0i',floor(d_time(1)),floor(d_time(end))));
		md.timestepping.start_time=d_time(1);
		md.timestepping.final_time=d_time(end);

		basalforcings_fname = sprintf('basalforcings_%s_%4.0i-%4.0i.mat',exp_name,floor(d_time(1)),floor(d_time(end)));
		disp(['   -- loading ' basalforcings_fname]);
		md.basalforcings=loadmodel(fullfile(org.repository,basalforcings_fname));

		% solve
		md.miscellaneous.name=sprintf('%sTransientRun_%4.0i-%4.0i',prefix,floor(d_time(1)),floor(d_time(end)));
		md=solve(md,'tr');

		% save results every 0.25 years
		fname = sprintf('./Models/%sresults_%4.0i-%4.0i',prefix,floor(d_time(1)),floor(d_time(end)));
		disp(['  saving results to ' fname])
		ind = mod([md.results.TransientSolution.time],0.25)==0;
		results = md.results.TransientSolution(ind);
      save(fname,'results');

      % reinitialize ISSM from results
      disp('  reinitializing ISSM from results.TransientSolution')
      md.geometry.base             = md.results.TransientSolution(end).Base;
      md.geometry.surface          = md.results.TransientSolution(end).Surface;
      md.geometry.thickness        = md.geometry.surface-md.geometry.base;
      md.initialization.vx         = md.results.TransientSolution(end).Vx;
      md.initialization.vy         = md.results.TransientSolution(end).Vy;
      md.initialization.vel        = md.results.TransientSolution(end).Vel;
      md.mask.ocean_levelset       = md.results.TransientSolution(end).MaskOceanLevelset;
      clear md.results;
	end
end % }}}
if perform(org,'VAF') % {{{	
	t=[];
	VAF=[];
	% run each decade
	decade_endyear = years(find(mod(years,10)==0));
   for d=1:5%length(decade_endyear)
		% time 
		d_time = time(time>=(decade_endyear(d)-10) &  time<=decade_endyear(d));
		disp(sprintf('YEARS: %4.0i--%4.0i',floor(d_time(1)),floor(d_time(end))));

		% load results
		fname = sprintf('./Models/%s_results_%4.0i-%4.0i',prefix,floor(d_time(1)),floor(d_time(end)));
		disp(['  loading ' fname])
      load(fname);

		t = [t results.time];
		VAF=[VAF results.IceVolumeAboveFloatation];
	end
	fname = sprintf('./Models/%s_VAF',exp_name);
	save(fname,'t','VAF');
end % }}}
