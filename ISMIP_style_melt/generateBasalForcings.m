steps=[2];
prophdir='/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET';
exp_name='Paris2C';
prefix=['PROPHET_ISMIPstyle_' exp_name '_'];
modeldir='./Models';
	
org=organizer('repository',modeldir,'prefix',prefix,'steps',steps);
if perform(org,'BasalForcings') % {{{
	% load mit structure
	mit=loadmodel(fullfile(prophdir,'experiments/MITgcm_initialization/Models/PROPHET_mitgcm_init_CompileMITgcm.mat')); % mit structure
	md=loadmodel(fullfile(prophdir,'experiments/ISSM_initialization/Models/PROPHET_issm_init_InversionC.mat')); % md structure

	processforcing=0;
	forcingfile=['./Models/' exp_name 'forcingdata.mat'];
	if (~isfile(forcingfile) | processforcing)
		% Initialize the time vector 
		years=2010:2100; % all available years (years)
		duration=3; % duration to load (years)
		nfiles = ceil(duration); % the number of files we need to load
		time=years(1) + ([1:duration*12]-1)./12; % the time vector with (years)

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
	bedmachinepath='/totten_1/ModelData/Antarctica/BedMachine/BedMachineAntarctica-v4.0.nc'; % path to dataset
	BED = interpBedmachineGreenland(Xq,Yq,'bed','linear',bedmachinepath);
	M   = interpBedmachineGreenland(Xq,Yq,'mask','nearest',bedmachinepath);

	disp('   -- Define open ocean pixels');
	[~,indx0] =min(abs(xq-mit.mesh.xc(1))); % find the nearest x in the ISMIP6 Grid to the OBW
	[~,indy0] =min(abs(yq-mit.mesh.yc(1))); % find the nearest y in the ISMIP6 Grid to the OBS
	[~,indxend] = min(abs(xq-max(X))); % find the nearest x in the ISMIP6 Grid to boundary end
	[~,indyend] = min(abs(yq-max(Y))); % find the nearest y in the ISMIP6 Grid to boundary end
	FAR = (Xq>=xq(indx0) & Xq<=xq(indxend) & Yq==yq(indy0)) | (Yq>=yq(indy0) & Yq<=yq(indyend) & Xq==xq(indx0));

	% define tf_depths, the levels that we will interpolate onto
	tf_depths = -(10:20:950)';

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

	disp('   -- Get temperature and salinity maps for each tf_depths');
	% visulatization of min depth indices for each level of tf_depths
	tf = cell(1,1,numel(tf_depths));
	for i=1:numel(tf_depths)
		disp(['   -- Depth = ' num2str(tf_depths(i)) ' ' num2str(i) '/' num2str(numel(tf_depths))]);

		% min depth indexing
		levelID=min(ID,i);
		% initialize the state fields for this level	
		temperature  =single(NaN([numel(yq),numel(xq)]));
		salinity     =single(NaN([numel(yq),numel(xq)]));
		% initialize the md.basalforcings.tf field time-series for this level
		tf{1,1,i}=zeros((md.mesh.numberofvertices+1),numel(time));

		[posi posj]=find(levelID~=0);
		ind = sub2ind(size(levelID),posi,posj);
		for t=1:numel(time)
			temperature(ind) = Tq(levelID(ind),t); % in-situ temp at min depth (deg C)
			salinity(ind)    = Sq(levelID(ind),t); % in-situ salt at min depth
			freezingpoint = a.*salinity + b + c.*tf_depths(i); % deg C
			theta = temperature - freezingpoint; % thermal forcing (deg C)
			tf{1,1,i}(1:end-1,t) = InterpFromGridToMesh(xq,yq,theta,md.mesh.x,md.mesh.y,0); % tf at md vert. for this time (deg C)
			tf{1,1,i}(end,t)     = time(t); % the time stamp for this column (y)
		end
		%	% plot level ID
		%	figure(1);clf;
		%	imagesc(xq,yq,levelID)
		%	set(gca,'ydir','normal');axis equal tight
		%	colorbar; caxis([0,length(tf_depths)]);
		%	title(['depth level: ' num2str(i) '/' num2str(numel(tf_depths))])
		%	drawnow;
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

	md.basalforcings=basalforcings;

	savemodel(org,md);
end % }}}
if perform(org,'TransientRun') % {{{
	md=loadmodel(org,'BasalForcings');

	%Find elements that don't have any ice
	pos = find(min(md.mask.ice_levelset(md.mesh.elements),[],2)<0);
	md.mask.ice_levelset(md.mesh.elements(pos,:))= -1;
	md.mask.ice_levelset   = reinitializelevelset(md, md.mask.ice_levelset);

	% solve stress

	%Set parameters
	md.inversion.iscontrol=0;
	md.transient.ismovingfront=0;
	md.transient.isthermal=0;
	md.transient.isstressbalance=1;
	md.transient.ismasstransport=1;
	md.transient.isgroundingline=1;
	md.groundingline.migration = 'SubelementMigration';
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','BasalforcingsFloatingiceMeltingRate','Thickness','MaskOceanLevelset'};
	md.settings.output_frequency = 1;
	md.timestepping=timesteppingadaptive();
	md.timestepping.time_step_max=0.05;
	md.timestepping.time_step_min=0.0005;

	% miscelleneous
	md.miscellaneous.name=[prefix 'run'];

	% timestepping
	md.timestepping.start_time=2010;
	md.timestepping.final_time=2013;

	% cluster
	md.cluster=generic('name',oshostname(),'np',75);

	% output
	md.verbose.solution=true;
	md.verbose.solver  =false;

	% solve
	md=solve(md,'tr');

	% save
	results=md.results;
	save([prefix 'results'],'results')
	savemodel(org,md);
end % }}}
