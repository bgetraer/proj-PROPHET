steps = [5];

analysisdir = '/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/analysis_Dan';
datadir = fullfile(analysisdir,'data');
expdir = '/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/exp';
uncoupleddir = '/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/ISMIP_style_melt/Models';
org=organizer('repository',analysisdir,'prefix','PROPHET_analysis_','steps',steps);

if perform(org,'basin_partitions') % {{{
	% LOAD SHAPE FILES AND INTERPOLATION MESH FROM DAN
	filename = 'IceBoundaries_Y2014-2016_Antarctica/IceBoundaries_Y2014-2016_Antarctica.shp';
	A = shaperead(fullfile(analysisdir, 'dan_analysis_files', filename));
	is_grounded = strcmp({A.TYPE},"GR"); % index of grounded ice domains
	is_floating = strcmp({A.TYPE},"FL"); % index of floating ice domains

	filename = 'for_averaging.mat';
	B = load(fullfile(analysisdir, 'dan_analysis_files', filename));
	[B.X, B.Y] = meshgrid(B.x_mesh,B.y_mesh); % mesh grid of the x and y cell boundaries
	B.XC = B.X(1:end-1,1:end-1) + diff(B.X(1:end-1,:),[],2)/2; % mesh grid of the x cell centers
	B.YC = B.Y(1:end-1,1:end-1) + diff(B.Y(:,1:end-1),[],1)/2; % mesh grid of the y cell centers
	B.mask = B.mask'==1;
	B.Areas = diff(B.X(1:end-1,:),[],2) .* diff(B.Y(:,1:end-1),[],1);

	% PIG — combination of grounded Pine_Island and floating Pine_Island
	name = "Pine_Island";
	PIG_GR = A(strcmp({A.NAME},name) & is_grounded);
	PIG_FL = A(strcmp({A.NAME},name) & is_floating);

	PIG_polyvec_gr = [polyshape(PIG_GR.X,PIG_GR.Y)]; % region grounded polyvec
	PIG_polyvec_fl = [polyshape(PIG_FL.X,PIG_FL.Y)]; % region floating polyvec
	[PIG_mask, PIG_shape] = mask_region(PIG_polyvec_gr,PIG_polyvec_fl,B.XC,B.YC,B.mask); % extract region shape and mask

	% THW — combination of grounded Thwaites, floating Thwaites, and grounded Haynes
	name = "Thwaites";
	THWAITES_GR = A(strcmp({A.NAME},name) & is_grounded);
	THWAITES_FL = A(strcmp({A.NAME},name) & is_floating);
	name = "Haynes";
	HAYNES_GR = A(strcmp({A.NAME},name) & is_grounded);

	THW_polyvec_gr = [polyshape(THWAITES_GR.X,THWAITES_GR.Y),...
		polyshape(HAYNES_GR.X,HAYNES_GR.Y)]; % region grounded polyvec
	THW_polyvec_fl = [polyshape(THWAITES_FL.X,THWAITES_FL.Y)]; % region floating polyvec

	[THW_mask, THW_shape] = mask_region(THW_polyvec_gr,THW_polyvec_fl,B.XC,B.YC,B.mask); % extract region shape and mask

	% SMITH - combination of grounded Kohler, grounded Pope, grounded Smith, floating Dotson, and floating Crosson.
	name = "Kohler";
	KOHLER_GR = A(strcmp({A.NAME},name) & is_grounded);
	name = "Pope";
	POPE_GR = A(strcmp({A.NAME},name) & is_grounded);
	name = "Smith";
	SMITH_GR = A(strcmp({A.NAME},name) & is_grounded);
	name = "Dotson";
	DOTSON_FL = A(strcmp({A.NAME},name) & is_floating);
	name = "Crosson";
	CROSSON_FL = A(strcmp({A.NAME},name) & is_floating);

	SMITH_polyvec_gr = [polyshape(KOHLER_GR.X,KOHLER_GR.Y),...
		polyshape(POPE_GR.X,POPE_GR.Y),...
		polyshape(SMITH_GR.X,SMITH_GR.Y)]; % region grounded polyvec
	SMITH_polyvec_fl = [polyshape(DOTSON_FL.X,DOTSON_FL.Y),...
		polyshape(CROSSON_FL.X,CROSSON_FL.Y)];  % region floating polyvec
	[SMITH_mask, SMITH_shape] = mask_region(SMITH_polyvec_gr,SMITH_polyvec_fl,B.XC,B.YC,B.mask); % extract region shape and mask

	% ALL - all grounded and floating areas
	ALL_mask = (PIG_mask | THW_mask | SMITH_mask);
	ALL_shape = union([PIG_shape, THW_shape, SMITH_shape]);

	% GET INTERPOLATION INDICES AND POINTS
	X_q = B.XC(ALL_mask);
	Y_q = B.YC(ALL_mask);
	Areas_q = B.Areas(ALL_mask);

	PIG_ind = find(PIG_mask(ALL_mask));
	THW_ind = find(THW_mask(ALL_mask));
	SMITH_ind = find(SMITH_mask(ALL_mask));

	% SAVE
	fname = 'interpPoints.mat';
	save(fname, 'X_q', 'Y_q', 'Areas_q', 'PIG_ind', 'THW_ind', 'SMITH_ind', ...
		'PIG_shape', 'THW_shape', 'SMITH_shape', 'ALL_shape');

	% PLOT
	figure(1);clf;hold on;
	% plot(PIG_shape,'FaceColor','b','EdgeColor','none')
	plot(X_q(PIG_ind),Y_q(PIG_ind),'.b')
	contour(B.XC,B.YC,PIG_mask,[0.5,0.5],'EdgeColor','k','LineWidth',2)
	plot(PIG_shape,'FaceColor','none','EdgeColor','k','LineWidth',2)
	% plot(THW_shape,'FaceColor','y','EdgeColor','none')
	plot(X_q(THW_ind),Y_q(THW_ind),'.g')
	contour(B.XC,B.YC,THW_mask,[0.5,0.5],'EdgeColor','k','LineWidth',2)
	plot(THW_shape,'FaceColor','none','EdgeColor','k','LineWidth',2)
	% plot(SMITH_shape,'FaceColor','r','EdgeColor','none')
	plot(X_q(SMITH_ind),Y_q(SMITH_ind),'.r')
	contour(B.XC,B.YC,SMITH_mask,[0.5,0.5],'EdgeColor','k','LineWidth',2)
	plot(SMITH_shape,'FaceColor','none','EdgeColor','k','LineWidth',2)
	expdisp(fullfile(expdir,'domain.exp'));

	axis equal tight
	title('Subregions of the model domain')
	subtitle('Specified by Dan Goldberg')
end % }}}
if perform(org,'assemble_coupled_issm_results') % {{{
	% load initial model
	md = loadmodel('/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/ISSM_initialization/Models/PROPHET_issm_init_MeshParam.mat');
	% load interpolation points
	disp('Loading interpolation points');
	fname = 'interpPoints.mat';
	D = load(fname);
	% calculate flotation thickness
	H_fl = max(0, -(md.materials.rho_water/md.materials.rho_ice)* md.geometry.bed); % flotation thickness on ISSM mesh (m)
	H_fl_q = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, H_fl, D.X_q, D.Y_q); % flotation thickness at cell centers (m)

	% parallel
	nproc = 56;
	if ~exist(gcp("nocreate"))
		parpool(nproc);
	end

	% loop over each experiment
	experiments = ["Paris2C","RCP85"];
	for i = 1:numel(experiments)
		fprintf('Assembling ISSM results for experiment %s\n',experiments(i));
		resultsdir = fullfile('/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/',experiments(i),'/RUN02/results/');
		files = {dir(resultsdir).name};
		ind = find(contains(files,'issmDiag'));
		fprintf('Found %i files\n',numel(ind));
		% loop through all files
		clear results;
		%for (j = 1:3)%numel(ind))
		parfor (j = 1:numel(ind), nproc)
		fprintf('   file %i/%i...\n',j,numel(ind));
		fname = fullfile(resultsdir,files{ind(j)});
		% load model in parallel
		md = loadmodel('/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/ISSM_initialization/Models/PROPHET_issm_init_MeshParam.mat');
		% load results into model
		md.results.TransientSolution = loadmodel(fname);
		results(j).time = md.results.TransientSolution.time/3600/24/360 + 2010; % model time (seconds)
		results(j).H = md.results.TransientSolution.Thickness; % ice thickness (m)
		results(j).Melt = md.results.TransientSolution.BasalforcingsFloatingiceMeltingRate; % Melt Rate (m/year)
		results(j).masksub500 = 1-2*(md.results.TransientSolution.Base<-500); % index of ice base less than 500m depth (<0 is within mask)
	end
	%% save concatenated results (reading and writing may not be faster than just loading again in parallel...)
	%fname = sprintf('%s_ISSMresults',experiments(i));
	%fprintf('   saving results to %s\n',fullfile(datadir,fname));
	%save(fullfile(datadir,fname),'results','-v7.3');

	% interp thickness from ISSM mesh to interpolation points
	Hmatrix = [results(:).H]; % make one nelement x ntimesteps matrix (m)
	H_q = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x ,md.mesh.y , Hmatrix, D.X_q, D.Y_q); % thickness at cell centers (m)
	Meltmatrix = [results(:).Melt]; % make one nelement x ntimesteps matrix (m/year)
	Melt_q = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x ,md.mesh.y, Meltmatrix,  D.X_q, D.Y_q); % melt rate at cell centers (m/year)
	masksub500matrix = [results(:).masksub500]; % make one nelement x ntimesteps matrix (mask)
	masksub500_q = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x ,md.mesh.y, masksub500matrix, D.X_q, D.Y_q); % cell centers with base <500m (mask)
	masksub500_q = masksub500_q<0; % cell centers with base <500m (logical)
	% calculate thickness above flotation
	HAF_q = max((H_q-H_fl_q), 0); % height above flotation (m)
	% calculate volume above flotation, simple entire cell "in or out" mask 
	VAF_q = HAF_q.*D.Areas_q; % volume above flotation at cell centers (m^3)
	% calculate melt volume per year, simple entire cell "in or out" mask
	Melt_vol_q = Melt_q.*D.Areas_q; % melt volume per year (m^3/year)

	% calculate Dan's requested fields:
	%	1. annual (or per-timestep) VAF for each of PIG, THW and Smith, in km^3. (MAF in Gt.)
	%	2. annual (or per-coupled-timestep) total melt in Gt/a for PIG, THW and Smith.
	%	3. annual (or per-coupled-timestep) total melt in Gt/a for PIG, THW and Smith below 500m.
	time = [results.time]; % calendar years
	VAF.PIG		= 1E-9 * sum(VAF_q(D.PIG_ind,:),1); % volume above flotation for Pine Island Glacier catchment (km^3)
	VAF.THW		= 1E-9 * sum(VAF_q(D.THW_ind,:),1);   % volume above flotation for Thwaites Glacier cathcment (km^3)
	VAF.SMITH	= 1E-9 * sum(VAF_q(D.SMITH_ind,:),1); % volume above flotation for Smith Glacier cathment (km^3)
	VAF.units	= 'km^3';
	MELT.PIG		= 1E-9 * sum(Melt_vol_q(D.PIG_ind,:),1);   % total melt from Pine Island Glacier catchment (Gt/year)
	MELT.THW		= 1E-9 * sum(Melt_vol_q(D.THW_ind,:),1);   % total melt from Thwaites Glacier catchment (Gt/year)
	MELT.SMITH	= 1E-9 * sum(Melt_vol_q(D.SMITH_ind,:),1); % total melt from Smith Glacier catchment (Gt/year)
	MELT.units	= 'Gt/year';
	MELT_sub500m.PIG	= 1E-9 * sum(masksub500_q(D.PIG_ind,:) .* Melt_vol_q(D.PIG_ind,:),1);   % total sub-500m melt from Pine Island catchment (Gt/year)
	MELT_sub500m.THW	= 1E-9 * sum(masksub500_q(D.THW_ind,:) .* Melt_vol_q(D.THW_ind,:),1);   % total sub-500m melt from Thwaites catchment (Gt/year)
	MELT_sub500m.SMITH= 1E-9 * sum(masksub500_q(D.SMITH_ind,:) .* Melt_vol_q(D.SMITH_ind,:),1); % total sub-500m melt from Smith catchment (Gt/year)
	MELT_sub500m.units= 'Gt/year';

	% save results
	fname = sprintf('%s_results',experiments(i));
	fprintf('   saving results to %s\n',fullfile(datadir,fname));
	save(fullfile(datadir,fname),'time','VAF','MELT','MELT_sub500m','-v7.3');
end	
delete(gcp("nocreate"));
end % }}}
if perform(org,'assemble_uncoupled_issm_results') % {{{
	% load initial model
	md = loadmodel('/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/ISSM_initialization/Models/PROPHET_issm_init_MeshParam.mat');
	% basin flags
	disp('Flagging elements in sub-basins');
	flag_basin21=FlagElements(md,fullfile(expdir,'reg21_thwaites.exp')); % flags for basin 21
	flag_basin22=FlagElements(md,fullfile(expdir,'reg22_pineisland.exp')); % flags for basin 22

	% loop over each experiment
	experiments = ["Paris2C","RCP85"];
	% initialize time vector
	years = 2013:2100; % all available years (years)
	duration = 88; % duration to load (years)
	time=years(1) + ([1:duration*12]-1)./12; % the time vector with (years)

	% run each decade
	decade_endyear = years(find(mod(years,10)==0));

	for i = 1:numel(experiments)
		fprintf('Assembling ISSM results for experiment %s\n',experiments(i));
		ind = 0; % this is the index of the time step we are entering
		for d=1:length(decade_endyear)
			% time
			d_time = time(time>=(decade_endyear(d)-10) &  time<=decade_endyear(d));
			disp(sprintf('Processing years: %4.0i--%4.0i',floor(d_time(1)),floor(d_time(end))));

			% load results
			fname = sprintf('PROPHET_ISMIPstyle_%s_results_%4.0i-%4.0i',experiments(i),floor(d_time(1)),floor(d_time(end)));
			disp(['  loading ' fname])
			md.results.TransientSolution = loadmodel(fullfile(uncoupleddir,fname));

			% save vaf results into structure
			for j = 1:numel(md.results.TransientSolution)
				ind = ind + 1;
				results(ind).time = md.results.TransientSolution(j).time; % model time (seconds)
				results(ind).ice_vaf = md.results.TransientSolution(j).IceVolumeAboveFloatation; % total domain volume above flotation (m^3)
				% calculate VAF for sub-basins
				results(ind).ice_vaf_basin21 = VolumeAboveFloatation(md,j,flag_basin21); % basin volume above flotation at TransientSolution(1) (m^3)
				results(ind).ice_vaf_basin22 = VolumeAboveFloatation(md,j,flag_basin22); % basin volume above flotation at TransientSolution(1) (m^3)
				% calculate melt for sub-basins
				meltrate = mean(md.results.TransientSolution(j).BasalforcingsFloatingiceMeltingRate(md.mesh.elements),2);
				results(ind).melt = sum(meltrate .* GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y)); % total melt (m^3 / yr)
				results(ind).melt_basin21 = sum(meltrate .* GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y) .* flag_basin21); % basin melt m^3/yr
				results(ind).melt_basin22 = sum(meltrate .* GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y) .* flag_basin22); % basin melt m^3/yr
			end
		end
		% save concatenated results
		fname = sprintf('%s_uncoupled_ISSMresults',experiments(i));
		fprintf('   saving results to %s\n',fullfile(datadir,fname));
		save(fullfile(datadir,fname),'results','-v7.3');
		disp('   done');
	end
end % }}}
if perform(org,'assemble_control_issm_results') % {{{
	% load results
	disp('Loading ISSM results for control experiment.');
	filename='/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/Depth_dependent_melt/Models/PROPHET_DepthDep_Results.mat';
	md = loadmodel(filename);
	% load interpolation points
	disp('Loading interpolation points');
	fname = 'interpPoints.mat';
	D = load(fname);
	% calculate flotation thickness
	H_fl = max(0, -(md.materials.rho_water/md.materials.rho_ice)* md.geometry.bed); % flotation thickness on ISSM mesh (m)
	H_fl_q = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, H_fl, D.X_q, D.Y_q); % flotation thickness at cell centers (m)

	fprintf('Assembling ISSM results for control experiment\n');

	% interp thickness from ISSM mesh to interpolation points
	Hmatrix = [md.results.TransientSolution.Thickness]; % make one nelement x ntimesteps matrix (m)
	H_q = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x ,md.mesh.y , Hmatrix, D.X_q, D.Y_q); % thickness at cell centers (m)
	Meltmatrix = [md.results.TransientSolution.BasalforcingsFloatingiceMeltingRate]; % make one nelement x ntimesteps matrix (m/year)
	Melt_q = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x ,md.mesh.y, Meltmatrix,  D.X_q, D.Y_q); % melt rate at cell centers (m/year)
	masksub500matrix = 1-2*([md.results.TransientSolution.Base]<-500); % index of ice base less than 500m depth (<0 is within mask)];
	masksub500_q = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x ,md.mesh.y, masksub500matrix, D.X_q, D.Y_q); % cell centers with base <500m (mask)
	masksub500_q = masksub500_q<0; % cell centers with base <500m (logical)
	% calculate thickness above flotation
	HAF_q = max((H_q-H_fl_q), 0); % height above flotation (m)
	% calculate volume above flotation, simple entire cell "in or out" mask 
	VAF_q = HAF_q.*D.Areas_q; % volume above flotation at cell centers (m^3)
	% calculate melt volume per year, simple entire cell "in or out" mask
	Melt_vol_q = Melt_q.*D.Areas_q; % melt volume per year (m^3/year)

	% calculate Dan's requested fields:
	%	1. annual (or per-timestep) VAF for each of PIG, THW and Smith, in km^3. (MAF in Gt.)
	%	2. annual (or per-coupled-timestep) total melt in Gt/a for PIG, THW and Smith.
	%	3. annual (or per-coupled-timestep) total melt in Gt/a for PIG, THW and Smith below 500m.
	time = [md.results.TransientSolution.time]; % calendar years
	time = time-time(1)+2013; % set nominal start time to 2013  
	VAF.PIG		= 1E-9 * sum(VAF_q(D.PIG_ind,:),1); % volume above flotation for Pine Island Glacier catchment (km^3)
	VAF.THW		= 1E-9 * sum(VAF_q(D.THW_ind,:),1);   % volume above flotation for Thwaites Glacier cathcment (km^3)
	VAF.SMITH	= 1E-9 * sum(VAF_q(D.SMITH_ind,:),1); % volume above flotation for Smith Glacier cathment (km^3)
	VAF.units	= 'km^3';
	MELT.PIG		= 1E-9 * sum(Melt_vol_q(D.PIG_ind,:),1);   % total melt from Pine Island Glacier catchment (Gt/year)
	MELT.THW		= 1E-9 * sum(Melt_vol_q(D.THW_ind,:),1);   % total melt from Thwaites Glacier catchment (Gt/year)
	MELT.SMITH	= 1E-9 * sum(Melt_vol_q(D.SMITH_ind,:),1); % total melt from Smith Glacier catchment (Gt/year)
	MELT.units	= 'Gt/year';
	MELT_sub500m.PIG	= 1E-9 * sum(masksub500_q(D.PIG_ind,:) .* Melt_vol_q(D.PIG_ind,:),1);   % total sub-500m melt from Pine Island catchment (Gt/year)
	MELT_sub500m.THW	= 1E-9 * sum(masksub500_q(D.THW_ind,:) .* Melt_vol_q(D.THW_ind,:),1);   % total sub-500m melt from Thwaites catchment (Gt/year)
	MELT_sub500m.SMITH= 1E-9 * sum(masksub500_q(D.SMITH_ind,:) .* Melt_vol_q(D.SMITH_ind,:),1); % total sub-500m melt from Smith catchment (Gt/year)
	MELT_sub500m.units= 'Gt/year';

	% save results
	fname = 'Control_results';
	fprintf('   saving results to %s\n',fullfile(datadir,fname));
	save(fullfile(datadir,fname),'time','VAF','MELT','MELT_sub500m','-v7.3');
end % }}}
if perform(org,'plot_issm_results') % {{{
	% constants 
	rhoGt = 917E-3; % density of ice in Gt/km^3
	gt2mmslr = 1/361.8;  % 361.8 Gt of ice will raise global sea levels by ~1 mm

	% load results
	disp('Loading results...')
	experiments = ["Paris2C","RCP85"];
	fname = sprintf('%s_results',experiments(1));
	Paris2C_results = load(fullfile(datadir,fname));
	fname = sprintf('%s_results',experiments(2));
	RCP85_results = load(fullfile(datadir,fname));
	fname = 'Control_results';
	Control_results = load(fullfile(datadir,fname));

%	fname = sprintf('%s_uncoupled_ISSMresults',experiments(1));
%	Paris2C_UC_results = loadmodel(fullfile(datadir,fname));
%	fname = sprintf('%s_uncoupled_ISSMresults',experiments(2));
%	RCP85_UC_results = loadmodel(fullfile(datadir,fname));

	% plot
	figure(1);clf;
	%subplot(2,1,1);
	% VAF data
	yyaxis('left');hold on;
	% RCP85
	t = RCP85_results.time;
	ice_maf_THW = (RCP85_results.VAF.THW - RCP85_results.VAF.THW(1)).*rhoGt;
	hRCP85_THW = plot(t,ice_maf_THW,'-.r');

	% Paris2C
	t = Paris2C_results.time;
	ice_maf_THW = (Paris2C_results.VAF.THW - Paris2C_results.VAF.THW(1)).*rhoGt;
	hParis2C_THW = plot(t,ice_maf_THW,'--b');

	% CONTROL
	t = Control_results.time;
   ice_maf_THW = (Control_results.VAF.THW - Control_results.VAF.THW(1)).*rhoGt;
   hControl_THW = plot(t,ice_maf_THW,'-k');

%	% Uncoupled RCP85
%	t = [RCP85_UC_results.time];
%	ice_maf_basin21 = ([RCP85_UC_results.ice_vaf_basin21] - RCP85_UC_results(1).ice_vaf_basin21).*rhoGt;
%	ice_maf_basin22 = ([RCP85_UC_results.ice_vaf_basin22] - RCP85_UC_results(1).ice_vaf_basin22).*rhoGt;
%	ice_maf_total = ice_maf_basin21+ice_maf_basin22;
%	%hRCP85_UC_total = plot(t,ice_maf_total,'-r','linewidth',2);
%	hRCP85_UC_basin21 = plot(t,ice_maf_basin21,':r','linewidth',1);
%	hRCP85_UC_basin22 = plot(t,ice_maf_basin22,':r','linewidth',1);
%
%	% Uncoupled Paris2C
%	t = [Paris2C_UC_results.time];
%	ice_maf_basin21 = ([Paris2C_UC_results.ice_vaf_basin21] - Paris2C_UC_results(1).ice_vaf_basin21).*rhoGt;
%	ice_maf_basin22 = ([Paris2C_UC_results.ice_vaf_basin22] - Paris2C_UC_results(1).ice_vaf_basin22).*rhoGt;
%	ice_maf_total = ice_maf_basin21+ice_maf_basin22;
%	%hParis2C_UC_total = plot(t,ice_maf_total,'-r','linewidth',2);
%	hParis2C_UC_basin21 = plot(t,ice_maf_basin21,':b','linewidth',1);
%	hParis2C_UC_basin22 = plot(t,ice_maf_basin22,':b','linewidth',1);


	set([hRCP85_THW,hParis2C_THW,hControl_THW],'linewidth',1.5);
	legend([hParis2C_THW,hRCP85_THW,hControl_THW],'Paris 2C','RCP 8.5','CONTROL','location','sw')
	ylabel('\Delta mass above flotation (Gt)');
	ylimleft=ylim;
	yyaxis('right');hold on;
	ylimright=ylimleft.*gt2mmslr;
	ylim(ylimright);
	ylabel('sea level rise equivalence (mm)');
	xlim([2013,2100])
	set(gca,'fontsize',14)


	return;
	% FIGURE 2: MELT
	subplot(2,1,2); cla;hold on;
	% Paris2C
	t = [Paris2C_results.time];
	%hParis2C   = plot(t,[Paris2C_results.melt]*1000*1E-12,'b');
	hParis2C = plot(t,[Paris2C_results.melt_basin21]*1000*1E-12,'b');
	plot(t,[Paris2C_results.melt_basin22]*1000*1E-12,'b');
	% RCP85
	t = [RCP85_results.time];
	%hRCP85 = plot(t,[RCP85_results.melt]*1000*1E-12,'r');
	hRCP85 = plot(t,[RCP85_results.melt_basin21]*1000*1E-12,'r');
	plot(t,[RCP85_results.melt_basin22]*1000*1E-12,'r');
	%set([hRCP85,hParis2C],'linewidth',1);

%	% Uncoupled Paris2C
%	t = [Paris2C_UC_results.time];
%	%hParis2C   = plot(t,[Paris2C_UC_results.melt]*1000*1E-12,':b');
%	plot(t,[Paris2C_UC_results.melt_basin21]*1000*1E-12,':b');
%	plot(t,[Paris2C_UC_results.melt_basin22]*1000*1E-12,':b');
%	% Uncoupled RCP85
%	t = [RCP85_UC_results.time];
%	%hRCP85 = plot(t,[RCP85_UC_results.melt]*1000*1E-12,':r');
%	plot(t,[RCP85_UC_results.melt_basin21]*1000*1E-12,':r');
%	plot(t,[RCP85_UC_results.melt_basin22]*1000*1E-12,':r');

	ylabel('total melt (Gt/yr)');

	legend([hParis2C,hRCP85],'Paris 2C','RCP 8.5','location','nw')
	xlabel('time (y)')

	xlim([2013,2100])
	set(gca,'fontsize',14)
end % }}}

if perform(org,'plot_thermocline') % {{{
	disp('loading mit structure')
	mit=loadmodel('/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/MITgcm_initialization/Models/PROPHET_mitgcm_init_CompileMITgcm.mat');

	Paris2C=load(fullfile(datadir,'thermoclineParis2C.mat'));
	RCP85=load(fullfile(datadir,'thermoclineRCP85.mat'));
	% Paris2C
	t=Paris2C.t./3600/24/360+2010;
	ymin = min([mit.mesh.zp(find(isnan(Paris2C.T.obw(:,1)),1)), ...
		mit.mesh.zp(find(isnan(Paris2C.T.thwaites(:,1)),1)), ...
		mit.mesh.zp(find(isnan(RCP85.T.obw(:,1)),1)), ...
		mit.mesh.zp(find(isnan(RCP85.T.thwaites(:,1)),1))]);

	figure(1);clf;
	subplot(2,2,1);hold on;
	h=pcolor(t,mit.mesh.zc,Paris2C.T.obw);
	set(h,'edgecolor','flat');
	set(gca,'ydir','normal');axis tight;
	set(gca,'XTickLabel',[],'layer', 'top')
	ylim([ymin,0])
	subplot(2,2,2);hold on;
	h=pcolor(t,mit.mesh.zc,Paris2C.T.thwaites);
	set(h,'edgecolor','flat');
	set(gca,'ydir','normal');axis tight;
	set(gca,'XTickLabel',[],'YTickLabel',[],'layer', 'top')
	ylim([ymin,0])

	% RCP85
	t=RCP85.t./3600/24/360+2010;
	subplot(2,2,3);hold on;
	h=pcolor(t,mit.mesh.zc,RCP85.T.obw);
	set(h,'edgecolor','flat');
	set(gca,'ydir','normal');axis tight;
	set(gca,'XTick',2020:10:2100,'layer', 'top')
	ylim([ymin,0])
	subplot(2,2,4);hold on;
	h=pcolor(t,mit.mesh.zc,RCP85.T.thwaites);
	set(h,'edgecolor','flat');
	set(gca,'ydir','normal');axis tight;
	set(gca,'XTick',2020:10:2100,'YTickLabel',[],'layer', 'top')
	ylim([ymin,0])

	%colorbar('location','north')
	colormap(flip(brewermap(100,'RdYlBu')))
end % }}}

% FUNCTIONS
function [region_mask, region_shape] = mask_region(polyvec_gr,polyvec_fl,Xq,Yq,mask) % {{{
	region_shape_gr = union(polyvec_gr); % create union of the polyshape vector
	region_mask_gr = inpolygon(Xq,Yq,region_shape_gr.Vertices(:,1),region_shape_gr.Vertices(:,2)) & mask; % mask the mesh center points for grounded nunataks

	region_shape_fl = union(polyvec_fl); % create union of the polyshape vector
	region_shape_fl = rmholes(region_shape_fl); % remove holes from polyshape
	region_mask_fl = inpolygon(Xq,Yq,region_shape_fl.Vertices(:,1),region_shape_fl.Vertices(:,2)); % do not mask the mesh center points for floating ice

	region_shape = union([region_shape_gr,region_shape_fl]); % create union of the polyshape vector
	region_mask = (region_mask_gr | region_mask_fl); % mask of points in the shape
end % }}}
