steps = 4;

analysisdir = '/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/analysis';
datadir = '/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/data';
expdir = '/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/exp';
uncoupleddir = '/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/ISMIP_style_melt/Models';
org=organizer('repository',analysisdir,'prefix','PROPHET_analysis_','steps',steps);

if perform(org,'basin_partitions') % {{{
	% https://earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems
	% Data table with columns "Lat", "Lon", "Basin_ID" 
	fname = 'Antarctic_Drainage_System/ant_full_drainagesystem_polygons.txt';
	A = readmatrix(fullfile(expdir,fname),'NumHeaderLines',7); 

	% Thwaites Glacier: id 21; Pine Island Glacier: id 22
	basin_id = [21,22];
	basin_name = ["thwaites","pineisland"];
	fmtstr = 'reg%02i_%s.exp';

	for i = 1:numel(basin_id)
		fname = sprintf(fmtstr,basin_id(i),basin_name(i));
		[EXP.x,EXP.y] = ll2xy(A(A(:,3)==basin_id(i),1),A(A(:,3)==basin_id(i),2),-1); % Polar stereographic xy coordinates
		expwrite(EXP,fullfile(expdir,fname)); % write exp file
	end

	% load example data
	md = loadmodel('../experiments/ISSM_initialization/Models/PROPHET_issm_init_MeshParam.mat');

	% plot the regions
	plotmodel(md,'data',md.initialization.vel)
	expdisp(fullfile(expdir,'domain.exp'));
	expdisp(fullfile(expdir,'reg21_thwaites.exp'),'linestyle','y');
	expdisp(fullfile(expdir,'reg21_pineisland.exp'),'linestyle','g');
	hold on
	axis equal tight
	title('Subregions of the model domain')
	subtitle('Antarctic basins 21 and 22')
end % }}}
if perform(org,'assemble_coupled_issm_results') % {{{
	% load initial model
	md = loadmodel('/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/ISSM_initialization/Models/PROPHET_issm_init_MeshParam.mat');
	% basin flags
	disp('Flagging elements in sub-basins');
	flag_basin21=FlagElements(md,fullfile(expdir,'reg21_thwaites.exp')); % flags for basin 21
	flag_basin22=FlagElements(md,fullfile(expdir,'reg22_pineisland.exp')); % flags for basin 22

	% parallel
	nproc = 50;
	parpool(nproc);

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
		parfor (j = 1:numel(ind), nproc)
			fprintf('   file %i/%i...\n ',j,numel(ind));
			fname = fullfile(resultsdir,files{ind(j)});
			% load model in parallel
			md = loadmodel('/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/ISSM_initialization/Models/PROPHET_issm_init_MeshParam.mat');
			% load results into model
			md.results.TransientSolution = loadmodel(fname);
			results(j).time = md.results.TransientSolution.time/3600/24/360 + 2010;; % model time (seconds)
			results(j).ice_vaf = md.results.TransientSolution.IceVolumeAboveFloatation; % total domain volume above flotation (m^3)
			% calculate VAF for sub-basins
			results(j).ice_vaf_basin21 = VolumeAboveFloatation(md,1,flag_basin21); % basin volume above flotation at TransientSolution(1) (m^3)
			results(j).ice_vaf_basin22 = VolumeAboveFloatation(md,1,flag_basin22); % basin volume above flotation at TransientSolution(1) (m^3)
			% calculate melt for sub-basins
			meltrate = mean(md.results.TransientSolution.BasalforcingsFloatingiceMeltingRate(md.mesh.elements),2);
			results(j).melt =	sum(meltrate .* GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y)); % total melt (m^3 / yr)
			results(j).melt_basin21 = sum(meltrate .* GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y) .* flag_basin21); % basin melt m^3/yr
			results(j).melt_basin22 = sum(meltrate .* GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y) .* flag_basin22); % basin melt m^3/yr
		end
		% save concatenated results
		fname = sprintf('%s_ISSMresults',experiments(i));
		fprintf('   saving results to %s\n',fullfile(datadir,fname));
		save(fullfile(datadir,fname),'results','-v7.3');
		disp('   done');
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
if perform(org,'plot_issm_results') % {{{
	% constants 
	rhoGt = 917*1e-12; % density of ice in Gt/m^3
	gt2mmslr = 1/361.8;  % 361.8 Gt of ice will raise global sea levels by ~1 mm

	% load results
	experiments = ["Paris2C","RCP85"];
	fname = sprintf('%s_ISSMresults',experiments(1));
	Paris2C_results = loadmodel(fullfile(datadir,fname));
	fname = sprintf('%s_ISSMresults',experiments(2));
	RCP85_results = loadmodel(fullfile(datadir,fname));

	fname = sprintf('%s_uncoupled_ISSMresults',experiments(1));
	Paris2C_UC_results = loadmodel(fullfile(datadir,fname));
	fname = sprintf('%s_uncoupled_ISSMresults',experiments(2));
	RCP85_UC_results = loadmodel(fullfile(datadir,fname));

	% plot
	figure(1);clf;
	subplot(2,1,1);
	% VAF data
	yyaxis('left');hold on;
	% RCP85
	t = [RCP85_results.time];
	ice_maf_basin21 = ([RCP85_results.ice_vaf_basin21] - RCP85_results(1).ice_vaf_basin21).*rhoGt;
	ice_maf_basin22 = ([RCP85_results.ice_vaf_basin22] - RCP85_results(1).ice_vaf_basin22).*rhoGt;
	ice_maf_total = ice_maf_basin21+ice_maf_basin22;
	%hRCP85_total = plot(t,ice_maf_total,'-r','linewidth',2);
	hRCP85_basin21 = plot(t,ice_maf_basin21,'-r');
	hRCP85_basin22 = plot(t,ice_maf_basin22,'-r');

	% Paris2C
	t = [Paris2C_results.time];
	ice_maf_basin21 = ([Paris2C_results.ice_vaf_basin21] - Paris2C_results(1).ice_vaf_basin21).*rhoGt;
	ice_maf_basin22 = ([Paris2C_results.ice_vaf_basin22] - Paris2C_results(1).ice_vaf_basin22).*rhoGt;
	ice_maf_total = ice_maf_basin21+ice_maf_basin22;
	%hParis2C_total = plot(t,ice_maf_total,'-b','linewidth',2);
	hParis2C_basin21 = plot(t,ice_maf_basin21,'-b');
	hParis2C_basin22 = plot(t,ice_maf_basin22,'-b');

	% Uncoupled RCP85
	t = [RCP85_UC_results.time];
	ice_maf_basin21 = ([RCP85_UC_results.ice_vaf_basin21] - RCP85_UC_results(1).ice_vaf_basin21).*rhoGt;
	ice_maf_basin22 = ([RCP85_UC_results.ice_vaf_basin22] - RCP85_UC_results(1).ice_vaf_basin22).*rhoGt;
	ice_maf_total = ice_maf_basin21+ice_maf_basin22;
	%hRCP85_UC_total = plot(t,ice_maf_total,'-r','linewidth',2);
	hRCP85_UC_basin21 = plot(t,ice_maf_basin21,':r','linewidth',1);
	hRCP85_UC_basin22 = plot(t,ice_maf_basin22,':r','linewidth',1);

	% Uncoupled Paris2C
	t = [Paris2C_UC_results.time];
	ice_maf_basin21 = ([Paris2C_UC_results.ice_vaf_basin21] - Paris2C_UC_results(1).ice_vaf_basin21).*rhoGt;
	ice_maf_basin22 = ([Paris2C_UC_results.ice_vaf_basin22] - Paris2C_UC_results(1).ice_vaf_basin22).*rhoGt;
	ice_maf_total = ice_maf_basin21+ice_maf_basin22;
	%hParis2C_UC_total = plot(t,ice_maf_total,'-r','linewidth',2);
	hParis2C_UC_basin21 = plot(t,ice_maf_basin21,':b','linewidth',1);
	hParis2C_UC_basin22 = plot(t,ice_maf_basin22,':b','linewidth',1);


	%set([hRCP85_total,hParis2C_total],'linewidth',3);
	set([hRCP85_basin21,hRCP85_basin22,hParis2C_basin21,hParis2C_basin22],'linewidth',1.5);
	legend([hParis2C_basin21,hRCP85_basin21],'Paris 2C','RCP 8.5','location','sw')
	ylabel('\Delta mass above flotation (Gt)');
	ylimleft=ylim;
	yyaxis('right');hold on;
	ylimright=ylimleft.*gt2mmslr;
	ylim(ylimright);
	ylabel('sea level rise equivalence (mm)');
	xlim([2013,2100])
	set(gca,'fontsize',14)


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

	% Uncoupled Paris2C
   t = [Paris2C_UC_results.time];
   %hParis2C   = plot(t,[Paris2C_UC_results.melt]*1000*1E-12,':b');
   plot(t,[Paris2C_UC_results.melt_basin21]*1000*1E-12,':b');
   plot(t,[Paris2C_UC_results.melt_basin22]*1000*1E-12,':b');
   % Uncoupled RCP85
   t = [RCP85_UC_results.time];
   %hRCP85 = plot(t,[RCP85_UC_results.melt]*1000*1E-12,':r');
   plot(t,[RCP85_UC_results.melt_basin21]*1000*1E-12,':r');
   plot(t,[RCP85_UC_results.melt_basin22]*1000*1E-12,':r');

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
