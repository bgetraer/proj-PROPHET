steps = 1;

analysisdir = '/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/analysis';
datadir = '/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/data';
expdir = '/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/exp';
uncoupleddir = '/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/ISMIP_style_melt/Models';
org=organizer('repository',analysisdir,'prefix','PROPHET_analysis_','steps',steps);

experiment_names = { ...
	'se002_KN_constant_clim', ...
	'se003_KN_monthly_clim', ...
	'se004_PW700_constant', ...
	'se005_PW600_constant', ...
	'se006_PW800_constant', ...
	'se007_PW700_amp50_per2', ...
	'se008_PW700_amp50_per5', ...
	'se009_PW700_amp50_per10', ...
	'se010_PW700_amp100_per2', ...
	'se011_PW700_amp100_per5', ...
	'se012_PW700_amp100_per10', ...
	'Paris2C',...
	'RCP85',...
	};

if perform(org,'plot_issm_results') % {{{
	% constants 
	rhoGt = 917*1e-12; % density of ice in Gt/m^3
	gt2mmslr = 1/361.8;  % 361.8 Gt of ice will raise global sea levels by ~1 mm

	% plot
	figure(1);clf;
	nax_x = 3;
	nax_y = 2;
	nax = nax_x*nax_y;
	for i=1:nax
		ax{i}=subplot(2,3,i);
		yyaxis(ax{i},'left');hold on;
	end

	c = brewermap(numel(experiment_names),'Set1');

	for i=1:numel(experiment_names)
		% load results
		fname = sprintf('issmResults_%s.mat',experiment_names{i});
		results = loadmodel(fullfile(datadir,fname));

		% VAF data
		t = [results.time];
		ice_maf_basin21 = ([results.ice_vaf_basin21] - results(1).ice_vaf_basin21).*rhoGt;
		ice_maf_basin22 = ([results.ice_vaf_basin22] - results(1).ice_vaf_basin22).*rhoGt;
		ice_maf_total = ice_maf_basin21+ice_maf_basin22;
		h_total   = plot(ax{1},t,ice_maf_total,  '-','color',c(i,:),'linewidth',2);
	   h_basin21 = plot(ax{2},t,ice_maf_basin21,'-','color',c(i,:),'linewidth',2);
		h_basin22 = plot(ax{3},t,ice_maf_basin22,'-','color',c(i,:),'linewidth',2);

		plot(ax{1+nax_x},t,[results.melt]        *1000*1E-12,'-','color',c(i,:))
		plot(ax{2+nax_x},t,[results.melt_basin21]*1000*1E-12,'-','color',c(i,:))
		plot(ax{3+nax_x},t,[results.melt_basin22]*1000*1E-12,'-','color',c(i,:))
	end

	for i=1:nax_x
		%legend(ax{i},experiment_names,'location','eastoutside')
		ylabel(ax{i},'\Delta mass above flotation (Gt)');
		ylimleft=ylim(ax{i});
      yyaxis(ax{i},'right');hold on;
      ylimright=ylimleft.*gt2mmslr;
      ylim(ax{i},ylimright);
      ylabel(ax{i},'sea level rise equivalence (mm)');
      xlim([2013,2100]);
      set(ax{i},'fontsize',14);

		ylabel(ax{i+nax_x},'total melt (Gt/yr)');
		xlabel(ax{i+nax_x},'time (y)');
		set(ax{i+nax_x},'xlim',[2013,2100],'fontsize',14);
   end

	return

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
if perform(org,'plotThwaites') % {{{
	% constants 
	rhoGt = 917*1e-12; % density of ice in Gt/m^3
	gt2mmslr = 1/361.8;  % 361.8 Gt of ice will raise global sea levels by ~1 mm

	% plot
	figure(1);clf;
	nax_x = 1;
	nax_y = 2;
	nax = nax_x*nax_y;
	for i=1:nax
		ax{i}=subplot(nax_y,nax_x,i);
		yyaxis(ax{i},'left');hold on;
	end

	c = brewermap(numel(experiment_names),'Set1');

	for i=1:numel(experiment_names)
		% load results
		fname = sprintf('issmResults_%s.mat',experiment_names{i});
		results = loadmodel(fullfile(datadir,fname));

		% VAF data
		t = [results.time];
		ice_maf_basin21 = ([results.ice_vaf_basin21] - results(1).ice_vaf_basin21).*rhoGt;
		ice_maf_basin22 = ([results.ice_vaf_basin22] - results(1).ice_vaf_basin22).*rhoGt;
		ice_maf_total = ice_maf_basin21+ice_maf_basin22;
	   %h_basin21 = plot(ax{1},t,ice_maf_basin21,'-','color',c(i,:),'linewidth',2);
	   h_basin21 = plot(ax{1},t,[0 diff(ice_maf_basin21)],'-','color',c(i,:),'linewidth',2);

		plot(ax{1+nax_x},t,[results.melt_basin21]*1000*1E-12,'-','color',c(i,:))
	end

	for i=1:nax_x
		%legend(ax{i},experiment_names,'location','eastoutside')
		ylabel(ax{i},'\Delta mass above flotation (Gt)');
		ylimleft=ylim(ax{i});
      yyaxis(ax{i},'right');hold on;
      ylimright=ylimleft.*gt2mmslr;
      ylim(ax{i},ylimright);
      ylabel(ax{i},'sea level rise equivalence (mm)');
      set(ax{i},'xlim',[2013,2100],'fontsize',14);

		ylabel(ax{i+nax_x},'total melt (Gt/yr)');
		xlabel(ax{i+nax_x},'time (y)');
		set(ax{i+nax_x},'xlim',[2013,2100],'fontsize',14);
   end
end % }}}
if perform(org,'plotAnomaly') % {{{
	% constants 
	rhoGt = 917*1e-12; % density of ice in Gt/m^3
	gt2mmslr = 1/361.8;  % 361.8 Gt of ice will raise global sea levels by ~1 mm

	% Plot anomaly of KN scenarios {{{
	ind_ctrl = 1; % index of control experiment for this anomaly (constant climate experiment)
	ind_var  = [2,12,13]; % index of variable experiments for this anomaly (monthly climatology, Paris2C, RCP85)

	% load control
	fname = sprintf('issmResults_%s.mat',experiment_names{ind_ctrl});
   results_ctrl = loadmodel(fullfile(datadir,fname));

	figure(1);clf;
	nax_x = 1;
	nax_y = 2;
	nax = nax_x*nax_y;
	for i=1:nax
		ax{i}=subplot(nax_y,nax_x,i);
		yyaxis(ax{i},'left');hold on;
	end

	c = brewermap(numel(experiment_names),'Set1');
	line_str = {'-','-.',':'};
	line_str_ind = [...
		1,1,...
		0,0,0,...
		0,0,0,...
		0,0,0,...
		1,1,...
		];
	c_ind = [...
		0,3,...
		0,0,0,...
		0,0,0,...
		0,0,0,...
		2,1,...
		];

	for j=1:numel(ind_var)	
		i = ind_var(j);
		% load results
		fname = sprintf('issmResults_%s.mat',experiment_names{i});
		results_var = loadmodel(fullfile(datadir,fname));

		% plotting data
		end_ind = min(numel(results_ctrl),numel(results_var)); % only plot to the minimum shared year of the data
		t = [results_var(1:end_ind).time]; % time vectory for these data (y)
		maf_anom_basin21  = [results_var(1:end_ind).ice_vaf_basin21] - [results_ctrl(1:end_ind).ice_vaf_basin21]; % MAF anomaly (Gt)
		melt_anom_basin21 = ([results_var(1:end_ind).melt_basin21] - [results_ctrl(1:end_ind).melt_basin21]) *1000*1E-12; % melt rate anomaly (Gt/yr)

		% plot
		yline(ax{1},0,'--k');
		h_basin21(j) = plot(ax{1},t,maf_anom_basin21,...
			'Marker','none',...
         'LineStyle',line_str{line_str_ind(i)},...
         'Color',c(c_ind(i),:),...
         'LineWidth',2);
		yline(ax{1+nax_x},0,'--k');
		plot(ax{1+nax_x},t,melt_anom_basin21,...
			'Marker','none',...
			'LineStyle',line_str{line_str_ind(i)},...
			'Color',c(c_ind(i),:),...
			'LineWidth',2)
		%% plot
		%yline(ax{1},0,'--k');
		%h_basin21(j) = plot(ax{1},t,maf_anom_basin21,'-','color',c(i,:),'linewidth',2);
		%yline(ax{1+nax_x},0,'--k');
		%plot(ax{1+nax_x},t,melt_anom_basin21,'-','color',c(i,:))


		%ice_maf_basin21 = ([results.ice_vaf_basin21] - results(1).ice_vaf_basin21).*rhoGt;
		%ice_maf_basin22 = ([results.ice_vaf_basin22] - results(1).ice_vaf_basin22).*rhoGt;
		%ice_maf_total = ice_maf_basin21+ice_maf_basin22;
	   %h_basin21 = plot(ax{1},t,ice_maf_basin21,'-','color',c(i,:),'linewidth',2);

	   %h_basin21 = plot(ax{1},t,[0 diff(ice_maf_basin21)],'-','color',c(i,:),'linewidth',2);

		%plot(ax{1+nax_x},t,[results.melt_basin21]*1000*1E-12,'-','color',c(i,:))
	end

	legend(h_basin21,experiment_names(ind_var))

	for i=1:nax_x
		%legend(ax{i},experiment_names,'location','eastoutside')
		ylabel(ax{i},'\Delta mass above flotation (Gt)');
		ylimleft=ylim(ax{i});
      yyaxis(ax{i},'right');hold on;
      ylimright=ylimleft.*gt2mmslr;
      ylim(ax{i},ylimright);
      ylabel(ax{i},'sea level rise equivalence (mm)');
      set(ax{i},'xlim',[2013,2100],'fontsize',14);

		ylabel(ax{i+nax_x},'total melt (Gt/yr)');
		xlabel(ax{i+nax_x},'time (y)');
		set(ax{i+nax_x},'xlim',[2013,2100],'fontsize',14);
   end 
	% }}}
	% Plot anomaly of sensitivity scenarios {{{
	ind_ctrl = 3; % index of control experiment for this anomaly (constant climate experiment)
	ind_var  = [4,5,6,  8,9,  11]; % index of variable experiments for this anomaly (monthly climatology, Paris2C, RCP85)

	% load control
	fname = sprintf('issmResults_%s.mat',experiment_names{ind_ctrl});
   results_ctrl = loadmodel(fullfile(datadir,fname));

	figure(2);clf;
	nax_x = 1;
	nax_y = 2;
	nax = nax_x*nax_y;
	for i=1:nax
		ax{i}=subplot(nax_y,nax_x,i);
		yyaxis(ax{i},'left');hold on;
	end

	c = brewermap(numel(experiment_names),'Set1');
	line_str = {'-','-.',':'};
	line_str_ind = [...
		0,0,...
		1,1,1,...
		2,2,2,...
		3,3,3,...
		0,0,...
		];
	c_ind = [...
		0,0,...
		0,2,1,...
		3,4,5,...
		3,4,5,...
		0,0,...
		];

	for j=1:numel(ind_var)	
		i = ind_var(j);
		% load results
		fname = sprintf('issmResults_%s.mat',experiment_names{i});
		results_var = loadmodel(fullfile(datadir,fname));

		% plotting data
		end_ind = min(numel(results_ctrl),numel(results_var)); % only plot to the minimum shared year of the data
		t = [results_var(1:end_ind).time]; % time vectory for these data (y)
		maf_anom_basin21  = [results_var(1:end_ind).ice_vaf_basin21] - [results_ctrl(1:end_ind).ice_vaf_basin21]; % MAF anomaly (Gt)
		melt_anom_basin21 = ([results_var(1:end_ind).melt_basin21] - [results_ctrl(1:end_ind).melt_basin21]) *1000*1E-12; % melt rate anomaly (Gt/yr)

		% plot
		yline(ax{1},0,'--k');
		h_basin21(j) = plot(ax{1},t,maf_anom_basin21,...
			'Marker','none',...
         'LineStyle',line_str{line_str_ind(i)},...
         'Color',c(c_ind(i),:),...
         'LineWidth',2);
		yline(ax{1+nax_x},0,'--k');
		plot(ax{1+nax_x},t,melt_anom_basin21,...
			'Marker','none',...
			'LineStyle',line_str{line_str_ind(i)},...
			'Color',c(c_ind(i),:),...
			'LineWidth',2)


		%ice_maf_basin21 = ([results.ice_vaf_basin21] - results(1).ice_vaf_basin21).*rhoGt;
		%ice_maf_basin22 = ([results.ice_vaf_basin22] - results(1).ice_vaf_basin22).*rhoGt;
		%ice_maf_total = ice_maf_basin21+ice_maf_basin22;
	   %h_basin21 = plot(ax{1},t,ice_maf_basin21,'-','color',c(i,:),'linewidth',2);

	   %h_basin21 = plot(ax{1},t,[0 diff(ice_maf_basin21)],'-','color',c(i,:),'linewidth',2);

		%plot(ax{1+nax_x},t,[results.melt_basin21]*1000*1E-12,'-','color',c(i,:))
	end

	legend(h_basin21,experiment_names(ind_var),'Location','nw')

	for i=1:nax_x
		%legend(ax{i},experiment_names,'location','eastoutside')
		ylabel(ax{i},'\Delta mass above flotation (Gt)');
		ylimleft=ylim(ax{i});
      yyaxis(ax{i},'right');hold on;
      ylimright=ylimleft.*gt2mmslr;
      ylim(ax{i},ylimright);
      ylabel(ax{i},'sea level rise equivalence (mm)');
      set(ax{i},'xlim',[2013,2100],'fontsize',18);

		ylabel(ax{i+nax_x},'total melt (Gt/yr)');
		xlabel(ax{i+nax_x},'time (y)');
		set(ax{i+nax_x},'xlim',[2013,2100],'fontsize',18);
   end 
	% }}}
end % }}}
