steps = 1;

projdir = '/u/bgetraer/backup/proj-PROPHET';
analysisdir = fullfile(projdir,'analysis');
datadir = fullfile(projdir,'data');
expdir = fullfile(projdir,'exp');

experiment_dirs = { ...
%'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se002_KN_constant_clim/runcoupled/', ...
%'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se003_KN_monthly_clim/runcoupled/', ...
%'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se004_PW700_constant/runcoupled/', ...
%'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se005_PW600_constant/runcoupled/', ...
%'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se006_PW800_constant/runcoupled/', ...
%'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se007_PW700_amp50_per2/runcoupled/', ...
%'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se008_PW700_amp50_per5/runcoupled/', ...
%'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se009_PW700_amp50_per10/runcoupled/', ...
%'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se010_PW700_amp100_per2/runcoupled/', ...
%'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se011_PW700_amp100_per5/runcoupled/', ...
%'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se012_PW700_amp100_per10/runcoupled/', ...
'/u/bgetraer/backup/proj-PROPHET/experiments/Paris2C/RUN02/runcoupled/',...
'/u/bgetraer/backup/proj-PROPHET/experiments/RCP85/RUN02/runcoupled/',...
};
experiment_names = { ...
%'se002_KN_constant_clim', ...
%'se003_KN_monthly_clim', ...
%'se004_PW700_constant', ...
%'se005_PW600_constant', ...
%'se006_PW800_constant', ...
%'se007_PW700_amp50_per2', ...
%'se008_PW700_amp50_per5', ...
%'se009_PW700_amp50_per10', ...
%'se010_PW700_amp100_per2', ...
%'se011_PW700_amp100_per5', ...
%'se012_PW700_amp100_per10', ...
'Paris2C',...
'RCP85',...
};

% load initial model
mdpath = '/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/ISSM_initialization/Models/PROPHET_issm_init_MeshParam.mat';
md = loadmodel(mdpath);
% basin flags
flagfile = fullfile(expdir,'basin_flags.mat');
if exist(flagfile)
	disp('Loading sub-basin element flags');
	load(flagfile);
else
	disp('Flagging elements in sub-basins');
	flag_basin21=FlagElements(md,fullfile(expdir,'reg21_thwaites.exp')); % flags for basin 21
	flag_basin22=FlagElements(md,fullfile(expdir,'reg22_pineisland.exp')); % flags for basin 22
	fprintf('Saving flags to %s\n',flagfile);
	save(flagfile,'flag_basin21','flag_basin22');
end


% loop over each experiment
for i=1:numel(experiment_names)
	fprintf('Experiment: %s\n',experiment_names{i});
	% issmDiag files
	resultsdir = experiment_dirs{i};
	files = {dir(fullfile(resultsdir,'issmDiag.*')).name};
	modeltimes = sort(str2double(extractBetween(files, 'issmDiag.', '.mat')));
	% check for issmResults file
	issmResults_fname = sprintf('issmResults_%s.mat',experiment_names{i});
	if exist(fullfile(datadir,issmResults_fname))
		fprintf('   loading existing results from %s\n',issmResults_fname);
		results = loadmodel(fullfile(datadir,issmResults_fname));

		ind = find(~any(modeltimes'==[results.modeltime],2));
		j0 = numel(results);
	else
		clear results;

		ind = find(modeltimes');
		j0 = 0;
	end

	fprintf('Found %i saved timesteps\n',j0);
	fprintf('Found %i new timesteps\n',numel(ind));
	% loop through all files
	for j = 1:numel(ind)
		fprintf('   file %i/%i...\n',j,numel(ind));
		fname = fullfile(resultsdir,files{ind(j)});
		% load results into model
		md.results.TransientSolution = loadmodel(fname);
		results(j0+j).modeltime = md.results.TransientSolution.time; % model time (seconds)
		results(j0+j).time = md.results.TransientSolution.time/3600/24/360 + 2010;; % calendar time since 2010 (seconds)
		results(j0+j).ice_vaf = md.results.TransientSolution.IceVolumeAboveFloatation; % total domain volume above flotation (m^3)
		% calculate VAF for sub-basins
		results(j0+j).ice_vaf_basin21 = VolumeAboveFloatation(md,1,flag_basin21); % basin volume above flotation at TransientSolution(1) (m^3)
		results(j0+j).ice_vaf_basin22 = VolumeAboveFloatation(md,1,flag_basin22); % basin volume above flotation at TransientSolution(1) (m^3)
		% calculate melt for sub-basins
		meltrate = mean(md.results.TransientSolution.BasalforcingsFloatingiceMeltingRate(md.mesh.elements),2);
		results(j0+j).melt = sum(meltrate .* GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y)); % total melt (m^3 / yr)
		results(j0+j).melt_basin21 = sum(meltrate .* GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y) .* flag_basin21); % basin melt m^3/yr
		results(j0+j).melt_basin22 = sum(meltrate .* GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y) .* flag_basin22); % basin melt m^3/yr
	end
	% save concatenated results
	fprintf('   saving results to %s\n',fullfile(datadir,issmResults_fname));
	save(fullfile(datadir,issmResults_fname),'results','-v7.3');
	disp('   done');
end % }}}
