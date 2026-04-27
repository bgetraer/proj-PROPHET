steps = 1;

% directories
projdir = '/u/bgetraer/backup/proj-PROPHET';
analysisdir = fullfile(projdir,'analysis');
datadir = fullfile(projdir,'data');
expdir = fullfile(projdir,'exp');

experiment_dirs = { ...
'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se002_KN_constant_clim/runcoupled/', ...
'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se003_KN_monthly_clim/runcoupled/', ...
'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se004_PW700_constant/runcoupled/', ...
'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se005_PW600_constant/runcoupled/', ...
'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se006_PW800_constant/runcoupled/', ...
'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se007_PW700_amp50_per2/runcoupled/', ...
'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se008_PW700_amp50_per5/runcoupled/', ...
'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se009_PW700_amp50_per10/runcoupled/', ...
'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se010_PW700_amp100_per2/runcoupled/', ...
'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se011_PW700_amp100_per5/runcoupled/', ...
'/u/bgetraer/backup/proj-PROPHET/experiments/sensitivity_experiments/se012_PW700_amp100_per10/runcoupled/', ...
'/u/bgetraer/backup/proj-PROPHET/experiments/Paris2C/RUN02/runcoupled/',...
'/u/bgetraer/backup/proj-PROPHET/experiments/RCP85/RUN02/runcoupled/',...
};
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

% load initial model
mdpath = fullfile(projdir,'/experiments/ISSM_initialization/Models/PROPHET_issm_init_MeshParam.mat');
md = loadmodel(mdpath);

% loop over experiments
y2s = 360*24*3600; % seconds per year
dt = 25; % extract every dt years 
year_0 = 2010;
start_year = 2025;
end_year = 2300;
target_modeltimes=([start_year:dt:end_year]-year_0)*y2s;
exp_ind = 13;
for i=exp_ind
	fprintf('Experiment: %s\n',experiment_names{i});
	% issmDiag files
	resultsdir = experiment_dirs{i};
%	files = {dir(fullfile(resultsdir,'issmDiag.*')).name};
%	modeltimes = sort(str2double(extractBetween(files, 'issmDiag.', '.mat')));
	% check for issmGrLine file
	issmResults_fname = sprintf('issmGrLine_%s.mat',experiment_names{i});
%	if exist(fullfile(datadir,issmResults_fname))
%		fprintf('   loading existing results from %s\n',issmResults_fname);
%		results = loadmodel(fullfile(datadir,issmResults_fname));
%
%		ind = find(~any(modeltimes'==[results.modeltime],2));
%		j0 = numel(results);
%	else
%		clear results;
%
%		ind = find(modeltimes');
%		j0 = 0;
%	end

%	fprintf('Found %i saved timesteps\n',j0);
%	fprintf('Found %i new timesteps\n',numel(ind));
	
	j0=0;
	% loop through all target timesteps
	for j = 1:numel(target_modeltimes)
		fname=sprintf('issmDiag.%010i.mat',target_modeltimes(j));
		if exist(fullfile(datadir,fname))
			fprintf('   step %i/%i...\n',j,numel(target_modeltimes));
		else
			continue;
		end
		% load results into model
		md.results.TransientSolution = loadmodel(fname);
		results(j0+j).modeltime = md.results.TransientSolution.time; % model time (seconds)
		results(j0+j).time = md.results.TransientSolution.time/3600/24/360 + 2010;; % calendar time since 2010 (seconds)
		results(j0+j).ice_vaf = md.results.TransientSolution.IceVolumeAboveFloatation; % total domain volume above flotation (m^3)
		results(j0+j).groundingline = md.results.TransientSolution.
	end
	% save concatenated results
	fprintf('   saving results to %s\n',fullfile(datadir,issmResults_fname));
	save(fullfile(datadir,issmResults_fname),'results','-v7.3');
	disp('   done');
end % }}}
