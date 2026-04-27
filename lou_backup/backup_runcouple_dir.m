root_dir        = '/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/';
dest_root_dir   = '/u/bgetraer/backup/proj-proph/experiments';
scenario_source_dir = {'Paris2C/runcoupledRUN02_dt100_ct1296000','RCP85/runcoupledRUN02_dt100_ct1296000'};
scenario_destination_dir = {'Paris2C/RUN02/runcoupled','RCP85/RUN02/runcoupled'};

% sync scenarios
%for i=1:numel(scenario_source_dir)
%	disp(scenario_source_dir{i})
%	source_dir      = fullfile(root_dir,scenario_source_dir{i});
%	destination_dir = fullfile(dest_root_dir,scenario_destination_dir{i});
%	sync_dir(source_dir,destination_dir);
%end

% sync sensitivity experiments
sensitivity_dirs = fullfile('sensitivity_experiments',{dir(fullfile(root_dir,'sensitivity_experiments','se*')).name},'runcoupled');
for i=1:numel(sensitivity_dirs)
	disp(sensitivity_dirs{i})
	source_dir      = fullfile(root_dir,sensitivity_dirs{i});
	destination_dir = fullfile(dest_root_dir,sensitivity_dirs{i});
	sync_dir(source_dir,destination_dir);
end

function sync_dir(source_dir,destination_dir)
	sync_list = '/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET/lou_backup/sync_list.txt';

	A1 = dir(fullfile(source_dir,'issmDiag*.mat'));
	t1 = sort(str2double(extractBetween({A1.name}, 'issmDiag.', '.mat')));

	A2 = dir(fullfile(source_dir,'*.save.*.*'));
	t2 = str2double(extractBetween({A2.name}, '.save.', '.'));

	t = unique([t1,t2]);
	fprintf('%i timesteps found to sync\n',numel(t));

	fid = fopen(sync_list,'w');
	nfiles = 0;
	for i=1:(numel(t)-1)
		fmtstr = sprintf('*.%010i.*',t(i));
		B = dir(fullfile(source_dir,fmtstr));
		names = {B.name};	
		fprintf(fid,'%s\n',names{:});
		nfiles=nfiles+numel(names);
	end
	fclose(fid);

	fprintf('%i files found to sync',nfiles);

	% sync /nobackup to lou
	command = sprintf('rsync --remove-source-files -v --files-from=%s %s lou:%s', sync_list, source_dir, destination_dir);
	disp('Syncing now and removing source files');
	system(command);
end
