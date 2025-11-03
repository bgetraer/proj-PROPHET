root_dir        = '/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/';
dest_root_dir   = '/u/bgetraer/backup/proj-proph/experiments';
exp_directories = ['Paris2C','RCP85', fullfile('sensitivity_experiments',{dir(fullfile(root_dir,'sensitivity_experiments','se*')).name})];
%testref_dir = fullfile(root_dir,'test_ref');
%test_dir = fullfile(root_dir,'test');
%command = sprintf('cp -r %s %s',testref_dir,test_dir);
%system(command);
%exp_directories = {'test'};

shiftID = nan(numel(exp_directories),1);

for i=1:numel(exp_directories)
	% find the directory
	fprintf('Directory name: %s\n', exp_directories{i});
	source = fullfile(root_dir,exp_directories{i});
	%destination = fullfile(dest_root_dir,exp_directories{i});
	destination = dest_root_dir;
	assert(isdir(source),'Source directory not found')

	% delete .old folders
	old_dir = dir(fullfile(source,'*.old'));
	for j=1:numel(old_dir)
		folder_path = fullfile(source,old_dir(j).name);
		fprintf('Removing %s\n',folder_path);
		rmdir(folder_path,'s');
	end

	% sync /nobackup to lou
	command = sprintf('shiftc --sync -r %s lou:%s', source, destination);
	[status,cmdout] = system(command);
	shiftID(i) = str2num(extractBetween(string(cmdout),'Shift id is ',char(10)));
end

return
%% 
% save the log
log_dir='/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET/lou_scripts/shiftc_log/';
prefix = 'shiftc_ID_';
A = [dir(fullfile(log_dir,[prefix '*'])).name];
if isempty(A)
	nfile = 1;
else
	nfile = max(str2num(extractAfter(A,prefix)))+1;
end

% confirm sync is done
shiftState = ones(numel(exp_directories),1);
while any(shiftState==1)
	for i=1:numel(shiftID)
		if shiftState(i)==1
			command = sprintf('shiftc --status --id=%i',shiftID(i));
			[~,cmdout] = system(command);
			shiftStatus=textscan(cmdout,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','Delimiter','|','EndOfLine','\n','HeaderLines',3);
			if any(contains(shiftStatus{1},'run'))
				shiftState(i)=1;
			elseif any(contains(shiftStatus{1},'error'))
				shiftState(i)=999;
			else
				shiftState(i)=0;
			end
		end
	end
end

% identify files you do not need from /nobackup

% remove files from /nobackup
%for i=1:numel(shiftState)
%	assert(shiftState(i)==0,sprintf('Unsuccessful shiftc for shiftID=%i',shiftID(i)));
%	fprintf('Succesful shiftc for shiftID=%i, removing source\n',shiftID(i));
%	source = fullfile(root_dir,exp_directories{i});
%	=rmdir(source,'s');
%end
