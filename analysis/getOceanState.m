% getOceanState reads MITgcm pickup files and saves the output for the thermal forcing

% directories to pull from
proph_dir = '/nobackupp18/bgetraer/issmjpl/proj-getraer/proj-PROPHET';
mit=loadmodel(fullfile(proph_dir,'experiments/MITgcm_initialization/Models/PROPHET_mitgcm_init_RuntimeOptionsCoupled.mat'));
experiments = ["Paris2C","RCP85"];
experiment_dir = fullfile(proph_dir,'experiments',experiments,'runcoupledRUN02_dt100_ct1296000');

% query points at Thwaites Ice Shelf
xq = [-1578409, -1579897, -1584359, -1591052, -1597745]; % approximate query location in horizontal x (m)
yq = [-397306,  -410346,  -422014,  -431622,  -437799];  % approximate query location in horizontal y (m)
% find intersection with the mesh
for i=1:numel(xq)
	[~,xind(i)] = min(abs(mit.mesh.xc-xq(i))); % exact query index on the horizontal x grid
	[~,yind(i)] = min(abs(mit.mesh.yc-yq(i))); % exact query index on the horizontal y grid
end
% find the linear indices for the full columns
ind=sub2ind([numel(mit.mesh.xc),numel(mit.mesh.yc),numel(mit.mesh.zc)],... 
	repmat(xind,numel(mit.mesh.zc),1),...
	repmat(yind,numel(mit.mesh.zc),1),...
	repmat([1:numel(mit.mesh.zc)]',1,numel(xind)));

% time vector at which to extract
month_q = 10; % extract only the month of October
year_q = 0:90; % years to extract 
month2s = 30*24*3600; 
year2s = 360*24*3600;
ts = year_q*year2s + month_q*month2s; % timestep vector
% initialize data structures
T=struct();
T.domainavg=[];
T.obw=[];
T.thwaites=[];
t = [];
for i = 1:numel(experiments)
	disp(experiments(i));
	for j = 1:numel(ts)
		fname=sprintf('pickup.save.%010.0f',ts(j));
		D=dir(fullfile(experiment_dir(i),[fname '*']));
		if ~isempty(D)
			disp(D(1).name)
			PickupData=rdmds(fullfile(experiment_dir{i},fname));
			temp = PickupData(:,:,(1:mit.mesh.Nz)+2*mit.mesh.Nz); % Temperature state (deg C)
			temp(temp==0)=NaN;
			T.domainavg(:,end+1) = squeeze(nanmean(temp,[1,2]));
			T.obw(:,end+1)       = squeeze(nanmean(temp(1:10,:,:),[1,2]));
			T.thwaites (:,end+1) = squeeze(nanmean(temp(ind),[2]));
			t(end+1) = ts(j);
		end
	end
	fname = sprintf('thermocline%s',experiments(i));
	save(fullfile(proph_dir,'data',fname),'T','t');
end

