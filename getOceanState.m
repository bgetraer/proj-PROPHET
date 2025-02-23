mit=loadmodel('experiments/MITgcm_initialization/Models/PROPHET_mitgcm_init_RuntimeOptionsCoupled.mat');
expdir1='experiments/Paris2C/runcoupled_dt100_ct1296000/';
expdir2='/u/bgetraer/backup/proj-proph/experiments/Paris2C/runcoupled/';


xq = [-1578409, -1579897, -1584359, -1591052, -1597745];
yq = [-397306,  -410346,  -422014,  -431622,  -437799];

for i=1:numel(xq)
	[~,xind(i)] = min(abs(mit.mesh.xc-xq(i)));
	[~,yind(i)] = min(abs(mit.mesh.yc-yq(i)));
end

ind=sub2ind([numel(mit.mesh.xc),numel(mit.mesh.yc),numel(mit.mesh.zc)],...
	repmat(xind,numel(mit.mesh.zc),1),...
	repmat(yind,numel(mit.mesh.zc),1),...
	repmat([1:numel(mit.mesh.zc)]',1,numel(xind)));


ts=((0:70)*360*24*3600) + 10*30*24*3600;
T=struct();
T.domainavg=[];
T.obw=[];
T.thwaites=[];
t = [];
for i=1:numel(ts)
	fname=sprintf('pickup.save.%010.0f',ts(i));
	D=dir(fullfile(expdir1,[fname '*']));
	if ~isempty(D)
		disp(D(1).name)
		PickupData=rdmds(fullfile(expdir,fname));
		temp = PickupData(:,:,(1:mit.mesh.Nz)+2*mit.mesh.Nz); % Temperature state (deg C)
		temp(temp==0)=NaN;
		T.domainavg(:,end+1) = squeeze(nanmean(temp,[1,2]));
		T.obw(:,end+1)       = squeeze(nanmean(temp(1:10,:,:),[1,2]));
		T.thwaites (:,end+1) = squeeze(nanmean(temp(ind),[2]));
		t(end+1) = ts(i);
	else
		D=dir(fullfile(expdir2,[fname '*']));
		if ~isempty(D)
			disp(D(1).name)
			PickupData=rdmds(fullfile(expdir,fname));
			temp = PickupData(:,:,(1:mit.mesh.Nz)+2*mit.mesh.Nz); % Temperature state (deg C)
			temp(temp==0)=NaN;
			T.domainavg(:,end+1) = squeeze(nanmean(temp,[1,2]));
			T.obw(:,end+1)       = squeeze(nanmean(temp(1:10,:,:),[1,2]));
			T.thwaites (:,end+1) = squeeze(nanmean(temp(ind),[2]));
			t(end+1) = ts(i);
		end
	end
end


