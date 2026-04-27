runid = 'RUN02' 
fname = sprintf('%s.mat',runid);
if exist(fname)
	OldFile = load(fname);
	nParis2C = numel(OldFile.vafParis2C) + 1;
	nRCP85 = numel(OldFile.vafRCP85) + 1;
else
	nParis2C = 0 + 1;
	nRCP85 = 0 + 1;
end

% load data
formatstr = '/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/%s/%s/results/';
parisDir = sprintf(formatstr,'Paris2C',runid);
rcpDir = sprintf(formatstr,'RCP85',runid);

load('/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/ISSM_initialization/Models/PROPHET_issm_init_InversionC.mat');
areas  =GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y);

disp('loading Paris2C results');
S = dir(fullfile(parisDir,'issmDiag*.mat'));
filenames=strcat(parisDir,{S.name});
vafParis2C = zeros(numel(filenames),1);
meltParis2C = zeros(numel(filenames),1);
vafParis2C(1:numel(OldFile.vafParis2C)) = OldFile.vafParis2C;
meltParis2C(1:numel(OldFile.meltParis2C)) = OldFile.meltParis2C;
for i=nParis2C:numel(filenames)
	load(filenames{i});
	vafParis2C(i)=results.IceVolumeAboveFloatation;
	meltElement=mean(results.BasalforcingsFloatingiceMeltingRate(md.mesh.elements),2);
	meltParis2C(i)=sum(meltElement.*areas);
	clear results;
end

disp('loading RCP85 results');
S = dir(fullfile(rcpDir,'issmDiag*.mat'));
filenames=strcat(rcpDir,{S.name});
vafRCP85 = zeros(numel(filenames),1);
meltRCP85 = zeros(numel(filenames),1);
vafRCP85(1:numel(OldFile.vafRCP85)) = OldFile.vafRCP85;
meltRCP85(1:numel(OldFile.meltRCP85)) = OldFile.meltRCP85;
for i=nRCP85:numel(filenames)
	load(filenames{i});
	vafRCP85(i)=results.IceVolumeAboveFloatation;
	meltElement=mean(results.BasalforcingsFloatingiceMeltingRate(md.mesh.elements),2);
	meltRCP85(i)=sum(meltElement.*areas);
	clear results;
end
save(fname,'vafParis2C','meltParis2C','vafRCP85','meltRCP85');
