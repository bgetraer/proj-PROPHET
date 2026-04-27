runid = 'RUN02' 

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
for i=1:numel(filenames)
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
for i=1:numel(filenames)
	load(filenames{i});
	vafRCP85(i)=results.IceVolumeAboveFloatation;
	meltElement=mean(results.BasalforcingsFloatingiceMeltingRate(md.mesh.elements),2);
	meltRCP85(i)=sum(meltElement.*areas);
	clear results;
end

fname = sprintf('%s.mat',runid);
save(fname,'vafParis2C','meltParis2C','vafRCP85','meltRCP85');
