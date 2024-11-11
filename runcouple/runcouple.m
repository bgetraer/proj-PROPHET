function runcouple(mdfile,mitfile)
%RUNCOUPLE is intended as a standalone script to run a coupled ISSM-MITGCM model with MCC compilation/
%The inputs are mdfile which points to the location of the ISSM model file, and mitfile, which points
%to the location of the mit model file.
%These variables are declared explicitly and named in this script as the executable needs to expect them in order to load them.
%RUNCOUPLE loads the existing environment variables, loops through the time steps calling the models runs,
%and saves the output. Is is assumed that you are already located within the mitgcm "run" directory. 
dispMITxISSM();
disp('************************************************************************************');
disp('*   - beginning RUNCOUPLE mcc deployable');
disp(['*   - current directory is ' pwd])

%declare all variables and classes we need to load from input
mit=struct();
md=model();
md.friction=frictionschoof();
md.timestepping=timesteppingadaptive();
md.inversion=m1qn3inversion();

%load model structures
load(mdfile); % ISSM model
load(mitfile); % MITgcm model

% Set parameters outside of the loop
npMIT=mit.build.SZ.nPx*mit.build.SZ.nPy; % number of processors
md.cluster=generic('name',oshostname(),'np',npMIT); % set number of processors for ISSM. 'name' will be filled at runtime
md_prefix = 'runcouple';

%load initial draft and save mask of ice cover
bathy=binread(mit.fname.bathyfile,8,[mit.mesh.Nx,mit.mesh.Ny]);
draft=binread(mit.fname.draftfile,8,[mit.mesh.Nx,mit.mesh.Ny]);

%DEBUG CODE!!!!!
md.timestepping.final_time=mit.timestepping.coupledTimeStep/md.constants.yts;

%loop through each coupled step, run the models, save the ouput 
% n is the coupled step number we are STARTING FROM, from 0:nsteps-1
% niter is the MITgcm step number
% to start from an advanced state, all you need to do is set the 
% niter0 parameter to start at the niter you want.
% Example: if n=0, we are starting the first coupled step from niter0. After running the MITgcm model, niter advances by
% the appropriate number of MITgcm timesteps, and those results are loaded for the ice model.

niter=mit.inputdata.PARM{3}.nIter0; % starting niter
for n=0:(mit.timestepping.nsteps-1);
	display(['COUPLED STEP ' num2str(n+1) '/' num2str(mit.timestepping.nsteps)]);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Update draft
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('  reading ice shelf draft from ISSM');
	% interpolate ice draft and masks from ISSM
	newdraft = zeros(mit.mesh.Ny,mit.mesh.Nx);    % initialize matrix for ice draft depth (m)
	newdraft(:)=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.geometry.base,mit.mesh.hXC(:),mit.mesh.hYC(:),'default',0); % m
	mask_oce=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ocean_levelset,mit.mesh.hXC(:),mit.mesh.hYC(:),'default',-1); % -1 ocean, 1 grounded
	mask_ice=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ice_levelset,mit.mesh.hXC(:),mit.mesh.hYC(:),'default',1); % -1 ice, 1 no ice
	newdraft(mask_oce>0)=mit.geometry.bathy(mask_oce>0); % set all grounded ice to have a draft equal to the bathymetry (m)
   newdraft(mask_ice>0)=0; % set all open ocean to have zero draft (m)
	newdraft(1,:) = mit.geometry.draftOBS; % set draft at bottom boundary (m)
   newdraft(:,1) = mit.geometry.draftOBW; % set draft at left boundary (m)
	
	% WORK IN COL, ROW LIKE MITGCM
	newdraft=permute(newdraft,[2,1]); % permute indices

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Read ocean pickup file
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('  reading ocean pickup file');
	pickupSuff=getpickup(pwd,niter); % find the right pickup file suffix (ckptA or ckptB)
	pickup_fname=['pickup.' pickupSuff '.data']; % the data filename
	PickupData=binread(pickup_fname,8,[mit.mesh.Nx, mit.mesh.Ny, 6*mit.mesh.Nz+3]); % read the whole file
	T=PickupData(:,:,(1:mit.mesh.Nz)+2*mit.mesh.Nz); % Temperature state (deg C)
	S=PickupData(:,:,(1:mit.mesh.Nz)+3*mit.mesh.Nz); % Salinity state (g/kg)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Open new cells as necessary
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('  opening new melt cells');
	% find indices of locations where ice shelf retreated
	hFacC=binread('hFacC.data',4,[mit.mesh.Nx, mit.mesh.Ny, mit.mesh.Nz]); % hFacC (m)
	hCol=sum(hFacC,3); % water column height (m)
   [iw jw]=find(hCol>0); % horizontal indices where there is water
	[im jm] = find(newdraft>draft & newdraft>mit.mesh.zp(end)); % horizontal indices where there is melt

	disp(['  found ' num2str(numel(im)) ' melt cells']);
	%Extrapolate T/S to locations where ice shelf retreated
	for i=1:length(im)
		% first try vertical extrapolation
		kw=find(hFacC(im(i),jm(i),:)); % the vertical index where there is water
		if ~isempty(kw)
			% do we need to open a new cell?
			k_old=find(draft(im(i),jm(i))<=mit.mesh.zp(1:end),1,'last');
			k_new=find(newdraft(im(i),jm(i))<=mit.mesh.zp(1:end),1,'last');
			% new ocean cell where draft falls in different vertical cell
			if k_new~=k_old
				S(im(i),jm(i),1:min(kw)) = S(im(i),jm(i),min(kw));
				T(im(i),jm(i),1:min(kw)) = T(im(i),jm(i),min(kw));
			end
		else	%If not succesful, use closest neighbor horizontal extrapolation
			[~,ind]=min((iw-im(i)).^2+(jw-jm(i)).^2);
			salt_profile=squeeze(S(iw(ind),jw(ind),:)); % salinity profile of closest neighbor
			temp_profile=squeeze(T(iw(ind),jw(ind),:)); % temperature profile of closest neighbor
			kw=find(hFacC(iw(ind),jw(ind),:)); % the vertical index where there is water
			salt_profile(1:min(kw))=salt_profile(min(kw)); % extrapolate salinity profile to top
			temp_profile(1:min(kw))=temp_profile(min(kw)); % extrapolate temperature profile to top
			salt_profile(max(kw):end)=salt_profile(max(kw)); % extrapolate salinity profile to bottom
			temp_profile(max(kw):end)=temp_profile(max(kw)); % extrapolate temperature profile to bottom
			S(im(i),jm(i),:)=salt_profile; % set salinity for new ocean column
			T(im(i),jm(i),:)=temp_profile; % set salinity for new ocean column
		end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Write updated MITgcm files
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('  write updated ocean files');
	% updated pickup file
	PickupData(:,:,(1:mit.mesh.Nz)+2*mit.mesh.Nz)=T;
	PickupData(:,:,(1:mit.mesh.Nz)+3*mit.mesh.Nz)=S;
	binwrite(pickup_fname,PickupData,8);
	% updated draft file
   binwrite(mit.fname.draftfile,newdraft,8);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Update the data files
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% the cal start time, model start time, and obcs all stay the same. the only difference is 
	% the interation step number we pickup from.
	newline=['  nIter0=' num2str(niter) ','];
   command=['sed "s/.*nIter0.*/' newline '/" data > data.temp; mv data.temp data'];
	system(command);
	% update the pickup suffix to read the correct file
	newline=['  pickupSuff=' pickupSuff ','];
   command=['sed "s/.*pickupSuff.*/' newline '/" data > data.temp; mv data.temp data'];
   system(command);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Run MITgcm
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('  running MITgcm')
	system(['mpirun -np ' int2str(npMIT) ' ./mitgcmuv > out 2> err']);
	niter=mit.inputdata.PARM{3}.nIter0 + (n+1)*mit.timestepping.coupledTimeStep/mit.inputdata.PARM{3}.deltaT; % updated niter
	disp('  done MITgcm')

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Get melt from MITgcm
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('  read melt from MITgcm')
	melt_fname=sprintf('SHICE_fwFluxtave.%010i.data',niter); % melt file
   meltq_mitgcm = binread(melt_fname,4,[mit.mesh.Nx, mit.mesh.Ny]); % melt flux at cell centers (kg/m^2/s)
	meltq_mitgcm=permute(meltq_mitgcm,[2,1]);  % put in ROW COL order (kg/m^2/s)
	meltq_issm=InterpFromGridToMesh(mit.mesh.xc(:),mit.mesh.yc(:),meltq_mitgcm,md.mesh.x,md.mesh.y,0); % melt flux at vertices (kg/m^2/s)

   %Set basal melting rate fields
	md.basalforcings.floatingice_melting_rate=-meltq_issm*md.constants.yts/md.materials.rho_ice; % melt rate at element vertices (m/yr)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Run ISSM
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('  running ISSM')
   %Solve
   md.miscellaneous.name = sprintf('%s%010i',md_prefix,niter);
   md=solve(md,'transient');
	disp('  done ISSM')

	%Reset model
   md.geometry.base             = md.results.TransientSolution(end).Base;
   md.geometry.surface          = md.results.TransientSolution(end).Surface;
   md.geometry.thickness        = md.geometry.surface-md.geometry.base;
   md.initialization.vx         = md.results.TransientSolution(end).Vx;
   md.initialization.vy         = md.results.TransientSolution(end).Vy;
   md.initialization.vel        = md.results.TransientSolution(end).Vel;
   md.initialization.pressure   = md.results.TransientSolution(end).Pressure;
   md.mask.ocean_levelset       = md.results.TransientSolution(end).MaskOceanLevelset;

	%Save ISSM results
   fname = sprintf('issmDiag.%010i.mat', niter);
   results = md.results.TransientSolution(end);
	results.step = (n+1);
	results.time = (n+1)*mit.timestepping.coupledTimeStep;
   save(fname,'results');
   clear md.results;

	%Clear execution folder in ISSM to avoid going over quota
	command = ['rm -r ' md.cluster.executionpath '/' md_prefix '*'];
	system(command);
end

% subfunctions
function [pickupSuff]=getpickup(parentdir,nIter0) % {{{
%GETPICKUP finds the pickup file in the parentdir that matches the niter number
   pickupSuffs = {'ckptA','ckptB'}; % the .meta file names
	nit=[0,0];
   for i=1:numel(pickupSuffs)
      fid=fopen(fullfile(parentdir,['pickup.' pickupSuffs{i} '.meta']),'r');
      if fid~=-1
         tline=fgetl(fid); % read the next line
         while ischar(tline)
            if contains(tline, 'timeStepNumber')
               break;
            end
            tline = fgetl(fid); % read the next line
         end
         fclose(fid);
         nit(i)=str2num(extractBefore(extractAfter(tline,'['),']'));
		end
   end
   pickup_ind=find(nit==nIter0); % match the right pickup file
   if isempty(pickup_ind)
      error('No pickup file is found for nIter0!');
	end
	pickupSuff=pickupSuffs{pickup_ind};
end % }}}
function D=binread(fname,prec,arrsize) % {{{
% read data from binary file into a matlab array D.
% Assumes big-endian architecture, and given precision
% and array size.
%
% fname: filename or path (string)
% prec: 4 or 8 for number of bits
% arrsize: dimensions of D (array)
%
% D: array of requested dimension
	fid=fopen(fname,'r','b');
	switch prec
		case 8
			D=fread(fid,inf,'real*8');
		case 4
			D=fread(fid,inf,'real*4');
		otherwise
			error('give precision of data');
	end
	D=reshape(D,arrsize);
	fclose(fid);
end % }}}
function q=binwrite(fname,D,prec) % {{{
% write a matlab array D of arbitrary dimension to binary file
% using big-endian architecture and given precision.
%
% fname: filename or path (string)
% D: array of arbitrary dimension (storage is independent of dimension
% sizes)
% prec: 4 or 8 for number of bits
	fid=fopen(fname,'w','b');
	switch prec
		case 8
			q=fwrite(fid,D,'real*8');
		case 4
			q=fwrite(fid,D,'real*4');
		otherwise
			error('use valid precision');
	end
	fclose(fid);
end % }}}
function dispMITxISSM() % {{{
%DISPMITXISSM prints an ascii graphic to terminal output
   disp(' __      __   _   _______                    _     _  _____   _____ ____ __      __ ')
   disp('|   \  /   | | | |__   __|                   \ \ / / |_   _|/ ____/ ____|   \  /   |')
   disp('| |\ \/ /| | | |    | | __ _  ___ _ __ ___    \   /    | |  \___ \\___ \| |\ \/ /| |')
   disp('| | \__/ | | | |    | |/ _` |/ __| `_ ` _ \   /   \   _| |_ ____| |___| | | \__/ | |')
   disp('|_|      |_| |_|    |_| (_| | (__| | | | | | /_/ \_\ |_____|_____/_____/|_|      |_|')
   disp('                       \___ |\___|_| |_| |_| ')
   disp('                        __/ |                ')
   disp('                       |___/                 ')
end % }}}
end
