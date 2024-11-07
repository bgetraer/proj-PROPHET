%function runcouple(rundir,md,mit)
%RUNCOUPLE is intended as a standalone script to run a coupled ISSM-MITGCM model with MCC compilation/
%The only input is envfile.mat, which must contain all environment variables needed for this function,
%exported from the main runme.m. These variables are declared explicitly and named in this script as 
%the executable needs to expect them in order to load them.
%RUNCOUPLE loads the existing environment variables, loops through the time steps calling the models runs,
%and saves the output. Is is assumed that you are already located within the mitgcm "run" directory. 
dispMITxISSM();
disp('************************************************************************************');
disp('*   - beginning RUNCOUPLE mcc deployable');
disp(['*   - current directory is ' pwd])

Nx=mit.mesh.Nx;
Ny=mit.mesh.Ny;
Nz=mit.mesh.Nz;

% load initial draft and save mask of ice cover
fname=fullfile(rundir,mit.fname.bathyfile);
bathy=binread(fname,8,[Nx,Ny]);
fname=fullfile(rundir,mit.fname.draftfile);
draft=binread(fname,8,[Nx,Ny]);
mask_ice=zeros(size(draft));
mask_ice(draft>0)=-1;
mask_ice(draft==0)=1;

hFacC_fname=fullfile(rundir,'hFacC.data'); % hFacC file
npMIT=mit.build.SZ.nPx*mit.build.SZ.nPy;

md.cluster=generic('name',oshostname(),'np',npMIT); % set number of processors for ISSM. 'name' will be filled at runtime
md_prefix = 'runcouple';
%loop through each coupled step, run the models, save the ouput 
% n is the coupled step number we are STARTING FROM, from 0:nsteps-1
% niter is the MITgcm step number
% to start from an advanced state, all you need to do is set the 
% niter0 parameter to start at the niter you want.
% Example: if n=0, we are starting the first coupled step from niter0. After running the MITgcm model, niter advances by
% the appropriate number of MITgcm timesteps, and those results are loaded for the ice model.

niter=mit.inputdata.PARM{3}.nIter0; % starting niter
for n=0:(nsteps-1);
	display(['COUPLED STEP ' num2str(n+1) '/' num2str(nsteps)]);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Update draft
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('reading ice shelf draft from ISSM');
	% interpolate ice draft and masks from ISSM
	newdraft = zeros(mit.mesh.Ny,mit.mesh.Nx);    % initialize matrix for ice draft depth (m)
	newdraft(:)=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.geometry.base,mit.mesh.xc(:),mit.mesh.yc(:),0); % (def 0) (m)
	mask_oce=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ocean_levelset,mit.mesh.xc(:),mit.mesh.yc(:),-1); % -1 ocean, 1 grounded (def -1)
	newdraft(mask_oce>0)=bathy(mask_oce>0); % set all grounded ice to have a draft equal to the bathymetry (m)
	newdraft(draft==0)=0; % set all open ocean to have zero draft (m)
	newdraft(1,:) = mit.geometry.draftOBS; % set draft at bottom boundary (m)
	newdraft(:,1) = mit.geometry.draftOBW; % set draft at left boundary (m)
	newdraft=permute(newdraft,[2,1]); % WORK IN COL, ROW LIKE MITGCM

	newdraft=draft;
	newdraft(end-1,end-1)=draft(end-1,end-1)+10;
	newdraft(167,33)=0;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Read ocean pickup file
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[~,pickup_fname]=getpickup(rundir,niter); % find the right pickup file
	PickupData=binread(fullfile(rundir,pickup_fname),8,[Nx, Ny, 6*Nz+3]); % read the whole file
	whos PickupData
	T=PickupData(:,:,(1:Nz)+2*Nz); % Temperature state (deg C)
	S=PickupData(:,:,(1:Nz)+3*Nz); % Salinity state (g/kg)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Open new cells as necessary
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% find indices of locations where ice shelf retreated
	hFacC=binread(hFacC_fname,4,[Nx, Ny, Nz]); % hFacC (m)
	hCol=sum(hFacC,3); % water column height (m)
   [iw jw]=find(hCol>0); % horizontal indices where there is water
	[im jm] = find(newdraft>draft & newdraft>mit.mesh.zp(end)); % horizontal indices where there is melt

	%Extrapolate T/S to locations where ice shelf retreated
	for i=1:length(im)
		% first try vertical extrapolation
		kw=find(hFacC(im(i),jm(i),:)); % the vertical index where there is water
		if ~isempty(kw)
			% do we need to open a new cell?
			k_old=find(draft(im(i),jm(i))<=mit.mesh.zp(1:end),1,'last')
			k_new=find(newdraft(im(i),jm(i))<=mit.mesh.zp(1:end),1,'last')
			% new ocean cell where draft falls in different vertical cell
			if k_new~=k_old
				S(im(i),jm(i),1:min(kw)) = S(im(i),jm(i),min(kw));
				T(im(i),jm(i),1:min(kw)) = T(im(i),jm(i),min(kw));
				im(i)
				jm(i)
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
	% updated pickup file
	PickupData(:,:,(1:Nz)+2*Nz)=T;
	PickupData(:,:,(1:Nz)+3*Nz)=S;
	binwrite(fullfile(rundir,pickup_fname),PickupData,8);
	% updated draft file
   binwrite(fullfile(rundir,mit.fname.draftfile),newdraft,8);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Update niter in the data file
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% the cal start time, model start time, and obcs all stay the same. the only difference is 
	% the interation step number we pickup from.
	newline=['  nIter0=' num2str(niter) ','];
   command=['sed "s/.*nIter0.*/' newline '/" data > data.temp; mv data.temp data'];
	system(command);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Run MITgcm
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%system(['mpirun -np ' int2str(npMIT) ' ./mitgcmuv > out 2> err']);
	niter=mit.inputdata.PARM{3}.nIter0 + (n+1)*mit.timestepping.coupledTimeStep/mit.inputdata.PARM{3}.deltaT; % updated niter

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Get melt from MITgcm
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	melt_fname=sprintf('SHICE_fwFluxtave.%010i.data',niter); % melt file
   meltq_mitgcm = binread(fullfile(rundir,melt_fname),4,[Nx, Ny]); % melt flux at cell centers (kg/m^2/s)
	meltq_mitgcm=permute(melt_mitgcm,[2,1]);  % put in ROW COL order (kg/m^2/s)
	meltq_issm=InterpFromGridToMesh(mit.mesh.xc(:),mit.mesh.yc(:),melt_mitgcm,md.mesh.x,md.mesh.y,0); % melt flux at vertices (kg/m^2/s)

   %Set basal melting rate fields
	md.basalforcings.floatingice_melting_rate=-meltq_issm*md.constants.yts/md.materials.rho_ice; % melt rate at element vertices (m/yr)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Run ISSM
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Solve
   md.miscellaneous.name = sprintf('%s%010i',md_prefix,niter);
   md=solve(md,'transient');

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
   fname = sprintf('results/TransientSolution.%010i.mat', niter);
   results = md.results.TransientSolution(end);
	results.step = (n+1);
	results.time = (n+1)*mit.timestepping.coupledTimeStep;
   save(fname,'results');
   clear md.results;

	%Clear execution folder in ISSM to avoid going over quota
	command = ['rm -r ' md.cluster.executionpath '/' md_prefix '*'];
	system(command);
%end



% subfunctions
function [fname_meta,fname_data]=getpickup(parentdir,nIter0) % {{{
%GETPICKUP finds the pickup file in the parentdir that matches the niter number
   pickup_fnames = {'pickup.ckptA','pickup.ckptB'}; % the .meta file names
	nit=[0,0];
   for i=1:numel(pickup_fnames)
      fid=fopen(fullfile(parentdir,[pickup_fnames{i} '.meta']),'r');
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
	fname=pickup_fnames(pickup_ind);
	fname_meta=[fname{:} '.meta'];
	fname_data=[fname{:} '.data'];
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
%end
