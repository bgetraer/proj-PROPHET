function runcouple(mdfile,mitfile)
%RUNCOUPLE is a script to run a coupled ISSM-MITGCM model with MCC compilation/or on an interactive node.
%The inputs are mdfile which points to the location of the ISSM model file, and mitfile, which points
%to the location of the mit model file.
%These objexts (and some of their subobjects) are declared explicitly and named in this script as the 
%executable needs to know their class in order to load them.
%RUNCOUPLE loads the existing environment variables, loops through the time steps calling the models runs,
%and saves the output. Is is assumed that you are already located within the mitgcm "run" directory. 
%
% Example:
%    runcouple(mdfile,mitfile);

% opening display {{{
dispMITxISSM();
disp('************************************************************************************');
disp('*   - beginning RUNCOUPLE');
disp(['*   - current directory is ' pwd])
disp('************************************************************************************');
disp('');
% }}}
% parse inputs {{{
%declare all variables and classes we need to load from input
mit=struct();
md=model();
md.friction=frictionschoof();
md.timestepping=timesteppingadaptive();
md.inversion=m1qn3inversion();
%load model structures
load(mdfile); % ISSM model
load(mitfile); % MITgcm model
% }}}
% static parameters and fields {{{
% File names
draft_file = 'draft.bin';
bathy_file = 'bathy.bin';
uvel_file  = 'uvel.bin';
vvel_file  = 'vvel.bin';
theta_file = 'theta.bin';
salt_file  = 'salt.bin';
etan_file  = 'etan.bin';

% Model execution parameters
nprocs=mit.build.SZ.nPx*mit.build.SZ.nPy; % number of processors for MITgcm and ISSM
md.cluster=generic('name',oshostname(),'np',nprocs); % set number of processors for ISSM.
md.timestepping.final_time=mit.timestepping.coupledTimeStep./md.constants.yts; % how long to run ISSM for
md_prefix = 'runcouple'; % the ISSM model name prefix for execution files

% Static fields for opening and closing draft cells
bathy=binread(bathy_file,8,[mit.mesh.Nx,mit.mesh.Ny]); % the MITgcm bathymetry in col,row matrix (m)
mask_ice=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ice_levelset,mit.mesh.hXC(:),mit.mesh.hYC(:),'default',1); % ISSM ice mask (m)
mask_ice=permute(reshape(mask_ice,mit.mesh.Ny,mit.mesh.Nx),[2,1]); % ISSM ice mask in col,row matrix (m)

% Timekeeping parameters
% n is the coupled step number we are STARTING FROM, from 0:nsteps-1
% niter is the MITgcm step number
% modeltime is the MITgcm modeltime starting from the calendar start date (2010)
% coupled_basetime is the MITgcm modeltime that we start the coupling at (2013)
% basetime is the MITgcm modeltime that the current model starts at 
coupled_basetime=mit.inputdata.PARM{3}.startTime; % modeltime that we start coupling at
niter0=mit.inputdata.PARM{3}.nIter0; % starting niter
modelIterEnd=mit.timestepping.coupledTimeStep/mit.inputdata.PARM{3}.deltaT; % the final timestep number of each model run
% }}}
% initial parameters and fields {{{
if niter0==0
	deltaBase=zeros(size(md.geometry.base)); % the change in ice shelf draft from ISSM (initialize to zero) (m)
elseif niter0>0
	fname = sprintf('issmDiag.%010i.mat', niter0); % issm results to load into md
	load(fname); % load results structure
	md.geometry.base            = results.Base;
	md.geometry.surface         = results.Surface;
	md.geometry.thickness       = md.geometry.surface-md.geometry.base;
   md.initialization.vx        = results.Vx;
   md.initialization.vy        = results.Vy;
   md.initialization.vel       = results.Vel;
   md.mask.ocean_levelset      = results.MaskOceanLevelset;
	deltaBase                   = results.deltaBase;
end
% }}}

disp(['starting at niter=' num2str(niter0)]);
for n=0:(mit.timestepping.nsteps-1);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Update timekeeping
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	display(['COUPLED STEP ' num2str(n+1) '/' num2str(mit.timestepping.nsteps)]);
	modeltime = coupled_basetime+(n)*mit.timestepping.coupledTimeStep;                      % the current modeltime
	niter     = niter0 + (n)*mit.timestepping.coupledTimeStep/mit.inputdata.PARM{3}.deltaT; % the current niter

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Update draft
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	fname=sprintf('draft.save.%010i.bin',niter);
	disp(['  reading last MITgcm draft file ' fname]);
	draft=binread(fname,8,[mit.mesh.Nx,mit.mesh.Ny]); % existing MITgcm draft in col,row matrix (m)
	
	disp('  reading updated draft change from ISSM');
	deltaDraft=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,deltaBase,mit.mesh.hXC(:),mit.mesh.hYC(:),'default',0); % change in ISSM draft (m)
	deltaDraft=permute(reshape(deltaDraft,mit.mesh.Ny,mit.mesh.Nx),[2,1]); % change in ISSM draft in col,row matrix (m)
	newdraft=draft+deltaDraft; % updated MITgcm draft in col,row matrix (m)
	disp(['num deltaDraft cells = ' num2str(sum(deltaDraft(:)~=0))]);

	% interpolate ice draft and masks from ISSM
	%newdraft=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.geometry.base,mit.mesh.hXC(:),mit.mesh.hYC(:),'default',0); % m
	%newdraft=permute(reshape(newdraft,mit.mesh.Ny,mit.mesh.Nx),[2,1]); % ISSM draft in col,row matrix (m)

	mask_oce=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ocean_levelset,mit.mesh.hXC(:),mit.mesh.hYC(:),'default',-1); % -1 ocean, 1 grounded
	mask_oce=permute(reshape(mask_oce,mit.mesh.Ny,mit.mesh.Nx),[2,1]); % ISSM ocean mask in col,row matrix (m)

	newdraft(mask_oce>0)=bathy(mask_oce>0); % set all grounded ice to have a draft equal to the bathymetry (m)
   newdraft(mask_ice>0)=0;                 % set all open ocean to have zero draft (m)
	newdraft(:,1)=mit.geometry.draftOBS;  % set draft at bottom boundary (m)
   newdraft(1,:)=mit.geometry.draftOBW;  % set draft at left boundary (m)
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Read ocean pickup file
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	fname=sprintf('pickup.save.%010i.data',niter); % the data filename
	disp(['  reading ocean pickup file ' fname]);
	PickupData=binread(fname,8,[mit.mesh.Nx, mit.mesh.Ny, 6*mit.mesh.Nz+3]); % read the whole file
	U=PickupData(:,:,(1:mit.mesh.Nz)+0*mit.mesh.Nz); % x component of velocity (m/s)
	V=PickupData(:,:,(1:mit.mesh.Nz)+1*mit.mesh.Nz); % y component of velocity (m/s)
	T=PickupData(:,:,(1:mit.mesh.Nz)+2*mit.mesh.Nz); % Temperature state (deg C)
	S=PickupData(:,:,(1:mit.mesh.Nz)+3*mit.mesh.Nz); % Salinity state (g/kg)
	E=PickupData(:,:,(1)+6*mit.mesh.Nz); % free surface state (m)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Open new cells as necessary
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% find indices of locations where ice shelf retreated
	fname=sprintf('hFacC.save.%010i.data',niter);
	disp(['  reading ocean hFacC file ' fname]);
	hFacC=binread(fname,4,[mit.mesh.Nx, mit.mesh.Ny, mit.mesh.Nz]); % hFacC (m)
	hCol=sum(hFacC,3); % water column height (m)
   [iw jw]=find(hCol>0); % horizontal indices where there is water
	[im jm] = find(newdraft>draft & newdraft>mit.mesh.zp(end)); % horizontal indices where there is melt

	disp(['  found ' num2str(numel(im)) ' melt cells']);
	disp(['   - max diff draft = ' num2str(max((newdraft(:)-draft(:))))]);
	disp(['   - min diff draft = ' num2str(min((newdraft(:)-draft(:))))]);

	%Extrapolate T/S to locations where ice shelf retreated
	for i=1:length(im)
		% first try vertical extrapolation
		kw=find(hFacC(im(i),jm(i),:)); % the vertical index where there is water
		if ~isempty(kw)
			%% do we need to open a new cell?
			%k_old=find(draft(im(i),jm(i))<=mit.mesh.zp(1:end),1,'last');
			%k_new=find(newdraft(im(i),jm(i))<=mit.mesh.zp(1:end),1,'last');
			%% new ocean cell where draft falls in different vertical cell
			%if k_new~=k_old
			S(im(i),jm(i),1:min(kw)) = S(im(i),jm(i),min(kw));
			T(im(i),jm(i),1:min(kw)) = T(im(i),jm(i),min(kw));
			%end
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
	disp('  writing updated ocean files');
	% update restart files
   binwrite(draft_file,newdraft,8);
	binwrite(uvel_file ,U,8);
	binwrite(vvel_file ,V,8);
	binwrite(theta_file,T,8);
	binwrite(salt_file ,S,8);
	binwrite(etan_file ,E,8);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Update the data files
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% the cal start time, model start time, and obcs all stay the same. the niter restarts at 0
	% but we update the basetime so that the model knows where we are
	newline=['  startTime=' num2str(modeltime) ','];
   command=['sed "s/.*startTime.*/' newline '/" data > data.temp; mv data.temp data'];
	system(command);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Run MITgcm
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('  running MITgcm')
	tic
	system(['mpirun -np ' int2str(nprocs) ' ./mitgcmuv > out 2> err']);
	toc
	disp('  done MITgcm')
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Move files to time-corrected niter suffix
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('  saving output files to time-corrected niter')
	modeltime = coupled_basetime+(n+1)*mit.timestepping.coupledTimeStep; % updated modeltime
	niter     = niter0 + (n+1)*mit.timestepping.coupledTimeStep/mit.inputdata.PARM{3}.deltaT; % the updated niter
	movefile(sprintf('pickup.%010i.data',modelIterEnd), sprintf('pickup.save.%010i.data',niter)); % pickup.data
	movefile(sprintf('pickup.%010i.meta',modelIterEnd), sprintf('pickup.save.%010i.meta',niter)); % pickup.meta
	movefile(sprintf('SHICE_fwFluxtave.%010i.data',modelIterEnd), sprintf('SHICE_fwFluxtave.save.%010i.data',niter)); % SHICE_fwFluxtave.data
	movefile(sprintf('SHICE_fwFluxtave.%010i.meta',modelIterEnd), sprintf('SHICE_fwFluxtave.save.%010i.meta',niter)); % SHICE_fwFluxtave.meta
	% save the hFacC and draft files
	movefile('hFacC.data', sprintf('hFacC.save.%010i.data',niter)); % hFacC.data
	movefile('hFacC.meta', sprintf('hFacC.save.%010i.meta',niter)); % hFacC.meta
	movefile(draft_file, sprintf('draft.save.%010i.bin',niter)); % draft.bin

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Get melt from MITgcm
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('  reading melt from MITgcm')
	melt_fname=sprintf('SHICE_fwFluxtave.save.%010i.data',niter); % melt file
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

	%Get the change in ice shelf draft from ISSM
	deltaBase = md.results.TransientSolution(end).Base - md.geometry.base; % m

	%Save ISSM results
   fname = sprintf('issmDiag.%010i.mat', niter);
	disp(['  saving ISSM results to ' fname])
   results = md.results.TransientSolution(end);
	results.step = (n+1);
	results.time = (n+1)*mit.timestepping.coupledTimeStep;
	results.deltaBase = deltaBase;
   save(fname,'results');

	%Reset model
	disp('  reinitializing ISSM from results.TransientSolution')
   md.geometry.base             = md.results.TransientSolution(end).Base;
   md.geometry.surface          = md.results.TransientSolution(end).Surface;
   md.geometry.thickness        = md.geometry.surface-md.geometry.base;
   md.initialization.vx         = md.results.TransientSolution(end).Vx;
   md.initialization.vy         = md.results.TransientSolution(end).Vy;
   md.initialization.vel        = md.results.TransientSolution(end).Vel;
   md.mask.ocean_levelset       = md.results.TransientSolution(end).MaskOceanLevelset;
   clear md.results;

	%Clear execution folder in ISSM to avoid going over quota
	command = ['rm -r ' md.cluster.executionpath '/' md_prefix '*'];
	system(command);
end

% subfunctions
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
   disp(' __      __ _____ _______                    _     _  _____   _____ ____ __      __ ')
   disp('|   \  /   |_  __|__   __|                   \ \ / / |_   _|/ ____/ ____|   \  /   |')
   disp('| |\ \/ /| | | |    | | __ _  ___ _ __ ___    \   /    | |  \___ \\___ \| |\ \/ /| |')
   disp('| | \__/ | |_| |_   | |/ _` |/ __| `_ ` _ \   /   \   _| |_ ____| |___| | | \__/ | |')
   disp('|_|      |_|_____|  |_| (_| | (__| | | | | | /_/ \_\ |_____|_____/_____/|_|      |_|')
   disp('                       \___ |\___|_| |_| |_| ')
   disp('                        __/ |                ')
   disp('                       |___/                 ')
end % }}}
end
