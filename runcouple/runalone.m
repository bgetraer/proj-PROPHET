function runalone(mdfile,mitfile)
%RUNCOUPLE is intended as a standalone script to run a coupled ISSM-MITGCM model with MCC compilation/
%The inputs are mdfile which points to the location of the ISSM model file, and mitfile, which points
%to the location of the mit model file.
%These variables are declared explicitly and named in this script as the executable needs to expect them in order to load them.
%RUNCOUPLE loads the existing environment variables, loops through the time steps calling the models runs,
%and saves the output. Is is assumed that you are already located within the mitgcm "run" directory. 
dispMITxISSM();
disp('************************************************************************************');
disp('*   - beginning RUNALONE mcc deployable');
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
%% the appropriate number of MITgcm timesteps, and those results are loaded for the ice model.

% File names
draft_file = 'draft.bin';
bathy_file = 'bathy.bin';
uvel_file  = 'uvel.bin';
vvel_file  = 'uvel.bin';
theta_file = 'theta.bin';
salt_file  = 'salt.bin';
etan_file  = 'etan.bin';

coupled_basetime=mit.inputdata.PARM{3}.baseTime; % modeltime that we start coupling at
modelIterEnd=mit.timestepping.coupledTimeStep/mit.inputdata.PARM{3}.deltaT; % the final timestep number of each model run
for n=0:(mit.timestepping.nsteps-1);
	display(['COUPLED STEP ' num2str(n+1) '/' num2str(mit.timestepping.nsteps)]);
	modeltime = coupled_basetime+(n)*mit.timestepping.coupledTimeStep; % the current modeltime
   niter     = modeltime/mit.inputdata.PARM{3}.deltaT;                % the current niter
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Read ocean pickup file
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('  reading ocean pickup file');
   pickup_file=sprintf('pickup.%010i.data',niter); % the data filename
   PickupData=binread(pickup_file,8,[mit.mesh.Nx, mit.mesh.Ny, 6*mit.mesh.Nz+3]); % read the whole file
   U=PickupData(:,:,(1:mit.mesh.Nz)+0*mit.mesh.Nz); % x component of velocity (m/s)
   V=PickupData(:,:,(1:mit.mesh.Nz)+1*mit.mesh.Nz); % y component of velocity (m/s)
   T=PickupData(:,:,(1:mit.mesh.Nz)+2*mit.mesh.Nz); % Temperature state (deg C)
   S=PickupData(:,:,(1:mit.mesh.Nz)+3*mit.mesh.Nz); % Salinity state (g/kg)
   E=PickupData(:,:,(1)+7*mit.mesh.Nz); % free surface state (m)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Write updated MITgcm files
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp('  write updated ocean files');
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
   newline=['  baseTime=' num2str(modeltime) ','];
   command=['sed "s/.*baseTime.*/' newline '/" data > data.temp; mv data.temp data'];
   system(command);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Run MITgcm
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('  running MITgcm')
	tic
	system(['mpirun -np ' int2str(npMIT) ' ./mitgcmuv > out 2> err']);
	toc
	disp('  done MITgcm')
	modeltime = coupled_basetime+(n+1)*mit.timestepping.coupledTimeStep; % updated modeltime
   niter     = modeltime/mit.inputdata.PARM{3}.deltaT;                % updated niter
	% pickup
   source=sprintf('pickup.%010i.data',modelIterEnd); % where the pickup file is written
   destination=sprintf('pickup.%010i.data',niter);   % where we move the pickup file
   movefile(sprintf('pickup.%010i.data',modelIterEnd), sprintf('pickup.%010i.data',niter));
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
