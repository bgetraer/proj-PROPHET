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
	md.timestepping.final_time=mit.timestepping.deltaT_coupled./md.constants.yts; % how long to run ISSM for (y)
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
	relaxEndIter = mit.timestepping.relaxT/mit.timestepping.deltaT_relax;        % the final timestep iteration number of each relaxation run
	contEndIter  = mit.timestepping.contT/mit.timestepping.deltaT_cont;          % the final timestep iteration number of each relaxation run
	% }}}
	% initial parameters and fields {{{
	if mit.timestepping.ispickup==0
		if mit.timestepping.coupled_basetime~=mit.timestepping.startTime
			error('mit.timestepping.startTime is not equal to mit.timestepping.coupled_basetime, but md.timestepping.ispickup is FALSE');
		end
		disp(['Starting new coupled run from mit.timestepping.coupled_basetime = ' num2str(mit.timestepping.coupled_basetime)]);
		deltaBase=zeros(size(md.geometry.base)); % the change in ice shelf draft from ISSM (initialize to zero) (m)	
	elseif mit.timestepping.ispickup==1
		disp(['Picking up coupled run from mit.timestepping.startTime = ', num2str(mit.timestepping.startTime)]);
		fname = sprintf('issmDiag.%010i.mat', mit.timestepping.startTime); % issm results to load into md
		disp(['  reading ISSM results file ' fname]);
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
	% loop coupled steps {{{
	for n=-1:(mit.timestepping.nsteps-1)
	%for n=0:(mit.timestepping.nsteps-1)
		% update timekeeping {{{
		disp(['COUPLED STEP ' num2str(n+1) '/' num2str(mit.timestepping.nsteps)]);
		modeltime       = mit.timestepping.startTime + (n)*mit.timestepping.deltaT_coupled; % the start modeltime of this coupled step
		modeltime_cont  = modeltime + mit.timestepping.relaxT;                              % the modeltime after the relaxation run
		modeltime_next  = modeltime + mit.timestepping.deltaT_coupled;                      % the modeltime at the end of this coupled step
		elapse_y        = floor(modeltime/mit.timestepping.y2s); % elapsed years
		elapse_remsec   = mod(modeltime,mit.timestepping.y2s);   % elapsed remaining seconds
		disp(sprintf('modeltime: %010i',modeltime));
		disp(sprintf('elapsed time: %i yr, %s dd:mm:hh:ss',elapse_y,string(seconds(elapse_remsec),'dd:hh:mm:ss')));
		% }}}
		% ocean model {{{
		if modeltime>=mit.timestepping.coupled_basetime
			% update draft {{{
			fname=sprintf('draft.save.%010i.bin',modeltime);
			disp(['  reading previous draft file ' fname]);
			draft=binread(fname,8,[mit.mesh.Nx,mit.mesh.Ny]); % existing MITgcm draft in col,row matrix (m)

			disp('  interpolating updated draft change from ISSM');
			deltaDraft=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,deltaBase,mit.mesh.hXC(:),mit.mesh.hYC(:),'default',0); % change in ISSM draft (m)
			deltaDraft=permute(reshape(deltaDraft,mit.mesh.Ny,mit.mesh.Nx),[2,1]); % change in ISSM draft in col,row matrix (m)
			deltaDraft(mask_ice>0)=0; % mask the deltaDraft
			disp(['   - num deltaDraft cells = ' num2str(sum(deltaDraft(:)~=0))]);

			disp('  calculating newdraft');
			newdraft=draft+deltaDraft; % updated MITgcm draft in col,row matrix (m)
			disp(['   - max diff draft = ' num2str(max((newdraft(:)-draft(:))))]);
			disp(['   - min diff draft = ' num2str(min((newdraft(:)-draft(:))))]);

			disp('   applying ocean mask');
			mask_oce=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ocean_levelset,mit.mesh.hXC(:),mit.mesh.hYC(:),'default',-1); % -1 ocean, 1 grounded
			mask_oce=permute(reshape(mask_oce,mit.mesh.Ny,mit.mesh.Nx),[2,1]); % ISSM ocean mask in col,row matrix (m)
			newdraft(mask_oce>0 & mask_ice<0)=bathy(mask_oce>0 & mask_ice<0); % set all grounded ice to have a draft equal to the bathymetry (m)
			disp(['   - max diff draft = ' num2str(max((newdraft(:)-draft(:))))]);
			disp(['   - min diff draft = ' num2str(min((newdraft(:)-draft(:))))]);

			disp('   applying ice mask');
			newdraft(mask_ice>0)=0;                 % set all open ocean to have zero draft (m)
			disp(['   - max diff draft = ' num2str(max((newdraft(:)-draft(:))))]);
			disp(['   - min diff draft = ' num2str(min((newdraft(:)-draft(:))))]);

			disp('   applying OBCS mask');
			newdraft(:,1)=mit.geometry.draftOBS;  % set draft at bottom boundary (m)
			newdraft(1,:)=mit.geometry.draftOBW;  % set draft at left boundary (m)
			disp(['   - max diff draft = ' num2str(max((newdraft(:)-draft(:))))]);
			disp(['   - min diff draft = ' num2str(min((newdraft(:)-draft(:))))]);
			% }}}
			% read pickup file, open new cells, write updated init files {{{
			% read pickup file
			fname=sprintf('pickup.save.%010i.data',modeltime); % the data filename
			disp(['  reading ocean pickup file ' fname]);
         PickupData=binread(fname,8,[mit.mesh.Nx, mit.mesh.Ny, 6*mit.mesh.Nz+3]); % read the whole file
			U=PickupData(:,:,(1:mit.mesh.Nz)+0*mit.mesh.Nz); % x component of velocity (m/s)
			V=PickupData(:,:,(1:mit.mesh.Nz)+1*mit.mesh.Nz); % y component of velocity (m/s)
			T=PickupData(:,:,(1:mit.mesh.Nz)+2*mit.mesh.Nz); % Temperature state (deg C)
			S=PickupData(:,:,(1:mit.mesh.Nz)+3*mit.mesh.Nz); % Salinity state (g/kg)
			E=PickupData(:,:,(1)+6*mit.mesh.Nz); % free surface state (m)
			% open new cells as necessary 
			% find indices of locations where ice shelf retreated
			fname=sprintf('hFacC.save.%010i.data',modeltime);
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
				if numel(kw)>0
					S(im(i),jm(i),1:min(kw)) = S(im(i),jm(i),min(kw));
					T(im(i),jm(i),1:min(kw)) = T(im(i),jm(i),min(kw));
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
			% write updated MITgcm files
			disp('  writing updated ocean files');
			% update restart files
			binwrite(draft_file,newdraft,8);
			binwrite(uvel_file ,U,8);
			binwrite(vvel_file ,V,8);
			binwrite(theta_file,T,8);
			binwrite(salt_file ,S,8);
			binwrite(etan_file ,E,8);
			% }}}
			% update ./data files for fine deltaT relaxation run {{{
			disp('  setting runtime options for relaxation run');
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% I am defining these with startTime and nIter0 because I want to start from nIter0=0 but with
			% a modeltime of the correct calendar. The cal start time, and obcs all stay the same.

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% ./data
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			mit.inputdata.PARM{3} = struct();
			% structure information
			mit.inputdata.PARM{3}.header='PARM03';
			mit.inputdata.PARM{3}.description='Time stepping parameters';
			% Run Start and Duration
			mit.inputdata.PARM{3}.nIter0      = 0;                             % starting timestep iteration number
			mit.inputdata.PARM{3}.nEndIter    = relaxEndIter;                  % end timestep iteration number
			mit.inputdata.PARM{3}.deltaT      = mit.timestepping.deltaT_relax; % mitgcm deltaT (s)
			mit.inputdata.PARM{3}.startTime   = modeltime;                     % run start time for this integration (s)
			% Restart/Pickup Files
			mit.inputdata.PARM{3}.pChkptFreq  = modeltime_cont;                % permanent pickup checkpoint file write interval (s)
			mit.inputdata.PARM{3}.ChkptFreq   = 0;                             % temporary pickup checkpoint file write interval (s)
			% Frequency/Amount of Output
			mit.inputdata.PARM{3}.monitorFreq = modeltime_cont;       % interval to write monitor output - every coupled time step (s)
			mit.inputdata.PARM{3}.cAdjFreq        = -1;                       % frequency of convective adj. scheme
         mit.inputdata.PARM{3}.monitorSelect   = 1;                        % group of monitor variables to output
         mit.inputdata.PARM{3}.dumpInitAndLast = '.FALSE.';                % write out initial and last iteration model state

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% input/data.diagnostics.relaxation
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			mit.inputdata.DIAG{1}.N(1).frequency = 0;
			mit.inputdata.DIAG{1}.N(2).frequency = 0;
			mit.inputdata.DIAG{1}.N(3).frequency = 0;

			disp(mit.inputdata.PARM{3});

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% write data files
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			disp([' - Set runtime options in data file']);
			write_datafile('data', mit.inputdata.PARM, 'MODEL PARAMETERS');

			disp([' - Set runtime options in data.diagnostics file']);
			write_datafile('data.diagnostics', mit.inputdata.DIAG, 'DIAGNOSTICS RUNTIME PARAMETERS');


			% }}}
			% run MITgcm to end of relaxation time {{{
			disp('  running MITgcm relaxation period')
			tic
			system(['mpirun -np ' int2str(nprocs) ' ./mitgcmuv > out 2> err']);
			toc
			disp('  done MITgcm relaxation period')

			% check if bad solve
			[~, r]=system('grep " cg2d: Sum(rhs),rhsMax =                    NaN  0.00000000000000E+00" STDOUT.0000 | uniq -c');
			if ~isempty(r)
				error('MITgcm bad solve: NaN in STDOUT. Ending run!');
			end
			% }}}
			% update ./data files for coarse deltaT continuation run {{{
			disp('  setting runtime options for continuation run');
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% I am defining these with startTime and nIter0 because I want to start from nIter0=0 but with
			% a modeltime of the correct calendar. The cal start time, and obcs all stay the same.

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% ./data
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			mit.inputdata.PARM{3} = struct();
			% structure information
			mit.inputdata.PARM{3}.header='PARM03';
			mit.inputdata.PARM{3}.description='Time stepping parameters';
			% Run Start and Duration
			mit.inputdata.PARM{3}.nIter0      = 0;                            % starting timestep iteration number
			mit.inputdata.PARM{3}.nEndIter    = contEndIter;                  % end timestep iteration number
			mit.inputdata.PARM{3}.deltaT      = mit.timestepping.deltaT_cont; % mitgcm deltaT (s)
			mit.inputdata.PARM{3}.startTime   = modeltime_cont;               % run start time for this integration (s)
			% Restart/Pickup Files
			mit.inputdata.PARM{3}.pChkptFreq  = modeltime_next;       % permanent pickup checkpoint file write interval (s)
			mit.inputdata.PARM{3}.ChkptFreq   = 0;                            % temporary pickup checkpoint file write interval (s)
			mit.inputdata.PARM{3}.pickupSuff  = sprintf('%010i',relaxEndIter); % force run to use pickups and read files with this suffix
			% Frequency/Amount of Output
			mit.inputdata.PARM{3}.monitorFreq     = modeltime_next;           % interval to write monitor output - every coupled time step (s)
			mit.inputdata.PARM{3}.cAdjFreq        = -1;                       % frequency of convective adj. scheme                    
			mit.inputdata.PARM{3}.monitorSelect   = 1;                        % group of monitor variables to output
			mit.inputdata.PARM{3}.dumpInitAndLast = '.FALSE.';                % write out initial and last iteration model state 

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% input/data.diagnostics.relaxation
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			mit.inputdata.DIAG{1}.N(1).frequency = 0;
			mit.inputdata.DIAG{1}.N(2).frequency = 0;
			mit.inputdata.DIAG{1}.N(3).frequency = modeltime_next;

			disp(mit.inputdata.PARM{3});

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% write data file
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			disp([' - Set runtime options in data file']);
			write_datafile('data', mit.inputdata.PARM, 'MODEL PARAMETERS');

			disp([' - Set runtime options in data.diagnostics file']);
			write_datafile('data.diagnostics', mit.inputdata.DIAG, 'DIAGNOSTICS RUNTIME PARAMETERS');
			% }}}
			% run MITgcm until end of coupled time step {{{
			disp('  running MITgcm continuation')
			tic
			system(['mpirun -np ' int2str(nprocs) ' ./mitgcmuv > out 2> err']);
			toc
			disp('  done MITgcm continuation')

			% check if bad solve
			[~, r]=system('grep " cg2d: Sum(rhs),rhsMax =                    NaN  0.00000000000000E+00" STDOUT.0000 | uniq -c');
			if ~isempty(r)
				error('MITgcm bad solve: NaN in STDOUT. Ending run!');
			end
			% }}}
			% move files to modeltime suffix {{{
			disp('  saving output files to modeltime suffix')

			movefile(sprintf('pickup.%010i.data',contEndIter), sprintf('pickup.save.%010i.data',modeltime_next)); % pickup.data
			movefile(sprintf('pickup.%010i.meta',contEndIter), sprintf('pickup.save.%010i.meta',modeltime_next)); % pickup.meta
			movefile(sprintf('SHICE_fwFluxtave.%010i.data',contEndIter), sprintf('SHICE_fwFluxtave.save.%010i.data',modeltime_next)); % SHICE_fwFluxtave.data
			movefile(sprintf('SHICE_fwFluxtave.%010i.meta',contEndIter), sprintf('SHICE_fwFluxtave.save.%010i.meta',modeltime_next)); % SHICE_fwFluxtave.meta
			% save the hFacC and draft files
			movefile('hFacC.data', sprintf('hFacC.save.%010i.data',modeltime_next)); % hFacC.data
			movefile('hFacC.meta', sprintf('hFacC.save.%010i.meta',modeltime_next)); % hFacC.meta
			movefile(draft_file,   sprintf('draft.save.%010i.bin', modeltime_next)); % draft.bin
			% delete the relaxation pickup file
			delete(sprintf('pickup.%010i.data',relaxEndIter));
			delete(sprintf('pickup.%010i.meta',relaxEndIter));
			% }}}
		end
		% }}}
		% ice model {{{
		% get melt from MITgcm {{{
		melt_fname=sprintf('SHICE_fwFluxtave.save.%010i.data',modeltime_next); % melt file
		disp(['  reading melt from MITgcm file: ' melt_fname])
		meltq_mitgcm = binread(melt_fname,4,[mit.mesh.Nx, mit.mesh.Ny]); % melt flux at cell centers (kg/m^2/s)
		meltq_mitgcm=permute(meltq_mitgcm,[2,1]);  % put in ROW COL order (kg/m^2/s)
		meltq_issm=InterpFromGridToMesh(mit.mesh.xc(:),mit.mesh.yc(:),meltq_mitgcm,md.mesh.x,md.mesh.y,0); % melt flux at vertices (kg/m^2/s)

		%Set basal melting rate fields
		md.basalforcings.floatingice_melting_rate=-meltq_issm*md.constants.yts/md.materials.rho_ice; % melt rate at element vertices (m/yr)
		% }}}
		% run ISSM {{{
		disp('  running ISSM')
		%Solve
		md.miscellaneous.name = sprintf('%s%010i',md_prefix,modeltime_next);
		md=solve(md,'transient');
		disp('  done ISSM')
		% }}}
		% save ISSM results {{{
		% get the change in ice shelf draft from ISSM
		deltaBase = md.results.TransientSolution(end).Base - md.geometry.base; % m
		% save to file
		fname = sprintf('issmDiag.%010i.mat', modeltime_next);
		disp(['  saving ISSM results to ' fname])
		results = md.results.TransientSolution(end);
		results.step = (n+1);
		results.time = (n+1)*mit.timestepping.deltaT_coupled;
		results.deltaBase = deltaBase;
		save(fname,'results');
		% }}}
		% reinitialize ISSM from results {{{
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
		% }}}
		% }}}
	end
	% }}}
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
function write_datafile(fname,C,head) % {{{
	% WRITE_DATAFILE writes structures in C to fname {{{
	% INPUT: fname   string file to write
	%        C       cell array of structures (P) to write
	%        head    string of file header
	%
	% OUTPUT: writes to file with the following form:
	%
	% # head
	% # C{1}.description
	%  &C{1}.header
	%  C{1}.field1=value1,
	%  C{1}.field2=value2,
	%  ...
	%  &
	%
	% # C{2}.description
	%  &C{2}.header
	%  C{2}.field1=value1,
	%  C{2}.field2=value2,
	%  ...
	%  &
	%  ...
	% }}}
	disp(['    writing namelist file to ' fname])
	fileID = fopen(fname,'w');
	fprintf(fileID,'# %s\n',head); % write descriptive file header

	% loop through structures in C
	for i=1:length(C)
		% write a descriptive comment if it exists {{{
		if isfield(C{i},'description')
			if iscell(C{i}.description)
				fprintf(fileID,'# %s\n',C{i}.description{:}); % write multi-line description
			else
				fprintf(fileID,'# %s\n',C{i}.description); % write description
			end
			C{i}=rmfield(C{i},'description'); % do not write as field
		end % }}}
		% write PARM header {{{
		fprintf(fileID,' &%s\n',C{i}.header); % write header
		C{i}=rmfield(C{i},'header'); % do not write as field
		% }}}
		% write parameter fields and end section {{{
		if isfield(C{i},'N') % if diagnostic fields
			writediagfields(fileID,C{i}.N);
			C{i}=rmfield(C{i},'N'); % done with diag fields
		end
		writefields(fileID,C{i}); % write non-diagnostic fields
		fprintf(fileID,' &\n\n'); % end section
		% }}}
	end
	fclose(fileID);
end % }}}
function writefields(fileID,P) % {{{
	% WRITEFIELDS writes each field in P to fileID as a new line
	%  Each line takes the form ' fieldname=value,'
	fields=fieldnames(P);
	for i=1:length(fields)
		val=getfield(P,fields{i});
		fprintf(fileID,'  %s=%s,\n',fields{i},num2str(val));
	end
end  % }}}
function writediagfields(fileID,N) % {{{
	% WRITEDIAGFIELDS writes diagnostic fields for each output stream
	% INPUT: N is a struct. array with an element for each output stream n
	%  Each line takes the form fieldname(n)=value except for
	%  'fields' with take the form fields(1:length,n)='fields{1}', 'fields{2}', ...

	% loop over each output stream
	for n=1:length(N)
		subfields=fieldnames(N(n));
		% loop over each subfield
		for i=1:length(subfields)
			switch subfields{i}
				case 'fields'
					LHS=[subfields{i} '(1:' num2str(length(N(n).fields)) ',' num2str(n) ')'];
					dfields=getfield(N(n),subfields{i});
					dfields=strcat('''',dfields,''', ');
					RHS=strcat(dfields{:});
					fprintf(fileID,'  %s=%s\n',LHS,RHS); % write to file
				case 'levels'
					error('diag. levels not supported');
				otherwise
					LHS=[subfields{i} '(' num2str(n) ')'];
					RHS=num2str(getfield(N(n),subfields{i}));
					fprintf(fileID,'  %s=%s,\n',LHS,RHS); % write to file
			end
		end
		fprintf(fileID,'\n'); % line break
	end
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
