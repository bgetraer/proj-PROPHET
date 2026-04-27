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
	bathy_ref_file = 'bathy_ref.bin';
	bathy_file = 'bathy.bin';
	uvel_file  = 'uvel.bin';
	vvel_file  = 'vvel.bin';
	theta_file = 'theta.bin';
	salt_file  = 'salt.bin';
	etan_file  = 'etan.bin';

	% Model execution parameters
	nprocs=mit.build.SZ.nPx*mit.build.SZ.nPy; % number of processors for MITgcm and ISSM
	md.cluster=generic('name',oshostname(),'np',min(nprocs,70)); % set number of processors for ISSM.
	md.timestepping.final_time=mit.timestepping.deltaT_coupled./md.constants.yts; % how long to run ISSM for (y)
	md_prefix = 'runcouple'; % the ISSM model name prefix for execution files

	% Static fields for opening and closing draft cells
	bathy_ref=binread(bathy_ref_file,8,[mit.mesh.Nx,mit.mesh.Ny]); % the MITgcm bathymetry in col,row matrix (m)
	mask_ice=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ice_levelset,mit.mesh.hXC(:),mit.mesh.hYC(:),'default',1); % ISSM ice mask (m)
	mask_ice=permute(reshape(mask_ice,mit.mesh.Ny,mit.mesh.Nx),[2,1]); % ISSM ice mask in col,row matrix (m)

	% Timekeeping parameters
	% n is the coupled step number we are STARTING FROM, from 0:nsteps-1
	% niter is the MITgcm step number
	% modeltime is the MITgcm modeltime starting from the calendar start date (2010)
	% coupled_basetime is the MITgcm modeltime that we start the coupling at (2013)
	% basetime is the MITgcm modeltime that the current model starts at 
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
	%for n=-1:(mit.timestepping.nsteps-1)
	for n=0:(mit.timestepping.nsteps-1)
		% update timekeeping {{{
		disp(['COUPLED STEP ' num2str(n+1) '/' num2str(mit.timestepping.nsteps)]);
		modeltime       = mit.timestepping.startTime + (n)*mit.timestepping.deltaT_coupled; % the start modeltime of this coupled step
		modeltime_next  = modeltime + mit.timestepping.deltaT_coupled;                      % the modeltime at the end of this coupled step
		elapse_y        = floor(modeltime/mit.timestepping.y2s); % elapsed years
		elapse_remsec   = mod(modeltime,mit.timestepping.y2s);   % elapsed remaining seconds
		disp(sprintf('modeltime: %010i',modeltime));
		disp(sprintf('elapsed time: %i yr, %s dd:mm:hh:ss',elapse_y,string(seconds(elapse_remsec),'dd:hh:mm:ss')));
		% }}}
		% ocean model {{{
		%if modeltime>=mit.timestepping.coupled_basetime
		% update draft {{{
		fname=sprintf('draft.save.%010i.bin',modeltime);
		disp(['  reading previous draft file ' fname]);
		olddraft=binread(fname,8,[mit.mesh.Nx,mit.mesh.Ny]); % existing MITgcm draft in col,row matrix (m)

		% get the cavity height from ISSM
		disp('  interpolating updated cavity height from ISSM');
		cavityH_issm = md.geometry.base - md.geometry.bed; % cavity height on ISSM mesh (m)
		cavityH = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,cavityH_issm,mit.mesh.hXC(:),mit.mesh.hYC(:),'default',0); % cavity height from ISSM on MITgcm grid (m)
		cavityH = permute(reshape(cavityH,mit.mesh.Ny,mit.mesh.Nx),[2,1]); % cavity height in col,row matrix (m)

		disp('  calculating newdraft');
		newdraft=bathy_ref+cavityH; % set the MITgcm draft to have the same cavity height as ISSM (m)
		disp(['   - max diff draft = ' num2str(max((newdraft(:)-olddraft(:))))]);
		disp(['   - min diff draft = ' num2str(min((newdraft(:)-olddraft(:))))]);

		disp('   applying ocean mask');
		mask_oce=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ocean_levelset,mit.mesh.hXC(:),mit.mesh.hYC(:),'default',-1); % -1 ocean, 1 grounded
		mask_oce=permute(reshape(mask_oce,mit.mesh.Ny,mit.mesh.Nx),[2,1]); % ISSM ocean mask in col,row matrix (m)
		newdraft(mask_oce>0 & mask_ice<0)=bathy_ref(mask_oce>0 & mask_ice<0); % set all grounded ice to have a draft equal to the bathymetry (m)
		disp(['   - max diff draft = ' num2str(max((newdraft(:)-olddraft(:))))]);
		disp(['   - min diff draft = ' num2str(min((newdraft(:)-olddraft(:))))]);

		disp('   applying ice mask');
		newdraft(mask_ice>0)=0; % set all open ocean to have zero draft (m)
		disp(['   - max diff draft = ' num2str(max((newdraft(:)-olddraft(:))))]);
		disp(['   - min diff draft = ' num2str(min((newdraft(:)-olddraft(:))))]);

		disp('   applying OBCS mask');
		newdraft(:,1)=mit.geometry.draftOBS;  % set draft at bottom boundary (m)
		newdraft(1,:)=mit.geometry.draftOBW;  % set draft at left boundary (m)
		disp(['   - max diff draft = ' num2str(max((newdraft(:)-olddraft(:))))]);
		disp(['   - min diff draft = ' num2str(min((newdraft(:)-olddraft(:))))]);
		% }}}
		% read pickup file, open new cells, write updated init files {{{
		% read pickup file
		fname=sprintf('pickup.save.%010i.data',modeltime); % the data filename
		disp(['  reading ocean pickup file ' fname]);
		if n=1
			PickupData=binread(fname,8,[mit.mesh.Nx, mit.mesh.Ny, 6*mit.mesh.Nz+3]); % read the whole file
			U=PickupData(:,:,(1:mit.mesh.Nz)+0*mit.mesh.Nz); % x component of velocity (m/s)
			V=PickupData(:,:,(1:mit.mesh.Nz)+1*mit.mesh.Nz); % y component of velocity (m/s)
			T=PickupData(:,:,(1:mit.mesh.Nz)+2*mit.mesh.Nz); % Temperature state (deg C)
			S=PickupData(:,:,(1:mit.mesh.Nz)+3*mit.mesh.Nz); % Salinity state (g/kg)
			E=PickupData(:,:,(1)+6*mit.mesh.Nz); % free surface state (m)
		end

		U_old=U;
		V_old=V;
		T_old=T;
		S_old=S;
		E_old=E;

		% NEW CELL FILLING
		% 1) identify newly opened cells based on change in draft
		%   1.1) identify new hFacC
		%   1.2) find connectivity between all cells
		%   1.3) close "lakes" which are not connected to the main cavity (assumes that the ocean is ALWAYS the largest connected region)
		% 2) find the horizontal index of columns with new cell
		% 3) loop over each column with new cell
		%   3.1) if column has open cell, extrapolate upwards
		%   3.2) if column is new, extrapolate as average of neighboring columns
		%   3.3) repeat loop until all cells are filled

		hFacC_new = compute_hfac(mit.mesh.zp, bathy_ref, newdraft, mit.inputdata.PARM{1}.hFacMin); % only return hfacC
		fname=sprintf('hFacC.save.%010i.data',modeltime);
		hFacC_old = binread(fname,4,[mit.mesh.Nx, mit.mesh.Ny, mit.mesh.Nz]);

		% 1.2) find connectivity between all cells
		BW_new = hFacC_new>0; % binary 3D mask of all open cells
		conn = 6; % 3D connectivity kernel follows the MITgcm connectivity: adjacent faces are connected, edges and corners are not.
		CC = bwconncomp(BW_new,conn); % get the connected components structure
		[~,ind_oce] = max(cellfun(@numel,CC.PixelIdxList)); % index of the largest connected volume
		ind_lakes = (1:numel(CC.PixelIdxList))~=ind_oce; % indices of the disconnected volumes
		BW_new(:) = 0; % reset all connectivites to zero
		BW_new(CC.PixelIdxList{ind_oce}) = 1; % open only cells in the largest connected region
		% display size of largest connected region, and any disconnected regions to be closed
		n_all = sum(cellfun(@numel,CC.PixelIdxList));
		n_oce = cellfun(@numel,CC.PixelIdxList(ind_oce));
		fprintf('n connected cells: %i\n',n_oce);
		fprintf('n disconnected cells: %i\n',n_all-n_oce);

		% 1.3) close "lakes" which are not connected to the main cavity (assumes that the ocean is ALWAYS the largest connected region)
		hFacC_new(~BW_new)=0; % CLOSE all cells disconnected from the largest region
		bathy = bathy_ref; % reset the runtime bathymetry to the reference (m)
		bathy(~sum(BW_new,3)) = 0; % close all cells which are not connected to the ocean region by setting runtime bathy to 0 (m)

		% 2) find the horizontal index of columns with new cell
		% define the cells which are opening
		cell_open = (hFacC_new & ~hFacC_old);

		% open new cells as necessary 
		% make sure that S and T match the old hFacC
		S(~hFacC_old) = NaN;
		T(~hFacC_old) = NaN;

		% linear index
		ind_w = find((nansum(S,3))~=0); % horizontal indices where there is water
		ind_o = find(sum(cell_open,3)); % horizontal indices where there is a new cell opening
		% columns with open cells to open vertically
		ind_o_fillup = intersect(ind_o,ind_w); % linear index of new cells with existing water column - fill upwards
		ind_o_fillhor = setdiff(ind_o,ind_w);  % linear index of new water columns - fill horizontally


		disp(['  found ' num2str(numel(ind_o)) ' columns to open']);
		disp(['   - max diff draft = ' num2str(max((newdraft(:)-olddraft(:))))]);
		disp(['   - min diff draft = ' num2str(min((newdraft(:)-olddraft(:))))]);

		%Extrapolate T/S to locations where ice shelf retreated
		% 3) loop over each column with new cell
		% 3.1) if column has open cell, extrapolate upwards
		fprintf('EXPANDING %i EXISTING COLUMNS, VERTICAL EXTRAPOLATION\n',numel(ind_o_fillup));
		for i=1:length(ind_o_fillup)
			[io,jo] = ind2sub(size(E),ind_o_fillup(i)); % subscript indices for this open column
			kw_old=find(hFacC_old(io,jo,:)); % the old vertical indices where there is water
			kw_new=find(hFacC_new(io,jo,:)); % the new vertical indices where there is water
			S(io,jo,min(kw_new):min(kw_old)) = S(io,jo,min(kw_old)); % Salinity, fill upwards
			T(io,jo,min(kw_new):min(kw_old)) = T(io,jo,min(kw_old)); % Theta, fill upwards
		end
		% 3.2) if column is new, extrapolate as average of neighboring columns
		% 3.3) repeat loop until all cells are filled
		%
		% ** known issues: the weighted averaging is dependent on the order of filling:
		%    some cells will not "see" the newly opened neighbors if they have previous
		%    open neighbors. For now this is left as is.
		fprintf('OPENING %i NEW COLUMNS, WEIGHTED HORIZONTAL INTERPOLATION\n',numel(ind_o_fillhor));
		open_counter = zeros(size(ind_o_fillhor));
		i_adj = [0,-1,1,0]; % column subscript of adjacent column cells, centered at (0,0)
		j_adj = [-1,0,0,1]; % row subscript of adjacent column cells, centered at (0,0)
		while_iter = 0; % number of while loop iterations
		while any(open_counter==0)
			for i=1:length(ind_o_fillhor)
				if open_counter(i)==0
					% indexing
					[io,jo] = ind2sub(size(E),ind_o_fillhor(i)); % subscript indices for this open column
					ind_adj = sub2ind(size(cell_open),...
						io+repmat(i_adj,mit.mesh.Nz,1),... % column subscript of adjacent column cells
						jo+repmat(j_adj,mit.mesh.Nz,1),... % row subscript of adjacent column cells
						repmat(1:mit.mesh.Nz,4,1)'); % Nz by 4 matrix of adjacent linear cell indices
					ind_adjH = sub2ind(size(E),io+i_adj,jo+j_adj); % subscript horizontal indices

					% data from adjacent columns
					salt_profile = S(ind_adj); % salinity profiles of adjacent columns
					temp_profile = T(ind_adj); % salinity profiles of adjacent columns
					hfacc_profile = hFacC_new(ind_adj); % hFacC profiles of adjacent columns

					% if any connected cells in the adjacent columns are open, interpolate from those columns!
					kw_new=find(hFacC_new(io,jo,:)); % the new vertical indices where there is water
					if any(~isnan(salt_profile(kw_new,:)),'all')
						open_adj = any(~isnan(salt_profile),1); % index of open adjacent columns
						% weight S and T by hFacC for each column to generate a reference column
						salt_profile = nansum(salt_profile(:,open_adj).*hfacc_profile(:,open_adj),2)...
							./nansum(hfacc_profile(:,open_adj),2); % weighted mean of adjacent columns
						temp_profile = nansum(temp_profile(:,open_adj).*hfacc_profile(:,open_adj),2)...
							./nansum(hfacc_profile(:,open_adj),2); % weighted mean of adjacent columns

						% fill new cells from the reference column, filling any nan holes and extrapolating values up and down as needed.
						S(io,jo,kw_new) = fillmissing(salt_profile(kw_new),'linear','EndValues','nearest'); % set salinity for new ocean column
						T(io,jo,kw_new) = fillmissing(temp_profile(kw_new),'linear','EndValues','nearest'); % set temperature for new ocean column
						if any(E(ind_adjH)~=0)
							E(io,jo) = mean(E(ind_adjH(E(ind_adjH)~=0))); % set the free surface to the average of adjacent non-zero values
						end
						% if successfully filled, update the register to reflect opened column
						% otherwise, let it try again after filling the rest of the columns
						if all(~isnan(S(io,jo,kw_new)))
							open_counter(i) = 1; 
						end
					end
				end
			end
			while_iter = while_iter+1;
			if while_iter==25
				disp('************************************************************************************');
				fprintf('WHILE LOOP NOT CONVERGING, ITER=%i\n',while_iter);
				disp('   saving workspace to runcouple/workspace/runcouple_workspace.mat');
				save('/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET/runcouple/workspace/runcouple_workspace.mat');
				disp('Ending run now!');
				disp('************************************************************************************');
				return;
			end
		end

		% Force all closed cells/columns have 0 or NaN values
		U(~hFacC_new)=0;
		V(~hFacC_new)=0;
		T(~hFacC_new)=NaN;
		S(~hFacC_new)=NaN;
		E(~sum(hFacC_new,3))=0;

		% Check all open cells do not have NaN values
		if any(hFacC_new & (isnan(S) | isnan(T)),'all')
			disp('************************************************************************************');
			fprintf('INTERPOLATION ERROR: NaN found in %i open cells!!\n',sum(hFacC_new & (isnan(S) | isnan(T)),'all'));
			disp('   saving workspace to runcouple/workspace/runcouple_workspace.mat');
			save('/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET/runcouple/workspace/runcouple_workspace.mat');
			disp('Ending run now!');
			disp('************************************************************************************');
			return;
		end
		% Check if any newly opened cells have exactly zero Eta
		if any(E(ind_o_fillhor)==0)
			disp('INTERPOLATION WARNING: Eta is exactly zero in newly open column');
		end

		%disp('   saving workspace to runcouple/workspace/runcouple_workspace.mat');
		%save('/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET/runcouple/workspace/runcouple_workspace.mat');
		%return;

		% write updated MITgcm files
		disp('  writing updated ocean files');
		% update restart files
		binwrite(draft_file,newdraft,8);
		binwrite(bathy_file,bathy,8);
		binwrite(uvel_file ,U,8);
		binwrite(vvel_file ,V,8);
		binwrite(theta_file,T,8);
		binwrite(salt_file ,S,8);
		binwrite(etan_file ,E,8);
		% }}}
		% update ./data file {{{
		disp('  setting runtime options');
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% I am defining these with startTime and nIter0 because I want to start from nIter0=0 but with
		% a modeltime of the correct calendar. The cal start time, and obcs all stay the same.

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% ./data
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Update the startTime to the current modeltime
		mit.inputdata.PARM{3}.startTime   = modeltime; % run start time for this integration (s)
		disp(mit.inputdata.PARM{3});

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% input/data.diagnostics
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		mit.inputdata.DIAG{1}.N(3).frequency = modeltime_next;

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% write data files
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		disp([' - Set runtime options in data file']);
		write_datafile('data', mit.inputdata.PARM, 'MODEL PARAMETERS');

		disp([' - Set runtime options in data.diagnostics file']);
		write_datafile('data.diagnostics', mit.inputdata.DIAG, 'DIAGNOSTICS RUNTIME PARAMETERS');
		% }}}
		% run MITgcm until end of coupled time step {{{
		%return
		disp('  SKIPPING MITgcm!')
		%disp('  running MITgcm')
		%tic
		%system(['mpirun -np ' int2str(nprocs) ' ./mitgcmuv > out 2> err']);
		%toc
		%disp('  done MITgcm')

		% check if bad solve
		[~, r]=system('grep " cg2d: Sum(rhs),rhsMax =                    NaN  0.00000000000000E+00" STDOUT.0000 | uniq -c');
		if ~isempty(r)
			error('MITgcm bad solve: NaN in STDOUT. Ending run!');
		end
		% }}}
		% move files to modeltime suffix {{{
		disp('  saving output files to modeltime suffix')

		%copyfile(sprintf('pickup.%010i.data',modeltime), sprintf('pickup.save.%010i.data',modeltime_next)); % pickup.data
		%copyfile(sprintf('pickup.%010i.meta',modeltime), sprintf('pickup.save.%010i.meta',modeltime_next)); % pickup.meta
		%copyfile(sprintf('SHICE_fwFluxtave.%010i.data',modeltime), sprintf('SHICE_fwFluxtave.save.%010i.data',modeltime_next)); % SHICE_fwFluxtave.data
		%movefile(sprintf('SHICE_fwFluxtave.%010i.meta',modeltime), sprintf('SHICE_fwFluxtave.save.%010i.meta',modeltime_next)); % SHICE_fwFluxtave.meta
		% save the hFacC and draft files
		binwrite(sprintf('hFacC.save.%010i.data',modeltime_next),hFacC_new,4)
		movefile(draft_file,   sprintf('draft.save.%010i.bin', modeltime_next)); % draft.bin
		% }}}
		%end
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
		results.time = modeltime_next;
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
function varargout = compute_hfac(zp, bathy, draft, hFacMin) % {{{
	% compute_hfac - Calculates fractional layer thicknesses (hFacC) based on bathymetry and ice draft.
	% This function calculates the fractional layer thicknesses (hFacC) for a given
	% vertical grid (zp), bathymetry (bathy), ice draft (draft), and minimum
	% fractional cell thickness (hFacMin).
	% NOTE: This function assumes that bathymetry and draft matrices are in column,row format (y,x) as
	% taken from MITgcm files.
	%
	% Inputs:
	%   zp      - Vector of layer face locations in z dimension (m)
	%   bathy   - Matrix of bathymetry depths in col,row format (m)
	%   draft   - Matrix of ice draft depths in col,row format (m)
	%   hFacMin - Minimum non-dimensional fraction of cell thickness allowed
	%
	% Outputs:
	%   hfacC   - Matrix of fractional layer thicknesses at cell centers
	%   hfacW   - (Optional) Matrix of fractional layer thicknesses at west faces
	%   hfacS   - (Optional) Matrix of fractional layer thicknesses at south faces
	%
	% Useage:
	%   hfacC = compute_hfac(zp, bathy, draft, hFacMin); % only return hfacC
	%   [hfacC, hfacW, hfacS] = compute_hfac(zp, bathy, draft, hFacMin); % return all hfacs

	% check that bathy and draft are same size
	assert(all(size(bathy)==size(draft)),'bathy and draft matrices are not consistent sizes');
	% put bathy and draft into row,col order

	% get hFacC, using the draft and bathy at cell centers
	hFacC = get_hfac(zp, bathy, draft, hFacMin);
	varargout{1} = hFacC;
	if nargout>1
		% get values at the U point (western edge) of the C-grid
		bathyW = max(bathy, bathy([end,1:end-1],:)); % bathymetry at western edge of cell (m)
		draftW = min(draft, draft([end,1:end-1],:)); % draft at western edge of cell (m)
		hFacW = get_hfac(zp, bathyW, draftW, hFacMin);
		varargout{2} = hFacW;
	end
	if nargout>2
		% get values at the V point (southern edge) of the C-grid
		bathyS = max(bathy, bathy(:,[end,1:end-1])); % bathymetry at southern edge of cell (m)
		draftS = min(draft, draft(:,[end,1:end-1])); % draft at southern edge of cell (m)
		hFacS = get_hfac(zp, bathyS, draftS, hFacMin);
		varargout{3} = hFacS;
	end

	% subfunction for calculating the hFac
function hFac = get_hfac(zp, bathy, draft, hFacMin)
	% dimensions
	Nx = size(bathy,1);
	Ny = size(bathy,2);
	Nh = numel(bathy);
	Nz = numel(zp)-1;

	% make bathy and draft row vectors
	bathy_row = reshape(bathy,1,Nh);
	draft_row = reshape(draft,1,Nh);

	% make z edges into column vectors
	delzF = abs(diff(zp));
	delzF = delzF(:);
	zp_upper = reshape(zp(1:end-1),Nz,1);
	zp_lower = reshape(zp(2:end),Nz,1);
	% initialize to all closed
	hFac = zeros(Nz,Nh);

	% lower edge is at or above bathymety; upper edge is at or below ice draft
	ind = (zp_lower>=bathy_row) & (zp_upper<=draft_row); % index of full hFac
	hFac(ind) = 1.0; % set full hFacC

	% Find partial cells due to bathymetry alone
	ind = find((zp_lower<bathy_row) & (zp_upper>bathy_row) & (zp_upper<=draft_row));
	[i1,i2]=ind2sub(size(hFac),ind);
	hFac(ind) = (zp_upper(i1) - bathy_row(i2)')./delzF(i1);

	% Find partial cells due to ice shelf draft alone
	ind = find((zp_lower>=bathy_row) & (zp_lower<draft_row) & (zp_upper>draft_row));
	[i1,i2]=ind2sub(size(hFac),ind);
	hFac(ind) = (draft_row(i2)' - zp_lower(i1))./delzF(i1);

	% Find partial cells which are intersected by both
	ind = find((zp_lower<bathy_row) & (zp_lower<draft_row) & (zp_upper>bathy_row) & (zp_upper>draft_row));
	[i1,i2]=ind2sub(size(hFac),ind);
	hFac_bathy = (zp_upper(i1) - bathy_row(i2)')./delzF(i1);
	hFac_draft = (draft_row(i2)' - zp_lower(i1))./delzF(i1);
	% apply hFac to hFac_bathy first
	ind_bathy = round(hFac_bathy-(hFacMin/2),10) < 0;
	hFac_bathy(ind_bathy) = 0;
	ind_bathy = round(hFac_bathy-hFacMin,10) < 0 & round(hFac_bathy-(hFacMin/2),10) >= 0;
	hFac_bathy(ind_bathy) = hFacMin;
	% calculate the hFac
	hFac(ind) = hFac_bathy + hFac_draft - 1;

	% Apply hFacMin everywhere
	ind = round(hFac-(hFacMin/2),10) < 0;
	hFac(ind) = 0;
	ind = round(hFac-hFacMin,10) < 0 & ~ind;
	hFac(ind) = hFacMin;

	% reshape hfac
	hFac = reshape(hFac',Nx,Ny,Nz);
end
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
