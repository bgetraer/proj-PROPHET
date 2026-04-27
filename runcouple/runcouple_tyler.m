function runcouple_tyler(mdfile,mitfile)
	%RUNCOUPLE_TYLER
	disp('parse inputs');
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
	%Prepare runs (set all runtime variables in step-1)
	disp('parameters');
	% {{{ Parameters:
	%Set parameters
	Nx=mit.mesh.Nx;         %number of longitude cells
	Ny=mit.mesh.Ny;         %number of latitude cells
	Nz=mit.mesh.Nz;         %number of vertical cells
	nPx=mit.build.SZ.nPx;   %number of MITgcm processes to use in x direction
	nPy=mit.build.SZ.nPx;   %number of MITgcm processes to use in y direction
	xgOrigin=mit.mesh.xp(1);   %origin of longitude
	ygOrigin=mit.mesh.yp(1);   %origin of latitude
	dLong=mit.mesh.delxF; %longitude grid spacing
	dLat=mit.mesh.delyF;  %latitude grid spacing
	rho_ice=md.materials.rho_ice;    %density of ice
	y2s=mit.timestepping.y2s;   %years to seconds conversion factor
	c_int=mit.timestepping.coupledTimeStep;  %coupling interval in seconds (2-weeks)
	deltat_ocn=mit.inputdata.PARM{3}.deltaT; %ocean time step
	c_iter=sprintf('%010d',round(c_int/deltat_ocn)); %coupling iteration for ocean model

	%Set paths to present directory and directory where MTgcm is being run
	%	mitgcm_dir = '/nobackupp2/tpelle1/MITgcm/install/run690x280x70_highres/run_lowemission';
	%	pres_dir = '/nobackupp2/tpelle1/totten_mcc_highres_lowemission';

	%Set names of MITgcm parameter files
	draft_file = 'draft.bin';
	bed_file = 'bathy.bin';
	theta_file = 'theta.bin';
	salt_file  = 'salt.bin';
	uvel_file  = 'uvel.bin';
	vvel_file  = 'vvel.bin';
	etan_file  = 'etan.bin';

	%Set sizes of MITgcm files
	pickup_size = 350280000;  
	hFacC_size  = 28980000;
	fwflux_size = 420000;

	% Timekeeping parameters
   % n is the coupled step number we are STARTING FROM, from 0:nsteps-1
   % niter is the MITgcm step number
   % modeltime is the MITgcm modeltime starting from the calendar start date (2010)
   % coupled_basetime is the MITgcm modeltime that we start the coupling at (2013)
   % basetime is the MITgcm modeltime that the current model starts at
   coupled_basetime=mit.inputdata.PARM{3}.startTime; % modeltime that we start coupling at
   niter0=mit.inputdata.PARM{3}.nIter0; % starting niter
   modelIterEnd=mit.timestepping.coupledTimeStep/mit.inputdata.PARM{3}.deltaT; % the final timestep number of each model run

	%%Set how many years (yrs) you would like the coupling projection to run,
	%%then we set the number of iterations automatically assuming a 2-week timestep
	%yrs = 1; duration = yrs*12*2;

	%%Set dates for data.cal, loop over 20130101 to 20171215 (2 week time-step)
	%%Exf files span 1992 to 2017, so we must loop these years
	%startyr = 2013; endyr = 2100; dates=[];
	%for i=startyr:endyr  %years
	%	for ii=1:12       %months
	%		for iii=[1 15] %days (on the 1st and 15th of each month)
	%			yr = num2str(i);
	%			month = sprintf('%02d',ii);
	%			day = sprintf('%02d',iii);
	%			temp = [yr month day];
	%			dates = [dates; temp];
	%		end
	%	end
	%end
	%%Trim dates to correct duration
	%dates = repmat(dates,4,1); dates = dates(1:duration,:);

	%%Set dates of for naming files (2 week time-step), our projection starts at 2017
	%startyr = 2017; endyr = 2017+yrs-1; dates_name=[];
	%for i=startyr:endyr
	%	for ii=[1:12]
	%		for iii=[1 15] 
	%			yr = num2str(i);
	%			month = sprintf('%02d',ii);
	%			day = sprintf('%02d',iii);
	%			temp = [yr month day];
	%			dates_name = [dates_name; temp];
	%	   end
	%	end
	%end
	%}}}
	disp('mitgcm model');
	% {{{ MITgcmModel:
	%Create MITgcm grid in ISSM to use in data interpolation
	if exist('./mdm.mat')
		disp('  - loading mitgcm model');
		load mdm
	else
		disp('  - meshing mitgcm model');
		x = permute(mit.mesh.hXC,[2,1]);
		y = permute(mit.mesh.hYC,[2,1]);
		x=x(:);
		y=y(:);
		%Calculate the mesh connectivity
		index=[];
		%  C  D
		%  A  B
		for j=1:Ny-1,
			for i=1:Nx-1,
				A=(j-1)*Nx+i;
				B=(j-1)*Nx+i+1;
				C=j*Nx+i;
				D=j*Nx+i+1;
				index(end+1,:)=[A B C];
				index(end+1,:)=[C B D];
			end
		end
		%Set MITgcm model and fill mesh (basically a square mesh with triangular grid
		mdm = meshconvert(md,index,x,y);
		save mdm
	end
	%}}}
	disp('bathymetry and draft');
	% {{{ bathymetry and draft:
	%Get MITgcm bed and ISSM ice mask for to use in draft-interpolation
	bed = binread(bed_file,8,Nx,Ny,1); % matrix in ROW,COL order
	mask_ice = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ice_levelset,mdm.mesh.x,mdm.mesh.y,'default',-1); % vector in ROW,COL order
	% }}}

	%Run Coupled ISSM and MITgcm - start from initial state
	disp('run coupled models');
	% {{{ RunCoupledModels:
	%Begin looping over dates
	for n=-1:(mit.timestepping.nsteps-1)
		% Update timekeeping
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      display(['COUPLED STEP ' num2str(n+1) '/' num2str(mit.timestepping.nsteps)]);
      modeltime  = coupled_basetime+(n)*mit.timestepping.coupledTimeStep;                      % the current modeltime
      niter      = niter0 + (n)*mit.timestepping.coupledTimeStep/mit.inputdata.PARM{3}.deltaT; % the current niter
      niter_next = niter0 + (n+1)*mit.timestepping.coupledTimeStep/mit.inputdata.PARM{3}.deltaT; % the updated niter

		if n>=0
			% {{{Send draft and initial condition files to MITgcm:

			disp('update draft');
			%Get old draft, new draft from ISSM, and ice masks from ISSM - interpolate onto MITgcm grid
			fname=sprintf('draft.save.%010i.bin',niter);
         disp(['  reading previous draft file ' fname]);
			old_draft = binread(fname,8,Nx,Ny,1); % existing MITgcm draft in col,row matrix (m)
			new_draft = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.geometry.base,mdm.mesh.x,mdm.mesh.y,'default',999);
			mask_oce   = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ocean_levelset,mdm.mesh.x,mdm.mesh.y,'default',-1);

			%Take new draft from issm (>900) and sub it into mitgcm draft (<900)
			draft = zeros(size(new_draft));
			pos_mit = find(new_draft<900);       pos_issm = find(new_draft>=900);
			draft(pos_mit) = new_draft(pos_mit); draft(pos_issm) = old_draft(pos_issm);

			%Set draft equal to bed where grounded and equal to zero where no ice is present
			pos_gr = find(mask_oce>0);  draft(pos_gr) = bed(pos_gr); % grounded
			pos_0  = find(mask_ice>0); draft(pos_0)  = 0;            % no ice
			clear pos_mit; clear pos_issm; clear pos_gr; clear pos_0;

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

			% find indices of locations where ice shelf retreated
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %Open new cells as necessary
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % find indices of locations where ice shelf retreated
         fname=sprintf('hFacC.save.%010i.data',niter);
         disp(['  reading ocean hFacC file ' fname]);
         h=binread(fname,4,[mit.mesh.Nx, mit.mesh.Ny, mit.mesh.Nz]); % hFacC (m)
			msk=sum(h,3);
			msk(find(msk))=1;
			[iw jw]=find(msk); % horizontal indices where there is water
			tmp=reshape(draft,[Nx,Ny])-reshape(old_draft,[Nx Ny]);
			tmp(find(tmp<0))=0;
			[im jm]=find(tmp); % horizontal indices where there is melt

			%Extrapolate T/S to locations where ice shelf retreated
			for i=1:length(im)
				% first try vertical extrapolation
				in=find(h(im(i),jm(i),:));
				if length(in)>0;
					S(im(i),jm(i),1:min(in)) = S(im(i),jm(i),min(in));
					T(im(i),jm(i),1:min(in)) = T(im(i),jm(i),min(in));
					continue
				end

				%If not succesful, use closest neighbor horizontal extrapolation
				[y c]=min((iw-im(i)).^2+(jw-jm(i)).^2);
				salt=squeeze(S(iw(c),jw(c),:)); % salinity profile of closest neighbor
				temp=squeeze(T(iw(c),jw(c),:)); % temperature profile of closest neighbor
				in=find(h(iw(c),jw(c),:));
				salt(1:min(in))=salt(min(in));
				temp(1:min(in))=temp(min(in));
				salt(max(in):end)=salt(max(in));
				temp(max(in):end)=temp(max(in));
				S(im(i),jm(i),:)=salt;
				T(im(i),jm(i),:)=temp;
			end

			%Make sure theta and salinity are within reasonable bounds (only impacts a few grid cells)
			posT = find(T(:)>10 | T(:)==0); T(posT) = -2;
			posS = find(S(:)>38 | S(:)==0); S(posS) = 34;
			clear posT; clear posS;

			%Write ocean initial conditions to MITgcm directory
			writebin([mitgcm_dir salt_file] ,S);
			writebin([mitgcm_dir theta_file],T);
			writebin([mitgcm_dir uvel_file] ,U);
			writebin([mitgcm_dir vvel_file] ,V);
			writebin([mitgcm_dir etan_file] ,E);
			writebin([mitgcm_dir draft_file],reshape(draft,[Nx Ny]));

			%Read seaice pickup file
			fnm=['run_big_200wc/pickup_seaice.' dates_name(t-1,:) '.data'];
			Area_si = readbin(fnm,[Nx Ny],1,'real*8',7);
			Heff_si = readbin(fnm,[Nx Ny],1,'real*8',8);
			Snow_si = readbin(fnm,[Nx Ny],1,'real*8',9);
			Salt_si = readbin(fnm,[Nx Ny],1,'real*8',10);

			%Write seaice initial conditions to MITgcm directory
			writebin([mitgcm_dir Area_file] ,Area_si);
			writebin([mitgcm_dir Heff_file] ,Heff_si);
			writebin([mitgcm_dir Hsnow_file] ,Snow_si);
			writebin([mitgcm_dir Hsalt_file] ,Salt_si);

			% }}}
		end		 
		% {{{ Edit calander start date (data.cal MITgcm file):

		%Load data.cal lines into cell A
		fidi = fopen([mitgcm_dir '/data.cal'],'r');
		tline = fgetl(fidi);
		i = 1; B = {}; B{i} = tline;
		while ischar(tline)
			i = i+1;
			tline = fgetl(fidi);
			B{i} = tline;
		end
		fclose(fidi);

		%Change calandar start date
		B{7} = [' startDate_1=' dates(t,:) ','];

		%Write new data.cal file
		fido = fopen([mitgcm_dir '/data.cal'],'w');
		for i=1:numel(B)
			if B{i+1} == -1
				fprintf(fido,'%s',B{i});
				break
			else
				fprintf(fido,'%s\n',B{i});
			end
		end
		fclose(fido);

		%Load data.obcs lines into cell A
		fidi = fopen([mitgcm_dir '/data.obcs'],'r');
		tline = fgetl(fidi);
		A={}; i = 1; A{i} = tline;
		while ischar(tline)
			i = i+1;
			tline = fgetl(fidi);
			A{i} = tline;
		end
		fclose(fidi);

		%Set BC file based on timestep number
		%I split up the ocean boundary condition files into 24-yr chuncks because
		%I was looping the original atmospheric forcing, which was provided for 24 yrs (1992-2016)
		%while using CMIP6-forced ocean boundary conditions. Looking back, I should have just
		%created new atospheric boundary forcing datasets that looped through 2100.
		if t<(628*1)+1
			num = 1;
		elseif t>=(628*1)+1 & t<(628*2)+1
			num = 2;
		elseif t>=(628*2)+1 & t<(628*3)+1
			num=3;
		else
			num=4;
		end

		%Change file names of BCs based on timestep
		A{16} = [' OBNsFile=''FOBNs_m_' num2str(num) '.bin'','];
		A{17} = [' OBNtFile=''FOBNt_m_' num2str(num) '.bin'','];
		A{18} = [' OBNuFile=''FOBNu_m_' num2str(num) '.bin'','];
		A{19} = [' OBNvFile=''FOBNv_m_' num2str(num) '.bin'','];
		A{26} = [' OBEsFile=''FOBEs_m_' num2str(num) '.bin'','];
		A{27} = [' OBEtFile=''FOBEt_m_' num2str(num) '.bin'','];
		A{28} = [' OBEuFile=''FOBEu_m_' num2str(num) '.bin'','];
		A{29} = [' OBEvFile=''FOBEv_m_' num2str(num) '.bin'','];
		A{31} = [' OBWsFile=''FOBWs_m_' num2str(num) '.bin'','];
		A{32} = [' OBWtFile=''FOBWt_m_' num2str(num) '.bin'','];
		A{33} = [' OBWuFile=''FOBWu_m_' num2str(num) '.bin'','];
		A{34} = [' OBWvFile=''FOBWv_m_' num2str(num) '.bin'','];

		%Write new data.cal file
		fido = fopen([mitgcm_dir '/data.obcs'],'w');
		for i=1:numel(A)
			if A{i+1} == -1
				fprintf(fido,'%s',A{i});
				break
			else
				fprintf(fido,'%s\n',A{i});
			end
		end
		fclose(fido);

		%Copy files
		copyfile([mitgcm_dir '/data.obcs'],[mitgcm_dir '/data.obcs~'])
		copyfile([mitgcm_dir '/data.cal'],[mitgcm_dir '/data.cal~'])
		%}}}
		% {{{ System Call to run MITgcm:
		cd(mitgcm_dir)
		!ln -sf /nobackup/hzhang1/forcing/era_xx/
		eval(['!mpiexec -np ' int2str(nPx*nPy) ' ./mitgcmuv']);
		pause('on');

		%Check for correct file size
		%Many times, MATLAB would collect the files before they were finished being written, which would crash the code. 
		%So here I make sure that the needed files are the correct size before progressing.
		%NOTE you will need to collect the correct sizes for the files and declare them in the first step.
		while 1
			f1=dir(['pickup.' c_iter '.data']);
			f2=dir(['SHICE_fwFluxtave.' c_iter '.data']);
			f3=dir(['hFacC.data']);
			f4=dir(['pickup_seaice.' c_iter '.data']);
			%If file size is correct, save files to run directory
			if (f1.bytes==pickup_size) && (f2.bytes==fwflux_size) && (f3.bytes==hFacC_size) && (f4.bytes==pickup_seaice_size)
				copyfile(['pickup.' c_iter '.data'],[pres_dir '/run_big_200wc/pickup.' dates_name(t,:) '.data']);
				copyfile(['pickup_seaice.' c_iter '.data'],[pres_dir '/run_big_200wc/pickup_seaice.' dates_name(t,:) '.data']);
				copyfile(['SHICE_fwFluxtave.' c_iter '.data'],[pres_dir '/run_big_200wc/SHICE_fwFluxtave.' dates_name(t,:) '.data']);
				copyfile('hFacC.data',[pres_dir '/run_big_200wc/hFacC.' dates_name(t,:) '.data']);
				break;
				%If file size is incorrect (not written yet), wait 2 sec and try again	
			else
				pause(2); continue;
			end
		end

		%clean up folder
		!rm *.log STD* *.data *.meta 
		% }}}
		% {{{ Run ISSM and save results:

		%Move to present directory
		cd(pres_dir);

		%Setting Controls
		md.inversion.iscontrol       = 0;
		md.transient.ismasstransport = 1;
		md.transient.isstressbalance = 1;
		md.transient.isgroundingline = 1;
		md.transient.ismovingfront   = 0;
		md.transient.isthermal       = 0;
		md.transient.isslr           = 0;
		md.timestepping.start_time   = (t/24)-(1/24);
		md.timestepping.final_time   = t/24;
		md.timestepping.time_step    = 1/24;
		md.levelset.kill_icebergs    = 1;

		%Load smb
		md.smb.mass_balance = SMB_CNRM_CM6_ssp126;

		%Interpolate mitgcm melt to issm mesh
		melting_rate = readbin(['run_big_200wc/SHICE_fwFluxtave.' dates_name(t,:) '.data'],[Nx Ny]);
		m_mit = -melting_rate(:)*y2s/rho_ice;
		m_tot = InterpFromMeshToMesh2d(mdm.mesh.elements,mdm.mesh.x,mdm.mesh.y,m_mit,md.mesh.x,md.mesh.y);

		%Set basal melting rate fields
		md.basalforcings = basalforcings();
		md.basalforcings.floatingice_melting_rate = m_tot;
		md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);

		%Solve
		md.miscellaneous.name = ['run_lowemission_' dates_name(t,:)];
		md.cluster=generic('name',oshostname(),'np',15);
		md.verbose.solution = 1;
		md=solve(md,'tr');

		%Reset model
		md.geometry.base             = md.results.TransientSolution(end).Base;
		md.geometry.surface          = md.results.TransientSolution(end).Surface;
		md.geometry.thickness        = md.geometry.surface-md.geometry.base;
		md.initialization.vx         = md.results.TransientSolution(end).Vx;
		md.initialization.vy         = md.results.TransientSolution(end).Vy;
		md.initialization.vel        = md.results.TransientSolution(end).Vel;
		md.initialization.pressure   = md.results.TransientSolution(end).Pressure;
		md.mask.groundedice_levelset = md.results.TransientSolution(end).MaskGroundediceLevelset;

		%Save results and story in results_200wc directory
		fname = ['results_200wc/results_' dates_name(t,:) '.mat'];
		results = md.results;
		save(fname,'results');
		clear md.results;

		%Clear execution folder in ISSM to avoid going over quota
		cd('/home1/tpelle1/trunk-jpl/execution');
		!rm -r run_lowemission_*
		cd(pres_dir);

		%Save
		% }}}
	end
	% }}}


%Concatenate ISSM results and save in results_200wc folder
% {{{ ConcatenateIssmResults:
if perform(org,'ConcatenateIssmResults'),
	%Move to present directory
	loaddata(org,'Parameters');
	cd(pres_dir);

	%Get all results files
	cd results_200wc;
	!rm results.mat;
	myfiles = dir(fullfile(pwd,'results*.mat'));
	cd ..;

	%Loop over result-files and build complete results matrix
	results_full = [];
	for i=1:length(myfiles)
		filename = ['results_200wc/' myfiles(i).name];
		load(filename);
		results_full = [results_full results];
		clear results;
	end

	%Save
	results=results_full;
	save('results_200wc/results.mat','results','-v7.3');
end
% }}}

%Run Coupled ISSM and MITgcm - start from advanced state
% {{{ RunCoupledModels_AdvStart:
if perform(org,'RunCoupledModels_AdvStart'),
	% {{{ General settings:
	%Load parameters
	loaddata(org,'Parameters');
	cd(pres_dir)

	%Load models (md=issm, msm=mitgcm) and results from unfinished run
	md=model; mdm=model;
	md=loadmodel(org,'MITgcmModel'); mdm=md;
	md=loadmodel(org,'TotModel');
	load results_200wc/results.mat;

	%Clear MITgcm
	cd(mitgcm_dir)
	!rm *
	!cp ../input_lowemission/* .
	cd(pres_dir)

	%Get MITgcm bed and ISSM ice mask for to use in draft-interpolation
	bed =readbin([mitgcm_dir bed_file],[Nx*Ny 1]);
	mask_ice = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ice_levelset,mdm.mesh.x,mdm.mesh.y,'default',-1);

	%Load smb
	load Models/SMB_CNRM_CM6_ssp126;

	% }}}
	% {{{ ISSM settings:
	%Reset model
	md.geometry.base             = results(end).TransientSolution(1).Base;
	md.geometry.surface          = results(end).TransientSolution(1).Surface;
	md.geometry.thickness        = md.geometry.surface-md.geometry.base;
	md.initialization.vx         = results(end).TransientSolution(1).Vx;
	md.initialization.vy         = results(end).TransientSolution(1).Vy;
	md.initialization.vel        = results(end).TransientSolution(1).Vel;
	md.initialization.pressure   = results(end).TransientSolution(1).Pressure;
	md.mask.groundedice_levelset = results(end).TransientSolution(1).MaskGroundediceLevelset;

	%SLR parameters (to satisfy model-consistency requirements, will not be used)
	md.slr.deltathickness = zeros(md.mesh.numberofelements,1);
	md.slr.sealevel       = zeros(md.mesh.numberofvertices,1);
	md.slr.spcthickness   = zeros(md.mesh.numberofvertices,1);
	md.slr.hydro_rate     = zeros(md.mesh.numberofvertices,1);
	md.slr.Ugia           = zeros(md.mesh.numberofvertices,1);
	md.slr.Ngia           = zeros(md.mesh.numberofvertices,1);

	%Set floating ice parameters
	md.transient.requested_outputs={'default','BasalforcingsFloatingiceMeltingRate','MaskIceLevelset','MaskGroundediceLevelset','IceVolume','IceVolumeAboveFloatation','GroundedArea','FloatingArea','SmbMassBalance','TotalFloatingBmb'};
	pos=find(md.mesh.vertexonboundary);
	md.masstransport.spcthickness=NaN(md.mesh.numberofvertices,1);
	md.masstransport.spcthickness(pos)=md.geometry.thickness(pos);

	%Set GL and Fr interpolation schemes
	md.groundingline.migration              = 'SubelementMigration';
	md.groundingline.friction_interpolation = 'SubelementFriction1';
	md.groundingline.melt_interpolation     = 'SubelementMelt1';

	%Set inversion
	md.inversion=m1qn3inversion(md.inversion);
	% }}}

	%Begin looping over dates
	for t=length(results)+1:size(dates,1)
		 disp(['------------------ Date of projection: ' dates_name(t,:)])
       if t>1
% {{{Send draft and initial condition files to MITgcm:

	%Get old draft, new draft from ISSM, and ice masks from ISSM - interpolate onto MITgcm grid
	old_draft = readbin([mitgcm_dir draft_file],[Nx*Ny 1]);
	new_draft = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.geometry.base,mdm.mesh.x,mdm.mesh.y,'default',999);
	mask_gr   = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.groundedice_levelset,mdm.mesh.x,mdm.mesh.y,'default',-1);

	%Take new draft from issm (>900) and sub it into mitgcm draft (<900)
	draft = zeros(size(new_draft));
	pos_mit = find(new_draft<900);       pos_issm = find(new_draft>=900);
	draft(pos_mit) = new_draft(pos_mit); draft(pos_issm) = old_draft(pos_issm);

	%Set draft equal to bed where grounded and equal to zero where no ice is present
	pos_gr = find(mask_gr>0);  draft(pos_gr) = bed(pos_gr);
   pos_0  = find(mask_ice>0); draft(pos_0)  = 0;
	clear pos_mit; clear pos_issm; clear pos_gr; clear pos_0;

	%Read ocean pickup file 
	fnm=['run_big_200wc/pickup.' dates_name(t-1,:) '.data'];
	U=readbin(fnm,[Nx Ny Nz],1,'real*8',0);
	V=readbin(fnm,[Nx Ny Nz],1,'real*8',1);
	T=readbin(fnm,[Nx Ny Nz],1,'real*8',2);
	S=readbin(fnm,[Nx Ny Nz],1,'real*8',3);
	E=readbin(fnm,[Nx Ny],1,'real*8',8);

	%Find indices of locations where ice shelf retreated
	h=readbin(['run_big_200wc/hFacC.' dates_name(t-1,:) '.data'],[Nx Ny Nz]);
	msk=sum(h,3);
	msk(find(msk))=1;
	[iw jw]=find(msk); % horizontal indices where there is water
	tmp=reshape(draft,[Nx,Ny])-reshape(old_draft,[Nx Ny]);
	tmp(find(tmp<0))=0;
	[im jm]=find(tmp); % horizontal indices where there is melt

	%Extrapolate T/S to locations where ice shelf retreated
	for i=1:length(im)

		%First try vertical extrapolation
		in=find(h(im(i),jm(i),:));
		if length(in)>0;
			S(im(i),jm(i),1:min(in)) = S(im(i),jm(i),min(in));
			T(im(i),jm(i),1:min(in)) = T(im(i),jm(i),min(in));
			continue
		end

		%If not succesful, use closest neighbor horizontal extrapolation
		[y c]=min((iw-im(i)).^2+(jw-jm(i)).^2);
		salt=squeeze(S(iw(c),jw(c),:)); %Salinity profile of closest neighbor
		temp=squeeze(T(iw(c),jw(c),:)); %Temperature profile of closest neighbor
		in=find(h(iw(c),jw(c),:));
		salt(1:min(in))=salt(min(in));
		temp(1:min(in))=temp(min(in));
		salt(max(in):end)=salt(max(in));
		temp(max(in):end)=temp(max(in));
		S(im(i),jm(i),:)=salt;
		T(im(i),jm(i),:)=temp;
	end

	%Make sure theta and salinity are within reasonable bounds
	posT = find(T(:)>10 | T(:)==0); T(posT) = -2;
	posS = find(S(:)>38 | S(:)==0); S(posS) = 34;
	clear posT; clear posS;

	%Write ocean initial conditions to MITgcm directory
	writebin([mitgcm_dir salt_file] ,S);
	writebin([mitgcm_dir theta_file],T);
	writebin([mitgcm_dir uvel_file] ,U);
	writebin([mitgcm_dir vvel_file] ,V);
	writebin([mitgcm_dir etan_file] ,E);
	writebin([mitgcm_dir draft_file],reshape(draft,[Nx Ny]));

   %Read seaice pickup file
   fnm=['run_big_200wc/pickup_seaice.' dates_name(t-1,:) '.data'];
   Area_si = readbin(fnm,[Nx Ny],1,'real*8',7);
	Heff_si = readbin(fnm,[Nx Ny],1,'real*8',8);
	Snow_si = readbin(fnm,[Nx Ny],1,'real*8',9);
	Salt_si = readbin(fnm,[Nx Ny],1,'real*8',10);

	%Write seaice initial condition files to MITgcm directory
	writebin([mitgcm_dir Area_file] ,Area_si);
	writebin([mitgcm_dir Heff_file] ,Heff_si);
	writebin([mitgcm_dir Hsnow_file] ,Snow_si);
	writebin([mitgcm_dir Hsalt_file] ,Salt_si);

% }}}
		 end		 
% {{{ Edit calander start date (data.cal MITgcm file):
   
	%Load data.cal lines into cell A
	fidi = fopen([mitgcm_dir '/data.cal'],'r');
	tline = fgetl(fidi);
	i = 1; B = {}; B{i} = tline;
	while ischar(tline)
		i = i+1;
		tline = fgetl(fidi);
		B{i} = tline;
	end
	fclose(fidi);

	%Change calandar start date
	B{7} = [' startDate_1=' dates(t,:) ','];

	%Write new data.cal file
	fido = fopen([mitgcm_dir '/data.cal'],'w');
	for i=1:numel(B)
		if B{i+1} == -1
			fprintf(fido,'%s',B{i});
			break
		else
			fprintf(fido,'%s\n',B{i});
		end
	end
	fclose(fido);

	%Load data.obcs lines into cell A
	fidi = fopen([mitgcm_dir '/data.obcs'],'r');
	tline = fgetl(fidi);
	A={}; i = 1; A{i} = tline;
	while ischar(tline)
		i = i+1;
		tline = fgetl(fidi);
		A{i} = tline;
	end
	fclose(fidi);

	%Set BC file based on timestep number
	if t<(628*1)+1
		num = 1;
	elseif t>=(628*1)+1 & t<(628*2)+1
		num = 2;
	elseif t>=(628*2)+1 & t<(628*3)+1
		num=3;
	else
		num=4;
	end

	%Change file names of BCs based on timestep
	A{16} = [' OBNsFile=''FOBNs_m_' num2str(num) '.bin'','];
	A{17} = [' OBNtFile=''FOBNt_m_' num2str(num) '.bin'','];
	A{18} = [' OBNuFile=''FOBNu_m_' num2str(num) '.bin'','];
	A{19} = [' OBNvFile=''FOBNv_m_' num2str(num) '.bin'','];
	A{26} = [' OBEsFile=''FOBEs_m_' num2str(num) '.bin'','];
	A{27} = [' OBEtFile=''FOBEt_m_' num2str(num) '.bin'','];
	A{28} = [' OBEuFile=''FOBEu_m_' num2str(num) '.bin'','];
	A{29} = [' OBEvFile=''FOBEv_m_' num2str(num) '.bin'','];
	A{31} = [' OBWsFile=''FOBWs_m_' num2str(num) '.bin'','];
	A{32} = [' OBWtFile=''FOBWt_m_' num2str(num) '.bin'','];
	A{33} = [' OBWuFile=''FOBWu_m_' num2str(num) '.bin'','];
	A{34} = [' OBWvFile=''FOBWv_m_' num2str(num) '.bin'','];

	%Write new data.cal file
	fido = fopen([mitgcm_dir '/data.obcs'],'w');
	for i=1:numel(A)
		if A{i+1} == -1
			fprintf(fido,'%s',A{i});
	      break
	   else
		   fprintf(fido,'%s\n',A{i});
	   end
   end
	fclose(fido);

   %Copy files
	copyfile([mitgcm_dir '/data.obcs'],[mitgcm_dir '/data.obcs~'])
	copyfile([mitgcm_dir '/data.cal'],[mitgcm_dir '/data.cal~'])
%}}}
% {{{ System Call to run MITgcm:
	cd(mitgcm_dir)
	!ln -sf /nobackup/hzhang1/forcing/era_xx/
   eval(['!mpiexec -np ' int2str(nPx*nPy) ' ./mitgcmuv']);
	pause('on');

	%Check for correct file size
	while 1
		f1=dir(['pickup.' c_iter '.data']);
		f2=dir(['SHICE_fwFluxtave.' c_iter '.data']);
		f3=dir(['hFacC.data']);
		f4=dir(['pickup_seaice.' c_iter '.data']);
		%If file size is correct, save files to run directory
		if (f1.bytes==pickup_size) && (f2.bytes==fwflux_size) && (f3.bytes==hFacC_size) && (f4.bytes==pickup_seaice_size)
			copyfile(['pickup.' c_iter '.data'],[pres_dir '/run_big_200wc/pickup.' dates_name(t,:) '.data']);
			copyfile(['pickup_seaice.' c_iter '.data'],[pres_dir '/run_big_200wc/pickup_seaice.' dates_name(t,:) '.data']);
			copyfile(['SHICE_fwFluxtave.' c_iter '.data'],[pres_dir '/run_big_200wc/SHICE_fwFluxtave.' dates_name(t,:) '.data']);
			copyfile('hFacC.data',[pres_dir '/run_big_200wc/hFacC.' dates_name(t,:) '.data']);
			break;
			%If file size is incorrect (not written yet), wait 2 sec and try again	
		else
			pause(2); continue;
		end
	end

	%clean up folder
   !rm *.log STD* *.data *.meta 
% }}}
	% {{{ Run ISSM and save results:

	%Move to present directory
	cd(pres_dir);

	%Setting Controls
	md.inversion.iscontrol       = 0;
	md.transient.ismasstransport = 1;
	md.transient.isstressbalance = 1;
	md.transient.isgroundingline = 1;
	md.transient.ismovingfront   = 0;
	md.transient.isthermal       = 0;
	md.transient.isslr           = 0;
	md.timestepping.start_time   = (t/24)-(1/24);
	md.timestepping.final_time   = t/24;
	md.timestepping.time_step    = 1/24;
	md.levelset.kill_icebergs    = 1;

	%Load smb
	md.smb.mass_balance = SMB_CNRM_CM6_ssp126;

	%Interpolate mitgcm melt to issm mesh
	melting_rate = readbin(['run_big_200wc/SHICE_fwFluxtave.' dates_name(t,:) '.data'],[Nx Ny]);
	m_mit = -melting_rate(:)*y2s/rho_ice;
	m_tot = InterpFromMeshToMesh2d(mdm.mesh.elements,mdm.mesh.x,mdm.mesh.y,m_mit,md.mesh.x,md.mesh.y);

	%Set basal melting rate fields
	md.basalforcings = basalforcings();
	md.basalforcings.floatingice_melting_rate = m_tot;
	md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);

	%Solve
	md.miscellaneous.name = ['run_lowemission_' dates_name(t,:)];
	md.cluster=generic('name',oshostname(),'np',15);
	md.verbose.solution = 1;
	md=solve(md,'tr');

	%Reset model
	md.geometry.base             = md.results.TransientSolution(end).Base;
	md.geometry.surface          = md.results.TransientSolution(end).Surface;
	md.geometry.thickness        = md.geometry.surface-md.geometry.base;
	md.initialization.vx         = md.results.TransientSolution(end).Vx;
	md.initialization.vy         = md.results.TransientSolution(end).Vy;
	md.initialization.vel        = md.results.TransientSolution(end).Vel;
	md.initialization.pressure   = md.results.TransientSolution(end).Pressure;
	md.mask.groundedice_levelset = md.results.TransientSolution(end).MaskGroundediceLevelset;

	%Save results and story in results_200wc directory
	fname = ['results_200wc/results_' dates_name(t,:) '.mat'];
	results = md.results;
	save(fname,'results');
	clear md.results;

	%Clear execution folder in ISSM to avoid going over quota
	cd('/home1/tpelle1/trunk-jpl/execution');
	!rm -r run_lowemission_*
	cd(pres_dir);

	%Save
	% }}}
end
end
% }}}

%These are scripts I used to test certain parts of the code and to perform data analysis
%(e.g. computing mean ocean temperatures near the ice fronts of select glaciers or mean
%speeds of currents in certain areas of the model grid). I figured these could be helpful
%for you when you get further into the Thwaites ice-ocean research. 
% {{{ RunIssm:
if perform(org,'RunIssm'),

	%Set ISSM
	md=model; mdm=model;
	md=loadmodel(org,'MITgcmModel'); mdm=md;
	md=loadmodel(org,'TotModel');
	load results_200wc/results.mat;
   t = length(results)+1;

	md.geometry.base             = results(end).TransientSolution(1).Base;
	md.geometry.surface          = results(end).TransientSolution(1).Surface;
	md.geometry.thickness        = md.geometry.surface-md.geometry.base;
	md.initialization.vx         = results(end).TransientSolution(1).Vx;
	md.initialization.vy         = results(end).TransientSolution(1).Vy;
	md.initialization.vel        = results(end).TransientSolution(1).Vel;
	md.initialization.pressure   = results(end).TransientSolution(1).Pressure;
	md.mask.groundedice_levelset = results(end).TransientSolution(1).MaskGroundediceLevelset;

	%Setting Controls
	md.inversion.iscontrol       = 0;
	md.transient.ismasstransport = 1;
	md.transient.isstressbalance = 1;
	md.transient.isgroundingline = 1;
	md.transient.ismovingfront   = 0;
	md.transient.isthermal       = 0;
	md.transient.isslr           = 0;
	md.timestepping.final_time   = 1/24;
	md.timestepping.time_step    = 1/24;
	md.levelset.kill_icebergs    = 1;

	%SLR parameters (to satisfy model-consistency requirements, will not be used)
	md.slr.deltathickness = zeros(md.mesh.numberofelements,1);
	md.slr.sealevel       = zeros(md.mesh.numberofvertices,1);
	md.slr.spcthickness   = zeros(md.mesh.numberofvertices,1);
	md.slr.hydro_rate     = zeros(md.mesh.numberofvertices,1);
	md.slr.Ugia           = zeros(md.mesh.numberofvertices,1);
	md.slr.Ngia           = zeros(md.mesh.numberofvertices,1);

	%Set floating ice parameters
	md.transient.requested_outputs={'default','BasalforcingsFloatingiceMeltingRate','MaskIceLevelset','MaskGroundediceLevelset','IceVolume','IceVolumeAboveFloatation','GroundedArea','FloatingArea','SmbMassBalance','TotalFloatingBmb'};
	pos=find(md.mesh.vertexonboundary);
	md.masstransport.spcthickness=NaN(md.mesh.numberofvertices,1);
	md.masstransport.spcthickness(pos)=md.geometry.thickness(pos);

	%Set GL and Fr interpolation schemes
	md.groundingline.migration              = 'SubelementMigration';
	md.groundingline.friction_interpolation = 'SubelementFriction1';
	md.groundingline.melt_interpolation     = 'SubelementMelt1';

	%Set inversion
	md.inversion=m1qn3inversion(md.inversion);

	%Interpolate mitgcm melt to issm mesh
	loaddata(org,'Parameters');
	melting_rate = readbin(['run_big_200wc/SHICE_fwFluxtave.' dates_name(t,:) '.data'],[Nx Ny]);
	m_mit = -melting_rate(:)*y2s/rho_ice;
	m_tot = InterpFromMeshToMesh2d(mdm.mesh.elements,mdm.mesh.x,mdm.mesh.y,m_mit,md.mesh.x,md.mesh.y);

	%Set basal melting rate fields
	md.basalforcings = basalforcings();
	md.basalforcings.floatingice_melting_rate = m_tot;
	md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);

	%Solve
	md.miscellaneous.name = ['run_lowemission_' dates_name(t,:)];
	md.cluster=generic('name',oshostname(),'np',15);
	md.verbose.solution = 1;
	md=solve(md,'tr');

	%Reset model
	md.geometry.base             = md.results.TransientSolution(end).Base;
	md.geometry.surface          = md.results.TransientSolution(end).Surface;
	md.geometry.thickness        = md.geometry.surface-md.geometry.base;
	md.initialization.vx         = md.results.TransientSolution(end).Vx;
	md.initialization.vy         = md.results.TransientSolution(end).Vy;
	md.initialization.vel        = md.results.TransientSolution(end).Vel;
	md.initialization.pressure   = md.results.TransientSolution(end).Pressure;
	md.mask.groundedice_levelset = md.results.TransientSolution(end).MaskGroundediceLevelset;

	%Save results and story in results_200wc directory
	fname = ['results_200wc/results_' dates_name(t,:) '.mat'];
	results = md.results;
	save(fname,'results');
	clear md.results;
end
% }}}
% {{{ ComputeBottomTempTimeSeriesTot:
if perform(org,'ComputeBottomTempTimeSeriesTot'),

	%Load model
	load 'Models/IceOcean_MITgcmModel.mat';

	%Set domain size
	Nx=690; Ny=280; Nz=70;

	%Get all melt files and store names in myfiles
	cd run_big_200wc;
	myfiles = dir(fullfile(pwd,'pickup.*.data'));
	cd ..

	%get areas of ISSM mesh for averaging
	areas = GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y);

	%Loop over melt data (1992 - 2017)
	for i=1:1992
		%Get 3d temperature field
		filename = ['run_big_200wc/' myfiles(i).name];
		T = readbin(filename,[Nx Ny Nz],1,'real*8',2);
		S = readbin(filename,[Nx Ny Nz],1,'real*8',3);

		%Compute bottom temperature
		bottom_T = NaN([Nx Ny]);
		bottom_S = NaN([Nx Ny]);
		for ii=1:Nx
			for j=1:Ny
				Tcol = T(ii,j,:);
				Scol = S(ii,j,:);
				pos = max(find(Tcol~=0));
				if isempty(pos)
					 bottom_T(ii,j) = NaN;
					 bottom_S(ii,j) = NaN;
				 else
					bottom_T(ii,j) = Tcol(pos);
					bottom_S(ii,j) = Scol(pos);
				end
				clear pos; clear Tcol; clear Scol;
			end
		end

		%Average bottom T over elements
		bottomT_el = mean(bottom_T(md.mesh.elements),2);
		bottomS_el = mean(bottom_S(md.mesh.elements),2);

		%Get bottom t in front of Totten ice shelf
		pos = find(ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,'Exp/TotContShelf.exp','element',2));
		bottomT_tot = bottomT_el(pos);
		bottomS_tot = bottomS_el(pos);
		area_tot = areas(pos);

		%Get weighted average of bottom t continental shelf
		topT = sum(area_tot.*bottomT_tot);
		topS = sum(area_tot.*bottomS_tot);
		bttm = sum(area_tot);
		T_picop_lowemission_Tot(i) = topT/bttm;
		S_picop_lowemission_Tot(i) = topS/bttm;

		%Clear variables
		clear pos; clear bottom_T; clear bottomT_el; clear bottomT_tot; clear bottom_S;
		clear bottomS_el; clear bottomS_tot; clear area_tot; clear topT; clear topS; clear bttm;
	end

	%Save
	save('results_200wc/T_picop_lowemission_Tot.mat','T_picop_lowemission_Tot');
	save('results_200wc/S_picop_lowemission_Tot.mat','S_picop_lowemission_Tot');
end
% }}}
% {{{ ComputeBottomTempTimeSeriesMuis:
if perform(org,'ComputeBottomTempTimeSeriesMuis'),

	%Load model
	load 'Models/IceOcean_MITgcmModel.mat';

	%Set domain size
	Nx=690; Ny=280; Nz=70;

	%Get all melt files and store names in myfiles
	cd run_big_200wc;
	myfiles = dir(fullfile(pwd,'pickup.*.data'));
	cd ..

	%get areas of ISSM mesh for averaging
	areas = GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y);

	%Loop over melt data (1992 - 2017)
	for i=1:1992
		%Get 3d temperature field
		filename = ['run_big_200wc/' myfiles(i).name];
		T = readbin(filename,[Nx Ny Nz],1,'real*8',2);
		S = readbin(filename,[Nx Ny Nz],1,'real*8',3);

		%Compute bottom temperature
		bottom_T = NaN([Nx Ny]);
		bottom_S = NaN([Nx Ny]);
		for ii=1:Nx
			for j=1:Ny
				Tcol = T(ii,j,:);
				Scol = S(ii,j,:);
				pos = max(find(Tcol~=0));
				if isempty(pos)
					 bottom_T(ii,j) = NaN;
					 bottom_S(ii,j) = NaN;
				 else
					bottom_T(ii,j) = Tcol(pos);
					bottom_S(ii,j) = Scol(pos);
				end
				clear pos; clear Tcol; clear Scol;
			end
		end

		%Average bottom T over elements
		bottomT_el = mean(bottom_T(md.mesh.elements),2);
		bottomS_el = mean(bottom_S(md.mesh.elements),2);

		%Get bottom t in front of Totten ice shelf
		pos = find(ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,'Exp/MuisContShelf.exp','element',2));
		bottomT_tot = bottomT_el(pos);
		bottomS_tot = bottomS_el(pos);
		area_tot = areas(pos);

		%Get weighted average of bottom t continental shelf
		topT = sum(area_tot.*bottomT_tot);
		topS = sum(area_tot.*bottomS_tot);
		bttm = sum(area_tot);
		T_picop_lowemission_Muis(i) = topT/bttm;
		S_picop_lowemission_Muis(i) = topS/bttm;

		%Clear variables
		clear pos; clear bottom_T; clear bottomT_el; clear bottomT_tot; clear bottom_S;
		clear bottomS_el; clear bottomS_tot; clear area_tot; clear topT; clear topS; clear bttm;
	end

	%Save
	save('results_200wc/T_picop_lowemission_Muis.mat','T_picop_lowemission_Muis');
	save('results_200wc/S_picop_lowemission_Muis.mat','S_picop_lowemission_Muis');
end
% }}}
% {{{ ComputeBottomTempTimeSeriesDen:
if perform(org,'ComputeBottomTempTimeSeriesDen'),

	%Load model
	load 'Models/IceOcean_MITgcmModel.mat';

	%Set domain size
	Nx=690; Ny=280; Nz=70;

	%Get all melt files and store names in myfiles
	cd run_big_200wc;
	myfiles = dir(fullfile(pwd,'pickup.*.data'));
	cd ..

	%get areas of ISSM mesh for averaging
	areas = GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y);

	%Loop over melt data (1992 - 2017)
	for i=1:1992
		%Get 3d temperature field
		filename = ['run_big_200wc/' myfiles(i).name];
		T = readbin(filename,[Nx Ny Nz],1,'real*8',2);
		S = readbin(filename,[Nx Ny Nz],1,'real*8',3);

		%Compute bottom temperature
		bottom_T = NaN([Nx Ny]);
		bottom_S = NaN([Nx Ny]);
		for ii=1:Nx
			for j=1:Ny
				Tcol = T(ii,j,:);
				Scol = S(ii,j,:);
				pos = max(find(Tcol~=0));
				if isempty(pos)
					 bottom_T(ii,j) = NaN;
					 bottom_S(ii,j) = NaN;
				 else
					bottom_T(ii,j) = Tcol(pos);
					bottom_S(ii,j) = Scol(pos);
				end
				clear pos; clear Tcol; clear Scol;
			end
		end

		%Average bottom T over elements
		bottomT_el = mean(bottom_T(md.mesh.elements),2);
		bottomS_el = mean(bottom_S(md.mesh.elements),2);

		%Get bottom t in front of Totten ice shelf
		pos = find(ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,'Exp/DenmanContShelf.exp','element',2));
		bottomT_tot = bottomT_el(pos);
		bottomS_tot = bottomS_el(pos);
		area_tot = areas(pos);

		%Get weighted average of bottom t continental shelf
		topT = sum(area_tot.*bottomT_tot);
		topS = sum(area_tot.*bottomS_tot);
		bttm = sum(area_tot);
		T_picop_Den_lowemission(i) = topT/bttm;
		S_picop_Den_lowemission(i) = topS/bttm;

		%Clear variables
		clear pos; clear bottom_T; clear bottomT_el; clear bottomT_tot; clear bottom_S;
		clear bottomS_el; clear bottomS_tot; clear area_tot; clear topT; clear topS; clear bttm;
	end

	%Save
	save('results_200wc/T_picop_Den_lowemission.mat','T_picop_Den_lowemission');
	save('results_200wc/S_picop_Den_lowemission.mat','S_picop_Den_lowemission');
end
% }}}
% {{{ ComputeDaTempTimeSeriesDen:
if perform(org,'ComputeDaTempTimeSeriesDen'),

	%Load model
	load 'Models/IceOcean_MITgcmModel.mat';

	%Set domain size
	Nx=690; Ny=280; Nz=70;

	%Get all melt files and store names in myfiles
	cd run_big_200wc;
	myfiles = dir(fullfile(pwd,'pickup.*.data'));
	cd ..

	%get areas of ISSM mesh for averaging
	areas = GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y);

	%Loop over melt data (1992 - 2017)
	for i=1:1992
		%Get 3d temperature field
		filename = ['run_big_200wc/' myfiles(i).name];
		T = readbin(filename,[Nx Ny Nz],1,'real*8',2);
		S = readbin(filename,[Nx Ny Nz],1,'real*8',3);

		%Compute bottom temperature
		bottom_T = NaN([Nx Ny]);
		bottom_S = NaN([Nx Ny]);
		for ii=1:Nx
			for j=1:Ny
				Tcol = T(ii,j,:);
				Scol = S(ii,j,:);
				pos = max(find(Tcol~=0));
				if isempty(pos)
					 bottom_T(ii,j) = NaN;
					 bottom_S(ii,j) = NaN;
				 else
					 Tcol(find(Tcol==0))=NaN;
					 Scol(find(Scol==0))=NaN;
					 bottom_T(ii,j) = nanmean(Tcol);
					 bottom_S(ii,j) = nanmean(Scol);
				end
				clear pos; clear Tcol; clear Scol;
			end
		end

		%Average bottom T over elements
		bottomT_el = mean(bottom_T(md.mesh.elements),2);
		bottomS_el = mean(bottom_S(md.mesh.elements),2);

		%Get bottom t in front of Totten ice shelf
		pos = find(ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,'Exp/DenmanContShelf.exp','element',2));
		bottomT_tot = bottomT_el(pos);
		bottomS_tot = bottomS_el(pos);
		area_tot = areas(pos);

		%Get weighted average of bottom t continental shelf
		topT = sum(area_tot.*bottomT_tot);
		topS = sum(area_tot.*bottomS_tot);
		bttm = sum(area_tot);
		T_picop_Den_lowemission_DA(i) = topT/bttm;
		S_picop_Den_lowemission_DA(i) = topS/bttm;

		%Clear variables
		clear pos; clear bottom_T; clear bottomT_el; clear bottomT_tot; clear bottom_S;
		clear bottomS_el; clear bottomS_tot; clear area_tot; clear topT; clear topS; clear bttm;
	end

	%Save
	save('results_200wc/T_picop_Den_lowemission_DA.mat','T_picop_Den_lowemission_DA');
	save('results_200wc/S_picop_Den_lowemission_DA.mat','S_picop_Den_lowemission_DA');
end
% }}}
% {{{ ComputeDaTempTimeSeriesTot:
if perform(org,'ComputeDaTempTimeSeriesTot'),

	%Load model
	load 'Models/IceOcean_MITgcmModel.mat';

	%Set domain size
	Nx=690; Ny=280; Nz=70;

	%Get all melt files and store names in myfiles
	cd run_big_200wc;
	myfiles = dir(fullfile(pwd,'pickup.*.data'));
	cd ..

	%get areas of ISSM mesh for averaging
	areas = GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y);

	%Loop over melt data (1992 - 2017)
	for i=1:1992
		%Get 3d temperature field
		filename = ['run_big_200wc/' myfiles(i).name];
		T = readbin(filename,[Nx Ny Nz],1,'real*8',2);
		S = readbin(filename,[Nx Ny Nz],1,'real*8',3);

		%Compute bottom temperature
		bottom_T = NaN([Nx Ny]);
		bottom_S = NaN([Nx Ny]);
		for ii=1:Nx
			for j=1:Ny
				Tcol = T(ii,j,:);
				Scol = S(ii,j,:);
				pos = max(find(Tcol~=0));
				if isempty(pos)
					 bottom_T(ii,j) = NaN;
					 bottom_S(ii,j) = NaN;
				 else
					 Tcol(find(Tcol==0))=NaN;
					 Scol(find(Scol==0))=NaN;
					 bottom_T(ii,j) = nanmean(Tcol);
					 bottom_S(ii,j) = nanmean(Scol);
				end
				clear pos; clear Tcol; clear Scol;
			end
		end

		%Average bottom T over elements
		bottomT_el = mean(bottom_T(md.mesh.elements),2);
		bottomS_el = mean(bottom_S(md.mesh.elements),2);

		%Get bottom t in front of Totten ice shelf
		pos = find(ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,'Exp/TotContShelf.exp','element',2));
		bottomT_tot = bottomT_el(pos);
		bottomS_tot = bottomS_el(pos);
		area_tot = areas(pos);

		%Get weighted average of bottom t continental shelf
		topT = sum(area_tot.*bottomT_tot);
		topS = sum(area_tot.*bottomS_tot);
		bttm = sum(area_tot);
		T_picop_Tot_lowemission_DA(i) = topT/bttm;
		S_picop_Tot_lowemission_DA(i) = topS/bttm;

		%Clear variables
		clear pos; clear bottom_T; clear bottomT_el; clear bottomT_tot; clear bottom_S;
		clear bottomS_el; clear bottomS_tot; clear area_tot; clear topT; clear topS; clear bttm;
	end

	%Save
	save('results_200wc/T_picop_Tot_lowemission_DA.mat','T_picop_Tot_lowemission_DA');
	save('results_200wc/S_picop_Tot_lowemission_DA.mat','S_picop_Tot_lowemission_DA');
end
% }}}
% {{{ ComputeDaTempTimeSeriesMuis:
if perform(org,'ComputeDaTempTimeSeriesMuis'),

	%Load model
	load 'Models/IceOcean_MITgcmModel.mat';

	%Set domain size
	Nx=690; Ny=280; Nz=70;

	%Get all melt files and store names in myfiles
	cd run_big_200wc;
	myfiles = dir(fullfile(pwd,'pickup.*.data'));
	cd ..

	%get areas of ISSM mesh for averaging
	areas = GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y);

	%Loop over melt data (1992 - 2017)
	for i=1:1992
		%Get 3d temperature field
		filename = ['run_big_200wc/' myfiles(i).name];
		T = readbin(filename,[Nx Ny Nz],1,'real*8',2);
		S = readbin(filename,[Nx Ny Nz],1,'real*8',3);

		%Compute bottom temperature
		bottom_T = NaN([Nx Ny]);
		bottom_S = NaN([Nx Ny]);
		for ii=1:Nx
			for j=1:Ny
				Tcol = T(ii,j,:);
				Scol = S(ii,j,:);
				pos = max(find(Tcol~=0));
				if isempty(pos)
					 bottom_T(ii,j) = NaN;
					 bottom_S(ii,j) = NaN;
				 else
					 Tcol(find(Tcol==0))=NaN;
					 Scol(find(Scol==0))=NaN;
					 bottom_T(ii,j) = nanmean(Tcol);
					 bottom_S(ii,j) = nanmean(Scol);
				end
				clear pos; clear Tcol; clear Scol;
			end
		end

		%Average bottom T over elements
		bottomT_el = mean(bottom_T(md.mesh.elements),2);
		bottomS_el = mean(bottom_S(md.mesh.elements),2);

		%Get bottom t in front of Totten ice shelf
		pos = find(ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,'Exp/MuisContShelf.exp','element',2));
		bottomT_tot = bottomT_el(pos);
		bottomS_tot = bottomS_el(pos);
		area_tot = areas(pos);

		%Get weighted average of bottom t continental shelf
		topT = sum(area_tot.*bottomT_tot);
		topS = sum(area_tot.*bottomS_tot);
		bttm = sum(area_tot);
		T_picop_Muis_lowemission_DA(i) = topT/bttm;
		S_picop_Muis_lowemission_DA(i) = topS/bttm;

		%Clear variables
		clear pos; clear bottom_T; clear bottomT_el; clear bottomT_tot; clear bottom_S;
		clear bottomS_el; clear bottomS_tot; clear area_tot; clear topT; clear topS; clear bttm;
	end

	%Save
	save('results_200wc/T_picop_Muis_lowemission_DA.mat','T_picop_Muis_lowemission_DA');
	save('results_200wc/S_picop_Muis_lowemission_DA.mat','S_picop_Muis_lowemission_DA');
end
% }}}
% {{{ ComputeUvelAscStrengthTimeSeries:
if perform(org,'ComputeUvelAscStrengthTimeSeries'),

		%Set domain size
		Nx=690; Ny=280; Nz=70;

		%Get all melt files and store names in myfiles
		cd run_big_200wc;
		myfiles = dir(fullfile(pwd,'pickup.*.data'));
		cd ..

		%Loop over melt data (1992 - 2017)
		for i=1:1992
			%Get 3d temperature field
			filename = ['run_big_200wc/' myfiles(i).name];
			U = readbin(filename,[Nx Ny Nz],1,'real*8',0);

			%Store mean ASC strength near Totten
			Uasc_lowemission(i) = mean(U(320,135,20:45));

			%Clear variables
			clear U; clear filename;
		end

		%Save
		save('results_200wc/Uasc_lowemission.mat','Uasc_lowemission');
	end
% }}}
% {{{ MakeNewBcs_TottenHighRes:
if perform(org,'MakeNewBcs_TottenHighRes'),

		%Set OLD MITgcm mesh (690 x 280)
		Nx1=690; %number of longitude cells
		Ny1=280; %number of latitude cells
		Nz=70;  %number of cells in the vertical
		xgOrigin=92.5; %origin of longitude
		ygOrigin=-70; %origin of latitude
		dLong=0.083333; %longitude grid spacing
		dLat=0.035225; %latitude grid spacing
		lat1=(ygOrigin+dLat/2):dLat:(ygOrigin+Ny1*dLat);
		long1=(xgOrigin+dLong/2):dLong:(xgOrigin+Nx1*dLong);

		%Set MITgcm mesh (700 x 300)
		Nx=700; %number of longitude cells
		Ny=300; %number of latitude cells
		Nz=70;  %number of cells in the vertical
		xgOrigin=105; %origin of longitude
		ygOrigin=-68.5; %origin of latitude
		dLong=25/Nx; %longitude grid spacing
		dLat=5/Ny; %latitude grid spacing
		lat=(ygOrigin+dLat/2):dLat:(ygOrigin+Ny*dLat);
		long=(xgOrigin+dLong/2):dLong:(xgOrigin+Nx*dLong);

		%Set vertical grid spacing
		delZ = [10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,...
		10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82,...
		22.00, 22.00, 22.00, 22.00, 22.00, 22.00, 22.00, 22.00,...
		22.00, 22.00, 22.00, 22.00, 22.00,   22.00, 22.00, 22.00,...
		22.00, 22.00, 22.00, 22.00, 22.00,   22.00, 22.00, 22.00,...
		22.00, 22.00, 22.00, 22.00, 22.00, 30.00, 45.00, 62.72,...
		98.25, 99.25,100.01, 101.33,104.56,111.33,122.83,139.09,...
		158.94,180.83,203.55,226.50,249.50,272.50,295.50,318.50,...
		341.50,364.50,387.50,410.50,433.50,456.50]; z(1) = 0;
		for i=2:length(delZ)
			z(i)=z(i-1)-delZ(i);
		end

		%Set boundary information
		south = 183; west = 151; east = 450; n=314;

		%Make meshgrids for interpolation
		[X1 Z1_ns] = meshgrid(z,long1);
		[Y1 Z1_ew] = meshgrid(z,lat1);
		[X Z_ns] = meshgrid(z,long);
		[Y Z_ew] = meshgrid(z,lat);

		%Get all pickup files and store names
		cd run_big_200wc;
		ocnfiles = dir(fullfile(pwd,'pickup.*.data'));
		icefiles = dir(fullfile(pwd,'pickup_seaice.*.data'));
		cd ..

		%%Preallocate new matrices
		%Un = zeros(Nx,Nz,n); Vn = zeros(Nx,Nz,n); Tn = zeros(Nx,Nz,n); Sn = zeros(Nx,Nz,n);
		Us = zeros(Nx,Nz,n); Vs = zeros(Nx,Nz,n); Ts = zeros(Nx,Nz,n); Ss = zeros(Nx,Nz,n);
		%Uw = zeros(Ny,Nz,n); Vw = zeros(Ny,Nz,n); Tw = zeros(Ny,Nz,n); Sw = zeros(Ny,Nz,n);
		%Ue = zeros(Ny,Nz,n); Ve = zeros(Ny,Nz,n); Te = zeros(Ny,Nz,n); Se = zeros(Ny,Nz,n);
		%Arean = zeros(Nx,n); Heffn = zeros(Nx,n); Snown = zeros(Nx,n); Saltn = zeros(Nx,n); Uicen = zeros(Nx,n); Vicen = zeros(Nx,n);
		Areas = zeros(Nx,n); Heffs = zeros(Nx,n); Snows = zeros(Nx,n); Salts = zeros(Nx,n); Uices = zeros(Nx,n); Vices = zeros(Nx,n);
		%Areaw = zeros(Ny,n); Heffw = zeros(Ny,n); Snoww = zeros(Ny,n); Saltw = zeros(Ny,n); Uicew = zeros(Ny,n); Vicew = zeros(Ny,n);
		%Areae = zeros(Ny,n); Heffe = zeros(Ny,n); Snowe = zeros(Ny,n); Salte = zeros(Ny,n); Uicee = zeros(Ny,n); Vicee = zeros(Ny,n);

		%Loop over melt data (1992 - 2017)
		for i=1:2:n*2

			%Progress message
			t = (i+1)/2;
			disp(['Interpolating file ' num2str(t) ' of ' num2str(n)])

			%Get file names
			filename_ocn = ['run_big_200wc/' ocnfiles(i).name];
			filename_ice = ['run_big_200wc/' icefiles(i).name];

			%Get ocean data
			U1=readbin(filename_ocn,[Nx1 Ny1 Nz],1,'real*8',0);
			V1=readbin(filename_ocn,[Nx1 Ny1 Nz],1,'real*8',1);
			T1=readbin(filename_ocn,[Nx1 Ny1 Nz],1,'real*8',2);
			S1=readbin(filename_ocn,[Nx1 Ny1 Nz],1,'real*8',3);
			E1=readbin(filename_ocn,[Nx1 Ny1],1,'real*8',8);

			%Read sea ice data
			Area1=readbin(filename_ice,[Nx1 Ny1],1,'float64',7);
			Heff1=readbin(filename_ice,[Nx1 Ny1],1,'float64',8);
			Snow1=readbin(filename_ice,[Nx1 Ny1],1,'float64',9);
			alt1=readbin(filename_ice,[Nx1 Ny1],1,'float64',10);
			Uice1=readbin(filename_ice,[Nx1 Ny1],1,'float64',11);
			Vice1=readbin(filename_ice,[Nx1 Ny1],1,'float64',12);

			%Interp ocean north
			Un(:,:,t) = interp2(X1,Z1_ns,squeeze(U1(:,183,:)),X,Z_ns);
			Vn(:,:,t) = interp2(X1,Z1_ns,squeeze(V1(:,183,:)),X,Z_ns);
			Tn(:,:,t) = interp2(X1,Z1_ns,squeeze(T1(:,183,:)),X,Z_ns);
			Sn(:,:,t) = interp2(X1,Z1_ns,squeeze(S1(:,183,:)),X,Z_ns);

			%Interp ocean east/west
			Uw(:,:,t) = interp2(Y1,Z1_ew,squeeze(U1(151,:,:)),Y,Z_ew);
			Vw(:,:,t) = interp2(Y1,Z1_ew,squeeze(V1(151,:,:)),Y,Z_ew);
			Tw(:,:,t) = interp2(Y1,Z1_ew,squeeze(T1(151,:,:)),Y,Z_ew);
			Sw(:,:,t) = interp2(Y1,Z1_ew,squeeze(S1(151,:,:)),Y,Z_ew);
			Ue(:,:,t) = interp2(Y1,Z1_ew,squeeze(U1(450,:,:)),Y,Z_ew);
			Ve(:,:,t) = interp2(Y1,Z1_ew,squeeze(V1(450,:,:)),Y,Z_ew);
			Te(:,:,t) = interp2(Y1,Z1_ew,squeeze(T1(450,:,:)),Y,Z_ew);
			Se(:,:,t) = interp2(Y1,Z1_ew,squeeze(S1(450,:,:)),Y,Z_ew);

			%Interp sea ice north
			Arean(:,t) = interp1(long1,squeeze(Area1(:,183)),long);
			Heffn(:,t) = interp1(long1,squeeze(Heff1(:,183)),long);
			Snown(:,t) = interp1(long1,squeeze(Snow1(:,183)),long);
			Saltn(:,t) = interp1(long1,squeeze(Salt1(:,183)),long);
			Uicen(:,t) = interp1(long1,squeeze(Uice1(:,183)),long);
			Vicen(:,t) = interp1(long1,squeeze(Vice1(:,183)),long);

			%Interp sea ice east/west
			Areaw(:,t) = interp1(lat1,squeeze(Area1(151,:)),lat);
			Heffw(:,t) = interp1(lat1,squeeze(Heff1(151,:)),lat);
			Snoww(:,t) = interp1(lat1,squeeze(Snow1(151,:)),lat);
			Saltw(:,t) = interp1(lat1,squeeze(Salt1(151,:)),lat);
			Uicew(:,t) = interp1(lat1,squeeze(Uice1(151,:)),lat);
			Vicew(:,t) = interp1(lat1,squeeze(Vice1(151,:)),lat);
			Areae(:,t) = interp1(lat1,squeeze(Area1(450,:)),lat);
			Heffe(:,t) = interp1(lat1,squeeze(Heff1(450,:)),lat);
			Snowe(:,t) = interp1(lat1,squeeze(Snow1(450,:)),lat);
			Salte(:,t) = interp1(lat1,squeeze(Salt1(450,:)),lat);
			Uicee(:,t) = interp1(lat1,squeeze(Uice1(450,:)),lat);
			Vicee(:,t) = interp1(lat1,squeeze(Vice1(450,:)),lat);

			%Clear variables
			clear filename_ocn; clear filename_ice; clear t;
		end

		%Save BC files
		writebin('BcFiles/FOBNu_m.bin',Un);
		writebin('BcFiles/FOBNv_m.bin',Vn);
		writebin('BcFiles/FOBNs_m.bin',Sn);
		writebin('BcFiles/FOBNt_m.bin',Tn);
		writebin('BcFiles/FOBSu_m.bin',Us);
		writebin('BcFiles/FOBSv_m.bin',Vs);
		writebin('BcFiles/FOBSs_m.bin',Ss);
		writebin('BcFiles/FOBSt_m.bin',Ts);
		writebin('BcFiles/FOBEu_m.bin',Ue);
		writebin('BcFiles/FOBEv_m.bin',Ve);
		writebin('BcFiles/FOBEs_m.bin',Se);
		writebin('BcFiles/FOBEt_m.bin',Te);
		writebin('BcFiles/FOBWu_m.bin',Uw);
		writebin('BcFiles/FOBWv_m.bin',Vw);
		writebin('BcFiles/FOBWs_m.bin',Sw);
		writebin('BcFiles/FOBWt_m.bin',Tw);
		writebin('BcFiles/OBNa_m.bin',   Arean);
		writebin('BcFiles/OBNh_m.bin',   Heffn);
		writebin('BcFiles/OBNsn_m.bin',  Snown);
		writebin('BcFiles/OBNsl_m.bin',  Saltn);
		writebin('BcFiles/OBNuice_m.bin',Uicen);
		writebin('BcFiles/OBNvice_m.bin',Vicen);
		writebin('BcFiles/OBSa_m.bin',   Areas);
		writebin('BcFiles/OBSh_m.bin',   Heffs);
		writebin('BcFiles/OBSsn_m.bin',  Snows);
		writebin('BcFiles/OBSsl_m.bin',  Salts);
		writebin('BcFiles/OBSuice_m.bin',Uices);
		writebin('BcFiles/OBSvice_m.bin',Vices);
		writebin('BcFiles/OBEa_m.bin',   Areae);
		writebin('BcFiles/OBEh_m.bin',   Heffe);
		writebin('BcFiles/OBEsn_m.bin',  Snowe);
		writebin('BcFiles/OBEsl_m.bin',  Salte);
		writebin('BcFiles/OBEuice_m.bin',Uicee);
		writebin('BcFiles/OBEvice_m.bin',Vicee);
		writebin('BcFiles/OBWa_m.bin',   Areaw);
		writebin('BcFiles/OBWh_m.bin',   Heffw);
		writebin('BcFiles/OBWsn_m.bin',  Snoww);
		writebin('BcFiles/OBWsl_m.bin',  Saltw);
		writebin('BcFiles/OBWuice_m.bin',Uicew);
		writebin('BcFiles/OBWvice_m.bin',Vicew);
	end
% }}}
% {{{ ChangeFilenames:
if perform(org,'ChangeFilenames'),

	%Load parameters
	loaddata(org,'Parameters');
   
	%Copy files to change start date from MITgcm directory
	copyfile([mitgcm_dir '/data.exf'],'.');
	copyfile([mitgcm_dir '/data.cal'],'.');

	%Load data.exf lines into cell A
	fidi = fopen('data.exf','r');
	tline = fgetl(fidi);
	A={}; B={};
	i = 1; A{i} = tline;
	while ischar(tline)
		i = i+1;
		tline = fgetl(fidi);
		A{i} = tline;
	end
	fclose(fidi);

	%Change start dates of exf
	A{39} = [' atempstartdate1   = ' dates(2,:) ','];
	A{43} = [' aqhstartdate1   = ' dates(2,:) ','];
	A{47} = [' precipstartdate1   = ' dates(2,:) ','];
	A{54} = [' uwindstartdate1   = ' dates(2,:) ','];
	A{58} = [' vwindstartdate1   = ' dates(2,:) ','];
	A{62} = [' swdownstartdate1   = ' dates(2,:) ','];
	A{66} = [' lwdownstartdate1   = ' dates(2,:) ','];

	%Change start dates of obsc
	A{151} = [' obcsNstartdate1   = ' dates(2,:) ','];
	A{155} = [' obcsSstartdate1   = ' dates(2,:) ','];
	A{159} = [' obcsEstartdate1   = ' dates(2,:) ','];
	A{163} = [' obcsWstartdate1   = ' dates(2,:) ','];
	A{167} = [' siobNstartdate1   = ' dates(2,:) ','];
	A{171} = [' siobSstartdate1   = ' dates(2,:) ','];
	A{175} = [' siobEstartdate1   = ' dates(2,:) ','];
	A{179} = [' siobWstartdate1   = ' dates(2,:) ','];

	%Write new data.exf file
	fido = fopen('run/data.exf','w');
	for i=1:numel(A)
		if A{i+1} == -1
			fprintf(fido,'%s',A{i});
			break
		else
			fprintf(fido,'%s\n',A{i});
		end
	end
	fclose(fido);

	%Load data.cal lines into cell A
	fidi = fopen('data.cal','r');
	tline = fgetl(fidi);
	i = 1; B{i} = tline;
	while ischar(tline)
		i = i+1;
		tline = fgetl(fidi);
		B{i} = tline;
	end
	fclose(fidi);

	%Change calandar start date
	B{7} = [' startDate_1=' dates(2,:) ','];

	%Write new data.cal file
	fido = fopen('run/data.cal','w');
	for i=1:numel(B)
		if B{i+1} == -1
			fprintf(fido,'%s',B{i});
			break
		else
			fprintf(fido,'%s\n',B{i});
		end
	end
	fclose(fido);
	
	end
%}}}
% {{{ ChangeFilenames2:
if perform(org,'ChangeFilenames2'),

	%Load parameters
	loaddata(org,'Parameters');
   
	%%Copy files to change start date from MITgcm directory
	%copyfile([mitgcm_dir '/data.obcs'],'.');
	%copyfile([mitgcm_dir '/data.cal'],'.');

	%Load data.exf lines into cell A
	fidi = fopen('data.obcs','r');
	tline = fgetl(fidi);
	A={}; B={};
	i = 1; A{i} = tline;
	while ischar(tline)
		i = i+1;
		tline = fgetl(fidi);
		A{i} = tline;
	end
	fclose(fidi);

	%Set BC file based on timestep number
	i=629
	if i<(628*1)+1
		num = 1;
	elseif i>=(628*1)+1 & i<(628*2)+1
		num = 2;
	elseif i>=(628*2)+1 & i<(628*3)+1
		num=3;
	else
		num=4;
	end

	%Change file names of BCs based on timestep
	A{16} = [' OBNsFile=''FOBNs_m_' num2str(num) '.bin'''];
	A{17} = [' OBNtFile=''FOBNt_m_' num2str(num) '.bin'''];
	A{26} = [' OBEsFile=''FOBEs_m_' num2str(num) '.bin'''];
	A{27} = [' OBEtFile=''FOBEt_m_' num2str(num) '.bin'''];
	A{31} = [' OBWsFile=''FOBWs_m_' num2str(num) '.bin'''];
	A{32} = [' OBWtFile=''FOBWt_m_' num2str(num) '.bin'''];

	%Write new data.exf file
	fido = fopen('data2.exf','w');
	for i=1:numel(A)
		if A{i+1} == -1
			fprintf(fido,'%s',A{i});
			break
		else
			fprintf(fido,'%s\n',A{i});
		end
	end
	fclose(fido);
	error

	%Load data.cal lines into cell A
	fidi = fopen('data.cal','r');
	tline = fgetl(fidi);
	i = 1; B{i} = tline;
	while ischar(tline)
		i = i+1;
		tline = fgetl(fidi);
		B{i} = tline;
	end
	fclose(fidi);

	%Change calandar start date
	B{7} = [' startDate_1=' dates(2,:) ','];

	%Write new data.cal file
	fido = fopen('run/data.cal','w');
	for i=1:numel(B)
		if B{i+1} == -1
			fprintf(fido,'%s',B{i});
			break
		else
			fprintf(fido,'%s\n',B{i});
		end
	end
	fclose(fido);
	
	end
%}}}
% {{{ RunUncoupledMITgcm:
if perform(org,'RunUncoupledMITgcm'),
	loaddata(org,'Parameters');
	cd(mitgcm_dir)
	!qsub run8_sandy_tracer_init.pbs
	pause('on');
	while 1
		if isfile('pickup.0000002928.data')
			copyfile('pickup.0000002928.data',[pres_dir '/run']);
			copyfile('SHICE_fwFluxtave.0000002928.data',[pres_dir '/run']);
			copyfile('hFacC.data',[pres_dir '/run']);
			break;
		else 
			pause(30); continue;
		end
	end

	%clean up folder
   !rm *.log STD* *.data *.meta
end
% }}}
% {{{ InterpBasalMelt:
if perform(org,'InterpBasalMelt'),

	%Load data and models (md=issm, mdm=MITgcm)
	loaddata(org,'Parameters');
	cd(pres_dir)
	md=model; mdm=model;
   md=loadmodel(org,'MITgcmModel'); mdm=md;
   md=loadmodel(org,'TottenModel');

	%Interpolate mitgcm melt to issm mesh
   melting_rate = readbin('run/SHICE_fwFluxtave.0000002928.data',[Nx Ny]);
	m_mit = -melting_rate(:)*y2s/rho_ice;
	m_tot = InterpFromMeshToMesh2d(mdm.mesh.elements,mdm.mesh.x,mdm.mesh.y,m_mit,md.mesh.x,md.mesh.y);
	end
% }}}
% {{{ RunUncoupledISSM:
if perform(org,'RunUncoupledISSM'),

	%Load parameters and models
	loaddata(org,'Parameters');
	cd(pres_dir)
	md=model; mdm=model;
   md=loadmodel(org,'MITgcmModel'); mdm=md;
   md=loadmodel(org,'TottenModel');
	md.inversion=m1qn3inversion(md.inversion);

	%Reset model
	md.geometry.base=md.results(end).TransientSolution(end).Base;
	md.geometry.surface=md.results(end).TransientSolution(end).Surface;
	md.geometry.thickness=md.geometry.surface-md.geometry.base;
	md.initialization.vx=md.results(end).TransientSolution(end).Vx;
	md.initialization.vy=md.results(end).TransientSolution(end).Vy;
	md.initialization.vel=md.results(end).TransientSolution(end).Vel;
	md.initialization.pressure=md.results(end).TransientSolution(end).Pressure;
	md.mask.groundedice_levelset=md.results(end).TransientSolution(end).MaskGroundediceLevelset;

	%timestepping:
	%Setting Controls
	md.inversion.iscontrol       = 0;
	md.transient.ismasstransport = 1;
	md.transient.isstressbalance = 1;
	md.transient.isgroundingline = 1;
	md.transient.ismovingfront   = 0;
	md.transient.isthermal       = 0;
	md.transient.isslr           = 0;
	md.timestepping.final_time   = 1/24;
	md.timestepping.time_step    = 1/24;
	md.levelset.kill_icebergs    = 1;

	%Set Pico Parameters
	melting_rate = readbin('run/SHICE_fwFluxtave.19920101.data',[Nx Ny]);
	m_mit = -melting_rate(:)*y2s/rho_ice;
	m_tot = InterpFromMeshToMesh2d(mdm.mesh.elements,mdm.mesh.x,mdm.mesh.y,m_mit,md.mesh.x,md.mesh.y);
   md.basalforcings = basalforcings();
	md.basalforcings.floatingice_melting_rate = m_tot;
	md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);

	%Set floating ice parameters
	md.transient.requested_outputs={'default','BasalforcingsFloatingiceMeltingRate','MaskIceLevelset','MaskGroundediceLevelset'};
	pos=find(md.mesh.vertexonboundary);
	md.masstransport.spcthickness=NaN(md.mesh.numberofvertices,1);
	md.masstransport.spcthickness(pos)=md.geometry.thickness(pos);

	%SLR parameters (to satisfy model-consistency requirements, will not be used)
	md.slr.deltathickness = zeros(md.mesh.numberofelements,1);
	md.slr.sealevel       = zeros(md.mesh.numberofvertices,1);
	md.slr.spcthickness   = zeros(md.mesh.numberofvertices,1);
	md.slr.hydro_rate     = zeros(md.mesh.numberofvertices,1);
	md.slr.Ugia           = zeros(md.mesh.numberofvertices,1);
	md.slr.Ngia           = zeros(md.mesh.numberofvertices,1);

	%Solve
	md.groundingline.migration = 'SubelementMigration';
	md.groundingline.friction_interpolation='SubelementFriction1';
	md.groundingline.melt_interpolation='SubelementMelt1';
	md.cluster=generic('name',oshostname(),'np',15);
	md.verbose.solution = 1;
	md=solve(md,'tr');

	%Reset model
	md.geometry.base=md.results(end).TransientSolution(end).Base;
	md.geometry.surface=md.results(end).TransientSolution(end).Surface;
	md.geometry.thickness=md.geometry.surface-md.geometry.base;
	md.initialization.vx=md.results(end).TransientSolution(end).Vx;
	md.initialization.vy=md.results(end).TransientSolution(end).Vy;
	md.initialization.vel=md.results(end).TransientSolution(end).Vel;
	md.initialization.pressure=md.results(end).TransientSolution(end).Pressure;
	md.mask.groundedice_levelset=md.results(end).TransientSolution(end).MaskGroundediceLevelset;

	%Save
	savemodel(org,md);
end
% }}}
% {{{ Interp_Draft_S_T:
if perform(org,'Interp_Draft_S_T'),

	%Load data and models (md=issm, mdm=MITgcm)
	loaddata(org,'Parameters');
	cd(pres_dir)
	%md=model; mdm=model;
	%md=loadmodel(org,'MITgcmModel'); mdm=md;
	%md=loadmodel(org,'RunUncoupledISSM');

	%Get drafts
	old_draft_temp=readbin([mitgcm_dir '/totten_draft_360x140_BMfinal.bin'],[Nx Ny]);
	old_draft = old_draft_temp(:);
	new_draft = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(1).Base,mdm.mesh.x,mdm.mesh.y,'default',999);

	%Take new draft from issm and sub it into mitgcm draft
	draft = zeros(size(new_draft));
	pos = find(new_draft<900); pos1 = find(new_draft>=900);
	draft(pos) = new_draft(pos);
	draft(pos1) = old_draft(pos1);

	%if t>1
	% Read pickup file
	fnm='run/pickup.0000002928.data';
	U=readbin(fnm,[Nx Ny Nz],1,'real*8',0);
	V=readbin(fnm,[Nx Ny Nz],1,'real*8',1);
	T=readbin(fnm,[Nx Ny Nz],1,'real*8',2);
	S=readbin(fnm,[Nx Ny Nz],1,'real*8',3);
	E=readbin(fnm,[Nx Ny],1,'real*8',8);

	% find indices of locations where ice shelf retreated
	h=readbin('run/hFacC.data',[Nx Ny Nz]);
	msk=sum(h,3);
	msk(find(msk))=1;
	[iw jw]=find(msk); % horizontal indices where there is water
	tmp=reshape(draft,[Nx,Ny])-reshape(old_draft,[Nx Ny]);
	tmp(find(tmp<0))=0;
	[im jm]=find(tmp); % horizontal indices where there is melt

	%Extrapolate T/S to locations where ice shelf retreated
	for i=1:length(im)

		% first try vertical extrapolation
		in=find(h(im(i),jm(i),:));
		if length(in)>0;
			S(im(i),jm(i),1:min(in)  ) = S(im(i),jm(i),min(in));
			T(im(i),jm(i),1:min(in)  ) = T(im(i),jm(i),min(in));
			continue
		end

		% if not succesful, use closest neighbor horizontal extrapolation
		[y c]=min((iw-im(i)).^2+(jw-jm(i)).^2);
		salt=squeeze(S(iw(c),jw(c),:)); % salinity profile of closest neighbor
		temp=squeeze(T(iw(c),jw(c),:)); % salinity profile of closest neighbor
		in=find(h(iw(c),jw(c),:));
		salt(1:min(in))=salt(min(in));
		temp(1:min(in))=temp(min(in));
		salt(max(in):end)=salt(max(in));
		temp(max(in):end)=temp(max(in));
		S(im(i),jm(i),:)=salt;
		T(im(i),jm(i),:)=temp;
	end

	% Write initial conditions
	writebin('run/Salt.bin' ,S);
	writebin('run/Theta.bin',T);
	writebin('run/Uvel.bin' ,U);
	writebin('run/Vvel.bin' ,V);
	writebin('run/Etan.bin' ,E);
	writebin('run/totten_draft_360x140_BMfinal.bin',reshape(draft,[Nx Ny]));
	%end
end
% }}}


