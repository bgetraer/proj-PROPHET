steps=[2];

prophdir='/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET';
prefix=['PROPHET_DepthDep_'];
modeldir='./Models';
if ~isdir(modeldir)
	mkdir(modeldir);
end

org=organizer('repository',modeldir,'prefix',prefix,'steps',steps);

if perform(org,'BasalForcings') % {{{
	% load md structure
	md=loadmodel(fullfile(prophdir,'experiments/ISSM_initialization/Models/PROPHET_issm_init_TransientPrep.mat')); % md structure

	disp('   -- Define basalforcings');
	%Basal melt rate
	md.basalforcings=linearbasalforcings();
	md.basalforcings.deepwater_melting_rate=50.; % m/yr ice equivalent
	md.basalforcings.deepwater_elevation=-500;
	md.basalforcings.upperwater_melting_rate=0; % no melting for zb>=0
	md.basalforcings.upperwater_elevation=0; % sea level
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1); % no melting on grounded ice
	md.basalforcings.geothermalflux=zeros(md.mesh.numberofvertices,1);
	savemodel(org,md);
end % }}}
if perform(org,'TransientRun') % {{{
   % load model
   md=loadmodel(org,'BasalForcings');

   % set options
   disp('setting transient options');
   md.cluster=generic('name',oshostname(),'np',75);
   md.verbose.solution=true;

   % run each decade
   md.timestepping.start_time=0;
   md.timestepping.final_time=100;
   md.timestepping.time_step_max=0.05;
   % solve
   md.toolkits.DefaultAnalysis.ksp_max_it=1000;
   md.miscellaneous.name='PROPHET_test_with_depth_dep_melt';
   md=solve(md,'tr');
end % }}}
if perform(org,'VAF') % {{{	
	t=[];
	VAF=[];
	% run each decade
	decade_endyear = years(find(mod(years,10)==0));
	for d=1:5%length(decade_endyear)
		% time 
		d_time = time(time>=(decade_endyear(d)-10) &  time<=decade_endyear(d));
		disp(sprintf('YEARS: %4.0i--%4.0i',floor(d_time(1)),floor(d_time(end))));

		% load results
		fname = sprintf('./Models/%s_results_%4.0i-%4.0i',prefix,floor(d_time(1)),floor(d_time(end)));
		disp(['  loading ' fname])
		load(fname);

		t = [t results.time];
		VAF=[VAF results.IceVolumeAboveFloatation];
	end
	fname = sprintf('./Models/%s_VAF',exp_name);
	save(fname,'t','VAF');
end % }}}
