steps=[4];

% ISSM alteration:
% file		src/c/classes/Elements/Element.cpp
% line		2600
% change		add code
%				if(base[i]>100.){
%					values[i]=0.;
%				}
% notes		alter behavior of linearbasalforcings() to set melt rate to 0 on floating ice with a draft less than 100 m. per Dan Goldberg requested control test.

prophdir='/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET';
prefix=['PROPHET_DepthDep_'];
modeldir='/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/Depth_dependent_melt/Models';

org=organizer('repository',modeldir,'prefix',prefix,'steps',steps);

if perform(org,'BasalForcings') % {{{
	% load md structure
	md=loadmodel(fullfile(prophdir,'experiments/ISSM_initialization/Models/PROPHET_issm_init_TransientPrep.mat')); % md structure

	disp('   -- Define basalforcings');
	%Basal melt rate

	md.basalforcings=linearbasalforcings();
	md.basalforcings.deepwater_melting_rate=60.; % m/yr ice equivalent
	md.basalforcings.deepwater_elevation=-700.;
	md.basalforcings.upperwater_melting_rate=0; % no melting for zb>=0
	md.basalforcings.upperwater_elevation=-100; % sea level
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1); % no melting on grounded ice
	md.basalforcings.geothermalflux=zeros(md.mesh.numberofvertices,1);
	savemodel(org,md);
end % }}}
if perform(org,'TransientRun') % {{{
	% load model
	md=loadmodel(org,'BasalForcings');

	% set options
	disp('setting transient options');
	md.cluster=generic('name',oshostname(),'np',90);
	md.verbose.solution=true;

	% run
	md.timestepping=timesteppingadaptive();
	md.timestepping.start_time=2012;
	md.timestepping.final_time=2100;
	md.timestepping.time_step_min=0.05;
	md.timestepping.time_step_max=0.50;
	md.settings.output_frequency = 2;
	md.cluster.interactive=0;
	% solve
	md.miscellaneous.name='PROPHET_control_with_depth_dep_melt';
	savemodel(org,md);
	%md=solve(md,'tr','loadonly',1);
end % }}}
if perform(org,'Results') % {{{
	disp('loading results from TransientRun');
	md=loadmodel(org,'TransientRun');
	filename='/totten_1/bgetraer/ISSM//execution/PROPHET_control_with_depth_dep_melt-10-27-2025-19-56-39-700988/PROPHET_control_with_depth_dep_melt.outbin';
	md=loadresultsfromdisk(md,filename);	
	savemodel(org,md);
end % }}}
if perform(org,'Plot') % {{{
	md=loadmodel(org,'Results') 
	oceanarea = 3.618E14; % m^2
	rho_ice = 917; % kg/m^3
	rho_freshwater = 1000; % kg/m^3

	t   = [md.results.TransientSolution.time];
	vaf = [md.results.TransientSolution.IceVolumeAboveFloatation];

	vaf_anom = vaf-vaf(1);

	maf_anom = vaf_anom*rho_ice; % kg
	sle_anom = maf_anom/rho_freshwater/oceanarea; % Sea level equivalence (m)

	yyaxis('left');
	plot(t,maf_anom);
	ylabel('\DeltaMAF (10^3 Gt)');
	ytick = get(gca,'ytick');
	set(gca,'yticklabels',ytick*1E-12*1E-3);
	yyaxis('right');
	ylabel('SLE (mm)');
	set(gca,'ytick',ytick,'yticklabels',round(-ytick/rho_freshwater/oceanarea*1E3),...
		'ylim',[ytick(1),ytick(end)]);
	
	title('Total ice loss (ASE)');
	xlabel('years');
end % }}}
