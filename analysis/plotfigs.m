mit=loadmodel('../experiments/MITgcm_initialization/Models/PROPHET_mitgcm_init_RuntimeOptions.mat');
md=loadmodel('../experiments/ISSM_initialization/Models/PROPHET_issm_init_MeshParam.mat');
%results{1}=loadmodel('../experiments/RCP85/results/issmDiag.1420416000.mat');

fig_on=[1];

plotfig=zeros(1,5);
plotfig(fig_on)=1;

if plotfig(1) 
	%500m Temp {{{
issmresults1=loadmodel('../experiments/Paris2C/RUN02/results/issmDiag.2021760000.mat');
A=rdmds('../experiments/Paris2C/RUN02/results/pickup.save.2021760000');
% extract fields
uvel = A(:,:,(1:mit.mesh.Nz)+0*mit.mesh.Nz);
vvel = A(:,:,(1:mit.mesh.Nz)+1*mit.mesh.Nz);
temp = A(:,:,(1:mit.mesh.Nz)+2*mit.mesh.Nz);
salt = A(:,:,(1:mit.mesh.Nz)+3*mit.mesh.Nz);
% put in row,col order
uvel = permute(uvel,[2,1,3]); 
vvel = permute(vvel,[2,1,3]); 
temp = permute(temp,[2,1,3]);
salt = permute(salt,[2,1,3]);
% make zero values NaN
uvel(uvel==0)=NaN;
vvel(vvel==0)=NaN;
temp(temp==0)=NaN;
salt(salt==0)=NaN;
% interpolate horizontal velocities onto cell centers
uvel = interp3(mit.mesh.xp(1:end-1),mit.mesh.yc,mit.mesh.zc,vvel,mit.mesh.XC,mit.mesh.YC,mit.mesh.ZC);
vvel = interp3(mit.mesh.xc,mit.mesh.yp(1:end-1),mit.mesh.zc,vvel,mit.mesh.XC,mit.mesh.YC,mit.mesh.ZC);
hvel = sqrt(uvel.^2+vvel.^2);

% 500 meter depth thermal forcing field
T500=interp3(mit.mesh.xc,mit.mesh.yc,mit.mesh.zc,temp,mit.mesh.hXC,mit.mesh.hYC,-500*ones(size(mit.mesh.hXC)));
S500=interp3(mit.mesh.xc,mit.mesh.yc,mit.mesh.zc,salt,mit.mesh.hXC,mit.mesh.hYC,-500*ones(size(mit.mesh.hXC)));
S500(S500<30)=min(S500(S500>30));
rho_0=1028; % approximate value for estimating pressure (kg/m^3)
T_0=10; % deg C
S_0=35; %
alpha=1.7E-4; % K^-1
beta=7.6E-4;
rho_w=rho_0.*(1-alpha.*(T500-T_0) + beta.*(S500-S_0)); % EOS (kg/m^3)
% Calculate in situ freezing point (from Holland, Jenkins, and Holland 2008)
g = 9.81;    % m/s^2
a = -0.0573; % deg C
b = 0.0832;  % deg C
c = 7.53E-3*1E-5*rho_w.*g; % deg C/Pa
freezingpoint = a.*S500 + b + c.*-500; % deg C
TF500=T500-freezingpoint; % deg C
T500alpha = ~(isnan(T500) | mit.geometry.bathy>=-500);

% 500 meter bathymetry mask
bathy500=mit.geometry.bathy>-500;
% contours
contBathy=mit.geometry.bathy;
contBathy(bwdist(~isnan(T500))>1)=NaN;
contBathysmoothed=imgaussfilt(mit.geometry.bathy,1);
contBathysmoothed(bwdist(~isnan(T500))>1)=NaN;

% }}}
%grounding lines {{{
	years = 2075;
	modeltime=(years-2010)*3600*24*360;

	% 2013
	contoursInit=isoline(md,md.mask.ocean_levelset,'output','matrix');

	fname=sprintf('../experiments/RCP85/RUN02/results/issmDiag.%010.0f.mat',modeltime);
	disp(['loading ' fname]);
	results=loadmodel(fname);
	contoursRCP85=isoline(md, results.MaskOceanLevelset,'output','matrix');


	fname=sprintf('../experiments/Paris2C/RUN02/results/issmDiag.%010.0f.mat',modeltime);
   disp(['loading ' fname]);
   results=loadmodel(fname);
   contoursParis2C=isoline(md, results.MaskOceanLevelset,'output','matrix');
	% }}}
	%thickness change {{{
	dH_md=results.Thickness-md.geometry.thickness;
	dH=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,dH_md,mit.mesh.xc,mit.mesh.yc,NaN);
	oceanmask=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,results.MaskOceanLevelset,mit.mesh.xc,mit.mesh.yc,-1);
	icemask=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ice_levelset,mit.mesh.xc,mit.mesh.yc,1);
	base=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.geometry.base,mit.mesh.xc,mit.mesh.yc,NaN);
	% }}}
	% ice flow {{{
	[L1,n]=bwlabel(oceanmask<0,4);
	L_gr = L1~=1;

	[L2,n]=bwlabel(oceanmask>0,4);
	B = groupcounts(L2(:));
	L_gr(L2~=0 & L2~=1 & L2~=5)=0;
	L_gr(L1==0)=1;

	L_gr=bwlabel(L_gr,4)>0;

	M=contour(mit.mesh.xc,mit.mesh.yc,imgaussfilt(single(L_gr),1),[0.5,0.5],'-k');
	pt = interparc(1000,M(1,2:end),M(2,2:end));

	xedge = [mit.mesh.hXC(find(L_gr(:,1),1):end,1); ... % W
		mit.mesh.hXC(end,2:end)'; ... % N
		mit.mesh.hXC(1:end-1,end); ... % E
		mit.mesh.hXC(1,find(L_gr(1,:),1):end-1)']; % S
	yedge = [mit.mesh.hYC(find(L_gr(:,1),1):end,1); ... % W
      mit.mesh.hYC(end,2:end)'; ... % N
      mit.mesh.hYC(1:end-1,end); ... % E
      mit.mesh.hYC(1,find(L_gr(1,:),1):end-1)']; % S


	vx=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,results.Vx,mit.mesh.xc,mit.mesh.yc,0);
	vy=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,results.Vy,mit.mesh.xc,mit.mesh.yc,0);
	vertEDGE = stream2(mit.mesh.hXC,mit.mesh.hYC,vx,vy,xedge(5:20:end),yedge(5:20:end));
	vertEDGE2 = stream2(mit.mesh.hXC,mit.mesh.hYC,vx,vy,xedge,yedge);

	vertGL = stream2(mit.mesh.hXC,mit.mesh.hYC,-vx,-vy,pt(:,1),pt(:,2));
	% }}}
	figure(1); clf;
	% T500
	ax1=axes;hold on;
	imagesc(mit.mesh.xc,mit.mesh.yc,TF500,'alphadata',double(T500alpha));
	cmap=[brewermap(120,'reds')];
	cmap=[[linspace(1,cmap(1,1),21)',linspace(1,cmap(1,2),21)',linspace(1,cmap(1,3),21)'];cmap(2:end,:)];
	set(ax1,'ydir','normal','colormap',cmap,'clim',[0 3.5]); 
	axis equal tight off;
	colorbar('Location','north')

	ax2=axes;hold on;
	imagesc(mit.mesh.xc,mit.mesh.yc,dH,'alphadata',double(~T500alpha & [base<-500 | oceanmask>0]));
	set(ax2,'ydir','normal','color','none','colormap',brewermap(100,'PRGn'))
	set(ax2,'clim',[-700,700])
	axis equal tight off;
	colorbar('Location','north')

	ax3=axes;hold on;
	imagesc(mit.mesh.xc,mit.mesh.yc,zeros(size(dH)),'alphadata',double(~T500alpha & mit.geometry.bathy>-500 & icemask>0));
	set(ax3,'ydir','normal','color','none','colormap',[1/3 1/3 1/3],'layer','top')
   axis equal tight;
	colorbar('Location','north')


	[~,hc500]=contour(mit.mesh.xc,mit.mesh.yc,contBathy,        [-500 -500],  'k','linewidth',1,'EdgeAlpha',1.0);
	[~,hc750]=contour(mit.mesh.xc,mit.mesh.yc,contBathysmoothed,[-750 -750],  'k','linewidth',1,'EdgeAlpha',0.2);
	[~,hc1000]=contour(mit.mesh.xc,mit.mesh.yc,contBathysmoothed,[-1000 -1000],'k','linewidth',1,'EdgeAlpha',0.5);
	hInit_GL=plot(contoursInit(:,1),contoursInit(:,2),'-b','linewidth',2);
	hRCP85_GL=  plot(contoursRCP85(:,1),contoursRCP85(:,2),'-g','linewidth',2);
	hParis2C_GL=plot(contoursParis2C(:,1),contoursParis2C(:,2),'-m','linewidth',2);
	hSL=streamline(vertEDGE);
	set(hSL,'linewidth',1,'LineStyle','--','color','k')

	legend([hc500,hc750,hc1000,hInit_GL,hRCP85_GL,hParis2C_GL,hSL(1)],'500m depth','750m depth','1000m depth',...
		'2013 grounding line','2050 grounding line (RCP8.5)','2050 grounding line (Paris2C)','ice flow streamlines')

	%% figure 2

	load('../../proj-n4/mat/domainplot.mat');
	figure(2);clf;hold on;
	contour(mit.mesh.xc,mit.mesh.yc,icemask,[0 0])
	contour(mit.mesh.xc,mit.mesh.yc,oceanmask,[0 0])
	plot(xrock,yrock,'-g','markersize',0.001);
	%h_shelf=patch(xshelf,yshelf,[1,1,1],'linestyle','--');
	%h_gl=plot(x1,y1,'-k','linewidth',0.75);
	%h_coast=plot(x2,y2,'-k','linewidth',2);
	axis equal tight;
	set(gca,'xlim',[min(mit.mesh.xp) max(mit.mesh.xp)],'ylim',[min(mit.mesh.yp) max(mit.mesh.yp)]);
	hSL=streamline(vertEDGE2);
	hSL_gl=streamline(vertGL);
	%save('coastline_data','x_coastline','y_coastline');

	%quiver(mIit.mesh.hXC,mit.mesh.hYC,U,V)
	%D = sym/divergence(mit.mesh.hXC,mit.mesh.hYC,U,V);

%% Bathy
%ax2=axes;
%imagesc(mit.mesh.xc,mit.mesh.yc,bathy500,'alphadata',double(bathy500));
%set(ax2,'color','none','ydir','normal','colormap',flip(gray),'clim',[0 1]);
%axis equal tight off;
%
%ax1=axes;
%surf(mit.mesh.xc,mit.mesh.yc,hVEL,'EdgeColor','none');
%view(2);
%axis equal tight;
%dT=issmresults1.Thickness-md.geometry.thickness;
%dT(md.mask.ice_levelset>0)=NaN;
%dT=mean(dT(md.mesh.elements),2);
%set(ax1,'colormap',jet,'clim',[-0.5 0.5]);
%
%ax2=axes;
%A=mean(issmresults1.BasalforcingsFloatingiceMeltingRate(md.mesh.elements),2);
%patch('Faces', md.mesh.elements, 'Vertices', [md.mesh.x md.mesh.y],'CData',A,'FaceColor','flat','EdgeColor','none')
%axis equal tight;
%set([ax1,ax2],'ydir','normal');
%set(ax2,'colormap',parula,'clim',[0 100]);
end
