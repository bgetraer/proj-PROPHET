plotfig=1;

load('/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/data/expData.mat');
UC=load('uncoupledResults.mat');

ts=93312000:1296000:2799360000;
t = ts/3600/24/360 +2010;
%t = (1:15*24*3600:90*360*24*3600)/3600/24/360 + 2010;

if plotfig == 1 % vaf and melt


	rhoGt = 917*1e-12; % density of ice in Gt/m^3
	gt2mmslr = 1/361.8;  % 361.8 Gt of ice will raise global sea levels by ~1 mm
	figure(1);clf;
	subplot(2,1,1);
	% VAF data
	yyaxis('left');hold on;
	dmafRCP85=(vafRCP85-vafRCP85(1)).*rhoGt;
	dmafParis2C=(vafParis2C-vafParis2C(1)).*rhoGt;
	%dmafParis2C_UC=(UC.VAF_P-UC.VAF_P(1)).*rhoGt;
	%dmafRCP85_UC=(UC.VAF_R-UC.VAF_R(1)).*rhoGt;
	hRCP85 = plot(t(1:numel(vafRCP85)),dmafRCP85,'-r');
	%hRCP85UC = plot(UC.t,dmafRCP85_UC,'--r');
	hParis2C   = plot(t(1:numel(vafParis2C)),dmafParis2C,'-b');
	%hParis2CUC = plot(UC.t,dmafParis2C_UC,'--b');
	set([hRCP85,hParis2C],'linewidth',3);
%	ind = 1:min([numel(vafRCP85),numel(vafParis2C)]);
%	hdiff = plot(t(ind),vafRCP85(ind) - vafParis2C(ind),'x-');
%  jump at index 512
%	set(hdiff,'linewidth',3);
	ylabel('\Delta mass above flotation (Gt)');
	ylimleft=ylim;
	yyaxis('right');hold on;
	ylimright=ylimleft.*gt2mmslr;
	ylim(ylimright);
	ylabel('sea level rise equivalence (mm)');
	xlim([2013,2100])
	set(gca,'fontsize',14)

	subplot(2,1,2); cla;hold on;
	hRCP85 = plot(t(1:numel(meltRCP85)),meltRCP85,'r:');
	hParis2C   = plot(t(1:numel(meltParis2C)),meltParis2C,'b:');
	set([hRCP85,hParis2C],'linewidth',3);
	ylabel('total melt (m^3/yr)');

	legend([hParis2C,hRCP85],'Paris 2C','RCP 8.5')
	xlabel('time (y)')

	xlim([2013,2100])
	set(gca,'fontsize',14)
end
	
