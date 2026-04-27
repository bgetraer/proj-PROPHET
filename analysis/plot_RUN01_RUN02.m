R1 = load('/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/data/RUN01.mat');
R2 = load('/totten_1/bgetraer/issmjpl/proj-getraer/proj-PROPHET/data/RUN02.mat');

ts=93312000:1296000:2799360000;
t = ts/3600/24/360 +2010;



rhoGt = 917*1e-12; % density of ice in Gt/m^3
gt2mmslr = 1/361.8;  % 361.8 Gt of ice will raise global sea levels by ~1 mm

figure(1);clf;
% PARIS2C
subplot(2,1,1);
% MAF
yyaxis('left');hold on;
mafParis2C_R1=(R1.vafParis2C-R1.vafParis2C(1)).*rhoGt;
hParis2C_R1   = plot(t(1:numel(R1.vafParis2C)),mafParis2C_R1,'--r');
mafParis2C_R2=(R2.vafParis2C-R2.vafParis2C(1)).*rhoGt;
hParis2C_R2   = plot(t(1:numel(R2.vafParis2C)),mafParis2C_R2,'-b');
set([hParis2C_R1,hParis2C_R2],'linewidth',3);
ylabel('\Delta mass above flotation (Gt)');
ylimleft=ylim;
% RUN 02
yyaxis('right');hold on;
ylimright=ylimleft.*gt2mmslr;
ylim(ylimright);
ylabel('sea level rise equivalence (mm)');
xlim([2013,2100])
set(gca,'fontsize',14)
title('Paris2C')
legend([hParis2C_R1,hParis2C_R2],'RUN01','RUN02')


% RCP85
subplot(2,1,2);
% MAF
yyaxis('left');hold on;
mafRCP85_R1=(R1.vafRCP85-R1.vafRCP85(1)).*rhoGt;
hRCP85_R1   = plot(t(1:numel(R1.vafRCP85)),mafRCP85_R1,'--r');
mafRCP85_R2=(R2.vafRCP85-R2.vafRCP85(1)).*rhoGt;
hRCP85_R2   = plot(t(1:numel(R2.vafRCP85)),mafRCP85_R2,'-b');
set([hRCP85_R1,hRCP85_R2],'linewidth',3);
ylabel('\Delta mass above flotation (Gt)');
ylimleft=ylim;
% RUN 02
yyaxis('right');hold on;
ylimright=ylimleft.*gt2mmslr;
ylim(ylimright);
ylabel('sea level rise equivalence (mm)');
xlim([2013,2100])
set(gca,'fontsize',14)
title('RCP85')
legend([hRCP85_R1,hRCP85_R2],'RUN01','RUN02')



figure(2);clf;
% PARIS2C
subplot(2,1,1);
% MAF
yyaxis('left');hold on;
dmafParis2C_R1=diff(R1.vafParis2C-R1.vafParis2C(1)).*rhoGt;
hParis2C_R1   = plot(t(1:numel(dmafParis2C_R1)),dmafParis2C_R1,'--r');
dmafParis2C_R2=diff(R2.vafParis2C-R2.vafParis2C(1)).*rhoGt;
hParis2C_R2   = plot(t(1:numel(dmafParis2C_R2)),dmafParis2C_R2,'-b');
set([hParis2C_R1,hParis2C_R2],'linewidth',3);
ylabel('\Delta mass above flotation (Gt)');
ylimleft=ylim;
% RUN 02
yyaxis('right');hold on;
ylimright=ylimleft.*gt2mmslr;
ylim(ylimright);
ylabel('sea level rise equivalence (mm)');
xlim([2013,2100])
set(gca,'fontsize',14)
title('Paris2C')
legend([hParis2C_R1,hParis2C_R2],'RUN01','RUN02')


% RCP85
subplot(2,1,2);
% MAF
yyaxis('left');hold on;
dmafRCP85_R1=diff(R1.vafRCP85-R1.vafRCP85(1)).*rhoGt;
hRCP85_R1   = plot(t(1:numel(dmafRCP85_R1)),dmafRCP85_R1,'--r');
dmafRCP85_R2=diff(R2.vafRCP85-R2.vafRCP85(1)).*rhoGt;
hRCP85_R2   = plot(t(1:numel(dmafRCP85_R2)),dmafRCP85_R2,'-b');
set([hRCP85_R1,hRCP85_R2],'linewidth',3);
ylabel('\Delta mass above flotation (Gt)');
ylimleft=ylim;
% RUN 02
yyaxis('right');hold on;
ylimright=ylimleft.*gt2mmslr;
ylim(ylimright);
ylabel('sea level rise equivalence (mm)');
xlim([2013,2100])
set(gca,'fontsize',14)
title('RCP85')
legend([hRCP85_R1,hRCP85_R2],'RUN01','RUN02')
