cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916
addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/retrieval/')
addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/')
%
date = '150916';
pathname = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/';
pathdate = [date,'/2015-09-16-landsat8-svc-spectra/'];

rho = 0.028;

%% IBAYN

% 150916_1225_R053_T054.sig dummy
% 150916_1226_R055_T056.sig water
% 150916_1226_R055_T057.sig water
% 150916_1227_R055_T058.sig water no picture!!!
% 150916_1227_R055_T059.sig sky
% 150916_1228_R055_T060.sig sky
% 150916_1228_R055_T061.sig dummy
% 150916_1306_R055_T062.sig dummy
% 150916_1307_R055_T063.sig dummy



filename = '150916_1226_R055_T056.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = '150916_1226_R055_T057.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = '150916_1227_R055_T058.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = '150916_1227_R055_T059.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskya,~] = extractSVC(filepath);

filename = '150916_1228_R055_T060.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskyb,~] = extractSVC(filepath);

Rrs1a = (Lt1-rho.*Lskya)./(pi.*Lg./0.99);
Rrs2a = (Lt2-rho.*Lskya)./(pi.*Lg./0.99);
Rrs3a = (Lt3-rho.*Lskya)./(pi.*Lg./0.99);

Rrs1b = (Lt1-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs2b = (Lt2-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs3b = (Lt3-rho.*Lskyb)./(pi.*Lg./0.99);



figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Rrs1a,'r')
hold on
plot(wavelengthSVC,Rrs2a,'b')
plot(wavelengthSVC,Rrs3a,'g')
plot(wavelengthSVC,Rrs1b,'--m')
plot(wavelengthSVC,Rrs2b,'--c')
plot(wavelengthSVC,Rrs3b,'--k')
legend('Lt1 Lskya','Lt2 Lskya','Lt3 Lskya',...
    'Lt1 Lskyb','Lt2 Lskyb','Lt3 Lskyb')
% plot(wavelengthSVC,r./100)

title('R_{rs} -- IBAYN','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([400 1000 0 0.03])
grid on

%% IBAYN average
RrsIBAYN150916=mean([...
Rrs1a,...
Rrs2a,...
Rrs3a,...
Rrs1b,...
Rrs2b,...
Rrs3b],2);

hold on
plot(wavelengthSVC,RrsIBAYN150916,'r','LineWidth',1.5)

%% Radiance Lsky
figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Lskya,'r')
hold on
plot(wavelengthSVC,Lskyb,'b')

legend('Lskya','Lskyb')

title('Lsky -- IBAYN','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.015])
grid on



%% ONTOS 1
% 150916_1308_R064_T065.sig water
% 150916_1308_R064_T066.sig water
% 150916_1308_R064_T067.sig water
% 150916_1309_R064_T068.sig sky
% 150916_1309_R064_T069.sig sky
% 150916_1311_R064_T070.sig dummy


filename = '150916_1308_R064_T065.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = '150916_1308_R064_T066.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = '150916_1308_R064_T067.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = '150916_1309_R064_T068.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskya,~] = extractSVC(filepath);

filename = '150916_1309_R064_T069.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskyb,~] = extractSVC(filepath);

Rrs1a_1 = (Lt1-rho.*Lskya)./(pi.*Lg./0.99);
Rrs2a_1 = (Lt2-rho.*Lskya)./(pi.*Lg./0.99);
Rrs3a_1 = (Lt3-rho.*Lskya)./(pi.*Lg./0.99);

Rrs1b_1 = (Lt1-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs2b_1 = (Lt2-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs3b_1 = (Lt3-rho.*Lskyb)./(pi.*Lg./0.99);

figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Rrs1a_1,'r')
hold on
plot(wavelengthSVC,Rrs2a_1,'b')
plot(wavelengthSVC,Rrs3a_1,'g')
plot(wavelengthSVC,Rrs1b_1,'--m')
plot(wavelengthSVC,Rrs2b_1,'--c')
plot(wavelengthSVC,Rrs3b_1,'--k')
legend('Lt1 Lskya','Lt2 Lskya','Lt3 Lskya',...
    'Lt1 Lskyb','Lt2 Lskyb','Lt3 Lskyb')
% plot(wavelengthSVC,r./100)

title('R_{rs} -- ONTOS 1','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([400 1000 0 0.03])
grid on

% Radiance Lsky
% figure
% fs = 15;
% set(gcf,'color','white')
% plot(wavelengthSVC,Lskya,'r')
% hold on
% plot(wavelengthSVC,Lskyb,'b')
% 
% legend('Lskya','Lskyb')
% 
% title('Lsky -- ONTOS 1','fontsize',fs)
% xlabel('wavelengthSVC [nm]','fontsize',fs)
% ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
% set(gca,'fontsize',fs)
% % axis([400 1000 0 0.015])
% grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ONTOS 2
% 150916_1312_R071_T072.sig water
% 150916_1312_R071_T073.sig water
% 150916_1313_R071_T074.sig water
% 150916_1313_R071_T075.sig sky
% 150916_1313_R071_T076.sig sky



filename = '150916_1312_R071_T072.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = '150916_1312_R071_T073.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = '150916_1313_R071_T074.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = '150916_1313_R071_T075.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskya,~] = extractSVC(filepath);

filename = '150916_1313_R071_T076.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskyb,~] = extractSVC(filepath);

Rrs1a_2 = (Lt1-rho.*Lskya)./(pi.*Lg./0.99);
Rrs2a_2 = (Lt2-rho.*Lskya)./(pi.*Lg./0.99);
Rrs3a_2 = (Lt3-rho.*Lskya)./(pi.*Lg./0.99);

Rrs1b_2 = (Lt1-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs2b_2 = (Lt2-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs3b_2 = (Lt3-rho.*Lskyb)./(pi.*Lg./0.99);


% figure
% fs = 15;
% set(gcf,'color','white')
plot(wavelengthSVC,Rrs1a_2,'r')
% hold on
plot(wavelengthSVC,Rrs2a_2,'b')
plot(wavelengthSVC,Rrs3a_2,'g')
plot(wavelengthSVC,Rrs1b_2,'--m')
plot(wavelengthSVC,Rrs2b_2,'--c')
plot(wavelengthSVC,Rrs3b_2,'--k')
legend('Lt1 Lskya','Lt2 Lskya','Lt3 Lskya',...
    'Lt1 Lskyb','Lt2 Lskyb','Lt3 Lskyb')
% plot(wavelengthSVC,r./100)

title('R_{rs} -- ONTOS 2','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([400 1000 0 0.03])
grid on

% %% Radiance Lsky
% figure
% fs = 15;
% set(gcf,'color','white')
% plot(wavelengthSVC,Lskya,'r')
% hold on
% plot(wavelengthSVC,Lskyb,'b')
% 
% legend('Lskya','Lskyb')
% 
% title('Lsky -- ONTOS 2','fontsize',fs)
% xlabel('wavelengthSVC [nm]','fontsize',fs)
% ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
% set(gca,'fontsize',fs)
% % axis([400 1000 0 0.015])
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% ONTOS 3
% 150916_1315_R077_T078.sig water
% 150916_1315_R077_T079.sig water
% 150916_1316_R077_T080.sig water
% 150916_1316_R077_T081.sig sky
% 150916_1316_R077_T082.sig sky



filename = '150916_1315_R077_T078.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = '150916_1315_R077_T079.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = '150916_1316_R077_T080.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = '150916_1316_R077_T081.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskya,~] = extractSVC(filepath);

filename = '150916_1316_R077_T082.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskyb,~] = extractSVC(filepath);

Rrs1a_3 = (Lt1-rho.*Lskya)./(pi.*Lg./0.99);
Rrs2a_3 = (Lt2-rho.*Lskya)./(pi.*Lg./0.99);
Rrs3a_3 = (Lt3-rho.*Lskya)./(pi.*Lg./0.99);

Rrs1b_3 = (Lt1-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs2b_3 = (Lt2-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs3b_3 = (Lt3-rho.*Lskyb)./(pi.*Lg./0.99);


% figure
% fs = 15,
% set(gcf,'color','white')
plot(wavelengthSVC,Rrs1a_3,'r')
% hold on
plot(wavelengthSVC,Rrs2a_3,'b')
plot(wavelengthSVC,Rrs3a_3,'g')
plot(wavelengthSVC,Rrs1b_3,'--m')
plot(wavelengthSVC,Rrs2b_3,'--c')
plot(wavelengthSVC,Rrs3b_3,'--k')
legend('Lt1 Lskya','Lt2 Lskya','Lt3 Lskya',...
    'Lt1 Lskyb','Lt2 Lskyb','Lt3 Lskyb')
% plot(wavelengthSVC,r./100)

title('R_{rs} -- ONTOS 3','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([400 1000 0 0.03])
grid on

% %% Radiance Lsky
% figure
% fs = 15;
% set(gcf,'color','white')
% plot(wavelengthSVC,Lskya,'r')
% hold on
% plot(wavelengthSVC,Lskyb,'b')
% 
% legend('Lskya','Lskyb')
% 
% title('Lsky -- ONTOS 3','fontsize',fs)
% xlabel('wavelengthSVC [nm]','fontsize',fs)
% ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
% set(gca,'fontsize',fs)
% % axis([400 1000 0 0.015])
% grid on
% 
% %% ONTOS average
RrsONTOS150916=mean([...
Rrs1a_1,...
Rrs2a_1,...
Rrs3a_1,...
Rrs1b_1,...
Rrs2b_1,...
Rrs3b_1,...
Rrs1a_2,...
Rrs2a_2,...
Rrs3a_2,...
Rrs1b_2,...
Rrs2b_2,...
Rrs3b_2,...
Rrs1a_3,...
Rrs2a_3,...
Rrs3a_3,...
Rrs1b_3,...
Rrs2b_3,...
Rrs3b_3],2);

hold on
plot(wavelengthSVC,RrsONTOS150916,'r','LineWidth',1.5)

clearvars -except wavelengthSVC RrsONTOS150916
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RVRPLM
% 150916_1345_R083_T084.sig water
% 150916_1345_R083_T085.sig water
% 150916_1345_R083_T086.sig water
% 150916_1346_R083_T087.sig sky
% 150916_1346_R083_T088.sig sky



filename = '150916_1345_R083_T084.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = '150916_1345_R083_T085.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = '150916_1345_R083_T086.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = '150916_1346_R083_T087.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskya,~] = extractSVC(filepath);

filename = '150916_1346_R083_T088.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskyb,~] = extractSVC(filepath);

Rrs1a = (Lt1-rho.*Lskya)./(pi.*Lg./0.99);
Rrs2a = (Lt2-rho.*Lskya)./(pi.*Lg./0.99);
Rrs3a = (Lt3-rho.*Lskya)./(pi.*Lg./0.99);

Rrs1b = (Lt1-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs2b = (Lt2-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs3b = (Lt3-rho.*Lskyb)./(pi.*Lg./0.99);

% RrsRVRPLM150916 = Rrs1a;

figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Rrs1a,'r')
hold on
plot(wavelengthSVC,Rrs2a,'b')
plot(wavelengthSVC,Rrs3a,'g')
plot(wavelengthSVC,Rrs1b,'--m')
plot(wavelengthSVC,Rrs2b,'--c')
plot(wavelengthSVC,Rrs3b,'--k')
legend('Lt1 Lskya','Lt2 Lskya','Lt3 Lskya',...
    'Lt1 Lskyb','Lt2 Lskyb','Lt3 Lskyb')
% plot(wavelengthSVC,r./100)

title('R_{rs} -- RVRPLM','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([400 1000 0 0.03])
grid on

%% RVRPL average
RrsRVRPL150916=mean([...
Rrs1a,...
Rrs2a,...
Rrs3a,...
Rrs1b,...
Rrs2b,...
Rrs3b],2);

hold on
plot(wavelengthSVC,RrsRVRPL150916,'r','LineWidth',1.5)

%% Radiance Lsky
figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Lskya,'r')
hold on
plot(wavelengthSVC,Lskyb,'b')

legend('Lskya','Lskyb')

title('Lsky -- RVRPLM','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.015])
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% CRANB
% 150916_1509_R093_T094.sig water
% 150916_1509_R093_T095.sig water
% 150916_1509_R093_T096.sig water
% 150916_1510_R093_T097.sig sky
% 150916_1510_R093_T098.sig sky
% 150916_1510_R093_T099.sig sky


filename = '150916_1509_R093_T094.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = '150916_1509_R093_T095.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = '150916_1509_R093_T096.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = '150916_1510_R093_T097.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskya,~] = extractSVC(filepath);

filename = '150916_1510_R093_T098.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskyb,~] = extractSVC(filepath);

filename = '150916_1510_R093_T099.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskyc,~] = extractSVC(filepath);

Rrs1a = (Lt1-rho.*Lskya)./(pi.*Lg./0.99);
Rrs2a = (Lt2-rho.*Lskya)./(pi.*Lg./0.99);
Rrs3a = (Lt3-rho.*Lskya)./(pi.*Lg./0.99);

Rrs1b = (Lt1-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs2b = (Lt2-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs3b = (Lt3-rho.*Lskyb)./(pi.*Lg./0.99);

Rrs1c = (Lt1-rho.*Lskyc)./(pi.*Lg./0.99);
Rrs2c = (Lt2-rho.*Lskyc)./(pi.*Lg./0.99);
Rrs3c = (Lt3-rho.*Lskyc)./(pi.*Lg./0.99);

% RrsLONGS150916 = Rrs2a;

figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Rrs1a,'r')
hold on
plot(wavelengthSVC,Rrs2a,'b')
plot(wavelengthSVC,Rrs3a,'g')
plot(wavelengthSVC,Rrs1b,'--m')
plot(wavelengthSVC,Rrs2b,'--c')
plot(wavelengthSVC,Rrs3b,'--k')
plot(wavelengthSVC,Rrs1c,'--r')
plot(wavelengthSVC,Rrs2c,'--b')
plot(wavelengthSVC,Rrs3c,'--g')
legend('Lt1 Lskya','Lt2 Lskya','Lt3 Lskya',...
    'Lt1 Lskyb','Lt2 Lskyb','Lt3 Lskyb',...
    'Lt1 Lskyc','Lt2 Lskyc','Lt3 Lskyc')

title('R_{rs} -- CRANB','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.03])
grid on
%%
% Ref = [wavelengthSVC*1E-3, RrsLONGS150916];
% save('LONGSRef_140919.txt','Ref','-ascii')

%% CRANB average
RrsCRANB150916=mean([...
Rrs1a,...
Rrs2a,...
Rrs3a,...
Rrs1b,...
Rrs2b,...
Rrs3b,...
Rrs1c,...
Rrs2c,...
Rrs3c],2);

hold on
plot(wavelengthSVC,RrsCRANB150916,'r','LineWidth',1.5)

%% Radiance Lsky
figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Lskya,'r')
hold on
plot(wavelengthSVC,Lskyb,'b')
plot(wavelengthSVC,Lskyc,'g')
legend('Lskya','Lskyb','Lskyc')

title('Lsky -- CRANB','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.015])
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LONGN
% 150916_1542_R101_T102.sig water
% 150916_1542_R101_T104.sig water
% 150916_1543_R101_T105.sig sky
% 150916_1543_R101_T106.sig sky
% 150916_1543_R101_T107.sig sky


filename = '150916_1542_R101_T102.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = '150916_1542_R101_T104.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = '150916_1543_R101_T105.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskya,~] = extractSVC(filepath);

filename = '150916_1543_R101_T106.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskyb,~] = extractSVC(filepath);

filename = '150916_1543_R101_T107.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskyc,~] = extractSVC(filepath);


Rrs1a = (Lt1-rho.*Lskya)./(pi.*Lg./0.99);
Rrs2a = (Lt2-rho.*Lskya)./(pi.*Lg./0.99);

Rrs1b = (Lt1-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs2b = (Lt2-rho.*Lskyb)./(pi.*Lg./0.99);

Rrs1c = (Lt1-rho.*Lskyc)./(pi.*Lg./0.99);
Rrs2c = (Lt2-rho.*Lskyc)./(pi.*Lg./0.99);
% RrsLONGS150916 = Rrs2a;

figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Rrs1a,'r')
hold on
plot(wavelengthSVC,Rrs2a,'b')
plot(wavelengthSVC,Rrs1b,'g')
plot(wavelengthSVC,Rrs2b,'--m')
plot(wavelengthSVC,Rrs1c,'--c')
plot(wavelengthSVC,Rrs1c,'--k')
legend('Lt1 Lskya','Lt2 Lskya','Lt1 Lskyb',...
    'Lt2 Lskyb','Lt1 Lskyc','Lt2 Lskyc')

title('R_{rs} -- LONGN','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([400 1000 0 0.03])
grid on

%% LONGN average
RrsLONGN150916=mean([...
Rrs1a,...
Rrs2a,...
Rrs1b,...
Rrs2b,...
Rrs1c,...
Rrs2c],2);

hold on
plot(wavelengthSVC,RrsLONGN150916,'r','LineWidth',1.5)

%% Radiance Lsky
figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Lskya,'r')
hold on
plot(wavelengthSVC,Lskyb,'g')
plot(wavelengthSVC,Lskyc,'b')
legend('Lskya','Lskyb','Lskyc')

title('Lsky -- LONGN','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.015])
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Rrs all
wavelengthSVC150916 = wavelengthSVC;

figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,RrsONTOS150916,'b')
hold on
plot(wavelengthSVC,RrsIBAYN150916,'r')
plot(wavelengthSVC,RrsRVRPL150916,'g')
plot(wavelengthSVC,RrsCRANB150916,'m')
plot(wavelengthSVC,RrsLONGN150916,'k')

legend('ONTOS','IBAYN','RVRPL','CRANB','LONGN')
title('R_{rs} -- 09/16/15 ','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.07])
grid on

save Rrs150916AllSites.mat wavelengthSVC150916 ...
    RrsIBAYN150916 RrsONTOS150916 RrsRVRPL150916 ...
    RrsCRANB150916 RrsLONGN150916

%% Spectrally sampled and save in text file
wlrange = wavelengthSVC>=400 & wavelengthSVC<=2500;
% wlzero = wavelengthSVC==2219.0;
wlavg = wavelengthSVC>=860 & wavelengthSVC<=870;

zeroavg = mean(RrsONTOS150916(wlavg));

RrsONTOScorr = RrsONTOS150916-zeroavg; % subtract NIR value for all bands
% RrsONTOScorr = RrsONTOS150916;%-zeroavg; % force NIR and SWIR bands to zero

RrsONTOSL8 = spect_sampL8(RrsONTOScorr(wlrange),wavelengthSVC(wlrange).*1E-3);

RrsONTOSL8corr = RrsONTOSL8;
RrsONTOSL8corr(5:7)=0;

L8bands = [0.4430,0.4826,0.5613,0.6546,0.8646,1.6090,2.2010];

figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,RrsONTOS150916,'r')
hold on
plot(wavelengthSVC,RrsONTOScorr,'--r')
% plot(wavelengthSVC(wlzero),RrsONTOS(wlzero),'.g')
plot(L8bands.*1E3,RrsONTOSL8,'.-b')
plot(L8bands.*1E3,RrsONTOSL8corr,'.-k')
plot(L8bands.*1E3,RrsONTOSL8corr*pi,'--k')
legend('RrsONTOS150916','RrsONTOScorr','RrsONTOSL8','RrsONTOSL8corr','RrsONTOSL8corr*pi')
title('R_{rs} -- 09/16/15 ','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([400 2500 0 0.03])
grid on
% 
% NSRef = [L8bands' RrsONTNSL8corr'];
% save([pathname,pathdate,'RrsONTNSL8.txt'],'NSRef','-ascii')

% Ref = [L8bands', RrsONTOSL8corr'];
% save('ONTOSL8_Rrs_150916corr.txt','Ref','-ascii')

%% Find best match in the HL LUT with 120 wl for ONTOS
% Run the LUT part of retrievalL8_150916.m first...
Rrs_SITE_test = RrsONTOS150916;

Rrs_SITE_test_HL = interp1(wavelengthSVC,Rrs_SITE_test,wavelength*1000);
Rrs_SITE_test_HL = Rrs_SITE_test_HL-Rrs_SITE_test_HL(end);

rule1 = strcmp(c{1}(:),'input150916ONTOS');
% rule1 = strcmp(c{1}(:),'input150916LONGS');
% rule1 = ~isnan(ones(size(Rrs,2),1));
c1_test = c{1}(rule1);
c2_test = c{2}(rule1);
c3_test = c{3}(rule1);
c4_test = c{4}(rule1);
c5_test = c{5}(rule1);


wl_nm = wavelength*1000;

% cond1 = wl_nm>500;
cond1 = wl_nm>0;

Rrs_test = Rrs(cond1,rule1);

figure
plot(wl_nm,Rrs_SITE_test_HL)
xlim([400 1000])
hold on

[Y,I2] = min(sqrt(mean((Rrs_test'-ones(size(Rrs_test,2),1)*Rrs_SITE_test_HL(cond1)').^2,2)));

plot(wl_nm(cond1),Rrs_test(:,I2),'g')
legend('Field','LUT')
str = sprintf('%s %f %f %f %s',char(c1_test(I2)),c2_test(I2),c3_test(I2),c4_test(I2),char(c5_test(I2)));
title(str)
grid on

%% Find best match in the HL LUT with 120 wl for LONGS
% Run the LUT part of retrievalL8_150916.m first...
Rrs_SITE_test = RrsLONGS150916;

Rrs_SITE_test_HL = interp1(wavelengthSVC,Rrs_SITE_test,wavelength*1000);
Rrs_SITE_test_HL = Rrs_SITE_test_HL-Rrs_SITE_test_HL(end);

rule1 = strcmp(c{1}(:),'input150916LONGS');
% rule1 = strcmp(c{1}(:),'input150916LONGS');
% rule1 = ~isnan(ones(size(Rrs,2),1));
c1_test = c{1}(rule1);
c2_test = c{2}(rule1);
c3_test = c{3}(rule1);
c4_test = c{4}(rule1);
c5_test = c{5}(rule1);


wl_nm = wavelength*1000;

cond1 = wl_nm>500;

Rrs_test = Rrs(cond1,rule1);

figure
plot(wl_nm,Rrs_SITE_test_HL)
xlim([400 1000])
hold on

[Y,I2] = min(sqrt(mean((Rrs_test'-ones(size(Rrs_test,2),1)*Rrs_SITE_test_HL(cond1)').^2,2)));

plot(wl_nm(cond1),Rrs_test(:,I2),'g')
legend('Field','LUT')
str = sprintf('%s %f %f %f %s',char(c1_test(I2)),c2_test(I2),c3_test(I2),c4_test(I2),char(c5_test(I2)));
title(str)
grid on

%%
figure
plot(wl_nm,Rrs)
hold on
plot(wl_nm(cond1),Rrs_test(:,I2),'g','linewidth',1.5)
%% Find best match in the HL LUT with 120 wl for CRANB
% Run the LUT part of retrievalL8_150916.m first...
Rrs_SITE_test = RrsCRANB150916;

Rrs_SITE_test_HL = interp1(wavelengthSVC,Rrs_SITE_test,wavelength*1000);
Rrs_SITE_test_HL = Rrs_SITE_test_HL-Rrs_SITE_test_HL(end);

rule1 = strcmp(c{1}(:),'input150916LONGS');
% rule1 = strcmp(c{1}(:),'input150916LONGS');
% rule1 = ~isnan(ones(size(Rrs,2),1));
c1_test = c{1}(rule1);
c2_test = c{2}(rule1);
c3_test = c{3}(rule1);
c4_test = c{4}(rule1);
c5_test = c{5}(rule1);


wl_nm = wavelength*1000;

cond1 = wl_nm>500;

Rrs_test = Rrs(cond1,rule1);

figure
plot(wl_nm,Rrs_SITE_test_HL)
xlim([400 1000])
hold on

[Y,I2] = min(sqrt(mean((Rrs_test'-ones(size(Rrs_test,2),1)*Rrs_SITE_test_HL(cond1)').^2,2)));

plot(wl_nm(cond1),Rrs_test(:,I2),'g')
legend('Field','LUT')
str = sprintf('%s %f %f %f %s',char(c1_test(I2)),c2_test(I2),c3_test(I2),c4_test(I2),char(c5_test(I2)));
title(str)
grid on

%%
figure
plot(wl_nm,Rrs)
hold on
plot(wl_nm(cond1),Rrs_test(:,I2),'g','linewidth',1.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAND 1 and SAND 2
% 150916_1706_R108_T109.sig
% 150916_1706_R108_T110.sig
% 150916_1710_R113_T114.sig
% 150916_1710_R113_T115.sig


filename = '150916_1706_R108_T109.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg1,Lt1,~] = extractSVC(filepath);

filename = '150916_1706_R108_T110.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = '150916_1710_R113_T114.sig';
filepath = [pathname,pathdate,filename];
[~,Lg2,Lt3,~] = extractSVC(filepath);

filename = '150916_1710_R113_T115.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt4,~] = extractSVC(filepath);


Rrs1 = (Lt1)./(pi.*Lg1./0.99);
Rrs2 = (Lt2)./(pi.*Lg1./0.99);
Rrs3 = (Lt3)./(pi.*Lg2./0.99);
Rrs4 = (Lt4)./(pi.*Lg2./0.99);

figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Rrs1,'r')
hold on
plot(wavelengthSVC,Rrs2,'--r')
plot(wavelengthSVC,Rrs3,'b')
plot(wavelengthSVC,Rrs4,'--b')
legend('Lt1 SAND 1','Lt2 SAND 1','Lt3 SAND 2',...
    'Lt4 SAND 2')

title('R_{rs} -- SAND 1 and SAND 2 ','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
% axis([0 3000 -0.01 0.01])
% axis([400 1000 0 0.015])
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %% Save for ENVI resampling to L8 response
% wavelengthSVC = CellRef{1,1}(:,1)*10^-3;
% cond = wavelengthSVC>=0.4 & wavelengthSVC<=1.0;
% wavelengthSVC = wavelengthSVC(cond);
% meanRef = mean(D,2)./100; % decimal
% ONTNSRef = meanRef(cond);
% 
% figure
% fs = 15;
% set(gcf,'color','white')
% set(gca,'fontsize',fs)
% plot(wavelengthSVC,ONTNSRef,'k')
% xlabel('wavelengthSVC [\mu m]','fontsize',fs)
% ylabel('reflectance','fontsize',fs)
% 
% 
% NSRef = [wavelengthSVC, ONTNSRef];
% save('ONTNSRef.txt','NSRef','-ascii')
