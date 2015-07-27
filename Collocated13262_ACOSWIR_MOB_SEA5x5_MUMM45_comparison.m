cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/
addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/gmregress/')
addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/retrieval/')
addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/130919')
%%
dir = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/';
tif13262 = 'LC80160302013262LGN00/SEADAS/Collocated13262_ACOSWIR_MOB_SEA5x5_MUMM45.tif';
% Image 2013262
filename = [dir tif13262];

im2013262 = imread(filename);

info = geotiffinfo(filename);
%1 lon_R_R_R_R_R
lon = double(im2013262(:,:,1));
%2 lat_R_R_R_R_R
lat = double(im2013262(:,:,2));
%3 rhoam_865_R_R_R_R_R
%4 rhow_443_ACO_R
%5 rhow_483_R_R_R_R_R
%6 rhow_561_R_R_R_R_R
%7 rhow_655_R_R_R_R_R
%8 rhow_865_R_R_R_R_R
%9 band_1_D_R_R_R_R
%10 band_2_D_R_R_R_R
%11 band_3_D_R_R_R_R
%12 band_4_D_R_R_R_R
%13 band_5_D_R_R_R_R
%14 band_6_D_R_R_R_R
%15 band_7_D_R_R_R_R
%16 lon_R_R_R_R
%17 lat_R_R_R_R
%18 Rrs_443_D_R_R_R
%19 chlor_a_D_R_R_R
%20 l2_flags_D_R_R_R
%21 Rrs_482_D_R_R_R
%22 Rrs_561_D_R_R_R
%23 Rrs_655_D_R_R_R
%24 longitude_R_R_R_R
%25 latitude_R_R_R_R
%26 lon_R_R_R
%27 lat_R_R_R
%28 chlor_a_MOB_D_R_R
CHL_MOB = double(im2013262(:,:,28)); % MOB-ELM
%29 lon_R_R
%30 lat_R_R
%31 Rrs_443_ACO_R_R
Rrs_443A = double(im2013262(:,:,31)); % Acolite
%32 Rrs_482_ACO_R_R
Rrs_483A = double(im2013262(:,:,32)); % Acolite
%33 Rrs_561_ACO_R_R
Rrs_561A = double(im2013262(:,:,33)); % Acolite
%34 Rrs_655_ACO_R_R
Rrs_655A = double(im2013262(:,:,34)); % Acolite
%35 Rrs_865_ACO_R_R
%36 Rrs_443_SEA5x5_R
Rrs_443S = double(im2013262(:,:,36))*2.0E-6+0.05; % SeaDAS, converted with scale_factor and add_offset from band attributes
%37 chlor_a_SEA5x5_R
CHL_SEA = double(im2013262(:,:,37)); % SeaDAS
%38 l2_flags_D_R
%39 Rrs_482_SEA5x5_R
Rrs_483S = double(im2013262(:,:,39))*2.0E-6+0.05; % SeaDAS, converted with scale_factor and add_offset from band attributes
%40 Rrs_561_SEA5x5_R
Rrs_561S = double(im2013262(:,:,40))*2.0E-6+0.05; % SeaDAS, converted with scale_factor and add_offset from band attributes
%41 Rrs_655_SEA5x5_R
Rrs_655S = double(im2013262(:,:,41))*2.0E-6+0.05; % SeaDAS, converted with scale_factor and add_offset from band attributes
%42 longitude_D_R
%43 latitude_D_R
%44 lon_R
%45 lat_R
%46 Rrs_443_MUMM45_R
Rrs_443M = double(im2013262(:,:,46))*2.0E-6+0.05; % MUMM (Ruddick);
%47 chlor_a_MUMM45_R
CHL_MUM = double(im2013262(:,:,47)); % MUMM
%48 l2_flags
%49 Rrs_482_MUMM45_R
Rrs_483M = double(im2013262(:,:,49))*2.0E-6+0.05; % MUMM (Ruddick);
%50 Rrs_561_MUMM45_R
Rrs_561M = double(im2013262(:,:,50))*2.0E-6+0.05; % MUMM (Ruddick);
%51 Rrs_655_MUMM45_R
Rrs_655M = double(im2013262(:,:,51))*2.0E-6+0.05; % MUMM (Ruddick);
%52 longitude
%53 latitude
%54 lat
%55 lon 
%56 Rrs_443_MOB
Rrs_443E = double(im2013262(:,:,56)); % MoB-ELM (E:ELM); for 2013262 image, not in Rrs 
%57 Rrs_482_MOB
Rrs_483E = double(im2013262(:,:,57)); % MoB-ELM (E:ELM); for 2013262 image, not in Rrs
%58 Rrs_561_MOB
Rrs_561E = double(im2013262(:,:,58)); % MoB-ELM (E:ELM); for 2013262 image, not in Rrs
%59 Rrs_655_MOB
Rrs_655E = double(im2013262(:,:,59)); % MoB-ELM (E:ELM); for 2013262 image, not in Rrs
%60 Rrs_865_MOB
%61 Rrs_1600_MOB
%62 Rrs_2100_MOB

tif13262 = 'LC80160302013262LGN00/ACOLITE/LC80160302013262LGN00_L2_SWIR_FranzAve.tif';
% Image 2013262
filename = [dir tif13262];

im2013262 = imread(filename);
%23 chlor_a_ACO_OC3
CHL_ACO = double(im2013262(:,:,23)); % Acolite
%% Chl ACOLITE
% R = log10(max(Rrs_443A,Rrs_483A)./Rrs_561A);
% CHL_ACO = 10.^(0.2412-2.0546.*R+1.1776.*R.^2-0.5538.*R.^3-0.4570.*R.^4);

%%
densityflag = 1;
regressiontype = 'RMA';
maxref = 0.030;
date = '2013262';

% Band 1: 443nm

% denscatplot(Rrs_443A,Rrs_443E,regressiontype,densityflag,'443',maxref,date,'Acolite','MoB-ELM')
% denscatplot(Rrs_443A,Rrs_443S,regressiontype,densityflag,'443',maxref,date,'Acolite','SeaDAS')
% denscatplot(Rrs_443A,Rrs_443M,regressiontype,densityflag,'443',maxref,date,'Acolite','MUMM')
% denscatplot(Rrs_443S,Rrs_443E,regressiontype,densityflag,'443',maxref,date,'SeaDAS','MoB-ELM')
% denscatplot(Rrs_443S,Rrs_443M,regressiontype,densityflag,'443',maxref,date,'SeaDAS','MUMM')
% denscatplot(Rrs_443M,Rrs_443E,regressiontype,densityflag,'443',maxref,date,'MUMM','MoB-ELM')

% Band 2: 483nm

% denscatplot(Rrs_483A,Rrs_483E,regressiontype,densityflag,'483',maxref,date,'Acolite','MoB-ELM')
% denscatplot(Rrs_483A,Rrs_483S,regressiontype,densityflag,'483',maxref,date,'Acolite','SeaDAS')
% denscatplot(Rrs_483A,Rrs_483M,regressiontype,densityflag,'483',maxref,date,'Acolite','MUMM')
% denscatplot(Rrs_483S,Rrs_483E,regressiontype,densityflag,'483',maxref,date,'SeaDAS','MoB-ELM')
% denscatplot(Rrs_483S,Rrs_483M,regressiontype,densityflag,'483',maxref,date,'SeaDAS','MUMM')
% denscatplot(Rrs_483M,Rrs_483E,regressiontype,densityflag,'483',maxref,date,'MUMM','MoB-ELM')

% Band 3: 561nm

% denscatplot(Rrs_561A,Rrs_561E,regressiontype,densityflag,'561',maxref,date,'Acolite','MoB-ELM')
% denscatplot(Rrs_561A,Rrs_561S,regressiontype,densityflag,'561',maxref,date,'Acolite','SeaDAS')
% denscatplot(Rrs_561A,Rrs_561M,regressiontype,densityflag,'561',maxref,date,'Acolite','MUMM')
% denscatplot(Rrs_561S,Rrs_561E,regressiontype,densityflag,'561',maxref,date,'SeaDAS','MoB-ELM')
% denscatplot(Rrs_561S,Rrs_561M,regressiontype,densityflag,'561',maxref,date,'SeaDAS','MUMM')
% denscatplot(Rrs_561M,Rrs_561E,regressiontype,densityflag,'561',maxref,date,'MUMM','MoB-ELM')

% Band 4: 655nm

% denscatplot(Rrs_655A,Rrs_655E,regressiontype,densityflag,'655',maxref,date,'Acolite','MoB-ELM')
% denscatplot(Rrs_655A,Rrs_655S,regressiontype,densityflag,'655',maxref,date,'Acolite','SeaDAS')
% denscatplot(Rrs_655A,Rrs_655M,regressiontype,densityflag,'655',maxref,date,'Acolite','MUMM')
% denscatplot(Rrs_655S,Rrs_655E,regressiontype,densityflag,'655',maxref,date,'SeaDAS','MoB-ELM')
% denscatplot(Rrs_655S,Rrs_655M,regressiontype,densityflag,'655',maxref,date,'SeaDAS','MUMM')
% denscatplot(Rrs_655M,Rrs_655E,regressiontype,densityflag,'655',maxref,date,'MUMM','MoB-ELM')

%% Chl
densityflag = 1;
regressiontype = 'RMA';
maxref = 150;
date = '2013262';

% sum(sum(CHL_ACO>150))
CHL_ACO(CHL_ACO>150)=NaN;
CHL_MOB(CHL_MOB>150)=NaN;
CHL_SEA(CHL_SEA>150)=NaN;
CHL_MUM(CHL_MUM>150)=NaN;
% sum(sum(CHL_ACO>150))

% denscatplot(CHL_ACO,CHL_MOB,regressiontype,densityflag,'C_a',maxref,date,'Acolite','MoB-ELM')
% 
% denscatplot(CHL_ACO,CHL_SEA,regressiontype,densityflag,'C_a',maxref,date,'Acolite','SeaDAS')
% 
% denscatplot(CHL_ACO,CHL_MUM,regressiontype,densityflag,'C_a',maxref,date,'Acolite','MUMM')
% 
% denscatplot(CHL_SEA,CHL_MOB,regressiontype,densityflag,'C_a',maxref,date,'SeaDAS','MoB-ELM')
% 
% denscatplot(CHL_SEA,CHL_MUM,regressiontype,densityflag,'C_a',maxref,date,'SeaDAS','MUMM')
% 
% denscatplot(CHL_MUM,CHL_MOB,regressiontype,densityflag,'C_a',maxref,date,'MUMM','MoB-ELM')

%if (rhos_469 != NaN and rhos_555 != NaN and rhos_645 != NaN)  then (       if (LAND)        then            (.091935692 + .61788 * atan(10*(rhos_645-.015)) )       else           (0.29319407 + 0.45585 * atan(50*(rhos_645-.015)) ) ) else NaN

%% Rrs Comparison with ground-truth
matObj = matfile('Rrs130919AllSites.mat');
whos(matObj)

load Rrs130919AllSites.mat

%% ONTNS from the field
sitename = 'ONTNS';
wlrange = wavelengthSVC130919>=400 & wavelengthSVC130919<=2500;
% wlzero = wavelengthSVC130919==2219.0;
wlavg = wavelengthSVC130919>=2000 & wavelengthSVC130919<=2350;

zeroavg = mean(RrsONTNS130919(wlavg));

RrsONTNS130919corr = RrsONTNS130919-zeroavg;

RrsONTNS130919L8 = spect_sampL8(RrsONTNS130919corr(wlrange),wavelengthSVC130919(wlrange).*1E-3);

RrsONTNS130919L8corr = RrsONTNS130919L8;
RrsONTNS130919L8corr(5:7)=0;

L8bands = [0.4430,0.4826,0.5613,0.6546,0.8646,1.6090,2.2010];

figure('name',sitename)
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC130919,RrsONTNS130919,'r')
hold on
plot(wavelengthSVC130919,RrsONTNS130919corr,'--r')
% plot(wavelengthSVC130919(wlzero),RrsONTNS130919(wlzero),'.g')
plot(L8bands.*1E3,RrsONTNS130919L8,'.-b')
plot(L8bands.*1E3,RrsONTNS130919L8corr,'.-k')
legend('RrsONTNS130919','RrsONTNS130919corr','RrsONTNS130919L8','RrsONTNS130919L8corr')
title('R_{rs} -- 09/19/13 ','fontsize',fs)
xlabel('wavelengthSVC130919 [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.03])
grid on
%% ONTNS Rrs comparison

ONTNSlat =	43.272159;
ONTNSlon = -77.538274;

dist2=sum(bsxfun(@minus, cat(3,ONTNSlat,ONTNSlon), cat(3,lat,lon)).^2,3);
[I,J]=find(dist2==min(dist2(:)));

% figure
% plot(lon(:),lat(:),'.')
% hold on
% plot(ONTNSlon,ONTNSlat,'r*')
% plot(lon(I,J),lat(I,J),'g*')

ONTNS_ACO = [Rrs_443A(I,J),Rrs_483A(I,J),Rrs_561A(I,J),Rrs_655A(I,J)];
ONTNS_SEA = [Rrs_443S(I,J),Rrs_483S(I,J),Rrs_561S(I,J),Rrs_655S(I,J)];
ONTNS_MOB = [Rrs_443E(I,J),Rrs_483E(I,J),Rrs_561E(I,J),Rrs_655E(I,J)];
ONTNS_MUM = [Rrs_443M(I,J),Rrs_483M(I,J),Rrs_561M(I,J),Rrs_655M(I,J)];

figure('name',sitename)
fs = 20;
lw = 1.5;
set(gcf,'color','white')
plot(L8bands(1:4).*1E3,RrsONTNS130919L8corr(1:4),'--b','LineWidth',lw)
hold on
plot(L8bands(1:4).*1E3,ONTNS_ACO(1:4),'r','LineWidth',lw)
plot(L8bands(1:4).*1E3,ONTNS_SEA(1:4),'g','LineWidth',lw)
plot(L8bands(1:4).*1E3,ONTNS_MOB(1:4),'b','LineWidth',lw)
plot(L8bands(1:4).*1E3,ONTNS_MUM(1:4),'k','LineWidth',lw)
% legend('Field','Acolite','SeaDAS','MoB-ELM','MUMM')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('R_{rs} [1/sr]','fontsize',fs)
set(gca,'fontsize',fs)
xlim([400 700])
ylim([0 0.016])
% axis([400 1000 0 0.03])
grid on

str3 = sprintf('RrsComp%s.eps',sitename);
dirname = '/Users/javier/Desktop/Javier/PHD_RIT/ConferencesAndApplications/2015_SPIE_SanDiego/Images/';
% print(h,[dirname str3],'-depsc')
saveas(gcf,[dirname str3],'epsc')

%% ONTOS from the field
sitename = 'ONTOS';
wlrange = wavelengthSVC130919>=400 & wavelengthSVC130919<=2500;
% wlzero = wavelengthSVC130919==2219.0;
wlavg = wavelengthSVC130919>=2000 & wavelengthSVC130919<=2350;

zeroavg = mean(RrsONTOS130919(wlavg));

RrsONTOS130919corr = RrsONTOS130919-zeroavg;

RrsONTOS130919L8 = spect_sampL8(RrsONTOS130919corr(wlrange),wavelengthSVC130919(wlrange).*1E-3);

RrsONTOS130919L8corr = RrsONTOS130919L8;
RrsONTOS130919L8corr(5:7)=0;

L8bands = [0.4430,0.4826,0.5613,0.6546,0.8646,1.6090,2.2010];

figure('name',sitename)
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC130919,RrsONTOS130919,'r')
hold on
plot(wavelengthSVC130919,RrsONTOS130919corr,'--r')
% plot(wavelengthSVC130919(wlzero),RrsONTOS130919(wlzero),'.g')
plot(L8bands.*1E3,RrsONTOS130919L8,'.-b')
plot(L8bands.*1E3,RrsONTOS130919L8corr,'.-k')
legend('RrsONTOS130919','RrsONTOS130919corr','RrsONTOS130919L8','RrsONTOS130919L8corr')
title('R_{rs} -- 09/19/13 ','fontsize',fs)
xlabel('wavelengthSVC130919 [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.03])
grid on
%% ONTOS Rrs comparison

ONTOSlat =	43.308923;
ONTOSlon = -77.540085;

dist2=sum(bsxfun(@minus, cat(3,ONTOSlat,ONTOSlon), cat(3,lat,lon)).^2,3);
[I,J]=find(dist2==min(dist2(:)));

% figure('name',sitename)
% plot(lon(:),lat(:),'.')
% hold on
% plot(ONTOSlon,ONTOSlat,'r*')
% plot(lon(I,J),lat(I,J),'g*')

ONTOS_ACO = [Rrs_443A(I,J),Rrs_483A(I,J),Rrs_561A(I,J),Rrs_655A(I,J)];
ONTOS_SEA = [Rrs_443S(I,J),Rrs_483S(I,J),Rrs_561S(I,J),Rrs_655S(I,J)];
ONTOS_MOB = [Rrs_443E(I,J),Rrs_483E(I,J),Rrs_561E(I,J),Rrs_655E(I,J)];
ONTOS_MUM = [Rrs_443M(I,J),Rrs_483M(I,J),Rrs_561M(I,J),Rrs_655M(I,J)];

figure('name',sitename)
fs = 20;
lw = 1.5;
set(gcf,'color','white')
plot(L8bands(1:4).*1E3,RrsONTOS130919L8corr(1:4),'--b','LineWidth',lw)
hold on
plot(L8bands(1:4).*1E3,ONTOS_ACO(1:4),'r','LineWidth',lw)
plot(L8bands(1:4).*1E3,ONTOS_SEA(1:4),'g','LineWidth',lw)
plot(L8bands(1:4).*1E3,ONTOS_MOB(1:4),'b','LineWidth',lw)
plot(L8bands(1:4).*1E3,ONTOS_MUM(1:4),'k','LineWidth',lw)
% legend('Field','Acolite','SeaDAS','MoB-ELM','MUMM')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('R_{rs} [1/sr]','fontsize',fs)
set(gca,'fontsize',fs)
xlim([400 700])
ylim([0 0.016])
% axis([400 1000 0 0.03])
grid on

str3 = sprintf('RrsComp%s.eps',sitename);
dirname = '/Users/javier/Desktop/Javier/PHD_RIT/ConferencesAndApplications/2015_SPIE_SanDiego/Images/';
% print(h,[dirname str3],'-depsc')
saveas(gcf,[dirname str3],'epsc')

%% ONTEX from the field
sitename = 'ONTEX';
wlrange = wavelengthSVC130919>=400 & wavelengthSVC130919<=2500;
% wlzero = wavelengthSVC130919==2219.0;
wlavg = wavelengthSVC130919>=2000 & wavelengthSVC130919<=2350;

zeroavg = mean(RrsONTEX130919(wlavg));

RrsONTEX130919corr = RrsONTEX130919-zeroavg;

RrsONTEX130919L8 = spect_sampL8(RrsONTEX130919corr(wlrange),wavelengthSVC130919(wlrange).*1E-3);

RrsONTEX130919L8corr = RrsONTEX130919L8;
RrsONTEX130919L8corr(5:7)=0;

L8bands = [0.4430,0.4826,0.5613,0.6546,0.8646,1.6090,2.2010];

figure('name',sitename)
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC130919,RrsONTEX130919,'r')
hold on
plot(wavelengthSVC130919,RrsONTEX130919corr,'--r')
% plot(wavelengthSVC130919(wlzero),RrsONTEX130919(wlzero),'.g')
plot(L8bands.*1E3,RrsONTEX130919L8,'.-b')
plot(L8bands.*1E3,RrsONTEX130919L8corr,'.-k')
legend('RrsONTEX130919','RrsONTEX130919corr','RrsONTEX130919L8','RrsONTEX130919L8corr')
title('R_{rs} -- 09/19/13 ','fontsize',fs)
xlabel('wavelengthSVC130919 [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.03])
grid on
%% ONTEX Rrs comparison

ONTEXlat =	43.244892;
ONTEXlon = -77.536671;

dist2=sum(bsxfun(@minus, cat(3,ONTEXlat,ONTEXlon), cat(3,lat,lon)).^2,3);
[I,J]=find(dist2==min(dist2(:)));

figure('name',sitename)
plot(lon(:),lat(:),'.')
hold on
plot(ONTEXlon,ONTEXlat,'r*')
plot(lon(I,J),lat(I,J),'g*')

ONTEX_ACO = [Rrs_443A(I,J),Rrs_483A(I,J),Rrs_561A(I,J),Rrs_655A(I,J)];
ONTEX_SEA = [Rrs_443S(I,J),Rrs_483S(I,J),Rrs_561S(I,J),Rrs_655S(I,J)];
ONTEX_MOB = [Rrs_443E(I,J),Rrs_483E(I,J),Rrs_561E(I,J),Rrs_655E(I,J)];
ONTEX_MUM = [Rrs_443M(I,J),Rrs_483M(I,J),Rrs_561M(I,J),Rrs_655M(I,J)];

figure('name',sitename)
fs = 20;
lw = 1.5;
set(gcf,'color','white')
plot(L8bands(1:4).*1E3,RrsONTEX130919L8corr(1:4),'--b','LineWidth',lw)
hold on
plot(L8bands(1:4).*1E3,ONTEX_ACO(1:4),'r','LineWidth',lw)
plot(L8bands(1:4).*1E3,ONTEX_SEA(1:4),'g','LineWidth',lw)
plot(L8bands(1:4).*1E3,ONTEX_MOB(1:4),'b','LineWidth',lw)
plot(L8bands(1:4).*1E3,ONTEX_MUM(1:4),'k','LineWidth',lw)
legend('Field','Acolite','SeaDAS','MoB-ELM','MUMM')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('R_{rs} [1/sr]','fontsize',fs)
set(gca,'fontsize',fs)
xlim([400 700])
ylim([0 0.016])
% axis([400 1000 0 0.03])
grid on

str3 = sprintf('RrsComp%s.eps',sitename);
dirname = '/Users/javier/Desktop/Javier/PHD_RIT/ConferencesAndApplications/2015_SPIE_SanDiego/Images/';
% print(h,[dirname str3],'-depsc')
saveas(gcf,[dirname str3],'epsc')

%% RVRPI from the field
sitename = 'RVRPI';
wlrange = wavelengthSVC130919>=400 & wavelengthSVC130919<=2500;
% wlzero = wavelengthSVC130919==2219.0;
wlavg = wavelengthSVC130919>=2000 & wavelengthSVC130919<=2350;

zeroavg = mean(RrsRVRPI130919(wlavg));

RrsRVRPI130919corr = RrsRVRPI130919-zeroavg;

RrsRVRPI130919L8 = spect_sampL8(RrsRVRPI130919corr(wlrange),wavelengthSVC130919(wlrange).*1E-3);

RrsRVRPI130919L8corr = RrsRVRPI130919L8;
RrsRVRPI130919L8corr(5:7)=0;

L8bands = [0.4430,0.4826,0.5613,0.6546,0.8646,1.6090,2.2010];

figure('name',sitename)
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC130919,RrsRVRPI130919,'r')
hold on
plot(wavelengthSVC130919,RrsRVRPI130919corr,'--r')
% plot(wavelengthSVC130919(wlzero),RrsRVRPI130919(wlzero),'.g')
plot(L8bands.*1E3,RrsRVRPI130919L8,'.-b')
plot(L8bands.*1E3,RrsRVRPI130919L8corr,'.-k')
legend('RrsRVRPI130919','RrsRVRPI130919corr','RrsRVRPI130919L8','RrsRVRPI130919L8corr')
title('R_{rs} -- 09/19/13 ','fontsize',fs)
xlabel('wavelengthSVC130919 [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.03])
grid on
%% RVRPI Rrs comparison

RVRPIlat =	43.259925;
RVRPIlon = -77.601587;

dist2=sum(bsxfun(@minus, cat(3,RVRPIlat,RVRPIlon), cat(3,lat,lon)).^2,3);
[I,J]=find(dist2==min(dist2(:)));

figure('name',sitename)
plot(lon(:),lat(:),'.')
hold on
plot(RVRPIlon,RVRPIlat,'r*')
plot(lon(I,J),lat(I,J),'g*')

RVRPI_ACO = [Rrs_443A(I,J),Rrs_483A(I,J),Rrs_561A(I,J),Rrs_655A(I,J)];
RVRPI_SEA = [Rrs_443S(I,J),Rrs_483S(I,J),Rrs_561S(I,J),Rrs_655S(I,J)];
RVRPI_MOB = [Rrs_443E(I,J),Rrs_483E(I,J),Rrs_561E(I,J),Rrs_655E(I,J)];
RVRPI_MUM = [Rrs_443M(I,J),Rrs_483M(I,J),Rrs_561M(I,J),Rrs_655M(I,J)];

figure('name',sitename)
fs = 20;
lw = 1.5;
set(gcf,'color','white')
plot(L8bands(1:4).*1E3,RrsRVRPI130919L8corr(1:4),'--b','LineWidth',lw)
hold on
plot(L8bands(1:4).*1E3,RVRPI_ACO(1:4),'r','LineWidth',lw)
plot(L8bands(1:4).*1E3,RVRPI_SEA(1:4),'g','LineWidth',lw)
plot(L8bands(1:4).*1E3,RVRPI_MOB(1:4),'b','LineWidth',lw)
plot(L8bands(1:4).*1E3,RVRPI_MUM(1:4),'k','LineWidth',lw)
% legend('Field','Acolite','SeaDAS','MoB-ELM','MUMM')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('R_{rs} [1/sr]','fontsize',fs)
set(gca,'fontsize',fs)
xlim([400 700])
ylim([0 0.016])
% axis([400 1000 0 0.03])
grid on

str3 = sprintf('RrsComp%s.eps',sitename);
dirname = '/Users/javier/Desktop/Javier/PHD_RIT/ConferencesAndApplications/2015_SPIE_SanDiego/Images/';
% print(h,[dirname str3],'-depsc')
saveas(gcf,[dirname str3],'epsc')

%% RVRPL from the field
sitename = 'RVRPL';
wlrange = wavelengthSVC130919>=400 & wavelengthSVC130919<=2500;
% wlzero = wavelengthSVC130919==2219.0;
wlavg = wavelengthSVC130919>=2000 & wavelengthSVC130919<=2350;

zeroavg = mean(RrsRVRPL130919(wlavg));

RrsRVRPL130919corr = RrsRVRPL130919-zeroavg;

RrsRVRPL130919L8 = spect_sampL8(RrsRVRPL130919corr(wlrange),wavelengthSVC130919(wlrange).*1E-3);

RrsRVRPL130919L8corr = RrsRVRPL130919L8;
RrsRVRPL130919L8corr(5:7)=0;

L8bands = [0.4430,0.4826,0.5613,0.6546,0.8646,1.6090,2.2010];

figure('name',sitename)
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC130919,RrsRVRPL130919,'r')
hold on
plot(wavelengthSVC130919,RrsRVRPL130919corr,'--r')
% plot(wavelengthSVC130919(wlzero),RrsRVRPL130919(wlzero),'.g')
plot(L8bands.*1E3,RrsRVRPL130919L8,'.-b')
plot(L8bands.*1E3,RrsRVRPL130919L8corr,'.-k')
legend('RrsRVRPL130919','RrsRVRPL130919corr','RrsRVRPL130919L8','RrsRVRPL130919L8corr')
title('R_{rs} -- 09/19/13 ','fontsize',fs)
xlabel('wavelengthSVC130919 [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.03])
grid on

%% RVRPL Rrs comparison

RVRPLlat =	43.270990;
RVRPLlon = -77.592282;

dist2=sum(bsxfun(@minus, cat(3,RVRPLlat,RVRPLlon), cat(3,lat,lon)).^2,3);
[I,J]=find(dist2==min(dist2(:)));

figure('name',sitename)
plot(lon(:),lat(:),'.')
hold on
plot(RVRPLlon,RVRPLlat,'r*')
plot(lon(I,J),lat(I,J),'g*')

RVRPL_ACO = [Rrs_443A(I,J),Rrs_483A(I,J),Rrs_561A(I,J),Rrs_655A(I,J)];
RVRPL_SEA = [Rrs_443S(I,J),Rrs_483S(I,J),Rrs_561S(I,J),Rrs_655S(I,J)];
RVRPL_MOB = [Rrs_443E(I,J),Rrs_483E(I,J),Rrs_561E(I,J),Rrs_655E(I,J)];
RVRPL_MUM = [Rrs_443M(I,J),Rrs_483M(I,J),Rrs_561M(I,J),Rrs_655M(I,J)];

figure('name',sitename)
fs = 20;
lw = 1.5;
set(gcf,'color','white')
plot(L8bands(1:4).*1E3,RrsRVRPL130919L8corr(1:4),'--b','LineWidth',lw)
hold on
plot(L8bands(1:4).*1E3,RVRPL_ACO(1:4),'r','LineWidth',lw)
plot(L8bands(1:4).*1E3,RVRPL_SEA(1:4),'g','LineWidth',lw)
plot(L8bands(1:4).*1E3,RVRPL_MOB(1:4),'b','LineWidth',lw)
plot(L8bands(1:4).*1E3,RVRPL_MUM(1:4),'k','LineWidth',lw)
% legend('Field','Acolite','SeaDAS','MoB-ELM','MUMM')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('R_{rs} [1/sr]','fontsize',fs)
set(gca,'fontsize',fs)
xlim([400 700])
ylim([0 0.016])
% axis([400 1000 0 0.03])
grid on

str3 = sprintf('RrsComp%s.eps',sitename);
dirname = '/Users/javier/Desktop/Javier/PHD_RIT/ConferencesAndApplications/2015_SPIE_SanDiego/Images/';
% print(h,[dirname str3],'-depsc')
saveas(gcf,[dirname str3],'epsc')

%% LONGN Rrs comparison

sitename = 'LONGN';

LONGNlat =	43.290836;
LONGNlon = -77.690662;

dist2=sum(bsxfun(@minus, cat(3,LONGNlat,LONGNlon), cat(3,lat,lon)).^2,3);
[I,J]=find(dist2==min(dist2(:)));

figure('name',sitename)
plot(lon(:),lat(:),'.')
hold on
plot(LONGNlon,LONGNlat,'r*')
plot(lon(I,J),lat(I,J),'g*')

LONGN_ACO = [Rrs_443A(I,J),Rrs_483A(I,J),Rrs_561A(I,J),Rrs_655A(I,J)];
LONGN_SEA = [Rrs_443S(I,J),Rrs_483S(I,J),Rrs_561S(I,J),Rrs_655S(I,J)];
LONGN_MOB = [Rrs_443E(I,J),Rrs_483E(I,J),Rrs_561E(I,J),Rrs_655E(I,J)];
LONGN_MUM = [Rrs_443M(I,J),Rrs_483M(I,J),Rrs_561M(I,J),Rrs_655M(I,J)];

figure('name',sitename)
fs = 20;
lw = 1.5;
set(gcf,'color','white')
% plot(L8bands(1:4).*1E3,RrsLONGN130919L8corr(1:4),'--b','LineWidth',lw)
hold on
plot(L8bands(1:4).*1E3,LONGN_ACO(1:4),'r','LineWidth',lw)
plot(L8bands(1:4).*1E3,LONGN_SEA(1:4),'g','LineWidth',lw)
plot(L8bands(1:4).*1E3,LONGN_MOB(1:4),'b','LineWidth',lw)
plot(L8bands(1:4).*1E3,LONGN_MUM(1:4),'k','LineWidth',lw)
% legend('Acolite','SeaDAS','MoB-ELM','MUMM')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('R_{rs} [1/sr]','fontsize',fs)
set(gca,'fontsize',fs)
xlim([400 700])
ylim([0 0.016])
% axis([400 1000 0 0.03])
grid on

str3 = sprintf('RrsComp%s.eps',sitename);
dirname = '/Users/javier/Desktop/Javier/PHD_RIT/ConferencesAndApplications/2015_SPIE_SanDiego/Images/';
% print(h,[dirname str3],'-depsc')
saveas(gcf,[dirname str3],'epsc')

%% LONGS Rrs comparison

sitename = 'LONGS';

LONGSlat =	43.289182;
LONGSlon = -77.696458;

dist2=sum(bsxfun(@minus, cat(3,LONGSlat,LONGSlon), cat(3,lat,lon)).^2,3);
[I,J]=find(dist2==min(dist2(:)));

figure('name',sitename)
plot(lon(:),lat(:),'.')
hold on
plot(LONGSlon,LONGSlat,'r*')
plot(lon(I,J),lat(I,J),'g*')

LONGS_ACO = [Rrs_443A(I,J),Rrs_483A(I,J),Rrs_561A(I,J),Rrs_655A(I,J)];
LONGS_SEA = [Rrs_443S(I,J),Rrs_483S(I,J),Rrs_561S(I,J),Rrs_655S(I,J)];
LONGS_MOB = [Rrs_443E(I,J),Rrs_483E(I,J),Rrs_561E(I,J),Rrs_655E(I,J)];
LONGS_MUM = [Rrs_443M(I,J),Rrs_483M(I,J),Rrs_561M(I,J),Rrs_655M(I,J)];

figure('name',sitename)
fs = 20;
lw = 1.5;
set(gcf,'color','white')
% plot(L8bands(1:4).*1E3,RrsLONGS130919L8corr(1:4),'--b','LineWidth',lw)
hold on
plot(L8bands(1:4).*1E3,LONGS_ACO(1:4),'r','LineWidth',lw)
plot(L8bands(1:4).*1E3,LONGS_SEA(1:4),'g','LineWidth',lw)
plot(L8bands(1:4).*1E3,LONGS_MOB(1:4),'b','LineWidth',lw)
plot(L8bands(1:4).*1E3,LONGS_MUM(1:4),'k','LineWidth',lw)
% legend('Acolite','SeaDAS','MoB-ELM','MUMM')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('R_{rs} [1/sr]','fontsize',fs)
set(gca,'fontsize',fs)
xlim([400 700])
ylim([0 0.016])
% axis([400 1000 0 0.03])
grid on

str3 = sprintf('RrsComp%s.eps',sitename);
dirname = '/Users/javier/Desktop/Javier/PHD_RIT/ConferencesAndApplications/2015_SPIE_SanDiego/Images/';
% print(h,[dirname str3],'-depsc')
saveas(gcf,[dirname str3],'epsc')

%% CRANB Rrs comparison

sitename = 'CRANB';

CRANBlat =	43.299938;
CRANBlon = -77.692915;

dist2=sum(bsxfun(@minus, cat(3,CRANBlat,CRANBlon), cat(3,lat,lon)).^2,3);
[I,J]=find(dist2==min(dist2(:)));

figure('name',sitename)
plot(lon(:),lat(:),'.')
hold on
plot(CRANBlon,CRANBlat,'r*')
plot(lon(I,J),lat(I,J),'g*')

CRANB_ACO = [Rrs_443A(I,J),Rrs_483A(I,J),Rrs_561A(I,J),Rrs_655A(I,J)];
CRANB_SEA = [Rrs_443S(I,J),Rrs_483S(I,J),Rrs_561S(I,J),Rrs_655S(I,J)];
CRANB_MOB = [Rrs_443E(I,J),Rrs_483E(I,J),Rrs_561E(I,J),Rrs_655E(I,J)];
CRANB_MUM = [Rrs_443M(I,J),Rrs_483M(I,J),Rrs_561M(I,J),Rrs_655M(I,J)];

figure('name',sitename)
fs = 20;
lw = 1.5;
set(gcf,'color','white')
% plot(L8bands(1:4).*1E3,RrsCRANB130919L8corr(1:4),'--b','LineWidth',lw)
hold on
plot(L8bands(1:4).*1E3,CRANB_ACO(1:4),'r','LineWidth',lw)
plot(L8bands(1:4).*1E3,CRANB_SEA(1:4),'g','LineWidth',lw)
plot(L8bands(1:4).*1E3,CRANB_MOB(1:4),'b','LineWidth',lw)
plot(L8bands(1:4).*1E3,CRANB_MUM(1:4),'k','LineWidth',lw)
% legend('Acolite','SeaDAS','MoB-ELM','MUMM')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('R_{rs} [1/sr]','fontsize',fs)
set(gca,'fontsize',fs)
xlim([400 700])
ylim([0 0.016])
% axis([400 1000 0 0.03])
grid on

str3 = sprintf('RrsComp%s.eps',sitename);
dirname = '/Users/javier/Desktop/Javier/PHD_RIT/ConferencesAndApplications/2015_SPIE_SanDiego/Images/';
% print(h,[dirname str3],'-depsc')
saveas(gcf,[dirname str3],'epsc')

%% BRADIN Rrs comparison

sitename = 'BRADIN';

BRADINlat =	43.313675;
BRADINlon = -77.717531;

dist2=sum(bsxfun(@minus, cat(3,BRADINlat,BRADINlon), cat(3,lat,lon)).^2,3);
[I,J]=find(dist2==min(dist2(:)));

figure('name',sitename)
plot(lon(:),lat(:),'.')
hold on
plot(BRADINlon,BRADINlat,'r*')
plot(lon(I,J),lat(I,J),'g*')

BRADIN_ACO = [Rrs_443A(I,J),Rrs_483A(I,J),Rrs_561A(I,J),Rrs_655A(I,J)];
BRADIN_SEA = [Rrs_443S(I,J),Rrs_483S(I,J),Rrs_561S(I,J),Rrs_655S(I,J)];
BRADIN_MOB = [Rrs_443E(I,J),Rrs_483E(I,J),Rrs_561E(I,J),Rrs_655E(I,J)];
BRADIN_MUM = [Rrs_443M(I,J),Rrs_483M(I,J),Rrs_561M(I,J),Rrs_655M(I,J)];

figure('name',sitename)
fs = 20;
lw = 1.5;
set(gcf,'color','white')
% plot(L8bands(1:4).*1E3,RrsBRADIN130919L8corr(1:4),'--b','LineWidth',lw)
hold on
plot(L8bands(1:4).*1E3,BRADIN_ACO(1:4),'r','LineWidth',lw)
plot(L8bands(1:4).*1E3,BRADIN_SEA(1:4),'g','LineWidth',lw)
plot(L8bands(1:4).*1E3,BRADIN_MOB(1:4),'b','LineWidth',lw)
plot(L8bands(1:4).*1E3,BRADIN_MUM(1:4),'k','LineWidth',lw)
% legend('Acolite','SeaDAS','MoB-ELM','MUMM')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('R_{rs} [1/sr]','fontsize',fs)
set(gca,'fontsize',fs)
xlim([400 700])
ylim([0 0.016])
% axis([400 1000 0 0.03])
grid on

str3 = sprintf('RrsComp%s.eps',sitename);
dirname = '/Users/javier/Desktop/Javier/PHD_RIT/ConferencesAndApplications/2015_SPIE_SanDiego/Images/';
% print(h,[dirname str3],'-depsc')
saveas(gcf,[dirname str3],'epsc')

%% BRADONT Rrs comparison

sitename = 'BRADONT';

BRADONTlat =	43.313675;
BRADONTlon = -77.717531;

dist2=sum(bsxfun(@minus, cat(3,BRADONTlat,BRADONTlon), cat(3,lat,lon)).^2,3);
[I,J]=find(dist2==min(dist2(:)));

figure('name',sitename)
plot(lon(:),lat(:),'.')
hold on
plot(BRADONTlon,BRADONTlat,'r*')
plot(lon(I,J),lat(I,J),'g*')

BRADONT_ACO = [Rrs_443A(I,J),Rrs_483A(I,J),Rrs_561A(I,J),Rrs_655A(I,J)];
BRADONT_SEA = [Rrs_443S(I,J),Rrs_483S(I,J),Rrs_561S(I,J),Rrs_655S(I,J)];
BRADONT_MOB = [Rrs_443E(I,J),Rrs_483E(I,J),Rrs_561E(I,J),Rrs_655E(I,J)];
BRADONT_MUM = [Rrs_443M(I,J),Rrs_483M(I,J),Rrs_561M(I,J),Rrs_655M(I,J)];

figure('name',sitename)
fs = 20;
lw = 1.5;
set(gcf,'color','white')
% plot(L8bands(1:4).*1E3,RrsBRADONT130919L8corr(1:4),'--b','LineWidth',lw)
hold on
plot(L8bands(1:4).*1E3,BRADONT_ACO(1:4),'r','LineWidth',lw)
plot(L8bands(1:4).*1E3,BRADONT_SEA(1:4),'g','LineWidth',lw)
plot(L8bands(1:4).*1E3,BRADONT_MOB(1:4),'b','LineWidth',lw)
plot(L8bands(1:4).*1E3,BRADONT_MUM(1:4),'k','LineWidth',lw)
% legend('Acolite','SeaDAS','MoB-ELM','MUMM')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('R_{rs} [1/sr]','fontsize',fs)
set(gca,'fontsize',fs)
xlim([400 700])
ylim([0 0.016])
% axis([400 1000 0 0.03])
grid on

str3 = sprintf('RrsComp%s.eps',sitename);
dirname = '/Users/javier/Desktop/Javier/PHD_RIT/ConferencesAndApplications/2015_SPIE_SanDiego/Images/';
% print(h,[dirname str3],'-depsc')
saveas(gcf,[dirname str3],'epsc')