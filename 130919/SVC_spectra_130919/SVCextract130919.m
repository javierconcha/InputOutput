cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/
addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/retrieval/')
%%
date = '130919';
pathname = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/';
pathdate = [date,'/SVC_spectra_130919/'];
%% Lsky from MODTRAN obtained from ../InputOutput/open_tape7scn.m
% method = 'rloess';
% SOL_SCATsmooth = smooth(SOL_SCAT,method);
% Lsky = interp1(WAVLEN_MCRN*1E3,SOL_SCATsmooth*1E7,wavelength);

Lsky = interp1(WAVLEN_MCRN*1E3,SOL_SCAT*1E7,wavelength); % sky was not measured
%% ONTNS 1
% L8_2013_09_19_R152_T153.sig
% L8_2013_09_19_R152_T154.sig
% L8_2013_09_19_R152_T155.sig
% L8_2013_09_19_R152_T156.sig
% L8_2013_09_19_R152_T157.sig
% L8_2013_09_19_R152_T158.sig

filename = 'L8_2013_09_19_R152_T153.sig';
filepath = [pathname,pathdate,filename];
[wavelength,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R152_T154.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R152_T155.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R152_T156.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt4,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R152_T157.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt5,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R152_T158.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt6,~] = extractSVC(filepath);

% Lt1smooth = smooth(Lt1,method);
% Lt2smooth = smooth(Lt2,method);
% Lt3smooth = smooth(Lt3,method);
% Lt4smooth = smooth(Lt4,method);
% Lt5smooth = smooth(Lt5,method);
% Lt6smooth = smooth(Lt6,method);

% Lgsmooth = smooth(Lg,method);

% Rrs1 = (Lt1smooth-0.028.*Lsky)./(pi.*Lgsmooth./0.99);
% Rrs2 = (Lt2smooth-0.028.*Lsky)./(pi.*Lgsmooth./0.99);
% Rrs3 = (Lt3smooth-0.028.*Lsky)./(pi.*Lgsmooth./0.99);
% Rrs4 = (Lt4smooth-0.028.*Lsky)./(pi.*Lgsmooth./0.99);
% Rrs5 = (Lt5smooth-0.028.*Lsky)./(pi.*Lgsmooth./0.99);
% Rrs6 = (Lt6smooth-0.028.*Lsky)./(pi.*Lgsmooth./0.99);

Rrs1a = (Lt1-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs2a = (Lt2-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs3a = (Lt3-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs4a = (Lt4-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs5a = (Lt5-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs6a = (Lt6-0.028.*Lsky)./(pi.*Lg./0.99);


figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,Rrs1a,'r')
hold on
plot(wavelength,Rrs2a,'b')
plot(wavelength,Rrs3a,'g')
plot(wavelength,Rrs4a,'--m')
plot(wavelength,Rrs5a,'--c')
plot(wavelength,Rrs6a,'--k')
legend('Lt1','Lt2','Lt3',...
    'Lt4','Lt5','Lt6')
% plot(wavelength,r./100)

title('R_{rs} -- ONTNS 1','fontsize',fs)
% title(method,'fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([0 3000 -0.01 0.01])
% axis([400 1000 0 0.015])
grid on


%% Radiance Lsky
figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,Lsky,'r')
% hold on
% plot(wavelength,Lskyb,'b')

% legend('Lskya','Lskyb')

title('Lsky -- ONTNS ','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.015])
grid on

%% ONTNS 2
% L8_2013_09_19_R159_T160.sig
% L8_2013_09_19_R159_T161.sig
% L8_2013_09_19_R159_T162.sig
% L8_2013_09_19_R159_T163.sig
% L8_2013_09_19_R159_T164.sig
% L8_2013_09_19_R159_T165.sig
% L8_2013_09_19_R159_T166.sig
% L8_2013_09_19_R159_T167.sig
% L8_2013_09_19_R159_T168.sig

filename = 'L8_2013_09_19_R159_T160.sig';
filepath = [pathname,pathdate,filename];
[wavelength,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R159_T161.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R159_T162.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R159_T163.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt4,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R159_T164.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt5,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R159_T165.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt6,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R159_T166.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt7,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R159_T167.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt8,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R159_T168.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt9,~] = extractSVC(filepath);

Rrs1b = (Lt1-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs2b = (Lt2-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs3b = (Lt3-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs4b = (Lt4-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs5b = (Lt5-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs6b = (Lt6-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs7b = (Lt7-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs8b = (Lt8-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs9b = (Lt9-0.028.*Lsky)./(pi.*Lg./0.99);

RrsONTNS = Rrs4b;

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,Rrs1b,'r')
hold on
plot(wavelength,Rrs2b,'b')
plot(wavelength,Rrs3b,'g')
plot(wavelength,Rrs4b,'--m')
plot(wavelength,Rrs5b,'--c')
plot(wavelength,Rrs6b,'--k')
plot(wavelength,Rrs7b,'-.m')
plot(wavelength,Rrs8b,'-.c')
plot(wavelength,Rrs9b,'-.k')
legend('Lt1','Lt2','Lt3',...
    'Lt4','Lt5','Lt6',...
    'Lt7','Lt8','Lt9')
% plot(wavelength,r./100)

title('R_{rs} -- ONTNS 2 ','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([0 3000 -0.01 0.01])
% axis([400 1000 0 0.015])
grid on

%% ONTOS

% L8_2013_09_19_R169_T170.sig
% L8_2013_09_19_R169_T171.sig
% L8_2013_09_19_R169_T172.sig
% L8_2013_09_19_R169_T173.sig
% L8_2013_09_19_R169_T174.sig

filename = 'L8_2013_09_19_R169_T170.sig';
filepath = [pathname,pathdate,filename];
[wavelength,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R169_T171.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,r] = extractSVC(filepath);

filename = 'L8_2013_09_19_R169_T172.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R169_T173.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt4,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R169_T174.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt5,~] = extractSVC(filepath);

Rrs1 = (Lt1-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs2 = (Lt2-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs3 = (Lt3-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs4 = (Lt4-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs5 = (Lt5-0.028.*Lsky)./(pi.*Lg./0.99);

RrsONTOS = Rrs1;

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,Rrs1,'r')
hold on
plot(wavelength,Rrs2,'b')
plot(wavelength,Rrs3,'g')
plot(wavelength,Rrs4,'--m')
plot(wavelength,Rrs5,'--c')


title('R_{rs} -- ONTOS ','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
% axis([0 3000 -0.01 0.01])
% axis([400 1000 0 0.015])
grid on

hold on
Rrs1old = (Lt2)./(pi.*Lg./0.99);
plot(wavelength,Rrs1old,'k')
legend('Lt1','Lt2','Lt3',...
    'Lt4','Lt5','w/o correction')

% %% ONTOS used for model-based ELM
% figure
% fs = 15;
% set(gcf,'color','white')
% 
% plot(wavelength,Rrs2,'b')
% hold on
% plot(wavelength,r/100/pi,'k')
% title('R_{rs} -- ONTOS ','fontsize',fs)
% xlabel('wavelength [nm]','fontsize',fs)
% ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
% set(gca,'fontsize',fs)
% % axis([400 2500 -0.001 0.015])
% grid on

%% RVRPLM

% L8_2013_09_19_R175_T176.sig
% L8_2013_09_19_R175_T177.sig
% L8_2013_09_19_R175_T178.sig
% L8_2013_09_19_R175_T179.sig
% L8_2013_09_19_R175_T180.sig
% L8_2013_09_19_R175_T181.sig

filename = 'L8_2013_09_19_R175_T176.sig';
filepath = [pathname,pathdate,filename];
[wavelength,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R175_T177.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R175_T178.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R175_T179.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt4,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R175_T180.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt5,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R175_T181.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt6,~] = extractSVC(filepath);

Rrs1 = (Lt1-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs2 = (Lt2-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs3 = (Lt3-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs4 = (Lt4-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs5 = (Lt5-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs6 = (Lt6-0.028.*Lsky)./(pi.*Lg./0.99);

RrsRVRPLM = Rrs1;

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,Rrs1,'r')
hold on
plot(wavelength,Rrs2,'b')
plot(wavelength,Rrs3,'g')
plot(wavelength,Rrs4,'--m')
plot(wavelength,Rrs5,'--c')
plot(wavelength,Rrs6,'--k')
legend('Lt1','Lt2','Lt3',...
    'Lt4','Lt5','Lt6')

title('R_{rs} -- RVRPLM ','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([0 3000 -0.01 0.01])
% axis([400 1000 0 0.015])
grid on

%% RVRPIER

% L8_2013_09_19_R182_T183.sig
% L8_2013_09_19_R182_T184.sig
% L8_2013_09_19_R182_T185.sig
% L8_2013_09_19_R182_T186.sig
% L8_2013_09_19_R182_T187.sig
% L8_2013_09_19_R182_T188.sig
% L8_2013_09_19_R182_T189.sig
% L8_2013_09_19_R182_T190.sig

filename = 'L8_2013_09_19_R182_T183.sig';
filepath = [pathname,pathdate,filename];
[wavelength,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R182_T184.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R182_T185.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R182_T186.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt4,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R182_T187.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt5,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R182_T188.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt6,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R182_T189.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt7,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R182_T190.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt8,~] = extractSVC(filepath);

Rrs1 = (Lt1-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs2 = (Lt2-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs3 = (Lt3-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs4 = (Lt4-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs5 = (Lt5-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs6 = (Lt6-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs7 = (Lt7-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs8 = (Lt8-0.028.*Lsky)./(pi.*Lg./0.99);

RrsRVRPIER = Rrs1;

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,Rrs1,'r')
hold on
plot(wavelength,Rrs2,'b')
plot(wavelength,Rrs3,'g')
plot(wavelength,Rrs4,'--m')
plot(wavelength,Rrs5,'--c')
plot(wavelength,Rrs6,'--k')
plot(wavelength,Rrs7,'-.m')
plot(wavelength,Rrs8,'-.c')
legend('Lt1','Lt2','Lt3',...
    'Lt4','Lt5','Lt6',...
    'Lt7','Lt8')

title('R_{rs} -- RVRPIER ','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
% axis([0 3000 -0.01 0.01])
% axis([400 1000 0 0.015])
grid on

%% ONTEX

% L8_2013_09_19_R191_T192.sig
% L8_2013_09_19_R191_T193.sig
% L8_2013_09_19_R191_T194.sig
% L8_2013_09_19_R191_T195.sig
% L8_2013_09_19_R191_T196.sig
% L8_2013_09_19_R191_T197.sig
% L8_2013_09_19_R191_T198.sig

filename = 'L8_2013_09_19_R191_T192.sig';
filepath = [pathname,pathdate,filename];
[wavelength,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R191_T193.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R191_T194.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R191_T195.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt4,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R191_T196.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt5,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R191_T197.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt6,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R191_T198.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt7,~] = extractSVC(filepath);

Rrs1 = (Lt1-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs2 = (Lt2-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs3 = (Lt3-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs4 = (Lt4-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs5 = (Lt5-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs6 = (Lt6-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs7 = (Lt7-0.028.*Lsky)./(pi.*Lg./0.99);

RrsONTEX = Rrs1;

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,Rrs1,'r')
hold on
plot(wavelength,Rrs2,'b')
plot(wavelength,Rrs3,'g')
plot(wavelength,Rrs4,'--m')
plot(wavelength,Rrs5,'--c')
plot(wavelength,Rrs6,'--k')
plot(wavelength,Rrs7,'-.m')
legend('Lt1','Lt2','Lt3',...
    'Lt4','Lt5','Lt6',...
    'Lt7')
% plot(wavelength,r./100)

title('R_{rs} -- ONTEX ','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([0 3000 -0.01 0.01])
% axis([400 1000 0 0.015])
grid on

%% IBAYN

% L8_2013_09_19_R199_T200.sig
% L8_2013_09_19_R199_T201.sig
% L8_2013_09_19_R199_T202.sig
% L8_2013_09_19_R199_T203.sig
% L8_2013_09_19_R199_T204.sig
% L8_2013_09_19_R199_T205.sig
% L8_2013_09_19_R199_T206.sig
% L8_2013_09_19_R199_T207.sig
% L8_2013_09_19_R199_T208.sig
% L8_2013_09_19_R199_T209.sig

filename = 'L8_2013_09_19_R199_T200.sig';
filepath = [pathname,pathdate,filename];
[wavelength,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R199_T201.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R199_T202.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R199_T203.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt4,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R199_T204.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt5,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R199_T205.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt6,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R199_T206.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt7,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R199_T207.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt8,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R199_T208.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt9,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R199_T209.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt10,~] = extractSVC(filepath);

Rrs1 = (Lt1-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs2 = (Lt2-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs3 = (Lt3-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs4 = (Lt4-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs5 = (Lt5-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs6 = (Lt6-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs7 = (Lt7-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs8 = (Lt8-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs9 = (Lt9-0.028.*Lsky)./(pi.*Lg./0.99);
Rrs10 = (Lt10-0.028.*Lsky)./(pi.*Lg./0.99);

RrsIBAYN = Rrs3;

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,Rrs1,'r')
hold on
plot(wavelength,Rrs2,'b')
plot(wavelength,Rrs3,'g')
plot(wavelength,Rrs4,'--m')
plot(wavelength,Rrs5,'--c')
plot(wavelength,Rrs6,'--k')
plot(wavelength,Rrs7,'-.m')
plot(wavelength,Rrs8,'-.c')
plot(wavelength,Rrs9,'-.k')
plot(wavelength,Rrs10,'k')
legend('Lt1','Lt2','Lt3',...
    'Lt4','Lt5','Lt6',...
    'Lt7','Lt8','Lt9','Lt10')

title('R_{rs} -- IBAYN ','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
% axis([0 3000 -0.01 0.01])
% axis([400 1000 0 0.015])
grid on

%% Rrs all
figure(21)
fs = 15;
set(gcf,'color','white')
plot(wavelength,RrsONTNS,'r')
hold on
plot(wavelength,RrsONTOS,'b')
plot(wavelength,RrsRVRPLM,'g')
plot(wavelength,RrsRVRPIER,'m')
plot(wavelength,RrsIBAYN,'k')
plot(wavelength,RrsONTEX,'--b')

legend('ONTNS','ONTOS','RVRPLM',...
    'RVRPIER','IBAYN','ONTEX')
title('R_{rs} -- 09/19/13 ','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.03])
grid on
%% Spectrally sampled and save in text file
wlrange = wavelength>=400 & wavelength<=2500;
% wlzero = wavelength==2219.0;
wlavg = wavelength>=2000 & wavelength<=2350;

zeroavg = mean(RrsONTNS(wlavg));

RrsONTNScorr = RrsONTNS-zeroavg;

RrsONTNSL8 = spect_sampL8(RrsONTNScorr(wlrange),wavelength(wlrange).*1E-3);

RrsONTNSL8corr = RrsONTNSL8;
RrsONTNSL8corr(5:7)=0;

L8bands = [0.4430,0.4826,0.5613,0.6546,0.8646,1.6090,2.2010];

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,RrsONTNS,'r')
hold on
plot(wavelength,RrsONTNScorr,'--r')
% plot(wavelength(wlzero),RrsONTNS(wlzero),'.g')
plot(L8bands.*1E3,RrsONTNSL8,'.-b')
plot(L8bands.*1E3,RrsONTNSL8corr,'.-k')
plot(L8bands.*1E3,RrsONTNSL8corr*pi,'--k')
legend('RrsONTNS','RrsONTNScorr','RrsONTNSL8','RrsONTNSL8corr','RrsONTNSL8corr*pi')
title('R_{rs} -- 09/19/13 ','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.03])
grid on
% 
% NSRef = [L8bands', RrsONTNSL8corr'];
% save([pathname,pathdate,'RrsONTNSL8.txt'],'NSRef','-ascii')

%% Sand spectra -- dried
% L8_2013_09_19_R210_T211.sig
% L8_2013_09_19_R210_T212.sig
% L8_2013_09_19_R210_T213.sig
% L8_2013_09_19_R210_T214.sig
% L8_2013_09_19_R210_T215.sig
% L8_2013_09_19_R210_T216.sig
% L8_2013_09_19_R210_T217.sig
% L8_2013_09_19_R210_T218.sig

filename = 'L8_2013_09_19_R210_T211.sig';
filepath = [pathname,pathdate,filename];
[wavelength,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R210_T212.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R210_T213.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R210_T214.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt4,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R210_T215.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt5,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R210_T216.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt6,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R210_T217.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt7,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R210_T218.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt8,~] = extractSVC(filepath);

Rrs1 = (Lt1)./(pi.*Lg./0.99);
Rrs2 = (Lt2)./(pi.*Lg./0.99);
Rrs3 = (Lt3)./(pi.*Lg./0.99);
Rrs4 = (Lt4)./(pi.*Lg./0.99);
Rrs5 = (Lt5)./(pi.*Lg./0.99);
Rrs6 = (Lt6)./(pi.*Lg./0.99);
Rrs7 = (Lt7)./(pi.*Lg./0.99);
Rrs8 = (Lt8)./(pi.*Lg./0.99);

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,Rrs1,'r')
hold on
plot(wavelength,Rrs2,'b')
plot(wavelength,Rrs3,'g')
plot(wavelength,Rrs4,'--m')
plot(wavelength,Rrs5,'--c')
plot(wavelength,Rrs6,'--k')
plot(wavelength,Rrs7,'-.m')
plot(wavelength,Rrs8,'-.c')
legend('Lt1','Lt2','Lt3',...
    'Lt4','Lt5','Lt6',...
    'Lt7','Lt8')

title('R_{rs} -- Dried Sand ','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
% axis([0 3000 -0.01 0.01])
% axis([400 1000 0 0.015])
grid on

%% Sand spectra -- wet

% L8_2013_09_19_R219_T220.sig
% L8_2013_09_19_R219_T221.sig
% L8_2013_09_19_R219_T222.sig
% L8_2013_09_19_R219_T223.sig
% L8_2013_09_19_R219_T224.sig

filename = 'L8_2013_09_19_R219_T220.sig';
filepath = [pathname,pathdate,filename];
[wavelength,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R219_T221.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R219_T222.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R219_T223.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt4,~] = extractSVC(filepath);

filename = 'L8_2013_09_19_R219_T224.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt5,~] = extractSVC(filepath);

Rrs1 = (Lt1)./(pi.*Lg./0.99);
Rrs2 = (Lt2)./(pi.*Lg./0.99);
Rrs3 = (Lt3)./(pi.*Lg./0.99);
Rrs4 = (Lt4)./(pi.*Lg./0.99);
Rrs5 = (Lt5)./(pi.*Lg./0.99);

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,Rrs1,'r')
hold on
plot(wavelength,Rrs2,'b')
plot(wavelength,Rrs3,'g')
plot(wavelength,Rrs4,'--m')
plot(wavelength,Rrs5,'--c')
legend('Lt1','Lt2','Lt3',...
    'Lt4','Lt5')

title('R_{rs} -- Wet Sand ','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
% axis([0 3000 -0.01 0.01])
% axis([400 1000 0 0.015])
grid on

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
% xlabel('wavelength [\mu m]','fontsize',fs)
% ylabel('reflectance','fontsize',fs)
% 
% 
% NSRef = [wavelengthSVC, ONTNSRef];
% save('ONTNSRef.txt','NSRef','-ascii')
