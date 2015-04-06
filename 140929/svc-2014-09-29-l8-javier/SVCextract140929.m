cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929
addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/retrieval/')
addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/')
%%
date = '140929';
pathname = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/';
pathdate = [date,'/svc-2014-09-29-l8-javier/'];

rho = 0.028;

%% ONTOS 1
% L8_2014_09_29_R185_T186.sig water
% L8_2014_09_29_R185_T187.sig water
% L8_2014_09_29_R185_T188.sig water
% L8_2014_09_29_R185_T189.sig sky
% L8_2014_09_29_R185_T190.sig sky


filename = 'L8_2014_09_29_R185_T186.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R185_T187.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R185_T188.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R185_T189.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskya,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R185_T190.sig';
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

title('R_{rs} -- ONTOS 1','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([400 1000 0 0.03])
grid on

%% Radiance Lsky
figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Lskya,'r')
hold on
plot(wavelengthSVC,Lskyb,'b')

legend('Lskya','Lskyb')

title('Lsky -- ONTOS 1','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.015])
grid on

%% ONTOS 2
% L8_2014_09_29_R191_T192.sig water
% L8_2014_09_29_R191_T193.sig water
% L8_2014_09_29_R191_T194.sig water
% L8_2014_09_29_R191_T195.sig sky
% L8_2014_09_29_R191_T196.sig sky


filename = 'L8_2014_09_29_R191_T192.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R191_T193.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R191_T194.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R191_T195.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskya,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R191_T196.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskyb,~] = extractSVC(filepath);

Rrs1a = (Lt1-rho.*Lskya)./(pi.*Lg./0.99);
Rrs2a = (Lt2-rho.*Lskya)./(pi.*Lg./0.99);
Rrs3a = (Lt3-rho.*Lskya)./(pi.*Lg./0.99);

Rrs1b = (Lt1-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs2b = (Lt2-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs3b = (Lt3-rho.*Lskyb)./(pi.*Lg./0.99);

RrsONTOS140929 = Rrs1a;

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

title('R_{rs} -- ONTOS 2','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([400 1000 0 0.03])
grid on

%% Radiance Lsky
figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Lskya,'r')
hold on
plot(wavelengthSVC,Lskyb,'b')

legend('Lskya','Lskyb')

title('Lsky -- ONTOS 2','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.015])
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LONGS 1
% L8_2014_09_29_R157_T158.sig water
% L8_2014_09_29_R157_T159.sig water
% L8_2014_09_29_R157_T160.sig water
% L8_2014_09_29_R157_T161.sig sky
% L8_2014_09_29_R157_T162.sig sky
% L8_2014_09_29_R157_T163.sig sky


filename = 'L8_2014_09_29_R157_T158.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R157_T159.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R157_T160.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R157_T161.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskya,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R157_T162.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskyb,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R157_T163.sig';
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

title('R_{rs} -- LONGS 1','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([400 1000 0 0.03])
grid on
%%
% Ref = [wavelengthSVC*1E-3, RrsLONGS140929];
% save('LONGSRef_140919.txt','Ref','-ascii')
%% Radiance Lsky
figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Lskya,'r')
hold on
plot(wavelengthSVC,Lskyb,'b')
plot(wavelengthSVC,Lskyc,'g')
legend('Lskya','Lskyb','Lskyc')

title('Lsky -- LONGS 1','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.015])
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LONGS 2
% L8_2014_09_29_R165_T166.sig water
% L8_2014_09_29_R165_T167.sig water
% L8_2014_09_29_R165_T168.sig water
% L8_2014_09_29_R165_T169.sig sky
% L8_2014_09_29_R165_T170.sig sky


filename = 'L8_2014_09_29_R165_T166.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R165_T167.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R165_T168.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R165_T169.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskya,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R165_T170.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskyb,~] = extractSVC(filepath);


Rrs1a = (Lt1-rho.*Lskya)./(pi.*Lg./0.99);
Rrs2a = (Lt2-rho.*Lskya)./(pi.*Lg./0.99);
Rrs3a = (Lt3-rho.*Lskya)./(pi.*Lg./0.99);

Rrs1b = (Lt1-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs2b = (Lt2-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs3b = (Lt3-rho.*Lskyb)./(pi.*Lg./0.99);

RrsLONGS140929 = Rrs2a;

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

title('R_{rs} -- LONGS 2','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([400 1000 0 0.03])
grid on

%% Radiance Lsky
figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Lskya,'r')
hold on
plot(wavelengthSVC,Lskyb,'b')
legend('Lskya','Lskyb')

title('Lsky -- LONGS 2','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.015])
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LONGN 1
% L8_2014_09_29_R144_T145.sig water
% L8_2014_09_29_R144_T146.sig water
% L8_2014_09_29_R144_T147.sig water
% L8_2014_09_29_R144_T148.sig sky
% L8_2014_09_29_R144_T149.sig sky


filename = 'L8_2014_09_29_R144_T145.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R144_T146.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R144_T147.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R144_T148.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskya,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R144_T149.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskyb,~] = extractSVC(filepath);


Rrs1a = (Lt1-rho.*Lskya)./(pi.*Lg./0.99);
Rrs2a = (Lt2-rho.*Lskya)./(pi.*Lg./0.99);
Rrs3a = (Lt3-rho.*Lskya)./(pi.*Lg./0.99);

Rrs1b = (Lt1-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs2b = (Lt2-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs3b = (Lt3-rho.*Lskyb)./(pi.*Lg./0.99);

RrsLONGN140929 = Rrs2a;

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
plot(wavelengthSVC,RrsLONGN140929,'g','linewidth',1.5)
legend('Lt1 Lskya','Lt2 Lskya','Lt3 Lskya',...
    'Lt1 Lskyb','Lt2 Lskyb','Lt3 Lskyb')

title('R_{rs} -- LONGN 1','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([400 1000 0 0.03])
grid on

%% Radiance Lsky
figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Lskya,'r')
hold on
plot(wavelengthSVC,Lskyb,'b')
legend('Lskya','Lskyb')

title('Lsky -- LONGN 1','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.015])
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LONGN 2
% L8_2014_09_29_R150_T151.sig water
% L8_2014_09_29_R150_T152.sig water
% L8_2014_09_29_R150_T153.sig water
% L8_2014_09_29_R150_T154.sig sky
% L8_2014_09_29_R150_T155.sig sky
% L8_2014_09_29_R150_T156.sig sky


filename = 'L8_2014_09_29_R150_T151.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R150_T152.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R150_T153.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R150_T154.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskya,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R150_T155.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskyb,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R150_T156.sig';
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

title('R_{rs} -- LONGN 2','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([400 1000 0 0.03])
grid on

%% Radiance Lsky
figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Lskya,'r')
hold on
plot(wavelengthSVC,Lskyb,'b')
plot(wavelengthSVC,Lskyb,'c')

legend('Lskya','Lskyb','Lskyc')

title('Lsky -- LONGN 2','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.015])
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CRANB 1
% L8_2014_09_29_R173_T174.sig water
% L8_2014_09_29_R173_T175.sig water
% L8_2014_09_29_R173_T176.sig water
% L8_2014_09_29_R173_T177.sig sky
% L8_2014_09_29_R173_T178.sig sky


filename = 'L8_2014_09_29_R173_T174.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R173_T175.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R173_T176.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R173_T177.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskya,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R173_T178.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskyb,~] = extractSVC(filepath);


Rrs1a = (Lt1-rho.*Lskya)./(pi.*Lg./0.99);
Rrs2a = (Lt2-rho.*Lskya)./(pi.*Lg./0.99);
Rrs3a = (Lt3-rho.*Lskya)./(pi.*Lg./0.99);

Rrs1b = (Lt1-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs2b = (Lt2-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs3b = (Lt3-rho.*Lskyb)./(pi.*Lg./0.99);

RrsCRANB140929 = Rrs3a;

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

title('R_{rs} -- CRANB 1','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([400 1000 0 0.03])
grid on

%% Radiance Lsky
figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Lskya,'r')
hold on
plot(wavelengthSVC,Lskyb,'b')
legend('Lskya','Lskyb')

title('Lsky -- CRANB 1','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.015])
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CRANB 2
% L8_2014_09_29_R179_T180.sig water
% L8_2014_09_29_R179_T181.sig water
% L8_2014_09_29_R179_T182.sig water
% L8_2014_09_29_R179_T183.sig sky
% L8_2014_09_29_R179_T184.sig sky


filename = 'L8_2014_09_29_R179_T180.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R179_T181.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R179_T182.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R179_T183.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskya,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R179_T184.sig';
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

title('R_{rs} -- CRANB 2','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([400 1000 0 0.03])
grid on

%% Radiance Lsky
figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Lskya,'r')
hold on
plot(wavelengthSVC,Lskyb,'b')
legend('Lskya','Lskyb')

title('Lsky -- CRANB 2','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.015])
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IBAYN 1
% L8_2014_09_29_R198_T199.sig water
% L8_2014_09_29_R198_T200.sig water
% L8_2014_09_29_R198_T201.sig water
% L8_2014_09_29_R198_T202.sig sky
% L8_2014_09_29_R198_T203.sig sky


filename = 'L8_2014_09_29_R198_T199.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R198_T200.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R198_T201.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R198_T202.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskya,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R198_T203.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskyb,~] = extractSVC(filepath);

Rrs1a = (Lt1-rho.*Lskya)./(pi.*Lg./0.99);
Rrs2a = (Lt2-rho.*Lskya)./(pi.*Lg./0.99);
Rrs3a = (Lt3-rho.*Lskya)./(pi.*Lg./0.99);

Rrs1b = (Lt1-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs2b = (Lt2-rho.*Lskyb)./(pi.*Lg./0.99);
Rrs3b = (Lt3-rho.*Lskyb)./(pi.*Lg./0.99);

RrsIBAYN140929 = Rrs1a;

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

title('R_{rs} -- IBAYN 1','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([400 1000 0 0.03])
grid on

Ref = [wavelengthSVC*1E-3, RrsIBAYN140929];
save('IBAYNRef_140919.txt','Ref','-ascii')
%% Radiance Lsky
figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Lskya,'r')
hold on
plot(wavelengthSVC,Lskyb,'b')

legend('Lskya','Lskyb')

title('Lsky -- IBAYN 1','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.015])
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IBAYN 2
% L8_2014_09_29_R204_T205.sig water
% L8_2014_09_29_R204_T206.sig water
% L8_2014_09_29_R204_T207.sig water
% L8_2014_09_29_R204_T208.sig sky
% L8_2014_09_29_R204_T209.sig sky
% L8_2014_09_29_R204_T210.sig sky


filename = 'L8_2014_09_29_R204_T205.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R204_T206.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R204_T207.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R204_T208.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskya,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R204_T209.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lskyb,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R204_T210.sig';
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

% RrsIBAYN140929 = Rrs1a;

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

title('R_{rs} -- IBAYN 2','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
axis([400 1000 0 0.03])
grid on

% Ref = [wavelengthSVC*1E-3, RrsIBAYN140929];
% save('IBAYNRef_140919.txt','Ref','-ascii')
%% Radiance Lsky
figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Lskya,'r')
hold on
plot(wavelengthSVC,Lskyb,'b')
plot(wavelengthSVC,Lskyb,'c')

legend('Lskya','Lskyb','Lskyc')

title('Lsky -- IBAYN 2','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.015])
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Rrs all
figure(21)
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,RrsONTOS140929,'r')
hold on
plot(wavelengthSVC,RrsLONGS140929,'b')
plot(wavelengthSVC,RrsLONGN140929,'g')
plot(wavelengthSVC,RrsCRANB140929,'m')
plot(wavelengthSVC,RrsIBAYN140929,'k')

legend('ONTOS','LONGS','LONGN','CRANB','IBAYN')
title('R_{rs} -- 09/29/14 ','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.03])
grid on
%% Spectrally sampled and save in text file
wlrange = wavelengthSVC>=400 & wavelengthSVC<=2500;
% wlzero = wavelengthSVC==2219.0;
wlavg = wavelengthSVC>=860 & wavelengthSVC<=870;

zeroavg = mean(RrsONTOS140929(wlavg));

RrsONTOScorr = RrsONTOS140929-zeroavg; % subtract NIR value for all bands
% RrsONTOScorr = RrsONTOS140929;%-zeroavg; % force NIR and SWIR bands to zero

RrsONTOSL8 = spect_sampL8(RrsONTOScorr(wlrange),wavelengthSVC(wlrange).*1E-3);

RrsONTOSL8corr = RrsONTOSL8;
RrsONTOSL8corr(5:7)=0;

L8bands = [0.4430,0.4826,0.5613,0.6546,0.8646,1.6090,2.2010];

figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,RrsONTOS140929,'r')
hold on
plot(wavelengthSVC,RrsONTOScorr,'--r')
% plot(wavelengthSVC(wlzero),RrsONTOS(wlzero),'.g')
plot(L8bands.*1E3,RrsONTOSL8,'.-b')
plot(L8bands.*1E3,RrsONTOSL8corr,'.-k')
plot(L8bands.*1E3,RrsONTOSL8corr*pi,'--k')
legend('RrsONTOS140929','RrsONTOScorr','RrsONTOSL8','RrsONTOSL8corr','RrsONTOSL8corr*pi')
title('R_{rs} -- 09/29/14 ','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.03])
grid on
% 
% NSRef = [L8bands' RrsONTNSL8corr'];
% save([pathname,pathdate,'RrsONTNSL8.txt'],'NSRef','-ascii')

Ref = [L8bands', RrsONTOSL8corr'];
save('ONTOSL8_Ref_140919corr.txt','Ref','-ascii')

%% Find best match in the HL LUT -- Run the LUT part of retrievalL8_140929.m first...
Rrs_SITE_test = RrsONTOS140929;

Rrs_SITE_test_HL = interp1(wavelengthSVC,Rrs_SITE_test,wavelength*1000);
Rrs_SITE_test_HL = Rrs_SITE_test_HL-Rrs_SITE_test_HL(end);

wl_um = wavelength*1000;

cond1 = wl_um>500;

Rrs_test = Rrs(cond1,:);

figure
plot(wl_um,Rrs_SITE_test_HL)
xlim([400 1000])

hold on

[Y,I1] = min(sqrt(mean((Rrs_test'-ones(size(Rrs_test,2),1)*Rrs_SITE_test_HL(cond1)').^2,2)));

plot(wl_um(cond1),Rrs_test(:,I1),'g')

legend('Field','LUT')

str = sprintf('%s %f %f %f %s',char(c{1}(I1)),c{2}(I1),c{3}(I1),c{4}(I1),char(c{5}(I1)));
title(str)
grid on
%%
figure
plot(wl_um,Rrs)
hold on
plot(wl_um(cond1),Rrs_test(:,I1),'g','linewidth',1.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAND 1
% L8_2014_09_29_R211_T212.sig
% L8_2014_09_29_R211_T213.sig
% L8_2014_09_29_R211_T214.sig
% L8_2014_09_29_R211_T215.sig


filename = 'L8_2014_09_29_R211_T212.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R211_T213.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R211_T214.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R211_T215.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt4,~] = extractSVC(filepath);


Rrs1 = (Lt1)./(pi.*Lg./0.99);
Rrs2 = (Lt2)./(pi.*Lg./0.99);
Rrs3 = (Lt3)./(pi.*Lg./0.99);
Rrs4 = (Lt4)./(pi.*Lg./0.99);

figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Rrs1,'r')
hold on
plot(wavelengthSVC,Rrs2,'b')
plot(wavelengthSVC,Rrs3,'g')
plot(wavelengthSVC,Rrs4,'--m')
legend('Lt1','Lt2','Lt3',...
    'Lt4')

title('R_{rs} -- SAND 1 ','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
% axis([0 3000 -0.01 0.01])
% axis([400 1000 0 0.015])
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAND 2
% L8_2014_09_29_R217_T218.sig
% L8_2014_09_29_R217_T219.sig
% L8_2014_09_29_R217_T220.sig
% L8_2014_09_29_R217_T221.sig



filename = 'L8_2014_09_29_R217_T218.sig';
filepath = [pathname,pathdate,filename];
[wavelengthSVC,Lg,Lt1,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R217_T219.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt2,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R217_T220.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt3,~] = extractSVC(filepath);

filename = 'L8_2014_09_29_R217_T221.sig';
filepath = [pathname,pathdate,filename];
[~,~,Lt4,~] = extractSVC(filepath);


Rrs1 = (Lt1)./(pi.*Lg./0.99);
Rrs2 = (Lt2)./(pi.*Lg./0.99);
Rrs3 = (Lt3)./(pi.*Lg./0.99);
Rrs4 = (Lt4)./(pi.*Lg./0.99);

figure
fs = 15;
set(gcf,'color','white')
plot(wavelengthSVC,Rrs1,'r')
hold on
plot(wavelengthSVC,Rrs2,'b')
plot(wavelengthSVC,Rrs3,'g')
plot(wavelengthSVC,Rrs4,'--m')
legend('Lt1','Lt2','Lt3',...
    'Lt4')

title('R_{rs} -- SAND 2 ','fontsize',fs)
xlabel('wavelengthSVC [nm]','fontsize',fs)
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
% xlabel('wavelengthSVC [\mu m]','fontsize',fs)
% ylabel('reflectance','fontsize',fs)
% 
% 
% NSRef = [wavelengthSVC, ONTNSRef];
% save('ONTNSRef.txt','NSRef','-ascii')
