addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/retrieval')
cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916
date = '09-16-2015';
%% Loading

% CHL
CH201051 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CH201051.ASC'); % CRANB
CH206251 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CH206251.ASC'); % ONTOS
CH217051 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CH217051.ASC'); % LONGN
CH219251 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CH219251.ASC'); % RVRPL
CH236091 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CH236091.ASC'); % IBAYN
% SM
SM201051 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/SM201051.ASC'); % CRANB
SM206251 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/SM206251.ASC'); % ONTOS
SM217051 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/SM217051.ASC'); % LONGN
SM219251 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/SM219251.ASC'); % RVRPL
SM236091 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/SM236091.ASC'); % IBAYN

SMBASE1 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/SMBASE1.ASC');

%% SM and Chl Absorption Coefficients

[a_SM_CRANB,a_CH_CRANB,wavelength] = ParticleExtraction(SM201051,CH201051,50); % CRANB
[a_SM_ONTOS,a_CH_ONTOS,~] = ParticleExtraction(SM206251,CH206251,250); % ONTOS
[a_SM_LONGN,a_CH_LONGN,~] = ParticleExtraction(SM217051,CH217051,50); % LONGN
[a_SM_RVRPL,a_CH_RVRPL,~] = ParticleExtraction(SM219251,CH219251,250); % RVRPL
[a_SM_IBAYN,a_CH_IBAYN,~] = ParticleExtraction(SM236091,CH236091,90); % IBAYN



figure('name',date)
fs = 15;
set(gcf,'color','white')
plot(wavelength,a_SM_LONGN)
hold on
plot(wavelength,a_SM_RVRPL,'k')
plot(wavelength,a_SM_CRANB,'r')
plot(wavelength,a_SM_ONTOS,'g')
plot(wavelength,a_SM_IBAYN,'m')
plot([wavelength(1) wavelength(end)],[0 0],'k')
title('SM Absorption Coefficient','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGN','RVRPL','CRANB','ONTOS','IBAYN');
xlim([400 900])
grid on

figure('name',date)
fs = 15;
set(gcf,'color','white')
plot(wavelength,a_CH_LONGN)
hold on
plot(wavelength,a_CH_RVRPL,'k')
plot(wavelength,a_CH_CRANB,'r')
plot(wavelength,a_CH_ONTOS,'g')
plot(wavelength,a_CH_IBAYN,'m')
plot([wavelength(1) wavelength(end)],[0 0],'k')
title('CH Absorption Coefficient','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGN','RVRPL','CRANB','ONTOS','IBAYN')
xlim([400 900])
grid on


%% CH mass-specific absorption coefficient
% astar = a/X where a is absorption coeff. and X is the component concentration
% From HE User's guide: For chl, astar is the chl- or mass-specific absorption coefficient 
% with units of m^2*(mg Chl)^{-1}, and chl concentration in (mg Chl) m^3.  
% For mineral particles , astar is in m^2*(g)^{-1} and the mineral particle
% concentration in g*m^3
% For our case:
% X_{Chl}: [ug/L]
% a: [m^{-1}]
% => astar: 0.001[m^2/ug] 

% Dr. Cristy Corrected
% X_CH_CRANB = ;%[ug/L] or [mg/m^3]
% X_CH_ONTOS = ;
% X_CH_LONGN = ;
% X_CH_RVRPL = ;
% X_CH_IBAYN = ;

% MONROE County lab Corrected
% X_CH_CRANB = 64.2;%[ug/L] or [mg/m^3]
% X_CH_ONTOS = 1.00; % <1.0
% X_CH_LONGN = 58.1;
% X_CH_RVRPL = 1.34;
% X_CH_IBAYN = 14.3;

% Me 10200HC Uncorrected
% X_CH_CRANB = ;%[ug/L] or [mg/m^3]
% X_CH_ONTOS = ;
% X_CH_LONGN = ;
% X_CH_RVRPL = ;
% X_CH_IBAYN = ;

% MONROE County lab Uncorrected
X_CH_CRANB = 99.3;% [ug/L] or [mg/m^3]
X_CH_ONTOS = 1.7 ;
X_CH_LONGN = 116 ;
X_CH_RVRPL = 3.91;
X_CH_IBAYN = 19.6;

astar_CH_CRANB = a_CH_CRANB/X_CH_CRANB;%[m^2/mg]
astar_CH_ONTOS = a_CH_ONTOS/X_CH_ONTOS;
astar_CH_LONGN = a_CH_LONGN/X_CH_LONGN; 
astar_CH_RVRPL = a_CH_RVRPL/X_CH_RVRPL;
astar_CH_IBAYN = a_CH_IBAYN/X_CH_IBAYN;

figure('name',date) % in a*, [m^2/mg]
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_CH_CRANB)
hold on
plot(wavelength,astar_CH_ONTOS,'k')
plot(wavelength,astar_CH_LONGN,'r')
plot(wavelength,astar_CH_RVRPL,'g')
plot(wavelength,astar_CH_IBAYN,'c')
plot([wavelength(1) wavelength(end)],[0 0],'k')
title('CH mass-specific absorption coeff.','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('a*, [m^2/mg]','fontsize',fs)
set(gca,'fontsize',fs)
legend('CRANB','ONTOS','LONGN','RVRPL','IBAYN');
xlim([400 900])
grid on

%% Correcting ONTOS and LONGS and CRANB

astar_CH_ONTOS_corrected = astar_CH_ONTOS;
astar_CH_ONTOS_corrected(wavelength>=750)=0;
astar_CH_ONTOS_corrected = sgolayfilt(astar_CH_ONTOS_corrected,5,41);
astar_CH_ONTOS_corrected(astar_CH_ONTOS_corrected<0)=0;

astar_CH_CRANB_corrected = astar_CH_CRANB;
astar_CH_CRANB_corrected(wavelength>=750)=0;
astar_CH_CRANB_corrected = sgolayfilt(astar_CH_CRANB_corrected,5,41);
astar_CH_CRANB_corrected(astar_CH_CRANB_corrected<0)=0;


figure('name',date) % in a*, [m^2/mg]
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_CH_CRANB_corrected,'r','linewidth',1.2)
hold on
plot(wavelength,astar_CH_ONTOS_corrected,'b','linewidth',1.2)
% plot(wavelength,astar_CH_CRANB,'k')
% hold on
% plot([wavelength(1) wavelength(end)],[0 0],'k')
grid on
title('CH mass-specific absorption coeff. -- CRANB','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('a*, [m^2/mg]','fontsize',fs)
set(gca,'fontsize',fs)
legend('CRANB','ONTOS')
% xlim([400 900])

%% saving to be used in HL
astar = [wavelength(end:-1:1) astar_CH_ONTOS_corrected(end:-1:1)];
save('astar_CH_ONTOS150916_CountyUncorr.txt','-ascii','astar')

astar = [wavelength(end:-1:1) astar_CH_CRANB_corrected(end:-1:1)];
save('astar_CH_CRANB150916_CountyUncorr.txt','-ascii','astar')

%% SM mass-specific absorption coefficient
% astar = a/X where a is absorption coeff. and X is the component concentration
% From HE User's guide: For chl, astar is the chl- or mass-specific absorption coefficient 
% with units of m^2*(mg Chl)^{-1}, and chl concentration in (mg Chl) m^3.  
% For mineral particles , astar is in m^2*(g)^{-1} and the mineral particle
% concentration in g*m^3
% For our case:
% X_{Chl}: [ug/L]
% a: [m^{-1}]
% => astar: 0.001[m^2/ug] 

% County
X_SM_CRANB = 28.3;
X_SM_ONTOS = 4.00;
X_SM_LONGN = 29.2;
X_SM_RVRPL = 4.0 ;
X_SM_IBAYN = 4.0 ;

astar_SM_CRANB = a_SM_CRANB/X_SM_CRANB;%[m^2/g]
astar_SM_ONTOS = a_SM_ONTOS/X_SM_ONTOS;
astar_SM_LONGN = a_SM_LONGN/X_SM_LONGN; 
astar_SM_RVRPL = a_SM_RVRPL/X_SM_RVRPL;
astar_SM_IBAYN = a_SM_IBAYN/X_SM_IBAYN;

figure('name',date) % in a*, [m^2/g]
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_SM_CRANB)
hold on
plot(wavelength,astar_SM_ONTOS,'k')
plot(wavelength,astar_SM_LONGN,'r')
plot(wavelength,astar_SM_RVRPL,'g')
plot(wavelength,astar_SM_IBAYN,'m')
plot([wavelength(1) wavelength(end)],[0 0],'k')
title('SM mass-specific absorption coeff. -- RIT','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('a*, [m^2/g]','fontsize',fs)
set(gca,'fontsize',fs)
legend('CRANB','ONTOS','LONGN','RVRPL','IBAYN')
xlim([400 900])
grid on

%% Fitting exp to SM (from Curve Fitting apps in Matlab)

astar_SM_CRANBoffset = astar_SM_CRANB - min(astar_SM_CRANB);
astar_SM_ONTOSoffset = astar_SM_ONTOS - min(astar_SM_ONTOS);

%% CRANB
% General model Exp2:
%      f(x) = a*exp(b*x) + c*exp(d*x)
% Coefficients (with 95% confidence bounds):
       a =      0.6075; %(0.3778, 0.8373)
       b =   -0.006162; %(-0.006671, -0.005653)
       c =        25.1; %(-6.872, 57.08)
       d =    -0.01702; %(-0.02093, -0.0131)

% Goodness of fit:
%   SSE: 0.001325
%   R-square: 0.9927
%   Adjusted R-square: 0.9927
%   RMSE: 0.001633

astar_SM_CRANB_fitted = a*exp(b*wavelength) + c*exp(d*wavelength);

astar_SM_CRANB_corrected = astar_SM_CRANB_fitted-astar_SM_CRANB_fitted(wavelength==900);
astar_SM_CRANB_corrected = sgolayfilt(astar_SM_CRANB_corrected,5,41);
astar_SM_CRANB_corrected(astar_SM_CRANB_corrected<0)=0;

% ONTOS
% General model Exp2:
%      f(x) = a*exp(b*x) + c*exp(d*x)
% Coefficients (with 95% confidence bounds):
       a =       1.901;  %(0.1859, 3.616)
       b =    -0.01134;  %(-0.01438, -0.008295)
       c =     0.04126;  %(-0.009542, 0.09206)
       d =   -0.003326;  %(-0.004764, -0.001889)

% Goodness of fit:
%   SSE: 0.001045
%   R-square: 0.9622
%   Adjusted R-square: 0.962
%   RMSE: 0.00145


astar_SM_ONTOS_fitted = a*exp(b*wavelength) + c*exp(d*wavelength);

astar_SM_ONTOS_corrected = astar_SM_ONTOS_fitted-astar_SM_ONTOS_fitted(wavelength==900);
astar_SM_ONTOS_corrected = sgolayfilt(astar_SM_ONTOS_corrected,5,41);
astar_SM_ONTOS_corrected(astar_SM_ONTOS_corrected<0)=0;

figure('name',date) % in a*, [m^2/g]
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_SM_CRANB)
hold on
plot(wavelength,astar_SM_CRANBoffset,'g')
plot(wavelength,astar_SM_CRANB_fitted,'r')
plot(wavelength,astar_SM_CRANB_corrected,'b','linewidth',1.5)
plot([wavelength(1) wavelength(end)],[0 0],'k')
title('SM mass-specific absorption coeff. -- RIT -- first fitted','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('a*, [m^2/g]','fontsize',fs)
set(gca,'fontsize',fs)
legend('CRANB','w/ offset','Fitted','corrected')
xlim([400 900])
grid on

figure('name',date) % in a*, [m^2/g]
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_SM_ONTOS)
hold on
plot(wavelength,astar_SM_ONTOSoffset,'g')
plot(wavelength,astar_SM_ONTOS_fitted,'r')
plot(wavelength,astar_SM_ONTOS_corrected,'b','linewidth',1.5)
plot([wavelength(1) wavelength(end)],[0 0],'k')
title('SM mass-specific absorption coeff. -- RIT -- first fitted','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('a*, [m^2/g]','fontsize',fs)
set(gca,'fontsize',fs)
legend('ONTOS','w/ offset','Fitted','corrected')
xlim([400 900])
grid on

figure('name',date) % in a*, [m^2/g]
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_SM_ONTOS_corrected,'r','linewidth',1.5)
hold on
plot(wavelength,astar_SM_CRANB_corrected,'b','linewidth',1.5)
plot([wavelength(1) wavelength(end)],[0 0],'k')
title('SM mass-specific absorption coeff. -- RIT corrected','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('a*, [m^2/g]','fontsize',fs)
set(gca,'fontsize',fs)
legend('ONTOS','CRANB')
xlim([400 900])
grid on


%% save the SM to be used in HL
% 
astar = [wavelength(end:-1:1) astar_SM_ONTOS_corrected(end:-1:1)];
save('astar_SM_ONTOS150916_County.txt','-ascii','astar')

astar = [wavelength(end:-1:1) astar_SM_CRANB_corrected(end:-1:1)];
save('astar_SM_CRANB150916_County.txt','-ascii','astar')

%% CDOM absorption coefficient


CD2012 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CD2012.ASC'); % CRANB
CD2061 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CD2061.ASC'); % ONTOS
CD2171 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CD2171.ASC'); % LONGN
CD2191 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CD2191.ASC'); % RVRPL
CD2361 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CD2361.ASC'); % IBAYN
CD2381 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CD2381.ASC'); % CRANB2

CDBASE1 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CDBASE1.ASC');
a_CD_CRAN2 =  CdomExtraction(CD2012);
a_CD_ONTOS =  CdomExtraction(CD2061);
a_CD_LONGN =  CdomExtraction(CD2171);
a_CD_RVRPL =  CdomExtraction(CD2191);
a_CD_IBAYN =  CdomExtraction(CD2361);
a_CD_CRANB =  CdomExtraction(CD2381);

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,a_CD_CRANB,'k');
hold on
plot(wavelength,a_CD_ONTOS,'r');
plot(wavelength,a_CD_LONGN,'g');
plot(wavelength,a_CD_RVRPL,'c');
plot(wavelength,a_CD_IBAYN,'--k');
plot(wavelength,a_CD_CRAN2,'--r');
title('CDOM Absorption Coefficient','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('CRANB','ONTOS','LONGN','RVRPL','IBAYN','CRAN2');
xlim([400 900])
grid on

a_CD_ONTOSoffset = a_CD_ONTOS-min(a_CD_ONTOS);
a_CD_LONGNoffset = a_CD_LONGN-min(a_CD_LONGN);
a_CD_RVRPLoffset = a_CD_RVRPL-min(a_CD_RVRPL);
a_CD_IBAYNoffset = a_CD_IBAYN-min(a_CD_IBAYN);
a_CD_CRANBoffset = a_CD_CRANB-min(a_CD_CRANB);

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,a_CD_CRANBoffset,'k');
hold on
plot(wavelength,a_CD_ONTOSoffset,'r');
plot(wavelength,a_CD_LONGNoffset,'g');
plot(wavelength,a_CD_RVRPLoffset,'c');
plot(wavelength,a_CD_IBAYNoffset,'--k');
title('CDOM Absorption Coefficient -- offset','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('CRANB','ONTOS','LONGN','RVRPL','IBAYN');
xlim([400 900])
grid on

%% fitting exp
% ONTOS
% General model Exp2:
%      f(x) = a*exp(b*x) + c*exp(d*x)
% Coefficients (with 95% confidence bounds):
       a =       779.3; %(511.8, 1047)
       b =    -0.01982; %(-0.0207, -0.01895)
       c =      0.1902; %(0.1644, 0.2159)
       d =   -0.001709; %(-0.001894, -0.001524)

% Goodness of fit:
%   SSE: 0.02882
%   R-square: 0.9883
%   Adjusted R-square: 0.9882
%   RMSE: 0.007615

a_CD_ONTOSfitted = a*exp(b*wavelength) + c*exp(d*wavelength);

a_CD_ONTOSzero = a_CD_ONTOSfitted-a_CD_ONTOSfitted(wavelength==900);

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,a_CD_ONTOS,'k');
hold on
plot(wavelength,a_CD_ONTOSoffset,'r');
plot(wavelength,a_CD_ONTOSfitted,'g');
plot(wavelength,a_CD_ONTOSzero,'b');
title('CDOM Absorption Coefficient -- ONTOS','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('ONTOS','offset','fitted','zero');
xlim([400 900])
grid on

% CRANB
% General model Exp2:
%      f(x) = a*exp(b*x) + c*exp(d*x)
% Coefficients (with 95% confidence bounds):
       a =        5509;  %(5038, 5980)
       b =    -0.01956;  %(-0.01979, -0.01933)
       c =       1.786;  %(1.535, 2.037)
       d =   -0.004192;  %(-0.00439, -0.003994)

% Goodness of fit:
%   SSE: 0.06265
%   R-square: 0.9995
%   Adjusted R-square: 0.9995
%   RMSE: 0.01123

a_CD_CRANBfitted = a*exp(b*wavelength) + c*exp(d*wavelength);

a_CD_CRANBzero = a_CD_CRANBfitted-a_CD_CRANBfitted(wavelength==900);

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,a_CD_CRANB,'k');
hold on
plot(wavelength,a_CD_CRANBoffset,'r');
plot(wavelength,a_CD_CRANBfitted,'g');
plot(wavelength,a_CD_CRANBzero,'b');
title('CDOM Absorption Coefficient -- CRANB','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('CRANB','offset','fitted','zero');
xlim([400 900])
grid on

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,a_CD_ONTOSzero,'k');
hold on
plot(wavelength,a_CD_CRANBzero,'b');
title('CDOM Absorption Coefficient -- zero','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('ONTOS','CRANB');
xlim([400 900])
grid on

%% CDOM absorption coefficient normalized at a(440)=1
astar_CD_ONTOS440 =  a_CD_ONTOSzero./a_CD_ONTOSzero(wavelength==440);
astar_CD_LONGN440 =  a_CD_LONGNzero./a_CD_LONGNzero(wavelength==440);
astar_CD_RVRPL440 =  a_CD_RVRPLzero./a_CD_RVRPLzero(wavelength==440);
astar_CD_IBAYN440 =  a_CD_IBAYNzero./a_CD_IBAYNzero(wavelength==440);
astar_CD_CRANB440 =  a_CD_CRANBzero./a_CD_CRANBzero(wavelength==440);


disp('a_CD_ONTOS(440):')
disp(a_CD_ONTOSzero(wavelength==440))

disp('a_CD_LONGN(440):')
disp(a_CD_LONGNzero(wavelength==440))

disp('a_CD_RVRPL(440):')
disp(a_CD_RVRPLzero(wavelength==440))

disp('a_CD_IBAYN(440):')
disp(a_CD_IBAYNzero(wavelength==440))

disp('a_CD_CRANB(440):')
disp(a_CD_CRANBzero(wavelength==440))


disp('astar_CD_ONTOS(440):')
disp(astar_CD_ONTOS440(wavelength==440))

% disp('astar_CD_LONGN(440):')
% disp(astar_CD_LONGN440(wavelength==440))

% disp('astar_CD_RVRPL(440):')
% disp(astar_CD_RVRPL440(wavelength==440))

% disp('astar_CD_IBAYN(440):')
% disp(astar_CD_IBAYN440(wavelength==440))

disp('astar_CD_CRANB(440):')
disp(astar_CD_CRANB440(wavelength==440))

figure('name',date)
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_CD_ONTOS440,'k');
hold on
% plot(wavelength,astar_CD_LONGN440,'r');
% plot(wavelength,astar_CD_RVRPL440,'g');
% plot(wavelength,astar_CD_IBAYN440,'c');
plot(wavelength,astar_CD_CRANB440,'b');
title('CDOM Absorption Coefficient','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('ONTOS','CRANB');
grid on

%% saving to be used in HL
astar = [wavelength(end:-1:1) astar_CD_CRANB440(end:-1:1)];
save('astar_CDOM_CRANB150916.txt','-ascii','astar')

astar = [wavelength(end:-1:1) astar_CD_ONTOS440(end:-1:1)];
save('astar_CDOM_ONTOS150916.txt','-ascii','astar')
