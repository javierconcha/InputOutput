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

%% CDOM

CD2012 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CD2012.ASC'); % CRANB
CD2061 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CD2061.ASC'); % ONTOS
CD2171 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CD2171.ASC'); % LONGN
CD2191 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CD2191.ASC'); % RVRPL
CD2361 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CD2361.ASC'); % IBAYN
CD2381 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CD2381.ASC'); % CRANB2

CDBASE1 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/150916/CDBASE1.ASC');
astar_CD_CRANB =  CdomExtraction(CD2012);
astar_CD_ONTOS =  CdomExtraction(CD2061);
astar_CD_LONGN =  CdomExtraction(CD2171);
astar_CD_RVRPL =  CdomExtraction(CD2191);
astar_CD_IBAYN =  CdomExtraction(CD2361);
astar_CD_CRAN2 =  CdomExtraction(CD2381);

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_CD_CRANB,'k');
hold on
plot(wavelength,astar_CD_ONTOS,'r');
plot(wavelength,astar_CD_LONGN,'g');
plot(wavelength,astar_CD_RVRPL,'c');
plot(wavelength,astar_CD_IBAYN,'--k');
plot(wavelength,astar_CD_CRAN2,'--r');
title('CDOM Absorption Coefficient','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('CRANB','ONTOS','LONGN','RVRPL','IBAYN','CRAN2');
xlim([400 900])
grid on
