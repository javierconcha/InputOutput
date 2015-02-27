addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/retrieval')
cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929
%% Loading
% CHL
CH219071 	= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/CH219071.ASC');
CH230251 	= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/CH230251.ASC');
CH240091 	= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/CH240091.ASC');
CH240092 	= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/CH240092.ASC');
CH246071 	= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/CH246071.ASC');
% SM
SM219071 	= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/SM219071.ASC');
SM229071 	= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/SM229071.ASC');
SM230251 	= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/SM230251.ASC');
SM240091 	= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/SM240091.ASC');
SM246071 	= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/SM246071.ASC');

SMBASE1 	= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/SMBASE1.ASC');
SMBASE2 	= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/SMBASE2.ASC');

%% SM and Chl Absorption Coefficients

[a_SM_LONGN,a_CH_LONGN,wavelength] = ParticleExtraction(SM219071,CH219071,70); % LONGN
[a_SM_LONGS,a_CH_LONGS,~] = ParticleExtraction(SM246071,CH246071,70); % LONGS
[a_SM_CRANB,a_CH_CRANB,~] = ParticleExtraction(SM229071,CH219071,70); % CRANB!!!! I do NOT have CH229071
[a_SM_ONTOS,a_CH_ONTOS,~] = ParticleExtraction(SM230251,CH230251,250); % ONTOS
[a_SM_IBAYN,a_CH_IBAYN,~] = ParticleExtraction(SM240091,CH240092,90); % IBAYN

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,a_SM_LONGN)
hold on
plot(wavelength,a_SM_LONGS,'k')
plot(wavelength,a_SM_CRANB,'r')
plot(wavelength,a_SM_ONTOS,'g')
plot(wavelength,a_SM_IBAYN,'m')
plot([wavelength(1) wavelength(end)],[0 0],'k')
title('SM Absorption Coefficient','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGN','LONGS','CRANB','ONTOS','IBAYN');
xlim([400 900])

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,a_CH_LONGN)
hold on
plot(wavelength,a_CH_LONGS,'k')
plot(wavelength,a_CH_CRANB,'r')
plot(wavelength,a_CH_ONTOS,'g')
plot(wavelength,a_CH_IBAYN,'m')
plot([wavelength(1) wavelength(end)],[0 0],'k')
title('CH Absorption Coefficient','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGN','LONGS','CRANB','ONTOS','IBAYN')
xlim([400 900])

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
% X_CH_LONGN = 41.65;%[ug/L] or [mg/m^3]
% X_CH_LONGS = 41.65;
% X_CH_CRANB = 49.66;
% X_CH_ONTOS = 2.40;
% X_CH_IBAYN = 28.20;

% MONROE County lab Corrected
% X_CH_LONGN = 28.10;%[ug/L] or [mg/m^3]
% X_CH_LONGS = 24.10;
% X_CH_CRANB = 20.10;
% X_CH_ONTOS = 1.00;
% X_CH_IBAYN = 20.10;;

% Me 10200HC Uncorrected
% X_CH_LONGN = 48.37;%[ug/L] or [mg/m^3]
% X_CH_LONGS = 48.37;
% X_CH_CRANB = 52.73;
% X_CH_ONTOS = 2.64;
% X_CH_IBAYN = 31.69;

% MONROE County lab Uncorrected
X_CH_LONGN = 47.90;% [ug/L] or [mg/m^3]
X_CH_LONGS = 46.10;
X_CH_CRANB = 58.30;
X_CH_ONTOS = 2.10;
X_CH_IBAYN = 28.30;

astar_CH_LONGN = a_CH_LONGN/X_CH_LONGN;%[m^2/mg]
astar_CH_LONGS = a_CH_LONGS/X_CH_LONGS;
astar_CH_CRANB = a_CH_CRANB/X_CH_CRANB; 
astar_CH_ONTOS = a_CH_ONTOS/X_CH_ONTOS;
astar_CH_IBAYN = a_CH_IBAYN/X_CH_IBAYN;

figure % in a*, [m^2/mg]
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_CH_LONGN)
hold on
plot(wavelength,astar_CH_LONGS,'k')
plot(wavelength,astar_CH_CRANB,'r')
plot(wavelength,astar_CH_ONTOS,'g')
plot(wavelength,astar_CH_IBAYN,'c')
plot([wavelength(1) wavelength(end)],[0 0],'k')
title('CH mass-specific absorption coeff.','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('a*, [m^2/mg]','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGN','LONGS','CRANB','ONTOS','IBAYN')
xlim([400 900])


%% Correcting ONTOS and LONGS

astar_CH_ONTOS_corrected = astar_CH_ONTOS;
astar_CH_ONTOS_corrected(wavelength>=750)=0;
astar_CH_ONTOS_corrected = sgolayfilt(astar_CH_ONTOS_corrected,5,41);
astar_CH_ONTOS_corrected(astar_CH_ONTOS_corrected<0)=0;

astar_CH_LONGS_corrected = astar_CH_LONGS;
astar_CH_LONGS_corrected(wavelength>=750)=0;
astar_CH_LONGS_corrected = sgolayfilt(astar_CH_LONGS_corrected,5,41);
astar_CH_LONGS_corrected(astar_CH_LONGS_corrected<0)=0;


figure % in a*, [m^2/mg]
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_CH_LONGS_corrected,'r','linewidth',1.2)
hold on
plot(wavelength,astar_CH_ONTOS_corrected,'b','linewidth',1.2)
% plot(wavelength,astar_CH_LONGS,'k')
% hold on
% plot([wavelength(1) wavelength(end)],[0 0],'k')
grid on
title('CH mass-specific absorption coeff. -- LONGS','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('a*, [m^2/mg]','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGS','ONTOS')
% xlim([400 900])

%% saving to be used in HL
astar = [wavelength(end:-1:1) astar_CH_ONTOS_corrected(end:-1:1)];
save('astar_CH_ONTOS140929_CountyUncorr.txt','-ascii','astar')

astar = [wavelength(end:-1:1) astar_CH_LONGS_corrected(end:-1:1)];
save('astar_CH_LONGS140929_CountyUncorr.txt','-ascii','astar')

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
X_SM_LONGN = 16.70;
X_SM_LONGS = 28.30;
X_SM_CRANB = 30.70;
X_SM_ONTOS = 1.40;
X_SM_IBAYN = 9.11;

astar_SM_LONGN = a_SM_LONGN/X_SM_LONGN;%[m^2/g]
astar_SM_LONGS = a_SM_LONGS/X_SM_LONGS;
astar_SM_CRANB = a_SM_CRANB/X_SM_CRANB; 
astar_SM_ONTOS = a_SM_ONTOS/X_SM_ONTOS;
astar_SM_IBAYN = a_SM_IBAYN/X_SM_IBAYN;

figure % in a*, [m^2/g]
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_SM_LONGN)
hold on
plot(wavelength,astar_SM_LONGS,'k')
plot(wavelength,astar_SM_CRANB,'r')
plot(wavelength,astar_SM_ONTOS,'g')
plot(wavelength,astar_SM_IBAYN,'m')
plot([wavelength(1) wavelength(end)],[0 0],'k')
title('SM mass-specific absorption coeff. -- RIT','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('a*, [m^2/g]','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGN','LONGS','CRANB','ONTOS','IBAYN')
xlim([400 900])


%% Correcting SM 

astar_SM_ONTOS_corrected = astar_SM_ONTOS-astar_SM_ONTOS(wavelength==753);
astar_SM_ONTOS_corrected(wavelength>=814)=0;
astar_SM_ONTOS_corrected = sgolayfilt(astar_SM_ONTOS_corrected,5,41);
astar_SM_ONTOS_corrected(astar_SM_ONTOS_corrected<0)=0;

astar_SM_LONGS_corrected = astar_SM_LONGS-astar_SM_LONGS(wavelength==814);
astar_SM_LONGS_corrected(wavelength>=814)=0;
astar_SM_LONGS_corrected = sgolayfilt(astar_SM_LONGS_corrected,5,41);
astar_SM_LONGS_corrected(astar_SM_LONGS_corrected<0)=0;

%% Fitting exp to SM (from Curve Fitting apps in Matlab)
a =       8.219;
b =    -0.01428;

SM_ONTOS_fitted = a*exp(b.*wavelength);

a =       3.905;
b =    -0.01036;

SM_LONGS_fitted = a*exp(b.*wavelength);

figure % in a*, [m^2/mg]
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_SM_LONGS_corrected,'--r','linewidth',1)
hold on
plot(wavelength,astar_SM_ONTOS_corrected,'--b','linewidth',1)
plot(wavelength,SM_LONGS_fitted,'r','linewidth',1)
plot(wavelength,SM_ONTOS_fitted,'b','linewidth',1)
grid on
title('SM mass-specific absorption coeff. -- LONGS','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('a*, [m^2/g]','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGS','ONTOS','LONGS fit','ONTOS fit')
% xlim([400 900])

%% save the SM to be used in HL

astar = [wavelength(end:-1:1) SM_ONTOS_fitted(end:-1:1)];
save('astar_SM_ONTOS140929_County.txt','-ascii','astar')

astar = [wavelength(end:-1:1) SM_LONGS_fitted(end:-1:1)];
save('astar_SM_LONGS140929_County.txt','-ascii','astar')
%% CDOM absorption coefficient

% Loading
CD2001		= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/CD2001.ASC'); % for 140828 
CD2131		= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/CD2131.ASC'); % ONTOS
CD2191		= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/CD2191.ASC'); % LONGN
CD2251		= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/CD2251.ASC'); %for 140828
CD2291		= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/CD2291.ASC'); % CRANB
CD2401		= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/CD2401.ASC'); % IBAYN
CD2461		= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/CD2461.ASC'); % LONGS
CDBASE1 	= load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/CDBASE1.ASC');

astar_CD_LONGN =  CdomExtraction(CD2191);
astar_CD_LONGS =  CdomExtraction(CD2461);
astar_CD_CRANB =  CdomExtraction(CD2291);
astar_CD_ONTOS =  CdomExtraction(CD2131);
astar_CD_IBAYN =  CdomExtraction(CD2401);

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_CD_LONGN,'k');
hold on
plot(wavelength,astar_CD_LONGS,'r');
plot(wavelength,astar_CD_CRANB,'g');
plot(wavelength,astar_CD_ONTOS,'c');
plot(wavelength,astar_CD_IBAYN,'--k');


plot([wavelength(1) wavelength(end)],[0 0],'k')
title('CDOM Absorption Coefficient','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGN','LONGS','CRANB','ONTOS','IBAYN');
xlim([400 900])
grid on
%% Correcting ONTOS and LONGS

astar_CD_ONTOS_corrected = astar_CD_ONTOS-astar_CD_ONTOS(wavelength==726);
astar_CD_ONTOS_corrected(wavelength>=726)=0;
astar_CD_ONTOS_corrected = sgolayfilt(astar_CD_ONTOS_corrected,5,41);
astar_CD_ONTOS_corrected(astar_CD_ONTOS_corrected<0)=0;

astar_CD_LONGS_corrected = astar_CD_LONGS-astar_CD_LONGS(wavelength==750);
astar_CD_LONGS_corrected(wavelength>=750)=0;
astar_CD_LONGS_corrected = sgolayfilt(astar_CD_LONGS_corrected,5,41);
astar_CD_LONGS_corrected(astar_CD_LONGS_corrected<0)=0;

%% fitting exp to ONTOS

a =       178.2;
b =    -0.01712;

CD_ONTOS_fitted = a*exp(b.*wavelength);


a =       850.9;
b =    -0.01545;

CD_LONGS_fitted = a*exp(b.*wavelength);

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_CD_ONTOS ,'k');
hold on
plot(wavelength,astar_CD_ONTOS_corrected ,'b');
plot(wavelength,CD_ONTOS_fitted ,'r');
% plot([wavelength(1) wavelength(end)],[0 0],'k')
title('CDOM Absorption Coefficient ONTOS','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
xlim([400 900])
grid on

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_CD_LONGS ,'k');
hold on
plot(wavelength,astar_CD_LONGS_corrected ,'b');
plot(wavelength,CD_LONGS_fitted ,'r');
% plot([wavelength(1) wavelength(end)],[0 0],'k')
title('CDOM Absorption Coefficient LONGS','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
xlim([400 900])
grid on

%% CDOM absorption coefficient normalized at a(440)=1

astar_CD_LONGS440 =  CD_LONGS_fitted./CD_LONGS_fitted(wavelength==440);
astar_CD_ONTOS440 =  CD_ONTOS_fitted./CD_ONTOS_fitted(wavelength==440);

disp('LONGS(440):')
CD_LONGS_fitted(wavelength==440)

disp('ONTOS(440):')
CD_ONTOS_fitted(wavelength==440)


figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_CD_LONGS440,'k');
hold on
plot(wavelength,astar_CD_ONTOS440,'c');
title('CDOM Absorption Coefficient','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGS','ONTOS');
grid on

%% saving to be used in HL
astar = [wavelength(end:-1:1) astar_CD_LONGS440(end:-1:1)];
save('astar_CDOM_LONGS140929.txt','-ascii','astar')

astar = [wavelength(end:-1:1) astar_CD_ONTOS440(end:-1:1)];
save('astar_CDOM_ONTOS140929.txt','-ascii','astar')
