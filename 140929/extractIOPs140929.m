addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/retrieval')
cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929
date = '09-29-2014';
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

figure('name',date)
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
grid on

figure('name',date)
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

figure('name',date) % in a*, [m^2/mg]
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
grid on


%% Correcting ONTOS and LONGS and CRANB

astar_CH_ONTOS_corrected = astar_CH_ONTOS;
astar_CH_ONTOS_corrected(wavelength>=750)=0;
astar_CH_ONTOS_corrected = sgolayfilt(astar_CH_ONTOS_corrected,5,41);
astar_CH_ONTOS_corrected(astar_CH_ONTOS_corrected<0)=0;

astar_CH_LONGS_corrected = astar_CH_LONGS;
astar_CH_LONGS_corrected(wavelength>=750)=0;
astar_CH_LONGS_corrected = sgolayfilt(astar_CH_LONGS_corrected,5,41);
astar_CH_LONGS_corrected(astar_CH_LONGS_corrected<0)=0;


figure('name',date) % in a*, [m^2/mg]
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
% astar = [wavelength(end:-1:1) astar_CH_ONTOS_corrected(end:-1:1)];
% save('astar_CH_ONTOS140929_CountyUncorr.txt','-ascii','astar')
% 
% astar = [wavelength(end:-1:1) astar_CH_LONGS_corrected(end:-1:1)];
% save('astar_CH_LONGS140929_CountyUncorr.txt','-ascii','astar')

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

figure('name',date) % in a*, [m^2/g]
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
grid on


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

figure('name',date) % in a*, [m^2/mg]
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
% 
% astar = [wavelength(end:-1:1) SM_ONTOS_fitted(end:-1:1)];
% save('astar_SM_ONTOS140929_County.txt','-ascii','astar')
% 
% astar = [wavelength(end:-1:1) SM_LONGS_fitted(end:-1:1)];
% save('astar_SM_LONGS140929_County.txt','-ascii','astar')
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

figure('name',date)
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_CD_LONGN,'k');
hold on
plot(wavelength,astar_CD_LONGS,'r');
plot(wavelength,astar_CD_CRANB,'g');
plot(wavelength,astar_CD_ONTOS,'c');
plot(wavelength,astar_CD_IBAYN,'--k');
title('CDOM Absorption Coefficient','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGN','LONGS','CRANB','ONTOS','IBAYN');
xlim([400 900])
grid on

%% Correcting for NIR signal

astar_CD_ONTOS_corrected = astar_CD_ONTOS-astar_CD_ONTOS(wavelength==726);
astar_CD_ONTOS_corrected(wavelength>=726)=0;
astar_CD_ONTOS_corrected = sgolayfilt(astar_CD_ONTOS_corrected,5,41);
astar_CD_ONTOS_corrected(astar_CD_ONTOS_corrected<0)=0;

astar_CD_LONGS_corrected = astar_CD_LONGS-astar_CD_LONGS(wavelength==887);
astar_CD_LONGS_corrected(wavelength>=887)=0;
astar_CD_LONGS_corrected = sgolayfilt(astar_CD_LONGS_corrected,5,41);
astar_CD_LONGS_corrected(astar_CD_LONGS_corrected<0)=0;

astar_CD_LONGN_corrected = astar_CD_LONGN-astar_CD_LONGN(wavelength==895);
astar_CD_LONGN_corrected(wavelength>=895)=0;
astar_CD_LONGN_corrected = sgolayfilt(astar_CD_LONGN_corrected,5,41);
astar_CD_LONGN_corrected(astar_CD_LONGN_corrected<0)=0;

astar_CD_CRANB_corrected = astar_CD_CRANB-astar_CD_CRANB(wavelength==877);
astar_CD_CRANB_corrected(wavelength>=877)=0;
astar_CD_CRANB_corrected = sgolayfilt(astar_CD_CRANB_corrected,5,41);
astar_CD_CRANB_corrected(astar_CD_CRANB_corrected<0)=0;

astar_CD_IBAYN_corrected = astar_CD_IBAYN-astar_CD_IBAYN(wavelength==836);
astar_CD_IBAYN_corrected(wavelength>=836)=0;
astar_CD_IBAYN_corrected = sgolayfilt(astar_CD_IBAYN_corrected,5,41);
astar_CD_IBAYN_corrected(astar_CD_IBAYN_corrected<0)=0;

figure('name',date)
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_CD_LONGN_corrected,'k');
hold on
plot(wavelength,astar_CD_LONGS_corrected,'r');
plot(wavelength,astar_CD_CRANB_corrected,'g');
plot(wavelength,astar_CD_ONTOS_corrected,'c');
plot(wavelength,astar_CD_IBAYN_corrected,'--k');
title('CDOM Absorption Coefficient; corrected','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGN','LONGS','CRANB','ONTOS','IBAYN');
xlim([400 900])
grid on

%% fitting exp
% ONTOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a =      0.3976;
% b =   -0.007702;
% c =       378.2;
% d =    -0.01921;

% CD_ONTOS_fitted = a*exp(b*wavelength) + c*exp(d*wavelength);

a =       178.2;
b =    -0.01712;

CD_ONTOS_fitted = a*exp(b*wavelength);
% General model Exp1:
%      f(x) = a*exp(b*x)
% Coefficients (with 95% confidence bounds):
%        a =       178.2  (169, 187.5)
%        b =    -0.01712  (-0.01724, -0.01699)

% Goodness of fit:
%   SSE: 0.00172
%   R-square: 0.9979
%   Adjusted R-square: 0.9979
%   RMSE: 0.001856

% General model Exp2:
%      f(x) = a*exp(b*x) + c*exp(d*x)
% Coefficients (with 95% confidence bounds):
%        a =      0.3976  (-0.1467, 0.9419)
%        b =   -0.007702  (-0.009626, -0.005778)
%        c =       378.2  (285.1, 471.2)
%        d =    -0.01921  (-0.01996, -0.01845)

% Goodness of fit:
%   SSE: 0.001136
%   R-square: 0.9986
%   Adjusted R-square: 0.9986
%   RMSE: 0.001512

% LONGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a =   2.412e+04;
% b =    -0.02488;
% c =       23.95;
% d =    -0.00874;

% CD_LONGS_fitted = a*exp(b*wavelength) + c*exp(d*wavelength);

a =       581.8;
b =    -0.01451;

CD_LONGS_fitted = a*exp(b*wavelength);


% General model Exp1:
%      f(x) = a*exp(b*x)
% Coefficients (with 95% confidence bounds):
%        a =       581.8  (538.9, 624.7)
%        b =    -0.01451  (-0.01468, -0.01434)

% Goodness of fit:
%   SSE: 0.474
%   R-square: 0.9937
%   Adjusted R-square: 0.9937
%   RMSE: 0.03082

% General model Exp2:
%      f(x) = a*exp(b*x) + c*exp(d*x)
% Coefficients (with 95% confidence bounds):
%        a =   2.412e+04  (2.068e+04, 2.757e+04)
%        b =    -0.02488  (-0.02528, -0.02448)
%        c =       23.95  (21.86, 26.04)
%        d =    -0.00874  (-0.008877, -0.008604)

% Goodness of fit:
%   SSE: 0.01142
%   R-square: 0.9998
%   Adjusted R-square: 0.9998
%   RMSE: 0.004794

% LONGN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a =   1.601e+04;
% b =    -0.02347;
% c =       9.859;
% d =   -0.007053;

% CD_LONGN_fitted = a*exp(b*wavelength) + c*exp(d*wavelength);

a =       472.1;
b =    -0.01395;

CD_LONGN_fitted = a*exp(b*wavelength);
% General model Exp1:
%      f(x) = a*exp(b*x)
% Coefficients (with 95% confidence bounds):
%        a =       472.1  (429.7, 514.5)
%        b =    -0.01395  (-0.01416, -0.01375)

% Goodness of fit:
%   SSE: 0.8069
%   R-square: 0.9896
%   Adjusted R-square: 0.9896
%   RMSE: 0.04021

% General model Exp2:
%      f(x) = a*exp(b*x) + c*exp(d*x)
% Coefficients (with 95% confidence bounds):
%        a =   1.601e+04  (1.485e+04, 1.716e+04)
%        b =    -0.02347  (-0.02367, -0.02328)
%        c =       9.859  (9.347, 10.37)
%        d =   -0.007053  (-0.007131, -0.006974)

% Goodness of fit:
%   SSE: 0.006542
%   R-square: 0.9999
%   Adjusted R-square: 0.9999
% 	RMSE: 0.003628

% CRANB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a =   3.773e+08;
% b =    -0.05318;
% c =        1185;
% d =    -0.01636;

% CD_CRANB_fitted = a*exp(b*wavelength) + c*exp(d*wavelength);

a =        1982;
b =    -0.01742;

CD_CRANB_fitted = a*exp(b*wavelength);
% General model Exp1:
%      f(x) = a*exp(b*x)
% Coefficients (with 95% confidence bounds):
%        a =        1982  (1926, 2039)
%        b =    -0.01742  (-0.01748, -0.01735)

% Goodness of fit:
%   SSE: 0.0485
%   R-square: 0.9994
%   Adjusted R-square: 0.9994
%   RMSE: 0.009859

% General model Exp2:
%      f(x) = a*exp(b*x) + c*exp(d*x)
% Coefficients (with 95% confidence bounds):
%        a =   3.773e+08  (-4.562e+08, 1.211e+09)
%        b =    -0.05318  (-0.05894, -0.04741)
%        c =        1185  (1107, 1262)
%        d =    -0.01636  (-0.01648, -0.01623)

% Goodness of fit:
%   SSE: 0.01401
%   R-square: 0.9998
%   Adjusted R-square: 0.9998
%   RMSE: 0.005309

% IBAYN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a =   5.309e+04;
% b =     -0.0268;
% c =       75.98;
% d =    -0.01114;

% CD_IBAYN_fitted = a*exp(b*wavelength) + c*exp(d*wavelength);

a =        1557;
b =     -0.0167;

CD_IBAYN_fitted = a*exp(b*wavelength);
% General model Exp1:
%      f(x) = a*exp(b*x)
% Coefficients (with 95% confidence bounds):
%        a =        1557  (1463, 1650)
%        b =     -0.0167  (-0.01683, -0.01656)

% Goodness of fit:
%   SSE: 0.2641
%   R-square: 0.997
%   Adjusted R-square: 0.9969
%   RMSE: 0.02301

% General model Exp2:
%      f(x) = a*exp(b*x) + c*exp(d*x)
% Coefficients (with 95% confidence bounds):
%        a =   5.309e+04  (4.225e+04, 6.393e+04)
%        b =     -0.0268  (-0.0274, -0.02621)
%        c =       75.98  (65.77, 86.18)
%        d =    -0.01114  (-0.01136, -0.01093)

% Goodness of fit:
%   SSE: 0.01187
%   R-square: 0.9999
%   Adjusted R-square: 0.9999
%   RMSE: 0.004887


figure('name',date)
fs = 15;
set(gcf,'color','white')
plot(wavelength,CD_LONGN_fitted,'k');
hold on
plot(wavelength,CD_LONGS_fitted,'r');
plot(wavelength,CD_CRANB_fitted,'g');
plot(wavelength,CD_ONTOS_fitted,'c');
plot(wavelength,CD_IBAYN_fitted,'--k');
title('CDOM Absorption Coefficient; fitted','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGN','LONGS','CRANB','ONTOS','IBAYN');
xlim([400 900])
grid on

% figure('name',date)
% fs = 15;
% set(gcf,'color','white')
% plot(wavelength,astar_CD_ONTOS ,'k');
% hold on
% plot(wavelength,astar_CD_ONTOS_corrected ,'b');
% plot(wavelength,CD_ONTOS_fitted ,'r');
% % plot([wavelength(1) wavelength(end)],[0 0],'k')
% title('CDOM Absorption Coefficient ONTOS','fontsize',fs)
% xlabel('wavelength [nm]','fontsize',fs)
% ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
% legend('original','corrected','fitted')
% set(gca,'fontsize',fs)
% xlim([400 900])
% grid on
% 
% figure('name',date)
% fs = 15;
% set(gcf,'color','white')
% plot(wavelength,astar_CD_LONGS ,'k');
% hold on
% plot(wavelength,astar_CD_LONGS_corrected ,'b');
% plot(wavelength,CD_LONGS_fitted ,'r');
% % plot([wavelength(1) wavelength(end)],[0 0],'k')
% title('CDOM Absorption Coefficient LONGS','fontsize',fs)
% xlabel('wavelength [nm]','fontsize',fs)
% ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
% legend('original','corrected','fitted')
% set(gca,'fontsize',fs)
% xlim([400 900])
% grid on
% 
% figure('name',date)
% fs = 15;
% set(gcf,'color','white')
% plot(wavelength,astar_CD_LONGN ,'k');
% hold on
% plot(wavelength,astar_CD_LONGN_corrected ,'b');
% plot(wavelength,CD_LONGN_fitted ,'r');
% % plot([wavelength(1) wavelength(end)],[0 0],'k')
% title('CDOM Absorption Coefficient LONGN','fontsize',fs)
% xlabel('wavelength [nm]','fontsize',fs)
% ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
% legend('original','corrected','fitted')
% set(gca,'fontsize',fs)
% xlim([400 900])
% grid on
% 
% figure('name',date)
% fs = 15;
% set(gcf,'color','white')
% plot(wavelength,astar_CD_CRANB ,'k');
% hold on
% plot(wavelength,astar_CD_CRANB_corrected ,'b');
% plot(wavelength,CD_CRANB_fitted ,'r');
% % plot([wavelength(1) wavelength(end)],[0 0],'k')
% title('CDOM Absorption Coefficient CRANB','fontsize',fs)
% xlabel('wavelength [nm]','fontsize',fs)
% ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
% legend('original','corrected','fitted')
% set(gca,'fontsize',fs)
% xlim([400 900])
% grid on
% 
% figure('name',date)
% fs = 15;
% set(gcf,'color','white')
% plot(wavelength,astar_CD_IBAYN ,'k');
% hold on
% plot(wavelength,astar_CD_IBAYN_corrected ,'b');
% plot(wavelength,CD_IBAYN_fitted ,'r');
% % plot([wavelength(1) wavelength(end)],[0 0],'k')
% title('CDOM Absorption Coefficient IBAYN','fontsize',fs)
% xlabel('wavelength [nm]','fontsize',fs)
% ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
% legend('original','corrected','fitted')
% set(gca,'fontsize',fs)
% xlim([400 900])
% grid on

%% CDOM absorption coefficient normalized at a(440)=1
astar_CD_ONTOS440 =  CD_ONTOS_fitted./CD_ONTOS_fitted(wavelength==440);
astar_CD_LONGS440 =  CD_LONGS_fitted./CD_LONGS_fitted(wavelength==440);
astar_CD_LONGN440 =  CD_LONGN_fitted./CD_LONGN_fitted(wavelength==440);
astar_CD_CRANB440 =  CD_CRANB_fitted./CD_CRANB_fitted(wavelength==440);
astar_CD_IBAYN440 =  CD_IBAYN_fitted./CD_IBAYN_fitted(wavelength==440);


disp('a_ONTOS(440):')
disp(CD_ONTOS_fitted(wavelength==440))

disp('a_LONGS(440):')
disp(CD_LONGS_fitted(wavelength==440))

disp('a_LONGN(440):')
disp(CD_LONGN_fitted(wavelength==440))

disp('a_CRANB(440):')
disp(CD_CRANB_fitted(wavelength==440))

disp('a_IBAYN(440):')
disp(CD_IBAYN_fitted(wavelength==440))

% ONTOS a(440)=1
% General model Exp1:
%      f(x) = a*exp(b*x)
% Coefficients (with 95% confidence bounds):
%        a =        1868  (1868, 1868)
%        b =    -0.01712  (-0.01712, -0.01712)
% 
% Goodness of fit:
%   SSE: 2.089e-28
%   R-square: 1
%   Adjusted R-square: 1
%   RMSE: 6.47e-16

% LONGS a(440)=1
% General model Exp1:
%      f(x) = a*exp(b*x)
% Coefficients (with 95% confidence bounds):
%        a =       592.5  (592.5, 592.5)
%        b =    -0.01451  (-0.01451, -0.01451)
% 
% Goodness of fit:
%   SSE: 1.806e-28
%   R-square: 1
%   Adjusted R-square: 1
%   RMSE: 6.016e-16


figure('name',date)
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_CD_LONGN440,'k');
hold on
plot(wavelength,astar_CD_LONGS440,'r');
plot(wavelength,astar_CD_CRANB440,'g');
plot(wavelength,astar_CD_ONTOS440,'c');
plot(wavelength,astar_CD_IBAYN440,'--k');
title('CDOM Absorption Coefficient','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGN','LONGS','CRANB','ONTOS','IBAYN');
grid on

%% saving to be used in HL
% astar = [wavelength(end:-1:1) astar_CD_LONGS440(end:-1:1)];
% save('astar_CDOM_LONGS140929.txt','-ascii','astar')
% 
% astar = [wavelength(end:-1:1) astar_CD_ONTOS440(end:-1:1)];
% save('astar_CDOM_ONTOS140929.txt','-ascii','astar')
