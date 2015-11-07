cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/retrieval
%% Optical Densities
CH230103 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919_2/CH230103.ASC');
CH230251 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919_2/CH230251.ASC');
CH230254 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919_2/CH230254.ASC');
SM230103 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919_2/SM230103.ASC');
SM230251 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919_2/SM230251.ASC');
CH222152 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919_2/CH222152.ASC');
CH226502 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919_2/CH226502.ASC');
CH240052 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919_2/CH240052.ASC');
CH249502 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919_2/CH249502.ASC');
SM222152 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919_2/SM222152.ASC');
SM226502 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919_2/SM226502.ASC');
SM240052 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919_2/SM240052.ASC');
SM249502 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919_2/SM249502.ASC');
CH242042 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919_2/CH242042.ASC');
CH244042 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919_2/CH244042.ASC');
SM242042 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919_2/SM242042.ASC');
SM244042 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919_2/SM244042.ASC');
SMBASE1 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919_2/SMBASE1.ASC');

% CDOM
CDIBAYN  = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919/CDIBAYN.ASC');
CDONTEX  = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919/CDONTEX.ASC');
CDONTNS  = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919/CDONTNS.ASC');
CDONTOS  = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919/CDONTOS.ASC');
CDRVRPIE = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919/CDRVRPIE.ASC');
CDRVRPLM = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919/CDRVRPLM.ASC');
CDBASE1  = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919/CDBASE1.ASC');
CDBRADIN = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919/CDBRADIN.ASC');
CDBRADON = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919/CDBRADON.ASC');
CDCRANB  = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919/CDCRANB.ASC');
CDLONGN  = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919/CDLONGN.ASC');
CDLONGS  = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/130919/CDLONGS.ASC');
%% SM and Chl Absorption Coefficients

[a_SM_LONGN,a_CH_LONGN,wavelength] = ParticleExtraction(SM244042,CH244042,40); % LONGN
[a_SM_LONGS,a_CH_LONGS,~] = ParticleExtraction(SM242042,CH242042,40); % LONGS
[a_SM_CRANB,a_CH_CRANB,~] = ParticleExtraction(SM240052,CH240052,50); % CRANB
[a_SM_ONTNS,a_CH_ONTNS,~] = ParticleExtraction(SM249502,CH249502,500); % ONTNS
[a_SM_IBAYN,a_CH_IBAYN,~] = ParticleExtraction(SM230103,CH230103,100); % ONTNS

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,a_SM_LONGN)
hold on
plot(wavelength,a_SM_LONGS,'k')
plot(wavelength,a_SM_CRANB,'r')
plot(wavelength,a_SM_ONTNS,'g')
plot(wavelength,a_SM_IBAYN,'m')
plot([wavelength(1) wavelength(end)],[0 0],'k')
title('SM Absorption Coefficient','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGN','LONGS','CRANB','ONTNS','IBAYN');
xlim([400 900])
grid on

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,a_CH_LONGN)
hold on
plot(wavelength,a_CH_LONGS,'k')
plot(wavelength,a_CH_CRANB,'r')
plot(wavelength,a_CH_ONTNS,'g')
plot(wavelength,a_CH_IBAYN,'m')
plot([wavelength(1) wavelength(end)],[0 0],'k')
title('CH Absorption Coefficient','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGN','LONGS','CRANB','ONTNS','IBAYN')
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
% X_CH_LONGN = 136.79;%[ug/L] or [mg/m^3]
% X_CH_LONGS = 112.76;
% X_CH_CRANB = 64.08;
% X_CH_ONTNS = .48;
% X_CH_IBAYN = 43.25;

% MONROE County lab Corrected
% X_CH_LONGN = 120;%[ug/L] or [mg/m^3]
% X_CH_LONGS = 104;
% X_CH_CRANB = 72.2;
% X_CH_ONTNS = 3.0;
% X_CH_IBAYN = 48.1;

% Me 10200HC Uncorrected
X_CH_LONGN = 148.47;%[ug/L] or [mg/m^3]
X_CH_LONGS = 124.03;
X_CH_CRANB = 65.77;
X_CH_ONTNS = 1.01;
X_CH_IBAYN = 46.81;

% MONROE County lab Uncorrected
% X_CH_LONGN = 140;%[ug/L] or [mg/m^3]
% X_CH_LONGS = 124.03;
% X_CH_CRANB = 75.1;
% X_CH_ONTNS = 1.72;
% X_CH_IBAYN = 51.7;

astar_CH_LONGN = a_CH_LONGN/X_CH_LONGN;%[m^2/mg]
astar_CH_LONGS = a_CH_LONGS/X_CH_LONGS;
astar_CH_CRANB = a_CH_CRANB/X_CH_CRANB; 
astar_CH_ONTNS = a_CH_ONTNS/X_CH_ONTNS;

figure % in a*, [m^2/mg]
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_CH_LONGN)
hold on
plot(wavelength,astar_CH_LONGS,'k')
plot(wavelength,astar_CH_CRANB,'r')
plot(wavelength,astar_CH_ONTNS,'g')
plot([wavelength(1) wavelength(end)],[0 0],'k')
title('CH mass-specific absorption coeff.','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('a*, [m^2/mg]','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGN','LONGS','CRANB','ONTNS')
xlim([400 900])
grid on

% astar = [wavelength(end:-1:1) astar_CH_CRANB(end:-1:1)];
% save('astarCHL091913.txt','-ascii','astar')


%% Save CHL ONTNS for using in HL

astar_CH_ONTNS_corrected = astar_CH_ONTNS;
astar_CH_ONTNS_corrected(wavelength>=750)=0;
astar_CH_ONTNS_corrected = sgolayfilt(astar_CH_ONTNS_corrected,5,41);
astar_CH_ONTNS_corrected(astar_CH_ONTNS_corrected<0)=0;


figure % in a*, [m^2/mg]
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_CH_ONTNS_corrected,'r')
hold on
% plot(wavelength,astar_CH_ONTNS,'k')
% hold on
plot([wavelength(1) wavelength(end)],[0 0],'k')
title('chlorophyll mass-specific absorption coefficient -- ONTNS','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('a*, [m^2/mg]','fontsize',fs)
set(gca,'fontsize',fs)
% legend('ONTNS')
% xlim([400 900])


% astar = [wavelength(end:-1:1) astar_CH_ONTNS_corrected(end:-1:1)];
% save('astar_CH_ONTNS091913_10200HC.txt','-ascii','astar')
%% Save CHL LONGS for using in HL

astar_CH_LONGS_corrected = astar_CH_LONGS;
astar_CH_LONGS_corrected(wavelength>=750)=0;
astar_CH_LONGS_corrected = sgolayfilt(astar_CH_LONGS_corrected,5,41);
astar_CH_LONGS_corrected(astar_CH_LONGS_corrected<0)=0;


figure % in a*, [m^2/mg]
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_CH_LONGS_corrected,'r','linewidth',1.5)
% hold on
% plot(wavelength,astar_CH_LONGS,'k')
% hold on
% plot([wavelength(1) wavelength(end)],[0 0],'k')
grid on
title('CH mass-specific absorption coeff. -- LONGS','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('a*, [m^2/mg]','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGS')
% xlim([400 900])


% astar = [wavelength(end:-1:1) astar_CH_LONGS_corrected(end:-1:1)];
% save('astar_CH_LONGS091913.txt','-ascii','astar')

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
% X_SM_LONGN = 48;
% X_SM_LONGS = 46;
% X_SM_CRANB = 26.7;
% X_SM_ONTNS = 1.6;
% X_SM_IBAYN = 10.3;

% RIT
X_SM_LONGN = 36.15;
X_SM_LONGS = 23.85;
X_SM_CRANB = 17.33;
X_SM_ONTNS = 0.10;
X_SM_IBAYN = 6.2;

astar_SM_LONGN = a_SM_LONGN/X_SM_LONGN;%[m^2/g]
astar_SM_LONGS = a_SM_LONGS/X_SM_LONGS;
astar_SM_CRANB = a_SM_CRANB/X_SM_CRANB; 
astar_SM_ONTNS = a_SM_ONTNS/X_SM_ONTNS;
astar_SM_IBAYN = a_SM_IBAYN/X_SM_IBAYN;

figure % in a*, [m^2/g]
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_SM_LONGN)
hold on
plot(wavelength,astar_SM_LONGS,'k')
plot(wavelength,astar_SM_CRANB,'r')
plot(wavelength,astar_SM_ONTNS,'g')
plot(wavelength,astar_SM_IBAYN,'m')
plot([wavelength(1) wavelength(end)],[0 0],'k')
title('SM mass-specific absorption coeff. -- RIT','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('a*, [m^2/g]','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGN','LONGS','CRANB','ONTNS','IBAYN')
xlim([400 900])
grid on
%  
% astar = [wavelength(end:-1:1) astar_CH_CRANB(end:-1:1)];
% save('astarCHL091913.txt','-ascii','astar')

%% Save SM ONTNS for using in HL

astar_SM_ONTNS_corrected = astar_SM_ONTNS-astar_SM_ONTNS(wavelength==814);
astar_SM_ONTNS_corrected(wavelength>=814)=0;
astar_SM_ONTNS_corrected = sgolayfilt(astar_SM_ONTNS_corrected,5,41);
astar_SM_ONTNS_corrected(astar_SM_ONTNS_corrected<0)=0;


figure % in a*, [m^2/mg]
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_SM_ONTNS_corrected,'r')
hold on
% plot(wavelength,astar_SM_ONTNS,'k')
% hold on
plot([wavelength(1) wavelength(end)],[0 0],'k')
title('mineral particle mass-specific absorption coefficient -- ONTNS','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('a*, [m^2/g]','fontsize',fs)
set(gca,'fontsize',fs)
% legend('ONTNS')
% xlim([400 900])


% astar = [wavelength(end:-1:1) astar_SM_ONTNS_corrected(end:-1:1)];
% save('astar_SM_ONTNS091913_RIT.txt','-ascii','astar')

%% Save SM LONGS for using in HL

astar_SM_LONGS_corrected = astar_SM_LONGS-astar_SM_LONGS(wavelength==814);
astar_SM_LONGS_corrected(wavelength>=814)=0;
astar_SM_LONGS_corrected = sgolayfilt(astar_SM_LONGS_corrected,5,41);
astar_SM_LONGS_corrected(astar_SM_LONGS_corrected<0)=0;


figure % in a*, [m^2/mg]
fs = 15;
set(gcf,'color','white')
plot(wavelength,astar_SM_LONGS_corrected,'r','linewidth',1.2)
% hold on
% plot(wavelength,astar_SM_LONGS,'k')
% % hold on
% plot([wavelength(1) wavelength(end)],[0 0],'k')
grid on
title('CH mass-specific absorption coeff. -- LONGS','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('a*, [m^2/g]','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGS')
% xlim([400 900])


% astar = [wavelength(end:-1:1) astar_SM_LONGS_corrected(end:-1:1)];
% save('astar_SM_LONGS091913.txt','-ascii','astar')
%% CDOM absorption coefficient

a_CDIBAYN  = CdomExtraction(CDIBAYN);
a_CDONTEX  = CdomExtraction(CDONTEX);
a_CDONTNS  = CdomExtraction(CDONTNS);
a_CDONTOS  = CdomExtraction(CDONTOS);
a_CDRVRPIE = CdomExtraction(CDRVRPIE);
a_CDRVRPLM = CdomExtraction(CDRVRPLM);
a_CDBRADIN = CdomExtraction(CDBRADIN);
a_CDBRADON = CdomExtraction(CDBRADON);
a_CDCRANB  = CdomExtraction(CDCRANB);
a_CDLONGN  = CdomExtraction(CDLONGN);
a_CDLONGS  = CdomExtraction(CDLONGS);


figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,a_CDIBAYN ,'.k');
hold on
plot(wavelength,a_CDONTEX ,'r');
plot(wavelength,a_CDONTNS ,'g');
plot(wavelength,a_CDONTOS ,'c');
plot(wavelength,a_CDRVRPIE,'k');
plot(wavelength,a_CDRVRPLM,'--k');
plot(wavelength,a_CDBRADIN,'--r');
plot(wavelength,a_CDBRADON,'--g');
plot(wavelength,a_CDCRANB ,'--c');
plot(wavelength,a_CDLONGN ,'.-k');
plot(wavelength,a_CDLONGS ,'.-r');

plot([wavelength(1) wavelength(end)],[0 0],'k')
title('CDOM Absorption Coefficient','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('CDIBAYN ','CDONTEX ','CDONTNS ','CDONTOS ','CDRVRPIE','CDRVRPLM','CDBRADIN','CDBRADON','CDCRANB ','CDLONGN ','CDLONGS ');
xlim([400 900])
%% PLOT ONTNS
figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,a_CDONTNS ,'k');
% plot([wavelength(1) wavelength(end)],[0 0],'k')
title('CDOM Absorption Coefficient ONTNS','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
xlim([400 900])
grid on



% ONTNS Chl-a [ug/L] = 0.48; TSS [mg/L] = 1.60

%% exponential fitting
a =       30.29 ; % (28.38, 32.19)
b =     -0.0126 ; % (-0.01275, -0.01246)

a440 = a_CDONTNS(wavelength==440);
a_CDONTNS440 = a_CDONTNS./a440;
hold on
aasterisk = a440*exp(b*(wavelength-440));

y = a*exp(b*wavelength)/0.1185;
figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,y ,'k');
hold on
plot(wavelength,aasterisk,'r')
plot(wavelength,a_CDONTNS440,'b');
plot([wavelength(1) wavelength(end)],[0 0],'k')
title('CDOM Absorption Coefficient ONTNS','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('y','a*','a_{CDONTNS}(440)')
grid on
% xlim([400 900])
%% CDOM absorption coefficient normalized at a(440)=1

a_CDIBAYN440 =  a_CDIBAYN./a_CDIBAYN(wavelength==440);
a_CDONTEX440 =  a_CDONTEX./a_CDONTEX(wavelength==440);
a_CDONTNS440 =  a_CDONTNS./a_CDONTNS(wavelength==440);
a_CDONTOS440 =  a_CDONTOS./a_CDONTOS(wavelength==440);
a_CDRVRPIE440 = a_CDRVRPIE./a_CDRVRPIE(wavelength==440);
a_CDRVRPLM440 = a_CDRVRPLM./a_CDRVRPLM(wavelength==440);
a_CDBRADIN440 = a_CDBRADIN./a_CDBRADIN(wavelength==440);
a_CDBRADON440 = a_CDBRADON./a_CDBRADON(wavelength==440);
a_CDCRANB440 =  a_CDCRANB./a_CDCRANB(wavelength==440);
a_CDLONGN440 =  a_CDLONGN./a_CDLONGN(wavelength==440);
a_CDLONGS440 =  a_CDLONGS./a_CDLONGS(wavelength==440);


disp('a_CDIBAYN(440):')
disp(a_CDIBAYN(wavelength==440))

disp('a_CDONTEX(440):')
disp(a_CDONTEX(wavelength==440))

disp('a_CDONTNS(440):')
disp(a_CDONTNS(wavelength==440))

disp('a_CDONTOS(440):')
disp(a_CDONTOS(wavelength==440))

disp('a_CDRVRPIE(440):')
disp(a_CDRVRPIE(wavelength==440))

disp('a_CDRVRPLM(440):')
disp(a_CDRVRPLM(wavelength==440))

disp('a_CDBRADIN(440):')
disp(a_CDBRADIN(wavelength==440))

disp('a_CDBRADON(440):')
disp(a_CDBRADON(wavelength==440))

disp('a_CDCRANB(440):')
disp(a_CDCRANB(wavelength==440))

disp('a_CDLONGN(440):')
disp(a_CDLONGN(wavelength==440))

disp('a_CDLONGS(440):')
disp(a_CDLONGS(wavelength==440))


figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,a_CDIBAYN440,'k');
hold on
plot(wavelength,a_CDONTEX440,'r');
plot(wavelength,a_CDONTNS440,'g');
plot(wavelength,a_CDONTOS440,'c');
plot(wavelength,a_CDRVRPIE440,'k');
plot(wavelength,a_CDRVRPLM440,'--k');
plot(wavelength,a_CDBRADIN440,'--r');
plot(wavelength,a_CDBRADON440,'--g');
plot(wavelength,a_CDCRANB440,'--c');
plot(wavelength,a_CDLONGN440,'-.k');
plot(wavelength,a_CDLONGS440,'-.r');

plot([wavelength(1) wavelength(end)],[0 0],'k')
title('CDOM Absorption Coefficient','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('CDIBAYN ','CDONTEX ','CDONTNS ','CDONTOS ','CDRVRPIE','CDRVRPLM','CDBRADIN','CDBRADON','CDCRANB ','CDLONGN ','CDLONGS ');
% xlim([400 900])
%% PLOT CDOM LONGS

a_CDLONGS440corrected = a_CDLONGS440-a_CDLONGS440(440);
a_CDLONGS440corrected(wavelength>800)=0;
figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,a_CDLONGS440corrected,'k');
% plot([wavelength(1) wavelength(end)],[0 0],'k')
title('CDOM Absorption Coefficient LONGS','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
xlim([400 900])
grid on


astar = [wavelength(end:-1:1) a_CDLONGS440corrected(end:-1:1)];
save('astar_CDOM_LONGS091913.txt','-ascii','astar')
%%
Rrs = [
3.046E-04;
2.661E-03;
1.044E-02;
1.535E-02;
1.656E-02;
1.532E-02;
1.520E-02;
1.446E-02;
1.349E-02;
7.790E-03;
3.856E-03;
3.328E-03;
2.692E-03;
1.924E-03];
figure
wl = 410:20:670;
plot(wl',pi*Rrs*100)
