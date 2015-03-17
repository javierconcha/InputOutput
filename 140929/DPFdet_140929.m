%% Determination of dpf
addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/retrieval/')
cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929
%%
wavelength = [...
  4.02500E+02  4.07500E+02  4.12500E+02  4.17500E+02  4.22500E+02  4.27500E+02  4.32500E+02  4.37500E+02  4.42500E+02  4.47500E+02 ...
  4.52500E+02  4.57500E+02  4.62500E+02  4.67500E+02  4.72500E+02  4.77500E+02  4.82500E+02  4.87500E+02  4.92500E+02  4.97500E+02 ...
  5.02500E+02  5.07500E+02  5.12500E+02  5.17500E+02  5.22500E+02  5.27500E+02  5.32500E+02  5.37500E+02  5.42500E+02  5.47500E+02 ...
  5.52500E+02  5.57500E+02  5.62500E+02  5.67500E+02  5.72500E+02  5.77500E+02  5.82500E+02  5.87500E+02  5.92500E+02  5.97500E+02 ...
  6.02500E+02  6.07500E+02  6.12500E+02  6.17500E+02  6.22500E+02  6.27500E+02  6.32500E+02  6.37500E+02  6.42500E+02  6.47500E+02 ...
  6.52500E+02  6.57500E+02  6.62500E+02  6.67500E+02  6.72500E+02  6.77500E+02  6.82500E+02  6.87500E+02  6.92500E+02  6.97500E+02 ...
  7.02500E+02  7.07500E+02  7.12500E+02  7.17500E+02  7.22500E+02  7.27500E+02  7.32500E+02  7.37500E+02  7.42500E+02  7.47500E+02 ...
  7.52500E+02  7.57500E+02  7.62500E+02  7.67500E+02  7.72500E+02  7.77500E+02  7.82500E+02  7.87500E+02  7.92500E+02  7.97500E+02 ...
  8.02500E+02  8.07500E+02  8.12500E+02  8.17500E+02  8.22500E+02  8.27500E+02  8.32500E+02  8.37500E+02  8.42500E+02  8.47500E+02 ...
  8.52500E+02  8.57500E+02  8.62500E+02  8.67500E+02  8.72500E+02  8.77500E+02  8.82500E+02  8.87500E+02  8.92500E+02  8.97500E+02 ...
  9.02500E+02  9.07500E+02  9.12500E+02  9.17500E+02  9.22500E+02  9.27500E+02  9.32500E+02  9.37500E+02  9.42500E+02  9.47500E+02 ...
  9.52500E+02  9.57500E+02  9.62500E+02  9.67500E+02  9.72500E+02  9.77500E+02  9.82500E+02  9.87500E+02  9.92500E+02  9.97500E+02];

wavelength = wavelength'*0.001;

L8bands = [0.4430,0.4826,0.5613,0.6546,0.8646,1.6090,2.2010];



%% LUT from HL for same concentrations but different DPFs
rr = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/RvectorDPFdet_140929_150305.txt');

nruns = size(rr,1)/size(wavelength,1);
LUT_DPFdet_ = reshape(rr(:,1),size(wavelength,1),nruns);
% LUT_DPFdet_ = LUT_DPFdet_*pi;

LUT_DPFdet_L8 = spect_sampL8(LUT_DPFdet_,wavelength);

% Concentrations

filename = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/concentration_listDPFdet140929_150306.txt';

fid = fopen(filename);
c1 = textscan(fid,'%s %s %s');
fclose all;

c = {[c1{1}] [c1{2}] [c1{3}]};
LUTconc = [c{1}(:) c{2}(:) c{3}(:)];

rule0 = strcmp(c{2}(:),c{3}(:));
rule1 = strcmp(c{1}(:),'input140929ONTOS')&strcmp(c{2}(:),c{3}(:));

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,LUT_DPFdet_L8(rule1,:))
title('DPF det. -- LUT','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('R_{rs} [1/sr]','fontsize',fs)
set(gca,'fontsize',fs)
%% ONTOS
pp = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/ONTOSL8_Ref_140919corr.txt');
wavelengthSVC = pp(:,1);
ONTOS_SVC = pp(:,2);
ONTOS_L8 = spect_sampL8(ONTOS_SVC,wavelengthSVC);

% ONTNS = [ 0.003355  0.004671  0.004345  0.000874 ... % from 140919
%           0.000042 -0.000003  0.000006 ]; % in Rrs



LUTused = LUT_DPFdet_L8(rule1,1:5);
LUTconcused = LUTconc(rule1,:);

[Y,I1] = min(sqrt(mean((LUTused(:,1:5)-ones(size(LUTused,1),1)*ONTOS_L8(:,1:5)).^2,2)));


disp(LUTconcused(I1,:))

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,ONTOS_L8,'.-k','linewidth',1.5)
hold on
plot(L8bands(1:5),LUTused(I1,:),'.-r','linewidth',1.5)
plot(L8bands,LUT_DPFdet_L8(rule1,:))
title('DPF det. -- ONTOS','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('R_{rs} [1/sr]','fontsize',fs)
set(gca,'fontsize',fs)
legend('ONTOS field',char(LUTconcused(I1,:)))
grid on

Ref = [L8bands', [LUTused(I1,:)';0;0]];
save('ONTOSL8_Ref_140919_HL.txt','Ref','-ascii')
%% LONGS
pp = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/LONGSRef_140919.txt');
wavelengthSVC = pp(:,1);
LONGS_SVC = pp(:,2);
LONGS_L8 = spect_sampL8(LONGS_SVC,wavelengthSVC);

figure
plot(L8bands,LONGS_L8)


rule2 = strcmp(c{1}(:),'input140929LONGS')&strcmp(c{2}(:),c{3}(:));

LUTused = LUT_DPFdet_L8(rule2,1:5);
LUTconcused = LUTconc(rule2,:);

[Y,I1] = min(sqrt(mean((LUTused(:,1:5)-ones(size(LUTused,1),1)*LONGS_L8(:,1:5)).^2,2)));


disp(LUTconcused(I1,:))

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,LONGS_L8,'.-k')
hold on
plot(L8bands(1:5),LUTused(I1,:),'.-r')
title('DPF det. -- LONGS','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('R_{rs} [1/sr]','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGS field',char(LUTconcused(I1,:)))

%% Cranb
rr = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/RvectorCranb.txt');


nruns = size(rr,1)/size(wavelength,1);
RrsCranb = reshape(rr(:,1),size(wavelength,1),nruns);
% RrsCranb = RrsCranb*pi;

RrsCranbL8 = spect_sampL8(RrsCranb,wavelength);

% figure
% fs = 15;
% set(gcf,'color','white')
% plot(L8bands,RrsCranbL8)
% % hold on
% % plot(wavelength,ONTNSRefinterp,'.-r')
% title('Rrs for Cranb','fontsize',fs)
% xlabel('wavelength [nm]','fontsize',fs)
% ylabel('reflectance','fontsize',fs)
% set(gca,'fontsize',fs)
% grid on
% xlim([400 2200])
% ylim([0 .3])

% Cranb = [0.010979 0.017012 0.044860 ...
%     0.025397 0.007961 0.001144 0.000468];
Cranb = [ 0.003640 0.005506 0.014292 0.008114 ...
           0.002373 0.000267 0.000190]; % in Rrs
% find Cranb in waterpixels with index I
[Y,I2] = min(sqrt(mean((RrsCranbL8(:,1:5)-ones(size(RrsCranbL8,1),1)*Cranb(:,1:5)).^2,2)));

C2 = textscan(c{1}{I2},'%s');
disp(C2{:})

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,Cranb,'.-k')
hold on
plot(L8bands,RrsCranbL8(I2,:),'.-r')
title('DPF det. -- Cranb','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
legend('Cranb field',char(C2{1}))

%% ONTNS
rr = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/RvectorONTNS.txt');


nruns = size(rr,1)/size(wavelength,1);
RrsONTNS = reshape(rr(:,1),size(wavelength,1),nruns);
% RrsONTNS = RrsONTNS*pi;

RrsONTNSL8LUT = spect_sampL8(RrsONTNS,wavelength);

% figure
% fs = 15;
% set(gcf,'color','white')
% plot(L8bands,RrsONTNSL8LUT)
% % hold on
% % plot(wavelength,ONTNSRefinterp,'.-r')
% title('Rrs for ONTNS','fontsize',fs)
% xlabel('wavelength [nm]','fontsize',fs)
% ylabel('reflectance','fontsize',fs)
% set(gca,'fontsize',fs)
% grid on
% % xlim([400 2200])
% % ylim([0 .3])

% ONTNS = [0.012930 0.019153 0.021296 0.004983 ...
%          0.000745 0.000520 0.000194]; %  before in reflectance

ONTNS = [ 0.003355  0.004671  0.004345  0.000874 ...
          0.000042 -0.000003  0.000006 ]; % in Rrs

% ONTNS = RrsONTNSL8corr;% for ELM from SVCextract130919.m

% find ONTNS in waterpixels with index I
[~,I3] = min(sqrt(mean((RrsONTNSL8LUT(:,1:5)-ones(size(RrsONTNSL8LUT,1),1)*ONTNS(:,1:5)).^2,2)));

C3 = textscan(c{1}{I3},'%s');
disp(C3{:})

RrsONTNSL8HL = RrsONTNSL8LUT(I3,:);
RrsONTNSL8HL(6:7) = 0; % because HL is not > 1000nm, so SWIR bands are NaN

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,ONTNS,'.-k')
hold on
plot(L8bands,RrsONTNSL8HL,'.-r')
title('DPF det. -- ONTNS','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
legend('ONTNS field',char(C3{1}))
%% Save to use in ELM
NSRef = [L8bands', RrsONTNSL8HL'];
save([pathname,pathdate,'RrsONTNSL8HL.txt'],'NSRef','-ascii')

%% Plot all
figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,LongS,'-k')
hold on
plot(L8bands,LUT_DPFdet_L8(I1,:),'.-k')
plot(L8bands,Cranb,'-r')
plot(L8bands,RrsCranbL8(I2,:),'.-r')
plot(L8bands,ONTNS,'-b')
plot(L8bands,RrsONTNSL8LUT(I3,:),'.-b')
title('DPF det. -- LongS','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('rem-sens reflectance R_{rs} (sr^{-1})','fontsize',fs)
set(gca,'fontsize',fs)
legend('LONGS field',char(C1{1}),'Cranb field',char(C2{1}),'ONTNS field',char(C3{1}))
%%
figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,LUT_DPFdet_L8(I1-1:I1+1,:),'k')
hold on
plot(L8bands,RrsCranbL8(I2-1:I2+1,:),'r')
plot(L8bands,RrsONTNSL8(I3-1:I3+1,:),'b')
title('Rrs for LONGS','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
grid on
legend('LONGS','Cranb','ONTNS')
