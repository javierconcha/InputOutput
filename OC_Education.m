addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/retrieval/')
addpath('/Users/javier/Downloads/mtit')
cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/
%% Landsat 8 image in Rrs
folderpath = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/';
filename = 'LC80160302013262LGN00/MOBELM/LC80160302013262LGN00_ONelm140629.tif';

filepath = [folderpath filename];

[RGB,map] = imread(filepath);

% RGB_R = double(RGB(:,:,4));
% RGB_R(RGB_R<0)=0;
% 
% RGB_G = double(RGB(:,:,3));
% RGB_G(RGB_G<0)=0;
% 
% RGB_B = double(RGB(:,:,2));
% RGB_B(RGB_B<0)=0;

% LUTs from HydroLight

clear c c1 Rrs LUT LUTconc LUTconcDPF LUTused InputType% if other retrieval's variables are in Workspace

% % new 03/02/15, spectrally sampling made in matlab
% LUTfilename1 = 'Rvector140929_150305_2.txt';
% LUTfilename1 = 'Rvector140929_150317.txt'; % with FFbb018.dpf for ONTOS
% LUTfilename1 = 'Rvector140929_150407.txt'; % with FFbb012.dpf for ONTOS
% LUTfilename1 = 'Rvector140929_150409.txt'; % with FFbb012.dpf for ONTOS
LUTfilename1 = 'Rvector140929_150420.txt'; % with more dpfs


% LUTconcfilename1 = 'concentration_list140929_150305_2.txt';
% LUTconcfilename1 = 'concentration_list140929_150317.txt';
% LUTconcfilename1 = 'concentration_list140929_150406.txt';
% LUTconcfilename1 = 'concentration_list140929_150409.txt';
LUTconcfilename1 = 'concentration_list140929_150420.txt';

filepath = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/';

LUTpath1 = [filepath LUTfilename1];
rr1 = load(LUTpath1); % Created for 09/29/14 image!!!

LUTconpath1 = [filepath LUTconcfilename1];
fid = fopen(LUTconpath1);
c1 = textscan(fid,'%s %f %f %f %s');
fclose all;



c = {[c1{1}] [c1{2}] [c1{3}] [c1{4}] [c1{5}]};


LUTconc = [c{2}(:) c{3}(:) c{4}(:)];

% LUT from HL with 120 wavelength

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

wavelength = wavelength'*0.001; % um

L8bands = [0.4430,0.4826,0.5613,0.6546,0.8646,1.6090,2.2010];

nruns = size(rr1,1)/size(wavelength,1);
Rrs = reshape(rr1(:,1),size(wavelength,1),nruns);



LUT = spect_sampL8(Rrs,wavelength);
LUTconcDPF = nan(size(LUT,1),1);
for index = 1:size(LUT,1) 
    if strcmp(c1{5}(index),'FFbb005.dpf')
        LUTconcDPF(index)= 0.5;
    elseif strcmp(c1{5}(index),'FFbb006.dpf')
        LUTconcDPF(index)= 0.6; 
    elseif strcmp(c1{5}(index),'FFbb007.dpf')
        LUTconcDPF(index)= 0.7;
    elseif strcmp(c1{5}(index),'FFbb008.dpf')
        LUTconcDPF(index)= 0.8;
    elseif strcmp(c1{5}(index),'FFbb009.dpf')
        LUTconcDPF(index)= 0.9;
    elseif strcmp(c1{5}(index),'FFbb010.dpf')
        LUTconcDPF(index)= 1.0; 
    elseif strcmp(c1{5}(index),'FFbb012.dpf')
        LUTconcDPF(index)= 1.2;
    elseif strcmp(c1{5}(index),'FFbb014.dpf')
        LUTconcDPF(index)= 1.4;
    elseif strcmp(c1{5}(index),'FFbb016.dpf')
        LUTconcDPF(index)= 1.6;
    elseif strcmp(c1{5}(index),'FFbb018.dpf')
        LUTconcDPF(index)= 1.8;  
    elseif strcmp(c1{5}(index),'FFbb020.dpf')
        LUTconcDPF(index)= 2.0;     
    elseif strcmp(c1{5}(index),'FFbb022.dpf')
        LUTconcDPF(index)= 2.2;      
    elseif strcmp(c1{5}(index),'FFbb024.dpf')
        LUTconcDPF(index)= 2.4;  
    elseif strcmp(c1{5}(index),'FFbb026.dpf')
        LUTconcDPF(index)= 2.6;
    end
end
%% Fixed CDOM
tilex = 200;
tiley = 50;
input = 'input140929LONGS';
DPFconc = 0.5;
CDconc = 0.0954;
CTE = 2.0; % constant for the color adjustment


rule = strcmp(c{1}(:),input) & LUTconcDPF(:)==DPFconc;
LUTused = LUT(rule,:);

% s0101
CHconc = 0.01; SMconc = 0.01;  

[M_s0101_RGB,Rrs_s0101] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0102
CHconc = 0.01; SMconc = 9.11;

[M_s0102_RGB,Rrs_s0102] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0103
CHconc = 0.01; SMconc = 30.7;

[M_s0103_RGB,Rrs_s0103] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0104
CHconc = 0.01; SMconc = 50.0;

[M_s0104_RGB,Rrs_s0104] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0201
CHconc = 25.0; SMconc = 0.01;

[M_s0201_RGB,Rrs_s0201] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0202
CHconc = 25.0; SMconc = 9.11;

[M_s0202_RGB,Rrs_s0202] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0203
CHconc = 25.0; SMconc = 30.7;

[M_s0203_RGB,Rrs_s0203] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0204
CHconc = 25.0; SMconc = 50.0;

[M_s0204_RGB,Rrs_s0204] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0301
CHconc = 60.0; SMconc = 0.01;

[M_s0301_RGB,Rrs_s0301] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0302
CHconc = 60.0; SMconc = 9.11;

[M_s0302_RGB,Rrs_s0302] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0303
CHconc = 60.0; SMconc = 30.7;

[M_s0303_RGB,Rrs_s0303] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0304
CHconc = 60.0; SMconc = 50.0;

[M_s0304_RGB,Rrs_s0304] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
%
M = [ ...
  M_s0101_RGB,M_s0102_RGB,M_s0103_RGB,M_s0104_RGB; ...
  M_s0201_RGB,M_s0202_RGB,M_s0203_RGB,M_s0204_RGB; ...
  M_s0301_RGB,M_s0302_RGB,M_s0303_RGB,M_s0304_RGB];

M_R = M(:,:,1);
M_G = M(:,:,2);
M_B = M(:,:,3);  

RGB_R = LUTused(:,4);
RGB_G = LUTused(:,3);
RGB_B = LUTused(:,2);

% Adjusting threshold for display


threshold_Bpos = mean(RGB_B(:))+CTE*std(RGB_B(:));
M_B(M_B>threshold_Bpos)=threshold_Bpos;
threshold_Bneg = mean(RGB_B(:))-CTE*std(RGB_B(:));
M_B(M_B<threshold_Bneg)=threshold_Bneg;
M_B(M_B<0)=0;

threshold_Gpos = mean(RGB_G(:))+CTE*std(RGB_G(:));
M_G(M_G>threshold_Gpos)=threshold_Gpos;
threshold_Gneg = mean(RGB_G(:))-CTE*std(RGB_G(:));
M_G(M_G<threshold_Gneg)=threshold_Gneg;
M_G(M_G<0)=0;

threshold_Rpos = mean(RGB_R(:))+CTE*std(RGB_R(:));
M_R(M_R>threshold_Rpos)=threshold_Rpos;
threshold_Rneg = mean(RGB_R(:))-CTE*std(RGB_R(:));
M_R(M_R<threshold_Rneg)=threshold_Rneg;
M_R(M_R<0)=0;

% Convert values to [0 1] for display
RGBdisplay_B = (M_B)/(max(RGB_B(:))-min(RGB_B(:)));
RGBdisplay_B(RGBdisplay_B<0)=0;
RGBdisplay_G = (M_G)/(max(RGB_G(:))-min(RGB_G(:)));
RGBdisplay_G(RGBdisplay_G<0)=0;
RGBdisplay_R = (M_R)/(max(RGB_R(:))-min(RGB_R(:)));
RGBdisplay_R(RGBdisplay_R<0)=0;

RGBdisplay(:,:,1) = RGBdisplay_R;
RGBdisplay(:,:,2) = RGBdisplay_G;
RGBdisplay(:,:,3) = RGBdisplay_B;

% figure
% imshow(RGBdisplay)

figure
image(RGBdisplay)

% figure
% imagesc(RGBdisplay)
%%
figure
plot(L8bands,Rrs_s0101)
hold on
plot(L8bands,Rrs_s0201)
plot(L8bands,Rrs_s0301)
plot(L8bands,Rrs_s0102)
plot(L8bands,Rrs_s0202)
plot(L8bands,Rrs_s0302)

%% FIXED SM
tilex = 200;
tiley = 50;
input = 'input140929LONGS';
DPFconc = 1.0;
SMconc =1.4;
CTE = 2.0; % constant for the color adjustment


rule = strcmp(c{1}(:),input) & LUTconcDPF(:)==DPFconc;
LUTused = LUT(rule,:);

% s0101
CHconc = 0.01; CDconc = 0.0954;  

[M_s0101_RGB,Rrs_s0101] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0102
CHconc = 0.01; CDconc = 0.9297;

[M_s0102_RGB,Rrs_s0102] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0103
CHconc = 0.01; CDconc = 1.0025;

[M_s0103_RGB,Rrs_s0103] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0104
CHconc = 0.01; CDconc = 1.0194;

[M_s0104_RGB,Rrs_s0104] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0201
CHconc = 25.0; CDconc = 0.0954;

[M_s0201_RGB,Rrs_s0201] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0202
CHconc = 25.0; CDconc = 0.9297;

[M_s0202_RGB,Rrs_s0202] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0203
CHconc = 25.0; CDconc = 1.0025;

[M_s0203_RGB,Rrs_s0203] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0204
CHconc = 25.0; CDconc = 1.0194;

[M_s0204_RGB,Rrs_s0204] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0301
CHconc = 60.0; CDconc = 0.0954;

[M_s0301_RGB,Rrs_s0301] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0302
CHconc = 60.0; CDconc = 0.9297;

[M_s0302_RGB,Rrs_s0302] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0303
CHconc = 60.0; CDconc = 1.0025;

[M_s0303_RGB,Rrs_s0303] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
% s0304
CHconc = 60.0; CDconc = 1.0194;

[M_s0304_RGB,Rrs_s0304] = createtile(input,CHconc,SMconc,CDconc,DPFconc,c,LUTconcDPF,LUTconc,LUT,tilex,tiley);
%
M = [ ...
  M_s0101_RGB,M_s0102_RGB,M_s0103_RGB,M_s0104_RGB; ...
  M_s0201_RGB,M_s0202_RGB,M_s0203_RGB,M_s0204_RGB; ...
  M_s0301_RGB,M_s0302_RGB,M_s0303_RGB,M_s0304_RGB];

M_R = M(:,:,1);
M_G = M(:,:,2);
M_B = M(:,:,3);  

RGB_R = LUTused(:,4);
RGB_G = LUTused(:,3);
RGB_B = LUTused(:,2);

% Adjusting threshold for display


threshold_Bpos = mean(RGB_B(:))+CTE*std(RGB_B(:));
M_B(M_B>threshold_Bpos)=threshold_Bpos;
threshold_Bneg = mean(RGB_B(:))-CTE*std(RGB_B(:));
M_B(M_B<threshold_Bneg)=threshold_Bneg;
M_B(M_B<0)=0;

threshold_Gpos = mean(RGB_G(:))+CTE*std(RGB_G(:));
M_G(M_G>threshold_Gpos)=threshold_Gpos;
threshold_Gneg = mean(RGB_G(:))-CTE*std(RGB_G(:));
M_G(M_G<threshold_Gneg)=threshold_Gneg;
M_G(M_G<0)=0;

threshold_Rpos = mean(RGB_R(:))+CTE*std(RGB_R(:));
M_R(M_R>threshold_Rpos)=threshold_Rpos;
threshold_Rneg = mean(RGB_R(:))-CTE*std(RGB_R(:));
M_R(M_R<threshold_Rneg)=threshold_Rneg;
M_R(M_R<0)=0;

% Convert values to [0 1] for display
RGBdisplay_B = (M_B)/(max(RGB_B(:))-min(RGB_B(:)));
RGBdisplay_B(RGBdisplay_B<0)=0;
RGBdisplay_G = (M_G)/(max(RGB_G(:))-min(RGB_G(:)));
RGBdisplay_G(RGBdisplay_G<0)=0;
RGBdisplay_R = (M_R)/(max(RGB_R(:))-min(RGB_R(:)));
RGBdisplay_R(RGBdisplay_R<0)=0;

RGBdisplay(:,:,1) = RGBdisplay_R;
RGBdisplay(:,:,2) = RGBdisplay_G;
RGBdisplay(:,:,3) = RGBdisplay_B;

% figure
% imshow(RGBdisplay)

figure
image(RGBdisplay)

% figure
% imagesc(RGBdisplay)
%%
figure
plot(L8bands,Rrs_s0101)
hold on
plot(L8bands,Rrs_s0201)
plot(L8bands,Rrs_s0301)
plot(L8bands,Rrs_s0102)
plot(L8bands,Rrs_s0202)
plot(L8bands,Rrs_s0302)

%% Curves with fixed CH
input = 'input140929ONTOS';
DPFconc = 1.0;
CHconc = 0.01;
CHconc = 25.0;
CHconc = 60.0;
rule = strcmp(c{1}(:),input) & LUTconcDPF(:)==DPFconc & ...
    LUTconc(:,1)==CHconc;
Rrsused = Rrs(:,rule);

figure
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(wavelength,Rrsused')
str1 = sprintf('C_a = %2.2f [ug/L]',CHconc);
title(str1,'fontsize',fs)
xlabel('wavelength [um]')
ylabel('remote-sensing reflectance [1/sr]')
%% Standard Algorithms, no fixed

% input = 'input140929LONGS';
% DPFconc = 0.5;


% rule = strcmp(c{1}(:),input) & LUTconcDPF(:)==DPFconc;
rule = ~isnan(LUTconc(:,1));

LUTused = LUT(rule,:);

for index=1:size(LUTused,1)
    if LUTused(:,1)>LUTused(:,2)
        RdivG = LUTused(:,1)./LUTused(:,3);
    else
        RdivG = LUTused(:,2)./LUTused(:,3);
    end

end

% Standard algorithms
R = 0.1:0.01:15;
Rlog10 = log10(R);

% OC2v4, O'Reilly et al. (2000)
a = [0.319,-2.336,0.879,-0.135,-0.071];
chl_oc2 = 10.0.^(a(1)+a(2)*Rlog10+a(3)*Rlog10.^2+a(4)*Rlog10.^3) + a(5);

% OC3
a = [0.283,-2.753,1.457,0.659,-1.403];
chl_oc3 = 10.0.^(a(1)+a(2).*Rlog10+a(3)*Rlog10.^2+a(4)*Rlog10.^3+a(5)*Rlog10.^4);

% OC4v4, O'Reilly et al. (2000)
a = [0.366,-3.067,1.930,0.649,-1.532];
chl_oc4 = 10.0.^(a(1)+a(2).*Rlog10+a(3)*Rlog10.^2+a(4)*Rlog10.^3+a(5)*Rlog10.^4);

% CI, Hu et al. (2012)

CI = LUTused(:,3)-(LUTused(:,1)+((561-443)/(655-443))*(LUTused(:,4)-LUTused(:,1)));

figure
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
loglog(RdivG,LUTconc(rule,1),'.')
hold on
loglog(R,chl_oc3,'Color',[0 0.5 0] ,'LineWidth',2.0)
loglog(R,chl_oc2,'Color',[0.5 0 0] ,'LineWidth',2.0)
loglog(R,chl_oc4,'Color',[0 0 0] ,'LineWidth',2.0)
legend('From Hydrolight','OC3 Model','OC2v4 Model','OC4v4 Model')
grid on
% str1 = sprintf('C_a = %2.2f [mg/L]',CHconc);
% title(str1,'fontsize',fs)
xlabel('R_{blue}/R_{green}')
ylabel('C_a [ug/L]')
xlim([0.1 12])
%% Standard Algorithms

input = 'input140929ONTOS';
DPFconc = 0.5;
SMconc = 0.0100;
CDconc = 0.0954;


rule = strcmp(c{1}(:),input) & LUTconcDPF(:)==DPFconc & ...
    LUTconc(:,2)==SMconc & LUTconc(:,3)==CDconc;
% rule = ~isnan(LUTconc(:,1));

LUTused = LUT(rule,:);

for index=1:size(LUTused,1)
    if LUTused(:,1)>LUTused(:,2)
        RdivG = LUTused(:,1)./LUTused(:,3);
    else
        RdivG = LUTused(:,2)./LUTused(:,3);
    end

end

figure
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
loglog(RdivG,LUTconc(rule,1),'.')
hold on
loglog(R,chl_oc3,'Color',[0 0.5 0] ,'LineWidth',2.0)
loglog(R,chl_oc2,'Color',[0.5 0 0] ,'LineWidth',2.0)
loglog(R,chl_oc4,'Color',[0 0 0] ,'LineWidth',2.0)
legend('Hydrolight','OC3 Model','OC2v4 Model','OC4v4 Model')
grid on
% str1 = sprintf('C_a = %2.2f [mg/L]',CHconc);
% title(str1,'fontsize',fs)
xlabel('R_{blue}/R_{green}')
ylabel('C_a [mg/m^3]')
xlim([0.1 12])