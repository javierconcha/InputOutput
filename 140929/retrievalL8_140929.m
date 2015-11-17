% Retrieve the concentration from a image using a LUT from HL for L8 image of 140929.
% Version 3.0
% Created by Javier A. Concha
% 03/02/15
addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/retrieval/')
addpath('/Users/javier/Downloads/mtit')
cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929
%% L8 image cropped
% From ELM using new L8 reflectance product
% folderpath = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/LC80170302014272LGN00/LC80170302014272LGN00_ROI_Rrs_150408.tif';
folderpath = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/LC80170302014272LGN00/MOBELM/';
filename = 'LC80170302014272LGN00_ROI_Rrs_150418CRANB_ONTOS.tif';
date = '140929';

filepath = [folderpath filename];
clear imL8crop imL8cropRGB maskRGB; % if other retrieval's variables are in Workspace

[imL8crop, cmap] = imread(filepath);
INFO = imfinfo(filepath);

%%%% Mask
imL8cropmask = imread(...
    '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/LC80170302014272LGN00/MOBELM/LC80170302014272LGN00_ROImask.tif');

imL8cropmask(imL8cropmask>0)=1;


L8bands = [0.4430,0.4826,0.5613,0.6546,0.8646,1.6090,2.2010];

%%%% water pixels. Convert each band in columns.

imnew = reshape(imL8crop,[size(imL8crop,1)*size(imL8crop,2) size(imL8crop,3)]);
masknew = reshape(imL8cropmask,[size(imL8cropmask,1)*size(imL8cropmask,2) size(imL8cropmask,3)]);



waterpixels = imnew(masknew==1,:);
waterpixels = double(waterpixels);

waterpixels(isnan(waterpixels))=0;

% added 01-11-14. plot radiance curves
% radiance image
imL8radcrop = imread(...
    '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images//LC80170302014272LGN00/MOBELM/LC80170302014272LGN00_ROItif.tif');
imradnew = reshape(imL8radcrop,[size(imL8radcrop,1)*size(imL8radcrop,2) size(imL8radcrop,3)]);

waterradpixels = imradnew(masknew==1,:);
waterradpixels = double(waterradpixels);

% for displaying

imL8cropRGB(:,:,1)=imadjust(imL8crop(:,:,4));
imL8cropRGB(:,:,2)=imadjust(imL8crop(:,:,3));
imL8cropRGB(:,:,3)=imadjust(imL8crop(:,:,2));


impos = double(imL8cropRGB);
% impos(impos<0)=0;% only positive values

maskRGB(:,:,1)=double(imL8cropmask);
maskRGB(:,:,2)=double(imL8cropmask);
maskRGB(:,:,3)=double(imL8cropmask);
 

impos = impos.*maskRGB;

%% RGB display
figure
set(gcf,'color','white')
imagesc(impos)
axis equal

%% mask display
figure
set(gcf,'color','white')
imshow(imadjust(imL8cropmask))



%% Stats water pixels
meanwp = mean(waterpixels,1);
stdwp = std(waterpixels,1);
maxwp = max(waterpixels,[],1);
minwp = min(waterpixels,[],1);


figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,meanwp,'k')
title('Reflectance water L8 image','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
hold on
plot(L8bands,meanwp+stdwp,'g')
plot(L8bands,meanwp-stdwp,'g')
plot(L8bands,maxwp,'r')
plot(L8bands,minwp,'r')
xlim([min(L8bands) max(L8bands)])

format short
disp('---------------------------------------------------')
disp('Basic Stats      Min       Max       Mean     Stdev  ') 
disp('---------------------------------------------------')

bands = [1 2 3 4 5 6 7];

for i = 1:size(minwp,2)
    str = sprintf('     band %i  %2.6f  %2.6f  %2.6f  %2.6f  %2.6f  %2.6f',bands(i),minwp(i), maxwp(i), meanwp(i), stdwp(i));
    disp(str)
end

nbins = 100;

figure
subplot(2,4,1)
fs = 15;
set(gcf,'color','white')
hist(waterpixels(:,1),nbins)
title('band 1','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,4,2)
fs = 15;
set(gcf,'color','white')
hist(waterpixels(:,2),nbins)
title('band 2','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,4,3)
fs = 15;
set(gcf,'color','white')
hist(waterpixels(:,3),nbins)
title('band 3','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,4,4)
fs = 15;
set(gcf,'color','white')
hist(waterpixels(:,4),nbins)
title('band 4','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,4,5)
fs = 15;
set(gcf,'color','white')
hist(waterpixels(:,5),nbins)
title('band 5','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,4,6)
fs = 15;
set(gcf,'color','white')
hist(waterpixels(:,6),nbins)
title('band 6','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,4,7)
fs = 15;
set(gcf,'color','white')
hist(waterpixels(:,7),nbins)
title('band 7','fontsize',fs)
set(gca,'fontsize',fs)

p=mtit('Histogram Water Pixels per Band',...
 	     'fontsize',fs+1,'xoff',0,'yoff',.025);
     
%% display All water pixels

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,waterpixels')
title('Reflectance water L8 image','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
% ylim([0 0.18])
 hold on
 plot(L8bands,meanwp,'g','linewidth',2)
 grid on
 
%% display All water pixels RADIANCE values

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,waterradpixels')
title('Radiance water L8 image','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('Radiance [W/m^2/sr]','fontsize',fs)
set(gca,'fontsize',fs)
% ylim([0 0.18])
%  hold on
%  plot(L8bands,meanwp+stdwp,'g','linewidth',2) 

%% negative values
im = double(imL8crop);
imneg = zeros(size(im));
imneg(im<0)=im(im<0);% only negatives

imnegmask = zeros(size(imL8crop));% for displaying
imnegmask(im<0)=1; % negative values are white
imnegmask(im>=0)=0; % positive values are black
imnegmask = imnegmask+0.5*repmat(~double(imL8cropmask),[1 1 size(imL8crop,3)]); % for the land appear gray

figure
subplot(2,2,1)
fs = 15;
set(gcf,'color','white')
imshow(imnegmask(:,:,5))
title('band 5','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,2,2)
fs = 15;
set(gcf,'color','white')
imshow(imnegmask(:,:,6))
title('band 6','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,2,3)
fs = 15;
set(gcf,'color','white')
imshow(imnegmask(:,:,7))
title('band 7','fontsize',fs)
set(gca,'fontsize',fs)


p=mtit('Negative Values',...
 	     'fontsize',fs+1,'xoff',0,'yoff',.025);
     
     
format short
disp('----------------------------------------------------------')
disp('Basic Stats      Min       Max       Mean     Stdev    N  ') 
disp('----------------------------------------------------------')

bands = [5 6 7];

for i = 1:size(bands,2)
    data = imneg(:,:,bands(i));
    str = sprintf('     band %i  %2.6f  %2.6f  %2.6f  %2.6f  %2.6f  %2.6f %6.0i'...
        ,bands(i),min(data(:)), max(data(:)), mean(data(:)), std(data(:)),sum(sum(data<0)));
    disp(str)
end
     
%% Display negative values
bn = 7;
waterpixels_neg = waterpixels(waterpixels(:,bn)<0,:);

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,waterpixels_neg')
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
str = sprintf('Curves with negatives values in Band %i',bn);
title(str,'fontsize',fs)
% ylim([0 0.18])

%% Display high values in the different bands
im = double(imL8crop);
bn = 5;
imhighNIR = zeros(size(im,1),size(im,2));
cond2 = im(:,:,bn)> (meanwp(bn)+2*stdwp(bn)) & imL8cropmask ~= 0;
imhighNIR(cond2)=1;% only high NIR

imhighNIRmask = imhighNIR+0.5*repmat(~double(imL8cropmask),[1 1 1]); % for the land appear gray


figure
fs = 15;
set(gcf,'color','white')
imshow(imhighNIRmask)
str = sprintf('High values for band %i',bn);
title(str,'fontsize',fs)
set(gca,'fontsize',fs)

%% display All water pixels with no high values in NIR

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,waterpixels(waterpixels(:,bn)<meanwp(bn)+2*stdwp(bn),:)')
str = sprintf('Reflectance water L8 image with low values band %i',bn);
title(str,'fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
% ylim([0 0.18])
hold on
plot(L8bands,meanwp+stdwp,'g','linewidth',2)

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,waterpixels(waterpixels(:,bn)>meanwp(bn)+2*stdwp(bn),:)')
str = sprintf('Reflectance water L8 image with high values band %i',bn);
title(str,'fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
% ylim([0 0.18])
hold on
plot(L8bands,meanwp+stdwp,'g','linewidth',2)

%% Pixels that B5 > B3
im = double(imL8crop);
imB5greaterthanB3 = zeros(size(im,1),size(im,2));
cond3 = (im(:,:,5)> im(:,:,3) )& (imL8cropmask ~= 0);
imB5greaterthanB3(cond3)=1;% only high NIR

imB5greaterthanB3mask = imB5greaterthanB3+0.5*repmat(~double(imL8cropmask),[1 1 1]); % for the land appear gray


figure
fs = 15;
set(gcf,'color','white')
imshow(imB5greaterthanB3mask)
str = sprintf('High values for band %i',bn);
title(str,'fontsize',fs)
set(gca,'fontsize',fs)

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,waterpixels(waterpixels(:,5)>waterpixels(:,3),:)')
str = sprintf('Reflectance water L8 image with high values band %i',bn);
title(str,'fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
% ylim([0 0.18])
hold on
plot(L8bands,meanwp+stdwp,'g','linewidth',2)

%% Pixels para incluir en el IGARSS14 abstract
% water pixels reflectance
waterpixelsamples = waterpixels(1:30:end,:);
figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,waterpixelsamples(waterpixelsamples(:,5)<waterpixelsamples(:,3),:)')
str = sprintf('Reflectance water pixels L8');
title(str,'fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
grid on
% ylim([0 .25])

% water pixels Radiance
% figure
% fs = 15;
% set(gcf,'color','white')
% plot(L8bands,waterradpixels(waterpixels(:,5)<waterpixels(:,3),:)')
% str = sprintf('Radiance water pixels L8');
% title(str,'fontsize',fs)
% xlabel('wavelength [\mu m]','fontsize',fs)
% ylabel('Radiance [W/m^2/sr]','fontsize',fs)
% set(gca,'fontsize',fs)



%% LUTs from HydroLight

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

% figure
% fs = 15;
% set(gcf,'color','white')
% plot(wavelength,Rrs)
% title('Rrs','fontsize',fs)
% xlabel('wavelength [\mum]','fontsize',fs)
% ylabel('reflectance','fontsize',fs)
% set(gca,'fontsize',fs)
% xlim([0.4 2.2])
% ylim([0 .2])


LUT = spect_sampL8(Rrs,wavelength);

LUTconcInput = zeros(size(LUT,1),1);
LUTconcDPF = zeros(size(LUT,1),1);

for index = 1:size(LUT,1) 
    
    if strcmp(c{1}(index),'input140929ONTOS')
        LUTconcInput(index)= 1;
    elseif strcmp(c{1}(index),'input140929LONGS')
        LUTconcInput(index)= 2;     
    end
    
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

% Cleaning LUT from repeated values
LUT_All = unique([LUT LUTconcInput LUTconc LUTconcDPF],'rows','stable');
if size(LUT_All,1)~=size(LUT,1)
    LUT             = LUT_All(:,1:7);
    LUTconc         = LUT_All(:,9:11);
    LUTconcInput    = LUT_All(:,8);
    LUTconcDPF      = LUT_All(:,12);
    disp('Repeated Values in LUT. They were deleted!')
end

% rule5 = strcmp(c{1}(:),'input140408ONTNS')& LUTconc(:,1)<10&LUTconc(:,2)<10&LUTconc(:,3)<0.9 &...
%     (strcmp(c{5}(:),'FFbb010.dpf')|strcmp(c{5}(:),'FFbb012.dpf'));
% rule2 = strcmp(c{1}(:),'input140408LONGS')& LUTconc(:,1)>=10&LUTconc(:,2)>=10&LUTconc(:,3)>=0.9&...
%     (strcmp(c{5}(:),'FFbb005.dpf')|strcmp(c{5}(:),'FFbb006.dpf')|strcmp(c{5}(:),'FFbb007.dpf')...
%     |strcmp(c{5}(:),'FFbb008.dpf')|strcmp(c{5}(:),'FFbb009.dpf'));

CHlimit = 28.30;
SMlimit = 9.11;
CDlimit = 0.9819;
DPFlimit = 1.2;

rule5 = strcmp(c{1}(:),'input140929ONTOS')& ...
    LUTconc(:,1)<CHlimit & LUTconc(:,2)<SMlimit & LUTconc(:,3)<CDlimit & ...
    LUTconcDPF(:)<DPFlimit;
rule2 = strcmp(c{1}(:),'input140929LONGS')& ...
    LUTconc(:,1)>=CHlimit & LUTconc(:,2)>=SMlimit & LUTconc(:,3)>=CDlimit & ...
    LUTconcDPF(:)>=DPFlimit;

LUTsmart = LUT(rule5|rule2,:);
LUTconcsmart = LUTconc(rule5|rule2,:);
Inputsmart = LUTconcInput(rule5|rule2);
DPFsmart = LUTconcDPF(rule5|rule2);

LUTlake = LUT(rule5,:);
LUTconclake = LUTconc(rule5,:);
Inputlake = LUTconcInput(rule5);
DPFlake = LUTconcDPF(rule5);

LUTpond = LUT(rule2,:);
LUTconcpond = LUTconc(rule2,:);
Inputpond = LUTconcInput(rule2);
DPFpond = LUTconcDPF(rule2);

WhichLUT =1;

switch WhichLUT
    case 0
        LUTused = LUT;
        LUTconcused = LUTconc;
        LUTInputused = LUTconcInput;
        LUTDPFused = LUTconcDPF;
        fprintf('Using full LUT\n');
        
    case 1
        LUTused = LUTsmart;
        LUTconcused = LUTconcsmart;
        LUTInputused = Inputsmart;
        LUTDPFused = DPFsmart;     
        fprintf('Using smart LUT\n');
        
    case 2
        LUTused = LUTlake;
        LUTconcused = LUTconclake;
        LUTInputused = Inputlake;
        LUTDPFused = DPFlake;     
        fprintf('Using lake LUT\n');     
        
    case 3
        LUTused = LUTpond;
        LUTconcused = LUTconcpond;
        LUTInputused = Inputpond;
        LUTDPFused = DPFpond;     
        fprintf('Using pond LUT\n');         
end



% rule3 =  strcmp(c{1}(:),'input140408LONGS');
% LUTLONGS = LUT(rule3,:);
% LUTconcLONGS = LUTconc(rule3,:);
% 
% rule4 =  strcmp(c{1}(:),'input140408ONTNS');
% LUTONTNS = LUT(rule4,:);
% LUTconcONTNS = LUTconc(rule4,:);

%% LUT Used
figure
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(L8bands,LUTused)
% str1 = sprintf('Reflectance LUT from HydroLight -- %s',LUTname);
% title(str1,'fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('R_{rs} [1/sr]','fontsize',fs)
xlim([0.4 2.5])
grid on

%% Retrieval Best Match %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------------------------------------------------------------')
disp('Running Best Match Routine')
    [XResultsLS,IMatrixLS] = BestMatchRetrieval(waterpixels(:,1:5),LUTused(:,1:5),LUTconcused);
disp('Routine finished Successfully')

XResults = XResultsLS;
IMatrix = IMatrixLS;

%% Retrieval using optimization

format short
CDconc = unique(LUTconcused(:,3))
SMconc = unique(LUTconcused(:,2))
CHconc = unique(LUTconcused(:,1))
%
disp('--------------------------------------------------------------------------')
disp('Running Optimization Routine')
    [XResultsOpt,residualOpt,IMatrixOpt] = opt(waterpixels(:,1:5),LUTused(:,1:5),LUTconcused,LUTInputused,LUTDPFused);
disp('Optimization Routine finished Successfully')

save LSQNONLIN_results.mat XResultsOpt residualOpt IMatrixOpt
%%
load LSQNONLIN_results.mat

XResults = XResultsOpt;
IMatrix = IMatrixOpt;
%% Maps

% to see what kind of input (ONTNS or LONGS) and DPFs were retrieved
InputRet = LUTInputused(IMatrix);
DPFRet = LUTDPFused(IMatrix);

ConcRet = nan(size(masknew,1),5);
ConcRet(masknew==1,:) = [XResults InputRet DPFRet]; % Concentration Retrieved

CHmap  = reshape(ConcRet(:,1),...
    [size(imL8cropmask,1) size(imL8cropmask,2) size(imL8cropmask,3)]);
SMmap   = reshape(ConcRet(:,2),...
    [size(imL8cropmask,1) size(imL8cropmask,2) size(imL8cropmask,3)]);
CDmap = reshape(ConcRet(:,3),...
    [size(imL8cropmask,1) size(imL8cropmask,2) size(imL8cropmask,3)]);

INPUTmap = reshape(ConcRet(:,4),...
    [size(imL8cropmask,1) size(imL8cropmask,2) size(imL8cropmask,3)]);

DPFmap = reshape(ConcRet(:,5),...
    [size(imL8cropmask,1) size(imL8cropmask,2) size(imL8cropmask,3)]);

CHmaplog10 = log10(CHmap);
CHmaplog10(CHmaplog10==-Inf)=-4;
% CHmaplog10masked = bsxfun(@times, CHmaplog10, landmask);

SMmaplog10 = log10(SMmap);
SMmaplog10(SMmaplog10==-Inf)=-4;
% SMmaplog10masked = bsxfun(@times, SMmaplog10, landmask);

CDmaplog10 = log10(CDmap);
CDmaplog10(CDmaplog10==-Inf)=-4;
% CDmaplog10masked = bsxfun(@times, CDmaplog10, landmask);

%% DISSERTACION FIGURES -- log10 scale
% define mask
mask = logical(imL8cropmask);
gray = cat(3,   0.7*ones(size(imL8cropmask)), ...
                0.7*ones(size(imL8cropmask)),...
                0.7*ones(size(imL8cropmask)));

%% ROI RGB image
figure('name',date)
fs = 16;
ms = 16;
set(gcf,'color','white')
imagesc(gray); % display color of the mask first
hold on
h0 = imagesc(impos); % display Map imaage second
set(h0, 'AlphaData', mask) % Apply transparency to the mask
set(gca,'fontsize',fs)
axis equal
axis image
axis off
set(gca, 'Units', 'normalized', 'Position', [0 0.1 1 1])
%% CHL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('name',date)
fs = 16;
ms = 16;
set(gcf,'color','white')
imagesc(gray); % display color of the mask first
hold on
h0 = imagesc(CHmaplog10); % display Map imaage second
set(gca, 'CLim', [min(CHmaplog10(:)), max(CHmaplog10(:))]);
set(h0, 'AlphaData', mask) % Apply transparency to the mask
set(gca,'fontsize',fs)
axis equal
axis image
axis off
h = colorbar;
set(h,'fontsize',fs,'Location','southoutside')
set(h,'Position',[.2 .08 .6 .05])
title(h,'L8 retrieved C_a [mg m^{-3}]','FontSize',fs)
set(gca, 'Units', 'normalized', 'Position', [0 0.1 1 1])
y = get(h,'XTick');
[xmin,xmax] = caxis;
x = [10^xmin 10.^y];
set(h,'XTick',log10(x));
set(h,'XTickLabel',x)

%% SM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('name',date)
fs = 16;
ms = 16;
set(gcf,'color','white')
imagesc(gray); % display color of the mask first
hold on
h0 = imagesc(SMmaplog10); % display Map imaage second
set(gca, 'CLim', [min(SMmaplog10(:)), max(SMmaplog10(:))]);
set(h0, 'AlphaData', mask) % Apply transparency to the mas
set(gca,'fontsize',fs)
axis equal
axis image
axis off
h = colorbar;
set(h,'fontsize',fs,'Location','southoutside')
set(h,'Position',[.2 .1 .6 .05])
title(h,'L8 retrieved TSS [g m^{-3}]','FontSize',fs)
set(gca, 'Units', 'normalized', 'Position', [0 .1 1 1])
y = get(h,'XTick')
[xmin,xmax] = caxis
x = [10^(xmin) 0.1 1 10 10^(xmax)];
set(h,'XTick',log10(x));
set(h,'XTickLabel',x)


%% CDOM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('name',date)
fs = 16;
ms = 16;
set(gcf,'color','white')
imagesc(gray); % display color of the mask first
hold on
h0 = imagesc(CDmaplog10); % display Map imaage second
set(gca, 'CLim', [min(CDmaplog10(:)), max(CDmaplog10(:))]);
set(h0, 'AlphaData', mask) % Apply transparency to the mask
set(gca,'fontsize',fs)
axis equal
axis image
axis off
h = colorbar;
set(h,'fontsize',fs,'Location','southoutside')
set(h,'Position',[.2 .07 .6 .05])
title(h,'L8 retrieved a_{CDOM}(440nm) [1/m]','FontSize',fs)
set(gca, 'Units', 'normalized', 'Position', [0 0.1 1 1])
y = get(h,'XTick') % values in x from 10^x
[xmin,xmax] = caxis
x = [10^xmin 0.2 0.5 10^xmax];
% x = 10.^y;
set(h,'XTick',log10(x));
set(h,'XTickLabel',x)
%% Maps Linear Scale
fs = 30; % font size
cbfs = 15; % colorbar font size

figure('name',date,'Position',get(0,'ScreenSize'))
set(gcf,'color','white')
subplot(2,2,1)
imagesc(impos)
title('RGB image ','fontsize',fs)
set(gca,'fontsize',fs)
axis equal
axis image
axis off

subplot(2,2,2)
imagesc(CHmap)
% caxis([min(XResults(:,1)) max(XResults(:,1))])
title('<CHL>, [\mug/L] ','fontsize',fs)
set(gca,'fontsize',fs)
h = colorbar;
set(h,'fontsize',cbfs)
axis equal
axis image
axis off

subplot(2,2,3)
imagesc(SMmap)
% caxis([min(XResults(:,2)) max(XResults(:,2))])
title('<TSS>, [mg/L]','fontsize',fs)
set(gca,'fontsize',fs)
h = colorbar;
set(h,'fontsize',cbfs)
axis equal
axis image
axis off

subplot(2,2,4)
imagesc(CDmap)
% caxis([min(XResults(:,3)) max(XResults(:,3))])
title('a_{CDOM}(440nm), [1/m] ','fontsize',fs)
set(gca,'fontsize',fs)
h = colorbar;
set(h,'fontsize',cbfs)
axis equal
axis image
axis off

%% Plot Input (ONTNS or LONGS) and DPFs retrieved
figure('name',date,'Position',get(0,'ScreenSize'))
subplot(1,2,1)
set(gcf,'color','white')
imagesc(INPUTmap)
title('INPUT map linear scale','fontsize',fs)
set(gca,'fontsize',fs)
axis equal
axis image
axis off
h = colorbar;
set(h,'fontsize',cbfs)

subplot(1,2,2)
set(gcf,'color','white')
imagesc(DPFmap)
title('DPF map linear scale','fontsize',fs)
set(gca,'fontsize',fs)
axis equal
axis image
axis off
h = colorbar;
set(h,'fontsize',cbfs)

%% Read ROIs from ENVI text Stat file 

pathname = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929/';
% statfilename = 'ROIstats_140929_150408.txt';
statfilename = 'ROIstats_140929_150418CRANB_ONTOS.txt';

fid = fopen([pathname statfilename]);
s = textscan(fid, '%s', 'delimiter', '\n');
% search for "data=":
% idx1 = find(strcmp('Column 7: Mean: SAND~~7',s{1}),1);
idx1 = find(strcmp('Column 7: Mean: IBAYN~~5',s{1}),1);% for MoB-ELM with bright: CRANB and dark: ONTOS

% and read from s{1}(idx1+1:idx2-1) the values using textscan again ...
data = s{1}(idx1+1:size(s{1},1));
fclose(fid);

for index = 1:size(data,1)
    statdata(index,:) = str2num(cell2mat(data(index)));
end

L8bands_ENVI = statdata(:,1);
% Cranb = statdata(:,2);
% LongS = statdata(:,3);
% LongN = statdata(:,4);
% IBayN = statdata(:,5);
% OntOS = statdata(:,6);
% Sand1 = statdata(:,7);

% for MoB-ELM with bright: CRANB and dark: ONTOS
Cranb = statdata(:,2); Cranb(isnan(Cranb(:)))=0;
LongS = statdata(:,3); LongS(isnan(LongS(:)))=0;
LongN = statdata(:,4); LongN(isnan(LongN(:)))=0;
OntOS = statdata(:,5); OntOS(isnan(OntOS(:)))=0;
Sand1 = statdata(:,6); Sand1(isnan(Sand1(:)))=0;
IBayN = statdata(:,7); IBayN(isnan(IBayN(:)))=0;

figure('name',date)
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(L8bands_ENVI,Cranb,'r')
hold on
plot(L8bands_ENVI,LongS,'g')
plot(L8bands_ENVI,LongN,'b')
plot(L8bands_ENVI,OntOS,'c')
plot(L8bands_ENVI,IBayN,'y')
plot(L8bands_ENVI,Sand1,'m')
legend('Cranb','LongS','LongN','OntOS','IBayN','Sand1')
grid on
axis([.4 2.5 0 0.025])


%% Find LONGS
rule6 = strcmp(Inputused(:),'input140929LONGS')&...
    LUTconcused(:,1)==45 & LUTconcused(:,2)==28.00 &...
    LUTconcused(:,3)==1.0025;


% find LongS in waterpixels with index I
[Y,I] = min(sqrt(mean((waterpixels-ones(size(waterpixels,1),1)*LongS').^2,2)));

Inputused(IMatrix(I))
DPFused(IMatrix(I))
LUTconcused(IMatrix(I),:)

LongSconc140929 = [46.10 28.30 0.9819];
LongSconc140929ret = XResults(I,:);

figure('name',date)
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(L8bands,LongS,'.-g')
hold on
plot(L8bands,waterpixels(I,:),'.-r')
plot(L8bands,LUTused(IMatrix(I),:),'.-b')
plot(L8bands,LUTused(rule6,:)','k')
title('LongS','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('R_{rs} [1/m]','fontsize',fs)
legend('ROI','Closest L8','ret. from HL','DPF LUT')
grid on
xlim([0.4 2.5])

%% Find LONGN
rule6 = strcmp(Inputused(:),'input140929LONGS')&...
    LUTconcused(:,1)==50 & LUTconcused(:,2)==16.7 &...
    LUTconcused(:,3)==1.0025;


% find LongN in waterpixels with index I
[Y,I] = min(sqrt(mean((waterpixels-ones(size(waterpixels,1),1)*LongN').^2,2)));

Inputused(IMatrix(I))
DPFused(IMatrix(I))
LUTconcused(IMatrix(I),:)

LongNconc140929 = [47.90 16.7 1.0194];
LongNconc140929ret = XResults(I,:);

figure('name',date) 
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(L8bands,LongN,'.-g')
hold on
plot(L8bands,waterpixels(I,:),'.-r')
plot(L8bands,LUTused(IMatrix(I),:),'.-b')
plot(L8bands,LUTused(rule6,:)','k')
title('LongN','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('R_{rs} [1/m]','fontsize',fs)
legend('ROI','Closest L8','ret. from HL','DPF LUT')
grid on
xlim([0.4 2.5])

%% Find Cranb
rule6 = strcmp(Inputused(:),'input140929LONGS')&...
    LUTconcused(:,1)==60 & LUTconcused(:,2)==30.7 &...
    LUTconcused(:,3)==1.0025;

% find Cranb in waterpixels with index I
[Y,I] = min(sqrt(mean((waterpixels-ones(size(waterpixels,1),1)*Cranb').^2,2)));

Inputused(IMatrix(I))
DPFused(IMatrix(I))
LUTconcused(IMatrix(I),:)

Cranbconc140929 = [58.3 30.70 0.9297];
Cranbconc140929ret = XResults(I,:);

figure('name',date)
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(L8bands,Cranb,'.-g')
hold on
plot(L8bands,waterpixels(I,:),'.-r')
plot(L8bands,LUTused(IMatrix(I),:),'.-b')
plot(L8bands,LUTused(rule6,:)','k')
title('Cranb','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('R_{rs} [1/m]','fontsize',fs)
legend('ROI','Closest L8','ret. from HL','DPF LUT')
grid on
xlim([0.4 2.5])

%% Find IBayN
rule6 = strcmp(Inputused(:),'input140929LONGS')&...
    LUTconcused(:,1)==30.0 & LUTconcused(:,2)==9.1100 &...
    LUTconcused(:,3)==1.0025;

% find IBayN in waterpixels with index I
[Y,I] = min(sqrt(mean((waterpixels-ones(size(waterpixels,1),1)*IBayN').^2,2)));

Inputused(IMatrix(I))
DPFused(IMatrix(I))
LUTconcused(IMatrix(I),:)

IBayNconc140929 = [28.30 9.11 1.0025];
IBayNconc140929ret = XResults(I,:);

figure('name',date)
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(L8bands,IBayN,'.-g')
hold on
plot(L8bands,waterpixels(I,:),'.-r')
plot(L8bands,LUTused(IMatrix(I),:),'.-b')
plot(L8bands,LUTused(rule6,:)','k')
title('IbayN','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('R_{rs} [1/m]','fontsize',fs)
legend('ROI','Closest L8','ret. from HL','DPF LUT')
grid on
xlim([0.4 2.5])

%% Find OntOS
rule6 = strcmp(Inputused(:),'input140929ONTOS')&...
    LUTconcused(:,1)==2.1 & LUTconcused(:,2)==1.4 &...
    LUTconcused(:,3)==0.0954;

% find OntOS in waterpixels with index I
[Y,I] = min(sqrt(mean((waterpixels-ones(size(waterpixels,1),1)*OntOS').^2,2)));

Inputused(IMatrix(I))
DPFused(IMatrix(I))
LUTconcused(IMatrix(I),:)

OntOSconc140929 = [2.10 1.4 0.0954];
OntOSconc140929ret = XResults(I,:);

figure('name',date)
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(L8bands,OntOS,'-g')
hold on
plot(L8bands,waterpixels(I,:),'.-r')
plot(L8bands,LUTused(IMatrix(I),:),'.-b')
plot(L8bands,LUTused(rule6,:)','k')
title('OntOS','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('R_{rs} [1/m]','fontsize',fs)
legend('ROI','Closest L8','ret. from HL','DPF LUT')
grid on
xlim([0.4 2.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scatter plot

fs = 25;
ms = 25; %marker size

figure('name',date,'Position',get(0,'ScreenSize'))
subplot(1,3,1)
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(LongSconc140929(1),LongSconc140929ret(1),'.r','MarkerSize', ms);
hold on
plot(LongNconc140929(1),LongNconc140929ret(1),'.k','MarkerSize', ms);
plot(Cranbconc140929(1),Cranbconc140929ret(1),'.b','MarkerSize', ms);
plot(IBayNconc140929(1),IBayNconc140929ret(1),'.g','MarkerSize', ms);
plot(OntOSconc140929(1),OntOSconc140929ret(1),'.m','MarkerSize', ms);
maxconcChl = 160;
plot([0 maxconcChl],[0 maxconcChl],'k')
axis equal
ylim([0 maxconcChl])
xlim([0 maxconcChl])
title('<Chl>, [\mug/L]','fontsize',fs)
xlabel('measured <Chl> [\mug/L] ','fontsize',fs)
ylabel('L8 retrieved <Chl> [\mug/L]','fontsize',fs)
legend('LONGS','LONGN','CRANB','IBAYN','ONTOS')

% figure
subplot(1,3,2)
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(LongSconc140929(2),LongSconc140929ret(2),'.r','MarkerSize', ms);
hold on
plot(LongNconc140929(2),LongNconc140929ret(2),'.k','MarkerSize', ms);
plot(Cranbconc140929(2),Cranbconc140929ret(2),'.b','MarkerSize', ms);
plot(IBayNconc140929(2),IBayNconc140929ret(2),'.g','MarkerSize', ms);
plot(OntOSconc140929(2),OntOSconc140929ret(2),'.m','MarkerSize', ms);
maxconcTSS = 60;
plot([0 maxconcTSS],[0 maxconcTSS],'k')
axis equal
ylim([0 maxconcTSS])
xlim([0 maxconcTSS])
title('<TSS>, [mg/L]','fontsize',fs)
xlabel('measured <TSS> [mg/L] ','fontsize',fs)
ylabel('L8 retrieved <TSS> [m/L]','fontsize',fs)
% legend('LONGS','LONGN','CRANB','IBAYN','ONTOS')

subplot(1,3,3)
% figure
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(LongSconc140929(3),LongSconc140929ret(3),'.r','MarkerSize', ms);
hold on
plot(LongNconc140929(3),LongNconc140929ret(3),'.k','MarkerSize', ms);
plot(Cranbconc140929(3),Cranbconc140929ret(3),'.b','MarkerSize', ms);
plot(IBayNconc140929(3),IBayNconc140929ret(3),'.g','MarkerSize', ms);
plot(OntOSconc140929(3),OntOSconc140929ret(3),'.m','MarkerSize', ms);
maxconcCDOM = 1.5;
plot([0 maxconcCDOM],[0 maxconcCDOM],'k')
axis equal
ylim([0 maxconcCDOM])
xlim([0 maxconcCDOM])
title('a_{CDOM}(440nm), [1/m]','fontsize',fs)
xlabel('measured a_{CDOM}(440nm) [1/m]','fontsize',fs)
ylabel('retrieved a_{CDOM}(440nm) [1/m]','fontsize',fs)
% legend('LONGS','LONGN','CRANB','IBAYN','ONTOS')



%% RS of ENVIRONMENT PAPER FIGURES
%% CHL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('name',date)
fs = 16;
ms = 16;
set(gcf,'color','white')
imagesc(CHmap)
set(gca,'fontsize',fs)
axis equal
axis image
axis off
h = colorbar;
set(h,'fontsize',fs,'Location','southoutside')
set(h,'Position',[.2 .05 .6 .05])
title(h,'L8 retrieved C_a [mg m^{-3}]','FontSize',fs)
set(gca, 'Units', 'normalized', 'Position', [0 0.07 1 1])
%%
figure('name',date)
fs = 20;
ms = 25;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(LongSconc140929(1),LongSconc140929ret(1),'.r','MarkerSize', ms);
hold on
plot(LongNconc140929(1),LongNconc140929ret(1),'.k','MarkerSize', ms);
plot(Cranbconc140929(1),Cranbconc140929ret(1),'.b','MarkerSize', ms);
plot(IBayNconc140929(1),IBayNconc140929ret(1),'.g','MarkerSize', ms);
plot(OntOSconc140929(1),OntOSconc140929ret(1),'.m','MarkerSize', ms);
maxconcChl = 100;
plot([0 maxconcChl],[0 maxconcChl],'--k')
axis equal
ylim([0 maxconcChl])
xlim([0 maxconcChl])
xlabel('measured C_a [mg m^{-3}] ','fontsize',fs,'Position',[55 -10])
ylabel('L8 retrieved C_a [mg m^{-3}]','fontsize',fs)
legend('LONGS','LONGN','CRANB','IBAYN','ONTOS','Location','best')
set(gca,'OuterPosition',[0 0.05 1 1])
% save('CHL.txt','-ascii','-double','-tabs','CHLmap')
%% SM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('name',date)
fs = 16;
ms = 16;
set(gcf,'color','white')
imagesc(SMmap)
set(gca,'fontsize',fs)
axis equal
axis image
axis off
h = colorbar;
set(h,'fontsize',fs,'Location','southoutside')
set(h,'Position',[.2 .055 .6 .05])
title(h,'L8 retrieved TSS [g m^{-3}]','FontSize',fs)
set(gca, 'Units', 'normalized', 'Position', [0 0.07 1 1])
%%
figure('name',date)
fs = 20;
ms = 25;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(LongSconc140929(2),LongSconc140929ret(2),'.r','MarkerSize', ms);
hold on
plot(LongNconc140929(2),LongNconc140929ret(2),'.k','MarkerSize', ms);
plot(Cranbconc140929(2),Cranbconc140929ret(2),'.b','MarkerSize', ms);
plot(IBayNconc140929(2),IBayNconc140929ret(2),'.g','MarkerSize', ms);
plot(OntOSconc140929(2),OntOSconc140929ret(2),'.m','MarkerSize', ms);
maxconcTSS = 40;
plot([0 maxconcTSS],[0 maxconcTSS],'--k')
axis equal
ylim([0 maxconcTSS])
xlim([0 maxconcTSS])
xlabel('measured TSS [g m^{-3}] ','fontsize',fs)
ylabel('L8 retrieved TSS [g m^{-3}]','fontsize',fs)
legend('LONGS','LONGN','CRANB','IBAYN','ONTOS','Location','best')
% save('TSS.txt','-ascii','-double','-tabs','SMmap')
%% CDOM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('name',date)
fs = 16;
ms = 16;
set(gcf,'color','white')
imagesc(CDmap)
set(gca,'fontsize',fs)
axis equal
axis image
axis off
h = colorbar;
set(h,'fontsize',fs,'Location','southoutside')
set(h,'Position',[.2 .06 .65 .05])
title(h,'L8 retrieved a_{CDOM}(440nm) [1/m]','FontSize',fs)
set(gca, 'Units', 'normalized', 'Position', [0 0.09 1 1])
%%
figure('name',date)
fs = 20;
ms = 25;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(LongSconc140929(3),LongSconc140929ret(3),'.r','MarkerSize', ms);
hold on
plot(LongNconc140929(3),LongNconc140929ret(3),'.k','MarkerSize', ms);
plot(Cranbconc140929(3),Cranbconc140929ret(3),'.b','MarkerSize', ms);
plot(IBayNconc140929(3),IBayNconc140929ret(3),'.g','MarkerSize', ms);
plot(OntOSconc140929(3),OntOSconc140929ret(3),'.m','MarkerSize', ms);
maxconcCDOM = 1.5;
plot([0 maxconcCDOM],[0 maxconcCDOM],'--k')
axis equal
ylim([0 maxconcCDOM])
xlim([0 maxconcCDOM])
xlabel('measured a_{CDOM}(440nm) [1/m]','fontsize',fs)
ylabel('retrieved a_{CDOM}(440nm) [1/m]','fontsize',fs)
legend('LONGS','LONGN','CRANB','IBAYN','ONTOS','Location','best')
% save('CDOM.txt','-ascii','-double','-tabs','CDOMmap')
%% Plot Input (ONTNS or LONGS) and DPFs retrieved
figure('name',date)
subplot(1,2,1)
set(gcf,'color','white')
imagesc(INPUTmap)
title('INPUT map linear scale','fontsize',fs)
set(gca,'fontsize',fs)
axis equal
axis image
axis off
h = colorbar;
set(h,'fontsize',cbfs)

subplot(1,2,2)
set(gcf,'color','white')
imagesc(DPFmap)
title('DPF map linear scale','fontsize',fs)
set(gca,'fontsize',fs)
axis equal
axis image
axis off
h = colorbar;
set(h,'fontsize',cbfs)


%% Mapping Concentrations linear scale
fs = 30; % font size
cbfs = 15; % colorbar font size

figure('name',date)
set(gcf,'color','white')
imagesc(impos)
title('RGB image ','fontsize',fs)
set(gca,'fontsize',fs)
axis equal
axis image
axis off

figure('name',date)
set(gcf,'color','white')
imagesc(CHLmap)
title('<CHL>, \mug/L','fontsize',fs)
set(gca,'fontsize',fs)
axis equal
axis image
axis off
colormap('gray')
h = colorbar;
set(h,'fontsize',cbfs)

% save('CHL.txt','-ascii','-double','-tabs','CHLmap')

figure('name',date)
set(gcf,'color','white')
imagesc(SMmap)
title('<TSS>, mg/L','fontsize',fs)
set(gca,'fontsize',fs)
axis equal
axis image
axis off
colormap('gray')
h = colorbar;
set(h,'fontsize',cbfs)

% save('TSS.txt','-ascii','-double','-tabs','SMmap')

figure('name',date)
set(gcf,'color','white')
imagesc(CDOMmap)
title('a_{CDOM}(440), 1/m','fontsize',fs)
set(gca,'fontsize',fs)
axis equal
axis image
axis off
colormap('gray')
h = colorbar;
set(h,'fontsize',cbfs)

% save('CDOM.txt','-ascii','-double','-tabs','CDOMmap')

%% Mapping Concentrations log scale


figure('name',date)
fs = 15;
set(gcf,'color','white')
imagesc(CHLmaplog10)
title('CHL map log scale','fontsize',fs)
set(gca,'fontsize',fs)
axis('equal')
colorbar

figure('name',date)
fs = 15;
set(gcf,'color','white')
imagesc(SMmaplog10)
title('SM map log scale','fontsize',fs)
set(gca,'fontsize',fs)
axis('equal')
colorbar

figure('name',date)
fs = 15;
set(gcf,'color','white')
imagesc(CDOMmaplog10)
title('CDOM map log scale','fontsize',fs)
set(gca,'fontsize',fs)
axis('equal')
colorbar


%% Save maps as TIFF

t1 = Tiff('LC80160302013262LGN00CHLmap.tif','w');
t2 = Tiff('LC80160302013262LGN00SMmap.tif','w');
t3 = Tiff('LC80160302013262LGN00CDOMmap.tif','w');

tagstruct.ImageLength = size(imL8crop,1);
tagstruct.ImageWidth = size(imL8crop,2);
tagstruct.Photometric = Tiff.Photometric.LinearRaw;
tagstruct.SampleFormat = 3; %'IEEEFP'
tagstruct.BitsPerSample = 64;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB';

t1.setTag(tagstruct)
t1.write(CHLmap)
t1.close();

t2.setTag(tagstruct)
t2.write(SMmap)
t2.close();

t3.setTag(tagstruct)
t3.write(CDOMmap)
t3.close();

% imagesc(imread('LC80160302013262LGN00CHLmap.tif'));
% axis equal

%% Histogram of concentrations log scale
nbins = 50;
figure
subplot(1,3,1)
fs = 15;
set(gcf,'color','white')
hist(CHLmaplog10(CHLmaplog10~=-Inf),nbins)
title('CHL','fontsize',fs)
set(gca,'fontsize',fs)

subplot(1,3,2)
fs = 15;
set(gcf,'color','white')
hist(SMmaplog10(SMmaplog10~=-Inf),nbins)
title('SM','fontsize',fs)
set(gca,'fontsize',fs)

subplot(1,3,3)
fs = 15;
set(gcf,'color','white')
hist(CDOMmaplog10(CDOMmaplog10~=-Inf),nbins)
title('CDOM','fontsize',fs)
set(gca,'fontsize',fs)
%% Histogram of concentrations linear scale
nbins = 50;
figure
subplot(1,3,1)
fs = 15;
set(gcf,'color','white')
hist(XResults(:,1),nbins)
title('CHL','fontsize',fs)
set(gca,'fontsize',fs)

subplot(1,3,2)
fs = 15;
set(gcf,'color','white')
hist(XResults(:,2),nbins)
title('SM','fontsize',fs)
set(gca,'fontsize',fs)

subplot(1,3,3)
fs = 15;
set(gcf,'color','white')
hist(XResults(:,3),nbins)
title('CDOM','fontsize',fs)
set(gca,'fontsize',fs)



%% Figure for Dr. John - 09/23/13 Log Scale
figure
fs = 20;
set(gcf,'color','white')
subplot(2,2,1)
imagesc(impos)
title('RGB image ','fontsize',fs)
set(gca,'fontsize',fs)
axis equal
axis off

subplot(2,2,2)
imagesc(CHLmaplog10)
title('CHL map ','fontsize',fs)
set(gca,'fontsize',fs)
colorbar
axis equal
axis off

subplot(2,2,3)
imagesc(SMmaplog10)
title('SM map ','fontsize',fs)
set(gca,'fontsize',fs)
colorbar
axis equal
axis off

subplot(2,2,4)
imagesc(CDOMmaplog10)
title('CDOM map ','fontsize',fs)
set(gca,'fontsize',fs)
colorbar
axis equal
axis off

%% Create for ENVI library
% % C = {c{1}(rule1|rule2) cellstr(num2str(LUTconc(rule1|rule2,:))) c{5}(rule1|rule2)};
% C = {char(c{1}(rule5|rule2)) ...
%     num2str(LUTconc(rule5|rule2,:)) ...
%     char(c{5}(rule5|rule2))};
% 
% C = cellstr(C);
% 
% filename = 'file.txt';
% fid = fopen(filename,'wt');
% 
% for index = 1:size(C{1},1)
% fprintf('%s %s %s\r\n',C{1}(index,:),C{2}(index,:),C{3}(index,:));
% end
% fclose(fid);
% 
% % tt = [L8bands LUTsmart];
% save('LUTtest.txt','tt','-ascii')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test the optimization algorithm

format short
CDconc = unique(LUTconc(:,3))
SMconc = unique(LUTconc(:,2))
CHconc = unique(LUTconc(:,1))
%
disp('--------------------------------------------------------------------------')
disp('Running Optimization Routine')
    [XResultstest,residual] = opt(LUTused(:,1:5),LUTused(:,1:5),LUTconcused,LUTInputused,LUTDPFused);
disp('Optimization Routine finished Successfully')

%% E_RMS
disp('--------------------------------------------------')
E_CH = sqrt(sum((XResultstest(:,1)-LUTconcused(:,1)).^2)/size(XResultstest,1));
E_CH = E_CH*100/max(CHconc);
str = sprintf('E_CH = %2.2f %%',E_CH);
disp(str)

E_SM = sqrt(sum((XResultstest(:,2)-LUTconcused(:,2)).^2)/size(XResultstest,1));
E_SM = E_SM*100/max(SMconc);
str = sprintf('E_SM = %2.2f %%',E_SM);
disp(str)

E_CD = sqrt(sum((XResultstest(:,3)-LUTconcused(:,3)).^2)/size(XResultstest,1));
E_CD = E_CD*100/max(CDconc);
str = sprintf('E_CD = %2.2f %%',E_CD);
disp(str)

%% Residual Histogram
figure
fs = 15;
set(gcf,'color','white')
subplot(2,3,1)
set(gca,'fontsize',fs)
hist(residual(:,1),20)
title('band 1')

subplot(2,3,2)
set(gca,'fontsize',fs)
hist(residual(:,2),20)
title('band 2')

subplot(2,3,3)
set(gca,'fontsize',fs)
hist(residual(:,3),20)
title('band 3')

subplot(2,3,4)
set(gca,'fontsize',fs)
hist(residual(:,4),20)
title('band 4')

subplot(2,3,5)
set(gca,'fontsize',fs)
hist(residual(:,5),20)
title('band 5')

%% Display data vs retrieved
% CHL
figure
fs = 15;
set(gcf,'color','white')
plot(LUTconcused(:,1),XResultstest(:,1),'.')
hold on
maxconcCHL = 150;
plot([0 maxconcCHL],[0 maxconcCHL],'--k')
axis equal
ylim([0 maxconcCHL])
xlim([0 maxconcCHL])
title('CHL Real vs retrieved','fontsize',fs)
xlabel('real','fontsize',fs)
ylabel('retrieved','fontsize',fs)
set(gca,'fontsize',fs)


% SM
figure
fs = 15;
set(gcf,'color','white')
plot(LUTconcused(:,2),XResultstest(:,2),'.')
hold on
maxconcSM = 60;
plot([0 maxconcSM],[0 maxconcSM],'--k')
axis equal
ylim([0 maxconcSM])
xlim([0 maxconcSM])
title('SM Real vs retrieved','fontsize',fs)
xlabel('real','fontsize',fs)
ylabel('retrieved','fontsize',fs)
set(gca,'fontsize',fs)

% CDOM
figure
fs = 15;
set(gcf,'color','white')
plot(LUTconcused(:,3),XResultstest(:,3),'.')
hold on
maxconcCD = 1.4;
plot([0 maxconcCD],[0 maxconcCD],'--k')
axis equal
ylim([0 maxconcCD])
xlim([0 maxconcCD])
title('CDOM Real vs retrieved','fontsize',fs)
xlabel('real','fontsize',fs)
ylabel('retrieved','fontsize',fs)
set(gca,'fontsize',fs)
%% Display data vs retrieved
% CHL
figure
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(LUTconcused(:,1))
hold on
plot(XResultstest(:,1),'--r')
title('CHL')
legend('real','retrieved')
xlabel('ith element')
ylabel('C_a')
grid on

% SM
figure
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(LUTconcused(:,2))
hold on
plot(XResultstest(:,2),'--r')
title('SM')
legend('real','retrieved')
xlabel('ith element')
ylabel('TSS')
grid on

% CDOM
figure
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(LUTconcused(:,3))
hold on
plot(XResultstest(:,3),'--r')
title('CDOM')
legend('real','retrieved')
xlabel('ith element')
ylabel('a_{CDOM}(440)')
grid on
%%
figure
bn = 1;
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(XResultsOpt(:,bn),XResultsLS(:,bn),'.')
axis equal
xlabel('C_a NonLin')
ylabel('C_a Least Squared')
Concmax = max([XResultsOpt(:,bn);XResultsLS(:,bn)]);
hold on 
plot([0 Concmax*1.1],[0 Concmax*1.1],'k')
axis([0 Concmax*1.1 0 Concmax*1.1])
%%
figure
bn = 2;
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(XResultsOpt(:,bn),XResultsLS(:,bn),'.')
axis equal
xlabel('TSS NonLin')
ylabel('TSS Least Squared')
Concmax = max([XResultsOpt(:,bn);XResultsLS(:,bn)]);
hold on 
plot([0 Concmax*1.1],[0 Concmax*1.1],'k')
axis([0 Concmax*1.1 0 Concmax*1.1])
%%
figure
bn = 3;
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(XResultsOpt(:,bn),XResultsLS(:,bn),'.')
axis equal
xlabel('CDOM NonLin')
ylabel('CDOM Least Squared')
Concmax = max([XResultsOpt(:,bn);XResultsLS(:,bn)]);
hold on 
plot([0 Concmax*1.1],[0 Concmax*1.1],'k')
axis([0 Concmax*1.1 0 Concmax*1.1])
