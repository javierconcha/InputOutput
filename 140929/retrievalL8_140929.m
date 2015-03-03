% Retrieve the concentration from a image using a LUT from HL for L8 image of 140929.
% Version 3.0
% Created by Javier A. Concha
% 03/02/15
addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/retrieval/')
addpath('/Users/javier/Downloads/mtit')
cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/140929
%% L8 image cropped
% From ELM using new L8 reflectance product
folderpath = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/LC80170302014272LGN00/LC80170302014272LGN00_ROI_Rrstif.tif';
filename = '';

filepath = [folderpath filename];

[imL8crop, cmap] = imread(filepath);
INFO = imfinfo(filepath);

%%%% Mask
imL8cropmask = imread(...
    '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/LC80170302014272LGN00/LC80170302014272LGN00_ROImask.tif');

imL8cropmask(imL8cropmask>0)=1;


L8bands = [0.4430,0.4826,0.5613,0.6546,0.8646,1.6090,2.2010];

%%%% water pixels. Convert each band in columns.

imnew = reshape(imL8crop,[size(imL8crop,1)*size(imL8crop,2) size(imL8crop,3)]);
masknew = reshape(imL8cropmask,[size(imL8cropmask,1)*size(imL8cropmask,2) size(imL8cropmask,3)]);



waterpixels = imnew(masknew==1,:);
waterpixels = double(waterpixels);

%% added 01-11-14. plot radiance curves
% radiance image
imL8radcrop = imread(...
    '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images//LC80170302014272LGN00/LC80170302014272LGN00_ROItif.tif');
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

%%
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
bn = 7;
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

% display All water pixels with no high values in NIR

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
% % new 03/02/15, spectrally sampling made in matlab
LUTfilename1 = 'Rvector140929_150301.txt'; 
LUTconcfilename1 = 'concentration_list140929_150301.txt';

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

% rule5 = strcmp(c{1}(:),'input140408ONTNS')& LUTconc(:,1)<10&LUTconc(:,2)<10&LUTconc(:,3)<0.9 &...
%     (strcmp(c{5}(:),'FFbb010.dpf')|strcmp(c{5}(:),'FFbb012.dpf'));
% rule2 = strcmp(c{1}(:),'input140408LONGS')& LUTconc(:,1)>=10&LUTconc(:,2)>=10&LUTconc(:,3)>=0.9&...
%     (strcmp(c{5}(:),'FFbb005.dpf')|strcmp(c{5}(:),'FFbb006.dpf')|strcmp(c{5}(:),'FFbb007.dpf')...
%     |strcmp(c{5}(:),'FFbb008.dpf')|strcmp(c{5}(:),'FFbb009.dpf'));

rule5 = strcmp(c{1}(:),'input140929ONTNS')& LUTconc(:,1)<10&LUTconc(:,2)<10&LUTconc(:,3)<0.9;
rule2 = strcmp(c{1}(:),'input140929LONGS')& LUTconc(:,1)>=10&LUTconc(:,2)>=10&LUTconc(:,3)>=0.9;

LUTsmart = LUT(rule5|rule2,:);
LUTconcsmart = LUTconc(rule5|rule2,:);
Inputsmart = c{1}(rule5|rule2);
DPFsmart = c{5}(rule5|rule2);

LUTlake = LUT(rule5,:);
LUTconclake = LUTconc(rule5,:);
Inputlake = c{1}(rule5);
DPFlake = c{5}(rule5);

LUTpond = LUT(rule2,:);
LUTconcpond = LUTconc(rule2,:);
Inputpond = c{1}(rule2);
DPFpond = c{5}(rule2);

WhichLUT =0;

switch WhichLUT
    case 0
        LUTused = LUT;
        LUTconcused = LUTconc;
        Inputused = c{1};
        DPFused = c{5};
        fprintf('Using full LUT\n');
        
    case 1
        LUTused = LUTsmart;
        LUTconcused = LUTconcsmart;
        Inputused = Inputsmart;
        DPFused = DPFsmart;     
        fprintf('Using smart LUT\n');
        
    case 2
        LUTused = LUTlake;
        LUTconcused = LUTconclake;
        Inputused = Inputlake;
        DPFused = DPFlake;     
        fprintf('Using lake LUT\n');     
        
    case 3
        LUTused = LUTpond;
        LUTconcused = LUTconcpond;
        Inputused = Inputpond;
        DPFused = DPFpond;     
        fprintf('Using pond LUT\n');         
        
end



% rule3 =  strcmp(c{1}(:),'input140408LONGS');
% LUTLONGS = LUT(rule3,:);
% LUTconcLONGS = LUTconc(rule3,:);
% 
% rule4 =  strcmp(c{1}(:),'input140408ONTNS');
% LUTONTNS = LUT(rule4,:);
% LUTconcONTNS = LUTconc(rule4,:);

%% Display LUTsmart
figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,LUTused)
title('Reflectance LUT from HydroLight','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
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



%%
% % %% Test the optimization algorhythm
% % LUTconc = LUTconc;
% % CDOMconc = unique(LUTconc(:,3))
% % SMconc   = unique(LUTconc(:,2))
% % CHLconc  = unique(LUTconc(:,1))
% % 
% % disp('--------------------------------------------------------------------------')
% % disp('Running Optimization Routine')
% %     [XResultstest,residual] = opt(LUT(:,1:5),LUT(:,1:5),LUTconc);
% % disp('Optimization Routine finished Successfully')
% % 
% % 
% % 
% % 
% % % E_RMS
% % disp('--------------------------------------------------')
% % E_Chl = sqrt(sum((XResultstest(:,1)-LUTconc(:,1)).^2)/size(XResultstest,1));
% % E_Chl = E_Chl*100/68;
% % str = sprintf('E_Chl  = %2.2f %%',E_Chl);
% % disp(str)
% % 
% % E_SM = sqrt(sum((XResultstest(:,2)-LUTconc(:,2)).^2)/size(XResultstest,1));
% % E_SM = E_SM*100/24;
% % str = sprintf('E_SM   = %2.2f %%',E_SM);
% % disp(str)
% % 
% % E_CDOM = sqrt(sum((XResultstest(:,3)-LUTconc(:,3)).^2)/size(XResultstest,1));
% % E_CDOM = E_CDOM*100/14;
% % str = sprintf('E_CDOM = %2.2f %%',E_CDOM);
% % disp(str)
% % 
% % %% Residual Histogram
% % figure
% % set(gcf,'color','white')
% % subplot(2,3,1)
% % hist(residual(:,1))
% % title('band 1')
% % 
% % subplot(2,3,2)
% % hist(residual(:,2))
% % title('band 2')
% % 
% % subplot(2,3,3)
% % hist(residual(:,3))
% % title('band 3')
% % 
% % subplot(2,3,4)
% % hist(residual(:,4))
% % title('band 4')
% % 
% % subplot(2,3,5)
% % hist(residual(:,5))
% % title('band 5')
% % %% Display data vs retrieved
% % 
% % figure
% % fs = 15;
% % set(gcf,'color','white')
% % plot(LUTconc(:,1),XResultstest(:,1),'.')
% % xLimits = get(gca,'XLim');  %# Get the range of the x axis
% % yLimits = get(gca,'YLim');  %# Get the range of the y axis
% % hold on
% % plot(xLimits,xLimits,'k')
% % ylim(xLimits)
% % xlim(xLimits)
% % title('CHL Real vs retrieved','fontsize',fs)
% % xlabel('real','fontsize',fs)
% % ylabel('retrieved','fontsize',fs)
% % set(gca,'fontsize',fs)
% % axis equal
% % 
% % figure
% % fs = 15;
% % set(gcf,'color','white')
% % plot(LUTconc(:,2),XResultstest(:,2),'.')
% % xLimits = get(gca,'XLim');  %# Get the range of the x axis
% % yLimits = get(gca,'YLim');  %# Get the range of the y axis
% % hold on
% % plot(xLimits,xLimits,'k')
% % ylim(xLimits)
% % xlim(xLimits)
% % title('SM Real vs retrieved','fontsize',fs)
% % xlabel('real','fontsize',fs)
% % ylabel('retrieved','fontsize',fs)
% % set(gca,'fontsize',fs)
% % axis equal
% % 
% % figure
% % fs = 15;
% % set(gcf,'color','white')
% % plot(LUTconc(:,3),XResultstest(:,3),'.')
% % xLimits = get(gca,'XLim');  %# Get the range of the x axis
% % yLimits = get(gca,'YLim');  %# Get the range of the y axis
% % hold on
% % plot(xLimits,xLimits,'k')
% % ylim(xLimits)
% % xlim(xLimits)
% % title('CDOM Real vs retrieved','fontsize',fs)
% % xlabel('real','fontsize',fs)
% % ylabel('retrieved','fontsize',fs)
% % set(gca,'fontsize',fs)
% % axis equal
% % 
% % %% Retrieval Opt, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % disp('--------------------------------------------------------------------------')
% % disp('Running Optimization Routine')
% %     XResults = opt(waterpixels(:,1:5),LUT(:,1:5),LUTconc);
% % disp('Optimization Routine finished Successfully')

%% Retrieval Best Match %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------------------------------------------------------------')
disp('Running Best Match Routine')
    [XResults,IMatrix] = BestMatchRetrieval(waterpixels(:,1:5),LUTused(:,1:5),LUTconcused);
disp('Routine finished Successfully')
% to see what kind of input (ONTNS or LONGS) and DPFs were retrieved

InputType = zeros(size(IMatrix,1),1);
DPFType = zeros(size(IMatrix,1),1);

for index = 1:size(IMatrix,1) 
    if strcmp(Inputused(IMatrix(index)),'input140929ONTNS')
        InputType(index)= 1;
    elseif strcmp(Inputused(IMatrix(index)),'input140929LONGS')
        InputType(index)= 2;     
    end
        
    if strcmp(DPFused(IMatrix(index)),'FFbb005.dpf')
        DPFType(index)= 0.5;
    elseif strcmp(DPFused(IMatrix(index)),'FFbb006.dpf')
        DPFType(index)= 0.6; 
    elseif strcmp(DPFused(IMatrix(index)),'FFbb007.dpf')
        DPFType(index)= 0.7;
    elseif strcmp(DPFused(IMatrix(index)),'FFbb008.dpf')
        DPFType(index)= 0.8;
    elseif strcmp(DPFused(IMatrix(index)),'FFbb009.dpf')
        DPFType(index)= 0.9;
    elseif strcmp(DPFused(IMatrix(index)),'FFbb010.dpf')
        DPFType(index)= 1.0; 
    elseif strcmp(DPFused(IMatrix(index)),'FFbb012.dpf')
        DPFType(index)= 1.2;
    elseif strcmp(DPFused(IMatrix(index)),'FFbb014.dpf')
        DPFType(index)= 1.4;
    elseif strcmp(DPFused(IMatrix(index)),'FFbb016.dpf')
        DPFType(index)= 1.6;
    end
end

% Maps
ConcRet = zeros(size(masknew,1),5);
ConcRet(masknew==1,:) = [XResults InputType DPFType]; % Concentration Retrieved

CHLmap  = reshape(ConcRet(:,1),...
    [size(imL8cropmask,1) size(imL8cropmask,2) size(imL8cropmask,3)]);
SMmap   = reshape(ConcRet(:,2),...
    [size(imL8cropmask,1) size(imL8cropmask,2) size(imL8cropmask,3)]);
CDOMmap = reshape(ConcRet(:,3),...
    [size(imL8cropmask,1) size(imL8cropmask,2) size(imL8cropmask,3)]);

INPUTmap = reshape(ConcRet(:,4),...
    [size(imL8cropmask,1) size(imL8cropmask,2) size(imL8cropmask,3)]);

DPFmap = reshape(ConcRet(:,5),...
    [size(imL8cropmask,1) size(imL8cropmask,2) size(imL8cropmask,3)]);

CHLmaplog10 = log10(CHLmap);
CHLmaplog10(CHLmaplog10==-Inf)=-4;
% CHLmaplog10masked = bsxfun(@times, CHLmaplog10, landmask);

SMmaplog10 = log10(SMmap);
SMmaplog10(SMmaplog10==-Inf)=-4;
% SMmaplog10masked = bsxfun(@times, SMmaplog10, landmask);

CDOMmaplog10 = log10(CDOMmap);
CDOMmaplog10(CDOMmaplog10==-Inf)=-4;
% CDOMmaplog10masked = bsxfun(@times, CDOMmaplog10, landmask);

%% Find LONGS
rule6 = strcmp(Inputused(:),'input140408LONGS')&...
    LUTconcused(:,1)==110 & LUTconcused(:,2)==45 &...
    LUTconcused(:,3)==1.2;

LongS = [0.009494 0.014197 0.032298 0.020612 ...
    0.009200 0.000870 0.000533];

% LongS = [0.010427 0.015314 0.032556 0.020807 ...
%     0.009461 0.001035 0.000651];

% LongS = [0.003286 0.004843 0.010526 0.006787 ...
%           0.003129 0.000323 0.000214]; % LONGN Rrs
      
% LongS = [0.002991 0.004547 0.010292 0.006556 ...
%             0.002935 0.000281 0.000176];  % LONGS Rrs

% find LongS in waterpixels with index I
[Y,I] = min(sqrt(mean((waterpixels-ones(size(waterpixels,1),1)*LongS).^2,2)));

Inputused(IMatrix(I))
DPFused(IMatrix(I))
LUTconcused(IMatrix(I),:)

LongSconc = [112.76 46 1.1953];
LongSconcret = XResults(I,:);

figure 
fs = 15;
set(gcf,'color','white')
plot(L8bands,LongS,'.-g')
hold on
plot(L8bands,waterpixels(I,:),'.-r')
plot(L8bands,LUTused(IMatrix(I),:),'.-b')
plot(L8bands,LUTused(rule6,:)','k')
title('LongS','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
legend('Field','Rrs','ret. from HL','DPF LUT')

%% Find Cranb
rule6 = strcmp(Inputused(:),'input140408LONGS')&...
    LUTconcused(:,1)==60 & LUTconcused(:,2)==25 &...
    LUTconcused(:,3)==1.00;

Cranb = [0.010979 0.017012 0.044860 ...
    0.025397 0.007961 0.001144 0.000468];

% Cranb = [ 0.003640 0.005506 0.014292 0.008114 ...
%            0.002373 0.000267 0.000190]; % in Rrs

% find LongS in waterpixels with index I
[Y,I] = min(sqrt(mean((waterpixels-ones(size(waterpixels,1),1)*Cranb).^2,2)));

Inputused(IMatrix(I))
DPFused(IMatrix(I))
LUTconcused(IMatrix(I),:)

Cranbconc = [64.08 26.7 1.0433];
Cranbconcret = XResults(I,:);

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,Cranb,'.-g')
hold on
plot(L8bands,waterpixels(I,:),'.-r')
plot(L8bands,LUTused(IMatrix(I),:),'.-b')
plot(L8bands,LUTused(rule6,:)','k')
title('Cranb','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
legend('Field','Rrs','ret. from HL','DPF LUT')

%% Find OntOS
% rule6 = strcmp(c{1}(:),'input140408ONTNS')&...
%     LUTconc(:,1)==1.0 & LUTconc(:,2)==1.0 &...
%     LUTconc(:,3)==0.21;

rule6 = strcmp(Inputused(:),'input140408ONTNS')&...
    LUTconcused(:,1)==1.0 & LUTconcused(:,2)==1.0 &...
    LUTconcused(:,3)==0.21;

% OntOS = [0.010809  0.014465  0.013691  ...
%     0.003076 -0.000265  0.000113 -0.000050];

OntOS = [ 0.007588 0.010373 0.008942 0.002004 ...
			 0.000197 0.000113 0.000057]; % ONTOS  Rrs


% find LongS in waterpixels with index I
[Y,I] = min(sqrt(mean((waterpixels-ones(size(waterpixels,1),1)*OntOS).^2,2)));

Inputused(IMatrix(I))
DPFused(IMatrix(I))
LUTconcused(IMatrix(I),:)

OntOSconc = [0.96 1.0 0.2188];
OntOSconcret = XResults(I,:);

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,OntOS,'-g')
hold on
plot(L8bands,waterpixels(I,:),'.-r')
plot(L8bands,LUTused(IMatrix(I),:),'.-b')
plot(L8bands,LUTused(rule6,:)','k')
title('OntOS','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
legend('Field','Rrs','ret. from HL','DPF LUT')


%% Find OntNS
rule6 = strcmp(Inputused,'input140408ONTNS')&...
    LUTconcused(:,1)==0.5 & LUTconcused(:,2)==2.0 &...
    LUTconcused(:,3)==0.11;

OntNS = [0.012930 0.019153 0.021296 0.004983 ...
    0.000745 0.000520 0.000194];
% 
% OntNS = [ 0.003355  0.004671  0.004345  0.000874 ...
%           0.000042 -0.000003  0.000006 ]; % in Rrs


% find LongS in waterpixels with index I
[~,I] = min(sqrt(mean((waterpixels-ones(size(waterpixels,1),1)*OntNS).^2,2)));

Inputused(IMatrix(I))
DPFused(IMatrix(I))
LUTconcused(IMatrix(I),:)

OntNSconc = [0.48 1.6 0.1152];
OntNSconcret = XResults(I,:);

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,OntNS,'-g')
hold on
plot(L8bands,waterpixels(I,:),'.-r')
plot(L8bands,LUTused(IMatrix(I),:),'.-b')
plot(L8bands,LUTused(rule6,:)','k')
title('OntNS','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
legend('Field','Rrs','ret. from HL','DPF LUT')


%% Scatter plot

fs = 25;
ms = 25; %marker size

figure
subplot(1,3,1)
set(gcf,'color','white')
plot(LongSconc(1),LongSconcret(1),'.r','MarkerSize', ms);
hold on
plot(Cranbconc(1),Cranbconcret(1),'.k','MarkerSize', ms);
plot(OntOSconc(1),OntOSconcret(1),'.b','MarkerSize', ms);
plot(OntNSconc(1),OntNSconcret(1),'.g','MarkerSize', ms);
maxconcChl = 160;
plot([0 maxconcChl],[0 maxconcChl],'k')
ylim([0 maxconcChl])
xlim([0 maxconcChl])
axis equal
title('<Chl>, [\mug/L]','fontsize',fs)
set(gca,'fontsize',fs)
xlabel('measured','fontsize',fs)
ylabel('retrieved','fontsize',fs)
legend('LongS','Cranb','OntOS','OntNS')

% figure
subplot(1,3,2)
set(gcf,'color','white')
plot(LongSconc(2),LongSconcret(2),'.r','MarkerSize', ms);
hold on
plot(Cranbconc(2),Cranbconcret(2),'.k','MarkerSize', ms);
plot(OntOSconc(2),OntOSconcret(2),'.b','MarkerSize', ms);
plot(OntNSconc(2),OntNSconcret(2),'.g','MarkerSize', ms);
maxconcTSS = 60;
plot([0 maxconcTSS],[0 maxconcTSS],'k')
ylim([0 maxconcTSS])
xlim([0 maxconcTSS])
axis equal
title('<TSS>, [mg/L]','fontsize',fs)
set(gca,'fontsize',fs)
xlabel('measured','fontsize',fs)
ylabel('retrieved','fontsize',fs)
% legend('LongS','Cranb','OntOS','OntNS')

subplot(1,3,3)
% figure
set(gcf,'color','white')
plot(LongSconc(3),LongSconcret(3),'.r','MarkerSize', ms);
hold on
plot(Cranbconc(3),Cranbconcret(3),'.k','MarkerSize', ms);
plot(OntOSconc(3),OntOSconcret(3),'.b','MarkerSize', ms);
plot(OntNSconc(3),OntNSconcret(3),'.g','MarkerSize', ms);
maxconcCDOM = 1.5;
plot([0 maxconcCDOM],[0 maxconcCDOM],'k')
ylim([0 maxconcCDOM])
xlim([0 maxconcCDOM])
axis equal
title('a_{CDOM}(440nm), [1/m]','fontsize',fs)
set(gca,'fontsize',fs)
xlabel('measured','fontsize',fs)
ylabel('retrieved','fontsize',fs)

% legend('LongS','Cranb','OntOS','OntNS')



%% Figure for Dr. John - 09/23/13 Linear Scale
fs = 30; % font size
cbfs = 15; % colorbar font size

figure
set(gcf,'color','white')
subplot(2,2,1)
imagesc(impos)
title('RGB image ','fontsize',fs)
set(gca,'fontsize',fs)
axis equal
axis image
axis off

subplot(2,2,2)
imagesc(CHLmap)
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
imagesc(CDOMmap)
% caxis([min(XResults(:,3)) max(XResults(:,3))])
title('a_{CDOM}(440nm), [1/m] ','fontsize',fs)
set(gca,'fontsize',fs)
h = colorbar;
set(h,'fontsize',cbfs)
axis equal
axis image
axis off
%% Mapping Concentrations linear scale
fs = 30; % font size
cbfs = 15; % colorbar font size

figure
set(gcf,'color','white')
imagesc(impos)
title('RGB image ','fontsize',fs)
set(gca,'fontsize',fs)
axis equal
axis image
axis off

figure
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

save('CHL.txt','-ascii','-double','-tabs','CHLmap')

figure
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

save('TSS.txt','-ascii','-double','-tabs','SMmap')

figure
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

save('CDOM.txt','-ascii','-double','-tabs','CDOMmap')
%% Plot Input (ONTNS or LONGS) and DPFs retrieved
figure
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



%% Mapping Concentrations log scale


figure
fs = 15;
set(gcf,'color','white')
imagesc(CHLmaplog10)
title('CHL map log scale','fontsize',fs)
set(gca,'fontsize',fs)
axis('equal')
colorbar

figure
fs = 15;
set(gcf,'color','white')
imagesc(SMmaplog10)
title('SM map log scale','fontsize',fs)
set(gca,'fontsize',fs)
axis('equal')
colorbar

figure
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
