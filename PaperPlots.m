%% Plot ROI 
clear
clc

folderpath = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/';
filename = 'LC80170302014272LGN00/LC80170302014272LGN00_ROI_RGB.tif';

filepath = [folderpath filename];
%%
% from http://www.mathworks.com/help/map/ref/geotiff2mstruct.html?refresh=true
% Compare inverse transform of points using projinv and minvtran.
% Obtain the projection structure of 'boston.tif'.
% filename = 'boston.tif';
proj = geotiffinfo(filepath);

% Convert the corner map coordinates to latitude and longitude.
x = proj.CornerCoords.X;
y = proj.CornerCoords.Y;
[latProj, lonProj] = projinv(proj, x, y);

% Obtain the mstruct from the GeoTIFF projection.
mstruct = geotiff2mstruct(proj);

% Convert the units of x and y to meter to match projection units.
x = unitsratio('meter','sf') * x;
y = unitsratio('meter','sf') * y;

% Convert the corner map coordinates to latitude and longitude.
[latMstruct, lonMstruct] = minvtran(mstruct, x, y);

% % Verify the values are within a tolerance of each other.
% abs(latProj - latMstruct) <= 1e-7
% abs(lonProj - lonMstruct) <= 1e-7


R = georasterref('RasterInterpretation','cells');
R.RasterSize = [proj.Height proj.Width];
R.Latlim  = [min(latProj(:)) max(latProj(:))];
R.Lonlim = [min(lonProj(:)) max(lonProj(:))];
R.ColumnsStartFrom = 'north';
R.RowsStartFrom = 'west';

folderpath = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/';
filename = 'LC80170302014272LGN00/LC80170302014272LGN00_ROI_RGB.tif';
filepath = [folderpath filename];

[RGB,map] = imread(filepath);

RGB_R = double(RGB(:,:,3));
RGB_G = double(RGB(:,:,2));
RGB_B = double(RGB(:,:,1));

% Adjusting threshold for display
CTE = 1.5;

threshold_B = mean(RGB_B(:))+CTE*std(RGB_B(:));
RGB_B(RGB_B>threshold_B)=threshold_B;

threshold_G = mean(RGB_G(:))+CTE*std(RGB_G(:));
RGB_G(RGB_G>threshold_G)=threshold_G;

threshold_R = mean(RGB_R(:))+CTE*std(RGB_R(:));
RGB_R(RGB_R>threshold_R)=threshold_R;

% Convert values to [0 1] for display
RGBdisplay_B = (RGB_B-min(RGB_B(:)))/(max(RGB_B(:))-min(RGB_B(:)));
RGBdisplay_G = (RGB_G-min(RGB_G(:)))/(max(RGB_G(:))-min(RGB_G(:)));
RGBdisplay_R = (RGB_R-min(RGB_R(:)))/(max(RGB_R(:))-min(RGB_R(:)));

RGBdisplay(:,:,3) = RGBdisplay_B;
RGBdisplay(:,:,2) = RGBdisplay_G;
RGBdisplay(:,:,1) = RGBdisplay_R;

figure('Renderer', 'opengl')
fs = 20;
set(gcf,'color','white')
set(gca,'fontsize',fs)
% ax = worldmap(R.Latlim,R.Lonlim);

ax = axesm(mstruct, 'Grid', 'on',...
    'GColor', [.9 .9 .9], ...
    'MapProjection','utm',...
    'MapLatlimit', R.Latlim, 'MapLonLimit', R.Lonlim, ...
    'ParallelLabel', 'on', 'PLabelLocation', .025, 'PlabelMeridian', 'west', ...
    'MeridianLabel', 'on', 'MlabelLocation', .05, 'MLabelParallel', 'south', ...
    'MLabelRound', -2, 'PLabelRound', -2, ...
    'PLineVisible', 'on', 'PLineLocation', .025, ...
    'MLineVisible', 'on', 'MlineLocation', .05, ...
    'Grid','off', ...
    'Frame','on', ...
    'FontWeight','bold', ...
    'FontSize',14, ...
    'LabelUnits','degrees');

geoshow(RGBdisplay, R)

axis off
tightmap
%%
% Read in a gray scale image and make it double 
I1 = imread('moon.tif');
I1(1:5,1:5,1)

I2 = imread(filepath);
I2(1:5,1:5,1)
% in the range 0-255.
doubleGrayImage1 = double(I1);
subplot(2,2,1);
imshow(doubleGrayImage1, []) % Doesn't look right.
doubleGrayImage1(1:5,1:5,1)

% Read in a gray scale image and make it double 
% in the range 0-1.
doubleGrayImage2 = im2double(I1);
subplot(2,2,2);
imshow(doubleGrayImage2, []) % Doesn't look right.
doubleGrayImage2(1:5,1:5,1)

% Read in an image and make it double 
% in the range 0-255.
doubleRGBImage3 = double(I2);
subplot(2,2,3);
imshow(doubleRGBImage3, []) % Doesn't look right.
doubleRGBImage3(1:5,1:5,1)

% Read in an image and make it double 
% in the range 0-1.
doubleRGBImage4 = im2double(I2);
subplot(2,2,4);
imshow(doubleRGBImage4) % Now looks right, but values are changed.
doubleRGBImage4(1:5,1:5,1)

%% ELM parameters

% ENVI ASCII Plot File [Wed Feb 25 11:13:20 2015]
% Column 1: Wavelength
% Column 2: Darkrad140929.txt:C2~~2##255,0,0
% Column 3: ONTOSRef_140919.txt:C2~~22##0,128,0
% Column 4: PIFrad140929.txt:C2~~4##0,0,255
% Column 5: RrsDTROC4272.txt:C2~~6##0,255,255
% Column 6: RrsDTROC4272.txt:C2~~7##255,0,255
ELMpar = [... 
 0.443000  47.715229  0.009902  77.846487  0.033414  0.033414;
 0.482600  38.305432  0.009372  75.344770  0.038789  0.038789;
 0.561300  21.273056  0.007824  65.115974  0.047262  0.047262;
 0.654600  10.047512  0.003416  57.065069  0.052816  0.052816;
 0.864600   2.932883  0.001856  40.646813  0.061532  0.061532;
 1.609000   0.161773  0.000688  10.359676  0.067213  0.067213;
 2.201000   0.029518  0.000066   2.766715  0.057581  0.057581];


Wavelength      = ELMpar(:,1);
Darkrad140929 	= ELMpar(:,2);
ONTOSRef_140919 = ELMpar(:,3);
PIFrad140929   	= ELMpar(:,4);
RrsDTROC4272 	= ELMpar(:,5);

figure
fs = 15;
lw = 1.5;
set(gcf,'color','white')
plot(ELMpar(:,1),Darkrad140929,'k','LineWidth',lw)
hold on 
plot(ELMpar(:,1),PIFrad140929 ,'LineWidth',lw)
legend('Dark: water ROI','Bright: PIF from L8')
title('Radiance values for ELM-based method','fontsize',fs)
xlabel('wavelength [\mum]','fontsize',fs)
ylabel('Radiance [W/m^2/sr/\mum]','fontsize',fs)
set(gca,'fontsize',fs)
xlim([.4 2.5])
grid on

figure
fs = 15;
lw = 1.5;
set(gcf,'color','white')
plot(ELMpar(:,1),ONTOSRef_140919,'k','LineWidth',lw)
hold on 
plot(ELMpar(:,1),RrsDTROC4272,'LineWidth',lw)
legend('Dark: ONTOS field','Bright: PIF L8 refl. product')
title('Reflectance values for ELM-based method','fontsize',fs)
xlabel('wavelength [\mum]','fontsize',fs)
ylabel('Rrs [1/\pi]','fontsize',fs)
set(gca,'fontsize',fs)
xlim([.4 2.5])
grid on

%% ROIs 
ROIstat = [...
  0.443000  0.009901  0.010249  0.010188  0.010691  0.010247  0.032758;
  0.482600  0.009343  0.009598  0.009578  0.010056  0.010751  0.040786;
  0.561300  0.007739  0.013285  0.013020  0.019612  0.015439  0.059070;
  0.654600  0.003437  0.008682  0.008326  0.012264  0.008745  0.072488;
  0.864600  0.001848  0.003361  0.003281  0.003740  0.003287  0.092126;
  1.609000  0.000699  0.000749  0.000777  0.000946  0.000744  0.115534;
  2.201000  0.000065  0.000056  0.000089  0.000254  0.000077  0.117337];

wl    = ROIstat(:,1);
ONTOS = ROIstat(:,2);
LONGS = ROIstat(:,3);
LONGN = ROIstat(:,4);
CRANB = ROIstat(:,5);
IBAYN = ROIstat(:,6);
Sand  = ROIstat(:,7);


figure
fs = 15;
lw = 1.5;
set(gcf,'color','white')
plot(ROIstat(:,1),ROIstat(:,2:end),'LineWidth',lw)
legend('ONTOS','LONGS','LONGN','CRANB','IBAYN','Sand')
xlabel('wavelength [\mum]','fontsize',fs)
ylabel('Rrs [1/\pi]','fontsize',fs)
set(gca,'fontsize',fs)
xlim([.4 2.5]) 
grid on