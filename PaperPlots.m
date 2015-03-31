cd /Users/javier/Desktop/Javier/PHD_RIT/ConferencesAndApplications/2015_Landsat_Special_Issue/Images
%% Plot ROI 
clear
clc

folderpath = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/';
% filename = 'LC80170302014272LGN00/LC80170302014272LGN00_ROI_RGB.tif';
filename = 'LC80160302013262LGN00/LC80160302013262LGN00_ONelm140629.tif';

% filename = 'LC80160302013262LGN00/LC80160302013262LGN00_WaterAndDTROC.tif';
% filename = 'LC80160302013262LGN00/LC80160302013262LGN00_WholeImage.tif';



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

% folderpath = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/';
% filename = 'LC80170302014272LGN00/LC80170302014272LGN00_ROI_RGB.tif';
% filepath = [folderpath filename];

[RGB,map] = imread(filepath);

RGB_R = double(RGB(:,:,4));
RGB_R(RGB_R<0)=0;

RGB_G = double(RGB(:,:,3));
RGB_G(RGB_G<0)=0;

RGB_B = double(RGB(:,:,2));
RGB_B(RGB_B<0)=0;

% Adjusting threshold for display
CTE = 1.5; % constant

threshold_Bpos = mean(RGB_B(:))+CTE*std(RGB_B(:));
RGB_B(RGB_B>threshold_Bpos)=threshold_Bpos;
threshold_Bneg = mean(RGB_B(:))-CTE*std(RGB_B(:));
RGB_B(RGB_B<threshold_Bneg)=threshold_Bneg;

threshold_Gpos = mean(RGB_G(:))+CTE*std(RGB_G(:));
RGB_G(RGB_G>threshold_Gpos)=threshold_Gpos;
threshold_Gneg = mean(RGB_G(:))-CTE*std(RGB_G(:));
RGB_G(RGB_G<threshold_Gneg)=threshold_Gneg;

threshold_Rpos = mean(RGB_R(:))+CTE*std(RGB_R(:));
RGB_R(RGB_R>threshold_Rpos)=threshold_Rpos;
threshold_Rneg = mean(RGB_R(:))-CTE*std(RGB_R(:));
RGB_R(RGB_R<threshold_Rneg)=threshold_Rneg;

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

sep = 0.05;
% sep = 0.5;

ax = axesm(mstruct, 'Grid', 'on',...
    'GColor', [.9 .9 .9], ...
    'MapProjection','utm',...
    'MapLatlimit', R.Latlim, 'MapLonLimit', R.Lonlim, ...
    'ParallelLabel', 'on', 'PLabelLocation', sep, 'PlabelMeridian', 'west', ...
    'MeridianLabel', 'on', 'MlabelLocation', sep, 'MLabelParallel', 'south', ...
    'MLabelRound', -2, 'PLabelRound', -2, ...
    'PLineVisible', 'on', 'PLineLocation', sep, ...
    'MLineVisible', 'on', 'MlineLocation', sep, ...
    'Grid','off', ...
    'Frame','on', ...
    'FontWeight','bold', ...
    'FontSize',24, ...
    'LabelUnits','degrees');
%
geoshow(RGBdisplay, R)

axis off
tightmap

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
% subplot(1,2,1)
fs = 15;
lw = 1.5;
set(gcf,'color','white')
plot(ELMpar(:,1),Darkrad140929,'--k','LineWidth',lw)
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
% subplot(1,2,2)
fs = 15;
lw = 1.5;
set(gcf,'color','white')
plot(ELMpar(:,1),ONTOSRef_140919,'--k','LineWidth',lw)
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
  0.443000   0.005023   0.005443   0.005369   0.005977   0.005441   0.032622;
  0.482600   0.007257   0.007529   0.007508   0.008020   0.008765   0.040927;
  0.561300   0.008571   0.014000   0.013740   0.020194   0.016108   0.058821;
  0.654600   0.001785   0.007206   0.006838   0.010907   0.007271   0.073146;
  0.864600   0.000069   0.001626   0.001545   0.002017   0.001550   0.093038;
  1.609000   0.000011   0.000062   0.000090   0.000261   0.000057   0.116034;
  2.201000  -0.000001  -0.000011   0.000023   0.000188   0.000011   0.117406];
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
plot(ROIstat(:,1),ROIstat(:,2:end-1),'LineWidth',lw)
legend('ONTOS','LONGS','LONGN','CRANB','IBAYN')
xlabel('wavelength [\mum]','fontsize',fs)
ylabel('Rrs [1/sr]','fontsize',fs)
set(gca,'fontsize',fs)
xlim([.4 2.5]) 
grid on