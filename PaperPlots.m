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

%% Retrieved vs Measured
% CHL
figure
fs = 20;
ms = 10;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(LongSconc130919(1),LongSconc130919ret(1),'*k','MarkerSize', ms);
hold on
plot(Cranbconc130919(1),Cranbconc130919ret(1),'*k','MarkerSize', ms);
plot(OntOSconc130919(1),OntOSconc130919ret(1),'*k','MarkerSize', ms);
plot(OntNSconc130919(1),OntNSconc130919ret(1),'*k','MarkerSize', ms);
plot(LongSconc140929(1),LongSconc140929ret(1),'^k','MarkerSize', ms);
plot(LongNconc140929(1),LongNconc140929ret(1),'^k','MarkerSize', ms);
plot(Cranbconc140929(1),Cranbconc140929ret(1),'^k','MarkerSize', ms);
plot(IBayNconc140929(1),IBayNconc140929ret(1),'^k','MarkerSize', ms);
plot(OntOSconc140929(1),OntOSconc140929ret(1),'^k','MarkerSize', ms);
maxconcChl = 150;
plot([0 maxconcChl],[0 maxconcChl],'--k')
axis equal
ylim([0 maxconcChl])
xlim([0 maxconcChl])
xlabel('measured C_a [mg m^{-3}] ','fontsize',fs,'Position',[80 -15])
ylabel('L8 retrieved C_a [mg m^{-3}]','fontsize',fs)
% legend('LONGS','LONGN','CRANB','IBAYN','ONTOS','Location','best')
set(gca,'OuterPosition',[0 0.05 1 1])

% TSS
figure
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(LongSconc130919(2),LongSconc130919ret(2),'*k','MarkerSize', ms);
hold on
plot(Cranbconc130919(2),Cranbconc130919ret(2),'*k','MarkerSize', ms);
plot(OntOSconc130919(2),OntOSconc130919ret(2),'*k','MarkerSize', ms);
plot(OntNSconc130919(2),OntNSconc130919ret(2),'*k','MarkerSize', ms);
plot(LongSconc140929(2),LongSconc140929ret(2),'^k','MarkerSize', ms);
plot(LongNconc140929(2),LongNconc140929ret(2),'^k','MarkerSize', ms);
plot(Cranbconc140929(2),Cranbconc140929ret(2),'^k','MarkerSize', ms);
plot(IBayNconc140929(2),IBayNconc140929ret(2),'^k','MarkerSize', ms);
plot(OntOSconc140929(2),OntOSconc140929ret(2),'^k','MarkerSize', ms);
maxconcTSS = 60;
plot([0 maxconcTSS],[0 maxconcTSS],'--k')
axis equal
ylim([0 maxconcTSS])
xlim([0 maxconcTSS])
xlabel('measured TSS [g m^{-3}] ','fontsize',fs)
ylabel('L8 retrieved TSS [g m^{-3}]','fontsize',fs)
set(gca,'OuterPosition',[0 0.05 1 1])
%% CDOM
figure
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(LongSconc130919(3),LongSconc130919ret(3),'*k','MarkerSize', ms);
hold on
plot(LongSconc140929(3),LongSconc140929ret(3),'^k','MarkerSize', ms); % not correct position!!!
legend('09-19-13','09-29-14','Location','best')
plot(Cranbconc130919(3),Cranbconc130919ret(3),'*k','MarkerSize', ms);
plot(OntOSconc130919(3),OntOSconc130919ret(3),'*k','MarkerSize', ms);
plot(OntNSconc130919(3),OntNSconc130919ret(3),'*k','MarkerSize', ms);
plot(LongSconc140929(3),LongSconc140929ret(3),'^k','MarkerSize', ms);
plot(LongNconc140929(3),LongNconc140929ret(3),'^k','MarkerSize', ms);
plot(Cranbconc140929(3),Cranbconc140929ret(3),'^k','MarkerSize', ms);
plot(IBayNconc140929(3),IBayNconc140929ret(3),'^k','MarkerSize', ms);
plot(OntOSconc140929(3),OntOSconc140929ret(3),'^k','MarkerSize', ms);
maxconcCDOM = 1.5;
plot([0 maxconcCDOM],[0 maxconcCDOM],'--k')
axis equal
ylim([0 maxconcCDOM])
xlim([0 maxconcCDOM])
xlabel('measured a_{CDOM}(440nm) [1/m]','fontsize',fs,'Position',[0.8 -0.16])
ylabel('retrieved a_{CDOM}(440nm) [1/m]','fontsize',fs)
set(gca,'OuterPosition',[0 0.05 1 1])

%%
CHL_data = [...
LongSconc130919(1),LongSconc130919ret(1);
Cranbconc130919(1),Cranbconc130919ret(1);
OntOSconc130919(1),OntOSconc130919ret(1);
OntNSconc130919(1),OntNSconc130919ret(1);
LongSconc140929(1),LongSconc140929ret(1);
LongNconc140929(1),LongNconc140929ret(1);
Cranbconc140929(1),Cranbconc140929ret(1);
IBayNconc140929(1),IBayNconc140929ret(1);
OntOSconc140929(1),OntOSconc140929ret(1)];

CHL_RMSE = 100*sqrt(mean((CHL_data(:,1)-CHL_data(:,2)).^2))/(max(CHL_data(:))-min(CHL_data(:)));

CHL_std = 100*std(abs(CHL_data(:,1)-CHL_data(:,2)))/max(CHL_data(:));

TSS_data = [...
LongSconc130919(2),LongSconc130919ret(2);
Cranbconc130919(2),Cranbconc130919ret(2);
OntOSconc130919(2),OntOSconc130919ret(2);
OntNSconc130919(2),OntNSconc130919ret(2);
LongSconc140929(2),LongSconc140929ret(2);
LongNconc140929(2),LongNconc140929ret(2);
Cranbconc140929(2),Cranbconc140929ret(2);
IBayNconc140929(2),IBayNconc140929ret(2);
OntOSconc140929(2),OntOSconc140929ret(2)];

TSS_RMSE = 100*sqrt(mean((TSS_data(:,1)-TSS_data(:,2)).^2))/(max(TSS_data(:))-min(TSS_data(:)));

TSS_std = 100*std(abs(TSS_data(:,1)-TSS_data(:,2)))/max(TSS_data(:));

CDO_data = [...
LongSconc130919(3),LongSconc130919ret(3);
Cranbconc130919(3),Cranbconc130919ret(3);
OntOSconc130919(3),OntOSconc130919ret(3);
OntNSconc130919(3),OntNSconc130919ret(3);
LongSconc140929(3),LongSconc140929ret(3);
LongNconc140929(3),LongNconc140929ret(3);
Cranbconc140929(3),Cranbconc140929ret(3);
IBayNconc140929(3),IBayNconc140929ret(3);
OntOSconc140929(3),OntOSconc140929ret(3)];

CDO_RMSE = 100*sqrt(mean((CDO_data(:,1)-CDO_data(:,2)).^2))/(max(CDO_data(:))-min(CDO_data(:)));

CDO_std = 100*std(abs(CDO_data(:,1)-CDO_data(:,2)))/max(CDO_data(:));
%


error = [CHL_RMSE    TSS_RMSE    CDO_RMSE];

figure
set(gcf,'color','white')
fs = 20;
bar(error,0.5)
% hold on
% errorbar(error,[CHL_std TSS_std CDO_std],'kx')

Labels = {'Chl-a','TSS','CDOM'};
set(gca, 'XTick', 1:size(Labels,2), 'XTickLabel',Labels,'FontSize',fs);
ylabel('Percentage of RMSE [%]','FontSize',fs)
barmap=[0.7 0.7 0.7];
colormap(barmap)
ylim([0 15])
grid on