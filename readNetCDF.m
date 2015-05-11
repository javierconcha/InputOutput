cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/

%%
dir = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/';
ncfile = 'LC80160302013262LGN00/ACOLITE/LC80160302013262LGN00_L2_SWIR_Unity.nc';
% ncfile = 'LC80160302013262LGN00/ACOLITE/LC80160302013262LGN00_L2_NIR.nc';

%%
ncdisp([dir ncfile])

rhow_443 = ncread([dir ncfile],'rhow_443');
rhow_443 = flipud(rhow_443');
rhow_443_display = (rhow_443-min(rhow_443(:)))/(max(rhow_443(:))-min(rhow_443(:)));

CTE = 1.5;
rhow_483 = ncread([dir ncfile],'rhow_483');
rhow_483 = flipud(rhow_483');
threshold_Bpos = nanmean(rhow_483(:))+CTE*nanstd(rhow_483(:));
rhow_483(rhow_483>threshold_Bpos)=threshold_Bpos;
threshold_Bneg = nanmean(rhow_483(:))-CTE*nanstd(rhow_483(:));
rhow_483(rhow_483<threshold_Bneg)=threshold_Bneg;
rhow_483_display = (rhow_483-min(rhow_483(:)))/(max(rhow_483(:))-min(rhow_483(:)));


rhow_561 = ncread([dir ncfile],'rhow_561');
rhow_561 = flipud(rhow_561');
threshold_Gpos = nanmean(rhow_561(:))+CTE*nanstd(rhow_561(:));
rhow_561(rhow_561>threshold_Gpos)=threshold_Gpos;
threshold_Gneg = nanmean(rhow_561(:))-CTE*nanstd(rhow_561(:));
rhow_561(rhow_561<threshold_Gneg)=threshold_Gneg;
rhow_561_display = (rhow_561-min(rhow_561(:)))/(max(rhow_561(:))-min(rhow_561(:)));

rhow_655 = ncread([dir ncfile],'rhow_655');
rhow_655 = flipud(rhow_655');
threshold_Rpos = nanmean(rhow_655(:))+CTE*nanstd(rhow_655(:));
rhow_655(rhow_655>threshold_Rpos)=threshold_Rpos;
threshold_Rneg = nanmean(rhow_655(:))-CTE*nanstd(rhow_655(:));
rhow_655(rhow_655<threshold_Rneg)=threshold_Rneg;
rhow_655_display = (rhow_655-min(rhow_655(:)))/(max(rhow_655(:))-min(rhow_655(:)));

rhow_865 = ncread([dir ncfile],'rhow_865');
rhow_865 = flipud(rhow_865');
rhow_865_display = (rhow_865-min(rhow_865(:)))/(max(rhow_865(:))-min(rhow_865(:)));

rhoam_865 = ncread([dir ncfile],'rhoam_865');
rhoam_865 = flipud(rhoam_865)';
rhoam_865_display = (rhoam_865-min(rhoam_865(:)))/(max(rhoam_865(:))-min(rhoam_865(:)));

lon = ncread([dir ncfile],'lon');
lon = flipud(lon');
lat = ncread([dir ncfile],'lat');
lat = flipud(lat');

% ncid = netcdf.open ([dir ncfile],'NC_NOWRITE');

RGB_display(:,:,3) = double(rhow_483_display);
RGB_display(:,:,2) = double(rhow_561_display);
RGB_display(:,:,1) = double(rhow_655_display);



figure
set(gcf,'color','white')
h = geoshow(lat,lon,RGB_display,'DisplayType','texturemap');
% set(h,'ZData',zeros(size(rhow_443)))
% colorbar
% axis off
%%
R = georasterref('RasterSize',size(rhow_443),...
    'Latlim',[double(min(lat(:))) double(max(lat(:)))],...
    'Lonlim',[double(min(lon(:))) double(max(lon(:)))]);
% figure
% set(gcf,'color','white')
% worldmap(double(rhow_443),R)
% geoshow(double(rhow_443),R,'DisplayType','surface',...
%     'ZData',zeros(size(rhow_443)),'CData',double(rhow_443))
% % demcmap(double(rhow_443))

figure
set(gcf,'color','white')

worldmap((RGB_display),R)
geoshow(RGB_display,R,'DisplayType','texturemap',...
    'ZData',zeros(size(rhow_443)),'CData',RGB_display)

%%
figure
set(gcf,'color','white')
imagesc(flipud(rhow_483_display))
colorbar
axis equal
% axis off



%% histogram

figure('name','NIR')
fs = 15;
set(gcf,'color','white')
subplot(2,3,1)
hist(rhow_443(:)/pi,50)
xlabel('R_{rs} 443','FontSize',fs)

subplot(2,3,2)
hist(rhow_483(:)/pi,50)
xlabel('R_{rs} 483','FontSize',fs)

subplot(2,3,3)
hist(rhow_561(:)/pi,50)
xlabel('R_{rs} 561','FontSize',fs)

subplot(2,3,4)
hist(rhow_655(:)/pi,50)
xlabel('R_{rs} 655','FontSize',fs)

subplot(2,3,5)
hist(rhow_865(:)/pi,50)
xlabel('R_{rs} 865','FontSize',fs)

subplot(2,3,6)
hist(rhoam_865(:),50)
xlabel('\rho_{am} 865','FontSize',fs)

set(gca,'fontsize',fs)
%% Negative values
bandname = '443'; 

rhow_443_zero = zeros(size(rhow_443));

rhow_443_zero(rhow_443<0) = 1;
rhow_443_zero(rhow_443>=0) = 0;

figure
set(gcf,'color','white')
imagesc(flipud(rhow_443_zero))
colormap('gray')
axis equal
axis off
title('rhow\_443 < 0','FontSize',fs)


%% Example
korea = load('korea');
Z = korea.map;
Z(80:100,:) = NaN;
R = georasterref('RasterSize',size(Z),...
    'Latlim',[30 45], 'Lonlim', [115 135]);
figure;
set(gcf,'color','white')
worldmap(Z,R)
geoshow(Z,R,'DisplayType','surface','ZData',zeros(size(Z)),'CData',Z)
demcmap(Z)

