dir = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/';
ncfile = 'LC80160302013262LGN00/ACOLITE/LC80160302013262LGN00_L2.nc';

ncdisp([dir ncfile])

rhow_443 = ncread([dir ncfile],'rhow_443');
rhow_443_display = (rhow_443-min(rhow_443(:)))/(max(rhow_443(:))-min(rhow_443(:)));

CTE = 1.5;
rhow_483 = ncread([dir ncfile],'rhow_483');
threshold_Bpos = nanmean(rhow_483(:))+CTE*nanstd(rhow_483(:));
rhow_483(rhow_483>threshold_Bpos)=threshold_Bpos;
threshold_Bneg = nanmean(rhow_483(:))-CTE*nanstd(rhow_483(:));
rhow_483(rhow_483<threshold_Bneg)=threshold_Bneg;
rhow_483_display = (rhow_483-min(rhow_483(:)))/(max(rhow_483(:))-min(rhow_483(:)));


rhow_561 = ncread([dir ncfile],'rhow_561');
threshold_Gpos = nanmean(rhow_561(:))+CTE*nanstd(rhow_561(:));
rhow_561(rhow_561>threshold_Gpos)=threshold_Gpos;
threshold_Gneg = nanmean(rhow_561(:))-CTE*nanstd(rhow_561(:));
rhow_561(rhow_561<threshold_Gneg)=threshold_Gneg;
rhow_561_display = (rhow_561-min(rhow_561(:)))/(max(rhow_561(:))-min(rhow_561(:)));

rhow_655 = ncread([dir ncfile],'rhow_655');
threshold_Rpos = nanmean(rhow_655(:))+CTE*nanstd(rhow_655(:));
rhow_655(rhow_655>threshold_Rpos)=threshold_Rpos;
threshold_Rneg = nanmean(rhow_655(:))-CTE*nanstd(rhow_655(:));
rhow_655(rhow_655<threshold_Rneg)=threshold_Rneg;
rhow_655_display = (rhow_655-min(rhow_655(:)))/(max(rhow_655(:))-min(rhow_655(:)));

rhow_865 = ncread([dir ncfile],'rhow_865');
rhow_865_display = (rhow_865-min(rhow_865(:)))/(max(rhow_865(:))-min(rhow_865(:)));

rhoam_865 = ncread([dir ncfile],'rhoam_865');
rhoam_865_display = (rhoam_865-min(rhoam_865(:)))/(max(rhoam_865(:))-min(rhoam_865(:)));

lon = ncread([dir ncfile],'lon');
lat = ncread([dir ncfile],'lat');

ncid = netcdf.open ([dir ncfile],'NC_NOWRITE');

RGB_display(:,:,3) = double(rhow_483_display);
RGB_display(:,:,2) = double(rhow_561_display);
RGB_display(:,:,1) = double(rhow_655_display);

figure
set(gcf,'color','white')
geoshow(lat,lon,RGB_display)
% colorbar
% axis off

figure
set(gcf,'color','white')
imagesc(rhow_443')
colorbar
axis equal
% axis off

