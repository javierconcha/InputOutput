cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/
addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/gmregress/')
%%
dir = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/';
tif_file = 'LC80160302013262LGN00/Collocated2013262.tif';
filename = [dir tif_file];

im2013262 = imread(filename);
proj = geotiffinfo(filename);
info = imfinfo(filename);

%% Band 1: 443nm
rho_443 = im2013262(:,:,4);
Rrs_443A = double(rho_443./pi);

Rrs_443E = double(im2013262(:,:,9));

denscatplot(Rrs_443A,Rrs_443E)