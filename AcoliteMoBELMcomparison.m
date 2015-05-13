cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/
addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/gmregress/')
%%
dir = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/';
tif_file = 'LC80160302013262LGN00/Collocated2013262.tif';
filename = [dir tif_file];

im2013262 = imread(filename);
proj = geotiffinfo(filename);
info = imfinfo(filename);

densityflag = 1;
regressiontype = 'RMA';
maxref = 0.030;
date = '2013262';
comptype = 'SWIR';
%% Band 1: 443nm
rho_443 = im2013262(:,:,4);
Rrs_443A = double(rho_443./pi);

Rrs_443E = double(im2013262(:,:,9));

denscatplot(Rrs_443A,Rrs_443E,regressiontype,densityflag,'443',maxref,date,comptype)

%% Band 2: 483nm
rho_483 = im2013262(:,:,5);
Rrs_483A = double(rho_483./pi);

Rrs_483E = double(im2013262(:,:,10));

denscatplot(Rrs_483A,Rrs_483E,regressiontype,densityflag,'483',maxref,date,comptype)

%% Band 3: 561nm
rho_561 = im2013262(:,:,6);
Rrs_561A = double(rho_561./pi);

Rrs_561E = double(im2013262(:,:,11));

denscatplot(Rrs_561A,Rrs_561E,regressiontype,densityflag,'561',maxref,date,comptype)

%% Band 4: 655nm
rho_655 = im2013262(:,:,7);
Rrs_655A = double(rho_655./pi);

Rrs_655E = double(im2013262(:,:,12));

denscatplot(Rrs_655A,Rrs_655E,regressiontype,densityflag,'655',maxref,date,comptype)

%% Band 5: 865nm
rho_865 = im2013262(:,:,8);
Rrs_865A = double(rho_865./pi);

Rrs_865E = double(im2013262(:,:,13));

denscatplot(Rrs_865A,Rrs_865E,regressiontype,densityflag,'865',maxref,date,comptype)