cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/
addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/gmregress/')
%%
dir = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/';


 tif13262_NF = 'LC80160302013262LGN00/Collocated2013262_140317bigger500_NIR_FranzAve.tif';
 tif13262_NP = 'LC80160302013262LGN00/Collocated2013262_140317bigger500_NIR_PahlevanL.tif';
 tif13262_NU = 'LC80160302013262LGN00/Collocated2013262_140317bigger500_NIR_Unity.tif';
 tif13262_SF = 'LC80160302013262LGN00/Collocated2013262_140317bigger500_SWIR_FranzAve.tif';
 tif13262_SP = 'LC80160302013262LGN00/Collocated2013262_140317bigger500_SWIR_PahlevalL.tif';
 tif13262_SU = 'LC80160302013262LGN00/Collocated2013262_140317bigger500_SWIR_Unity.tif';
 tif14272_NF = 'LC80170302014272LGN00/Collocated2014272_150418CRANB_ONTOS_NIR_Franzave.tif';
 tif14272_NP = 'LC80170302014272LGN00/Collocated2014272_150418CRANB_ONTOS_NIR_PahlevanL.tif';
 tif14272_NU = 'LC80170302014272LGN00/Collocated2014272_150418CRANB_ONTOS_NIR_Unity.tif';
 tif14272_SF = 'LC80170302014272LGN00/Collocated2014272_150418CRANB_ONTOS_SWIR_FranzAve.tif';
 tif14272_SP = 'LC80170302014272LGN00/Collocated2014272_150418CRANB_ONTOS_SWIR_PahlevanL.tif';
 tif14272_SU = 'LC80170302014272LGN00/Collocated2014272_150418CRANB_ONTOS_SWIR_Unity.tif';

% 1: band_1_R
% 2: band_2_R
% 3: band_3_R
% 4: band_4_R
% 5: band_5_R
% 6: band_6_R
% 7: band_7_R
% 8: lon_D
% 9: lat_D
% 10: rhoam_865_D
% 11: rho_443_D
% 12: rho_483_D
% 13: rho_561_D
% 14: rho_655_D
% 15: rho_865_D

% tif_file = 'LC80160302013262LGN00/Collocated2013262.tif';
%% Image 2013262
filename = [dir tif13262_SU];

im2013262 = imread(filename);
proj = geotiffinfo(filename);
info = imfinfo(filename);

densityflag = 0;
regressiontype = 'RMA';
maxref = 0.030;
date = '2013262';
comptype = 'SWIR_Unity';

% Band 1: 443nm
Rrs_443A = double(im2013262(:,:,11))/pi;
Rrs_443E = double(im2013262(:,:,1))/pi; % for 2013262 image, not in Rrs
denscatplot(Rrs_443A,Rrs_443E,regressiontype,densityflag,'443',maxref,date,comptype)

% Band 2: 483nm
Rrs_483A = double(im2013262(:,:,12))/pi;
Rrs_483E = double(im2013262(:,:,2))/pi; % for 2013262 image, not in Rrs
denscatplot(Rrs_483A,Rrs_483E,regressiontype,densityflag,'483',maxref,date,comptype)

% Band 3: 561nm
Rrs_561A = double(im2013262(:,:,13))/pi;
Rrs_561E = double(im2013262(:,:,3))/pi; % for 2013262 image, not in Rrs
denscatplot(Rrs_561A,Rrs_561E,regressiontype,densityflag,'561',maxref,date,comptype)

% Band 4: 655nm
Rrs_655A = double(im2013262(:,:,14))/pi;
Rrs_655E = double(im2013262(:,:,4))/pi; % for 2013262 image, not in Rrs
denscatplot(Rrs_655A,Rrs_655E,regressiontype,densityflag,'655',maxref,date,comptype)

% Band 5: 865nm
Rrs_865A = double(im2013262(:,:,15))/pi;
Rrs_865E = double(im2013262(:,:,5))/pi; % for 2013262 image, not in Rrs
denscatplot(Rrs_865A,Rrs_865E,regressiontype,densityflag,'865',maxref,date,comptype)

%% Image 2014272
filename = [dir tif14272_SU];

im2014272 = imread(filename);
proj = geotiffinfo(filename);
info = imfinfo(filename);

densityflag = 1; % 0 or 1
regressiontype = 'RMA'; % RMA or OLS
maxref = 0.030;
date = '2014272';
comptype = 'SWIR_Unity'

% Band 1: 443nm
Rrs_443A = double(im2014272(:,:,11)/pi);
Rrs_443E = double(im2014272(:,:,1));
denscatplot(Rrs_443A,Rrs_443E,regressiontype,densityflag,'443',maxref,date,comptype)

% Band 2: 483nm
Rrs_483A = double(im2014272(:,:,12))/pi;
Rrs_483E = double(im2014272(:,:,2));
denscatplot(Rrs_483A,Rrs_483E,regressiontype,densityflag,'483',maxref,date,comptype)

% Band 3: 561nm
Rrs_561A = double(im2014272(:,:,13))/pi;
Rrs_561E = double(im2014272(:,:,3));
denscatplot(Rrs_561A,Rrs_561E,regressiontype,densityflag,'561',maxref,date,comptype)

% Band 4: 655nm
Rrs_655A = double(im2014272(:,:,14))/pi;
Rrs_655E = double(im2014272(:,:,4));
denscatplot(Rrs_655A,Rrs_655E,regressiontype,densityflag,'655',maxref,date,comptype)

% Band 5: 865nm
Rrs_865A = double(im2014272(:,:,15))/pi;
Rrs_865E = double(im2014272(:,:,5));
denscatplot(Rrs_865A,Rrs_865E,regressiontype,densityflag,'865',maxref,date,comptype)