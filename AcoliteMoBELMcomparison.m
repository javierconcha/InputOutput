cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/
addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/gmregress/')
%%
dir = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/';


 tif13262_NF = 'LC80160302013262LGN00/SEADAS/Collocated2013262_140317bigger500_NIR_FranzAve.tif';
 tif13262_NP = 'LC80160302013262LGN00/SEADAS/Collocated2013262_140317bigger500_NIR_PahlevanL.tif';
 tif13262_NU = 'LC80160302013262LGN00/SEADAS/Collocated2013262_140317bigger500_NIR_Unity.tif';
 tif13262_SF = 'LC80160302013262LGN00/SEADAS/Collocated2013262_140317bigger500_SWIR_FranzAve.tif';
 tif13262_SP = 'LC80160302013262LGN00/SEADAS/Collocated2013262_140317bigger500_SWIR_PahlevalL.tif';
 tif13262_SU = 'LC80160302013262LGN00/SEADAS/Collocated2013262_140317bigger500_SWIR_Unity.tif';
 tif14272_NF = 'LC80170302014272LGN00/SEADAS/Collocated2014272_150418CRANB_ONTOS_NIR_Franzave.tif';
 tif14272_NP = 'LC80170302014272LGN00/SEADAS/Collocated2014272_150418CRANB_ONTOS_NIR_PahlevanL.tif';
 tif14272_NU = 'LC80170302014272LGN00/SEADAS/Collocated2014272_150418CRANB_ONTOS_NIR_Unity.tif';
 tif14272_SF = 'LC80170302014272LGN00/SEADAS/Collocated2014272_150418CRANB_ONTOS_SWIR_FranzAve.tif';
 tif14272_SP = 'LC80170302014272LGN00/SEADAS/Collocated2014272_150418CRANB_ONTOS_SWIR_PahlevanL.tif';
 tif14272_SU = 'LC80170302014272LGN00/SEADAS/Collocated2014272_150418CRANB_ONTOS_SWIR_Unity.tif';

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

densityflag = 0;
regressiontype = 'RMA';
maxref = 0.030;
date = '2013262';
% comptype = 'SWIR_Unity';

% Band 1: 443nm
Rrs_443A = double(im2013262(:,:,11))/pi;
Rrs_443E = double(im2013262(:,:,1))/pi; % for 2013262 image, not in Rrs
denscatplot(Rrs_443A,Rrs_443E,regressiontype,densityflag,'443',maxref,date,'Acolite','MOB-ELM')

% Band 2: 483nm
Rrs_483A = double(im2013262(:,:,12))/pi;
Rrs_483E = double(im2013262(:,:,2))/pi; % for 2013262 image, not in Rrs
denscatplot(Rrs_483A,Rrs_483E,regressiontype,densityflag,'483',maxref,date,'Acolite','MOB-ELM')

% Band 3: 561nm
Rrs_561A = double(im2013262(:,:,13))/pi;
Rrs_561E = double(im2013262(:,:,3))/pi; % for 2013262 image, not in Rrs
denscatplot(Rrs_561A,Rrs_561E,regressiontype,densityflag,'561',maxref,date,'Acolite','MOB-ELM')

% Band 4: 655nm
Rrs_655A = double(im2013262(:,:,14))/pi;
Rrs_655E = double(im2013262(:,:,4))/pi; % for 2013262 image, not in Rrs
denscatplot(Rrs_655A,Rrs_655E,regressiontype,densityflag,'655',maxref,date,'Acolite','MOB-ELM')

% Band 5: 865nm
Rrs_865A = double(im2013262(:,:,15))/pi;
Rrs_865E = double(im2013262(:,:,5))/pi; % for 2013262 image, not in Rrs
denscatplot(Rrs_865A,Rrs_865E,regressiontype,densityflag,'865',maxref,date,'Acolite','MOB-ELM')

%% Image 2014272 SWIR_Unity
filename = [dir tif14272_SU];

im2014272 = imread(filename);

densityflag = 0; % 0 or 1
regressiontype = 'RMA'; % RMA or OLS
maxref = 0.030;
date = '2014272';
% comptype = 'SWIR_Unity'

% Band 1: 443nm
Rrs_443A_SU = double(im2014272(:,:,11)/pi);
Rrs_443E_SU = double(im2014272(:,:,1));
denscatplot(Rrs_443A_SU,Rrs_443E_SU,regressiontype,densityflag,'443',maxref,date,'Acolite','MOB-ELM')

% Band 2: 483nm
Rrs_483A_SU = double(im2014272(:,:,12))/pi;
Rrs_483E_SU = double(im2014272(:,:,2));
denscatplot(Rrs_483A_SU,Rrs_483E_SU,regressiontype,densityflag,'483',maxref,date,'Acolite','MOB-ELM')

% Band 3: 561nm
Rrs_561A_SU = double(im2014272(:,:,13))/pi;
Rrs_561E_SU = double(im2014272(:,:,3));
denscatplot(Rrs_561A_SU,Rrs_561E_SU,regressiontype,densityflag,'561',maxref,date,'Acolite','MOB-ELM')

% Band 4: 655nm
Rrs_655A_SU = double(im2014272(:,:,14))/pi;
Rrs_655E_SU = double(im2014272(:,:,4));
denscatplot(Rrs_655A_SU,Rrs_655E_SU,regressiontype,densityflag,'655',maxref,date,'Acolite','MOB-ELM')

% Band 5: 865nm
Rrs_865A_SU = double(im2014272(:,:,15))/pi;
Rrs_865E_SU = double(im2014272(:,:,5));
denscatplot(Rrs_865A_SU,Rrs_865E_SU,regressiontype,densityflag,'865',maxref,date,'Acolite','MOB-ELM')

%% Image 2014272 SWIR_Pahlevan
filename = [dir tif14272_SP];

im2014272 = imread(filename);

densityflag = 0; % 0 or 1
regressiontype = 'RMA'; % RMA or OLS
maxref = 0.030;
date = '2014272';
% comptype = 'SWIR_Pahlevan'

% Band 1: 443nm
Rrs_443A_SP = double(im2014272(:,:,11)/pi);
Rrs_443E_SP = double(im2014272(:,:,1));
denscatplot(Rrs_443A_SP,Rrs_443E_SP,regressiontype,densityflag,'443',maxref,date,'Acolite','MOB-ELM')

% Band 2: 483nm
Rrs_483A_SP = double(im2014272(:,:,12))/pi;
Rrs_483E_SP = double(im2014272(:,:,2));
denscatplot(Rrs_483A_SP,Rrs_483E_SP,regressiontype,densityflag,'483',maxref,date,'Acolite','MOB-ELM')

% Band 3: 561nm
Rrs_561A_SP = double(im2014272(:,:,13))/pi;
Rrs_561E_SP = double(im2014272(:,:,3));
denscatplot(Rrs_561A_SP,Rrs_561E_SP,regressiontype,densityflag,'561',maxref,date,'Acolite','MOB-ELM')

% Band 4: 655nm
Rrs_655A_SP = double(im2014272(:,:,14))/pi;
Rrs_655E_SP = double(im2014272(:,:,4));
denscatplot(Rrs_655A_SP,Rrs_655E_SP,regressiontype,densityflag,'655',maxref,date,'Acolite','MOB-ELM')

% Band 5: 865nm
Rrs_865A_SP = double(im2014272(:,:,15))/pi;
Rrs_865E_SP = double(im2014272(:,:,5));
denscatplot(Rrs_865A_SP,Rrs_865E_SP,regressiontype,densityflag,'865',maxref,date,'Acolite','MOB-ELM')

%% Image 2014272 SWIR_Franz
filename = [dir tif14272_SF];

im2014272 = imread(filename);
proj = geotiffinfo(filename);
info = imfinfo(filename);

densityflag = 0; % 0 or 1
regressiontype = 'RMA'; % RMA or OLS
maxref = 0.030;
date = '2014272';
comptype = 'SWIR_Franz'

% Band 1: 443nm
Rrs_443A_SF = double(im2014272(:,:,11)/pi);
Rrs_443E_SF = double(im2014272(:,:,1));
denscatplot(Rrs_443A_SF,Rrs_443E_SF,regressiontype,densityflag,'443',maxref,date,'Acolite','MOB-ELM')

% Band 2: 483nm
Rrs_483A_SF = double(im2014272(:,:,12))/pi;
Rrs_483E_SF = double(im2014272(:,:,2));
denscatplot(Rrs_483A_SF,Rrs_483E_SF,regressiontype,densityflag,'483',maxref,date,'Acolite','MOB-ELM')

% Band 3: 561nm
Rrs_561A_SF = double(im2014272(:,:,13))/pi;
Rrs_561E_SF = double(im2014272(:,:,3));
denscatplot(Rrs_561A_SF,Rrs_561E_SF,regressiontype,densityflag,'561',maxref,date,'Acolite','MOB-ELM')

% Band 4: 655nm
Rrs_655A_SF = double(im2014272(:,:,14))/pi;
Rrs_655E_SF = double(im2014272(:,:,4));
denscatplot(Rrs_655A_SF,Rrs_655E_SF,regressiontype,densityflag,'655',maxref,date,'Acolite','MOB-ELM')

% Band 5: 865nm
Rrs_865A_SF = double(im2014272(:,:,15))/pi;
Rrs_865E_SF = double(im2014272(:,:,5));
denscatplot(Rrs_865A_SF,Rrs_865E_SF,regressiontype,densityflag,'865',maxref,date,'Acolite','MOB-ELM')

%% Image 2014272 SWIR_Unity VS SWIR_PahlevanL
densityflag = 0; % 0 or 1
regressiontype = 'RMA'; % RMA or OLS
denscatplot(Rrs_443A_SU,Rrs_443A_SF,regressiontype,densityflag,'443',maxref,date,'Acolite','MOB-ELM')