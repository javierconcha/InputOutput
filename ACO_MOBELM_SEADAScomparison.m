cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/
addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/gmregress/')
%%
dir = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/';


tif13262 = 'LC80160302013262LGN00/Collocated13262_SWIRFranz_L2LAC_140315bigger500.tif';

% 1: lon_R_R
% 2: lat_R_R
% 3: rhoam_865
% 4: rhow_443_R_R
% 5: rhow_483_R_R
% 6: rhow_561_R_R
% 7: rhow_655_R_R
% 8: rhow_865_R_R
% 9: aot_865_D_R
% 10: angstrom_D_R
% 11: Rrs_443_D_R
% 12: Rrs_483_D_R
% 13: Rrs_561_D_R
% 14: Rrs_655_D_R
% 15: chlor_a_D_R
% 16: kd_490_D_R
% 17: l2_flags_D_R
% 18: longitude_D_R
% 19: latitude_D_R
% 20: lat_R
% 21: lon_R
% 22: band_1_D
% 23: band_2_D
% 24: band_3_D
% 25: band_4_D
% 26: band_5_D
% 27: band_6_D
% 28: band_7_D


% tif_file = 'LC80160302013262LGN00/Collocated2013262.tif';
%% Image 2013262
filename = [dir tif13262];

im2013262 = imread(filename);

densityflag = 1;
regressiontype = 'RMA';
maxref = 0.030;
date = '2013262';

% Band 1: 443nm
Rrs_443A = double(im2013262(:,:, 4))/pi; % Acolite
Rrs_443S = (double(im2013262(:,:,12))*2.0E-6+0.05); % SeaDAS, converted with scale_factor and add_offset from band attributes
Rrs_443E = double(im2013262(:,:,22))/pi; % MoBELM; for 2013262 image, not in Rrs
denscatplot(Rrs_443A,Rrs_443E,regressiontype,densityflag,'443',maxref,date,'Acolite','MoB-ELM')
denscatplot(Rrs_443A,Rrs_443S,regressiontype,densityflag,'443',maxref,date,'Acolite','SeaDAS')
denscatplot(Rrs_443S,Rrs_443E,regressiontype,densityflag,'443',maxref,date,'SeaDAS','MoB-ELM')
%% Plot Chl
% Chl_aS = double(im2013262(:,:,15));
% figure
% image(Chl_aS)

%% Band 2: 483nm
Rrs_483A = double(im2013262(:,:, 5))/pi; % Acolite
Rrs_483S = (double(im2013262(:,:,13))*2.0E-6+0.05); % SeaDAS, converted with scale_factor and add_offset from band attributes
Rrs_483E = double(im2013262(:,:,23))/pi; % MoBELM; for 2013262 image, not in Rrs
denscatplot(Rrs_483A,Rrs_483E,regressiontype,densityflag,'483',maxref,date,'Acolite','MoB-ELM')
denscatplot(Rrs_483A,Rrs_483S,regressiontype,densityflag,'483',maxref,date,'Acolite','SeaDAS')
denscatplot(Rrs_483S,Rrs_483E,regressiontype,densityflag,'483',maxref,date,'SeaDAS','MoB-ELM')

%% Band 3: 561nm
Rrs_561A = double(im2013262(:,:, 6))/pi; % Acolite
Rrs_561S = (double(im2013262(:,:,14))*2.0E-6+0.05); % SeaDAS, converted with scale_factor and add_offset from band attributes
Rrs_561E = double(im2013262(:,:,24))/pi; % MoBELM; for 2013262 image, not in Rrs
denscatplot(Rrs_561A,Rrs_561E,regressiontype,densityflag,'561',maxref,date,'Acolite','MoB-ELM')
denscatplot(Rrs_561A,Rrs_561S,regressiontype,densityflag,'561',maxref,date,'Acolite','SeaDAS')
denscatplot(Rrs_561S,Rrs_561E,regressiontype,densityflag,'561',maxref,date,'SeaDAS','MoB-ELM')

%% Band 4: 655nm
Rrs_655A = double(im2013262(:,:, 7))/pi; % Acolite
Rrs_655S = double(im2013262(:,:,15))*2.0E-6+0.05; % SeaDAS, converted with scale_factor and add_offset from band attributes
Rrs_655E = double(im2013262(:,:,25))/pi; % MoBELM; for 2013262 image, not in Rrs
denscatplot(Rrs_655A,Rrs_655E,regressiontype,densityflag,'655',maxref,date,'Acolite','MoB-ELM')
denscatplot(Rrs_655A,Rrs_655S,regressiontype,densityflag,'655',maxref,date,'Acolite','SeaDAS')
denscatplot(Rrs_655S,Rrs_655E,regressiontype,densityflag,'655',maxref,date,'SeaDAS','MoB-ELM')

%% Band 5: 865nm
Rrs_865A = double(im2013262(:,:, 7))/pi; % Acolite
Rrs_865E = double(im2013262(:,:,26))/pi; % MoBELM; for 2013262 image, not in Rrs
denscatplot(Rrs_865A,Rrs_865E,regressiontype,densityflag,'865',maxref,date,'Acolite','MoB-ELM')