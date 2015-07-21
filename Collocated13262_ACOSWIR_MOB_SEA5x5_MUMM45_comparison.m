cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/
addpath('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/gmregress/')
%%
dir = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/';
tif13262 = 'LC80160302013262LGN00/SEADAS/Collocated13262_ACOSWIR_MOB_SEA5x5_MUMM45.tif';
% Image 2013262
filename = [dir tif13262];

im2013262 = imread(filename);

info = geotiffinfo(filename);
%1 lon_R_R_R_R_R
%2 lat_R_R_R_R_R
%3 rhoam_865_R_R_R_R_R
%4 rhow_443_ACO_R
%5 rhow_483_R_R_R_R_R
%6 rhow_561_R_R_R_R_R
%7 rhow_655_R_R_R_R_R
%8 rhow_865_R_R_R_R_R
%9 band_1_D_R_R_R_R
%10 band_2_D_R_R_R_R
%11 band_3_D_R_R_R_R
%12 band_4_D_R_R_R_R
%13 band_5_D_R_R_R_R
%14 band_6_D_R_R_R_R
%15 band_7_D_R_R_R_R
%16 lon_R_R_R_R
%17 lat_R_R_R_R
%18 Rrs_443_D_R_R_R
%19 chlor_a_D_R_R_R
%20 l2_flags_D_R_R_R
%21 Rrs_482_D_R_R_R
%22 Rrs_561_D_R_R_R
%23 Rrs_655_D_R_R_R
%24 longitude_R_R_R_R
%25 latitude_R_R_R_R
%26 lon_R_R_R
%27 lat_R_R_R
%28 chlor_a_MOB_D_R_R
CHL_MOB = double(im2013262(:,:,28)); % MOB-ELM
%29 lon_R_R
%30 lat_R_R
%31 Rrs_443_ACO_R_R
Rrs_443A = double(im2013262(:,:,31)); % Acolite
%32 Rrs_482_ACO_R_R
Rrs_483A = double(im2013262(:,:,32)); % Acolite
%33 Rrs_561_ACO_R_R
Rrs_561A = double(im2013262(:,:,33)); % Acolite
%34 Rrs_655_ACO_R_R
Rrs_655A = double(im2013262(:,:,34)); % Acolite
%35 Rrs_865_ACO_R_R
%36 Rrs_443_SEA5x5_R
Rrs_443S = double(im2013262(:,:,36))*2.0E-6+0.05; % SeaDAS, converted with scale_factor and add_offset from band attributes
%37 chlor_a_SEA5x5_R
CHL_SEA = double(im2013262(:,:,37)); % SeaDAS
%38 l2_flags_D_R
%39 Rrs_482_SEA5x5_R
Rrs_483S = double(im2013262(:,:,39))*2.0E-6+0.05; % SeaDAS, converted with scale_factor and add_offset from band attributes
%40 Rrs_561_SEA5x5_R
Rrs_561S = double(im2013262(:,:,40))*2.0E-6+0.05; % SeaDAS, converted with scale_factor and add_offset from band attributes
%41 Rrs_655_SEA5x5_R
Rrs_655S = double(im2013262(:,:,41))*2.0E-6+0.05; % SeaDAS, converted with scale_factor and add_offset from band attributes
%42 longitude_D_R
%43 latitude_D_R
%44 lon_R
%45 lat_R
%46 Rrs_443_MUMM45_R
Rrs_443M = double(im2013262(:,:,46))*2.0E-6+0.05; % MUMM (Ruddick);
%47 chlor_a_MUMM45_R
CHL_MUM = double(im2013262(:,:,47)); % MUMM
%48 l2_flags
%49 Rrs_482_MUMM45_R
Rrs_483M = double(im2013262(:,:,49))*2.0E-6+0.05; % MUMM (Ruddick);
%50 Rrs_561_MUMM45_R
Rrs_561M = double(im2013262(:,:,50))*2.0E-6+0.05; % MUMM (Ruddick);
%51 Rrs_655_MUMM45_R
Rrs_655M = double(im2013262(:,:,51))*2.0E-6+0.05; % MUMM (Ruddick);
%52 longitude
%53 latitude
%54 lat
%55 lon 
%56 Rrs_443_MOB
Rrs_443E = double(im2013262(:,:,56)); % MoB-ELM (E:ELM); for 2013262 image, not in Rrs 
%57 Rrs_482_MOB
Rrs_483E = double(im2013262(:,:,57)); % MoB-ELM (E:ELM); for 2013262 image, not in Rrs
%58 Rrs_561_MOB
Rrs_561E = double(im2013262(:,:,58)); % MoB-ELM (E:ELM); for 2013262 image, not in Rrs
%59 Rrs_655_MOB
Rrs_655E = double(im2013262(:,:,59)); % MoB-ELM (E:ELM); for 2013262 image, not in Rrs
%60 Rrs_865_MOB
%61 Rrs_1600_MOB
%62 Rrs_2100_MOB
%%
densityflag = 1;
regressiontype = 'RMA';
maxref = 0.030;
date = '2013262';

% Band 1: 443nm

denscatplot(Rrs_443A,Rrs_443E,regressiontype,densityflag,'443',maxref,date,'Acolite','MoB-ELM')
denscatplot(Rrs_443A,Rrs_443S,regressiontype,densityflag,'443',maxref,date,'Acolite','SeaDAS')
denscatplot(Rrs_443A,Rrs_443M,regressiontype,densityflag,'443',maxref,date,'Acolite','MUMM')
denscatplot(Rrs_443S,Rrs_443E,regressiontype,densityflag,'443',maxref,date,'SeaDAS','MoB-ELM')
denscatplot(Rrs_443S,Rrs_443M,regressiontype,densityflag,'443',maxref,date,'SeaDAS','MUMM')
denscatplot(Rrs_443M,Rrs_443E,regressiontype,densityflag,'443',maxref,date,'MUMM','MoB-ELM')

% Band 2: 483nm

denscatplot(Rrs_483A,Rrs_483E,regressiontype,densityflag,'483',maxref,date,'Acolite','MoB-ELM')
denscatplot(Rrs_483A,Rrs_483S,regressiontype,densityflag,'483',maxref,date,'Acolite','SeaDAS')
denscatplot(Rrs_483A,Rrs_483M,regressiontype,densityflag,'483',maxref,date,'Acolite','MUMM')
denscatplot(Rrs_483S,Rrs_483E,regressiontype,densityflag,'483',maxref,date,'SeaDAS','MoB-ELM')
denscatplot(Rrs_483S,Rrs_483M,regressiontype,densityflag,'483',maxref,date,'SeaDAS','MUMM')
denscatplot(Rrs_483M,Rrs_483E,regressiontype,densityflag,'483',maxref,date,'MUMM','MoB-ELM')

% Band 3: 561nm

denscatplot(Rrs_561A,Rrs_561E,regressiontype,densityflag,'561',maxref,date,'Acolite','MoB-ELM')
denscatplot(Rrs_561A,Rrs_561S,regressiontype,densityflag,'561',maxref,date,'Acolite','SeaDAS')
denscatplot(Rrs_561A,Rrs_561M,regressiontype,densityflag,'561',maxref,date,'Acolite','MUMM')
denscatplot(Rrs_561S,Rrs_561E,regressiontype,densityflag,'561',maxref,date,'SeaDAS','MoB-ELM')
denscatplot(Rrs_561S,Rrs_561M,regressiontype,densityflag,'561',maxref,date,'SeaDAS','MUMM')
denscatplot(Rrs_561M,Rrs_561E,regressiontype,densityflag,'561',maxref,date,'MUMM','MoB-ELM')

% Band 4: 655nm

denscatplot(Rrs_655A,Rrs_655E,regressiontype,densityflag,'655',maxref,date,'Acolite','MoB-ELM')
denscatplot(Rrs_655A,Rrs_655S,regressiontype,densityflag,'655',maxref,date,'Acolite','SeaDAS')
denscatplot(Rrs_655A,Rrs_655M,regressiontype,densityflag,'655',maxref,date,'Acolite','MUMM')
denscatplot(Rrs_655S,Rrs_655E,regressiontype,densityflag,'655',maxref,date,'SeaDAS','MoB-ELM')
denscatplot(Rrs_655S,Rrs_655M,regressiontype,densityflag,'655',maxref,date,'SeaDAS','MUMM')
denscatplot(Rrs_655M,Rrs_655E,regressiontype,densityflag,'655',maxref,date,'MUMM','MoB-ELM')

%% Chl
densityflag = 0;
regressiontype = 'RMA';
maxref = 100;
date = '2013262';

% denscatplot(CHL_ACO,CHL_MOB,regressiontype,densityflag,'',maxref,date,'Acolite','MoB-ELM')
% denscatplot(CHL_ACO,CHL_SEA,regressiontype,densityflag,'',maxref,date,'Acolite','SeaDAS')
% denscatplot(CHL_ACO,CHL_MUM,regressiontype,densityflag,'',maxref,date,'Acolite','MUMM')
denscatplot(CHL_SEA,CHL_MOB,regressiontype,densityflag,'',maxref,date,'SeaDAS','MoB-ELM')
denscatplot(CHL_SEA,CHL_MUM,regressiontype,densityflag,'',maxref,date,'SeaDAS','MUMM')
denscatplot(CHL_MUM,CHL_MOB,regressiontype,densityflag,'',maxref,date,'MUMM','MoB-ELM')
%if (rhos_469 != NaN and rhos_555 != NaN and rhos_645 != NaN)  then (       if (LAND)        then            (.091935692 + .61788 * atan(10*(rhos_645-.015)) )       else           (0.29319407 + 0.45585 * atan(50*(rhos_645-.015)) ) ) else NaN