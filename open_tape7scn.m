date = '130919';
% date = '140602';
pathname = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/';

filename = 'tape7.scn';
filepath = [pathname,date,'/',filename];
fid = fopen(filepath);
s = textscan(fid, '%s', 'delimiter', '\n');
% search for "WAVLEN":
S = strfind(s{1},'WAVLEN');
idx1 = find(~cellfun('isempty', S));
% and read from s{1}(idx1+1:idx2-1) the values using textscan again ...
data = s{1}(idx1(1)+1:size(s{1},1)-1);
fclose(fid);

for index = 1:size(data,1)
    tape7scnDATA(index,:) = str2num(cell2mat(data(index)));
end

% WAVLEN MCRN  TRANS  PTH THRML  THRML SCT  SURF EMIS   SOL SCAT  SING SCAT  GRND RFLT  DRCT RFLT  TOTAL RAD  REF SOL  SOL@OBS   DEPTH
%    0.400000 0.3731 0.0000E+00            0.0000E+00 1.0970E-02 4.9229E-03 0.0000E+00 0.0000E+00 1.0970E-02 6.22E-02 7.25E-02

WAVLEN_MCRN = tape7scnDATA(:,1);
TRANS       = tape7scnDATA(:,2);
PTH_THRML   = tape7scnDATA(:,3);
SURF_EMIS   = tape7scnDATA(:,4);
SOL_SCAT    = tape7scnDATA(:,5);
SING_SCAT   = tape7scnDATA(:,6);
GRND_RFLT   = tape7scnDATA(:,7);
DRCT_RFLT   = tape7scnDATA(:,8);
TOTAL_RAD   = tape7scnDATA(:,9);
REF_SOL     = tape7scnDATA(:,10);
SOL_OBS     = tape7scnDATA(:,11);

% SOL_SCATsmooth = sgolayfilt(SOL_SCAT,5,41);


figure
fs = 15;
set(gcf,'color','white')
plot(WAVLEN_MCRN*1E3,SOL_SCAT*1E7,'b')
hold on
% plot(WAVLEN_MCRN*1E3,TOTAL_RAD*1E7,'r')
% plot(WAVLEN_MCRN*1E3,SOL_SCATsmooth*1E7,'k')
legend('SOL\_SCAT')
% plot(wavelength,r./100)

title('tape7.scn -- SOL\_SCAT ','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('radiance [\muW/sr/cm^2/n m]','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.015])
grid on

% hold on
% plot(wavelength/1000,Lskya./max(Lskya),'b')

%% Comparison between mea. and MODTRAN
figure(20)
fs = 15;
set(gcf,'color','white')
plot(WAVLEN_MCRN*1E3,SOL_SCAT*1E7,'k')
hold on
plot(wavelength,Lskya,'r')

legend('MODTRAN 06/02/14','SVC 06/02/14','MODTRAN 09/19/13')

title('L_{sky} ','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('radiance (10^{-10}*W/(cm^2*nm*sr))','fontsize',fs)
set(gca,'fontsize',fs)
% axis([400 1000 0 0.015])
grid on
