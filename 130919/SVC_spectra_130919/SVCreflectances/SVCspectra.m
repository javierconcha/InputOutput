cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/
%%
date = '130919';
pathname = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/';
pathdate = [date,'/SVC_spectra_',date,'/SVCreflectances/'];
listname = 'SVCSpectraList.txt';
listpath = [pathname,pathdate,listname];
%%
fid = fopen(listpath);
c = textscan(fid,'%s','delimiter','\n');
fclose all;

CellRef = cell([size(c{:},1),2]);

fignum = 102;
figure(fignum)
fs = 15;
set(gcf,'color','white')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Reflectance [%]','fontsize',fs)
set(gca,'fontsize',fs)

for row = 1:size(c{:},1)
    filename = char(c{1}(row));
    
    filepath = [pathname,pathdate,filename];
    SVCdata = load(filepath);
    
    CellRef(row,2) = c{1}(row);
    CellRef(row,1) = {SVCdata};
    
    disp(filename)
    disp(row)
    
    wl = CellRef{row,1}(:,1);
    refl = CellRef{row,1}(:,4);
    
    figure(fignum)
    hold on
    plot(wl,refl)
    title(filename)
%     pause(1)
    
end

% SpectraList = open(''.join((dirname,'spectralist.txt')))
% 
% for line in SpectraList:
%% Dark pixel - measurements taken in the NS Ontario - ONTNS
figure
fs = 15;
set(gcf,'color','white')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('reflectance [%]','fontsize',fs)
set(gca,'fontsize',fs)
D = zeros(size(CellRef{1,1}(:,4)));
for j = 1:6

    hold on
    plot(CellRef{1,1}(:,1),CellRef{j,1}(:,4))
    
    if j == 1
        D = CellRef{j,1}(:,4);
    else
        D = [D CellRef{j,1}(:,4)];
    end
    
    xlim([400 1000])
    disp(j)
%   pause(0.5)
    
end

%% Dark pixels stats
figure
fs = 15;
set(gcf,'color','white')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('reflectance [%]','fontsize',fs)
set(gca,'fontsize',fs)
plot(CellRef{1,1}(:,1),mean(D,2),'k')
hold on
plot(CellRef{1,1}(:,1),mean(D,2)-1*std(D,0,2),'g');
plot(CellRef{1,1}(:,1),mean(D,2)+1*std(D,0,2),'g');
plot(CellRef{1,1}(:,1),max(D,[],2),'r')
plot(CellRef{1,1}(:,1),min(D,[],2),'r')
legend('mean','mean-std','mean+std','max','min')
xlim([min(CellRef{1,1}(:,1)) 1000])


%% Save for ENVI resampling to L8 response
wavelengthSVC = CellRef{1,1}(:,1)*10^-3;
cond = wavelengthSVC>=0.4 & wavelengthSVC<=1.0;
wavelengthSVC = wavelengthSVC(cond);
meanRef = mean(D,2)./100; % decimal
ONTNSRef = meanRef(cond);

cond2 = wavelengthSVC>=0.8 & wavelengthSVC<=0.9;


ONTNSRefcorr = ONTNSRef-mean(ONTNSRef(cond2));


figure
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
% plot(wavelengthSVC.*1E3,ONTNSRefcorr,'k')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('reflectance','fontsize',fs)

ONTSRefL8 = spect_sampL8(ONTNSRefcorr,wavelengthSVC);
L8bands = [0.4430,0.4826,0.5613,0.6546,0.8646,1.6090,2.2010];

hold on
plot(L8bands*1E3,ONTSRefL8,'.-r')
plot(L8bands.*1E3,RrsONTNSL8corr*pi,'.-b')
legend('old','new')
title('ONTNS 09/19/13')
grid on


NSRef = [wavelengthSVC, ONTNSRef];
save('ONTNSRef.txt','NSRef','-ascii')
%% Bright pixel - measurements taken in the Charlotte beach
figure
fs = 15;
set(gcf,'color','white')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('reflectance [%]','fontsize',fs)
set(gca,'fontsize',fs)

j = 58-6:58+1;
D = zeros(size(CellRef{1,1}(:,4)));

for t = 1:size(j,2)

    hold on
    plot(CellRef{1,1}(:,1),CellRef{j(t),1}(:,4))
    
    if t == 1
        D = CellRef{j(t),1}(:,4);
    else
        D = [D CellRef{j(t),1}(:,4)];
    end
    
%     xlim([400 1000])
    disp(j(t))
%     pause(0.5)  
end
figure(gcf)
hold on
% 427 478 546 608 659 724 831 908x
m = get(gca,'ylim');
lw = 1.0;
line([L8bands(1)*1000 L8bands(1)*1000],m,'Color','c','LineWidth',lw)
line([L8bands(2)*1000 L8bands(2)*1000],m,'Color','b','LineWidth',lw)
line([L8bands(3)*1000 L8bands(3)*1000],m,'Color','g','LineWidth',lw)
line([L8bands(4)*1000 L8bands(4)*1000],m,'Color','r','LineWidth',lw)
line([L8bands(5)*1000 L8bands(5)*1000],m,'Color','m','LineWidth',lw)
line([L8bands(6)*1000 L8bands(6)*1000],m,'Color','k','LineWidth',lw)
line([L8bands(7)*1000 L8bands(7)*1000],m,'Color','k','LineWidth',lw)
%% Bright pixels stats
figure
fs = 15;
set(gcf,'color','white')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('reflectance [%]','fontsize',fs)
set(gca,'fontsize',fs)
plot(CellRef{1,1}(:,1),mean(D,2),'k')
hold on
plot(CellRef{1,1}(:,1),mean(D,2)-1*std(D,0,2),'g');
plot(CellRef{1,1}(:,1),mean(D,2)+1*std(D,0,2),'g');
plot(CellRef{1,1}(:,1),max(D,[],2),'r')
plot(CellRef{1,1}(:,1),min(D,[],2),'r')
legend('mean','mean-std','mean+std','max','min')
xlim([min(CellRef{1,1}(:,1)) 1000])
%% Save for ENVI resampling to L8 response
wavelength = CellRef{1,1}(:,1)*10^-3;
cond = wavelength>=0.4 & wavelength<=1.2;
wavelength = wavelength(cond);
meanRef = mean(D,2)./100; % decimal
Ref = meanRef(cond);

figure
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(wavelength,Ref,'k')
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)


CharlotteRef = [wavelength, Ref];
save('Charlotte.txt','CharlotteRef','-ascii')
%% Ibay pixel - measurements taken in Irondequoit bay IBAYN
figure
fs = 15;
set(gcf,'color','white')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('reflectance [%]','fontsize',fs)
set(gca,'fontsize',fs)

j = 42:51;
D = zeros(size(CellRef{1,1}(:,4)));

for t = 1:size(j,2)

    hold on
    plot(CellRef{1,1}(:,1),CellRef{j(t),1}(:,4))
    
    if t == 1
        D = CellRef{j(t),1}(:,4);
    else
        D = [D CellRef{j(t),1}(:,4)];
    end
    
    xlim([400 1000])
    disp(j(t))
    pause(0.5)
    
end
%% Ibay pixel stats
figure
fs = 15;
set(gcf,'color','white')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('reflectance [%]','fontsize',fs)
set(gca,'fontsize',fs)
plot(CellRef{1,1}(:,1),mean(D,2),'k')
hold on
plot(CellRef{1,1}(:,1),mean(D,2)-1*std(D,0,2),'g');
plot(CellRef{1,1}(:,1),mean(D,2)+1*std(D,0,2),'g');
plot(CellRef{1,1}(:,1),max(D,[],2),'r')
plot(CellRef{1,1}(:,1),min(D,[],2),'r')
legend('mean','mean-std','mean+std','max','min')
xlim([min(CellRef{1,1}(:,1)) 1000])
%% Save for ENVI resampling to L8 response
wavelength = CellRef{1,1}(:,1)*10^-3;
cond = wavelength>=0.4 & wavelength<=1.2;
wavelength = wavelength(cond);
meanRef = mean(D,2)./100; % decimal
Ref = meanRef(cond);

figure
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(wavelength,Ref,'k')
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)


Ibay = [wavelength, Ref];
save('IBAYNRef.txt','Ibay','-ascii')
%% ONTOS pixel - measurements taken in Lake Ontario off shore
figure
fs = 15;
set(gcf,'color','white')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('reflectance [%]','fontsize',fs)
set(gca,'fontsize',fs)

j = 16:20;
D = zeros(size(CellRef{1,1}(:,4)));

for t = 1:size(j,2)

    hold on
    plot(CellRef{1,1}(:,1),CellRef{j(t),1}(:,4))
    
    if t == 1
        D = CellRef{j(t),1}(:,4);
    else
        D = [D CellRef{j(t),1}(:,4)];
    end
    
    xlim([400 1000])
    disp(j(t))
    pause(0.5)
    
end
%% ONTOS pixel stats
figure
fs = 15;
set(gcf,'color','white')
xlabel('wavelength [nm]','fontsize',fs)
ylabel('reflectance [%]','fontsize',fs)
set(gca,'fontsize',fs)
plot(CellRef{1,1}(:,1),mean(D,2),'k')
hold on
plot(CellRef{1,1}(:,1),mean(D,2)-1*std(D,0,2),'g');
plot(CellRef{1,1}(:,1),mean(D,2)+1*std(D,0,2),'g');
plot(CellRef{1,1}(:,1),max(D,[],2),'r')
plot(CellRef{1,1}(:,1),min(D,[],2),'r')
legend('mean','mean-std','mean+std','max','min')
xlim([min(CellRef{1,1}(:,1)) 1000])
%% Save for ENVI resampling to L8 response
wavelength = CellRef{1,1}(:,1)*10^-3;
cond = wavelength>=0.4 & wavelength<=1.2;
wavelength = wavelength(cond);
meanRef = mean(D,2)./100; % decimal
Ref = meanRef(cond);

figure
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(wavelength,Ref,'k')
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)


ONTOS = [wavelength, Ref];
save('ONTOS.txt','ONTOS','-ascii')



