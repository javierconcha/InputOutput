cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/
%%
pathname = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/';
pathdate = '130919/SVC_spectra_130919/SVCreflectances/';
listname = 'SVCSpectraList.txt';
listpath = [pathname,pathdate,listname];
%%
fid = fopen(listpath);
c = textscan(fid,'%s','delimiter','\n');
fclose all;

CellRef = cell([size(c{:},1),2]);

fignum = 101;
figure(fignum)
fs = 15;
set(gcf,'color','white')
% xlabel('wavelength [nm]','fontsize',fs)
% ylabel('Radiance [W/m^2/sr/nm]','fontsize',fs)
set(gca,'fontsize',fs)

for row = 1:size(c{:},1)
    filename = char(c{1}(row));
    
    filepath = [pathname,pathdate,filename];
    SVCdata = load(filepath);
    
    CellRef(row,2) = c{1}(row);
    CellRef(row,1) = {SVCdata};
    
    disp(filename)
    
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
