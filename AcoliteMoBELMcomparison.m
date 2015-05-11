tif_file = 'LC80160302013262LGN00/Collocated2013262.tif';
filename = [dir tif_file];

im2013262 = imread(filename);
proj = geotiffinfo(filename);
info = imfinfo(filename);

%% Band 1: 443nm
rho_443 = im2013262(:,:,4);
Rrs_443A = double(rho_443./pi);


Rrs_443E = double(im2013262(:,:,9));

cond1 = Rrs_443A(:)>0;
cond2 = Rrs_443E(:)>0;
cond3 = cond1&cond2;
Rrs_443A_used = Rrs_443A(cond3);
Rrs_443E_used = Rrs_443E(cond3);

disp('Band 443nm: Acolite Pos. [%]:')
disp(100*sum(cond1)/size(Rrs_443A(:),1))
disp('Band 443nm: MoB-ELM Pos. [%]:')
disp(100*sum(cond2)/size(Rrs_443E(:),1))

figure
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(Rrs_443A_used,Rrs_443E_used,'.')
xlabel('Rrs 443 nm (Acolite) [1/sr]')
ylabel('Rrs 443 nm (MoB-ELM) [1/sr]')
grid on
axis equal
maxref = 0.015;
axis([0 maxref 0 maxref])
hold on
plot([0 maxref],[0 maxref],'--k')

Rrs_443_RMSE = sqrt(mean((Rrs_443A_used-Rrs_443E_used).^2));
hold on % regression
[a,b] = polyfit(Rrs_443A_used,Rrs_443E_used,1);
x1=[0 maxref];
y1=a(1).*x1+a(2);
plot(x1,y1,'r-','LineWidth',1)
C = corrcoef([Rrs_443A_used,Rrs_443E_used]);
r2 = C(1,2)^2;
str1 = sprintf('y: %2.4f x + %2.4f \n R^2: %2.4f; N: %i \n RMSE: %2.4f',a(1),a(2),r2,size(Rrs_443A_used,1),Rrs_443_RMSE);
xLimits = get(gca,'XLim');
yLimits = get(gca,'YLim');
xLoc = xLimits(1)+0.1*(xLimits(2)-xLimits(1));
yLoc = yLimits(1)+0.85*(yLimits(2)-yLimits(1));
h = text(xLoc,yLoc,str1,'FontSize',fs,'FontWeight','normal');


SStot = sum((Rrs_443E_used-mean(Rrs_443E_used)).^2);
SSres = sum((Rrs_443E_used-polyval(a,Rrs_443A_used)).^2);
rsq = 1-(SSres/SStot)




% %%%%%%%% RMA Regression %%%%%%%%%%%%%
% [[b1 b0],bintr,bintjm] = gmregress(Rrs_443A_used,Rrs_443E_used);
b1 = std(Rrs_443E_used)/std(Rrs_443A_used); % slope
b0 = mean(Rrs_443E_used)-mean(Rrs_443A_used)*b1; % y intercept
y3=b1.*x1+b0;
plot(x1,y3,'g-','LineWidth',1)

SStot = sum((Rrs_443E_used-mean(Rrs_443E_used)).^2);
SSres = sum((Rrs_443E_used-polyval([b1 b0],Rrs_443A_used)).^2);
rsq = 1-(SSres/SStot)
%% Density
po = 3;
method = 'squares';
radius = max([Rrs_443A_used;Rrs_443E_used])/512/2;
% radius = [];
N = 512;
n = [];
ms = 5;
figure
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
out = scatplot(Rrs_443A_used,Rrs_443E_used,method,radius,N,n,po,ms);
xlabel('Rrs 443 nm (Acolite) [1/sr]')
ylabel('Rrs 443 nm (MoB-ELM) [1/sr]')
colorbar('southoutside')
grid on
axis equal
maxref = 0.015;
axis([0 maxref 0 maxref])
hold on
plot([0 maxref],[0 maxref],'--k')

Rrs_443_RMSE = sqrt(mean((Rrs_443A_used-Rrs_443E_used).^2)); 

% regression
[a,b] = polyfit(Rrs_443A_used,Rrs_443E_used,1);
x1=[0 maxref];
y1=a(1).*x1+a(2);
plot(x1,y1,'r-','LineWidth',1)
C = corrcoef([Rrs_443A_used,Rrs_443E_used]);
r2 = C(1,2)^2;
str1 = sprintf('y: %2.4f x + %2.4f \n R^2: %2.4f; N: %i \n RMSE: %2.4f',a(1),a(2),r2,size(Rrs_443A_used,1),Rrs_443_RMSE);
xLimits = get(gca,'XLim');
yLimits = get(gca,'YLim');
xLoc = xLimits(1)+0.1*(xLimits(2)-xLimits(1));
yLoc = yLimits(1)+0.85*(yLimits(2)-yLimits(1));
h = text(xLoc,yLoc,str1,'FontSize',16,'FontWeight','bold');
%%
figure
hAxes = dscatter(Rrs_443A_used,Rrs_443E_used,'MARKER','.','BINS',[512,512] )
colorbar
%%
figure
values = hist3([Rrs_443E_used Rrs_443A_used],[512 512]);
imagesc(values)
colorbar
axis equal
axis xy
