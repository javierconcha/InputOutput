function denscatplot(x,y)% density scatterplot

cond1 = x(:)>0;
cond2 = y(:)>0;
cond3 = cond1&cond2;
x_used = x(cond3);
y_used = y(cond3);

disp('Band 443nm: Acolite Pos. [%]:')
disp(100*sum(cond1)/size(x(:),1))
disp('Band 443nm: MoB-ELM Pos. [%]:')
disp(100*sum(cond2)/size(y(:),1))

figure
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(x_used,y_used,'.')
xlabel('Rrs 443 nm (Acolite) [1/sr]')
ylabel('Rrs 443 nm (MoB-ELM) [1/sr]')
grid on
axis equal
maxref = 0.015;
axis([0 maxref 0 maxref])
hold on
plot([0 maxref],[0 maxref],'--k')

RMSE = sqrt(mean((x_used-y_used).^2));
hold on % regression
[a,b] = polyfit(x_used,y_used,1);
x1=[0 maxref];
y1=a(1).*x1+a(2);
plot(x1,y1,'r-','LineWidth',1)
C = corrcoef([x_used,y_used]);
r2 = C(1,2)^2;
str1 = sprintf('y: %2.4f x + %2.4f \n R^2: %2.4f; N: %i \n RMSE: %2.4f',a(1),a(2),r2,size(x_used,1),RMSE);
xLimits = get(gca,'XLim');
yLimits = get(gca,'YLim');
xLoc = xLimits(1)+0.1*(xLimits(2)-xLimits(1));
yLoc = yLimits(1)+0.85*(yLimits(2)-yLimits(1));
h = text(xLoc,yLoc,str1,'FontSize',fs,'FontWeight','normal');


SStot = sum((y_used-mean(y_used)).^2);
SSres = sum((y_used-polyval(a,x_used)).^2);
rsq = 1-(SSres/SStot)




% %%%%%%%% RMA Regression %%%%%%%%%%%%%
% [[b1 b0],bintr,bintjm] = gmregress(x_used,y_used);
b1 = std(y_used)/std(x_used); % slope
b0 = mean(y_used)-mean(x_used)*b1; % y intercept
y3=b1.*x1+b0;
plot(x1,y3,'g-','LineWidth',1)

SStot = sum((y_used-mean(y_used)).^2);
SSres = sum((y_used-polyval([b1 b0],x_used)).^2);
rsq = 1-(SSres/SStot)
%% Density
po = 3;
method = 'squares';
radius = max([x_used;y_used])/512/2;
% radius = [];
N = 512;
n = [];
ms = 5;
figure
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
out = scatplot(x_used,y_used,method,radius,N,n,po,ms);
xlabel('Rrs 443 nm (Acolite) [1/sr]')
ylabel('Rrs 443 nm (MoB-ELM) [1/sr]')
colorbar('southoutside')
grid on
axis equal
maxref = 0.015;
axis([0 maxref 0 maxref])
hold on
plot([0 maxref],[0 maxref],'--k')

RMSE = sqrt(mean((x_used-y_used).^2)); 

% regression
[a,b] = polyfit(x_used,y_used,1);
x1=[0 maxref];
y1=a(1).*x1+a(2);
plot(x1,y1,'r-','LineWidth',1)
C = corrcoef([x_used,y_used]);
r2 = C(1,2)^2;
str1 = sprintf('y: %2.4f x + %2.4f \n R^2: %2.4f; N: %i \n RMSE: %2.4f',a(1),a(2),r2,size(x_used,1),RMSE);
xLimits = get(gca,'XLim');
yLimits = get(gca,'YLim');
xLoc = xLimits(1)+0.1*(xLimits(2)-xLimits(1));
yLoc = yLimits(1)+0.85*(yLimits(2)-yLimits(1));
h = text(xLoc,yLoc,str1,'FontSize',16,'FontWeight','bold');
% %%
% figure
% hAxes = dscatter(x_used,y_used,'MARKER','.','BINS',[512,512] )
% colorbar
% %%
% figure
% values = hist3([y_used x_used],[512 512]);
% imagesc(values)
% colorbar
% axis equal
% axis xy
