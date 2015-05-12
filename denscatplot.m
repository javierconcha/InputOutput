function denscatplot(x,y,regressiontype,densityflag,bandname,maxref)% density scatterplot
% Preparing data

disp('Band:')
disp(bandname)
cond1 = x(:)>0;
cond2 = y(:)>0;
cond3 = cond1&cond2;
x_used = x(cond3);
y_used = y(cond3);

disp('x Pos. [%]:')
disp(100*sum(cond1)/sum(isfinite(x(:))))
disp('y Pos. [%]:')
disp(100*sum(cond2)/sum(isfinite(y(:))))

%% Density 
if  densityflag == 0
    h = figure;
    plot(x_used,y_used,'.')
elseif densityflag == 1
    po = 3;
    method = 'squares';
    radius = max([x_used;y_used])/512/2;
    % radius = [];
    N = 512;
    n = [];
    ms = 5;
    h = figure;
    [~] = scatplot(x_used,y_used,method,radius,N,n,po,ms);
    colorbar off
    colorbar('southoutside')
end

figure(h)
fs = 15;
set(gcf,'color','white')
set(gcf,'name',regressiontype)
set(gca,'fontsize',fs)
str1 = sprintf('Rrs %s nm (Acolite) [1/sr]',bandname);
str2 = sprintf('Rrs %s nm (MoB-ELM) [1/sr]',bandname);
xlabel(str1)
ylabel(str2)
grid on
axis equal
axis([0 maxref 0 maxref])

figure(h)
hold on
plot([0 maxref],[0 maxref],'--k')

RMSE = sqrt(mean((x_used-y_used).^2));

figure(h)
hold on 

% regression
if strcmp(regressiontype,'OLS') 
    [a,~] = polyfit(x_used,y_used,1);
    x1=[0 maxref];
    y1=a(1).*x1+a(2);
    plot(x1,y1,'r-','LineWidth',1)
elseif strcmp(regressiontype,'RMA')
    % %%%%%%%% RMA Regression %%%%%%%%%%%%%
    % [[b1 b0],bintr,bintjm] = gmregress(x_used,y_used);
    a(1) = std(y_used)/std(x_used); % slope
    a(2) = mean(y_used)-mean(x_used)*a(1); % y intercept
    x1=[0 maxref];
    y3=a(1).*x1+a(2);
    plot(x1,y3,'r-','LineWidth',1)
    
end

% r-squared or coefficient of determination
SStot = sum((y_used-mean(y_used)).^2);
SSres = sum((y_used-polyval(a,x_used)).^2);
rsq = 1-(SSres/SStot);

str1 = sprintf('y: %2.4f x + %2.4f \n R^2: %2.4f; N: %i \n RMSE: %2.4f',a(1),a(2),rsq,size(x_used,1),RMSE);
xLimits = get(gca,'XLim');
yLimits = get(gca,'YLim');
xLoc = xLimits(1)+0.1*(xLimits(2)-xLimits(1));
yLoc = yLimits(1)+0.85*(yLimits(2)-yLimits(1));
h = text(xLoc,yLoc,str1,'FontSize',fs,'FontWeight','normal');

