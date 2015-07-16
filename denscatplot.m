function denscatplot(x,y,regressiontype,densityflag,bandname,maxref,date,labelx,labely)% density scatterplot
%% Preparing data

cond1 = x(:)>0;
cond2 = y(:)>0;
cond3 = cond1&cond2;
cond4 = x(:)<0;
cond5 = y(:)<0;
x_used = x(cond3);
y_used = y(cond3);

str0 = sprintf('Band %s: x Neg: %2.2f%%; y Neg: %2.2f%%;usable: %2.0f%%',...
    bandname,100*sum(cond4)/sum(isfinite(x(:))),100*sum(cond5)/sum(isfinite(y(:))),...
    100*sum(cond3)/sum(isfinite(y(:))));
disp(str0)

h = figure;
fs = 15;
set(gcf,'color','white')
set(gcf,'name',[date ' ' regressiontype ' ' labelx ' ' labely])
set(gca,...
    'fontsize'  , fs    , ...
    'XMinorTick', 'on'  , ...
    'YMinorTick', 'on'  );

%% Density 
if  densityflag == 0
    figure(h)
    hold on
    plot(x_used,y_used,'.')
elseif densityflag == 1
    po = 3;
    method = 'squares';
    radius = max([x_used;y_used])/512/2;
    % radius = [];
    N = 512;
    n = [];
    ms = 5;
    figure(h)
    hold on
    [~] = scatplot(x_used,y_used,method,radius,N,n,po,ms);
    colorbar off
    hc = colorbar('southoutside','FontSize',14);
    title(hc, 'Pixel Density','FontSize',14)
    set(gca,'position',[0.15 0.28 0.75 0.7]);
    set(hc,'position',[0.2 0.07 0.6 0.05]);
end

figure(h)
hold on
str1 = sprintf('R_{rs} %s nm (%s) [1/sr]',bandname,labelx);
str2 = sprintf('R_{rs} %s nm (%s) [1/sr]',bandname,labely);
xlabel(str1)
ylabel(str2)
% grid on
axis equal
% axis([0 maxref 0 maxref])
box on

figure(h)
hold on
plot([0 maxref],[0 maxref],'--k','LineWidth',1.2)
%% regression
if strcmp(regressiontype,'OLS') 
    [a,~] = polyfit(x_used,y_used,1);
elseif strcmp(regressiontype,'RMA')
    % %%%%%%%% RMA Regression %%%%%%%%%%%%%
    % [[b1 b0],bintr,bintjm] = gmregress(x_used,y_used);
    a(1) = std(y_used)/std(x_used); % slope
    
    if corr(x_used,y_used)<0
        a(1) = -abs(a(1));
    elseif corr(x_used,y_used)>=0
        a(1) = abs(a(1));
    end
    
    a(2) = mean(y_used)-mean(x_used)*a(1); % y intercept
    
end

x1=[0 maxref];
y1=a(1).*x1+a(2);

figure(h)
plot(x1,y1,'r-','LineWidth',1.2)

%% Statistics: r-squared (or coefficient of determination) 
SStot = sum((y_used-mean(y_used)).^2);
SSres = sum((y_used-polyval(a,x_used)).^2);
rsq = 1-(SSres/SStot)
corr(x_used,y_used)^2

RMSE = sqrt(mean((x_used-y_used).^2));

if a(2)>=0
    str1 = sprintf('y: %2.4f x + %2.4f \n R^2: %2.4f; N: %i \n RMSE: %2.4f',...
        a(1),abs(a(2)),rsq,size(x_used,1),RMSE);
else
    str1 = sprintf('y: %2.4f x - %2.4f \n R^2: %2.4f; N: %i \n RMSE: %2.4f',...
        a(1),abs(a(2)),rsq,size(x_used,1),RMSE);
end

axis([0 maxref 0 maxref])

xLimits = get(gca,'XLim');
yLimits = get(gca,'YLim');
xLoc = xLimits(1)+0.1*(xLimits(2)-xLimits(1));
yLoc = yLimits(1)+0.85*(yLimits(2)-yLimits(1));
figure(h)
hold on
text(xLoc,yLoc,str1,'FontSize',fs,'FontWeight','normal');
disp(str1)
% %% cell mean
% N = 101;
% x_interval = linspace(0,max(x_used),N);
% y_interval = linspace(0,max(y_used),2);
% 
% x_cellmean = nan(size(x_interval,2),size(y_interval,2));
% y_cellmean = nan(size(x_interval,2),size(y_interval,2));
% cell_density = nan(size(x_interval,2),size(y_interval,2));
% 
% for x_index = 1:size(x_interval,2)-1
% %     for y_index = 1:size(y_interval,2)-1
%         y_index = 1;
%         x_rule = x_used>x_interval(x_index)&x_used<=x_interval(x_index+1);
%         x_cell = x_used(x_rule);
%         y_cell = y_used(x_rule);
%         y_rule = y_cell>y_interval(y_index)&y_cell<=y_interval(y_index+1);
%         
%         x_cellmean(x_index,y_index) = mean(x_cell(y_rule));
%         y_cellmean(x_index,y_index) = mean(y_cell(y_rule));
%         cell_density(x_index,y_index) = sum(y_rule);
%         
% %     end
% end
% 
% % plot(x_cellmean(:),y_cellmean(:),'g.')
% % regression
% if strcmp(regressiontype,'OLS')
%     [a,~] = polyfit(isfinite(x_cellmean(:)),isfinite(y_cellmean(:)),1);
% elseif strcmp(regressiontype,'RMA')
%     a(1) = nanstd(y_cellmean(:))/nanstd(x_cellmean(:)); % slope
%     
%     if corr(x_cellmean(isfinite(x_cellmean(:))),y_cellmean(isfinite(y_cellmean(:))))<0
%         a(1) = -abs(a(1));
%     elseif corr(x_cellmean(isfinite(x_cellmean(:))),y_cellmean(isfinite(y_cellmean(:))))>=0
%         a(1) = abs(a(1));
%     end
% end
% 
% a(2) = nanmean(y_cellmean(:))-nanmean(x_cellmean(:))*a(1); % y intercept
% 
% 
% 
% x2=[0 maxref];
% y2=a(1).*x2+a(2);
% 
% plot(x2,y2,'c-','LineWidth',2)


%% Save figure
% str3 = sprintf('%s_AcoliteMoBELMcomp_%s_%s',date,bandname,comptype);
% dirname = '/Users/javier/Desktop/Javier/PHD_RIT/ConferencesAndApplications/2015_SPIE_SanDiego/Images/';
% print(h,[dirname str3],'-depsc')

