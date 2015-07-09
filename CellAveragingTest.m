
figure
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
% plot(x_used,y_used,'*')
hold on
plot(x1,y1,'r','LineWidth',2')
plot(ylim,ylim,'k')
axis equal
% axis([0 0.02 0 0.02])

N = 11;
x_interval = linspace(0,max(x_used),N);
y_interval = linspace(0,max(y_used),N);

x_cellmean = nan(size(x_interval,2),size(y_interval,2));
y_cellmean = nan(size(x_interval,2),size(y_interval,2));
cell_density = nan(size(x_interval,2),size(y_interval,2));

for x_index = 1:size(x_interval,2)-1
    for y_index = 1:size(y_interval,2)-1
        x_rule = x_used>x_interval(x_index)&x_used<=x_interval(x_index+1);
        x_cell = x_used(x_rule);
        y_cell = y_used(x_rule);
        y_rule = y_cell>y_interval(y_index)&y_cell<=y_interval(y_index+1);
        
        x_cellmean(x_index,y_index) = mean(x_cell(y_rule));
        y_cellmean(x_index,y_index) = mean(y_cell(y_rule));
        cell_density(x_index,y_index) = sum(y_rule);
        
    end
end

plot(x_cellmean(:),y_cellmean(:),'g.')

% OLS
[a3,~] = polyfit(x_cellmean(isfinite(x_cellmean(:))),y_cellmean(isfinite(y_cellmean(:))),1);
x3=[0 maxref];
y3=a3(1).*x3+a3(2);

plot(x3,y3,'b-','LineWidth',2)

% RMA
a(1) = nanstd(y_cellmean(:))/nanstd(x_cellmean(:)); % slope

if corr(x_cellmean(isfinite(x_cellmean(:))),y_cellmean(isfinite(y_cellmean(:))))<0
    a(1) = -abs(a(1));
elseif corr(x_cellmean(isfinite(x_cellmean(:))),y_cellmean(isfinite(y_cellmean(:))))>=0
    a(1) = abs(a(1));
end

a(2) = nanmean(y_cellmean(:))-nanmean(x_cellmean(:))*a(1); % y intercept



x2=[0 maxref];
y2=a(1).*x2+a(2);

plot(x2,y2,'c-','LineWidth',2)

[xx,yy] = meshgrid(x_interval,y_interval);

% plot(xx(:),yy(:),'.k')
% legend('data','Linear Regression','1:1','cell mean')
%%
figure
surf(xx,yy,cell_density)
view(2)
colorbar