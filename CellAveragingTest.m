
figure
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(x_used,y_used,'*')
hold on
plot(x1,y1,'r')
plot(ylim,ylim,'k')
axis equal
axis([0 0.02 0 0.02])

N = 101;
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

[xx,yy] = meshgrid(x_interval,y_interval);

plot(xx(:),yy(:),'.k')
legend('data','Linear Regression','1:1','cell mean','grid')
%%
figure
surf(xx,yy,cell_density)
view(2)
colorbar