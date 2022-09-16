%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots the non-optimized $D^{(1,2,4)}\left(u_r,r\right)$ with respect to scale $r$ 
% and velocity increment $u_r$.
%
% Arguments IN
% evaluated = evaluated struct array from the function 'KM_Calculation'
% Fs = Acquisition/Sampling Frequency in Hz
% taylor_L = Taylor length scale in meters
% point = Multipoint condition 1=YES or 2=NO
% condition = This input is for multipoint statistics
% norm_ur = normalization of the data using $\sigma_\infty$? data 1=Yes, 0=No
% norm_r = normalization of the scale using $\lambda$? data 1=Yes, 0=No
% save_path = path for saving figures and files
% save_name = name for saving files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function KM_plot_raw(evaluated,Fs,taylor_L,point,condition,norm_ur,norm_r,save_path,save_name)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');

if point==1
    k=askInput({sprintf('Plot: Multi-point number of point condition:')}, {num2str(ceil(length(condition)/2))});
    my_field = sprintf('point_eval_%d', k);  
    evaluated=evaluated.(my_field);
end

aufl    = size(evaluated(1).D1,2);
xx      = nan(length(evaluated),aufl);
y       = nan(length(evaluated),aufl);
z1      = nan(length(evaluated),aufl);
z2      = nan(length(evaluated),aufl);
z3      = nan(length(evaluated),aufl);
z4      = nan(length(evaluated),aufl);

% for i=find((evaluated.r)/taylor_L>4):length(evaluated) 
for i=1:length(evaluated)
    xx(i,1:sum(evaluated(i).x_bin_not_nan))=evaluated(i).y_mean_bin(evaluated(i).x_bin_not_nan);
    if norm_r==1
    y(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).r./taylor_L);
    else
    y(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).r);   
    end
    z1(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).D1(evaluated(i).x_bin_not_nan));
    z2(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).D2(evaluated(i).x_bin_not_nan));
    z3(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).D3(evaluated(i).x_bin_not_nan));
    z4(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).D4(evaluated(i).x_bin_not_nan));
end


[xData, yData, zData] = prepareSurfaceData( xx, y, (z2.^2)./z4);
h(6) = figure;
scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
hold on
set(gca,'yScale','log') 
% set(gca,'zScale','log')  
set(gcf, 'Color', 'w')
set(gca, 'FontSize',18)
if norm_ur==1
    xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
else
    xlabel('$u_r\ (m/s)$', 'interpreter','latex')   
end
if norm_r==1
    ylabel('$r / \lambda$', 'interpreter','latex')
else
    ylabel('$r\ (m)$', 'interpreter','latex')
end
zlabel('$D^{(2)} (u_r,r)^2 / D^{(4)} (u_r,r)$', 'interpreter','latex')
colormap (parula(200))
% colorbar
grid off
xlim([min(xData) max(xData)])
ylim([min(y(:,1)) max(y(:,1))])
% zlim([0 2*max(zData)])
% caxis([0 max(zData)])
axis square
set(gca,'YDir','reverse');
rotate3d on



[xData, yData, zData] = prepareSurfaceData( xx, y, z4./(z2.^2));
h(5) = figure;
scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
hold on
set(gca,'yScale','log') 
% set(gca,'zScale','log')  
set(gcf, 'Color', 'w')
set(gca, 'FontSize',18)
if norm_ur==1
    xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
else
    xlabel('$u_r\ (m/s)$', 'interpreter','latex')   
end
if norm_r==1
    ylabel('$r / \lambda$', 'interpreter','latex')
else
    ylabel('$r\ (m)$', 'interpreter','latex')
end
zlabel('$D^{(4)} (u_r,r)/D^{(2)} (u_r,r)^2$', 'interpreter','latex')
colormap (parula(200))
% colorbar
grid off
xlim([min(xData) max(xData)])
ylim([min(y(:,1)) max(y(:,1))])
% zlim([0 2*max(zData)])
% caxis([0 max(zData)])
axis square
set(gca,'YDir','reverse');
rotate3d on
fig_setup
legend off

[xData, yData, zData] = prepareSurfaceData( xx, y, z4);
h(4) = figure;
scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
set(gca,'yScale','log') 
% set(gca,'zScale','log')  
set(gcf, 'Color', 'w')
set(gca, 'FontSize',18)
if norm_ur==1
    xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
else
    xlabel('$u_r\ (m/s)$', 'interpreter','latex')   
end
if norm_r==1
    ylabel('$r / \lambda$', 'interpreter','latex')
else
    ylabel('$r\ (m)$', 'interpreter','latex')
end
zlabel('$D^{(4)} (u_r,r)$', 'interpreter','latex')
colormap (parula(200))
% colorbar
grid off
xlim([min(xData) max(xData)])
ylim([min(y(:,1)) max(y(:,1))])
% zlim([0 2*max(zData)])
% caxis([0 max(zData)])
axis square
% close
set(gca,'YDir','reverse');
rotate3d on
fig_setup
legend off

[xData, yData, zData] = prepareSurfaceData( xx, y, z3);
h(3) = figure;
scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
set(gca,'yScale','log') 
% set(gca,'zScale','log')  
set(gcf, 'Color', 'w')
set(gca, 'FontSize',18)
if norm_ur==1
    xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
else
    xlabel('$u_r\ (m/s)$', 'interpreter','latex')   
end
if norm_r==1
    ylabel('$r / \lambda$', 'interpreter','latex')
else
    ylabel('$r\ (m)$', 'interpreter','latex')
end
zlabel('$D^{(3)} (u_r,r)$', 'interpreter','latex')
colormap (parula(200))
% colorbar
grid off
xlim([min(xData) max(xData)])
ylim([min(y(:,1)) max(y(:,1))])
% zlim([0 2*max(zData)])
% caxis([0 max(zData)])
axis square
% close
set(gca,'YDir','reverse');
rotate3d on
fig_setup
legend off

[xData, yData, zData] = prepareSurfaceData( xx, y, z2);
h(2) = figure;
scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
hold on
set(gca,'yScale','log') 
% set(gca,'zScale','log')  
set(gcf, 'Color', 'w')
set(gca, 'FontSize',18)
if norm_ur==1
    xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
else
    xlabel('$u_r\ (m/s)$', 'interpreter','latex')   
end
if norm_r==1
    ylabel('$r / \lambda$', 'interpreter','latex')
else
    ylabel('$r\ (m)$', 'interpreter','latex')
end
zlabel('$D^{(2)} (u_r,r)$', 'interpreter','latex')
colormap (parula(200))
% colorbar
grid off
xlim([min(xData) max(xData)])
ylim([min(y(:,1)) max(y(:,1))])
% zlim([0 2*max(zData)])
% caxis([0 max(zData)])
axis square
set(gca,'YDir','reverse');
rotate3d on
fig_setup
legend off
txt = {'(b)'};
b=text(4,0.5,txt);
b.Units='normalized';

[xData, yData, zData] = prepareSurfaceData( xx, y, z1);
h(1) = figure;
scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca, 'FontSize',18)
if norm_ur==1
    xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
else
    xlabel('$u_r\ (m/s)$', 'interpreter','latex')   
end
if norm_r==1
    ylabel('$r / \lambda$', 'interpreter','latex')
else
    ylabel('$r\ (m)$', 'interpreter','latex')
end
zlabel('$D^{(1)} (u_r,r)$', 'interpreter','latex')
colormap (parula(200))
% colorbar
grid off
xlim([min(xData) max(xData)])
ylim([min(y(:,1)) max(y(:,1))])
% zlim([-0.005 0.005])
axis square
set(gca,'YDir','reverse');
rotate3d on
fig_setup
legend off
txt = {'(a)'};
a=text(4,0.5,txt);
a.Units='normalized';

font=22;
a.FontSize = font;
b.FontSize = font;

pos_txt=[-0.23   0.9];
a.Position=pos_txt;
b.Position=pos_txt;

if ischar(save_path)
    savefig(h,fullfile(save_path,append(save_name,'_','KM_raw.fig')),'compact')
    for a = 1:length(h)
        exportgraphics(h(a),fullfile(save_path,append(save_name,'_',sprintf('KM_raw_%d.png', a))))
    end
end

for i=1:6
    tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
    uiwait(gcf);
    delete(tmp_ui);
    close
end
end