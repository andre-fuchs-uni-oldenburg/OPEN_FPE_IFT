%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots the optimized $D^{(1)}\left(u_r,r\right)$ and $D^{(2)}\left(u_r,r\right)$
% and the surface fit.
%
% Arguments IN
% co = Coefficients $d_{ij}(r)$ of surface fits with a linear function for 
%      $D^{(1)}\left(u_r,r\right)$ and a parabolic function for $D^{(2)}\left(u_r,r\right)$
% evaluated = evaluated struct array from the function 'KM_Calculation'
% Fs = Acquisition/Sampling Frequency in Hz
% taylor_L = Taylor length scale in meters
% multi_point = Multipoint condition 1=YES or 2=NO
% condition = This input is for multipoint statistics
% norm_ur = normalization of the data using $\sigma_\infty$? data 1=Yes, 0=No
% norm_r = normalization of the scale using $\lambda$? data 1=Yes, 0=No
% save_path = path for saving figures and files
% save_name = name for saving files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function KM_plot(co,evaluated,Fs,taylor_L,int_L,multi_point,condition,norm_ur,norm_r,save_path,save_name,varargin)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');

if multi_point==1
    k           = askInput({sprintf('Plot: Multi-point number of point condition:')}, {num2str(ceil(length(condition)/2))});
    my_field    = sprintf('point_eval_%d', k);  
    evaluated   = evaluated.(my_field);
    co          = co.(my_field);
end

aufl        = size(evaluated(1).D1,2);
xx          = nan(length(evaluated),aufl);
y           = nan(length(evaluated),aufl);
z1          = nan(length(evaluated),aufl);
% z1e         = nan(length(evaluated),aufl);
z2          = nan(length(evaluated),aufl);
% z2e         = nan(length(evaluated),aufl);
z1_opti     = nan(length(evaluated),aufl);
z2_opti     = nan(length(evaluated),aufl);


% for i=find((evaluated.r)/taylor_L>4):length(evaluated) 
for i=1:length(evaluated)
    xx(i,1:sum(evaluated(i).x_bin_not_nan))=evaluated(i).y_mean_bin(evaluated(i).x_bin_not_nan);
    if norm_r==1
    y(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).r./taylor_L);
    else
    y(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).r);   
    end
    % y(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).r./(m_data*markov/Fs));

    z1(i,:)=(evaluated(i).D1);
    % z1e(i,1:sum(evaluated(i).x_bin_not_nan))=abs(1./(evaluated(i).eD1(evaluated(i).x_bin_not_nan)));
    z2(i,:)=(evaluated(i).D2);
    % z2e(i,1:sum(evaluated(i).x_bin_not_nan))=abs(1./(evaluated(i).eD2(evaluated(i).x_bin_not_nan)));

    z1_opti(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).D1_opti);
    z2_opti(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).D2_opti);
    
    if nargin >= 12  
        z1(i,:)=z1(i,:)./(evaluated(i).r);
        z2(i,:)=z2(i,:)./(evaluated(i).r);
        z1_opti(i,:)=z1_opti(i,:)./(evaluated(i).r);
        z2_opti(i,:)=z2_opti(i,:)./(evaluated(i).r);
    end
end

% z1_opti=z1;
% z2_opti=z2;

% z1_opti=z1_opti.*y;
% z2_opti=z2_opti.*y;


% %% d_ij(r) plots
% h(1) = figure;
% subplot(2,2,1);
% plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%     co.a(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co.ea(1)+co.a(2),'k','LineWidth',2)
% hold on
% % plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
% %     -0.5.*(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100)).^(-0.5))
% x_fit=logspace(log10(min(y(:,1))),log10(max(y(:,1))),100);
% y_fit=co.a(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co.ea(1)+co.a(2);
% [xData, yData] = prepareCurveData( x_fit, y_fit );
% ft = fittype( 'power1' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% [fitresult_d11, gof] = fit( xData, yData, ft, opts );
% plot(linspace(0,max(y(:,1)),100),feval(fitresult_d11,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);
% % plot(linspace(0,15,100),feval(fitresult,linspace(0,15,100)),'LineWidth',2 );
% axis square
% ylabel('$d_{11}(r)$', 'interpreter','latex')
% xlabel('$r / \lambda$', 'interpreter','latex')
% xlim([min(y(:,1)) max(y(:,1))])
% set(gca, 'FontSize',24)
% % set(gca,'xScale','log') 
% % set(gca,'yScale','log') 
% %
% subplot(2,2,2);
% plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%     co.b(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co.eb(1)+co.b(2),'k','LineWidth',2)
% hold on
% x_fit=logspace(log10(min(y(:,1))),log10(max(y(:,1))),100);
% y_fit=co.b(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co.eb(1)+co.b(2);
% [xData, yData] = prepareCurveData( x_fit, y_fit );
% ft = fittype( 'power1' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% [fitresult_d20, gof] = fit( xData, yData, ft, opts );
% plot(linspace(0,max(y(:,1)),100),feval(fitresult_d20,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);
% axis square
% ylabel('$d_{20}(r)$', 'interpreter','latex')
% xlabel('$r / \lambda$', 'interpreter','latex')
% xlim([min(y(:,1)) max(y(:,1))])
% set(gca, 'FontSize',24)
% % set(gca,'xScale','log') 
% %
% subplot(2,2,3);
% plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%     co.b(3).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co.eb(2)+co.b(4),'k','LineWidth',2)
% hold on
% x_fit=logspace(log10(min(y(:,1))),log10(max(y(:,1))),100);
% y_fit=co.b(3).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co.eb(2)+co.b(4);
% [xData, yData] = prepareCurveData( x_fit, y_fit );
% ft = fittype( 'power1' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% [fitresult_d21, gof] = fit( xData, yData, ft, opts );
% plot(linspace(0,max(y(:,1)),100),feval(fitresult_d21,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);
% axis square
% ylabel('$d_{21}(r)$', 'interpreter','latex')
% xlabel('$r / \lambda$', 'interpreter','latex')
% xlim([min(y(:,1)) max(y(:,1))])
% set(gca, 'FontSize',24)
% % set(gca,'xScale','log') 
% %
% subplot(2,2,4);
% plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%     co.b(5).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co.eb(3)+co.b(6),'k','LineWidth',2)
% hold on
% x_fit=logspace(log10(min(y(:,1))),log10(max(y(:,1))),100);
% y_fit=co.b(5).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co.eb(3)+co.b(6);
% [xData, yData] = prepareCurveData( x_fit, y_fit );
% ft = fittype( 'power1' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% [fitresult_d22, gof] = fit( xData, yData, ft, opts );
% plot(linspace(0,max(y(:,1)),100),feval(fitresult_d22,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);
% axis square
% ylabel('$d_{22}(r)$', 'interpreter','latex')
% xlabel('$r / \lambda$', 'interpreter','latex')
% xlim([min(y(:,1)) max(y(:,1))])
% set(gca, 'FontSize',24)
% set(gcf, 'Color', 'w')
% % set(gca,'yScale','log') 
% % set(gca,'xScale','log') 
% fitresult_d11
% fitresult_d20
% fitresult_d21
% fitresult_d22

%% 3D plots
[xData, yData, zData] = prepareSurfaceData( xx, y, z1_opti);
tmp_u=repmat(linspace(min(xData),max(xData),100),100,1);
tmp_r=repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100);


% %%
% h(1) = figure;
% scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
% hold on
% surf(tmp_u,tmp_r,(co.a(1).*tmp_r.^co.ea(1)+co.a(2)).*tmp_u,'EdgeColor','none')
% set(gca,'yScale','log') 
% set(gcf, 'Color', 'w')
% set(gca, 'FontSize',18)
% xlabel('$u_r / \sigma_\infty$', 'interpreter','latex',...
%     'Rotation',20,...
%     'VerticalAlignment','bottom',...
%     'HorizontalAlignment','left')
% % ylabel('$r / [m]$', 'interpreter','latex')
% ylabel('$r / \lambda$', 'interpreter','latex',...
%     'Rotation',325,...
%     'VerticalAlignment','middle',...
%     'HorizontalAlignment','right')
% zlabel('$D^{(1)} (u_r,r)$', 'interpreter','latex')
% % zlabel('$D1 (x,r) / [\frac{m}{s^2}]$', 'interpreter','latex')
% colormap (parula(200))
% % colorbar
% grid off
% xlim([min(xData) max(xData)])
% ylim([min(y(:,1)) max(y(:,1))])
% % zlim([-0.005 0.005])
% axis square
% set(gca,'YDir','reverse');
% rotate3d on
% 
% view(axes1,[-37.5 30]);
% set(get(gca,'xlabel'),'rotation',20);

% axislabel_rotation(h)
% 
% tmp_fig = rotate3d;
% % set(tmp_fig, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
% set(tmp_fig, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel_log)')
% set(tmp_fig, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
% align_axislabel([], gca)




%%
h(1) = figure;
scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
hold on
if nargin >= 12  
   surf(tmp_u,tmp_r,((co.a(1).*tmp_r.^co.ea(1)+co.a(2)).*tmp_u)./tmp_r,'EdgeColor','none')
else
   surf(tmp_u,tmp_r,(co.a(1).*tmp_r.^co.ea(1)+co.a(2)).*tmp_u,'EdgeColor','none')
end
surf(tmp_u,tmp_r,((-0.362.*tmp_r.^(-1)).*tmp_u./0.5399),'FaceColor',[0.501960813999176 0.501960813999176 0.501960813999176],'EdgeColor','none')
% surf(tmp_u,tmp_r,((-0.362.*(tmp_r.*taylor_L).^(-1)).*tmp_u./0.5399).*taylor_L,'FaceColor','r','EdgeColor','none')

%Renner
%         d10_R = -0.0015.*(tmp_r.*taylor_L).^0.6;
%         d11_R =   -0.61.*(tmp_r.*taylor_L).^(-0.67);
%         d12_R =  0.0096.*(tmp_r.*taylor_L).^(0);
%         d13_R = -0.0023.*(tmp_r.*taylor_L).^(0.3);      
%         d20_R =   0.033.*(tmp_r.*taylor_L).^(0.25);
%         d21_R =  -0.009.*(tmp_r.*taylor_L).^(0.2);
%         d22_R =   0.043.*(tmp_r.*taylor_L).^(-0.73);
%         d10_R = -0.0015.*(tmp_r).^0.6;
%         d11_R =   -0.61.*(tmp_r).^(-0.67);
%         d12_R =  0.0096.*(tmp_r).^(0);
%         d13_R = -0.0023.*(tmp_r).^(0.3);      
%         d20_R =   0.033.*(tmp_r).^(0.25);
%         d21_R =  -0.009.*(tmp_r).^(0.2);
%         d22_R =   0.043.*(tmp_r).^(-0.73);     
% surf(tmp_u,tmp_r,...
%     (d10_R+d11_R.*tmp_u+d12_R.*tmp_u.^2+d13_R.*tmp_u.^3)...
%     ,'FaceColor',[0.831372559070587 0.815686285495758 0.7843137383461])

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
% zlabel('$D1 (x,r) / [\frac{m}{s^2}]$', 'interpreter','latex')
colormap (parula(200))
% colorbar
grid off
xlim([min(xData) max(xData)])
ylim([min(y(:,1)) max(y(:,1))])
% zlim([-0.005 0.005])
axis square
set(gca,'YDir','reverse');
rotate3d on
colorbar('north','TickLabelInterpreter','latex','AxisLocation','out');
    fig_setup
    legend off
    txt = {'(a)'};
    a=text(4,0.5,txt);
    a.Units='normalized';

% level=round(logspace(log10(min(min((zData)))*10^6),log10(0.9*max(max((zData)))*10^6),10))./10^6;
level=round(linspace(min(min((zData))),0.9*max(max((zData))),10),2);

h(2) = figure;
% scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
% hold on
contour3(tmp_u,tmp_r,(co.a(1).*tmp_r.^co.ea(1)+co.a(2)).*tmp_u,level,'-','color','k','LineWidth',2,'ShowText','on')
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
% zlabel('$D1 (x,r) / [\frac{m}{s^2}]$', 'interpreter','latex')
grid off
box off
xlim([min(xData) max(xData)])
ylim([min(y(:,1)) max(y(:,1))])
% zlim([-0.005 0.005])
axis square
set(gca,'YDir','reverse');
rotate3d on

[xData, yData, zData] = prepareSurfaceData( xx, y, z2_opti);
tmp_u=repmat(linspace(min(xData),max(xData),100),100,1);
tmp_r=repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100);



h(3) = figure;
scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
hold on
surf(tmp_u,tmp_r,...
    ((co.b(1).*tmp_r.^co.eb(1)+co.b(2)).*tmp_u.^0)+...
    ((co.b(3).*tmp_r.^co.eb(2)+co.b(4)).*tmp_u.^1)+...
    ((co.b(5).*tmp_r.^co.eb(3)+co.b(6)).*tmp_u.^2),'EdgeColor','none')

surf(tmp_u,tmp_r,((0.0144.*tmp_r.^(-1)).*tmp_u.^2./(0.5399.^2)),'FaceColor',[0.501960813999176 0.501960813999176 0.501960813999176],'EdgeColor','none')


% surf(tmp_u,tmp_r,...
%     (d20_R+d21_R.*tmp_u+d22_R.*tmp_u.^2)...
%     ,'FaceColor',[0.831372559070587 0.815686285495758 0.7843137383461])

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
% zlabel('$D2 (x,r) / [\frac{m^2}{s^3}]$', 'interpreter','latex')
zlabel('$D^{(2)} (u_r,r)$', 'interpreter','latex')
colormap (parula(200))
% colorbar
grid off
xlim([min(xData) max(xData)])
ylim([min(y(:,1)) max(y(:,1))])
% zlim([0 2*max(zData)])
caxis([0 max(zData)])
axis square
% close
set(gca,'YDir','reverse');
rotate3d on
colorbar('north','TickLabelInterpreter','latex','AxisLocation','out');
    fig_setup
    legend off
    txt = {'(b)'};
    b=text(4,0.5,txt);
    b.Units='normalized';


level=round(logspace(log10(min(min((zData)))*10^6),log10(0.9*max(max((zData)))*10^6),8)./10^6,2);
h(4) = figure;
contour3(tmp_u,tmp_r,...
    ((co.b(1).*tmp_r.^co.eb(1)+co.b(2)).*tmp_u.^0)+...
    ((co.b(3).*tmp_r.^co.eb(2)+co.b(4)).*tmp_u.^1)+...
    ((co.b(5).*tmp_r.^co.eb(3)+co.b(6)).*tmp_u.^2),level,'-','color','k','LineWidth',2,'ShowText','on')
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
% zlabel('$D2 (x,r) / [\frac{m^2}{s^3}]$', 'interpreter','latex')
zlabel('$D^{(2)} (u_r,r)$', 'interpreter','latex')
grid off
box off
xlim([min(xData) max(xData)])
ylim([min(y(:,1)) max(y(:,1))])
% zlim([0 2*max(zData)])
caxis([0 max(zData)])
axis square
% close
set(gca,'YDir','reverse');
rotate3d on



h(5) = figure;
surf(tmp_u,tmp_r,((co.a(1).*tmp_r.^co.ea(1)+co.a(2)).*tmp_u)./...
    (((co.b(1).*tmp_r.^co.eb(1)+co.b(2)).*tmp_u.^0)+...
    ((co.b(3).*tmp_r.^co.eb(2)+co.b(4)).*tmp_u.^1)+...
    ((co.b(5).*tmp_r.^co.eb(3)+co.b(6)).*tmp_u.^2)),'EdgeColor','none')
hold on
% surf(tmp_u,tmp_r,...
% (d10_R+d11_R.*tmp_u+d12_R.*tmp_u.^2+d13_R.*tmp_u.^3)./(d20_R+d21_R.*tmp_u+d22_R.*tmp_u.^2),'FaceColor',[0.831372559070587 0.815686285495758 0.7843137383461]);

% surf(tmp_u,tmp_r,((-0.36.*tmp_r.^(-1)).*tmp_u)./...
%     ((0.0144.*tmp_r.^(-1)).*tmp_u.^2),'FaceColor',[0.831372559070587 0.815686285495758 0.7843137383461])
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
% zlabel('$D2 (x,r) / [\frac{m^2}{s^3}]$', 'interpreter','latex')
zlabel('$D^{(1)} (u_r,r)/D^{(2)} (u_r,r)$', 'interpreter','latex')
colormap (parula(200))
% colorbar
grid off
xlim([min(xData) max(xData)])
ylim([min(y(:,1)) max(y(:,1))])
% zlim([-10 +10])
% caxis([-7 +7])
axis square
% close
set(gca,'YDir','reverse');
rotate3d on
colorbar('north','TickLabelInterpreter','latex','AxisLocation','out');
font=22;
    a.FontSize = font;
    b.FontSize = font;
    pos_txt=[-0.22   0.8];
    a.Position=pos_txt;
    b.Position=[-0.25  0.8];

if ischar(save_path)
    savefig(h,fullfile(save_path,append(save_name,'_','KM_opti.fig')),'compact')
    for a = 1:length(h)
        exportgraphics(h(a),fullfile(save_path,append(save_name,'_',sprintf('KM_opti_%d.png', a))))
    end
end
% set(gca,'ColorScale','log')

for i=1:5
    tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
    uiwait(gcf);
    delete(tmp_ui);
    close
end
end