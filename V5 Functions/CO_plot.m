%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots the parameters $d_{11}$, $d_{20}$, $d_{21}$ and $d_{22}$ as a function of $r$ 
% for optimized and non-optimized $D^{(1)}$ and $D^{(2)}$.
% Arguments IN
% co = Coefficients $d_{ij}(r)$ of surface fits with a linear function for 
%      $D^{(1)}\left(u_r,r\right)$ and a parabolic function for $D^{(2)}\left(u_r,r\right)$
% co_a = co_KM_opti
% co_b = co_KM_non_opti
% evaluated = evaluated struct array from the function 'KM_Calculation'
% taylor_L = Taylor length scale in meters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CO_plot(conv,co_a,co_b,evaluated,taylor_L,int_L,norm_ur,norm_r,siginf,varargin)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');

aufl    = size(evaluated(1).D1,2);
y       = nan(length(evaluated),aufl);
if norm_r==1
    for i=1:length(evaluated)
        y(i,1:length(evaluated(i).x_bin_not_nan))=(evaluated(i).r./taylor_L);
    end
else
    for i=1:length(evaluated)
        y(i,1:length(evaluated(i).x_bin_not_nan))=(evaluated(i).r);
    end  
end

r_tmp=logspace(log10(min(y(:,1))),log10(max(y(:,1))),100);
    
if conv==1
    conv_r=1;
else
    conv_r=0;
end

if isfield(co_b,'a') && ~isempty(co_b.a) 
    if nargin >= 10
        co_c       = varargin{1};
    end

%% d_ij(r) plots
    figure
    subplot(2,2,1);
    plot(r_tmp,...
        (co_a.a(1).*r_tmp.^co_a.ea(1)+co_a.a(2)).*r_tmp.^(conv_r),'LineWidth',2)
    hold on
    plot(r_tmp,...
        (co_b.a(1).*r_tmp.^co_b.ea(1)+co_b.a(2)).*r_tmp.^(conv_r),'LineWidth',2)
    if nargin >= 10
    plot(r_tmp,...
        (co_c.a(1).*r_tmp.^co_c.ea(1)+co_c.a(2)).*r_tmp.^(conv_r),'k','LineWidth',2)
    end
    plot(r_tmp,...
        (-0.362)*r_tmp.^(conv-1),'LineWidth',2,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
    plot(r_tmp,...
        (-0.362)*r_tmp.^(conv_r-1)/siginf,'LineWidth',2,'Color','r')
%     plot(r_tmp,...
%         (-0.362).*logspace(log10(min(y(:,1).*taylor_L)),log10(max(y(:,1))),100).^(0),'LineWidth',2,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
%    plot(r_tmp,...
%         (-1/3).*r_tmp.^(-1),'LineWidth',2,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])


    %Renner
%     plot(r_tmp,...
%         (-0.61).*r_tmp.^(-0.67).*r_tmp.^(conv_r),'g','LineWidth',2);

%     x_fit=r_tmp;
%     y_fit=co_a.a(1).*r_tmp.^co_a.ea(1)+co_a.a(2);
%     [xData, yData] = prepareCurveData( x_fit, y_fit );
%     ft = fittype( 'power1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     [fitresult_d11, gof] = fit( xData, yData, ft, opts );
%     plot(linspace(0,max(y(:,1)),100),feval(fitresult_d11,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);

    axis square
%     ylabel('$d_{11}(r)$', 'interpreter','latex')
    if conv==1
        ylabel('$d_{11} * (r / \lambda)$', 'interpreter','latex')
    else
        ylabel('$d_{11}$', 'interpreter','latex')
    end
    if norm_r==1
        xlabel('$r / \lambda$', 'interpreter','latex')
    else
        xlabel('$r\ (m)$', 'interpreter','latex')
    end
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',14)
    % set(gca,'xScale','log') 
    fig_setup
        legend off
        txt = {'(a)'};
        a=text(4,0.5,txt);
        a.Units='normalized';
    

    subplot(2,2,2);
    plot(r_tmp,...
        (co_a.b(1).*r_tmp.^co_a.eb(1)+co_a.b(2)).*r_tmp.^(conv_r),'LineWidth',2)
    hold on
    plot(r_tmp,...
        (co_b.b(1).*r_tmp.^co_b.eb(1)+co_b.b(2)).*r_tmp.^(conv_r),'LineWidth',2)
    if nargin >= 10
    plot(r_tmp,...
        (co_c.b(1).*r_tmp.^co_c.eb(1)+co_c.b(2)).*r_tmp.^(conv_r),'k','LineWidth',2)
    end
    %Renner
%     plot(r_tmp,...
%         (0.033).*r_tmp.^(0.25).*r_tmp.^(conv_r),'g','LineWidth',2);

%     y_fit=co_a.b(1).*r_tmp.^co_a.eb(1)+co_a.b(2);
%     [xData, yData] = prepareCurveData( x_fit, y_fit );
%     ft = fittype( 'power1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     [fitresult_d20, gof] = fit( xData, yData, ft, opts );
%     plot(linspace(0,max(y(:,1)),100),feval(fitresult_d20,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);

    axis square
%     ylabel('$d_{20}(r)$', 'interpreter','latex')
    if conv==1
    ylabel('$d_{20} * (r / \lambda)$', 'interpreter','latex')
    else
        ylabel('$d_{20}$', 'interpreter','latex')
    end
    if norm_r==1
        xlabel('$r / \lambda$', 'interpreter','latex')
    else
        xlabel('$r\ (m)$', 'interpreter','latex')
    end
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',14)
    % set(gca,'xScale','log') 
    fig_setup
    legend off
    txt = {'(b)'};
    b=text(4,0.5,txt);
    b.Units='normalized';
    
    subplot(2,2,3);
    plot(r_tmp,...
        (co_a.b(3).*r_tmp.^co_a.eb(2)+co_a.b(4)).*r_tmp.^(conv_r),'LineWidth',2)
    hold on
    plot(r_tmp,...
        (co_b.b(3).*r_tmp.^co_b.eb(2)+co_b.b(4)).*r_tmp.^(conv_r),'LineWidth',2)
    if nargin >= 10
    plot(r_tmp,...
        (co_c.b(3).*r_tmp.^co_c.eb(2)+co_c.b(4)).*r_tmp.^(conv_r),'k','LineWidth',2)
    end
    %Renner
%     plot(r_tmp,...
%         (-0.009).*r_tmp.^(0.2).*r_tmp.^(conv_r),'g','LineWidth',2);

%     y_fit=co_a.b(3).*r_tmp.^co_a.eb(2)+co_a.b(4);
%     [xData, yData] = prepareCurveData( x_fit, y_fit );
%     ft = fittype( 'power1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     [fitresult_d21, gof] = fit( xData, yData, ft, opts );
%     plot(linspace(0,max(y(:,1)),100),feval(fitresult_d21,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);

    axis square
%     ylabel('$d_{21}(r)$', 'interpreter','latex')
    if conv==1
        ylabel('$d_{21} * (r / \lambda)$', 'interpreter','latex')
    else
        ylabel('$d_{21}$', 'interpreter','latex')
    end
    if norm_r==1
        xlabel('$r / \lambda$', 'interpreter','latex')
    else
        xlabel('$r\ (m)$', 'interpreter','latex')
    end
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',14)
    % set(gca,'xScale','log') 
    fig_setup
    legend off
    txt = {'(c)'};
    c=text(4,0.5,txt);
    c.Units='normalized';
    
    subplot(2,2,4);
    plot(r_tmp,...
        (co_a.b(5).*r_tmp.^co_a.eb(3)+co_a.b(6)).*r_tmp.^(conv_r),'LineWidth',2)
    hold on
    plot(r_tmp,...
        (co_b.b(5).*r_tmp.^co_b.eb(3)+co_b.b(6)).*r_tmp.^(conv_r),'LineWidth',2)
    if nargin >= 10
    plot(r_tmp,...
        (co_c.b(5).*r_tmp.^co_c.eb(3)+co_c.b(6)).*r_tmp.^(conv_r),'k','LineWidth',2)
    end
    plot(r_tmp,...
          0.0144.*r_tmp.^(conv_r-1),'LineWidth',2,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
      plot(r_tmp,...
        0.0144.*r_tmp.^(conv_r-1)/(siginf^2),'LineWidth',2,'Color','r')
%      plot(r_tmp,...
%         0.0144.*r_tmp.^(0),'LineWidth',2,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])

    %Renner
%     plot(r_tmp,...
%         (0.043).*r_tmp.^(-0.73).*r_tmp.^(conv_r),'g','LineWidth',2);

%     y_fit=co_a.b(5).*r_tmp.^co_a.eb(3)+co_a.b(6);
%     [xData, yData] = prepareCurveData( x_fit, y_fit );
%     ft = fittype( 'power1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     [fitresult_d22, gof] = fit( xData, yData, ft, opts );
%     plot(linspace(0,max(y(:,1)),100),feval(fitresult_d22,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);

    axis square
%     ylabel('$d_{22}(r)$', 'interpreter','latex')
    if conv==1
        ylabel('$d_{22} * (r / \lambda)$', 'interpreter','latex')
    else
        ylabel('$d_{22}$', 'interpreter','latex')
    end
    if norm_r==1
        xlabel('$r / \lambda$', 'interpreter','latex')
    else
        xlabel('$r\ (m)$', 'interpreter','latex')
    end
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',14)
    set(gcf, 'Color', 'w')
    % set(gca,'xScale','log') 
    legend({'Optimized','Non-optimized','K62'},'Position',[0.45 0.5 0.1 0.05])
    if nargin >= 10
        legend({'IFT-optimized','con PDF-optimized','Non-optimized','K62'},'Position',[0.45 0.5 0.1 0.05],'FontSize',18)
    end
        fig_setup
% legend off
    txt = {'(d)'};
    d=text(4,0.5,txt);
    d.Units='normalized';

%     fitresult_d11
%     fitresult_d20
%     fitresult_d21
%     fitresult_d22
       
else
%% d_ij(r) plots
    figure
    subplot(2,2,1);
    plot(r_tmp,...
       (co_a.a(1).*r_tmp.^co_a.ea(1)+co_a.a(2)).*r_tmp.^(conv_r),'LineWidth',2)
%         co_a.a(1).*r_tmp.^co_a.ea(1)+co_a.a(2),'LineWidth',2)
    hold on
%     plot(r_tmp,...
%         (-0.36).*r_tmp.^(-1),'LineWidth',2,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
    %Renner
    % plot(r_tmp,...
    %     (-0.61).*r_tmp.^(-0.67),'g','LineWidth',2);

%     x_fit=r_tmp;
%     y_fit=co_a.a(1).*r_tmp.^co_a.ea(1)+co_a.a(2);
%     [xData, yData] = prepareCurveData( x_fit, y_fit );
%     ft = fittype( 'power1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     [fitresult_d11, gof] = fit( xData, yData, ft, opts );
%   plot(linspace(0,max(y(:,1)),100),feval(fitresult_d11,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);

    axis square
%     ylabel('$d_{11}(r)$', 'interpreter','latex')
    if conv==1
    ylabel('$d_{11} * (r / \lambda)$', 'interpreter','latex')
    else
    ylabel('$d_{11}$', 'interpreter','latex')
    end
    if norm_r==1
        xlabel('$r / \lambda$', 'interpreter','latex')
    else
        xlabel('$r\ (m)$', 'interpreter','latex')
    end
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',14)
    % set(gca,'xScale','log') 
fig_setup
    legend off
    txt = {'(a)'};
    a=text(4,0.5,txt);
    a.Units='normalized';
    
    subplot(2,2,2);
    plot(r_tmp,...
       (co_a.b(1).*r_tmp.^co_a.eb(1)+co_a.b(2)).*r_tmp.^(conv_r),'LineWidth',2)
%         co_a.b(1).*r_tmp.^co_a.eb(1)+co_a.b(2),'LineWidth',2)
    hold on
    %Renner
    % plot(r_tmp,...
    %     (0.033).*r_tmp.^(0.25),'g','LineWidth',2);

%     y_fit=co_a.b(1).*r_tmp.^co_a.eb(1)+co_a.b(2);
%     [xData, yData] = prepareCurveData( x_fit, y_fit );
%     ft = fittype( 'power1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     [fitresult_d20, gof] = fit( xData, yData, ft, opts );
%   plot(linspace(0,max(y(:,1)),100),feval(fitresult_d20,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);

    axis square
%     ylabel('$d_{20}(r)$', 'interpreter','latex')
%     ylabel('$d_{20}$', 'interpreter','latex')
    if conv==1
        ylabel('$d_{20} * (r / \lambda)$', 'interpreter','latex')
    else
        ylabel('$d_{20}$', 'interpreter','latex')
    end
    if norm_r==1
        xlabel('$r / \lambda$', 'interpreter','latex')
    else
        xlabel('$r\ (m)$', 'interpreter','latex')
    end
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',14)
    % set(gca,'xScale','log') 
fig_setup
    legend off
    txt = {'(b)'};
    b=text(4,0.5,txt);
    b.Units='normalized';
    
    subplot(2,2,3);
    plot(r_tmp,...
       (co_a.b(3).*r_tmp.^co_a.eb(2)+co_a.b(4)).*r_tmp.^(conv_r),'LineWidth',2)
%         co_a.b(3).*r_tmp.^co_a.eb(2)+co_a.b(4),'LineWidth',2)
    hold on
    %Renner
    % plot(r_tmp,...
    %     (-0.009).*r_tmp.^(0.2),'g','LineWidth',2);
    
%     y_fit=co_a.b(3).*r_tmp.^co_a.eb(2)+co_a.b(4);
%     [xData, yData] = prepareCurveData( x_fit, y_fit );
%     ft = fittype( 'power1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     [fitresult_d21, gof] = fit( xData, yData, ft, opts );
%     plot(linspace(0,max(y(:,1)),100),feval(fitresult_d21,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);

    axis square
%     ylabel('$d_{21}(r)$', 'interpreter','latex')
%     ylabel('$d_{21}$', 'interpreter','latex')
    if conv==1
        ylabel('$d_{21} * (r / \lambda)$', 'interpreter','latex')
    else
        ylabel('$d_{21}$', 'interpreter','latex')
    end
    if norm_r==1
        xlabel('$r / \lambda$', 'interpreter','latex')
    else
        xlabel('$r\ (m)$', 'interpreter','latex')
    end
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',14)
    % set(gca,'xScale','log') 
fig_setup
    legend off
    txt = {'(c)'};
    c=text(4,0.5,txt);
    c.Units='normalized';
    
    subplot(2,2,4);
    plot(r_tmp,...
       (co_a.b(5).*r_tmp.^co_a.eb(3)+co_a.b(6)).*r_tmp.^(conv_r),'LineWidth',2)
%         co_a.b(5).*r_tmp.^co_a.eb(3)+co_a.b(6),'LineWidth',2)
    hold on
%   plot(r_tmp,...
%         0.0144.*r_tmp.^(-1),'LineWidth',2,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
    %Renner
    % plot(r_tmp,...
    %     (0.043).*r_tmp.^(-0.73),'g','LineWidth',2);

%     y_fit=co_a.b(5).*r_tmp.^co_a.eb(3)+co_a.b(6);
%     [xData, yData] = prepareCurveData( x_fit, y_fit );
%     ft = fittype( 'power1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     [fitresult_d22, gof] = fit( xData, yData, ft, opts );
%    plot(linspace(0,max(y(:,1)),100),feval(fitresult_d22,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);

    axis square
%     ylabel('$d_{22}(r)$', 'interpreter','latex')
%     ylabel('$d_{22}$', 'interpreter','latex')
    if conv==1
        ylabel('$d_{22} * (r / \lambda)$', 'interpreter','latex')
    else
        ylabel('$d_{22}$', 'interpreter','latex')
    end
    if norm_r==1
        xlabel('$r / \lambda$', 'interpreter','latex')
    else
        xlabel('$r\ (m)$', 'interpreter','latex')
    end
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',14)
    set(gcf, 'Color', 'w')
    % set(gca,'xScale','log') 
    legend({'Optimized','K62'},'Position',[0.45 0.5 0.1 0.05])
    fig_setup
% legend off
    txt = {'(d)'};
    d=text(4,0.5,txt);
    d.Units='normalized';
%     fitresult_d11
%     fitresult_d20
%     fitresult_d21
%     fitresult_d22
end

% figure
% plot(r_tmp,...
%     co_b.b(1).*r_tmp.^co_b.eb(1))
% hold on
% plot(r_tmp,co_b.b(2)*ones(1,100));
% 
% figure
% plot(r_tmp,...
%     co_a.b(1).*r_tmp.^co_a.eb(1))
% hold on
% plot(r_tmp,co_a.b(2)*ones(1,100));

 font=22;
    a.FontSize = font;
    b.FontSize = font;
    c.FontSize = font;
    d.FontSize = font;

    pos_txt=[-0.68   0.9];
    a.Position=pos_txt;
    b.Position=pos_txt;
    c.Position=pos_txt;
    d.Position=pos_txt;

for i=1
    tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
    uiwait(gcf);
    delete(tmp_ui);
    close
end
end