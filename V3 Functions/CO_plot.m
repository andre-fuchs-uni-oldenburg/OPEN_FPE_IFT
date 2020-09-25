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
function CO_plot(co_a,co_b,evaluated,taylor_L,varargin)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');

aufl    = size(evaluated(1).D1,2);
y       = nan(length(evaluated),aufl);
for i=1:length(evaluated)
    y(i,1:length(evaluated(i).x_bin_not_nan))=(evaluated(i).r./taylor_L);
end

if isfield(co_b,'a') && ~isempty(co_b.a) 
    if nargin >= 5
        co_c       = varargin{1};
    end

%% d_ij(r) plots
    figure
    subplot(2,2,1);
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        co_a.a(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_a.ea(1)+co_a.a(2),'LineWidth',2)
    hold on
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        co_b.a(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_b.ea(1)+co_b.a(2),'LineWidth',2)
    if nargin >= 5
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        co_c.a(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_c.ea(1)+co_c.a(2),'k','LineWidth',2)
    end
%     plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%         (-0.36).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^(-1),'LineWidth',2,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])

    %Renner
    % plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
    %     (-0.61).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^(-0.67),'g','LineWidth',2);

%     x_fit=logspace(log10(min(y(:,1))),log10(max(y(:,1))),100);
%     y_fit=co_a.a(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_a.ea(1)+co_a.a(2);
%     [xData, yData] = prepareCurveData( x_fit, y_fit );
%     ft = fittype( 'power1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     [fitresult_d11, gof] = fit( xData, yData, ft, opts );
%     plot(linspace(0,max(y(:,1)),100),feval(fitresult_d11,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);

    axis square
    ylabel('$d_{11}(r)$', 'interpreter','latex')
    xlabel('$r / \lambda$', 'interpreter','latex')
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',24)
    % set(gca,'xScale','log') 

    
    subplot(2,2,2);
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        co_a.b(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_a.eb(1)+co_a.b(2),'LineWidth',2)
    hold on
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        co_b.b(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_b.eb(1)+co_b.b(2),'LineWidth',2)
    if nargin >= 5
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        co_c.b(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_c.eb(1)+co_c.b(2),'k','LineWidth',2)
    end
    %Renner
    % plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
    %     (0.033).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^(0.25),'g','LineWidth',2);

%     y_fit=co_a.b(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_a.eb(1)+co_a.b(2);
%     [xData, yData] = prepareCurveData( x_fit, y_fit );
%     ft = fittype( 'power1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     [fitresult_d20, gof] = fit( xData, yData, ft, opts );
%     plot(linspace(0,max(y(:,1)),100),feval(fitresult_d20,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);

    axis square
    ylabel('$d_{20}(r)$', 'interpreter','latex')
    xlabel('$r / \lambda$', 'interpreter','latex')
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',24)
    % set(gca,'xScale','log') 
    
    
    subplot(2,2,3);
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        co_a.b(3).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_a.eb(2)+co_a.b(4),'LineWidth',2)
    hold on
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        co_b.b(3).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_b.eb(2)+co_b.b(4),'LineWidth',2)
    if nargin >= 5
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        co_c.b(3).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_c.eb(2)+co_c.b(4),'k','LineWidth',2)
    end
    %Renner
    % plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
    %     (-0.009).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^(0.2),'g','LineWidth',2);

%     y_fit=co_a.b(3).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_a.eb(2)+co_a.b(4);
%     [xData, yData] = prepareCurveData( x_fit, y_fit );
%     ft = fittype( 'power1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     [fitresult_d21, gof] = fit( xData, yData, ft, opts );
%     plot(linspace(0,max(y(:,1)),100),feval(fitresult_d21,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);

    axis square
    ylabel('$d_{21}(r)$', 'interpreter','latex')
    xlabel('$r / \lambda$', 'interpreter','latex')
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',24)
    % set(gca,'xScale','log') 
    
    
    subplot(2,2,4);
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        co_a.b(5).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_a.eb(3)+co_a.b(6),'LineWidth',2)
    hold on
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        co_b.b(5).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_b.eb(3)+co_b.b(6),'LineWidth',2)
    if nargin >= 5
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        co_c.b(5).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_c.eb(3)+co_c.b(6),'k','LineWidth',2)
    end
%     plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%         0.0144.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^(-1),'LineWidth',2,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
    %Renner
    % plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
    %     (0.043).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^(-0.73),'g','LineWidth',2);

%     y_fit=co_a.b(5).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_a.eb(3)+co_a.b(6);
%     [xData, yData] = prepareCurveData( x_fit, y_fit );
%     ft = fittype( 'power1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     [fitresult_d22, gof] = fit( xData, yData, ft, opts );
%     plot(linspace(0,max(y(:,1)),100),feval(fitresult_d22,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);

    axis square
    ylabel('$d_{22}(r)$', 'interpreter','latex')
    xlabel('$r / \lambda$', 'interpreter','latex')
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',24)
    set(gcf, 'Color', 'w')
    % set(gca,'xScale','log') 
    legend({'Optimized','Non-optimized','K62'},'Position',[0.45 0.5 0.1 0.05])
    if nargin >= 5
        legend({'IFT-optimized','con PDF-optimized','Non-optimized','K62'},'Position',[0.45 0.5 0.1 0.05],'FontSize',18)
    end

%     fitresult_d11
%     fitresult_d20
%     fitresult_d21
%     fitresult_d22
       
else
%% d_ij(r) plots
    figure
    subplot(2,2,1);
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        co_a.a(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_a.ea(1)+co_a.a(2),'LineWidth',2)
    hold on
%     plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%         (-0.36).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^(-1),'LineWidth',2,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
    %Renner
    % plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
    %     (-0.61).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^(-0.67),'g','LineWidth',2);

%     x_fit=logspace(log10(min(y(:,1))),log10(max(y(:,1))),100);
%     y_fit=co_a.a(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_a.ea(1)+co_a.a(2);
%     [xData, yData] = prepareCurveData( x_fit, y_fit );
%     ft = fittype( 'power1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     [fitresult_d11, gof] = fit( xData, yData, ft, opts );
%   plot(linspace(0,max(y(:,1)),100),feval(fitresult_d11,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);

    axis square
    ylabel('$d_{11}(r)$', 'interpreter','latex')
    xlabel('$r / \lambda$', 'interpreter','latex')
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',24)
    % set(gca,'xScale','log') 

    
    subplot(2,2,2);
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        co_a.b(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_a.eb(1)+co_a.b(2),'LineWidth',2)
    hold on
    %Renner
    % plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
    %     (0.033).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^(0.25),'g','LineWidth',2);

%     y_fit=co_a.b(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_a.eb(1)+co_a.b(2);
%     [xData, yData] = prepareCurveData( x_fit, y_fit );
%     ft = fittype( 'power1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     [fitresult_d20, gof] = fit( xData, yData, ft, opts );
%   plot(linspace(0,max(y(:,1)),100),feval(fitresult_d20,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);

    axis square
    ylabel('$d_{20}(r)$', 'interpreter','latex')
    xlabel('$r / \lambda$', 'interpreter','latex')
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',24)
    % set(gca,'xScale','log') 

    
    subplot(2,2,3);
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        co_a.b(3).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_a.eb(2)+co_a.b(4),'LineWidth',2)
    hold on
    %Renner
    % plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
    %     (-0.009).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^(0.2),'g','LineWidth',2);
    
%     y_fit=co_a.b(3).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_a.eb(2)+co_a.b(4);
%     [xData, yData] = prepareCurveData( x_fit, y_fit );
%     ft = fittype( 'power1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     [fitresult_d21, gof] = fit( xData, yData, ft, opts );
%     plot(linspace(0,max(y(:,1)),100),feval(fitresult_d21,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);

    axis square
    ylabel('$d_{21}(r)$', 'interpreter','latex')
    xlabel('$r / \lambda$', 'interpreter','latex')
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',24)
    % set(gca,'xScale','log') 

    
    subplot(2,2,4);
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        co_a.b(5).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_a.eb(3)+co_a.b(6),'LineWidth',2)
    hold on
%   plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%         0.0144.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^(-1),'LineWidth',2,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
    %Renner
    % plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
    %     (0.043).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^(-0.73),'g','LineWidth',2);

%     y_fit=co_a.b(5).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_a.eb(3)+co_a.b(6);
%     [xData, yData] = prepareCurveData( x_fit, y_fit );
%     ft = fittype( 'power1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     [fitresult_d22, gof] = fit( xData, yData, ft, opts );
%    plot(linspace(0,max(y(:,1)),100),feval(fitresult_d22,linspace(0,max(y(:,1)),100)),'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176]);

    axis square
    ylabel('$d_{22}(r)$', 'interpreter','latex')
    xlabel('$r / \lambda$', 'interpreter','latex')
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',24)
    set(gcf, 'Color', 'w')
    % set(gca,'xScale','log') 
    legend({'Optimized','K62'},'Position',[0.45 0.5 0.1 0.05])

%     fitresult_d11
%     fitresult_d20
%     fitresult_d21
%     fitresult_d22
end

% figure
% plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%     co_b.b(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_b.eb(1))
% hold on
% plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),co_b.b(2)*ones(1,100));
% 
% figure
% plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%     co_a.b(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_a.eb(1))
% hold on
% plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),co_a.b(2)*ones(1,100));
end