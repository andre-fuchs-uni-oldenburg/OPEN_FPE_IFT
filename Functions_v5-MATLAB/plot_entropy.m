%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots the empirical average $\langle e^{\mathrm{-}\Delta S_{tot}} \rangle_N$ of 
% $\Delta S_{tot}$ as a function of the number, $N$ (sample size), of cascade trajectories 
% $\left[u(\cdot) \right]$ with errorbars. In addtion, the probability density function of the system, 
% medium and total entropy will be plotted while displaying the value of $\langle\Delta S_{tot}\rangle$ 
% which should be $\geq$ 0. The  integral fluctuation theorem (IFT) expresses the integral balance 
% between the entropy-consuming ($\Delta S_{tot}<0$) and the entropy-producing ($\Delta S_{tot}>0$) 
% cascade trajectories
% 
% Arguments IN
% Sm = Entropy of the medium
% Ds = Entropy of the system or the Shanon entropy
% DS = Total entropy production=Sm+Ds
% increment_bin = number of bins
% save_path = path for saving figures and files
% save_name = name for saving files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_entropy(Sm,Ds,DS,increment_bin,save_path,save_name)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');


N               = size(DS,1);   % number of trajectories u(r)
NN              = 20;           % steps for the calculation of mean and std value
Nvec            = round(logspace(1,log10(N),NN));
FTvec           = nan(1,NN);
FTvec_err       = nan(1,NN);
FTvec_med       = nan(1,NN);
FTvec_err_med   = nan(1,NN);
FTvec_sys       = nan(1,NN);
FTvec_err_sys   = nan(1,NN);

for n=1:NN
    FTvec(n)            = nanmean( exp( -DS(1:Nvec(n)) ) );
%   FTvec(n)            = median( exp( -DS(1:Nvec(n)) ) );
    FTvec_err(n)        = nanstd( exp( -DS(1:Nvec(n)) ) )/sqrt(Nvec(n));
    
    FTvec_med(n)        = nanmean( exp( -Sm(1:Nvec(n)) ) );
    FTvec_err_med(n)    = nanstd( exp( -Sm(1:Nvec(n)) ) )/sqrt(Nvec(n));
    
    FTvec_sys(n)        = nanmean( exp( -Ds(1:Nvec(n)) ) );
    FTvec_err_sys(n)    = nanstd( exp( -Ds(1:Nvec(n)) ) )/sqrt(Nvec(n));
end


h(1) = figure;
errorbar(Nvec,FTvec,FTvec_err,'MarkerSize',8,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',3,'LineStyle','none','Color',[0 0 0])
hold on
% According to the fluctuation theorem the average has to converge to the horizontal line (=1).
plot(Nvec,ones(numel(Nvec),1),'k','LineStyle','--','LineWidth',2);
axis square      
set(gca, 'FontSize',18)
set(gcf, 'Color', 'w')
xlabel('$N$','interpreter','latex')
ylabel('$\langle e^{\mathrm{-}\Delta S_{tot}} \rangle_N$','interpreter','latex')
set(gca,'xScale','log') 
% set(gca,'yScale','log') 
xlim([min(Nvec) max(Nvec)])
% set(gca,'XTick',[10 1000 100000])
% ylim([0.1 1.1])
title(['$\langle e^{\mathrm{-}\Delta S_{tot}} \rangle_{max(N)}=$',num2str(nanmean( exp( -DS )),'%1.2f')],'interpreter','latex')
fig_setup
legend off
txt = {'(a)'};
a=text(4,0.5,txt);
a.Units='normalized';


h(2) = figure;
subplot(1,2,1)
errorbar(Nvec,abs(1-FTvec),FTvec_err,'MarkerSize',8,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',3,'LineStyle','none','Color',[0 0 0])
% plot(Nvec,abs(1-FTvec),'MarkerSize',8,'Marker','*','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',4,'LineStyle','none')
hold on 
[xData, yData, weights] = prepareCurveData( Nvec, abs(1-FTvec), FTvec_err );
% Set up fittype and options.
ft = fittype( 'power1' );
ft_2 = fittype( 'a*x^-0.5', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Weights = weights;
% [fit_IFT, gof] = fit( xData, yData, ft, opts );
% plot(Nvec,fit_IFT.a*Nvec.^fit_IFT.b,'k','LineWidth',4)

[fit_IFT_2, gof] = fit( xData, yData, ft_2, opts );
plot(Nvec,fit_IFT_2.a*Nvec.^-0.5,'r','LineWidth',4)
% plot(Nvec,abs(1-FTvec),'MarkerSize',8,'Marker','*','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',4,'LineStyle','none')
axis square      
set(gca, 'FontSize',18)
set(gcf, 'Color', 'w')
xlabel('$N$','interpreter','latex')
ylabel('$|1-\langle e^{\mathrm{-}\Delta S_{tot}} \rangle_N|$','interpreter','latex')
set(gca,'xScale','log') 
% set(gca,'yScale','log') 
xlim([min(Nvec) max(Nvec)])
% set(gca,'XTick',[10 1000 100000])
% ylim([-0.3 1.2])
% title(['$\langle e^{\mathrm{-}\Delta S_{tot}} \rangle_{max(N)}=$',num2str(nanmean( exp( -DS )),'%1.4f')],'interpreter','latex')
plot(Nvec,0.*ones(numel(Nvec),1),'k','LineStyle','--','LineWidth',2);
legend({'exp','$aN^{-0.5}$'},'Interpreter','latex','Location','northeast')


% figure
subplot(1,2,2)
% errorbar(Nvec,abs(1-FTvec),FTvec_err,'MarkerSize',8,'Marker','*','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',4,'LineStyle','none')
plot(Nvec,abs(1-FTvec),'MarkerSize',8,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',3,'LineStyle','none')
hold on 
[xData, yData, weights] = prepareCurveData( Nvec, abs(1-FTvec), FTvec_err );
% Set up fittype and options.
ft = fittype( 'power1' );
ft_2 = fittype( 'a*x^-0.5', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Weights = weights;
% [fit_IFT, gof] = fit( xData, yData, ft, opts );
% plot(Nvec,fit_IFT.a*Nvec.^fit_IFT.b,'k','LineWidth',4)

[fit_IFT_2, gof] = fit( xData, yData, ft_2, opts );
plot(Nvec,fit_IFT_2.a*Nvec.^-0.5,'r','LineWidth',4)
% plot(Nvec,abs(1-FTvec),'MarkerSize',8,'Marker','*','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',4,'LineStyle','none')
axis square      
set(gca, 'FontSize',18)
set(gcf, 'Color', 'w')
xlabel('$N$','interpreter','latex')
ylabel('$|1-\langle e^{\mathrm{-}\Delta S_{tot}} \rangle_N|$','interpreter','latex')
set(gca,'xScale','log') 
set(gca,'yScale','log') 
xlim([min(Nvec) max(Nvec)])
% set(gca,'XTick',[10 1000 100000])
% ylim([-0.3 1.2])
% title(['$\langle e^{\mathrm{-}\Delta S_{tot}} \rangle_{max(N)}=$',num2str(nanmean( exp( -DS )),'%1.4f')],'interpreter','latex')
legend({'exp','$aN^{-0.5}$'},'Interpreter','latex','Location','northeast')


hist_min = min(min([Sm Ds DS]));
hist_max = max(max([Sm Ds DS]));
hist_tmp = linspace(hist_min,hist_max,increment_bin);

% [fuL uL] = hist(Sm,501); 
[fuL uL] = hist(Sm,hist_tmp); 
puL = fuL / ( sum(fuL) * mean(diff(uL)) ); 

h(3) = figure;
plot(uL,puL,'-','LineWidth',2)
% uL(find(puL==max(puL)))
% max(puL)
hold on

% [fuL uL] = hist(Ds,501); 
[fuL uL] = hist(Ds,hist_tmp); 
puL = fuL / ( sum(fuL) * mean(diff(uL)) ); 
plot(uL,puL,'-','LineWidth',2)
% uL(find(puL==max(puL)))
% max(puL)

% [fuL uL] = hist(DS,501); 
[fuL uL] = hist(DS,hist_tmp); 
puL = fuL / ( sum(fuL) * mean(diff(uL)) ); 
plot(uL,puL,'-','LineWidth',2,'Color','k')
% plot(uL,puL,'MarkerSize',8,'MarkerEdgeColor',[0 0 0],'Marker','o','LineWidth',1.5,'LineStyle','none','Color','k')
% uL(find(puL==max(puL)))
% max(puL)
axis square
set(gca, 'FontSize',18)
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
legend('$\Delta S_{med}$','$\Delta S_{sys}$','$\Delta S_{tot}$')
xlabel('$Entropy$','interpreter','latex')
ylabel('$PDF$','interpreter','latex')
title(['$\langle \Delta S_{tot}\rangle=$',num2str(nanmean(DS),'%1.2f')],'interpreter','latex')
% xlabel('$S_{med}, S_{sys}, \Delta S_{tot}$','interpreter','latex')
% ylabel('$pdf(S)$','interpreter','latex')
xlim([min(hist_tmp) max(hist_tmp)])
% set(gca,'yScale','lin') 
% xlim([-11 16])
fig_setup
vline(0,'k')
txt = {'(b)'};
b=text(4,0.5,txt);
b.Units='normalized';
font=22;
a.FontSize = font;
b.FontSize = font;
pos_txt=[-0.25   0.9];
a.Position=[-0.22   0.9];
b.Position=pos_txt;


h(4) = figure;
semilogy(uL,puL,'MarkerSize',8,'MarkerEdgeColor',[0 0 0],'Marker','o','LineWidth',1.5,'LineStyle','none')
hold on
ylabel('$p(\Delta S_{tot})$');
ylabel('$PDF$','interpreter','latex')
xlabel('$\Delta S_{tot}$');
axis square
set(gca, 'FontSize',18)
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
title(['$\langle \Delta S_{tot}\rangle=$',num2str(nanmean(DS),'%1.2f')],'interpreter','latex')
% hold on
% pd      = fitdist(DS,'GeneralizedExtremeValue')
% plot(uL,pdf(pd,uL),'Color',[0.501960813999176 0.501960813999176 0.501960813999176],'LineWidth',2)
% ylim([10^-6 2*max(puL)])
fig_setup
legend off
vline(0,'k')



h(5) = figure;
errorbar(Nvec,FTvec_med,FTvec_err_med,'MarkerSize',8,'Marker','o','LineWidth',3)
hold on
errorbar(Nvec,FTvec_sys,FTvec_err_sys,'MarkerSize',8,'Marker','o','LineWidth',3)
axis square      
set(gca, 'FontSize',18)
set(gcf, 'Color', 'w')
xlabel('$N$','interpreter','latex')
ylabel('$\langle e^{\mathrm{-}\Delta S_{med,sys}} \rangle_N$','interpreter','latex')
set(gca,'xScale','log') 
% set(gca,'yScale','log') 
xlim([min(Nvec) max(Nvec)])
% set(gca,'XTick',[10 1000 100000])
% ylim([-0.3 1.2])
title(['$\langle e^{\mathrm{-}\Delta S_{tot}} \rangle_{max(N)}=$',num2str(nanmean( exp( -DS )),'%1.2f')],'interpreter','latex')
legend('$\Delta S_{med}$','$\Delta S_{sys}$')
fig_setup


% figure
% scatter(Sm,Ds)
% axis square      
% set(gca, 'FontSize',18)
% set(gcf, 'Color', 'w')
% xlabel('$\Delta S_{med}$','interpreter','latex')
% ylabel('$\Delta S_{sys}$','interpreter','latex')
% box on
% fig_setup
% legend off
% 
% figure
% scatter(Sm,DS)
% axis square      
% set(gca, 'FontSize',18)
% set(gcf, 'Color', 'w')
% xlabel('$\Delta S_{med}$','interpreter','latex')
% ylabel('$\Delta S_{tot}$','interpreter','latex')
% box on
% fig_setup
% legend off
% 
% figure
% scatter(Ds,DS)
% axis square      
% set(gca, 'FontSize',18)
% set(gcf, 'Color', 'w')
% xlabel('$\Delta S_{sys}$','interpreter','latex')
% ylabel('$\Delta S_{tot}$','interpreter','latex')
% box on
% fig_setup
% legend off
% 
% 
% x=Sm;
% y=Ds;
% 
% x=Sm;
% y=DS;
% 
% x=Ds;
% y=DS;
% 
% 
% [c,lags] = xcorr(x,y);
% figure
% stem(lags,c)
% xlim([-100, 100])
% 
% % [xcf,lags,bounds] = crosscorr(x(y>0) ,y(y>0),'NumLags',100);
% figure
% crosscorr(x,y ,'NumLags',100)
% % vline(0,'k')
% axis square   
% set(gcf, 'Color', 'w')
% % lags(xcf==max(xcf))
% % xlim([-3*abs(lags(xcf==max(xcf))), 3*abs(lags(xcf==max(xcf)))])
% % xlim([-30, 30])
% 
% % if nargin >= 4
% %    save_path = varargin{1}; 
% %    savefig(h,fullfile(save_path,append(save_name,'_','entropy_plot.fig')),'compact')
% % end

if ischar(save_path)
    savefig(h,fullfile(save_path,append(save_name,'_','entropy_plot.fig')),'compact')
    for a = 1:length(h)
        exportgraphics(h(a),fullfile(save_path,append(save_name,'_',sprintf('entropy_plot_%d.png', a))))
    end
end

for i=1:5
    tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
    uiwait(gcf);
    delete(tmp_ui);
    close
end
end