%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots the $k$-th order structure function $S^{k}(r)$ with $k={1-7}$ for scales 
% $\lambda \leq r \leq L$. Dashed line represents -4/5 law. The 4/5 law states that 3rd-order 
% structure functions of longitudinal velocity increments should scale linearly with their 
% displacement distances. The 3rd-order structure function is one of the most fundamental of all 
% Navier-Stokes equation results; it is exact, with no adjustable constants, so any model that is 
% expected to produce accurate turbulence results must reasonably well duplicate this. 
%
% Arguments IN
% r_n = Number of seperated scales between Integral and Taylor (calculation of structure functions)
% Fs = Acquisition/Sampling Frequency in Hz
% data = it is the filtered data 
% m_data = mean of the data
% int_L = Integral length scale in meters
% taylor_L = Taylors length scale in meters
%
% Arguments OUT
% r = scales between Integral and Taylor (calculation of structure functions)
% S_exp = $k$-th order structure function $S^{k}(r)$ with $k={1-7}$
% Lam_sq_L_1 = shape parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [r,S_exp_2,S_exp_3,S_exp_4,S_exp_5,S_exp_6,S_exp_7,Lam_sq_L_1,Lam_sq_L_2,Lam_sq_L_3]=plot_struc_function(r_n,Fs,data,m_data,int_L,taylor_L,save_path)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');

clear r tau step
r               = unique(round(logspace(log10(taylor_L*10^6),log10(int_L*10^6),r_n))./10^6,'stable');
% r                 = unique(round(logspace(log10(taylor_L*10^6),log10(int_L*10^6),r_n))./10^6,'stable');
% r                 = logspace(-3,0,r_n);
% r                 = unique(round(logspace(log10(evaluated(1).r*10^6),log10(evaluated(end).r*10^6),r_n))./10^6,'stable');

tau(1:length(r))    = nan;
step(1:length(r))   = nan;

% m_data = mean(data);

for k=1:length(r)
    tau(k)  = r(k)/m_data;
    step(k) = round(tau(k)*Fs);  
end

stepp = unique(step);

if stepp(1) ==0
    stepp  = stepp(2:end);
end

clear r r_2 S_exp_2 S_exp_3 S_exp_4 S_exp_5 S_exp_6

r(1:length(stepp))                     = nan;
S_exp_2(1:length(stepp))               = nan; 
S_exp_3(1:length(stepp))               = nan;
S_exp_4(1:length(stepp))               = nan; 
S_exp_5(1:length(stepp))               = nan; 
S_exp_6(1:length(stepp))               = nan;  
S_exp_7(1:length(stepp))               = nan;  

Lam_sq_L_1(1:length(stepp))            = nan;  
Lam_sq_L_2(1:length(stepp))            = nan;  
Lam_sq_L_3(1:length(stepp))            = nan;  


tmp = length(data);

parfor k=1:length(stepp)
    k/length(stepp)   
    tau1            = stepp(k); 
    r(k)            = (tau1./Fs).*m_data;
    
    incr1   =   (data(tau1+1:tmp)     - data(1:tmp-tau1)).';
    
%% EXP    
    Lam_sq_L_1(k)  = log(kurtosis(incr1)/3)/4;
    Lam_sq_L_2(k)  = log(kurtosis(incr1))/4;

    S_exp_2(k)  = mean(incr1.^2);
    S_exp_3(k)  = mean(incr1.^3);
    S_exp_4(k)  = mean(incr1.^4);
    S_exp_5(k)  = mean(incr1.^5);
    S_exp_6(k)  = mean(incr1.^6);
    S_exp_7(k)  = mean(incr1.^7);
    
    Lam_sq_L_3(k)   = log(S_exp_4(k)/(3*S_exp_2(k)*S_exp_2(k)));
end

%% Plot structure functions
in_n_nan   = ~isnan(r);
r          = r(in_n_nan);
r          = r./taylor_L;
S_exp_2    = S_exp_2(in_n_nan); 
S_exp_3    = S_exp_3(in_n_nan);
S_exp_4    = S_exp_4(in_n_nan); 
S_exp_5    = S_exp_5(in_n_nan);
S_exp_6    = S_exp_6(in_n_nan);
S_exp_7    = S_exp_7(in_n_nan);

fontsize=20;
% line=2;

h(1) = figure('WindowState', 'maximized');
set(gcf, 'Color', 'w')
subplot(3,2,1)
plot(r,S_exp_2,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
set(gca,'xScale','log') 
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
xlabel('$r/m$', 'interpreter','latex')
xlabel('$r / \lambda$', 'interpreter','latex')
ylabel('$S^2(r)$', 'interpreter','latex')
xlim([min(r)*0.5 max(r)*2])
ylim([min(S_exp_2) max(S_exp_2)])
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
axis square


subplot(3,2,2)
plot(r,-S_exp_3,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
hold on
[xData, yData] = prepareCurveData( r,-S_exp_3 );
% Set up fittype and options.
ft = fittype( '(4/5)*a*x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
fitresult.a=fitresult.a*1.6;
plot(r,feval(fitresult,r),'LineWidth',2,'LineStyle','--','Color','r')
set(gca,'xScale','log')
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
% xlabel('$r$', 'interpreter','latex')
xlabel('$r/m$', 'interpreter','latex')
xlabel('$r / \lambda$', 'interpreter','latex')
ylabel('$-S^3(r)$', 'interpreter','latex')
xlim([min(r)*0.5 max(r)*2])
ylim([min(-S_exp_3) max(-S_exp_3)])
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
axis square


subplot(3,2,3)
plot(r,S_exp_4,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
set(gca,'xScale','log')
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
% xlabel('$r$', 'interpreter','latex')
xlabel('$r/m$', 'interpreter','latex')
xlabel('$r / \lambda$', 'interpreter','latex')
ylabel('$S^4(r)$', 'interpreter','latex')
xlim([min(r)*0.5 max(r)*2])
ylim([min(S_exp_4) max(S_exp_4)])
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
axis square

subplot(3,2,4)
plot(r,-S_exp_5,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
set(gca,'xScale','log')
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
% xlabel('$r$', 'interpreter','latex')
xlabel('$r/m$', 'interpreter','latex')
xlabel('$r / \lambda$', 'interpreter','latex')
ylabel('$-S^5(r)$', 'interpreter','latex')
xlim([min(r)*0.5 max(r)*2])
ylim([min(-S_exp_5) max(-S_exp_5)])
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
axis square

subplot(3,2,5)
plot(r,S_exp_6,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
set(gca,'xScale','log')
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
% xlabel('$r$', 'interpreter','latex')
xlabel('$r/m$', 'interpreter','latex')
xlabel('$r / \lambda$', 'interpreter','latex')
ylabel('$S^6(r)$', 'interpreter','latex')
xlim([min(r)*0.5 max(r)*2])
ylim([min(S_exp_6) max(S_exp_6)])
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
axis square 

subplot(3,2,6)
plot(r,-S_exp_7,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
set(gca,'xScale','log')
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
xlabel('$r / \lambda$', 'interpreter','latex')
xlabel('$r / \lambda$', 'interpreter','latex')
ylabel('$-S^7(r)$', 'interpreter','latex')
xlim([min(r)*0.5 max(r)*2])
ylim([min(-S_exp_7) max(-S_exp_7)])
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
axis square
% set(h(1),'Position',[1 46 835 1299]);
% fig_setup
% legend off



%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
h(2) = figure;
subplot(1,2,1)
p1=plot(r,S_exp_2,'MarkerSize',8,'Marker','+','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
hold on
p2=plot(r,S_exp_4,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
p3=plot(r,S_exp_6,'MarkerSize',8,'Marker','*','LineWidth',2,'LineStyle','none','Color',[0 0 0]); 
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
set(gca,'xScale','log')
set(gca,'yScale','log') 
axis square
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
xlabel('$r/m$', 'interpreter','latex')
xlabel('$r / \lambda$', 'interpreter','latex')
xlim([min(r)*0.5 max(r)*2])
ylim([0.5*min([S_exp_2 S_exp_4 S_exp_6]) 2*max([S_exp_2 S_exp_4 S_exp_6])])
legend([p1,p2,p3],{'$S^2(r)$','$S^4(r)$','$S^6(r)$'},'Interpreter','latex','Location','northwest','FontSize',18)  

subplot(1,2,2)
p4=plot(r,-S_exp_3,'MarkerSize',8,'Marker','+','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
hold on
p5=plot(r,-S_exp_5,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
p6=plot(r,-S_exp_7,'MarkerSize',8,'Marker','*','LineWidth',2,'LineStyle','none','Color',[0 0 0]); 
plot(r,feval(fitresult,r),'LineWidth',2,'LineStyle','--','Color','r')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
set(gca,'xScale','log')
set(gca,'yScale','log') 
axis square
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
xlabel('$r/m$', 'interpreter','latex')
xlabel('$r / \lambda$', 'interpreter','latex')
xlim([min(r)*0.5 max(r)*2])
ylim([0.5*min([-S_exp_3 -S_exp_5 -S_exp_7]) 2*max([-S_exp_3 -S_exp_5 -S_exp_7])])
legend([p4,p5,p6],{'$-S^3(r)$','$-S^5(r)$','$-S^7(r)$'},'Interpreter','latex','Location','northwest','FontSize',18)  
% vline(taylor_L,'k')
% vline(int_L,'k')
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')

h(3) = figure;
subplot(1,3,1)
plot(r,Lam_sq_L_1,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none');
hold on
plot(r,Lam_sq_L_2,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none');
plot(r,Lam_sq_L_3,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
% set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
% set(gca,'xScale','log')
% set(gca,'yScale','log') 
axis square
set(gcf, 'Color', 'w')
set(gca,'FontSize',20)
xlabel('$r/m$', 'interpreter','latex')
xlabel('$r / \lambda$', 'interpreter','latex')
% xlim([min(r)*0.5 max(r)*2])

legend({'$\frac{ln\left(K(u_r)/3\right)}{4}$',...
    '$\frac{ln\left(K(u_r)\right)}{4}$',...
    '$ln\left(\frac{S^4(r)}{3\ \left(S^2(r)\right)^2}\right)$'},'Interpreter','latex','Location','northeast','FontSize',18)  

subplot(1,3,2)
plot(r,Lam_sq_L_1,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none');
hold on
plot(r,Lam_sq_L_2,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none');
plot(r,Lam_sq_L_3,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
set(gca,'xScale','log')
% set(gca,'yScale','log') 
axis square
set(gcf, 'Color', 'w')
set(gca,'FontSize',20)
xlabel('$r/m$', 'interpreter','latex')
xlabel('$r / \lambda$', 'interpreter','latex')
xlim([min(r)*0.5 max(r)*2])

legend({'$\frac{ln\left(K(u_r)/3\right)}{4}$',...
    '$\frac{ln\left(K(u_r)\right)}{4}$',...
    '$ln\left(\frac{S^4(r)}{3\ \left(S^2(r)\right)^2}\right)$'},'Interpreter','latex','Location','northeast','FontSize',18)  

subplot(1,3,3)
plot(r,Lam_sq_L_1,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none');
hold on
plot(r,Lam_sq_L_2,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none');
plot(r,Lam_sq_L_3,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
set(gca,'xScale','log')
set(gca,'yScale','log') 
axis square
set(gcf, 'Color', 'w')
set(gca,'FontSize',20)
xlabel('$r/m$', 'interpreter','latex')
xlabel('$r / \lambda$', 'interpreter','latex')
xlim([min(r)*0.5 max(r)*2])

legend({'$\frac{ln\left(K(u_r)/3\right)}{4}$',...
    '$\frac{ln\left(K(u_r)\right)}{4}$',...
    '$ln\left(\frac{S^4(r)}{3\ \left(S^2(r)\right)^2}\right)$'},'Interpreter','latex','Location','northeast','FontSize',18)  
fig_setup

if ischar(save_path)
   savefig(h,fullfile(save_path,'structure_function.fig'),'compact')
   for a = 1:length(h)
        exportgraphics(h(a),fullfile(save_path,sprintf('structure_function_%d.png', a)))
   end
end

for i=1:3
    tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
    uiwait(gcf);
    delete(tmp_ui);
    close
end
end