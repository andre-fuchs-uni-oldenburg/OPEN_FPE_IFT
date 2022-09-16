%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The conditional PDF's will be plotted using different representations (shown here are only two see 
% Fig. \ref{fig:optimi_a}) to see the differences between optimized and non-optimized and 
% experimental conditional PDF's.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_pointwise_KM(y,x,y_mean_bin,x_mean_bin,x_bin_not_nan,x0,x1,P_B,P_AnB,P_AIB,counter_A,counter_B,min_events,eD1,eD2,markov,fval_start,fval,tau1,tau2,m_data,Fs,taylor_L,tol,dy,increment_bin,ub,lb,norm_ur,norm_r,save_path,save_name)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');


y_input     = y_mean_bin(counter_B>=min_events);
x_input     = x_mean_bin(counter_A>=min_events);

P_AIB       = P_AIB(counter_A>=min_events,counter_B>=min_events);
P_AnB       = P_AnB(counter_A>=min_events,counter_B>=min_events);
% P_A       = P_A(counter_A>=min_events);
P_B         = P_B(counter_B>=min_events);

D1_poly     = x0(1:sum(x_bin_not_nan));
D2_poly     = abs(x0(sum(x_bin_not_nan)+1:end));

D1_poly_optimiert     = x1(1:sum(x_bin_not_nan));
D2_poly_optimiert     = x1(sum(x_bin_not_nan)+1:2*sum(x_bin_not_nan));

%% non-optimized D1, D2
P_moment = nan;
P_moment = ShortTimeProp(y_input,x_input,tau1,tau2,D1_poly,D2_poly,taylor_L,m_data,Fs,norm_ur,norm_r);

%% optimized D1, D2
P_optimiert = nan;
P_optimiert = ShortTimeProp(y_input,x_input,tau1,tau2,D1_poly_optimiert,D2_poly_optimiert,taylor_L,m_data,Fs,norm_ur,norm_r);%an dieser Stelle auch x_mean_bin berücksichtigen!

%% Joint PDF
h(1) = figure;
surf(y(counter_B>=min_events),x(counter_A>=min_events),P_AnB)
hold on
% mesh(y_input,x_input,P_moment.*P_B)
mesh(y_input,x_input,P_optimiert.*P_B)
set(gca,'FontSize',18)
% set(gca,'zScale','log') 
axis([min(y_input) max(y_input) min(x_input) max(x_input)])
set(gcf, 'Color', 'w')
legend('$p_{exp}$','$p_{stp,opti}$','Location','Northeast')
% title(['directly=' num2str(fval_start,2) '; ' 'optimized=' num2str(fval,2)])
% xlabel(['$u_{r_2=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
% ylabel(['$u_{r_1=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
if norm_ur==1 
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
ylabel(['$u_{r''=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')
ylabel(['$u_{r''=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex') 
end
grid off
axis vis3d
fig_setup
rotate3d on
% close

%% conditional PDF
h(2) = figure;
surf(y(counter_B>=min_events),x(counter_A>=min_events),P_AIB)
hold on
% surf(y(counter_B>=min_events),x(counter_A>=min_events),real(P_AIB))
% mesh(y_input,x_input,real(P_moment))
%  set(gca,'zScale','log') 
%  axis([min(y_input) max(y_input) min(x_input) max(x_input) 0.001 10])
mesh(y_input,x_input,(P_optimiert))
set(gca,'FontSize',18)
axis([min(y_input) max(y_input) min(x_input) max(x_input)])
axis vis3d
set(gcf, 'Color', 'w')
legend('$p_{exp}$','$p_{stp,opti}$','Location','Northeast')
% title(['directly=' num2str(fval_start,2) '; ' 'optimized=' num2str(fval,2)])
% xlabel(['$u_{r_2=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
% ylabel(['$u_{r_1=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
if norm_ur==1 
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
ylabel(['$u_{r''=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')
ylabel(['$u_{r''=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')
end   
grid off
% close
fig_setup
rotate3d on


h(3) = figure;
surf(y(counter_B>=min_events),x(counter_A>=min_events),P_AIB)
hold on
% surf(y(counter_B>=min_events),x(counter_A>=min_events),real(P_AIB))
% mesh(y_input,x_input,real(P_moment))
 set(gca,'zScale','log') 
axis([min(y_input) max(y_input) min(x_input) max(x_input) 0.001 10])
mesh(y_input,x_input,(P_optimiert))
set(gca,'FontSize',18)
% axis([min(y_input) max(y_input) min(x_input) max(x_input)])
axis vis3d
set(gcf, 'Color', 'w')
legend('$p_{exp}$','$p_{stp,opti}$','Location','NorthWest')
% title(['directly=' num2str(fval_start,2) '; ' 'optimized=' num2str(fval,2)])
% xlabel(['$u_{r_2=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
% ylabel(['$u_{r_1=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
if norm_ur==1 
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
ylabel(['$u_{r''=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')
ylabel(['$u_{r''=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')  
end
grid off
% close
fig_setup
rotate3d on


%% difference: non-optimized and optimized joint PDF
h(4) = figure;
surf(y_input,x_input,P_optimiert.*P_B-P_AnB)
set(gca,'FontSize',18)
axis([min(y_input) max(y_input) min(x_input) max(x_input)])
axis vis3d
set(gcf, 'Color', 'w')
legend('$p_{stp,opti} - p_{exp}$')
% title(['directly=' num2str(fval_start,2) '; ' 'optimized=' num2str(fval,2)])
% xlabel(['$u_{r_2=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
% ylabel(['$u_{r_1=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
if norm_ur==1 
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
ylabel(['$u_{r''=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')
ylabel(['$u_{r''=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')  
end
grid off
colormap (parula(200))
colorbar
% close
fig_setup
rotate3d on


%% difference: non-optimized and optimized conditional PDF
h(5) = figure;
surf(y_input,x_input,P_optimiert-P_AIB)
set(gca,'FontSize',18)
axis([min(y_input) max(y_input) min(x_input) max(x_input)])
axis vis3d
set(gcf, 'Color', 'w')
legend('$p_{stp,opti} - p_{exp}$')
% title(['directly=' num2str(fval_start,2) '; ' 'optimized=' num2str(fval,2)])
% xlabel(['$u_{r_2=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
% ylabel(['$u_{r_1=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
if norm_ur==1 
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
ylabel(['$u_{r''=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')
ylabel(['$u_{r''=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')  
end
grid off
% caxis([min(z1(1,:)) max(z1(1,:))])
colormap (parula(200))
colorbar
% close
fig_setup
rotate3d on


%% D1
h(6) = figure;
errorbar(y_input,D1_poly,(eD1(x_bin_not_nan)),'MarkerSize',8,'Marker','o','LineStyle','none','LineWidth',1.5)
hold on
plot(y_input,D1_poly_optimiert,'MarkerSize',8,'Marker','o','LineStyle','none','LineWidth',1.5)
axis square
xlim([min(y_input) max(y_input)])
xPos = 0.0;
yPos = 0.0;
set(gca,'FontSize',18)
set(gcf, 'Color', 'w')
title(['$r=' num2str(tau2*m_data/(Fs*taylor_L),'%2.1f') ' ' '\lambda$'])
% xlabel(['$u_{r_2=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
% ylabel('$D^{(1)} (u_{r_2},r_2)$', 'interpreter','latex')
if norm_ur==1 
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')
end
ylabel('$D^{(1)} (u_{r},r)$', 'interpreter','latex')
% ylabel('$D1 (x,r) / [\frac{m}{s^2}]$', 'interpreter','latex')
fig_setup
plot([xPos xPos], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
plot(get(gca,'xlim'),[yPos yPos],'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
legend('non-optimized','optimized','Location','Northeast')
txt = {'(c)'};
c=text(4,0.5,txt);
c.Units='normalized';

%% difference: non-optimized and optimized D1
h(7) = figure;
plot(y_input,(D1_poly_optimiert-D1_poly),'MarkerSize',8,'Marker','o','LineStyle','none','LineWidth',1.5)
hold on
plot(y_input,(repmat(tol*(range(D1_poly)),sum(x_bin_not_nan))),'k','LineWidth',2,'LineStyle','--')
plot(y_input,-(repmat(tol*(range(D1_poly)),sum(x_bin_not_nan))),'k','LineWidth',2,'LineStyle','--')
% plot(y_input,ub(1:sum(x_bin_not_nan)),'LineWidth',2,'LineStyle','--','Color','k');
% plot(y_input,lb(1:sum(x_bin_not_nan)),'LineWidth',2,'LineStyle','--','Color','k');
set(gca,'FontSize',18)
xlim([min(y_input) max(y_input)])
axis square
set(gcf, 'Color', 'w')
% xlabel(['$u_{r_2=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
% ylabel('$D^{(1)} (u_{r_2},r_2)$', 'interpreter','latex')
if norm_ur==1 
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')
end
ylabel('$D^{(1)} (u_{r},r)$', 'interpreter','latex')
title('optimization difference')
fig_setup
plot([xPos xPos], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
plot(get(gca,'xlim'),[yPos yPos],'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
legend off


%% Potential of D1
h(8) = figure;
plot(y_input,-cumsum(D1_poly_optimiert.*dy),'MarkerSize',8,'Marker','o','LineStyle','none','LineWidth',1.5)
hold on
set(gca,'FontSize',18)
xlim([min(y_input) max(y_input)])
axis square
set(gcf, 'Color', 'w')
% xlabel(['$u_{r_2=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
ylabel('$-\int_{-\infty}^{u_{r_2}}D^{(1)} (u_{r_2},r_2)*du_{r_2}$', 'interpreter','latex')
if norm_ur==1 
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')   
end
ylabel('$-\int_{-\infty}^{u_{r}}D^{(1)} (u_{r},r)*du_{r}$', 'interpreter','latex')
xPos = 0.0;
yPos = 0.0;
% plot(get(gca,'xlim'),[yPos yPos],'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
fig_setup
plot([xPos xPos], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
legend off
title('potential')

%% D2
h(9) = figure;
errorbar(y_input,abs(D2_poly),(eD2(x_bin_not_nan)),'MarkerSize',8,'Marker','o','LineStyle','none','LineWidth',1.5)
hold on
plot(y_input,D2_poly_optimiert,'MarkerSize',8,'Marker','o','LineStyle','none','LineWidth',1.5)
hold on
a_grenze=y_input(find(D2_poly_optimiert==min(D2_poly_optimiert)));
b_grenze=y_input(find(y_input>=a_grenze,1,'first'));
plot(b_grenze-(y_input(y_input>=a_grenze)-b_grenze)...
    ,(D2_poly_optimiert(y_input>=a_grenze)),'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','LineWidth',1.5)
title(['$r=' num2str(tau2*m_data/(Fs*taylor_L),'%2.1f') ' ' '\lambda$'])
set(gca,'FontSize',18)
xlim([min(y_input) max(y_input)])
ylim([0 max([abs(D2_poly) D2_poly_optimiert])])
axis square
% xlabel(['$u_{r_2=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
% ylabel('$D^{(2)} (u_{r_2},r_2)$', 'interpreter','latex')
if norm_ur==1 
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')   
end
ylabel('$D^{(2)} (u_{r},r)$', 'interpreter','latex')
% ylabel('$D2 (x,r) / [\frac{m^2}{s^3}]$', 'interpreter','latex')
set(gcf, 'Color', 'w')
xPos = 0.0;
yPos = 0.0;
fig_setup
plot([xPos xPos], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
% plot(get(gca,'xlim'),[yPos yPos],'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
legend('non-optimized','optimized','Location','Northeast')
txt = {'(d)'};
d=text(4,0.5,txt);
d.Units='normalized';


%% difference: non-optimized and optimized D2
h(10) = figure;
plot(y_input,(D2_poly_optimiert-(D2_poly)),'MarkerSize',8,'Marker','o','LineStyle','none','LineWidth',1.5)
hold on
plot(y_input,(repmat(tol*(range(D2_poly)),sum(x_bin_not_nan))),'k','LineWidth',2,'LineStyle','--')
plot(y_input,-(repmat(tol*(range(D2_poly))+max(D2_poly),sum(x_bin_not_nan))),'k','LineWidth',2,'LineStyle','--')
% plot(y_input,ub(x0(sum(x_bin_not_nan)+1:end)),'LineWidth',2,'LineStyle','--','Color','k');
% plot(y_input,lb(x0(sum(x_bin_not_nan)+1:end)),'LineWidth',2,'LineStyle','--','Color','k');
set(gca,'FontSize',18)
xlim([min(y_input) max(y_input)])
axis square
set(gcf, 'Color', 'w')
% xlabel(['$u_r (r_2=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}) / \sigma_\infty$'], 'interpreter','latex')
% ylabel('$D^{(2)} (u_{r_2},r_2)$', 'interpreter','latex')
if norm_ur==1 
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')   
end
ylabel('$D^{(2)} (u_{r},r)$', 'interpreter','latex')
title('optimization difference')
fig_setup
plot([xPos xPos], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
plot(get(gca,'xlim'),[yPos yPos],'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
legend off




%% contour plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [C,h]=contour(y,x,(P_AIB),level,'-','color','k','LineWidth',2,'ShowText','on');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cut_a=round(3*increment_bin/10);
% cut_0=round(increment_bin/2);
% cut_b=round(7*increment_bin/10);

% a_tmp=max(max(P_AIB));
% level = [0.005*a_tmp 0.05*a_tmp 0.1*a_tmp 0.20*a_tmp 0.45*a_tmp 0.65*a_tmp 0.85*a_tmp 0.95*a_tmp];
level=round(logspace(log10(0.01*max(max((P_AIB)))*10^6),log10(0.9*max(max((P_AIB)))*10^6),7))./10^6;

h(11) = figure;
contour(y(counter_B>=min_events),x(counter_A>=min_events),real(P_AIB),level,'-','color','k','LineWidth',2)
hold on
contour(y_input,x_input,real(P_moment),level,'-','color',[0    0.4470    0.7410],'LineWidth',2)
contour(y_input,x_input,real(P_optimiert),level,'-','color','r','LineWidth',2)
set(gca,'FontSize',18)
% axis([-5 5 -5 5])
axis square
set(gcf, 'Color', 'w')
% xPos = 0.0;
% yPos = 0.0;
% plot([xPos xPos], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
% plot(get(gca,'xlim'),[yPos yPos],'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
% plot(get(gca,'ylim'),get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','r');
% hold off
% xlabel(['$u_{r_2=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
% ylabel(['$u_{r_1=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
if norm_ur==1 
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
ylabel(['$u_{r''=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')
ylabel(['$u_{r''=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')
end   
fig_setup
title(['$\xi_{stp}$=' num2str(fval_start,2) '; ' '$\xi_{stp,opti}$=' num2str(fval,2)])
xPos = 0.0;
yPos = 0.0;
plot([xPos xPos], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
plot(get(gca,'xlim'),[yPos yPos],'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
plot(get(gca,'ylim'),get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k');
% legend('experiment','non-optimized','optimized','Location','NorthWest')
legend('$p_{exp}$','$p_{stp}$','$p_{stp,opti}$','Location','NorthWest')
txt = {'(a)'};
a=text(4,0.5,txt);
a.Units='normalized';

%vertical cut
% plot([y(cut_0) y(cut_0)], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
% plot([y(cut_a) y(cut_a)], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color',[0    0.4470    0.7410]); % Adapts to y limits of current axes
% plot([y(cut_b) y(cut_b)], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','r'); % Adapts to y limits of current axes
% horizontal cut
% plot(get(gca,'xlim'),[x(cut_0) x(cut_0)],'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
% % plot(get(gca,'ylim'),get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','r');
% plot(get(gca,'xlim'),[x(cut_a) x(cut_a)],'LineWidth',2,'LineStyle','--','Color',[0    0.4470    0.7410]); % Adapts to y limits of current axes
% plot(get(gca,'xlim'),[x(cut_b) x(cut_b)],'LineWidth',2,'LineStyle','--','Color','r'); % Adapts to y limits of current axes

% xlabel(['$\xi (r_2=' num2str(tau2*m_data/Fs) 'm' ') / \sigma_\infty$'], 'interpreter','latex')
% ylabel(['$\xi (r_1=' num2str(tau1*m_data/Fs) 'm' ') / \sigma_\infty$'], 'interpreter','latex')
% xlabel(['$\xi (r_2=' num2str(tau2*m_data/(Fs*taylor_L),'%2.1f') 'r/\lambda' ') / \sigma_\infty$'], 'interpreter','latex')
% ylabel(['$\xi (r_1=' num2str(tau1*m_data/(Fs*taylor_L),'%2.1f') 'r/\lambda' ') / \sigma_\infty$'], 'interpreter','latex')



h(12) = figure;
contour(y(counter_B>=min_events),x(counter_A>=min_events),real(P_AIB),level,'-','color','k','LineWidth',2)
hold on
contour(y_input,x_input,real(P_optimiert),level,'-','color','r','LineWidth',2)
set(gca,'FontSize',18)
% axis([-5 5 -5 5])
axis square
set(gcf, 'Color', 'w')
title(['$\xi_{stp}$=' num2str(fval_start,2) '; ' '$\xi_{stp,opti}$=' num2str(fval,2)])
xPos = 0.0;
yPos = 0.0;
%vertical cut
% plot([y(cut_0) y(cut_0)], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
% plot([y(cut_a) y(cut_a)], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color',[0    0.4470    0.7410]); % Adapts to y limits of current axes
% plot([y(cut_b) y(cut_b)], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','r'); % Adapts to y limits of current axes
% horizontal cut
% plot(get(gca,'xlim'),[x(cut_0) x(cut_0)],'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
% % plot(get(gca,'ylim'),get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','r');
% plot(get(gca,'xlim'),[x(cut_a) x(cut_a)],'LineWidth',2,'LineStyle','--','Color',[0    0.4470    0.7410]); % Adapts to y limits of current axes
% plot(get(gca,'xlim'),[x(cut_b) x(cut_b)],'LineWidth',2,'LineStyle','--','Color','r'); % Adapts to y limits of current axes

% xlabel(['$\xi (r_2=' num2str(tau2*m_data/Fs) 'm' ') / \sigma_\infty$'], 'interpreter','latex')
% ylabel(['$\xi (r_1=' num2str(tau1*m_data/Fs) 'm' ') / \sigma_\infty$'], 'interpreter','latex')
% xlabel(['$\xi (r_2=' num2str(tau2*m_data/(Fs*taylor_L),'%2.1f') 'r/\lambda' ') / \sigma_\infty$'], 'interpreter','latex')
% ylabel(['$\xi (r_1=' num2str(tau1*m_data/(Fs*taylor_L),'%2.1f') 'r/\lambda' ') / \sigma_\infty$'], 'interpreter','latex')
fig_setup
plot([xPos xPos], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
plot(get(gca,'xlim'),[yPos yPos],'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
plot(get(gca,'ylim'),get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k');
% legend('experiment','optimized','Location','NorthWest')
legend('$p_{exp}$','$p_{stp,opti}$','Location','NorthWest')
% xlabel(['$u_{r_2=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
% ylabel(['$u_{r_1=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
if norm_ur==1 
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
ylabel(['$u_{r''=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')
ylabel(['$u_{r''=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')
end




h(13) = figure;
contour3(y(counter_B>=min_events),x(counter_A>=min_events),real(P_AIB),level,'-','color','k','LineWidth',2)
hold on
contour3(y_input,x_input,real(P_optimiert),level,'-','color','r','LineWidth',2)
set(gca,'FontSize',18)
set(gca,'zScale','log')
axis([min(y_input) max(y_input) min(x_input) max(x_input) min(level) 2*max(level)])
% axis([-5 5 -5 5])
axis square
set(gcf, 'Color', 'w')
% legend('experiment','optimized','Location','NorthWest')
legend('$p_{exp}$','$p_{stp,opti}$','Location','NorthWest','Position',[0.281785714285714 0.657642855417161 0.165338666098458 0.113571428117298])
title(['$\xi_{stp}$=' num2str(fval_start,2) '; ' '$\xi_{stp,opti}$=' num2str(fval,2)])
% title(['directly=' num2str(fval_start,2) '; ' 'optimized=' num2str(fval,2)])
grid off
% xlabel(['$u_{r_2=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
% ylabel(['$u_{r_1=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
if norm_ur==1 
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
ylabel(['$u_{r''=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
xlabel(['$u_{r=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')
ylabel(['$u_{r''=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}$'], 'interpreter','latex')
end
zlabel('PDF','interpreter','latex')
% xPos = 0.0;
% yPos = 0.0;
% plot([xPos xPos], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
% plot(get(gca,'xlim'),[yPos yPos],'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
% plot(get(gca,'ylim'),get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k');
%vertikaler Schnitt
% plot([y(cut_0) y(cut_0)], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
% plot([y(cut_a) y(cut_a)], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color',[0    0.4470    0.7410]); % Adapts to y limits of current axes
% plot([y(cut_b) y(cut_b)], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','r'); % Adapts to y limits of current axes
% %horizontaler Schnitt
% plot(get(gca,'xlim'),[x(cut_0) x(cut_0)],'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
% % plot(get(gca,'ylim'),get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','r');
% plot(get(gca,'xlim'),[x(cut_a) x(cut_a)],'LineWidth',2,'LineStyle','--','Color',[0    0.4470    0.7410]); % Adapts to y limits of current axes
% plot(get(gca,'xlim'),[x(cut_b) x(cut_b)],'LineWidth',2,'LineStyle','--','Color','r'); % Adapts to y limits of current axes

% xlabel(['$\xi (r_2=' num2str(tau2*m_data/Fs) 'm' ') / \sigma_\infty$'], 'interpreter','latex')
% ylabel(['$\xi (r_1=' num2str(tau1*m_data/Fs) 'm' ') / \sigma_\infty$'], 'interpreter','latex')
% xlabel(['$\xi (r_2=' num2str(tau2*m_data/(Fs*taylor_L),'%2.1f') 'r/\lambda' ') / \sigma_\infty$'], 'interpreter','latex')
% ylabel(['$\xi (r_1=' num2str(tau1*m_data/(Fs*taylor_L),'%2.1f') 'r/\lambda' ') / \sigma_\infty$'], 'interpreter','latex')
fig_setup
rotate3d on

txt = {'(b)'};
b=text(4,0.5,txt);
b.Units='normalized';

font=22;
a.FontSize = font;
b.FontSize = font;
c.FontSize = font;
d.FontSize = font;

pos_txt=[-0.22   0.9];
a.Position=pos_txt;
b.Position=[-0.3   0.9];
c.Position=pos_txt;
d.Position=pos_txt;

% gaus_bin = -5:(2*5/(increment_bin-1)):5; 
% % norm = normpdf(gaus_bin,0,1);
% norm1 = normpdf(gaus_bin,0,std(incr1)); 
% norm2 = normpdf(gaus_bin,0,std(incr2));   
% 
% figure; 
% plot(gaus_bin,norm1,'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176])
% hold on
%     plot(x_mean_bin,(hist(incr1,x_mean_bin)./nansum(hist(incr1,x_mean_bin))./dx),'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',3,'Marker','o', 'LineStyle','none');
%     plot(gaus_bin,norm2.*10,'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176])
%     plot(y_mean_bin,(hist(incr1,y_mean_bin)./nansum(hist(incr2,y_mean_bin))./dy).*10,'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],'MarkerSize',3,'Marker','o', 'LineStyle','none');
%    
% set(gca, 'FontSize',22,'PlotBoxAspectRatio',[1 1 1])
% set(gca,'yScale','log') 
% set(gcf, 'Color', 'w')
% xlabel('$\xi / \sigma_\infty$', 'interpreter','latex')
% ylabel('pdf($\xi,r$)(a.u.)', 'interpreter','latex')
% axis([-5 5 10^(-5) 0.5*10^(2)])


% optimized FIT of D1, D2
% co=co_KM_opti;
% r_2  = ((tau2./Fs).*m_data)/taylor_L;
% 
% D1_poly_optimiert_fit =    (co.a(1).*(r_2(k)).^co.ea(1)+co.a(2)).*y_input;
% D2_poly_optimiert_fit =   ((co.b(1).*(r_2(k)).^co.eb(1)+co.b(2)).*y_input.^0)+...
%             ((co.b(3).*(r_2(k)).^co.eb(2)+co.b(4)).*y_input.^1)+...
%             ((co.b(5).*(r_2(k)).^co.eb(3)+co.b(6)).*y_input.^2);
% 
% P_optimiert_fit = ShortTimeProp(y_input,x_input,tau1,tau2,D1_poly_optimiert_fit,D2_poly_optimiert_fit,taylor_L,m_data,Fs);
% 
% % plot(y_input,D1_poly_optimiert_fit,'-','color','k','LineWidth',2)
% % plot(y_input,D2_poly_optimiert_fit,'-','color','k','LineWidth',2)
% 
% % level=round(logspace(log10(0.01*max(max((P_AIB)))*10^6),log10(0.9*max(max((P_AIB)))*10^6),7))./10^6;
% 
% h(14) = figure;
% contour(y(counter_B>=min_events),x(counter_A>=min_events),real(P_AIB),level,'-','color','k','LineWidth',2)
% hold on
% % contour(y_input,x_input,real(P_moment),level,'-','color',[0    0.4470    0.7410],'LineWidth',2)
% contour(y_input,x_input,real(P_optimiert_fit),level,'-','color',[0    0.4470    0.7410],'LineWidth',2)
% contour(y_input,x_input,real(P_optimiert),level,'-','color','r','LineWidth',2)
% set(gca,'FontSize',18)
% % axis([-5 5 -5 5])
% axis square
% set(gcf, 'Color', 'w')
% % xPos = 0.0;
% % yPos = 0.0;
% % plot([xPos xPos], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
% % plot(get(gca,'xlim'),[yPos yPos],'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
% % plot(get(gca,'ylim'),get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','r');
% % hold off
% xlabel(['$u_{r_2=' num2str(tau2/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
% ylabel(['$u_{r_1=' num2str(tau1/markov,'%1.1f') '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
% fig_setup
% title(['directly=' num2str(fval_start,2) '; ' 'optimized=' num2str(fval,2)])
% xPos = 0.0;
% yPos = 0.0;
% plot([xPos xPos], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
% plot(get(gca,'xlim'),[yPos yPos],'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
% plot(get(gca,'ylim'),get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k');
% % legend('experiment','non-optimized','optimized','Location','NorthWest')



if ischar(save_path)
    savefig(h,fullfile(save_path,append(save_name,'_','KM_STP_optimization.fig')),'compact')
    for a = 1:length(h)
        exportgraphics(h(a),fullfile(save_path,append(save_name,'_',sprintf('KM_STP_optimization_%d.png', a))))
    end
end

for i=1:13
    tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
    uiwait(gcf);
    delete(tmp_ui);
    close
end
end