function [r,S_exp_2,S_exp_3,S_exp_6,S_KM_ST_opti_func_2,S_KM_ST_opti_func_3,S_KM_ST_opti_func_6]=...
    Struc_recon(r_n,increment_bin,co,m_data,Fs,data,evaluated,int_L,taylor_L,step_r,norm_ur,norm_r)
% data=data_filter(1:3*10^6);
% % co=co_IFT_opti;
% co=co_KM_opti;
% step_r=1;
% r_n=60;

% We begin with the 3rd-order structure function because, as we noted in Chap. 2, this is one of the
% most fundamental of all Navier?Stokes equation results; it is exact, with no adjustable constants, so
% any model that is expected to produce accurate turbulence results must reasonably well duplicate
% this. This is the so-called 4/5 law which states that 3rd-order structure functions of longitudinal
% velocity increments should scale linearly with their displacement distances.

%beim direkt kram muss noch dire Grenzen dx und so weiter angepasst werden
%-siehe unten bei marginal pdf habe ich schon mal angefangen aber nicht zu
%ende gebracht ----- es musss insbesondere auf das r geachtet werden  -
%bei der optiemierung ist heir erst bei 2 lambda lsogegangen....


%% short time propagator
clear r tau step
r               = unique(round(logspace(log10(0.8*taylor_L*10^6),log10(1.2*int_L*10^6),r_n))./10^6,'stable');
% r                 = unique(round(logspace(log10(taylor_L*10^6),log10(int_L*10^6),r_n))./10^6,'stable');
% r                 = logspace(-3,0,r_n);
% r                 = unique(round(logspace(log10(evaluated(1).r*10^6),log10(evaluated(end).r*10^6),r_n))./10^6,'stable');

tau(1:length(r))    = nan;
step(1:length(r))   = nan;

for k=1:length(r)
    tau(k)  = r(k)/m_data;
    step(k) = round(tau(k)*Fs);  
end

stepp = unique(step);
stepp = stepp+step_r;

if stepp(1) ==0
     stepp  = stepp(2:end);
end
stepp=stepp(stepp>step_r);

clear r r_2 S_exp_2 S_exp_3 S_exp_4 S_exp_5 S_exp_6 S_KM_ST_opti_func_2 S_KM_ST_opti_func_3 S_KM_ST_opti_func_4 S_KM_ST_opti_func_5 S_KM_ST_opti_func_6

r(1:length(stepp))                     = nan;
r_2(1:length(stepp))                   = nan; 
S_exp_2(1:length(stepp))               = nan; 
S_exp_3(1:length(stepp))               = nan;
S_exp_4(1:length(stepp))               = nan; 
S_exp_5(1:length(stepp))               = nan; 
S_exp_6(1:length(stepp))               = nan;  
S_exp_7(1:length(stepp))               = nan;  
S_KM_ST_opti_func_2(1:length(stepp))   = nan;
S_KM_ST_opti_func_3(1:length(stepp))   = nan; 
S_KM_ST_opti_func_4(1:length(stepp))   = nan;
S_KM_ST_opti_func_5(1:length(stepp))   = nan; 
S_KM_ST_opti_func_6(1:length(stepp))   = nan; 
S_KM_ST_opti_func_7(1:length(stepp))   = nan; 
x_A_STP(1:length(stepp),increment_bin) = nan;
P_A_STP(1:length(stepp),increment_bin) = nan;

%  S_pdf(1:length(stepp))                 = nan; 
%  S_P_AIB1(1:length(stepp))              = nan; 
%  S_P_AIB2(1:length(stepp))              = nan; 
%  S_KM_ST_opti_2(1:length(stepp))        = nan; 
%  S_KM_ST_opti_3(1:length(stepp))        = nan; 
%  S_KM_ST_opti_6(1:length(stepp))        = nan;   
%  S_KM_ST_direkt(1:length(stepp))        = nan; 
%  mean_bin(1:increment_bin)                       = nan;

parfor k=1:length(stepp)
    k/length(stepp)   
    tau1            = stepp(k)-step_r;
    tau2            = stepp(k);
    [incr1,incr2]   = Increment(tau1,tau2,data);
    [dx,dy,x,y,dA]  = limit(incr1,incr2,increment_bin);  
    r(k)            = (tau1./Fs).*m_data;
%     ((tau1./Fs).*m_data)./taylor_L
        
%% Min EVENTS 
%     min_events          = 0;
%     [fuL uL]            = hist(incr1,1000);
%     min_index           = fuL>=min_events;
%     incr1_tmp           = incr1>=(uL(find(min_index==1,1,'first'))) &  incr1 <=(uL(find(min_index==1,1,'last')));
%     %  incr1_tmp=( incr1 <=(uL(find(min_index==1,1,'last'))));
%     incr1               = incr1(incr1_tmp);
%     incr2               = incr2(incr1_tmp);

    [P_AIB,P_BIA,P_AnB,P_A,P_B,binP_AIB,x_mean_bin,y_mean_bin,events,counter_A,counter_B] = distribution(tau1,tau2,incr1,incr2,increment_bin,dx,dy,x,y,dA,1);
    
%% EXP    
    S_exp_2(k)  = mean(incr1.^2);
    S_exp_3(k)  = mean(incr1.^3);
    S_exp_4(k)  = mean(incr1.^4);
    S_exp_5(k)  = mean(incr1.^5);
    S_exp_6(k)  = mean(incr1.^6);
    S_exp_7(k)  = mean(incr1.^7);
%     S_exp_2(k)  = trapz(x,(x.^2.*P_A));
%     S_exp_3(k)  = trapz(x,(x.^3.*P_A));
%     S_exp_4(k)  = trapz(x,(x.^4.*P_A));
%     S_exp_5(k)  = trapz(x,(x.^5.*P_A));
%     S_exp_6(k)  = trapz(x,(x.^6.*P_A)); 
%     S_exp_7(k)  = trapz(x,(x.^7.*P_A)); 

%% short time propagator    
%   r_2     = r(k);    
    r_2(k)  = ((tau2./Fs).*m_data)/taylor_L;
            
    D1_poly =    (co.a(1).*(r_2(k)).^co.ea(1)+co.a(2)).*y;
    D2_poly =   ((co.b(1).*(r_2(k)).^co.eb(1)+co.b(2)).*y.^0)+...
                ((co.b(3).*(r_2(k)).^co.eb(2)+co.b(4)).*y.^1)+...
                ((co.b(5).*(r_2(k)).^co.eb(3)+co.b(6)).*y.^2);
            
    P_n_opti_func           = ShortTimeProp(y,x,tau1,tau2,D1_poly,D2_poly,taylor_L,m_data,Fs,norm_ur,norm_r);      
    
    x_A_STP(k,:)            = x;
    P_A_STP(k,:)            = (P_n_opti_func*P_B.'.*dy).';
    
%     if k==1
%         P_A_STP(k,:)            = (P_n_opti_func*P_B.'.*dy).';
%     else
%         
%     end

    figure
    plot(y,P_B.','LineWidth',2)
    hold on
    plot(x,P_A.','LineWidth',2)
%     P_B_STP(k,:)            = (P_BIA*P_A.'.*dx).';
%     plot(x,P_B_STP(k,:),'k','LineWidth',2)
    plot(x,P_A_STP(k,:),'k','LineWidth',2)
    set(gca,'yScale','log') 
    set(gcf, 'Color', 'w')
    set(gca,'FontSize',22)
    xlabel('$u_{r} / \sigma_\infty$', 'interpreter','latex')
    ylabel('$p(u_{r} / \sigma_\infty)$','interpreter','latex')
    axis square

    %an dieser Stelle vielleicht besser k=1 mit exp und dann später nur noch S_KM_ST_opti_func
    %nutzen um P_B zu berechenen 
    
    S_KM_ST_opti_func_2(k)  = trapz(x,(x.^2.*P_A_STP(k,:)));
    S_KM_ST_opti_func_3(k)  = trapz(x,(x.^3.*P_A_STP(k,:)));
    S_KM_ST_opti_func_4(k)  = trapz(x,(x.^4.*P_A_STP(k,:)));
    S_KM_ST_opti_func_5(k)  = trapz(x,(x.^5.*P_A_STP(k,:)));
    S_KM_ST_opti_func_6(k)  = trapz(x,(x.^6.*P_A_STP(k,:)));  
    S_KM_ST_opti_func_7(k)  = trapz(x,(x.^7.*P_A_STP(k,:)));  
    
%     S_KM_ST_opti_func_2(k)  = trapz(x,(x.^2.*(P_n_opti_func*P_B.'.*dy).'));
%     S_KM_ST_opti_func_3(k)  = trapz(x,(x.^3.*(P_n_opti_func*P_B.'.*dy).'));
%     S_KM_ST_opti_func_4(k)  = trapz(x,(x.^4.*(P_n_opti_func*P_B.'.*dy).'));
%     S_KM_ST_opti_func_5(k)  = trapz(x,(x.^5.*(P_n_opti_func*P_B.'.*dy).'));
%     S_KM_ST_opti_func_6(k)  = trapz(x,(x.^6.*(P_n_opti_func*P_B.'.*dy).'));  
%     S_KM_ST_opti_func_7(k)  = trapz(x,(x.^7.*(P_n_opti_func*P_B.'.*dy).'));  
%     S_KM_ST_opti_func_2(k)  = nansum(x.^2.*(P_n_opti_func*P_B.'.*dy).'.*dx);
%     S_KM_ST_opti_func_3(k)  = nansum(x.^3.*(P_n_opti_func*P_B.'.*dy).'.*dx);  
%     S_KM_ST_opti_func_6(k)  = nansum(x.^6.*(P_n_opti_func*P_B.'.*dy).'.*dx);

% close all
% Struc_recon_test_plot(P_AIB,P_BIA,P_AnB,P_A,P_B,binP_AIB,x_mean_bin,y_mean_bin,events,counter_A,counter_B,min_events,dx,dy,x,y,P_n_opti_func,tau1,tau2,incr1,incr2,step_r,increment_bin,m_data,Fs,taylor_L)
% -mean(incr1.^3)
% -trapz(x,(x.^3.*(P_AIB*P_B.'.*dy).'))
% -trapz(x,(x.^3.*(P_n_opti_func*P_B.'.*dy).'))
% 
% 
% figure
% plot(x,cumtrapz(x,(x.^3.*P_A)))
% hold on
% plot(x,cumtrapz(x,(x.^3.*(P_n_opti_func*P_B.'.*dy).')))
% % 
% % 
% % 
% figure
% plot(x,P_A,'LineWidth',2)
% hold on
% plot(y,P_B,'LineWidth',2)
% plot(x,(P_n_opti_func*P_B.'.*dy).','k','LineWidth',2)
% set(gca,'yScale','log') 
% set(gcf, 'Color', 'w')
% set(gca,'FontSize',22)
% xlabel('$u_{r} / \sigma_\infty$', 'interpreter','latex')
% ylabel('$p(u_{r} / \sigma_\infty)$','interpreter','latex')
% axis square
% legend({'$p(u_{r_1})$','$p(u_{r_2})$','$p_{stp}(u_{r_1})$'},'Interpreter','latex','Location','south','FontSize',18)
% xlim([min(y) max(y)])
% ylim([0.5*min(P_A) 2*max(P_A)])

end


    [incr1,incr2] = Increment(stepp(end)-step_r,stepp(1)-step_r,data);
    
    [fuL, uL]       = hist(incr1.',increment_bin); 
    puL             = fuL / ( nansum(fuL) * mean(diff(uL)) ); % relative frequency to approximate p(u,r=L)
    [ful, ul]       = hist(incr2.',uL); 
    pul             = ful / ( nansum(ful) * mean(diff(ul)) ); % relative frequency to approximate p(u,r=lambda)

figure
plot(uL,puL,'LineWidth',2)
hold on
plot(ul,pul,'LineWidth',2)
plot(x_A_STP(1,:),P_A_STP(1,:),'k','LineWidth',2)
plot(x_A_STP(end,:),P_A_STP(end,:),'k','LineWidth',2)
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca,'FontSize',22)
xlabel('$u_{r} / \sigma_\infty$', 'interpreter','latex')
ylabel('$p(u_{r} / \sigma_\infty)$','interpreter','latex')
axis square
% legend({'$p(u_{r_1})$','$p(u_{r_2})$','$p_{stp}(u_{r_1})$'},'Interpreter','latex','Location','south','FontSize',18)
xlim([min(uL) max(uL)])
ylim([10^-7 2*max(pul)])




%% short time propagator plot
%% Plot structure functions
in_n_nan   = ~isnan(r);
r          = r(in_n_nan);
S_exp_2    = S_exp_2(in_n_nan); 
S_exp_3    = S_exp_3(in_n_nan);
S_exp_4    = S_exp_4(in_n_nan); 
S_exp_5    = S_exp_5(in_n_nan);
S_exp_6    = S_exp_6(in_n_nan);
S_exp_7    = S_exp_7(in_n_nan);

%  S_pdf=S_pdf(in_n_nan); 
%  S_P_AIB1=S_P_AIB1(in_n_nan);
%  S_P_AIB2=S_P_AIB2(in_n_nan);
%  S_KM_ST_opti_2=S_KM_ST_opti_2(in_n_nan); 
%  S_KM_ST_opti_3=S_KM_ST_opti_3(in_n_nan);
%  S_KM_ST_opti_6=S_KM_ST_opti_6(in_n_nan);
%  S_KM_ST_direkt=S_KM_ST_direkt(in_n_nan);

S_KM_ST_opti_func_2    = S_KM_ST_opti_func_2(in_n_nan);
S_KM_ST_opti_func_3    = S_KM_ST_opti_func_3(in_n_nan);
S_KM_ST_opti_func_4    = S_KM_ST_opti_func_4(in_n_nan);
S_KM_ST_opti_func_5    = S_KM_ST_opti_func_5(in_n_nan);
S_KM_ST_opti_func_6    = S_KM_ST_opti_func_6(in_n_nan);
S_KM_ST_opti_func_7    = S_KM_ST_opti_func_7(in_n_nan);

fontsize=22;
line=2;

h(1) = figure;
set(gcf, 'Color', 'w')
subplot(3,2,1)
plot(r./taylor_L,S_exp_2,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
hold on
plot(r./taylor_L,S_KM_ST_opti_func_2,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 
set(gca,'xScale','log') 
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
xlabel('$r/\lambda$', 'interpreter','latex')
% xlabel('$r / \lambda$', 'interpreter','latex')
ylabel('$S^2(r)$', 'interpreter','latex')
xlim([min(r./taylor_L)*0.5 max(r./taylor_L)*2])
ylim([min(S_exp_2) max(S_exp_2)])
vline(evaluated(1).r/taylor_L,'r')
vline(taylor_L/taylor_L,'k')
vline(int_L/taylor_L,'k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0]);
axis square


subplot(3,2,2)
plot(r./taylor_L,-S_exp_3,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
hold on
plot(r./taylor_L,-S_KM_ST_opti_func_3,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 
set(gca,'xScale','log')
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
% xlabel('$r$', 'interpreter','latex')
xlabel('$r / \lambda$', 'interpreter','latex')
ylabel('$-S^3(r)$', 'interpreter','latex')
xlim([min(r./taylor_L)*0.5 max(r./taylor_L)*2])
ylim([min(-S_exp_3) max(-S_exp_3)])
vline(evaluated(1).r/taylor_L,'r')
vline(taylor_L/taylor_L,'k')
vline(int_L/taylor_L,'k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0]);
axis square


subplot(3,2,3)
plot(r./taylor_L,S_exp_4,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
hold on
plot(r./taylor_L,S_KM_ST_opti_func_4,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
set(gca,'xScale','log')
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
% xlabel('$r$', 'interpreter','latex')
xlabel('$r / \lambda$', 'interpreter','latex')
ylabel('$S^4(r)$', 'interpreter','latex')
xlim([min(r./taylor_L)*0.5 max(r./taylor_L)*2])
ylim([min(S_exp_4) max(S_exp_4)])
vline(evaluated(1).r/taylor_L,'r')
vline(taylor_L/taylor_L,'k')
vline(int_L/taylor_L,'k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0]);
axis square

subplot(3,2,4)
plot(r./taylor_L,-S_exp_5,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
hold on
plot(r./taylor_L,-S_KM_ST_opti_func_5,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 
set(gca,'xScale','log')
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
% xlabel('$r$', 'interpreter','latex')
xlabel('$r / \lambda$', 'interpreter','latex')
ylabel('$-S^5(r)$', 'interpreter','latex')
xlim([min(r./taylor_L)*0.5 max(r./taylor_L)*2])
ylim([min(-S_exp_5) max(-S_exp_5)])
vline(evaluated(1).r/taylor_L,'r')
vline(taylor_L/taylor_L,'k')
vline(int_L/taylor_L,'k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0]);
axis square

subplot(3,2,5)
plot(r./taylor_L,S_exp_6,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
hold on
plot(r./taylor_L,S_KM_ST_opti_func_6,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
set(gca,'xScale','log')
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
% xlabel('$r$', 'interpreter','latex')
xlabel('$r / \lambda$', 'interpreter','latex')
ylabel('$S^6(r)$', 'interpreter','latex')
xlim([min(r./taylor_L)*0.5 max(r./taylor_L)*2])
ylim([min(S_exp_6) max(S_exp_6)])
vline(evaluated(1).r/taylor_L,'r')
vline(taylor_L/taylor_L,'k')
vline(int_L/taylor_L,'k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0]);
axis square 

subplot(3,2,6)
plot(r./taylor_L,-S_exp_7,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
hold on
plot(r./taylor_L,-S_KM_ST_opti_func_7,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 
set(gca,'xScale','log')
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
% xlabel('$r$', 'interpreter','latex')
xlabel('$r / \lambda$', 'interpreter','latex')
ylabel('$-S^7(r)$', 'interpreter','latex')
xlim([min(r./taylor_L)*0.5 max(r./taylor_L)*2])
ylim([min(-S_exp_7) max(-S_exp_7)])
vline(evaluated(1).r/taylor_L,'r')
vline(taylor_L/taylor_L,'k')
vline(int_L/taylor_L,'k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0]);
axis square
set(h(1),'Position',[1 46 835 1299]);


% h(1) = figure;
% set(gcf, 'Color', 'w')
% subplot(1,5,1)
% plot(r,S_exp_2,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
% hold on
% plot(r,S_KM_ST_opti_func_2,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 
% set(gca,'xScale','log') 
% set(gca,'yScale','log') 
% set(gcf, 'Color', 'w')
% set(gca,'FontSize',fontsize)
% xlabel('$r / \lambda$', 'interpreter','latex')
% ylabel('$S^2(r)$', 'interpreter','latex')
% xlim([min(r)*0.5 max(r)*2])
% ylim([min(S_exp_2) max(S_exp_2)])
% vline(evaluated(1).r/taylor_L,'r')
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
% axis square
% 
% subplot(1,5,2)
% plot(r,-S_exp_3,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
% hold on
% plot(r,-S_KM_ST_opti_func_3,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 
% set(gca,'xScale','log')
% set(gca,'yScale','log') 
% set(gcf, 'Color', 'w')
% set(gca,'FontSize',fontsize)
% xlabel('$r / \lambda$', 'interpreter','latex')
% ylabel('$-S^3(r)$', 'interpreter','latex')
% xlim([min(r)*0.5 max(r)*2])
% ylim([min(-S_exp_3) max(-S_exp_3)])
% vline(evaluated(1).r/taylor_L,'r')
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
% axis square
% 
% subplot(1,5,3)
% plot(r,S_exp_4,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
% hold on
% plot(r,S_KM_ST_opti_func_4,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
% set(gca,'xScale','log')
% set(gca,'yScale','log') 
% set(gcf, 'Color', 'w')
% set(gca,'FontSize',fontsize)
% xlabel('$r / \lambda$', 'interpreter','latex')
% ylabel('$S^4(r)$', 'interpreter','latex')
% xlim([min(r)*0.5 max(r)*2])
% ylim([min(S_exp_4) max(S_exp_4)])
% vline(evaluated(1).r/taylor_L,'r')
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
% axis square
% 
% subplot(1,5,4)
% plot(r,-S_exp_5,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
% hold on
% plot(r,-S_KM_ST_opti_func_5,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 
% set(gca,'xScale','log')
% set(gca,'yScale','log') 
% set(gcf, 'Color', 'w')
% set(gca,'FontSize',fontsize)
% xlabel('$r / \lambda$', 'interpreter','latex')
% ylabel('$-S^5(r)$', 'interpreter','latex')
% xlim([min(r)*0.5 max(r)*2])
% ylim([min(-S_exp_5) max(-S_exp_5)])
% vline(evaluated(1).r/taylor_L,'r')
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
% axis square
% 
% subplot(1,5,5)
% plot(r,S_exp_6,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
% hold on
% plot(r,S_KM_ST_opti_func_6,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
% set(gca,'xScale','log')
% set(gca,'yScale','log') 
% set(gcf, 'Color', 'w')
% set(gca,'FontSize',fontsize)
% xlabel('$r / \lambda$', 'interpreter','latex')
% ylabel('$S^6(r)$', 'interpreter','latex')
% xlim([min(r)*0.5 max(r)*2])
% ylim([min(S_exp_6) max(S_exp_6)])
% vline(evaluated(1).r/taylor_L,'r')
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
% axis square 
% set(h(1),'Position',[1 923 2560 415]);




%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
h(1) = figure;
p11=plot(r,S_exp_2,'MarkerSize',8,'Marker','+','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
hold on
p1=plot(r,-S_exp_3,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
% p2=loglog(r,4/5*(r)*0.0175,'k','LineWidth',line,'LineStyle',':'); 
p3=plot(r,S_exp_6,'MarkerSize',8,'Marker','*','LineWidth',2,'LineStyle','none','Color',[0 0 0]); 
p33=plot(r,S_KM_ST_opti_func_2,'LineWidth',line,'Color','k');
p4=plot(r,-S_KM_ST_opti_func_3,'LineWidth',line,'Color','k');
p5=plot(r,S_KM_ST_opti_func_6,'LineWidth',line,'Color','k');
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0]);
set(gca,'xScale','log')
set(gca,'yScale','log') 
axis square
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
xlabel('$r$', 'interpreter','latex')
% xlabel('$r / \lambda$', 'interpreter','latex')
% ylabel('$S^k(r)$', 'interpreter','latex')
xlim([min(r)*0.5 max(r)*2])
ylim([0.5*min([S_exp_2 -S_exp_3 S_exp_6]) 2*max([S_exp_2 -S_exp_3 S_exp_6])])
legend([p11,p1,p3,p4],{'exp $S^2(r)$','exp $-S^3(r)$','exp $S^6(r)$','$D^{(1,2)}$ reconst'},'Interpreter','latex','Location','northwest','FontSize',18)  
% legend([p11,p1,p3,p2,p4],{'exp $S^2(r)$','exp $-S^3(r)$','exp $S^6(r)$','Kol 4/5-law','$D^{(1,2)}$ reconst'},'Interpreter','latex','Location','northwest','FontSize',18)  
% vline(evaluated(1).r,'r')
vline(taylor_L,'k')
vline(int_L,'k')



% h(2) = figure;
% p11=plot(r,S_exp_2,'MarkerSize',8,'Marker','+','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
% hold on
% p1=plot(r,-S_exp_3,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
% p2=loglog(r,S_exp_4,'MarkerSize',8,'Marker','x','LineWidth',2,'LineStyle','none','Color',[0 0 0]); 
% p22=plot(r,-S_exp_5,'MarkerSize',8,'Marker','v','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
% p3=plot(r,S_exp_6,'MarkerSize',8,'Marker','*','LineWidth',2,'LineStyle','none','Color',[0 0 0]); 
% 
% p33=plot(r,S_KM_ST_opti_func_2,'LineWidth',line,'Color','k');
% p4=plot(r,-S_KM_ST_opti_func_3,'LineWidth',line,'Color','k');
% p44=plot(r,S_KM_ST_opti_func_4,'LineWidth',line,'Color','k');
% p55=plot(r,-S_KM_ST_opti_func_5,'LineWidth',line,'Color','k');
% p5=plot(r,S_KM_ST_opti_func_6,'LineWidth',line,'Color','k');
% set(gca,'xScale','log')
% set(gca,'yScale','log') 
% axis square
% set(gcf, 'Color', 'w')
% set(gca,'FontSize',fontsize)
% xlabel('$r / \lambda$', 'interpreter','latex')
% % ylabel('$S^k(r)$', 'interpreter','latex')
% xlim([min(r/taylor_L) max(r/taylor_L)]);
% ylim([0.5*min([S_exp_2 -S_exp_3 S_exp_6]) 2*max([S_exp_2 -S_exp_3 S_exp_6])])
% 
% legend([p11,p1,p2,p22,p3,p4],{'exp $S^2(r)$','exp $-S^3(r)$','exp $S^4(r)$','exp $-S^5(r)$','exp $S^6(r)$','$D^{(1,2)}$ reconst'},'Interpreter','latex','Location','northwest','FontSize',18)  
% 
% vline(evaluated(1).r/taylor_L,'r')
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')



                   

% figure
% plot(r,S_exp_2./r,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
% hold on
% plot(r,(-2.*(d11+d22).*S_exp_2+d20),'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 

%  figure
%  plot(r,-S_exp_3,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
%  hold on
% %  plot(r,-S_KM_ST_opti_func_3./r,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 
%  plot(r,-(-3.*(d11+2.*d22).*S_exp_3-(6.d21.*S_exp_2)),'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 
 
% figure
% plot(r,S_exp_6./r,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
% hold on
% plot(r,...
%     (-6.*(d11+5.*d22).*S_exp_6)-(30.d21.*S_exp_5)),'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 





%% chapmann Kolmogorov 










%% Diff d_ij for very small steps
clear S_exp_1 S_exp_2 S_exp_3 S_exp_4 S_exp_5 S_exp_6 S_exp_7 S_exp_8
% r_diff      = round(r(1)/m_data*Fs):1:round(r(end)/m_data*Fs);
r_diff      = round(r(end)/m_data*Fs):-1:round(r(1)/m_data*Fs);
% r_diff      = round(evaluated(end).r/m_data*Fs):-1:round(evaluated(1).r/m_data*Fs);
% Renner
% r_diff      = round((1.2*int_L)/m_data*Fs):-1:round((0.8*taylor_L)/m_data*Fs);

parfor k=1:length(r_diff)
    k/length(r_diff)
    [incr1,~]   = Increment(r_diff(k),r_diff(k),data);
    
    S_exp_2(k)  = mean(incr1.^2);
    S_exp_3(k)  = mean(incr1.^3);
    S_exp_4(k)  = mean(incr1.^4);
    S_exp_5(k)  = mean(incr1.^5);
    S_exp_6(k)  = mean(incr1.^6);
    S_exp_7(k)  = mean(incr1.^7);
    S_exp_8(k)  = mean(incr1.^8);
end

r_diff_exp      = r_diff*m_data/Fs;
r_diff_exp      = r_diff_exp./taylor_L;
dr_diff_exp     = mean(diff(r_diff_exp));

clear r_diff
r_diff          = linspace(r_diff_exp(1),r_diff_exp(end),10^3);
% r_diff          = unique(round(logspace(log10(0.8*taylor_L*10^6),log10(1.2*int_L*10^6),r_n))./10^6,'stable');
% r_diff          = unique(round(logspace(log10(evaluated(end).r/taylor_L*10^6),log10(evaluated(1).r/taylor_L*10^6),r_n))./10^6,'stable');
dr_diff         = mean(diff(r_diff));

% r_diff          = flip(r_diff);
% r_diff_exp      = flip(r_diff_exp);
% S_exp_1         = flip(S_exp_1);
% S_exp_2         = flip(S_exp_2);
% S_exp_3         = flip(S_exp_3);
% S_exp_4         = flip(S_exp_4);
% S_exp_5         = flip(S_exp_5);
% S_exp_6         = flip(S_exp_6);
% dr_diff_exp     = mean(diff(r_diff_exp));
% dr_diff         = mean(diff(r_diff));

% d20             = (co.b(1).*(r_diff_exp).^co.eb(1)+co.b(2));        
% d21_R           =  -0.009.*(r_diff_exp(1:end-1)).^(0.2);
% d22_R           =   0.043.*(r_diff_exp(1:end-1)).^(-0.73);


%% Solve systems of linear equations AF = B
clear F
F(1:4,1:length(r_diff_exp)-1) = nan; 
for k=2:length(r_diff_exp)
    A = [2*S_exp_2(1,k-1)      2                                 0                      2*S_exp_2(1,k-1)   ;              
         3*S_exp_3(1,k-1)      0                           6*S_exp_2(1,k-1)             6*S_exp_3(1,k-1)   ;
         4*S_exp_4(1,k-1)     12*S_exp_2(1,k-1)           12*S_exp_3(1,k-1)            12*S_exp_4(1,k-1)   ;
         5*S_exp_5(1,k-1)     20*S_exp_3(1,k-1)           20*S_exp_4(1,k-1)            20*S_exp_5(1,k-1)   ;
         6*S_exp_6(1,k-1)     30*S_exp_4(1,k-1)           30*S_exp_5(1,k-1)            30*S_exp_6(1,k-1)   ;
         7*S_exp_7(1,k-1)     42*S_exp_5(1,k-1)           42*S_exp_6(1,k-1)            42*S_exp_7(1,k-1)   ;
         8*S_exp_8(1,k-1)     56*S_exp_6(1,k-1)           56*S_exp_7(1,k-1)            56*S_exp_8(1,k-1)]  ;

%% mit d20 aus dem Fit von D2
%     d20_R = 0.033.*(r_diff_exp(k-1)).^(0.25);
%     A = [2*S_exp_2(1,k-1)      d20_R                          0              2*S_exp_2(1,k-1)   ;              
%          3*S_exp_3(1,k-1)      d20_R            6*S_exp_2(1,k-1)             6*S_exp_3(1,k-1)   ;
%          4*S_exp_4(1,k-1)      d20_R           12*S_exp_3(1,k-1)            12*S_exp_4(1,k-1)   ;
%          5*S_exp_5(1,k-1)      d20_R           20*S_exp_4(1,k-1)            20*S_exp_5(1,k-1)   ;
%          6*S_exp_6(1,k-1)      d20_R           30*S_exp_5(1,k-1)            30*S_exp_6(1,k-1)]  ;

    b = [(S_exp_2(1,k)-S_exp_2(1,k-1))./dr_diff_exp;...
         (S_exp_3(1,k)-S_exp_3(1,k-1))./dr_diff_exp;...
         (S_exp_4(1,k)-S_exp_4(1,k-1))./dr_diff_exp;...
         (S_exp_5(1,k)-S_exp_5(1,k-1))./dr_diff_exp;...
         (S_exp_6(1,k)-S_exp_6(1,k-1))./dr_diff_exp;... 
         (S_exp_7(1,k)-S_exp_7(1,k-1))./dr_diff_exp;...
         (S_exp_8(1,k)-S_exp_8(1,k-1))./dr_diff_exp];
 
F(:,k-1)    = A \ b;
clear A b

% tmp = linsolve(A,b); 
%% orthonormal basis
% tmp=null(A);
% d = det(A(1:4,:));
%% If the matrix is full rank, then the rank is equal to the number of columns, size(A,2)
% rank(A) 
% size(A,2)
%% If you do not get back the original vector b, than Ax = b does not have an exact solution. In this case, pinv(A)*b returns a least-squares solution
% A*pinv(A)*b-b
%% You can determine whether Ax =b has an exact solution by finding the row reduced echelon form of the augmented matrix [A b].
% rref([A b])
%% Since the bottom row contains all zeros except for the last entry, the equation does not have a solution. In this case, pinv(A) returns a least-squares solution.
end

%%%%%%%%%%%%% Reconstruction of S^k
clear S_direct_1_b S_direct_3_b S_direct_4_b S_direct_5_b S_direct_6_b
for k=1:length(r_diff_exp)
    if k==1
        S_direct_2_b(1,k)    = S_exp_2(1);
        S_direct_3_b(1,k)    = S_exp_3(1);    
        S_direct_4_b(1,k)    = S_exp_4(1); 
        S_direct_5_b(1,k)    = S_exp_5(1);  
        S_direct_6_b(1,k)    = S_exp_6(1);  
    else
        d11 = F(1,k-1);
        d20 = F(2,k-1);
        d21 = F(3,k-1);
        d22 = F(4,k-1); 
       
       S_direct_2_b(1,k)   = S_direct_2_b(1,k-1) +... 
                             ((2*S_direct_2_b(1,k-1)*(d11+d22) + 2*d20)*dr_diff_exp);
      
       S_direct_3_b(1,k)   = S_direct_3_b(1,k-1) +... 
                             ((3*S_direct_3_b(1,k-1)*(d11+2*d22) + 6*S_direct_2_b(1,k-1)*d21)*dr_diff_exp);
       
       S_direct_4_b(1,k)   = S_direct_4_b(1,k-1) +... 
                             ((4*S_direct_4_b(1,k-1)*(d11+3*d22) + 12*S_direct_3_b(1,k-1)*d21 + 12*S_direct_2_b(1,k-1)*d20)*dr_diff_exp);
       
       S_direct_5_b(1,k)   = S_direct_5_b(1,k-1) +... 
                             ((5*S_direct_5_b(1,k-1)*(d11+4*d22) + 20*S_direct_4_b(1,k-1)*d21 + 20*S_direct_3_b(1,k-1)*d20)*dr_diff_exp);
       
       S_direct_6_b(1,k)   = S_direct_6_b(1,k-1) +... 
                             ((6*S_direct_6_b(1,k-1)*(d11+5*d22) + 30*S_direct_5_b(1,k-1)*d21 + 30*S_direct_4_b(1,k-1)*d20)*dr_diff_exp);       
    end
end

%% mit meinen optimierten D1,D2
% co    = co_IFT_opti;
% co    = co_renner;
% co    = co_KM_opti;
% co    = co_KM_opti_no_offset;

% K62
% co.a=[-0.36 0];
% co.ea=-1;
% % % co.b=[0 0 0 0 0 0];
% % % co.eb=[0 0 0];
% co.b=[0 0 0 0 0.0144 0];
% co.eb=[0 0 -1];

% idy     = linspace(0,2,10);
idy     = 0:0.2:2;
idy     = idy(2:end);
idy     = 1;

clear S_direct_2 S_direct_3 S_direct_4 S_direct_5 S_direct_6 a b c
S_direct_2(1:length(idy),1:length(r_diff))               = nan; 
S_direct_3(1:length(idy),1:length(r_diff))               = nan; 
S_direct_4(1:length(idy),1:length(r_diff))               = nan;
S_direct_5(1:length(idy),1:length(r_diff))               = nan; 
S_direct_6(1:length(idy),1:length(r_diff))               = nan;
% S_direct_2b(1:length(r_diff))                           = nan; 
% S_direct_3b(1:length(r_diff))                           = nan; 
% S_direct_6b(1:length(r_diff))                           = nan; 
% a(1:length(r_diff))                                     = nan; 
% b(1:length(r_diff))                                     = nan; 
% c(1:length(r_diff))                                     = nan; 

for idx=1:length(idy)
    for k=1:length(r_diff)
        if k==1
            i_tmp=find(abs(r_diff_exp-r_diff(1))==min(abs(r_diff_exp-r_diff(1))));
    %         S_direct_1(idx,k)    = S_exp_1(i_tmp)*idy(idx);
            S_direct_2(idx,k)    = S_exp_2(i_tmp)*idy(idx);
            S_direct_3(idx,k)    = S_exp_3(i_tmp)*idy(idx);    
            S_direct_4(idx,k)    = S_exp_4(i_tmp)*idy(idx); 
            S_direct_5(idx,k)    = S_exp_5(i_tmp)*idy(idx);  
            S_direct_6(idx,k)    = S_exp_6(i_tmp)*idy(idx);  
    %         S_direct_6b(k)   = S_exp_6(i_tmp);
    %         S_direct_2b(k)   = S_exp_2(i_tmp);
    %         S_direct_3b(k)   = S_exp_3(i_tmp);
        else

            d11 = (co.a(1).*(r_diff(k-1)).^co.ea(1)+co.a(2));
            d20 = (co.b(1).*(r_diff(k-1)).^co.eb(1)+co.b(2));
            d21 = (co.b(3).*(r_diff(k-1)).^co.eb(2)+co.b(4));
            d22 = (co.b(5).*(r_diff(k-1)).^co.eb(3)+co.b(6));

    %         d11 = -d11;
    %         d20 = d20;
    %         d21 = -d21;
    %         d22 = -d22;

    %         d11 = (co.a(1).*(r_diff(k)).^co.ea(1)+co.a(2));
    %         d20 = (co.b(1).*(r_diff(k)).^co.eb(1)+co.b(2));
    %         d21 = (co.b(3).*(r_diff(k)).^co.eb(2)+co.b(4));
    %         d22 = (co.b(5).*(r_diff(k)).^co.eb(3)+co.b(6)); 

    %Renner
    %         d10 = -0.0015.*(r_diff(k-1)).^0.6;
    %         d11 =   -0.61.*(r_diff(k-1)).^(-0.67);
    %         d12 =  0.0096.*(r_diff(k-1)).^(0);
    %         d13 = -0.0023.*(r_diff(k-1)).^(0.3);      
    %         d20 =   0.033.*(r_diff(k-1)).^(0.25);
    %         d21 =  -0.009.*(r_diff(k-1)).^(0.2);
    %         d22 =   0.043.*(r_diff(k-1)).^(-0.73);

           S_direct_2(idx,k)   = S_direct_2(idx,k-1) -... 
                                 ((2*S_direct_2(idx,k-1)*(d11+d22) + 2*d20)*dr_diff);

           S_direct_3(idx,k)   = S_direct_3(idx,k-1) -... 
                                 ((3*S_direct_3(idx,k-1)*(d11+2*d22) + 6*S_direct_2(idx,k-1)*d21)*dr_diff);

           S_direct_4(idx,k)   = S_direct_4(idx,k-1) -... 
                                 ((4*S_direct_4(idx,k-1)*(d11+3*d22) + 12*S_direct_3(idx,k-1)*d21 + 12*S_direct_2(idx,k-1)*d20)*dr_diff);

           S_direct_5(idx,k)   = S_direct_5(idx,k-1) -... 
                                 ((5*S_direct_5(idx,k-1)*(d11+4*d22) + 20*S_direct_4(idx,k-1)*d21 + 20*S_direct_3(idx,k-1)*d20)*dr_diff);

           S_direct_6(idx,k)   = S_direct_6(idx,k-1) -... 
                                 ((6*S_direct_6(idx,k-1)*(d11+5*d22) + 30*S_direct_5(idx,k-1)*d21 + 30*S_direct_4(idx,k-1)*d20)*dr_diff);


    %% LGS Lösung
    %     clear A b 
    %     A = [2*S_direct_2(idx,k-1)      2                           0                           2*S_direct_2(idx,k-1)   ;              
    %          3*S_direct_3(idx,k-1)      0                           6*S_direct_2(idx,k-1)       6*S_direct_3(idx,k-1)   ;
    %          4*S_direct_4(idx,k-1)     12*S_direct_2(idx,k-1)      12*S_direct_3(idx,k-1)      12*S_direct_4(idx,k-1)   ;
    %          5*S_direct_5(idx,k-1)     20*S_direct_3(idx,k-1)      20*S_direct_4(idx,k-1)      20*S_direct_5(idx,k-1)]   ;

    %     b = [(S_direct_2(idx,k)-S_direct_2(idx,k-1))./-mean(diff(r_diff));...
    %         (S_direct_3(idx,k)-S_direct_3(idx,k-1))./-mean(diff(r_diff));...
    %         (S_direct_4(idx,k)-S_direct_4(idx,k-1))./-mean(diff(r_diff));...
    %         (S_direct_5(idx,k)-S_direct_5(idx,k-1))./-mean(diff(r_diff))] ;

    %     F(:,k-1) = A \ b;       

    %     a(k)   =  -2*d11*S_direct_2(k-1)*dr_diff;
    %     b(k)   =  -2*d22*S_direct_2(k-1)*dr_diff;
    %     c(k)   =  -2*d20*dr_diff;

    %     a(k)   =  -3*d11*S_direct_3(k-1)*dr_diff;
    %     b(k)   =  -3*2*d22*S_direct_3(k-1)*dr_diff;
    %     c(k)   =  -6*d21*S_direct_2(k-1)*dr_diff; 

    % exp Stützstellen
    %     S_direct_2b(k)   = S_exp_2(k-1) +... 
    %                        ((-2*(d11+d22)*S_exp_2(k-1) - 2*d20)*dr_diff);

    %     S_direct_3b(k)   = S_exp_3(k-1) +... 
    %                        ((-3*(d11+2*d22)*S_exp_3(k-1) - 6*d21*S_exp_2(k-1))*dr_diff);

    %     S_direct_6b(k)   = S_exp_6(k-1) +... 
    %                        ((-6*(d11+5*d22)*S_exp_6(k-1) - 30*(d21*S_exp_5(k-1)+d20*S_exp_4(k-1)))*dr_diff); 

    % TEST
    %     S_direct_2(k)   = S_direct_2(k-1) -... 
    %                       ((2*S_direct_2(k-1)*(d11+d22) + 2*d20)*dr_diff);

    %     S_direct_3(k)   = S_direct_3(k-1) -... 
    %                       ((3*S_direct_3(k-1)*(d11+2*d22) + 6*S_direct_2(k-1)*d21)*dr_diff);

    %     S_direct_2(k)   = S_direct_2(k-1) -... 
    %                       ((2*S_direct_2(k-1)*(d11+d22)*r_diff(k-1) + 2*d20*r_diff(k-1))*dr_diff);

    %     S_direct_3(k)   = S_direct_3(k-1) -... 
    %                       ((3*S_direct_3(k-1)*(d11+2*d22)*r_diff(k-1) + 6*S_direct_2(k-1)*d21*r_diff(k-1))*dr_diff);
        end
    end
end

% close all
Struc_recon_plot


 



        


% -----------------------------------------------------------------------------------------------------   
%% tmp Diff d_ij for very small steps
% mit S1
% A = [1*S_exp_1(1,k-1)      0                           0                            0   ; 
%      2*S_exp_2(1,k-1)      2                           2*S_exp_1(1,k-1)             2*S_exp_2(1,k-1)   ;              
%      3*S_exp_3(1,k-1)      6*S_exp_1(1,k-1)            6*S_exp_2(1,k-1)             6*S_exp_3(1,k-1)   ;
%      4*S_exp_4(1,k-1)     12*S_exp_2(1,k-1)           12*S_exp_3(1,k-1)            12*S_exp_4(1,k-1)   ;
%      5*S_exp_5(1,k-1)     20*S_exp_3(1,k-1)           20*S_exp_4(1,k-1)            20*S_exp_5(1,k-1)   ;
%      6*S_exp_6(1,k-1)     30*S_exp_4(1,k-1)           30*S_exp_5(1,k-1)            30*S_exp_6(1,k-1)]   ;
%  
% b = [(S_exp_1(1,k)-S_exp_1(1,k-1))./dr_diff_exp;...
%      (S_exp_2(1,k)-S_exp_2(1,k-1))./dr_diff_exp;...
%      (S_exp_3(1,k)-S_exp_3(1,k-1))./dr_diff_exp;...
%      (S_exp_4(1,k)-S_exp_4(1,k-1))./dr_diff_exp;...
%      (S_exp_5(1,k)-S_exp_5(1,k-1))./dr_diff_exp;...
%      (S_exp_6(1,k)-S_exp_6(1,k-1))./dr_diff_exp] ;

%d20=0 und d21=0 
% A = [2*S_exp_2(1,k-1)      0                            0                             2*S_exp_2(1,k-1)   ;              
%      3*S_exp_3(1,k-1)      0                            0                             6*S_exp_3(1,k-1)   ;
%      4*S_exp_4(1,k-1)      0                            0                            12*S_exp_4(1,k-1)   ;
%      5*S_exp_5(1,k-1)      0                            0                            20*S_exp_5(1,k-1)   ;
%      6*S_exp_6(1,k-1)      0                            0                            30*S_exp_6(1,k-1)]   ;
% 
% F(:,k-1) = A \ b; 
% clear A b
% end
% 
% F(2,:)=0;
% F(3,:)=0;


%        S_direct_1_b(1,k)   = S_direct_1_b(1,k-1) +... 
%            (S_direct_1_b(1,k-1)*d11*dr_diff_exp); 
%        
%        S_direct_2_b(1,k)   = S_direct_2_b(1,k-1) +... 
%            ((2*S_direct_2_b(1,k-1)*(d11+d22) + 2*S_direct_1_b(1,k-1)*d21 + 2*d20)*dr_diff_exp);
%       
%        S_direct_3_b(1,k)   = S_direct_3_b(1,k-1) +... 
%            ((3*S_direct_3_b(1,k-1)*(d11+2*d22) + 6*S_direct_2_b(1,k-1)*d21 + 6*S_direct_1_b(1,k-1)*d20)*dr_diff_exp);
%        
%        S_direct_4_b(1,k)   = S_direct_4_b(1,k-1) +... 
%            ((4*S_direct_4_b(1,k-1)*(d11+3*d22) + 12*S_direct_3_b(1,k-1)*d21 + 12*S_direct_2_b(1,k-1)*d20)*dr_diff_exp);
%        
%        S_direct_5_b(1,k)   = S_direct_5_b(1,k-1) +... 
%            ((5*S_direct_5_b(1,k-1)*(d11+4*d22) + 20*S_direct_4_b(1,k-1)*d21 + 20*S_direct_3_b(1,k-1)*d20)*dr_diff_exp);
%        
%        S_direct_6_b(1,k)   = S_direct_6_b(1,k-1) +... 
%            ((6*S_direct_6_b(1,k-1)*(d11+5*d22) + 30*S_direct_5_b(1,k-1)*d21 + 30*S_direct_4_b(1,k-1)*d20)*dr_diff_exp);


% mit S1
%        S_direct_1(1,k)   = S_direct_1(1,k-1) -... 
%            (S_direct_1(1,k-1)*d11*dr_diff); 
%        
%        S_direct_2(idx,k)   = S_direct_2(idx,k-1) -... 
%            ((2*S_direct_2(idx,k-1)*(d11+d22) + 2*S_direct_1(1,k-1)*d21 + 2*d20)*dr_diff);
%       
%        S_direct_3(idx,k)   = S_direct_3(idx,k-1) -... 
%            ((3*S_direct_3(idx,k-1)*(d11+2*d22) + 6*S_direct_2(idx,k-1)*d21 + 6*S_direct_1(1,k-1)*d20)*dr_diff);
%        
%        S_direct_4(idx,k)   = S_direct_4(idx,k-1) -... 
%            ((4*S_direct_4(idx,k-1)*(d11+3*d22) + 12*S_direct_3(idx,k-1)*d21 + 12*S_direct_2(idx,k-1)*d20)*dr_diff);
%        
%        S_direct_5(idx,k)   = S_direct_5(idx,k-1) -... 
%            ((5*S_direct_5(idx,k-1)*(d11+4*d22) + 20*S_direct_4(idx,k-1)*d21 + 20*S_direct_3(idx,k-1)*d20)*dr_diff);
%        
%        S_direct_6(idx,k)   = S_direct_6(idx,k-1) -... 
%            ((6*S_direct_6(idx,k-1)*(d11+5*d22) + 30*S_direct_5(idx,k-1)*d21 + 30*S_direct_4(idx,k-1)*d20)*dr_diff);

% -----------------------------------------------------------------------------------------------------










%%
% h(2) = figure;
% set(gcf, 'Color', 'w')
% subplot(1,5,1)
% plot(r_diff_exp(1:end-1),diff(S_exp_2)./-mean(diff(r_diff_exp)),'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
% hold on
% plot(r_diff,S_direct_2,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 
% % plot(r_diff,S_direct_2,'LineWidth',line) 
% set(gca,'xScale','log') 
%  set(gca,'yScale','log') 
%  set(gcf, 'Color', 'w')
%  set(gca,'FontSize',fontsize)
%  xlabel('$r / \lambda$', 'interpreter','latex')
%  ylabel('$S^2(r)$', 'interpreter','latex')
%  xlim([min(r_diff_exp) max(r_diff_exp)])
% %  ylim([min(S_exp_2) max(S_exp_2)])
%  vline(evaluated(1).r/taylor_L,'r')
%  vline(taylor_L/taylor_L,'k')
%  vline(int_L/taylor_L,'k')
%   axis square
% 
% subplot(1,5,2)
% plot(r_diff_exp,-S_exp_3,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
%  hold on
%  plot(r_diff,-S_direct_3,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 
% %  plot(r_diff,-S_direct_3,'LineWidth',line) 
%  set(gca,'xScale','log')
%  set(gca,'yScale','log') 
%  set(gcf, 'Color', 'w')
%  set(gca,'FontSize',fontsize)
%  xlabel('$r / \lambda$', 'interpreter','latex')
%  ylabel('$-S^3(r)$', 'interpreter','latex')
%  xlim([min(r_diff_exp) max(r_diff_exp)])
% %  ylim([min(-S_exp_3) max(-S_exp_3)])
%  vline(evaluated(1).r/taylor_L,'r')
%  vline(taylor_L/taylor_L,'k')
%  vline(int_L/taylor_L,'k')
%  axis square
% 
%  subplot(1,5,3)
% plot(r_diff_exp,S_exp_4,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
%  hold on
%  plot(r_diff,S_direct_4,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
% %  plot(r_diff,S_direct_4,'LineWidth',line)
%  set(gca,'xScale','log')
%  set(gca,'yScale','log') 
%  set(gcf, 'Color', 'w')
%  set(gca,'FontSize',fontsize)
%  xlabel('$r / \lambda$', 'interpreter','latex')
%  ylabel('$S^4(r)$', 'interpreter','latex')
%  xlim([min(r_diff_exp) max(r_diff_exp)])
% %  ylim([min(S_exp_4) max(S_exp_4)])
%  vline(evaluated(1).r/taylor_L,'r')
%  vline(taylor_L/taylor_L,'k')
%  vline(int_L/taylor_L,'k')
%  axis square
%  
%  subplot(1,5,4)
% plot(r_diff_exp,-S_exp_5,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
%  hold on
%  plot(r_diff,-S_direct_5,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 
% %   plot(r_diff,-S_direct_5,'LineWidth',line)
%  set(gca,'xScale','log')
%  set(gca,'yScale','log') 
%  set(gcf, 'Color', 'w')
%  set(gca,'FontSize',fontsize)
%  xlabel('$r / \lambda$', 'interpreter','latex')
%  ylabel('$-S^5(r)$', 'interpreter','latex')
%  xlim([min(r_diff_exp) max(r_diff_exp)])
% %  ylim([min(-S_exp_5) max(-S_exp_5)])
%  vline(evaluated(1).r/taylor_L,'r')
%  vline(taylor_L/taylor_L,'k')
%  vline(int_L/taylor_L,'k')
%  axis square
%  
% subplot(1,5,5)
% plot(r_diff_exp,S_exp_6,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
%  hold on
%  plot(r_diff,S_direct_6,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
% %  plot(r_diff,S_direct_6,'LineWidth',line)
%  set(gca,'xScale','log')
%  set(gca,'yScale','log') 
%  set(gcf, 'Color', 'w')
%  set(gca,'FontSize',fontsize)
%  xlabel('$r / \lambda$', 'interpreter','latex')
%  ylabel('$S^6(r)$', 'interpreter','latex')
%  xlim([min(r_diff_exp) max(r_diff_exp)])
% %  ylim([min(S_exp_6) max(S_exp_6)])
%  vline(evaluated(1).r/taylor_L,'r')
%  vline(taylor_L/taylor_L,'k')
%  vline(int_L/taylor_L,'k')
%  axis square 
%  set(h(2),'Position',[1 923 2560 415]); 


















% ylabel('$d_{11}(r)*(r / \lambda)$', 'interpreter','latex')
% 
% 
% 
% plot(r_diff,S_direct_3,'LineWidth',line)
% 
% 
% 
% figure
% plot(r_diff,-S_exp_3,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
% hold on
% plot(r_diff,-S_direct_3,'LineWidth',line,'Color',[0    0.4470    0.7410]);
% plot(r_diff,-S_direct_3b,'LineWidth',line,'Color',[ 0.8500    0.3250    0.0980]);
% 
% plot(r_diff,a,'LineWidth',line)
% plot(r_diff,b,'LineWidth',line)
% plot(r_diff,c,'LineWidth',line)
% set(gca,'xScale','log')
% axis square
% set(gcf, 'Color', 'w')
% set(gca,'FontSize',fontsize)
% xlabel('$r / \lambda$', 'interpreter','latex')
% xlim([min(r/taylor_L) max(r/taylor_L)]);
% vline(evaluated(1).r/taylor_L,'r')
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
% 
% figure
% plot(r_diff,a,'LineWidth',line)
% hold on
% plot(r_diff,b,'LineWidth',line)
% plot(r_diff,c,'LineWidth',line)
% set(gca,'xScale','log')
% axis square
% set(gcf, 'Color', 'w')
% set(gca,'FontSize',fontsize)
% xlabel('$r / \lambda$', 'interpreter','latex')
% xlim([min(r/taylor_L) max(r/taylor_L)]);
% vline(evaluated(1).r/taylor_L,'r')
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
% 
% 
% 
% 
% 
% 
% 
% % Renner
%         d10 = -0.0015.*(r_diff(1:end-1)).^0.6;
%         d11 =   -0.61.*(r_diff(1:end-1)).^(-0.67);
%         d12 =  0.0096.*(r_diff(1:end-1)).^(0);
%         d13 = -0.0023.*(r_diff(1:end-1)).^(0.3);
%         
%         d20 =   0.033.*(r_diff(1:end-1)).^(0.25);
%         d21 =  -0.009.*(r_diff(1:end-1)).^(0.2);
%         d22 =   0.043.*(r_diff(1:end-1)).^(-0.73);
% 
% figure
% plot(r_diff(1:end-1).*taylor_L,d11.*r_diff(1:end-1),'LineWidth',line)
% hold on
% 
%         
% % co=co_IFT_opti;
% % co=co_KM_opti;
% % co=co_KM_opti_no_offset;
% co.a=[-0.33 0];
% co.ea=-1;
% co.b=[0 0 0 0 0 0];
% co.eb=[0 0 0];
% co.b=[0 0 0 0 0.0167 0];
% co.eb=[0 0 -1];        
%         
% d11 = (co.a(1).*(r_diff(1:end-1)).^co.ea(1)+co.a(2));
% d20 = (co.b(1).*(r_diff(1:end-1)).^co.eb(1)+co.b(2));
% d21 = (co.b(3).*(r_diff(1:end-1)).^co.eb(2)+co.b(4));
% d22 = (co.b(5).*(r_diff(1:end-1)).^co.eb(3)+co.b(6)); 
% 
% plot(r_diff(1:end-1),d11.*r_diff(1:end-1),'LineWidth',line)
%       
%       
% figure
% plot(r_diff(1:end-1),d22.*r_diff(1:end-1),'LineWidth',line)
% hold on
%    
% 
% 
% plot(r_diff(1:end-1),d20,'LineWidth',line)
% plot(r_diff(1:end-1),d21,'LineWidth',line)
% plot(r_diff(1:end-1),d22,'LineWidth',line)
% set(gca,'xScale','log')
% axis square
% set(gcf, 'Color', 'w')
% set(gca,'FontSize',fontsize)
% xlabel('$r / \lambda$', 'interpreter','latex')
% xlim([min(r/taylor_L) max(r/taylor_L)]);
% vline(evaluated(1).r/taylor_L,'r')
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
% 
% figure
% plot(r_diff(1:end-1),d11./(r_diff(1:end-1)),'LineWidth',line)
% hold on
% plot(r_diff(1:end-1),d20./(r_diff(1:end-1)),'LineWidth',line)
% plot(r_diff(1:end-1),d21./(r_diff(1:end-1)),'LineWidth',line)
% plot(r_diff(1:end-1),d22./(r_diff(1:end-1)),'LineWidth',line)
% set(gca,'xScale','log')
% axis square
% set(gcf, 'Color', 'w')
% set(gca,'FontSize',fontsize)
% xlabel('$r / \lambda$', 'interpreter','latex')
% xlim([min(r/taylor_L) max(r/taylor_L)]);
% vline(evaluated(1).r/taylor_L,'r')
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
% 
% 
% figure
% plot(r_diff,S_exp_6,'MarkerSize',8,'Marker','+','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
% hold on
% % plot(r_diff,S_direct_6,'LineWidth',line,'Color',[0    0.4470    0.7410]);
% plot(r_diff,S_direct_6b,'LineWidth',line,'Color',[ 0.8500    0.3250    0.0980]);
% set(gca,'xScale','log')
% set(gca,'yScale','log') 
% axis square
% set(gcf, 'Color', 'w')
% set(gca,'FontSize',fontsize)
% xlabel('$r / \lambda$', 'interpreter','latex')
% ylabel('$S^6(r)$', 'interpreter','latex')
% xlim([min(r/taylor_L) max(r/taylor_L)]);
% ylim([min(S_exp_6) max(S_exp_6)])
% vline(evaluated(1).r/taylor_L,'r')
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
end