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
function [r,S_exp_2,S_exp_3,S_exp_4,S_exp_5,S_exp_6,S_exp_7,Lam_sq_L_1,Lam_sq_L_2,Lam_sq_L_3,T_exp_3,T_exp_5,T_exp_7,exponent]=plot_struc_function(r_n,Fs,data,m_data,int_L,taylor_L,mu,D,norm_r,save_path,save_name)
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

T_exp_3(1:length(stepp))               = nan; 
T_exp_5(1:length(stepp))               = nan;  
T_exp_7(1:length(stepp))               = nan;  

K(1:length(stepp))            = nan;  
Lam_sq_L_1(1:length(stepp))            = nan;  
Lam_sq_L_2(1:length(stepp))            = nan;  
Lam_sq_L_3(1:length(stepp))            = nan;  

Skew(1:length(stepp))               = nan; 

tmp = length(data);

parfor k=1:length(stepp)
    k/length(stepp)   
    tau1            = stepp(k); 
    r(k)            = (tau1./Fs).*m_data;
    
    incr1   =   (data(tau1+1:tmp)     - data(1:tmp-tau1)).';
    
%% EXP    
    K(k)                    =         kurtosis(incr1);         
    Lam_sq_L_1(k)  = log(kurtosis(incr1)/3)/4;
    Lam_sq_L_2(k)  = log(kurtosis(incr1))/4;

    Skew(k)  = skewness(incr1);
    
    S_exp_2(k)  = mean(incr1.^2);
    S_exp_3(k)  = mean(incr1.^3);
    S_exp_4(k)  = mean(incr1.^4);
    S_exp_5(k)  = mean(incr1.^5);
    S_exp_6(k)  = mean(incr1.^6);
    S_exp_7(k)  = mean(incr1.^7);
    
    T_exp_3(k)  = mean(abs(incr1).^3);
    T_exp_5(k)  = mean(abs(incr1).^5);
    T_exp_7(k)  = mean(abs(incr1).^7);
       
    Lam_sq_L_3(k)   = log(S_exp_4(k)/(3*S_exp_2(k)*S_exp_2(k)));
end

%% Plot structure functions
in_n_nan   = ~isnan(r);
r          = r(in_n_nan);
if norm_r==1
    r          = r./taylor_L;
end
S_exp_2    = S_exp_2(in_n_nan); 
S_exp_3    = S_exp_3(in_n_nan);
S_exp_4    = S_exp_4(in_n_nan); 
S_exp_5    = S_exp_5(in_n_nan);
S_exp_6    = S_exp_6(in_n_nan);
S_exp_7    = S_exp_7(in_n_nan);

T_exp_3    = T_exp_3(in_n_nan);
T_exp_5    = T_exp_5(in_n_nan);
T_exp_7    = T_exp_7(in_n_nan);

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
if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
else
    xlabel('$r\ (m)$', 'interpreter','latex')
end
ylabel('$S^2$', 'interpreter','latex')
xlim([min(r)*0.5 max(r)*2])
ylim([min(S_exp_2) max(S_exp_2)])
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
axis square
txt = {'(a)'};
a=text(4,0.5,txt);
a.Units='normalized';
fig_setup
legend off

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
if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
else
    xlabel('$r\ (m)$', 'interpreter','latex')
end
ylabel('$-S^3$', 'interpreter','latex')
xlim([min(r)*0.5 max(r)*2])
ylim([min(-S_exp_3) max(-S_exp_3)])
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
axis square
txt = {'(b)'};
b=text(4,0.5,txt);
b.Units='normalized';
fig_setup
legend off

subplot(3,2,3)
plot(r,S_exp_4,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
set(gca,'xScale','log')
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
% xlabel('$r$', 'interpreter','latex')
xlabel('$r/m$', 'interpreter','latex')
if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
else
    xlabel('$r\ (m)$', 'interpreter','latex')
end
ylabel('$S^4$', 'interpreter','latex')
xlim([min(r)*0.5 max(r)*2])
ylim([min(S_exp_4) max(S_exp_4)])
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
axis square
txt = {'(c)'};
c=text(4,0.5,txt);
c.Units='normalized';
fig_setup
legend off

subplot(3,2,4)
plot(r,-S_exp_5,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
set(gca,'xScale','log')
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
% xlabel('$r$', 'interpreter','latex')
xlabel('$r/m$', 'interpreter','latex')
if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
else
    xlabel('$r\ (m)$', 'interpreter','latex')
end
ylabel('$-S^5$', 'interpreter','latex')
xlim([min(r)*0.5 max(r)*2])
ylim([min(-S_exp_5) max(-S_exp_5)])
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
axis square
txt = {'(d)'};
d=text(4,0.5,txt);
d.Units='normalized';
fig_setup
legend off

subplot(3,2,5)
plot(r,S_exp_6,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
set(gca,'xScale','log')
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
% xlabel('$r$', 'interpreter','latex')
xlabel('$r/m$', 'interpreter','latex')
if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
else
    xlabel('$r\ (m)$', 'interpreter','latex')
end
ylabel('$S^6$', 'interpreter','latex')
xlim([min(r)*0.5 max(r)*2])
ylim([min(S_exp_6) max(S_exp_6)])
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
axis square 
txt = {'(e)'};
e=text(4,0.5,txt);
e.Units='normalized';
fig_setup
legend off

subplot(3,2,6)
plot(r,-S_exp_7,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
set(gca,'xScale','log')
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
xlabel('$r / \lambda$', 'interpreter','latex')
if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
else
    xlabel('$r\ (m)$', 'interpreter','latex')
end
ylabel('$-S^7$', 'interpreter','latex')
xlim([min(r)*0.5 max(r)*2])
ylim([min(-S_exp_7) max(-S_exp_7)])
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
axis square
% set(h(1),'Position',[1 46 835 1299]);
% fig_setup
% legend off
txt = {'(f)'};
f=text(4,0.5,txt);
f.Units='normalized';
fig_setup
legend off

font=22;
a.FontSize = font;
b.FontSize = font;
c.FontSize = font;
d.FontSize = font;
e.FontSize = font;
f.FontSize = font;

pos_txt=[-0.29   0.9];
a.Position=pos_txt;
b.Position=pos_txt;
c.Position=pos_txt;
d.Position=pos_txt;
e.Position=pos_txt;
f.Position=pos_txt;


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
if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
else
    xlabel('$r\ (m)$', 'interpreter','latex')
end
xlim([min(r)*0.5 max(r)*2])
ylim([0.5*min([S_exp_2 S_exp_4 S_exp_6]) 2*max([S_exp_2 S_exp_4 S_exp_6])])
legend([p1,p2,p3],{'$S^2(r)$','$S^4(r)$','$S^6(r)$'},'Interpreter','latex','Location','northwest','FontSize',18)  
legend([p1,p2,p3],{'$S^2$','$S^4$','$S^6$'},'Interpreter','latex','Location','northwest','FontSize',18)  
txt = {'(a)'};
a=text(4,0.5,txt);
a.Units='normalized';
fig_setup


subplot(1,2,2)
p4=plot(r,-S_exp_3,'MarkerSize',8,'Marker','+','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
hold on
p5=plot(r,-S_exp_5,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
p6=plot(r,-S_exp_7,'MarkerSize',8,'Marker','*','LineWidth',2,'LineStyle','none','Color',[0 0 0]); 

p7=plot(r,T_exp_3,'MarkerSize',8,'Marker','+','LineWidth',2,'LineStyle','none');
p8=plot(r,T_exp_5,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none');
p9=plot(r,T_exp_7,'MarkerSize',8,'Marker','*','LineWidth',2,'LineStyle','none'); 

plot(r,feval(fitresult,r),'LineWidth',2,'LineStyle','--','Color','r')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
set(gca,'xScale','log')
set(gca,'yScale','log') 
axis square
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
xlabel('$r/m$', 'interpreter','latex')
if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
else
    xlabel('$r\ (m)$', 'interpreter','latex')
end
xlim([min(r)*0.5 max(r)*2])
ylim([0.5*min([-S_exp_3 -S_exp_5 -S_exp_7]) 2*max([-S_exp_3 -S_exp_5 -S_exp_7])])
legend([p4,p5,p6,p7,p8,p9],{'$-S^3(r)$','$-S^5(r)$','$-S^7(r)$','$T^3(r)$','$T^5(r)$','$T^7(r)$'},'Interpreter','latex','Location','northwest','FontSize',18)  
legend([p4,p5,p6,p7,p8,p9],{'$-S^3$','$-S^5$','$-S^7$','$T^3$','$T^5$','$T^7$'},'Interpreter','latex','Location','northwest','FontSize',18)  
% vline(taylor_L,'k')
% vline(int_L,'k')
% vline(taylor_L/taylor_L,'k')
% vline(int_L/taylor_L,'k')
txt = {'(b)'};
b=text(4,0.5,txt);
b.Units='normalized';
fig_setup


pos_txt=[-0.25   0.9];
a.Position=pos_txt;
b.Position=pos_txt;
font=22;
a.FontSize = font;
b.FontSize = font;


exponent=nan(1,7);
% x=log(-S_exp_3);
x=log(T_exp_3);


h(3) = figure;
subplot(1,2,1)
y=log(S_exp_2);
p1=plot(x,y,'MarkerSize',8,'Marker','+','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
hold on
[xData, yData] = prepareCurveData( x, y );
ft = fittype( 'poly1' );
opts = fitoptions( ft );
[fitresult, gof] = fit( xData, yData, ft, opts );
exponent(2)=fitresult.p1;
plot(x,feval(fitresult,x),'LineWidth',2,'Color',[0 0 0])


y=log(S_exp_4);
p2=plot(x,y,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
[xData, yData] = prepareCurveData( x, y );
ft = fittype( 'poly1' );
opts = fitoptions( ft );
[fitresult, gof] = fit( xData, yData, ft, opts );
exponent(4)=fitresult.p1;
plot(x,feval(fitresult,x),'LineWidth',2,'Color',[0 0 0])


y=log(S_exp_6);
p3=plot(x,y,'MarkerSize',8,'Marker','*','LineWidth',2,'LineStyle','none','Color',[0 0 0]); 
[xData, yData] = prepareCurveData( x, y );
ft = fittype( 'poly1' );
opts = fitoptions( ft );
[fitresult, gof] = fit( xData, yData, ft, opts );
exponent(6)=fitresult.p1;
plot(x,feval(fitresult,x),'LineWidth',2,'Color',[0 0 0])
axis square
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
xlabel('$log(T^3)$', 'interpreter','latex')
legend('$log(S^2)$','$log(S^4)$','$log(S^6)$','Interpreter','latex','Location','northwest','FontSize',18)  
txt = {'(a)'};
a=text(4,0.5,txt);
a.Units='normalized';
fig_setup



subplot(1,2,2)
% y=log(-S_exp_5);
y=log(T_exp_5);
p1=plot(x,y,'MarkerSize',8,'Marker','+','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
hold on
[xData, yData] = prepareCurveData( x, y );
ft = fittype( 'poly1' );
opts = fitoptions( ft );
[fitresult, gof] = fit( xData, yData, ft, opts );
exponent(5)=fitresult.p1;
plot(x,feval(fitresult,x),'LineWidth',2,'Color',[0 0 0])

% y=log(-S_exp_7);
y=log(T_exp_7);
p2=plot(x,y,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
[xData, yData] = prepareCurveData( x, y );
[fitresult, gof] = fit( xData, yData, ft, opts );
exponent(7)=fitresult.p1;
plot(x,feval(fitresult,x),'LineWidth',2,'Color',[0 0 0])
axis square
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
xlabel('$log(T^3)$', 'interpreter','latex')
legend('$log(T^5)$','$log(T^6)$','Interpreter','latex','Location','northwest','FontSize',18)  
txt = {'(b)'};
b=text(4,0.5,txt);
b.Units='normalized';
fig_setup


pos_txt=[-0.25   0.9];
a.Position=pos_txt;
b.Position=pos_txt;
font=22;
a.FontSize = font;
b.FontSize = font;




h(4) = figure;
x_tmp=0:8;
plot(exponent,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
hold on
plot(x_tmp,x_tmp./3,'LineWidth',2,'LineStyle','--')
plot(x_tmp,x_tmp./3-((mu*(x_tmp.^2-3*x_tmp))/18),'LineWidth',2,'LineStyle','--')
plot(x_tmp,x_tmp./9+2*(1-(2/3).^(x_tmp./3)),'LineWidth',2,'LineStyle','--')
plot(x_tmp,x_tmp./3+(3-D).*(1-(x_tmp./3)),'LineWidth',2,'LineStyle','--')
axis square
set(gcf, 'Color', 'w')
set(gca,'FontSize',fontsize)
xlabel('$k$', 'interpreter','latex')
ylabel('$\zeta_k$', 'interpreter','latex')
legend({'ESS: exp','$K41$','$K62$','She-Leveque','$\beta$-Model'},'Interpreter','latex','Location','northwest','FontSize',18)  




h(5) = figure;
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
if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
else
    xlabel('$r\ (m)$', 'interpreter','latex')
end
% xlim([min(r)*0.5 max(r)*2])
legend({'$\frac{ln\left(K(u_r)/3\right)}{4}$',...
    '$\frac{ln\left(K(u_r)\right)}{4}$',...
    '$ln\left(\frac{S^4(r)}{3\ \left(S^2(r)\right)^2}\right)$'},'Interpreter','latex','Location','northeast','FontSize',18)  
legend({'$\frac{ln\left(K(u_r)/3\right)}{4}$',...
    '$\frac{ln\left(K(u_r)\right)}{4}$',...
    '$ln\left(\frac{S^4}{3\ \left(S^2\right)^2}\right)$'},'Interpreter','latex','Location','northeast','FontSize',18)  
fig_setup
txt = {'(a)'};
a=text(4,0.5,txt);
a.Units='normalized';


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
if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
else
    xlabel('$r\ (m)$', 'interpreter','latex')
end
xlim([min(r)*0.5 max(r)*2])
legend({'$\frac{ln\left(K(u_r)/3\right)}{4}$',...
    '$\frac{ln\left(K(u_r)\right)}{4}$',...
    '$ln\left(\frac{S^4(r)}{3\ \left(S^2(r)\right)^2}\right)$'},'Interpreter','latex','Location','northeast','FontSize',18)  
legend({'$\frac{ln\left(K(u_r)/3\right)}{4}$',...
    '$\frac{ln\left(K(u_r)\right)}{4}$',...
    '$ln\left(\frac{S^4}{3\ \left(S^2\right)^2}\right)$'},'Interpreter','latex','Location','northeast','FontSize',18)  
fig_setup
txt = {'(b)'};
b=text(4,0.5,txt);
b.Units='normalized';


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
if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
else
    xlabel('$r\ (m)$', 'interpreter','latex')
end
xlim([min(r)*0.5 max(r)*2])
legend({'$\frac{ln\left(K(u_r)/3\right)}{4}$',...
    '$\frac{ln\left(K(u_r)\right)}{4}$',...
    '$ln\left(\frac{S^4(r)}{3\ \left(S^2(r)\right)^2}\right)$'},'Interpreter','latex','Location','northeast','FontSize',18)  
legend({'$\frac{ln\left(K(u_r)/3\right)}{4}$',...
    '$\frac{ln\left(K(u_r)\right)}{4}$',...
    '$ln\left(\frac{S^4}{3\ \left(S^2\right)^2}\right)$'},'Interpreter','latex','Location','northeast','FontSize',18)  
fig_setup
txt = {'(c)'};
c=text(4,0.5,txt);
c.Units='normalized';
font=22;
a.FontSize = font;
b.FontSize = font;
c.FontSize = font;
pos_txt=[-0.25   0.9];
a.Position=pos_txt;
b.Position=pos_txt;
c.Position=pos_txt;
% fig_setup






% 
% close all
% figure
% % r          = r./taylor_L;
% plot(stepp./markov,K,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none');
% axis square
% set(gcf, 'Color', 'w')
% set(gca,'FontSize',20)
% fig_setup
% legend off
% hold on
% set(gca,'xScale','log')
% set(gca,'yScale','log') 
% vline(1,'k')
% % vline(60,'k')
% plot(stepp./markov,3.*ones(numel(stepp),1),'k','LineStyle','--','LineWidth',2);
% xlabel('$\Delta \tau/\Delta_{EM}$', 'interpreter','latex')
% xlim([min(stepp./markov) max(stepp./markov)])
% 
% 
% 
% figure
% % r          = r./taylor_L;
% plot(stepp*0.16666,K,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none');
% axis square
% set(gcf, 'Color', 'w')
% set(gca,'FontSize',20)
% fig_setup
% legend off
% hold on
% set(gca,'xScale','log')
% set(gca,'yScale','log') 
% vline(markov*0.16666,'k')
% % vline(60*0.16666,'k')
% plot(stepp*0.16666,3.*ones(numel(stepp),1),'k','LineStyle','--','LineWidth',2);
% xlabel('$\Delta \tau/\tau_\eta$', 'interpreter','latex')
% xlim([min(stepp*0.16666) max(stepp*0.16666)])

if ischar(save_path)
   savefig(h,fullfile(save_path,append(save_name,'_','structure_function.fig')),'compact')
%    for a = 1:length(h)
%         exportgraphics(h(a),fullfile(save_path,append(save_name,'_',sprintf('structure_function_%d.png',
%         a))))
%    end
end

for i=5:-1:1
    tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
    uiwait(gcf);
    if ischar(save_path)
        exportgraphics(h(i),fullfile(save_path,append(save_name,'_',sprintf('structure_function_%d.png', i))))
    end
    delete(tmp_ui);
    close
end
end