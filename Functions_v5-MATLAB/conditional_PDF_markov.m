%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The main work of this function is to plot the single conditioned and double conditioned PDF's and 
% check visually if the markovian property is fulfilled or not. This is a qualitative check.
% This function calculates the single conditioned and double conditioned PDF of velocity increments 
% for a pair of two scales each of which is separated by the $\Delta_{EM}$. To do this, a pop-up 
% dialog box is generated to enter the conditioned value for large scale increment $u_{r_{3}}$, for 
% example $u_{r_{3}}=\pm 1$. Note, the condition $u_{r_{3}}=0$ corresponds to the maximum number of 
% statistics. This function also plots various representations of the single and double conditioned 
% PDF's. This is a qualitative/visual check for the validation of Markov property based on the 
% alignment or misalignment of these single and double conditioned PDF's. If there is not a good 
% agreement between the single conditioned and double conditioned PDF of velocity increments, it is 
% possible to modify the Einstein-Markov length and/or the minimum number of events and repeat this 
% qualitative check.
%
% Arguments IN
% data = it is the filtered data 
% Fs = Acquisition/Sampling Frequency in Hz
% m_data = mean of the filtered data
% markov = markov length in samples/steps NOT in meters
% second_condi = increment of the second condition
% num_bin = Number of bins 
% min_events = minimum number of events to consider in a single bin
% norm_ur = normalization of the scale using $\lambda$? data 1=Yes, 0=No
% save_path = path for saving figures and files
% save_name = name for saving files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function conditional_PDF_markov(data,Fs,m_data,markov,second_condi,num_bin,min_events,norm_ur,save_path,save_name)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');

tau2    = markov;       %in steps
tau1    = 2*tau2;       %in steps
tau0    = 3*tau2;       %in steps

r2      = tau2*m_data/Fs; %in meters
r1      = tau1*m_data/Fs; %in meters
r0      = tau0*m_data/Fs; %in meters

%tau2<tau1<tau0
%%r2<r1<r0
[incr1,incr2,incr3]    = Increment_double_condition(tau2,tau1,tau0,data);
[dx,dy,dz,x,y,z,dA,dV] = limit_double_condition(incr1,incr2,incr3,num_bin);


% Ideally, tol = 0 but after/during data processing you will never end up with increments(incr3) 
% which are exactly equal to zero. So in practice we must put some tolerance on incr3 tolerance(tol) 
% can be changes accordingly
% tol            = range(z)*0.01; 
tol                     = 0.1;


%% Plot PDF of increments
% nbins=1001;
% [fuL uL] = hist(incr1,nbins);                 % compute frequency of u-values at scale r=L
% puL = fuL / ( sum(fuL) * mean(diff(uL)) );    % relative frequency to approximate p(u,r=L)
% figure
% plot(uL,puL)
% uL(find(puL==max(puL)))
% max(puL)
% hold on
% [fuL uL] = hist(incr2,nbins); 
% puL = fuL / ( sum(fuL) * mean(diff(uL)) ); 
% plot(uL,puL)
% uL(find(puL==max(puL)))
% max(puL)
% [fuL uL] = hist(incr3,nbins);
% puL = fuL / ( sum(fuL) * mean(diff(uL)) );
% plot(uL,puL,'k')
% uL(find(puL==max(puL)))
% max(puL)
% axis square      
% set(gca, 'FontSize',16,'PlotBoxAspectRatio',[1 1 1])
% set(gca,'yScale','log') 
% set(gcf, 'Color', 'w')
% xPos = 0.0;
% plot([xPos xPos], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% single conditional probability
[P_AIB,P_BIA,P_AnB,P_A,P_B,binP_AIB,x_mean_bin,y_mean_bin,events_1,counter_A,counter_B] = distribution(tau2,tau1,incr1,incr2,num_bin,dx,dy,x,y,dA);

P_AIB_div	= P_AIB;
P_AIB_div(counter_A<min_events,counter_B<min_events)=0;

%cut out nan's and min_events
y_input     = y_mean_bin(counter_B>=min_events);
x_input     = x_mean_bin(counter_A>=min_events);

P_AIB      = P_AIB(counter_A>=min_events,counter_B>=min_events);
P_AnB     = P_AnB(counter_A>=min_events,counter_B>=min_events);
P_A         = P_A(counter_A>=min_events);
P_B         = P_B(counter_B>=min_events);
                         
level = round(logspace(log10(0.01*max(max((P_AIB)))*10^6),log10(0.9*max(max((P_AIB)))*10^6),7))./10^6;

fig(1)=figure;
[C,h]=contour(y_input,x_input,P_AIB,level,'-','color','k','LineWidth',2,'ShowText','off');
% [C,h]=contour3(y_input,x_input,P_AIB,level,'-','color','k','LineWidth',2,'ShowText','off');
% [C,h]=contour3(y_input,x_mean_bin,P_AIB,50,'-','color','k','LineWidth',2,'ShowText','off');
set(gca,'FontSize',18)
% xlabel(['$u_r (r_1=' num2str(r2) 'm) / \sigma_\infty$'], 'interpreter','latex')
% ylabel(['$u_r (r_2=' num2str(r1) 'm) / \sigma_\infty$'], 'interpreter','latex')
if norm_ur==1  
xlabel(['$u_{r_1=' num2str(tau1/markov) '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
ylabel(['$u_{r_2=' num2str(tau2/markov) '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
xlabel(['$u_{r_1=' num2str(tau1/markov) '\Delta_{EM}}$'], 'interpreter','latex')
ylabel(['$u_{r_2=' num2str(tau2/markov) '\Delta_{EM}}$'], 'interpreter','latex')
end
% axis([-5 5 -5 5])
axis square
set(gcf, 'Color', 'w')
hold on    
    



%% double_conditional Probability
% cut for the second conditon (incr3) 
index = find(incr3>-tol+second_condi & incr3<tol+second_condi);
[P_AIB_2,P_BIA_2,P_AnB_2,P_A_2,P_B_2,binP_AIB_2,x_mean_bin_2,y_mean_bin_2,events_2,counter_A_2,counter_B_2] = distribution_double(tau2,tau1,incr1(index),incr2(index),num_bin,dx,dy,x,y,dA);   
P_AIB_div_2=P_AIB_2;
P_AIB_div_2(counter_A_2<min_events,counter_B_2<min_events)=0;

% cut out nan's and min_events
y_input_2=y_mean_bin_2(counter_B_2>=min_events);
x_input_2=x_mean_bin_2(counter_A_2>=min_events);

P_AIB_2=P_AIB_2(counter_A_2>=min_events,counter_B_2>=min_events);
P_AnB_2=P_AnB_2(counter_A_2>=min_events,counter_B_2>=min_events);
P_A_2=P_A_2(counter_A_2>=min_events);
P_B_2=P_B_2(counter_B_2>=min_events);

% plot double_condition PDF
[C,h]=contour(y_input_2,x_input_2,P_AIB_2,level,'-','color','r','LineWidth',2,'ShowText','off'); 
% event_index_2=events_2>=min_events;
% [C,h]=contour(y_mean_bin_2(event_index_2),x_mean_bin_2,(P_AIB_2(:,event_index_2)),level,'-','color','r','LineWidth',2,'ShowText','off');  
% [C,h]=contour(y_mean_bin_2,x_mean_bin_2,(P_AIB_2),level,'-','color','r','LineWidth',2,'ShowText','off');  

% cut_a=find(y_mean_bin > -1-tol & y_mean_bin < -1+tol,1,'first');
cut_0=find(y_input > 0-tol & y_input < 0+tol,1,'first');
cut_b=find(y_input > 0.5-tol & y_input < 0.5+tol,1,'first');

%vertical cut
plot([y_input(cut_0) y_input(cut_0)], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
% plot([y(cut_a) y(cut_a)], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color',[0    0.4470    0.7410]); % Adapts to y limits of current axes
plot([y_input(cut_b) y_input(cut_b)], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes

legend({'$p(u_{r_2}|u_{r_1})$','$p(u_{r_2}|u_{r_1},u_{r_0})$'}, 'interpreter','latex','Location','northwest')




% Titel hier einfügen
%% weighted mean square error function in logarithmic space; conditional probability densities
d_M_1   = 0;
d_M_2   = 0;
for j =1:num_bin
    for i =1:num_bin
         if (P_AIB_div_2(i,j) > 0) && (P_AIB_div(i,j) > 0)     
            d_M_1 = d_M_1 + ((P_AIB_div_2(i,j) + P_AIB_div(i,j)).*...
                            (log(P_AIB_div_2(i,j))    - log(P_AIB_div(i,j))).^2);                        
            d_M_2 = d_M_2 + ((P_AIB_div_2(i,j) + P_AIB_div(i,j)).*...
                            (log(P_AIB_div_2(i,j)).^2 + log(P_AIB_div(i,j)).^2));
        end
%         if (P_verbund_n(i,j) < 0)
%             d_M_1 = d_M_1 + 10.^(10);                        
%         end
    end
end
fval     = d_M_1/d_M_2;

title(['$\xi_{\Delta_{EM}}$=' num2str(fval,2)])


fig_setup
txt = {'(a)'};
a=text(4,0.5,txt);
a.Units='normalized';
%horizontal cut
% plot(get(gca,'xlim'),[x(cut_0) x(cut_0)],'LineWidth',2,'LineStyle','--','Color','k'); % Adapts to y limits of current axes
% plot(get(gca,'ylim'),get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','r');
% plot(get(gca,'xlim'),[x(cut_a) x(cut_a)],'LineWidth',2,'LineStyle','--','Color',[0    0.4470    0.7410]); % Adapts to y limits of current axes
% plot(get(gca,'xlim'),[x(cut_b) x(cut_b)],'LineWidth',2,'LineStyle','--','Color','r'); % Adapts to y limits of current axes

%vertical cut
% figure
% h = axes;
% plot(x_mean_bin,P_AIB(:,find(y_mean_bin > -1-tol & y_mean_bin < -1+tol,1, 'first')),'Color','k','LineWidth',2) 
% hold on
% plot(x_mean_bin_2,P_AIB_2(:,find(y_mean_bin_2 > -1-tol & y_mean_bin_2 < -1+tol,1, 'first')),'Color','r','LineWidth',2)                 
% view([90 90]);
% set(h,'XDir','reverse');
% xlabel(['$\xi (r_1=' num2str(r1) 'm) / \sigma_\infty$'], 'interpreter','latex')
% set(gca,'FontSize',18)
% xlim([-5 5])
% axis square
% set(gcf, 'Color', 'w')
% set(gca,'yscale','log');
% ylim([10^-4 2])
% title(['$\xi (r_1=' num2str(r2) 'm) / \sigma_\infty$ = ',num2str(y(cut_a),'%1.2f')],'interpreter','latex')


%% Plot miscellaneous
fig(2)=figure;
h = axes;
plot(x_input,P_AIB(:,find(y_input > 0-tol & y_input < 0+tol,1,'first')),'Color','k','LineWidth',2) 
hold on
if ~isempty(find(y_input_2 > 0-tol & y_input_2 < 0+tol,1,'first'))
    plot(x_input_2,P_AIB_2(:,find(y_input_2 > 0-tol & y_input_2 < 0+tol,1,'first')),'Color','r','LineWidth',2)
end
% y_input(find(y_input > 0-tol & y_input < 0+tol,1,'first'))
% y_input_2(find(y_input_2 > 0-tol & y_input_2 < 0+tol,1,'first'))
view([90 90]);
set(h,'XDir','reverse');
if norm_ur==1  
    xlabel(['$u_{r_2=' num2str(tau2/markov) '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
    xlabel(['$u_{r_2=' num2str(tau2/markov) '\Delta_{EM}}$'], 'interpreter','latex')   
end
set(gca,'FontSize',18)
% xlim([-5 5])
axis square
set(gcf, 'Color', 'w')
% set(gca,'yscale','log');
% ylim([10^-4 2])
if norm_ur==1  
title(['$u_{r_1=' num2str(tau1/markov) '\Delta_{EM}}/ \sigma_\infty$ = ',num2str(y_input(cut_0),'%1.2f')],'interpreter','latex')
else
    title(['$u_{r_1=' num2str(tau1/markov) '\Delta_{EM}}$ = ',num2str(y_input(cut_0),'%1.2f')],'interpreter','latex')
end   
legend({'$p(u_{r_2}|u_{r_1})$','$p(u_{r_2}|u_{r_1},u_{r_0})$'},'Position',[0.541676093510218 0.77902142683665 0.257004955836705 0.113571428117298],'interpreter','latex')
fig_setup
ylabel('PDF','interpreter','latex')
txt = {'(c)'};
c=text(4,0.5,txt);
c.Units='normalized';


fig(3)=figure;
h = axes;
plot(x_input,P_AIB(:,find(y_input > 0.5-tol & y_input < 0.5+tol,1, 'first')),'Color','k','LineWidth',2) 
hold on
if ~isempty(find(y_input_2 > 0.5-tol & y_input_2 < 0.5+tol,1, 'first'))
plot(x_input_2,P_AIB_2(:,find(y_input_2 > 0.5-tol & y_input_2 < 0.5+tol,1, 'first')),'Color','r','LineWidth',2) 
end
% y_input(find(y_input > 0.5-tol & y_input < 0.5+tol,1,'first'))
% y_input_2(find(y_input_2 > 0.5-tol & y_input_2 < 0.5+tol,1,'first'))
view([90 90]);
set(h,'XDir','reverse');
if norm_ur==1  
    xlabel(['$u_{r_2=' num2str(tau2/markov) '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
    xlabel(['$u_{r_2=' num2str(tau2/markov) '\Delta_{EM}}$'], 'interpreter','latex')    
end
set(gca,'FontSize',18)
% xlim([-5 5])
axis square
set(gcf, 'Color', 'w')
% set(gca,'yscale','log');
% ylim([10^-4 2])
if norm_ur==1  
    title(['$u_{r_1=' num2str(tau1/markov) '\Delta_{EM}}/ \sigma_\infty$ = ',num2str(y_input(cut_b),'%1.2f')],'interpreter','latex')
else
    title(['$u_{r_1=' num2str(tau1/markov) '\Delta_{EM}}$ = ',num2str(y_input(cut_b),'%1.2f')],'interpreter','latex')  
end
legend({'$p(u_{r_2}|u_{r_1})$','$p(u_{r_2}|u_{r_1},u_{r_0})$'},'Position',[0.541676093510218 0.77902142683665 0.257004955836705 0.113571428117298],'interpreter','latex')
fig_setup
ylabel('PDF','interpreter','latex')
txt = {'(d)'};
d=text(4,0.5,txt);
d.Units='normalized';



fig(4)=figure;
[C,h]=contour3(y_input,x_input,P_AIB,level,'-','color','k','LineWidth',2,'ShowText','off');
set(gca,'FontSize',18)
if norm_ur==1  
xlabel(['$u_{r_1=' num2str(tau1/markov) '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
ylabel(['$u_{r_2=' num2str(tau2/markov) '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
xlabel(['$u_{r_1=' num2str(tau1/markov) '\Delta_{EM}}$'], 'interpreter','latex')
ylabel(['$u_{r_2=' num2str(tau2/markov) '\Delta_{EM}}$'], 'interpreter','latex')    
end
axis square
set(gcf, 'Color', 'w')
hold on    
[C,h]=contour3(y_input_2,x_input_2,P_AIB_2,level,'-','color','r','LineWidth',2,'ShowText','off'); 
grid off
set(gca,'zscale','log');
zlim([min(level) max(level)])
fig_setup
set(gca,'zTick',[10^-2 10^-1  10^0]);
legend({'$p(u_{r_2}|u_{r_1})$','$p(u_{r_2}|u_{r_1},u_{r_0})$'}, 'interpreter','latex','Location','northwest')
rotate3d on
% tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
% uiwait(gcf);
% delete(tmp_ui);
% xlabel(['$u_{r_1=' num2str(tau1/markov) '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex',...
%     'Rotation',20,...
%     'VerticalAlignment','bottom',...
%     'HorizontalAlignment','left')
% ylabel(['$u_{r_2=' num2str(tau2/markov) '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex',...
%     'Rotation',325,...
%     'VerticalAlignment','middle',...
%     'HorizontalAlignment','right')
txt = {'(b)'};
b=text(4,0.5,txt);
b.Units='normalized';



fig(5)=figure;
h = axes;
plot(x_input,P_AIB(:,find(y_input > 0-tol & y_input < 0+tol,1,'first')),'Color','k','LineWidth',2) 
hold on
if ~isempty(find(y_input_2 > 0-tol & y_input_2 < 0+tol,1,'first'))
    plot(x_input_2,P_AIB_2(:,find(y_input_2 > 0-tol & y_input_2 < 0+tol,1,'first')),'Color','r','LineWidth',2)
end
% y_input(find(y_input > 0-tol & y_input < 0+tol,1,'first'))
% y_input_2(find(y_input_2 > 0-tol & y_input_2 < 0+tol,1,'first'))
view([90 90]);
set(h,'XDir','reverse');
if norm_ur==1  
    xlabel(['$u_{r_2=' num2str(tau2/markov) '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
    xlabel(['$u_{r_2=' num2str(tau2/markov) '\Delta_{EM}}$'], 'interpreter','latex')   
end
set(gca,'FontSize',18)
% xlim([-5 5])
axis square
set(gcf, 'Color', 'w')
set(gca,'yscale','log');
% ylim([10^-4 2])
if norm_ur==1  
title(['$u_{r_1=' num2str(tau1/markov) '\Delta_{EM}}/ \sigma_\infty$ = ',num2str(y_input(cut_0),'%1.2f')],'interpreter','latex')
else
    title(['$u_{r_1=' num2str(tau1/markov) '\Delta_{EM}}$ = ',num2str(y_input(cut_0),'%1.2f')],'interpreter','latex')
end   
legend({'$p(u_{r_2}|u_{r_1})$','$p(u_{r_2}|u_{r_1},u_{r_0})$'},'Position',[0.541676093510218 0.77902142683665 0.257004955836705 0.113571428117298],'interpreter','latex')
fig_setup
ylabel('PDF','interpreter','latex')
txt = {'(e)'};
e=text(4,0.5,txt);
e.Units='normalized';


fig(6)=figure;
h = axes;
plot(x_input,P_AIB(:,find(y_input > 0.5-tol & y_input < 0.5+tol,1, 'first')),'Color','k','LineWidth',2) 
hold on
if ~isempty(find(y_input_2 > 0.5-tol & y_input_2 < 0.5+tol,1, 'first'))
plot(x_input_2,P_AIB_2(:,find(y_input_2 > 0.5-tol & y_input_2 < 0.5+tol,1, 'first')),'Color','r','LineWidth',2) 
end
% y_input(find(y_input > 0.5-tol & y_input < 0.5+tol,1,'first'))
% y_input_2(find(y_input_2 > 0.5-tol & y_input_2 < 0.5+tol,1,'first'))
view([90 90]);
set(h,'XDir','reverse');
if norm_ur==1  
    xlabel(['$u_{r_2=' num2str(tau2/markov) '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
else
    xlabel(['$u_{r_2=' num2str(tau2/markov) '\Delta_{EM}}$'], 'interpreter','latex')    
end
set(gca,'FontSize',18)
% xlim([-5 5])
axis square
set(gcf, 'Color', 'w')
set(gca,'yscale','log');
% ylim([10^-4 2])
if norm_ur==1  
    title(['$u_{r_1=' num2str(tau1/markov) '\Delta_{EM}}/ \sigma_\infty$ = ',num2str(y_input(cut_b),'%1.2f')],'interpreter','latex')
else
    title(['$u_{r_1=' num2str(tau1/markov) '\Delta_{EM}}$ = ',num2str(y_input(cut_b),'%1.2f')],'interpreter','latex')  
end
legend({'$p(u_{r_2}|u_{r_1})$','$p(u_{r_2}|u_{r_1},u_{r_0})$'},'Position',[0.541676093510218 0.77902142683665 0.257004955836705 0.113571428117298],'interpreter','latex')
fig_setup
ylabel('PDF','interpreter','latex')
txt = {'(f)'};
f=text(4,0.5,txt);
f.Units='normalized';






font=22;
a.FontSize = font;
b.FontSize = font;
c.FontSize = font;
d.FontSize = font;
e.FontSize = font;
f.FontSize = font;

pos_txt=[-0.22   0.9];
a.Position=pos_txt;
b.Position=pos_txt;
c.Position=pos_txt;
d.Position=pos_txt;
e.Position=pos_txt;
f.Position=pos_txt;

fig_setup
if ischar(save_path)
    savefig(fig,fullfile(save_path,append(save_name,'_','con_PDF.fig')),'compact')
    for a = 1:length(fig)
        exportgraphics(fig(a),fullfile(save_path,append(save_name,'_',sprintf('con_PDF_%d.png', a))))
    end
end




%Joint PDF
figure
surf(y_input,x_input,P_AnB)
hold on
mesh(y_input_2,x_input_2,P_AnB_2)
set(gca,'FontSize',18)
% set(gca,'zScale','log') 
axis([min(y_input) max(y_input) min(x_input) max(x_input)])
set(gcf, 'Color', 'w')
% legend('experiment','optimized','Location','NorthWest')
xlabel(['$u_{r_1=' num2str(tau1/markov) '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
ylabel(['$u_{r_2=' num2str(tau2/markov) '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
grid off
% set(gca,'zscale','log');
axis square
close

%conditional PDF
figure
surf(y_input,x_input,P_AIB)
hold on
% surf(y_input,x_input,real(P_AIB))
% mesh(y_input,x_input,real(P_moment))
%  set(gca,'zScale','log') 
%  axis([min(y_input) max(y_input) min(x_input) max(x_input) 0.001 10])
mesh(y_input_2,x_input_2,P_AIB_2)
set(gca,'FontSize',18)
axis([min(y_input) max(y_input) min(x_input) max(x_input)])
axis vis3d
set(gcf, 'Color', 'w')
% legend('experiment','optimized','Location','NorthWest')
% title(['directly=' num2str(fval_start,2) '; ' 'optimized=' num2str(fval,2)])
xlabel(['$u_{r_1=' num2str(tau1/markov) '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
ylabel(['$u_{r_2=' num2str(tau2/markov) '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
grid off
close

figure
surf(y_input,x_input,P_AIB)
hold on
% surf(y_input,x_input,real(P_AIB))
% mesh(y_input,x_input,real(P_moment))
set(gca,'zScale','log') 
axis([min(y_input) max(y_input) min(x_input) max(x_input) 0.01*min(level) 2*max(level)])
mesh(y_input_2,x_input_2,P_AIB_2)
set(gca,'FontSize',18)
% axis([min(y_input) max(y_input) min(x_input) max(x_input)])
axis vis3d
set(gcf, 'Color', 'w')
% legend('experiment','optimized','Location','NorthWest')
% title(['directly=' num2str(fval_start,2) '; ' 'optimized=' num2str(fval,2)])
xlabel(['$u_{r_1=' num2str(tau1/markov) '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
ylabel(['$u_{r_2=' num2str(tau2/markov) '\Delta_{EM}}/ \sigma_\infty$'], 'interpreter','latex')
grid off
close

for i=1:6
    tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
    uiwait(gcf);
    delete(tmp_ui);
    close
end
end