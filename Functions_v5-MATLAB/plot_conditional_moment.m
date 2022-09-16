%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots the first and second conditional moments as a function of the scale 
% separation. For this purpose, a scale and the number of a bin (value of the increment) condition 
% must be specified. In addition, a linear extrapolation in $\Delta r$ (solid black line) of the 
% first and second order conditional moments is plotted (see Chapter: Estimation of Kramers-Moyal 
% coefficients).
%
% Arguments IN
% scal = Plot conditional moment number of Scale
% bin_num = Number of a bin (value of the velocity increment u_r)
% evaluated = struct array containing all the information about conditional moments for each scale and for each bin
% step_con_moment = the steps at which the value of conditional moments are calculated
% markov = markov length in number of samples
% multi_point = Multipoint condition 1=YES or 2=NO
% increment_bin = number of bins
% condition = This input is for multipoint statistics
% tol = This input is for multipoint statistics
% data_filter = filtered data 
% save_path = path for saving figures and files
% save_name = name for saving files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_conditional_moment(scal,bin_num,evaluated,step_con_moment,markov,multi_point,increment_bin,condition,tol,data,save_path,save_name)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');

if multi_point==1
    k=askInput({sprintf('Plot: Multi-point number of point condition:')}, {num2str(ceil(length(condition)/2))});
    my_field = sprintf('point_eval_%d', k);  
    evaluated=evaluated.(my_field);
end

if (evaluated(scal).r_samp-1)<3*markov
    step_con_moment_tmp     = unique(round(logspace(0,log10(evaluated(scal).r_samp-1),12)),'stable');
else  
    step_con_moment_tmp     = unique(round(logspace(0,log10(3*markov),12)),'stable');
% step_con_moment_tmp     = round(linspace(0,2*markov,5));
% step_con_moment_tmp      = round(linspace(markov,2*markov,5));
end

tau2                = evaluated(scal).r_samp;
bin_events          = NaN(increment_bin,length(step_con_moment_tmp));
M11                 = nan; 
M21                 = nan;
M31                 = nan;
M41                 = nan;
eM1                 = nan;
eM2                 = nan;

for jj=length(step_con_moment_tmp):-1:1
    if multi_point==1
        [incr1,incr2]  = Increment_point((tau2 - step_con_moment_tmp(jj)),tau2,data,condition(k),tol);
    elseif multi_point==0
        [incr1,incr2]  = Increment((tau2 - step_con_moment_tmp(jj)),tau2,data);
    end

    [dx,dy,x,y,dA] = limit(incr1,incr2,increment_bin);
    if length(x)==increment_bin & length(y)==increment_bin
       %     [M11_direct(:,jj),M21_direct(:,jj),M41_direct(:,jj),x_mean_bin_direct,y_mean_bin_direct,events] = conditional_moments_direct(tau1,tau2,incr1,incr2,aufl,min_events);
       %     bin_events_direct(:,jj)=events.';
        [P_AIB,P_BIA,P_AnB,P_A,P_B,binP_AIB,evaluated(scal).x_mean_bin,evaluated(scal).y_mean_bin,events,counter_A,counter_B] = distribution((tau2 - step_con_moment_tmp(jj)),tau2,incr1,incr2,increment_bin,dx,dy,x,y,dA);
        bin_events(:,jj)=events.';
        %----------------------------- conditional moments --------------------------------
        for i=1:increment_bin
            if (tau2 - step_con_moment_tmp(jj))>tau2
                   M11(i,jj)    = nansum(binP_AIB(:,i).^1.*P_AIB(:,i).*dy);
                   M21(i,jj)    = nansum(binP_AIB(:,i).^2.*P_AIB(:,i).*dy);
                   M31(i,jj)    = nansum(binP_AIB(:,i).^3.*P_AIB(:,i).*dy);
                   M41(i,jj)    = nansum(binP_AIB(:,i).^4.*P_AIB(:,i).*dy);
            else
                   M11(i,jj)    = nansum(binP_AIB(:,i).^1.*P_AIB(:,i).*dx);
                   M21(i,jj)    = nansum(binP_AIB(:,i).^2.*P_AIB(:,i).*dx);
                   M31(i,jj)    = nansum(binP_AIB(:,i).^3.*P_AIB(:,i).*dx);
                   M41(i,jj)    = nansum(binP_AIB(:,i).^4.*P_AIB(:,i).*dx);
            end
        end
        %Fehler von den Momenten
        eM1= real(sqrt(abs((M21 - (M11.^2))./bin_events)));
        eM2= real(sqrt(abs((M41 - (M21.^2))./bin_events)));
    end
end

h(1) = figure;
plot(step_con_moment,evaluated(scal).M11(bin_num,:),'Color',[0.00,0.45,0.74],'MarkerSize',8,'Marker','o', 'LineStyle','none','LineWidth',2)
hold
plot(step_con_moment,evaluated(scal).M21(bin_num,:),'Color',[0.85,0.33,0.10],'MarkerSize',8,'Marker','+', 'LineStyle','none','LineWidth',2)

plot(step_con_moment_tmp,M11(bin_num,:),'k','MarkerSize',8,'Marker','o', 'LineStyle','none','LineWidth',1.5)
plot(step_con_moment_tmp,M21(bin_num,:),'k','MarkerSize',8,'Marker','+', 'LineStyle','none','LineWidth',1.5)

xlabel('$\Delta r/sample$','Interpreter','latex')
% ylabel('Conditional Moment')
set(gca,'FontSize',18)
set(gcf, 'Color', 'w')
axis square
fig_setup
xlim([min(step_con_moment_tmp)-2 max(step_con_moment_tmp)+2])      
ylim([min([M11(bin_num,:) M21(bin_num,:)]) max([M11(bin_num,:) M21(bin_num,:)])])      


[xData, yData] = prepareCurveData( step_con_moment,evaluated(scal).M11(bin_num,:));
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
[fitresult, gof]        = fit( xData, yData, ft, opts );
plot(linspace(0,2*max(step_con_moment_tmp),100),feval(fitresult,linspace(0,2*max(step_con_moment_tmp),100)),'k','LineWidth',1.5,'LineStyle','-')

[xData, yData] = prepareCurveData( step_con_moment,evaluated(scal).M21(bin_num,:));
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
[fitresult, gof]        = fit( xData, yData, ft, opts );
plot(linspace(0,2*max(step_con_moment_tmp),100),feval(fitresult,linspace(0,2*max(step_con_moment_tmp),100)),'k','LineWidth',1.5,'LineStyle','-')

vline(min(step_con_moment),'k')
vline(max(step_con_moment),'k')

legend('$M^{(1)}\left(u_r,r\right)$ ','$M^{(2)}\left(u_r,r\right)$ ','Location','northwest')

if ischar(save_path)
    savefig(h,fullfile(save_path,append(save_name,'_','con_moment.fig')),'compact')
    for a = 1:length(h)
        exportgraphics(h(a),fullfile(save_path,append(save_name,'_',sprintf('con_moment_%d.png', a))))
    end
end

set(gca,'xScale','lin') 
set(gca,'yScale','lin') 
tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf);
delete(tmp_ui);
close
end