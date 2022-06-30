%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs the pointwise optimization of Kramers-Moyal coefficients 
% $D^{(1,2)}\left(u_r,r\right)$ at each scale and value of velocity increment to minimize possible 
% uncertainties in the absolute values of the Kramers-Moyal coefficients. The object of this 
% optimization is to find the best Fokker-Planck equation to reproduce the conditional PDF's as 
% these are the essential part of the Markov process. This optimization procedure is proposed in 
% \cite{kleinhans2005iterative,Nawroth2007,Reinke2018} and it includes the reconstruction of the 
% conditional probability density functions $p\left(u_{r'} | u_r \right)$ via the short time 
% propagator \cite{Risken}. The optimization procedure systematically changes 
% $D^{(1,2)}\left(u_r,r\right)$ until the error function is minimized.
%
% Arguments IN
% evaluated = struct array calculated in the function 'conditional_moment'
% scal = A scale number which needs to be optimized
% increment_bin = number of bins
% markov = markov length in number of samples
% data = filtered data
% Fs = Acquisition/Sampling Frequency in Hz
% tol = Optimization: Tolerance of the range of Kramers-Moyal coefficients in % 
%       It is the percentage(For Ex: 0.1 for 10 percent or 0.2 for 20 percent)of D1 or D2 
%       within these limit which you want to optimize these coeffcients D1 & D2
% m_data = mean of the data
% taylor_L = Taylor length scale in meters
% test_opti = weather to plot or not ==> 'Plot? 1=Yes, 0=No'
% multi_point = weather to do multi-point analysis or not 1=Yes, 0=No
% condition = condition for multi-point analysis
% tol_point = This input is for multipoint statistics
% min_events = minimum number of events
% k = Multi-point number of point condition
%
% Arguments OUT
% evaluated = a modified/updated struct 'evaluated' array after the optimization at each scale and at each bin for D1 & D2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [evaluated]=pointwise_KM_STP_optimization(evaluated,scal,increment_bin,markov,data,Fs,tol,m_data,taylor_L,test_opti,multi_point,condition,tol_point,min_events,k,norm_ur,norm_r,save_path,save_name) 
% Initial D1 and D2 before optimization 
clear x0
x0(1,1:2*sum(evaluated(scal).x_bin_not_nan))=[(evaluated(scal).D1(evaluated(scal).x_bin_not_nan)),(abs(evaluated(scal).D2(evaluated(scal).x_bin_not_nan)))];       

tau1            =   evaluated(scal).r_short_sample;
tau2            =   evaluated(scal).r_samp;
x_bin_not_nan   =   evaluated(scal).x_bin_not_nan;
D1              =   evaluated(scal).D1;
D2              =   evaluated(scal).D2;
eD1             =   evaluated(scal).eD1;
eD2             =   evaluated(scal).eD2;
x_mean_bin      =   evaluated(scal).x_mean_bin;
y_mean_bin      =   evaluated(scal).y_mean_bin;


if multi_point==1
    [incr1,incr2]  = Increment_point(tau1,tau2,data,condition(k),tol_point);
elseif multi_point==0
    [incr1,incr2]  = Increment(tau1,tau2,data);
end
          
[dx,dy,x,y,dA] = limit(incr1,incr2,increment_bin);
% [dx,dy,x,y,dA] = limit_STP(incr1,incr2,increment_bin);


%% calculation of PDF
[P_AIB,P_BIA,P_AnB,P_A,P_B,binP_AIB,~,~,~,counter_A,counter_B] = distribution(tau1,tau2,incr1,incr2,increment_bin,dx,dy,x,y,dA,1);
counter_B(isnan(D1)) = nan;


%% Tolerance of the range of Kramers-Moyal coefficients
% constant value
% ub       = [(D1(x_bin_not_nan))+tol,abs(D2(x_bin_not_nan))+tol]; %Upper bound
% lb       = [(D1(x_bin_not_nan))-tol,zeros(1,sum(x_bin_not_nan))]; %Lower bound
% percent of the range 
ub       = [(D1(x_bin_not_nan))+tol*(range(D1(x_bin_not_nan))),abs(D2(x_bin_not_nan))+tol*(range(D2(x_bin_not_nan)))]; %Upper bound
lb        = [(D1(x_bin_not_nan))-tol*(range(D1(x_bin_not_nan))),zeros(1,sum(x_bin_not_nan))];%Lower bound
% percent of the value itself
% tol=0.1; %in Prozent vom eigentlichen Wert
% ub       = [(D1(x_bin_not_nan))+(tol.*abs(D1(x_bin_not_nan))),abs(D2(x_bin_not_nan))+(tol.*abs(D2(x_bin_not_nan)))]; %Upper bound
% lb       = [(D1(x_bin_not_nan))-(tol.*abs(D1(x_bin_not_nan))),zeros(1,length(x_bin_not_nan))]; %Lower bound


%% Optimization statement
% >>opt=optimset(’TolFun’,1.0e-10) % Cancel, if function value difference
% less than or equal to 1.0e-10
x1(1:length(x0)) = nan;
% starting values without optimization 
fval_start  = calc_divergence(x0,y_mean_bin,x_mean_bin,P_B,P_AnB,P_AIB,markov,x_bin_not_nan,counter_A,counter_B,min_events,taylor_L,m_data,Fs,tau1,tau2,norm_ur,norm_r); 
% print details of optimzation steps 
% options     = optimset('LargeScale','off','GradObj','off','Hessian','off','TolX',1e-14,'TolFun',1e-14,'Display','iter-detailed',...
%              'MaxIter',400, 'MaxFunEvals',1e7,'FunValCheck','on','Algorithm','interior-point');
% no steps details
options     = optimset('LargeScale','off','GradObj','off','Hessian','off','TolX',1e-14,'TolFun',1e-14,...
             'MaxIter',400, 'MaxFunEvals',1e7,'FunValCheck','on','Algorithm','interior-point');
[x1,fval,exitflag,output] = ...
fmincon(@(x0)calc_divergence(x0,y_mean_bin,x_mean_bin,P_B,P_AnB,P_AIB,markov,x_bin_not_nan,counter_A,counter_B,min_events,taylor_L,m_data,Fs,tau1,tau2,norm_ur,norm_r),x0,[],[],[],[],lb,ub,[],options);
disp(fval_start)
disp(fval)

% figure
% plot(ub)
% hold on
% plot(lb)
% plot(x0,'k')
% plot(x1)

evaluated(scal).start_opti      = fval_start;
evaluated(scal).end_opti        = fval;
evaluated(scal).x0              = x0;
evaluated(scal).x1              = x1;
evaluated(scal).D1_opti         = x1(1:sum(x_bin_not_nan));
evaluated(scal).D2_opti         = x1(sum(x_bin_not_nan)+1:2*sum(x_bin_not_nan));
               
if test_opti==1  
    plot_pointwise_KM(y,x,y_mean_bin,x_mean_bin,x_bin_not_nan,x0,x1,P_B,P_AnB,P_AIB,counter_A,counter_B,min_events,eD1,eD2,markov,fval_start,fval,tau1,tau2,m_data,Fs,taylor_L,tol,dy,increment_bin,ub,lb,norm_ur,norm_r,save_path,save_name);
end
end