%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function tests whether the data have to be flipped or not. The decision of flipping of data 
% depends on a simple relation of $3^{rd}$ order structure function ($S^{3}$) with the dissipation 
% based on the assumption of homogeneous isotropic turbulence (HIT). The thumb rule is that the 
% quantity ($S^{3}$) must be negative. In the literature, the keyword that goes with this picture 
% is vortex stretching. To verify this, a plot of ($S^{3}$) as a function of the scale $r$ is 
% plotted, from which it is possible to decide whether it is essential to flip the data or not.
% The condition is S3 should always be negative! If we see S3 to be positve===> we flip the data
% Actual flipping will be done in next function called 'normalization'
%
% Arguments IN
% data = it is the filtered data 
% Fs = Acquisition/Sampling Frequency in Hz
% int_L = Integral length scale in meters
% taylor_L = Taylors length scale in meters
% save_path = path for saving figures and files
% save_name = name for saving files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Struc_flip_test(Fs,data,int_L,taylor_L,save_path,save_name)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');

m_data  =   mean(data);% mean of the data
clear r tau step S_exp_3
% r=logspace(-3,0,30);

r   = unique(round(logspace(log10(0.1*taylor_L*10^6),log10(2*int_L*10^6),40))./10^6,'stable');

tau(1:length(r))    = nan;
step(1:length(r))   = nan;

for k=1:length(r)
    tau(k)  = r(k)/m_data;
    step(k) = round(tau(k)*Fs);  
end

 stepp  = unique(step);
 stepp  = stepp(stepp>0);
 S_exp_3(1:length(stepp)) = nan; % Third order structure function
 Skew(1:length(stepp)) = nan;
 
 clear r
parfor k=1:length(stepp)
%     k/length(stepp)   
    tau         = stepp(k)
    r(k)        = (tau./Fs).*m_data;  
    incr        = (data(tau+1:length(data)) - data(1:length(data)-tau)).';    
    S_exp_3(k)  = mean(incr.^3); % Third order structure function
%     Skew(k)     = skewness(incr);
end

in_n_nan   = ~isnan(r);
r          = r(in_n_nan);
S_exp_3    = S_exp_3(in_n_nan);
% Skew       = Skew(in_n_nan); 

% Plotting r Vs. S_exp_3
% Remember the quantity 'S_exp_3' must be negative, if not you need to flip the data
h(1)=figure;
plot(r,S_exp_3,'MarkerSize',8,'Marker','o','LineWidth',1.5,'LineStyle','none','Color',[0 0 0])
% hold on
% plot(r,Skew,'MarkerSize',8,'Marker','o','LineWidth',1.5,'LineStyle','none')
set(gca,'xScale','log')
%  set(gca,'yScale','log') 
axis square
set(gcf, 'Color', 'w')
set(gca,'FontSize',24)
xlabel('$r/m$', 'interpreter','latex')
ylabel('$S^3(r)$', 'interpreter','latex')
xlabel('$r\ (m)$','interpreter','latex');
ylabel('$S^3$', 'interpreter','latex')
xlim([min(r)*0.5 max(r)*2])
%ylim([2e-3 0.2])
vline(taylor_L,'k')
vline(int_L,'k')
fig_setup
set(gca,'XTick',[10^-3 10^-2  10^-1  10^0]);
legend off

if ischar(save_path)
    savefig(h,fullfile(save_path,append(save_name,'_','struc_flip_test.fig')),'compact')
end
end