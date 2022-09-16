%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs the pointwise optimization of Kramers-Moyal coefficients towards the 
% integral fluctuation theorem will be done. Thereby the separation of scales/step increment 
% (in samples) referred to the sequence from large to small scales in the cascade trajectory is set 
% to a minimum of 1 sample and for the optimization, independent cascade trajectories (z=3) are used. 
% Note, we use here a separation that is less than or equal to the Einstein-Markov length of $\Delta_{EM}$.
% 
% This function performs the optimization of $D^{(1,2)}\left(u_r,r\right)$ at each scale and at each 
% value of velocity increment in order to satisfy the integral fluctuation theorem with minimum 
% possible error and plots the optimized $d_{ij}$ as a function of $r$. The optimization procedure 
% systematically changes $D^{(1,2)}\left(u_r,r\right)$ until the error function 
% \begin{eqnarray}
% \xi =|1-\langle e^{\mathrm{-}\Delta S_{tot}} \rangle_{max(N)}|
% \end{eqnarray}
% is minimized. 
% 
% Arguments IN
% tol_D1          = Tolerance of the range of Kramers-Moyal coefficients in %% for the Optimization
% tol_D2          = tol_D1;
% evaluated = a modified/updated struct 'evaluated' array with in function 'KM_STP_optimization'
% co_KM_opti = Coefficients $d_{ij}(r)$ of the optimized Kramers-Moyal coefficients towards the best 
%              Fokker-Planck equation to reproduce the conditional PDF's using the surface fits
% fitresult_D1_conf = Confidence intervals for fit coefficients of D1 (co_KM_opti)
% fitresult_D2_conf = Confidence intervals for fit coefficients of D2 (co_KM_opti)
% int_L = Integral length scale in meters
% taylor_L = Taylor length scale in meters
% m_data = mean of the data
% Fs = Acquisition/Sampling Frequency in Hz
% data = filtered data 
% z = This is the parameter which decides how to compute entropy using different methods using either 
% overlapping or independent trajectories
% z = 1 ==> Overlapping trajectories 
% z = 3 ==> Independent trajectories
% iter = The maximum number of iteration used for the optimization algorithms
% markov = markov length in number of samples
% trajec = 1 ==> The start/end of the cascade trajectory will be adjusted.
% dr_ind = Separation of scales/step increment (in samples) referred to the sequence from large to small scales in the cascade trajectory
% norm_ur = normalization of the data using $\sigma_\infty$? data 1=Yes, 0=No
% norm_r = normalization of the scale using $\lambda$? data 1=Yes, 0=No
% save_path = path for saving figures and files
% save_name = name for saving files
% 
% Arguments OUT
% co_IFT_opti = Coefficients $d_{ij}(r)$ of the optimized Kramers-Moyal coefficients towards the 
%               integral fluctuation theorem using the surface fits
% history = details of the optimization for every single iteration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [co_IFT_opti,history]=OPTI_IFT_dij(tol_D1,tol_D2,evaluated,co_KM_opti,fitresult_D1_conf,fitresult_D2_conf,int_L,taylor_L,m_data,Fs,data,z,iter,markov,trajec,dr_ind,norm_ur,norm_r,save_path,save_name)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');

% if z==1
%     [Sm,Ds,DS,~,~,~,~]=checkFT(evaluated,int_L,taylor_L,Fs,data,m_data,z,co_KM_opti,1,markov,trajec,dr_ind);
%     %Plot of Entropy
%     plot_entropy(Sm,Ds,DS)
%     set(gca, 'FontSize',10)
%     data_length=str2double(inputdlg({'Data length for optimization in samples:'}));
%     clear Sm Ds DS
% end


%% Initial D1 and D2 before optimization  
clear r_norm
r_norm           = nan(1,length(evaluated));

if norm_r==1
    for i=1:length(evaluated)
        r_norm(1,i)  = (evaluated(i).r./taylor_L);
    end
else
    for i=1:length(evaluated)
        r_norm(1,i)  = (evaluated(i).r);
    end
end


% fitresult_D1
d11 = co_KM_opti.a(1).*r_norm.^co_KM_opti.ea(1)+co_KM_opti.a(2);
% fitresult_D2
d20 = co_KM_opti.b(1).*r_norm.^co_KM_opti.eb(1)+co_KM_opti.b(2);
d21 = co_KM_opti.b(3).*r_norm.^co_KM_opti.eb(2)+co_KM_opti.b(4);
% d21=0.*d21;
d22 = co_KM_opti.b(5).*r_norm.^co_KM_opti.eb(3)+co_KM_opti.b(6);
x0  = [d11 d20 d21 d22];


%% Tolerance of the range of Kramers-Moyal coefficients
clear x1 lb ub
% x1(1:length(x0))    = nan;
lb(1:length(x0))    = nan;
ub(1:length(x0))    = nan;

% percent of the range 
% D1 Opti
t                       = length(r_norm);
lb(1:t)                 = x0(1:t)-abs(x0(1:t).*tol_D1); % Lower bound
ub(1:t)                 = x0(1:t)+abs(x0(1:t).*tol_D1); % Upper bound

% D2 Opti
for i=1:3
    lb(t*i+1:(t*i+1)+t-1)   = x0(t*i+1:(t*i+1)+t-1)-abs(x0(t*i+1:(t*i+1)+t-1).*tol_D2);
    ub(t*i+1:(t*i+1)+t-1)   = x0(t*i+1:(t*i+1)+t-1)+abs(x0(t*i+1:(t*i+1)+t-1).*tol_D2);
end


% constant value 10% of max/min value
% D1 Opti
% t                       = length(r_norm);
% lb(1:t)                 = x0(1:t)-(max(abs(x0(1:t))).*tol_D1);
% ub(1:t)                 = x0(1:t)+(max(abs(x0(1:t))).*tol_D1);
 
% D2 Opti
% for i=1:3
%     lb(t*i+1:(t*i+1)+t-1)   = x0(t*i+1:(t*i+1)+t-1)-(max(abs(x0(t*i+1:(t*i+1)+t-1))).*tol_D2);
%     ub(t*i+1:(t*i+1)+t-1)   = x0(t*i+1:(t*i+1)+t-1)+(max(abs(x0(t*i+1:(t*i+1)+t-1))).*tol_D2);
% end


% fix to min and max value
% d11
% lb(1:t)                 = min(lb(1:t));
% ub(1:t)                 = max(ub(1:t));

% d20
% lb(t*1+1:(t*1+1)+t-1)   = min(lb(t*1+1:(t*1+1)+t-1));
% ub(t*1+1:(t*1+1)+t-1)   = max(ub(t*1+1:(t*1+1)+t-1));
 
% d21
% lb(t*2+1:(t*2+1)+t-1)   = min(lb(t*2+1:(t*2+1)+t-1));
% ub(t*2+1:(t*2+1)+t-1)   = max(ub(t*2+1:(t*2+1)+t-1));

% d22
% lb(t*3+1:(t*3+1)+t-1)   = min(lb(t*3+1:(t*3+1)+t-1));
% ub(t*3+1:(t*3+1)+t-1)   = max(ub(t*3+1:(t*3+1)+t-1));


% using confidence intervals for fit coefficients of D1 and D2 (co_KM_opti)
% D1 Opti
% lb(1:2)    = [fitresult_D1_conf(1,1:2)];
% lb(9)      = [fitresult_D1_conf(1,3)];
% ub(1:2)    = [fitresult_D1_conf(2,1:2)];
% ub(9)      = [fitresult_D1_conf(2,3)];
 
% D2 Opti
% lb(3:8)       = [fitresult_D2_conf(1,1:6)];
% lb(10:end)    = [fitresult_D2_conf(1,7:9)];
% ub(3:8)       = [fitresult_D2_conf(2,1:6)];
% ub(10:end)    = [fitresult_D2_conf(2,7:9)];

% lb(isnan(lb))=0;
% ub(isnan(ub))=0;


%% The constraints were set in a physically and mathematically meaningful way: 
% (see for more details \textbf{\textit{FIT\_KM}})
% d11 <= 0
if ub(1:t)>0
    ub(ub(1:t)>0) = 0;
end
% d20 >= 0 
if lb(t*1+1:(t*1+1)+t-1)<0
    lb(lb(t*1+1:(t*1+1)+t-1)<0) = 0;
end
% d22 >= 0 
if lb(t*3+1:(t*3+1)+t-1)<0
    lb(lb(t*3+1:(t*3+1)+t-1)) = 0;
end



%% Selcet dij for the optimization
if askYesno('Do you want to perform the pointwise optimization of d_11?', 'Yes')
else
    lb(1:t)                 = x0(1:t); % Lower bound
    ub(1:t)                 = x0(1:t); % Upper bound
end

if askYesno('Do you want to perform the pointwise optimization of d_20?', 'Yes')
else
    lb(t*1+1:(t*1+1)+t-1)   = x0(t*1+1:(t*1+1)+t-1);
    ub(t*1+1:(t*1+1)+t-1)   = x0(t*1+1:(t*1+1)+t-1);
end

if askYesno('Do you want to perform the pointwise optimization of d_21?', 'Yes')
else
    lb(t*2+1:(t*2+1)+t-1)   = x0(t*2+1:(t*2+1)+t-1);
    ub(t*2+1:(t*2+1)+t-1)   = x0(t*2+1:(t*2+1)+t-1);
end

if askYesno('Do you want to perform the pointwise optimization of d_22?', 'Yes')
else
    lb(t*3+1:(t*3+1)+t-1)   = x0(t*3+1:(t*3+1)+t-1);
    ub(t*3+1:(t*3+1)+t-1)   = x0(t*3+1:(t*3+1)+t-1);
end


%% Plot dij(r) and optimizaiton range
figure
subplot(2,2,1);
plot(r_norm,...
    lb(1:t),'LineWidth',2,'LineStyle','--')
hold on
plot(r_norm,...
    d11,'k','LineWidth',2)
plot(r_norm,...
    ub(1:t),'LineWidth',2,'LineStyle','--')
axis square
ylabel('$d_{11}$', 'interpreter','latex')
if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
    else
    xlabel('$r/m $', 'interpreter','latex')
end
xlim([min(r_norm) max(r_norm)])
set(gca, 'FontSize',24)
%
subplot(2,2,2);
plot(r_norm,...
    lb(t*1+1:(t*1+1)+t-1),'LineWidth',2,'LineStyle','--')
hold on
plot(r_norm,...
    d20,'k','LineWidth',2)
plot(r_norm,...
    ub(t*1+1:(t*1+1)+t-1),'LineWidth',2,'LineStyle','--')
axis square
ylabel('$d_{d20}$', 'interpreter','latex')
if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
    else
    xlabel('$r/m $', 'interpreter','latex')
end
xlim([min(r_norm) max(r_norm)])
set(gca, 'FontSize',24)
%
subplot(2,2,3);
plot(r_norm,...
    lb(t*2+1:(t*2+1)+t-1),'LineWidth',2,'LineStyle','--')
hold on
plot(r_norm,...
    d21,'k','LineWidth',2)
plot(r_norm,...
    ub(t*2+1:(t*2+1)+t-1),'LineWidth',2,'LineStyle','--')
axis square
ylabel('$d_{d21}$', 'interpreter','latex')
if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
    else
    xlabel('$r/m $', 'interpreter','latex')
end
xlim([min(r_norm) max(r_norm)])
set(gca, 'FontSize',24)
%
subplot(2,2,4);
plot(r_norm,...
    lb(t*3+1:(t*3+1)+t-1),'LineWidth',2,'LineStyle','--')
hold on
plot(r_norm,...
    d22,'k','LineWidth',2)
plot(r_norm,...
    ub(t*3+1:(t*3+1)+t-1),'LineWidth',2,'LineStyle','--')
axis square
ylabel('$d_{d22}$', 'interpreter','latex')
if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
    else
    xlabel('$r/m $', 'interpreter','latex')
end
xlim([min(r_norm) max(r_norm)])
set(gca, 'FontSize',24)
set(gcf, 'Color', 'w')


%% Fix length of trajectories
if trajec==1
        L       = 0;
        % Generally each Overlapping/Non-Overlapping trajectory starts at L and
        % ends at Lambda for the calculation of the entropy; but its a choice
        % So we ask you at which biggest scale you want to start the trajectory
        tmp     = str2double(inputdlg({'Adjustment: Start of trajectory 1=Yes, 0=No'}));
        if tmp == 1
            % Here, if you want to start your trajectories at scale smaller than L then
            % you should enter a value between 0 & 1 
            % For ex: if your L=1 meters and you enter 0.9 your trajectory will start at 0.9 meters
            L   = str2double(inputdlg({'Start of trajectory in multiples of Integral length scale:'}))*int_L/taylor_L;
        else
            L   = int_L/taylor_L; % dimless integral scale
        end
        % So we ask you at which smallest scale you want to end the trajectory
        tmp     = str2double(inputdlg({'Adjustment: End of trajectory 1=Yes, 0=No'}));
        if tmp == 1
            % Here, if you want to end your trajectories at scale bigger than lambda then
            % you should enter a value between 1 and (L/lambda) 
            % For ex: if your lambda=1 cm and if you enter 5 your
            % trajectory will end at 5 cms
            l   = str2double(inputdlg({'End of trajectory in multiples of Taylor length scale:'}))*taylor_L/taylor_L;
        else
            l   = evaluated(1).r/taylor_L;% dimless smallest scale by D1 and D2
        end
else % start=L, end=lambda
    L   = int_L/taylor_L;
    l   = 1;
end

x       = m_data.*(1:ceil(1.1*L*Fs*taylor_L/m_data)).'./Fs;
x       = x/taylor_L; % dimless position
Lind    = find(abs(x(:,1)-L  )==min(abs(x(:,1)-L  )), 1, 'last' ); % index for x showing the LAST  entry where x is within range lambda..L
lind    = find(abs(x(:,1)-l  )==min(abs(x(:,1)-l  )), 1, 'first'); % index for x showing the FIRST entry where x is within range lambda..L

% Index for calculating the increment using data
Lind    = Lind+1;   
lind    = lind+1;

if lind<2
    lind = 2;
end

rindful = 1:Lind;                       % index-vector for all scales
rind    = fliplr(rindful(lind:Lind));   % index-vector only for scales from start to end of the cascade trajectory
% rind = round(linspace(Lind,lind,floor(L/l)+1)); % index-vector for scales in inertial range in steps of lambda
stp     = Lind;                         % index used to cut v(x) in non-overlapping pieces of length equal to integral scale
% R(k).r   = x(rind); % scales from index-vector 
if norm_r==1
    r       = x(rind); % scales from index-vector 
else
    r       = x(rind).*taylor_L; % scales from index-vector 
end
clear x
clear Sm Ds DS tmp_Sm tmp_Ds tmp_DS not_nan_index u x



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calc Entropy - overlapping trajectories
% if z==1
%     %% PDF of start/end of the cascade trajectory ==> Entropy of the system or the Shanon entropy
%     rind_tmp    = rind(1:dr_ind:end);
%     end_L       = size(data,1);
%     tmp_size    = size(data,1)-stp+1;
%     u_tmp       = nan( tmp_size, 2 );
% 
%     for i=0:end_L-stp
%         %% u
%         % u_tmp(i+1,:) = (data([rind_tmp(1)+i,rind_tmp(end)+i])- data(1+i)).';
%         %% u_s
%         u_tmp(i+1,:) = (data([rind_tmp(2)+i,rind_tmp(end-1)+i])- data(1+i)).';
%     end   
% 
%     nbins=301;
%     [fuL, uL]       = hist(u_tmp(:,1),nbins); 
%     puL             = fuL / ( nansum(fuL) * mean(diff(uL)) ); % relative frequency to approximate p(u,r=L)
% 
%     [ful, ul]       = hist(u_tmp(:,2),uL); 
%     pul             = ful / ( nansum(ful) * mean(diff(ul)) ); % relative frequency to approximate p(u,r=lambda)
%     clear u_tmp
% 
%     [history] = runfmincon_overlap(x0,lb,ub,data(1:data_length),stp,r, uL, puL, ul, pul,rind,iter,r_norm,t,co_KM_opti);
%     [history, x1] = runfmincon_overlap(x0,lb,ub,data(1:10^4),stp,r, uL, puL, ul, pul,rind,iter);
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calc Entropy - non-overlapping trajectories
if z==3
    tmp_size    = ceil(numel(data)/stp);
    [u]         = create_us(data,r,rind,stp);
%     [u]       = create_us(data(1:0.5*10^7),r,rind,stp);
    [history]   = runfmincon(x0,lb,ub,u,r,tmp_size,iter,r_norm,t,co_KM_opti,markov,taylor_L,Fs,m_data,trajec,dr_ind,z);
end


% Find minimized  error function 
index                       = find(history.fval==min(history.fval),1,'first');
[fit_d11]                   = FIT_d1j(r_norm,history.x1_iter(index,1:t));
[fit_d20, fit_d21, fit_d22] = FIT_d2j(r_norm,history.x1_iter(index,:),t);


%% IFT optimzed dij coefficients
co_IFT_opti.a  = [fit_d11.a1 fit_d11.a11];
co_IFT_opti.ea = [fit_d11.ra1];
co_IFT_opti.b  = [fit_d20.b0 fit_d20.b00 fit_d21.b1 fit_d21.b11 fit_d22.b2 fit_d22.b22];
co_IFT_opti.eb = [fit_d20.rb0 fit_d21.rb1 fit_d22.rb2]; 


%% Plot: comparision co_KM_opti and co_IFT_opti dij coefficients
OPTI_IFT_plot(data,stp,r,rind,co_KM_opti,co_IFT_opti,evaluated,Fs,taylor_L,lb,ub,t,history,index,d11,d20,d21,d22,m_data,markov,trajec,dr_ind,z,norm_ur,norm_r,save_path,save_name)
end