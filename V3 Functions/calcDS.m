%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the total entropy variation DS for trajectories in u_in along scale r_in.
%
% Arguments IN
% z = This is the parameter which decides how to compute entropy using different methods using either 
% overlapping or independent trajectories
% z = 1 ==> Overlapping trajectories 
% z = 3 ==> Independent trajectories
% u_in = cascade trajectories for which the total entropy shall be calculated 
% r_in = Normalized scale vector from the start to end of the cascade trajectory
% co = Coefficients $d_{ij}(r)$ of the optimized Kramers-Moyal coefficients, where first entry is zeroth power 
% markov = markov length in number of samples
% taylor_L = Taylor length scale in meters
% Fs = Acquisition/Sampling Frequency in Hz
% m_data = mean of the data
% dr_ind = Separation of scales/step increment (in samples) referred to the sequence from large to small scales in the cascade trajectory
% 
% Arguments OUT
% Sm = Entropy of the medium
% Ds = Entropy of the system or the Shanon entropy
% DS = Total entropy production=Sm+Ds
% dr = separation of scales/step increment (in normalized scale) referred to the sequence from 
%      large to small scales in the cascade trajectory
% r_s_tmp = mid-point scale vector (Stratonovich convention)
% The cascade process takes some initial probability s_L = $p(u_L, L)$ and changes to a different 
% final probability s_l = $p(u_\lambda, \lambda)$.
% A = action functional, pathintegral of Lagrangian
% Lag,p,H,tmpvar,H1,H2,ur_s_tmp (are not yet included in the current version)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sm,Ds,DS,dr,r_s_tmp,s_l,s_L,Lag,A,p,H,tmp,H1,H2,ur_s_tmp] = calcDS(z,u_in,r_in,co,markov,taylor_L,Fs,m_data,dr_ind,varargin)
if z==1
    if nargin >= 10
        uL          = varargin{1};
        puL         = varargin{2};
        ul          = varargin{3};
        pul         = varargin{4};
    end
end


%% initialise scale vector and trajectory matrix
r           = r_in(1:dr_ind:end);   % get r-vector in steps of dr, end is cut off
u           = u_in(:,1:dr_ind:end); % same for u
dr          = (r(:,3:end) - r(:,1:end-2))./2;
dr          = -repmat(dr,[length(u(:,1)) 1]); %same dimension as u
%dr          = dr(2:end-1);     
%dr          = mean(abs(diff(r_in)));
%dr          = (dr_ind/Fs*m_data)/taylor_L;     
%dr          = abs(diff(r_in(1:dr_ind:end)));   
%dr          = (r_in(:,3:end) - r_in(:,1:end-2));


%% prepare trajectories/definition of the discretisation rule/interpretation of the SDE
% Due to the stochastic variable u, the integrals are stochastic integrals for which the usual 
% Riemann definition of an integral is not suitable. The consequence is that, in the continuous 
% limit, the value of the stochastic integral is not unique, depending on whether we take the value
% of the integrand at the beginning (Ito: pre-point rule) of each discretisation interval,
% at the end (end-point rule) or to somewhere in-betwee (for example Statonovich: the mid-point rule). 


%% pre-point rule (Ito)
% u_s     = u(:,1:end-1);

%% mid-point rule (Stratonovich convention, usual calculus rules, symmetric wrt scale-reversal, is identical with the trapezoidal rule)
% Stratonovich stochastic integral takes the average of beginning and end of the discretisation intervals
% u_s     = (u(:,1:end-1) + u(:,2:end))/2;    

%% mid-point for 2*dr
% u_s     = (u(:,3:end) + u(:,1:end-2))/2;                    
u_s     = u(:,2:end-1);   

%% r-derivative of u(r)
% Forward difference quotient 
% ur_s    = (u(:,2:end) - u(:,1:end-1)) ./ dr;                
% ur_s    = reshape(smooth(ur_s,ceil(0.25*size(ur_s,2))),size(ur_s,1),size(ur_s,2));

% Central difference quotient
ur_s    = (u(:,3:end) - u(:,1:end-2))./ (2.*dr);  
% u_central = [u u_in(:,end)];
% ur_s    = (u_central(:,3:end) - u_central(:,1:end-2)) ./ (2*dr);          

%% pre-point rule (Ito): r-values at beginning of intervals, repeat as matrix to match u_s 
% r_s     = repmat(r(1:end-1),[length(u(:,1)) 1]);

%% mid-point 
% r_s   = repmat((r(1:end-1) + r(2:end))/2,[length(u(:,1)) 1]);
% r_s   = repmat((r(1:end-2) + r(3:end))/2,[length(u(:,1)) 1]);     

%% mid-point for 2*dr
r_s         = repmat(r(2:end-1),[length(u(:,1)) 1]);
ur_s_tmp    = ur_s;
r_s_tmp     = r_s;


%% PDF of start/end of the cascade trajectory ==> Entropy of the system or the Shanon entropy
if z==3
    if nargin < 10
        nbins       = 301;
        % [fuL, uL]    = hist(u(:,1),nbins); 
        [fuL, uL]    = hist(u_s(:,1),nbins); 
        puL          = fuL / ( nansum(fuL) * mean(diff(uL)) ); % relative frequency to approximate p(u,r=L)

        % [ful, ul]    = hist(u(:,end),uL); 
        [ful, ul]    = hist(u_s(:,end),uL); 
        pul          = ful / ( nansum(ful) * mean(diff(ul)) ); % relative frequency to approximate p(u,r=lambda)

% probability of start/end of the cascade trajectory
        % [fuL uL]    = hist(u(:,1),nbins);  
        % puL         = fuL / ( sum(fuL)); 
        % [ful ul]    = hist(u(:,end),uL);
        % % % % [ful ul]    = hist(u(:,end),nbins); % shows a significant offset
        % pul         = ful / ( sum(ful));

%         figure
%         plot(uL,puL,'LineWidth',2)
%         hold on
%         plot(ul,pul,'LineWidth',2)
%         % plot(ul,smooth(pul,floor(0.1*nbins)))
%         norm1 = pdf(fitdist(u_s(:,1),'Normal'),uL);
%         plot(uL,norm1,'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176])
%         norm1 = pdf(fitdist(u_s(:,end),'Normal'),uL);
%         plot(ul,norm1,'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176])
%         set(gca,'yScale','log') 
%         set(gcf, 'Color', 'w')
%         set(gca, 'FontSize',18)
%         xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
%         ylabel('$p(u_r / \sigma_\infty)$','interpreter','latex')
%         axis square
%         legend({'$u_L$','$u_\lambda$'},'Interpreter','latex')
%         xlim([min(uL) max(uL)])
%         ylim([0.0000001*max(puL) 2*max(pul)])
%         fig_setup 
%         tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
%         uiwait(gcf);
%         delete(tmp_ui);
%         close
        clear u_tmp   
    end
end


%% Entropy of the medium
[Sm_1,Sm_2,F,D,FD,Sm_tmp]     = calcSm(ur_s, u_s, r_s, dr, co); 
% plot_entropy_trajectory(r_s,u_s,ur_s,r_in,u_in,F,D,FD,Sm_tmp)

% Sm_1: rectangles to approximation, Sm_2: trapezoidal approximation
Sm              = Sm_1;


%% System entropy production: entropy difference between initial and final states
% [Ds,s_l,s_L]    = calcD_system_entropy(u(:,1),u(:,end),uL,puL,ul,pul);
[Ds,s_l,s_L]    = calcD_system_entropy(u_s(:,1),u_s(:,end),uL,puL,ul,pul);


%% Lagrangian, conjugate Varaible, Hamiltionian, Action = pathintegral of Lagrangian
[Lag,A,p,H,tmp,H1,H2] = calc_Path_integral(ur_s,u_s,r_s,dr,co,z);
% plotFreddy(Sm,Lag,p,H,u_s,r_s,co)

clear r_in u_in ur_s u_s r_s
clear r u uL puL ul pul co
DS              = Sm + Ds;
DS(isinf(DS))   = nan;

r_s_tmp = r_s_tmp(1,:);
dr      = dr(1,:);
end