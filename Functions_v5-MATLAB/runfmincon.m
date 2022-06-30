%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The optimization procedure systematically changes $D^{(1,2)}\left(u_r,r\right)$ until the error function 
% \begin{eqnarray}
% \xi =|1-\langle e^{\mathrm{-}\Delta S_{tot}} \rangle_{max(N)}|
% \end{eqnarray}
% is minimized. 
% 
% Arguments IN
% x0 = Initial D1 and D2 before optimization 
% lb = Lower bound
% ub = Upper bound
% u = cascade trajectories for which the total entropy shall be calculated 
% r = scale vector from the start to end of the cascade trajectory
% tmp_size = number of cascade trajectories
% iter = The maximum number of iteration used for the optimization algorithms
% r_norm = Normalized scale vector from the start to end of the cascade trajectory
% t = length(r_norm)
% co_KM_opti = Coefficients $d_{ij}(r)$ of the optimized Kramers-Moyal coefficients towards the best 
%              Fokker-Planck equation to reproduce the conditional PDF's using the surface fits
% markov = markov length in number of samples
% taylor_L = Taylor length scale in meters
% Fs = Acquisition/Sampling Frequency in Hz
% m_data = mean of the data
% trajec = 1 ==> The start/end of the cascade trajectory will be adjusted.
% dr_ind = Separation of scales/step increment (in samples) referred to the sequence from large to small scales in the cascade trajectory
% z = 3 ==> Independent trajectories
% 
% Arguments OUT
% co_IFT_opti = Coefficients $d_{ij}(r)$ of the optimized Kramers-Moyal coefficients towards the 
%               integral fluctuation theorem using the surface fits
% history = details of the optimization for every single iteration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [history] = runfmincon(x0,lb,ub,u,r,tmp_size,iter,r_norm,t,co_KM_opti,markov,taylor_L,Fs,m_data,trajec,dr_ind,z)
% Set up shared variables with OUTFUN
history.x1_iter     = [];
history.fval        = [];
history.IFT         = [];

%% Optimization statement
% call optimization
% 'PlotFcn','optimplotx', ...

options = optimoptions(@fmincon,'GradObj','off','Hessian','off','TolX',1e-10,'TolFun',1e-10,...
'OutputFcn',@outfun,...
'PlotFcn','optimplotfval', ...
'MaxIter',iter, 'MaxFunEvals',1e14,'FunValCheck','on','Algorithm','interior-point','Display','iter-detailed');

[x1,fval,exitflag,output] = ...
fmincon(@(x0)DIRECT_IFT_OPTIMIZATION_dij(x0,u,r,tmp_size,r_norm,t,co_KM_opti,markov,taylor_L,Fs,m_data,trajec,dr_ind,z),x0,[],[],[],[],lb,ub,[],options);

function stop = outfun(x,optimValues,state)
    stop = false;

     switch state
         case 'init'
             hold on
         case 'iter'
         % Concatenate current point and objective function value with history.
            history.fval = [history.fval; optimValues.fval];
            history.x1_iter = [history.x1_iter; x];
         % Concatenate current search direction with searchdir.
         % searchdir = [searchdir; optimValues.searchdirection'];

         
%% IFT Value concatenate current point
            [fit_d11]                   = FIT_d1j(r_norm,x(1:t),co_KM_opti);
            [fit_d20, fit_d21, fit_d22] = FIT_d2j(r_norm,x,t,co_KM_opti);

            co.a  = [fit_d11.a1 fit_d11.a11];
            co.ea = [fit_d11.ra1];
            co.b  = [fit_d20.b0 fit_d20.b00 fit_d21.b1 fit_d21.b11 fit_d22.b2 fit_d22.b22];
            co.eb = [fit_d20.rb0 fit_d21.rb1 fit_d22.rb2]; 

            DS                      = nan( tmp_size, 1 );
            [~,~,DS_tmp,~,~,~]      = calcDS(z,u,r.',co,markov,taylor_L,Fs,m_data,dr_ind);
            
            DS_tmp(isinf(DS_tmp))           = nan;
            DS(1:sum(~isnan(DS_tmp),1),1)   = DS_tmp(~isnan(DS_tmp));
            DS                              = reshape(DS,[],1);
            DS(isinf(DS))                   = nan;
            DS                              = DS(~isnan(DS));

            
            d_M_100     = nanmean( exp(-DS));
            history.IFT = [history.IFT; d_M_100];
         case 'done'
             hold off
         otherwise
     end
end
set(gca,'yScale','log')
history.iter = [0:output.iterations].';
end