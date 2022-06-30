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
% u = cascade trajectories for which the total entropy shall be calculated 
% r = scale vector from the start to end of the cascade trajectory
% tmp_size = number of cascade trajectories
% r_norm = Normalized scale vector from the start to end of the cascade trajectory
% t                       = length(r_norm)
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
% d_M = error function:\xi =|1-\langle e^{\mathrm{-}\Delta S_{tot}} \rangle_{max(N)}|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d_M=DIRECT_IFT_OPTIMIZATION_dij(x0,u,r,tmp_size,r_norm,t,co_KM_opti,markov,taylor_L,Fs,m_data,trajec,dr_ind,z)
[fit_d11]                   = FIT_d1j(r_norm,x0(1:t),co_KM_opti);
[fit_d20, fit_d21, fit_d22] = FIT_d2j(r_norm,x0,t,co_KM_opti);

co.a  = [fit_d11.a1 fit_d11.a11];
co.ea = [fit_d11.ra1];
co.b  = [fit_d20.b0 fit_d20.b00 fit_d21.b1 fit_d21.b11 fit_d22.b2 fit_d22.b22];
co.eb = [fit_d20.rb0 fit_d21.rb1 fit_d22.rb2]; 

clear DS
DS                      = nan( tmp_size, 1 );
[Sm_tmp,Ds_tmp,DS_tmp]  = calcDS(z,u,r.',co,markov,taylor_L,Fs,m_data,dr_ind);

Sm(1:sum(~isnan(DS_tmp),1),1)           = Sm_tmp(~isnan(DS_tmp));
Ds(1:sum(~isnan(DS_tmp),1),1)           = Ds_tmp(~isnan(DS_tmp));
DS(1:sum(~isnan(DS_tmp),1),1)           = DS_tmp(~isnan(DS_tmp));

clear Sm_tmp Ds_tmp DS_tmp

Sm      = reshape(Sm,[],1);
Ds      = reshape(Ds,[],1);
DS      = reshape(DS,[],1);

DS(isinf(DS))   = nan;
Sm              = Sm(~isnan(DS));
Ds              = Ds(~isnan(DS));
DS              = DS(~isnan(DS));

while mean( exp( -DS ) )== Inf 
%     | mean(Sm)
    quant   = 1;
    Sm      = Sm(exp(-DS)<=quantile(-DS,quant));
    Ds      = Ds(exp(-DS)<=quantile(-DS,quant));
    DS      = DS(exp(-DS)<=quantile(-DS,quant));
    quant   = quant-0.0001;
end

if ceil(0.6*length(DS))>5000
    d_M_60                  = abs(1-nanmean( exp( -DS(1:ceil(0.6*length(DS))) )));
else
    d_M_60=0;
end

if ceil(0.7*length(DS))>5000   
    d_M_70                  = abs(1-nanmean( exp( -DS(1:ceil(0.7*length(DS))) )));
else
    d_M_70=0;
end

if ceil(0.8*length(DS))>5000
    d_M_80                  = abs(1-nanmean( exp( -DS(1:ceil(0.8*length(DS))) )));
else
    d_M_80=0;
end

d_M_90                  = abs(1-nanmean( exp( -DS(1:ceil(0.9*length(DS))) )));
d_M_100                 = abs(1-nanmean( exp( -DS(1:(length(DS))) )));

d_M                     = d_M_60+d_M_70+d_M_80+d_M_90+d_M_100
% d_M                   = d_M_100;
% d_M                   = abs(1-nanmean( exp( -DS )));
end