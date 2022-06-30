%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction of the conditional probability density functions $p\left(u_{r'} | u_r \right)$
% via the short time propagator \cite{Risken}
%
% Arguments IN
% y = 1D array of centers of bins of incr1
% x = 1D array of centers of bins of incr1
% tau1 =   is r1 ===> (r2>r1) 
% tau2 =   is the r in number of samples==>is the scale in meters at which moments will be calculated and hence this r will
% be the same at which D1 & D2 will be calculated
% D1_poly = D1 for the reconstruction
% D2_poly = D2 for the reconstruction
% taylor_L = Taylor length scale in meters
% m_data = mean of the data
% Fs = Acquisition/Sampling Frequency in Hz
%
% Arguments OUT
% P = Reconstruction of the conditional probability density functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = ShortTimeProp(y,x,tau1,tau2,D1_poly,D2_poly,taylor_L,m_data,Fs,norm_ur,norm_r)
%     D1_poly=0.5.*D1_poly;
%     D2_poly=0.5.*D2_poly;
% delta          = tau/(tau2);
% delta          = markov;
if norm_r==1
step_taylor_L   = (taylor_L/(m_data/Fs));
else
step_taylor_L   = 1;    
end
% step_taylor_L=ceil(Fs*taylor_L/m_data);  
delta           = abs(tau2-tau1)/step_taylor_L;
%Renner S.37 Dissertation

% y_mean_bin(counter_B>=min_events),x_mean_bin(counter_A>=min_events),P_AIB(counter_A>=min_events,counter_B>=min_events)
P = zeros(length(x),length(y));
for s = 1:length(y)
    for q = 1:length(x)      
            P(q,s) = 1./(2.*sqrt(pi.*D2_poly(s).*delta)).*...
                     exp(-(((x(q)-y(s)-D1_poly(s).*delta).^2)./ (4.*D2_poly(s).*delta)));
    end
end
P=real(P);
end