%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the total Entropy of the system or the Shanon entropy for trajectories.
%
% Arguments IN
% uLs = mid-point increment at the start of the cascade trajectory
% uls = mid-point increment at the end of the cascade trajectory
% puL = Probability density function &  uL = is the center of the bin at start of the cascade trajectory
% pul = Probability density function &  ul = is the center of the bin at end of the cascade trajectory
% 
% Arguments OUT
% Ds = Entropy of the system or the Shanon entropy
% The cascade process takes some initial probability s_L = $p(u_L, L)$ and changes to a different 
% final probability s_l = $p(u_\lambda, \lambda)$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ds,s_l,s_L]  = calcD_system_entropy(uLs,uls,uL,puL,ul,pul)
index_l         = (pul==0);
ul              = ul(~index_l);
pul             = pul(~index_l);

index_L         = (puL==0);
uL              = uL(~index_L);
puL             = puL(~index_L);

pul             = (smooth(pul,floor(0.1*length(pul)))).';
puL             = (smooth(puL,floor(0.1*length(puL)))).';

redind_fun      = @(x)~(isnan(x)) & ~(isinf(x)) & real(x)==x; % returns 1 only if x is not a NaN, not Inf and real
uLs_ind         = redind_fun(uLs); 
uls_ind         = redind_fun(uls);
redind          = uLs_ind + uls_ind == 2; % index pointing where both u(L) and u(l) are normal floats
uLs(~redind)    = NaN; % cut out NaNs, Infs and complex valued entries
uls(~redind)    = NaN;

%% ... maybe better to use Castaing's formula here
s_L = interp1(uL, puL, uLs, 'linear'); % plug u(L) and u(l) into p_L(u) and p_l(u) via linear interpolation ...
s_l = interp1(ul, pul, uls, 'linear'); % 
% s_L = interp1(uL, puL, uLs, 'linear','extrap'); % 
% s_l = interp1(ul, pul, uls, 'linear','extrap'); % 
% s_L = interp1(uL, puL, uLs, 'spline','extrap'); % 
% s_l = interp1(ul, pul, uls, 'spline','extrap'); %

% Ds  = -log10( s_l./s_L );
Ds  = -log( s_l./s_L );

%  figure
%  plot(s_L)
%  hold on
%  plot(s_l)
%  plot( s_l./s_L )
%  plot(-log( s_l./s_L ))
end