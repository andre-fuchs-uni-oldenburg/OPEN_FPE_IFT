%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the total Entropy of the medium SM for trajectories.
% Comment: Ito convention is identical with the pre-point rule (rectangles to approximation), 
% and the Stratonovich convention is identical with the trapezoidal rule.
%
% Arguments IN
% ur_s = central difference quotient for the separation of scales/step increment 
%        (in normalized scale) referred to the sequence from large to small scales in the trajectory
% u_s = mid-point increment trajectories 
% r_s = mid-point scale vector (Stratonovich convention)
% dr = separation of scales/step increment (in normalized scale) referred to the sequence from 
%      large to small scales in the cascade trajectory
% co = Coefficients $d_{ij}(r)$ of the optimized Kramers-Moyal coefficients, where first entry is 
%      zeroth power and power laws for r-dependency with exponents co.ae and co.be, for each power 
%      in D1 and D2 respectively
% 
% Arguments OUT
% Sm_1 = rectangles to approximation
% Sm_2 = trapezoidal approximation
% F_tmp,D_tmp,FD,Sm_tmp (are not yet included/used in the current version)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sm_1,Sm_2,F_tmp,D_tmp,FD,Sm_tmp] = calcSm(ur_s, u_s, r_s, dr, co)
u_s(u_s==0)     = eps;
D1              = D1_fct(u_s,r_s,co);       % Drift 
D2              = D2_fct(u_s,r_s,co);       % Diffusion D2>0
D2_diff         = D2_fct_diff(u_s,r_s,co);  % D2'


% D1              = D1./r_s;       
% D2              = D2./r_s;       
% D2_diff         = D2_diff./r_s;  


F               = D1-(D2_diff./2);          % force  F = D1-D2'
F_tmp           = F;
D_tmp           = D2;
FD              = F./D2;
Sm_tmp          = ur_s.*F./D2.*dr;

%% rectangles to approximation
Sm_1            = nansum( (ur_s.*F./D2).* dr , 2); % stochastic integral, minus sign
% Sm_1 = nansum( smooth(ur_s.*F./D2,3).' , 2) * dr;
% Sm_1 = trapz(flip(r_s(~isnan(r_s))),ur_s(~isnan(ur_s)).*F(~isnan(F))./D2(~isnan(D2)));
% Sm_1 = real(Sm_1);
% parfor i=1:size(F,1)
%     Sm_1(i,1) = trapz(flip(r_s(i,:)),smooth(ur_s(i,:).*F(i,:)./D2(i,:),3).');
%     Sm_1(i,1) = trapz(flip(r_s(i,:)),ur_s(i,:).*F(i,:)./D2(i,:),2);
% end


%% trapezoidal approximation
Sm_2            = trapz(flip(r_s(1,:)),(ur_s.*F./D2),2);
% Sm_2          = cumtrapz((ur_s.*F./D2),2);

% figure
% plot(Sm_1(1:100))
% hold on
% plot(Sm_2(1:100))
% 
% nbins=301;
% [fuL, uL]     = hist(Sm_1,nbins);  
% puL           = fuL / ( sum(fuL)); 
% [ful, ul]     = hist(Sm_2,uL);
% pul           = ful / ( sum(ful));
% 
% figure
% plot(uL,puL)
% hold on
% plot(ul,pul)
% set(gca,'yScale','log')

clear ur_s u_s r_s dr co F D2
end


%% Drift
function D1 = D1_fct(u_s,r_s,co)
    D1   = zeros(size(u_s));   
    D1   = D1 + (r_s.^co.ea(1) .* co.a(1) + co.a(2)).*u_s.^1; % co.a(i) is a_(i-1), exponents accordingly
    %           - ( (r_s.^co.eb(3) .* co.b(5) + co.b(6)).*2.*u_s.^1 + ((r_s.^co.eb(2) .* co.b(3) + co.b(4)).*u_s.^0)); % substract u-derivative of D2
    D1              = real(D1);
    D1(isinf(D1))   = nan;       
end


%% Diffusion D2>0
function D2 = D2_fct(u_s,r_s,co)
    Nb  = numel(co.b); % number of coefficients in D2
    D2   = zeros(size(u_s));
    k=1;
    for i=1:2:Nb
        D2 = D2 + ( (r_s.^co.eb(k) .* co.b(i) + co.b(i+1)).*u_s.^(k-1) ); % co.d(i) is d_(i-1), exponents accordingly       end
        k = k+1;
    end
    D2           = real(D2);
    D2(isinf(D2)) = nan;
    % D2=D2+eps;
    D2(D2<=0)     = nan; % cut out zero and negative diffusion
    clear Nb i 
end


%% Derivative D2'
function D2_diff = D2_fct_diff(u_s,r_s,co)
    D2_diff = zeros(size(u_s));
    D2_diff = D2_diff + ((r_s.^co.eb(3) .* co.b(5) + co.b(6)).*2.*u_s.^1 + ((r_s.^co.eb(2) .* co.b(3) + co.b(4)).*u_s.^0)); 
    D2_diff = real(D2_diff);
    D2_diff(isinf(D2_diff))=nan;
end