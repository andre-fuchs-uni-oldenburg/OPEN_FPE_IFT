%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path integral formalism
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
% z = 1 ==> Overlapping trajectories 
% z = 3 ==> Independent trajectories
% 
% Arguments OUT
% A = action functional
% Lag,p,H,tmp,H1,H2 (are not yet included/used in the current version)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Lag,A,p,H,tmp,H1,H2] = calc_Path_integral(ur_s,u_s,r_s,dr,co,z)
D1              = D1_fct(u_s,r_s,co);           % Drift 
D2              = D2_fct(u_s,r_s,co);           % Diffusion D2>0

D1_diff           = D1_fct_diff(u_s,r_s,co);      % D1' 
% D1_diff_diff  = D1_fct_diff_diff(u_s,r_s,co); % D1''
D1_diff_diff    = 0;                            % D1''
D2_diff           = D2_fct_diff(u_s,r_s,co);      % D2'
D2_diff_diff    = D2_fct_diff_diff(u_s,r_s,co); % D2''


%% Lagrangian
% Lag       = ((ur_s-D1).^2)./(2.*D2);
% Lag       = ((ur_s-D1).^2)./(4.*D2);


%% Lagrangian 07/20
Lag         = (((ur_s-D1+(D2_diff./2)).^2)./(4.*D2))-(D1_diff./2);

%% Action: scale integral of the Lagrangian L
A           = nansum( Lag.* dr , 2) ;
% A       = trapz(flip(r_s(1,:)), Lag , 2);
% for i=1:size(Lag,1)
% %     A(i,1)      = -1.*nansum( log(D2(i,:))-cumtrapz(u_s(i,:),(D1(i,:)./D2(i,:))), 2) * dr;
%     A(i,1)      = trapz(r_s(1,:),log(D2(i,:))-cumtrapz(u_s(i,:),(D1(i,:)./D2(i,:))), 2);
% end

if z==3
    clear tmp
    tmp         = cumsum( Lag.* dr , 2) ;

    %% Conjugate variable
    % p       = (ur_s-D1)./(D2);
    % p       = (ur_s+D1)./(2.*D2);
    % p       = (ur_s+D1);

    %% Conjugate variable 07/20
    p           = (ur_s-D1+(D2_diff./2))./(2.*D2);

    %% Check
    % ur_s    = D2.*p+D1;
    % L       = D2.*p.^2;

    %% Hamiltonian: checken ob die Legendre Transformation so richitg ist
    % H       = D1.*p+D2.*p.^2;
    % % H       = p.*ur_s-Lag;
    % 
    % H1      = D1+(2.*D2.*p);                  % H1=d_r ur = ur_s 
    % H2      = (-p.*D1_diff)-(p.^2.*D2_diff);  % H2=d_r p


    %% Hamiltonian 07/20
    H           = (D1-(D2_diff./2)).*p + D2.*p.^2 + (D1_diff./2);

    % Hp
    H1          = D1+(2.*D2.*p)-(D2_diff./2);                                                              % H1 = d_r ur = ur_s 
    % Hu
    H2          = (-p.*(D1_diff-(D2_diff_diff./2)))-(p.^2.*D2_diff)-(D1_diff_diff./2);      % H2 = d_r p
%     H2          = - H2;
    % tmp     = D2.*p.^2;

    %% Hamiltonian 12/21 - nur minimaler Unterschied
    H           = D2.*p.^2 + (D1 - (D2_diff./2)).*p - (D1_diff./2);

    % Hp
    H1          = (2.*D2.*p) + D1 - (D2_diff./2);                                                              % H1 = d_r ur = ur_s 

    % Hu
    H2          = - ((p.^2.*D2_diff) + (D1_diff - (D2_diff_diff./2)).*p - (D1_diff_diff./2));  % H2 = d_r p
    H2          = - H2;

else
    Lag = nan;
    p   = nan;
    H   = nan;
    tmp = nan;
    H1  = nan;
    H2  = nan;  
end
clear  D1 D2 D1_diff D2_diff
end


%% Drift
function D1 = D1_fct(u_s,r_s,co)
    D1   = zeros(size(u_s));   
    D1   = D1 + (r_s.^co.ea(1) .* co.a(1) + co.a(2)).*u_s.^1; % co.a(i) is a_(i-1), exponents accordingly
    %           - ( (r_s.^co.eb(3) .* co.b(5) + co.b(6)).*2.*u_s.^1 + ((r_s.^co.eb(2) .* co.b(3) + co.b(4)).*u_s.^0)); % substract u-derivative of D2
    D1              = real(D1);
    D1(isinf(D1))   = nan;       
end


%% Derivative D1'
function D1_diff = D1_fct_diff(u_s,r_s,co)
    D1_diff = zeros(size(u_s)); 
    D1_diff = D1_diff + (r_s.^co.ea(1) .* co.a(1) + co.a(2)).*u_s.^0; % co.a(i) is a_(i-1), exponents accordingly
    D1_diff                     = real(D1_diff);
    D1_diff(isinf(D1_diff))     = nan;       
    clear Na Nb i 
end

% 2. Derivative D1''
% function D1_diff_diff = D1_fct_diff_diff(u_s,r_s,co)
%     D1_diff_diff = zeros(size(u_s)); 
%     D1_diff_diff = D1_diff_diff + (r_s.^co.ea(1) .* co.a(1) + co.a(2)).*u_s.^0; % co.a(i) is a_(i-1), exponents accordingly
% D1_diff_diff        = real(D1_diff_diff);
% D1_diff_diff(isinf(D1_diff_diff)) = nan;       
% clear Na Nb i 
% end


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


%% 2. Derivative D2''
function D2_diff_diff = D2_fct_diff_diff(u_s,r_s,co)
    D2_diff_diff = zeros(size(u_s));
    D2_diff_diff = D2_diff_diff + ((r_s.^co.eb(3) .* co.b(5) + co.b(6)).*2.*u_s.^0); 
    D2_diff_diff = real(D2_diff_diff);
    D2_diff_diff(isinf(D2_diff_diff))=nan;
end