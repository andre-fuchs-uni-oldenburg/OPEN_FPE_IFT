%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs the follwing fit to the coefficients $d_{ij}(r)$:
% $\alpha (r/\lambda)^{\beta}+\gamma$. 
% 
% Arguments IN
% r_norm = Normalized scale vector from the start to end of the cascade trajectory, 
% y = Coefficients $d_{ij}(r)$
% 
% Arguments OUT
% fit_dij = Coefficients $d_{ij}(r)$ of the optimized Kramers-Moyal coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit_d11] = FIT_d1j(r_norm, y, varargin )
%% d11-Fit
[xData, yData] = prepareCurveData(r_norm, y);
% Set up fittype and options.
ft = fittype( '(a1*x^ra1+a11)',  'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.DiffMinChange = 1e-12;
opts.Display = 'Off';
opts.MaxFunEvals = 1000;
opts.MaxIter = 1000;
opts.TolFun = 1e-12;
opts.TolX = 1e-12;
%%%%%%%%%%% opts.Weights = weights;
% %--------------------------------------------------
% %--------------------------------------------------
%%%%%  2019 %%%%%%%%
% d11, must be negative
% opts.StartPoint   = [-1E-1 1E-3 -1];
% opts.Lower        = [-Inf -Inf -Inf];
% opts.Upper        = [0     Inf   0];

if nargin >= 3  
   co_KM_opti       = varargin{1};
   opts.StartPoint  = [co_KM_opti.a(1) co_KM_opti.a(2) co_KM_opti.ea(1)];
else
    % opts.StartPoint = [-0.1 0.1 -1];
    opts.StartPoint = [0 0 0];
end

%% Lower bound and Upper bound
% (a1*r_norm^ra1+a11)
%                  a1  a11  ra1
% opts.Lower    = [-2   -2  -2]; 
% opts.Upper    = [ 0    2   0];
% opts.Lower    = [-1   -1 -Inf];
% opts.Upper    = [ 1    1   0];
opts.Lower      = [-2   -2 -Inf]; 
opts.Upper      = [2     2   0];


%% constant in scale
% opts.Lower      = [-2   0   0];
% opts.Upper      = [2    0   0];


% Fit model to data.
[fit_d11, gof] = fit(xData, yData, ft, opts );
end