%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs the surface fit with a linear function for $D^{(1)}\left(u_r,r\right)$.
% 
% Arguments IN
% xx = velocity increment vector (center of bin)
% r_norm = Normalized scale vector from the start to end of the cascade trajectory, 
% z = Drift coefficient
% ze = error associated with drift coefficient
% 
% Arguments OUT
% fitresult_D1 = Coefficients $d_{ij}(r)$ of the optimized Kramers-Moyal coefficients using the
%                surface fits: D1
% xData, yData, zData, weights = xx, r_norm, z, ze reshape into one single vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fitresult_D1, xData, yData, zData, weights] = FIT_D1(xx, r_norm, z, ze, varargin )
%% D1-Fit
[xData, yData, zData, weights] = prepareSurfaceData( xx, r_norm, z, ze );
% Set up fittype and options.
ft = fittype( '(a1*r_norm^ra1+a11)*x', 'independent', {'x', 'r_norm'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.DiffMinChange = 1e-12;
% opts.Display = 'Off';
% opts.MaxFunEvals = 1000;
% opts.MaxIter = 1000;
% opts.TolFun = 1e-12;
% opts.TolX = 1e-12;
opts.Weights = weights;

% %--------------------------------------------------
% %--------------------------------------------------
%%%%%  2019 %%%%%%%%
% d11, must be negative
% opts.StartPoint   = [-1E-1 1E-3 -1];
% opts.Lower        = [-Inf -Inf -Inf];
% opts.Upper        = [0     Inf   0];

if nargin >= 5  
   co_KM_opti       = varargin{1};
   opts.StartPoint  = [co_KM_opti.a(1) co_KM_opti.a(2) co_KM_opti.ea(1)];
else
    % opts.StartPoint = [-0.1 0.1 -1];
    opts.StartPoint = [0 0 0];
%     %K41
%     opts.StartPoint = [-0.33 0 -1];
%     %K62
%     opts.StartPoint = [-0.36 0 -1];
end

%% Lower bound and Upper bound
% (a1*r_norm^ra1+a11)*xx
%                  a1  a11  ra1
% opts.Lower    = [-2   -1  -2];
% opts.Upper    = [ 0    1   0];
% opts.Lower    = [-2   -2  -2];
% opts.Upper    = [ 0    2   0];
% opts.Lower    = [-1   -1  -Inf];
% opts.Upper    = [1     1   0];
opts.Lower      = [-2   -2  -Inf];
opts.Upper      = [2     2   0];

%% constant in scale
% opts.Lower      = [-2    0   0];
% opts.Upper      = [2     0   0];


% K41 
% opts.Lower = [-1 0 -Inf];
% opts.Upper = [1 0 0];

% Fit model to data.
[fitresult_D1, gof] = fit( [xData, yData], zData, ft, opts );
end