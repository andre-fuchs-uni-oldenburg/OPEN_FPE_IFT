%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs the follwing fit to the coefficients $d_{ij}(r)$:
% $\alpha (r/\lambda)^{\beta}+\gamma$. 
% 
% Arguments IN
% r_norm = Normalized scale vector from the start to end of the cascade trajectory, 
% y = Coefficients $d_{ij}(r)$
% t = length(r_norm)
% 
% Arguments OUT
% fit_dij = Coefficients $d_{ij}(r)$ of the optimized Kramers-Moyal coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit_d20, fit_d21, fit_d22 ] = FIT_d2j(r_norm, y, t, varargin )
%% D2-Fit
%% d20
[xData, yData] = prepareCurveData(r_norm,y(t*1+1:(t*1+1)+t-1));
ft = fittype( '(b0*x^rb0+b00)',  'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.DiffMinChange = 1e-12;
opts.Display = 'Off';
opts.MaxFunEvals = 1000;
opts.MaxIter = 1000;
opts.TolFun = 1e-12;
opts.TolX = 1e-12;
% b0 b00 rb0
if nargin >= 4  
    co_KM_opti       = varargin{1};
    opts.StartPoint = [co_KM_opti.b(1) co_KM_opti.b(2) co_KM_opti.eb(1)];
else
%     opts.StartPoint = [0.02 -0.02 0.001];
    opts.StartPoint = [0 0 0];
end  
% Lower bound and Upper bound
% (b0*r_norm^rb0+b00)
%                  b0  b00 rb0
% opts.Lower    = [ 0  -2   0];
% opts.Upper    = [ 2   2   2];
opts.Lower      = [-1  -1   0];
opts.Upper      = [1    1   Inf];

%% constant in scale
% opts.Lower      = [-1   0   0];
% opts.Upper      = [1    0   0];

[fit_d20, gof] = fit(xData, yData, ft, opts );


%% d21
[xData, yData] = prepareCurveData(r_norm,y(t*2+1:(t*2+1)+t-1));
ft = fittype( '(b1*x^rb1+b11)',  'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.DiffMinChange = 1e-12;
opts.Display = 'Off';
opts.MaxFunEvals = 1000;
opts.MaxIter = 1000;
opts.TolFun = 1e-12;
opts.TolX = 1e-12;
% b1 b11 rb1
% opts.StartPoint = [co.b(3) co.b(4) co.eb(2)];
if nargin >= 4  
    co_KM_opti       = varargin{1};
    opts.StartPoint = [co_KM_opti.b(3) co_KM_opti.b(4) co_KM_opti.eb(2)];
else
%     opts.StartPoint = [0.02 -0.02 -1];
    opts.StartPoint = [0 0 0];
end   
% Lower bound and Upper bound
% (b1*r_norm^rb1+b11)
%                  b1 b11 rb1
% opts.Lower    = [-1 -1  -2];
% opts.Upper    = [ 1  1   2];
opts.Lower      = [-1 -1 -Inf];
opts.Upper      = [ 1  1  Inf];

%% constant in scale
% opts.Lower      = [-1   0   0];
% opts.Upper      = [1    0   0];

[fit_d21, gof] = fit(xData, yData, ft, opts );


%% d22
[xData, yData] = prepareCurveData(r_norm,y(t*3+1:(t*3+1)+t-1));
ft = fittype( '(b2*x^rb2+b22)',  'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.DiffMinChange = 1e-12;
opts.Display = 'Off';
opts.MaxFunEvals = 1000;
opts.MaxIter = 1000;
opts.TolFun = 1e-12;
opts.TolX = 1e-12;
% b2 b22 rb2
% opts.StartPoint = [co.b(5) co.b(6) co.eb(3)];
 if nargin >= 4  
    co_KM_opti       = varargin{1};
    opts.StartPoint = [co_KM_opti.b(5) co_KM_opti.b(6) co_KM_opti.eb(3)];
else
%     opts.StartPoint = [0.007 -0.007 -0.3];
    opts.StartPoint = [0 0 0];
 end
% Lower bound and Upper bound
% (b2*r_norm^rb2+b22)
%                  b2 b22 rb2
% opts.Lower    = [0  -1  -2];
% opts.Upper    = [1   1   0];
opts.Lower      = [-1 -1 -Inf];
opts.Upper      = [1   1   0];

%% constant in scale
% opts.Lower      = [-1   0   0];
% opts.Upper      = [1    0   0];

[fit_d22, gof] = fit(xData, yData, ft, opts );
end