%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs the surface fit with a parabolic function for $D^{(2)}\left(u_r,r\right)$. 
% 
% Arguments IN
% r_norm = Normalized scale vector from the start to end of the cascade trajectory, 
% z = Diffusion coefficient coefficient
% ze = error associated with Diffusion coefficient coefficient
% 
% Arguments OUT
% fitresult_D2 = Coefficients $d_{ij}(r)$ of the optimized Kramers-Moyal coefficients using the
%                surface fits: D2
% xData, yData, zData, weights = xx, r_norm, z, ze reshape into one single vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fitresult_D2, xData, yData, zData, weights] = FIT_D2(xx, r_norm, z, ze, varargin )
%% D2-Fit
%D2(xi,r) Fit
[xData, yData, zData, weights] = prepareSurfaceData( xx, r_norm, z, ze );
% [xData, yData, zData] = prepareSurfaceData( xx, r_norm, z);

ft = fittype( '(b0*r_norm^rb0+b00) + (b1*r_norm^rb1+b11)*x + (b2*r_norm^rb2+b22)*x^2', 'independent', {'x', 'r_norm'}, 'dependent', 'z' );
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
if nargin >= 5  
   co_KM_opti       = varargin{1};
   opts.StartPoint  = [co_KM_opti.b(1) co_KM_opti.b(2) co_KM_opti.b(3) co_KM_opti.b(4) co_KM_opti.b(5) co_KM_opti.b(6) co_KM_opti.eb(1) co_KM_opti.eb(2) co_KM_opti.eb(3)];
else
%     opts.StartPoint = [0.02 -0.02 0.02 -0.02 0.007 -0.007 0.001 -1 -0.3];
    %                  b0  b00   b1   b11  b2  b22  rb0  rb1  rb2
    opts.StartPoint  = [0   0     0    0    0   0    0    0    0];   
    %K62
%     opts.StartPoint = [0 0 0 0 -0.33 0 0 0 -1];
end
% Lower bound and Upper bound
% (b0*r_norm^rb0+b00) + (b1*r_norm^rb1+b11)*x + (b2*r_norm^rb2+b22)*x^2
%              b0  b00   b1   b11   b2  b22   rb0    rb1   rb2
% opts.Lower = [0   -1   -1    -1   0   -1     0      -2   -2];
% opts.Upper = [1    1    1     1   1    1     1       2    0];
% opts.Lower = [0   -2   -1    -1   0   -1     0      -2   -2];
% opts.Upper = [2    2    1     1   1    1     2       2    0];

%             b0  b00   b1   b11   b2  b22   rb0    rb1   rb2
opts.Lower = [-1   -1   -1    -1   -1   -1     0   -Inf  -Inf];
opts.Upper = [ 1    1    1     1    1    1   Inf    Inf    0];

%d21<0
%             b0  b00   b1   b11   b2  b22   rb0    rb1   rb2
% opts.Lower = [-1   -1   -1    -1   -1   -1     0   -Inf  -Inf];
% opts.Upper = [ 1    1    0     0    1    1   Inf    Inf    0];


%d21=0
%             b0  b00   b1   b11   b2  b22   rb0    rb1   rb2
% opts.Lower = [-1   -1    0    0   -1   -1     0      0  -Inf];
% opts.Upper = [ 1    1    0    0    1    1   Inf      0    0];
% 
%d21=0, d22=0
% %             b0  b00   b1   b11   b2  b22   rb0    rb1   rb2
% opts.Lower = [-1   -1    0    0    0    0     0      0     0];
% opts.Upper = [ 1    1    0    0    0    0   Inf      0     0];

%K62
%             b0  b00  b1  b11    b2   b22   rb0   rb1   rb2
% opts.Lower = [0   0    0   0   -1    -1    0      0   -Inf];
% opts.Upper = [0   0    0   0    1     1    0      0     0];


% constant in scale
%             b0  b00   b1   b11   b2  b22   rb0    rb1   rb2
% opts.Lower = [-1   0   -1     0   -1    0     0      0     0];
% opts.Upper = [ 1   0    1     0    1    0     0      0     0];


[fitresult_D2, gof] = fit( [xData, yData], zData, ft, opts );

% a1=-2; a11=0; ra1=-2;
% figure
% plot(logspace(log10(min(r_norm(:,1))),log10(max(r_norm(:,1))),100),...
% a1.*logspace(log10(min(r_norm(:,1))),log10(max(r_norm(:,1))),100).^ra1+a11)
% % plot(logspace(log10(min(r_norm(:,1))),log10(max(r_norm(:,1))),100),...
% %     a1.*logspace(log10(min(r_norm(:,1))),log10(max(r_norm(:,1))),100).^ra1+a11)
% loglog(logspace(log10(min(r_norm(:,1))),log10(max(r_norm(:,1))),100),...
%     fitresult_D1.a1.*logspace(log10(min(r_norm(:,1))),log10(max(r_norm(:,1))),100).^fitresult_D1.ra1+fitresult_D1.a11)
% hold on
% loglog(logspace(log10(min(r_norm(:,1))),log10(max(r_norm(:,1))),100),...
%     -0.0818.*logspace(log10(min(r_norm(:,1))),log10(max(r_norm(:,1))),100).^-0.75)
% figure
% scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
% hold on
% surf(repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(r_norm(:,1))),log10(max(r_norm(:,1))),100).',1,100),...
% (a1.*repmat(logspace(log10(min(r_norm(:,1))),log10(max(r_norm(:,1))),100).',1,100).^ra1+a11).*repmat(linspace(min(xData),max(xData),100),100,1),'EdgeColor','none')
end