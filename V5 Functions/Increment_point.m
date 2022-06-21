%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This function calculates two increment time series for the two different lags tau1 and tau2,
% taking into account the multi-point analysis.
%
% Arguments IN
% tau1 = first increments in steps
% tau2 = second increments in steps
% d = data 
% condition = condition for multi-point analysis
% tol = This input is for multipoint statistics
%
% Arguments OUT
% incr1 = Velocity increment for step tau1
% incr2 = Velocity increment for step tau2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [incr1,incr2] = Increment_point(tau1,tau2,d,condition,tol)
incr1           = nan;
incr2           = nan;
if  tau1>tau2
    delta = tau1;

    index = find(d(1:end-delta)>-tol+condition & d(1:end-delta)<tol+condition);
    
    incr1   =   (d(tau1+index) - d(index)).';
    incr2   =   (d(tau2+index) - d(index)).';

else %     tau1<tau2;
    delta = tau2;
    
    index = find(d(1:end-delta)>-tol+condition & d(1:end-delta)<tol+condition);
    
    incr1   =   (d(tau1+index) - d(index)).';
    incr2   =   (d(tau2+index) - d(index)).';
end
end