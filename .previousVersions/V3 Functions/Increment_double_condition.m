%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This function calculates three increment time series for the different lags tau1, tau2 and tau3.
%
% Arguments IN
% tau1 = first increments in steps
% tau2 = second increments in steps
% tau3 = third increments in steps
% d = data 
%
% Arguments OUT
% incr1 = Velocity increment for step tau1
% incr2 = Velocity increment for step tau2
% incr3 = Velocity increment for step tau3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [incr1,incr2,incr3] = Increment_double_condition(tau1,tau2,tau3,d)
incr1           = nan;
incr2           = nan;
incr3           = nan;

%%%%%%
% delta=max([tau1 tau2 tau3]);
% tmp1=tau3-tau1;
% tmp2=tau3-tau2;
% 
% incr1   =   (d(tau1+1:length(d)-tmp1) - d(1:length(d)-delta)).';
% incr2   =   (d(tau2+1:length(d)-tmp2) - d(1:length(d)-delta)).';
% incr3   =   (d(tau3+1:length(d))      - d(1:length(d)-delta)).';
%%%%%

%tau1<tau2<tau3
delta=max([tau1 tau2 tau3]);
if tau1>tau3
    tmp2    = tau1-tau2;
    tmp3    = tau1-tau3;
    
    incr1   =   (d(tau1+1:length(d))        - d(1:length(d)-delta)).';
    incr2   =   (d(tau2+1:length(d)-tmp2)   - d(1:length(d)-delta)).';
    incr3   =   (d(tau3+1:length(d)-tmp3)   - d(1:length(d)-delta)).';

else
    tmp2    = tau3-tau2;
    tmp3    = tau3-tau1;
    
    incr1   =   (d(tau1+1:length(d)-tmp3)   - d(1:length(d)-delta)).';
    incr2   =   (d(tau2+1:length(d)-tmp2)   - d(1:length(d)-delta)).';
    incr3   =   (d(tau3+1:length(d))        - d(1:length(d)-delta)).';
end 

% if  tau1>tau2&tau3;
%     delta=tau1;
%     tmp=abs(max([tau1 tau2 tau3])-min([tau1 tau2 tau3]));
%     incr1   =   (d(tau1+1:length(d))     - d(1:length(d)-delta)).';
%     incr2   =   (d(tau2+1:length(d)-tmp) - d(1:length(d)-delta)).';
%     incr3   =   (d(tau3+1:length(d)-tmp) - d(1:length(d)-delta)).';
% else %     tau1<tau2;
%     delta=tau2;
%     tmp=abs(max([tau1 tau2 tau3])-min([tau1 tau2 tau3]));
%     incr1   =   (d(tau1+1:length(d)-tmp) - d(1:length(d)-delta)).';
%     incr2   =   (d(tau2+1:length(d))     - d(1:length(d)-delta)).';
end