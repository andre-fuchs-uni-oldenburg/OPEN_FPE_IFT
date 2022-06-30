%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This function calculates two increment time series for the two different lags tau1 and tau2.
%
% Arguments IN
% tau1 = first increments in steps
% tau2 = second increments in steps
% d = data 
%
% Arguments OUT
% incr1 = Velocity increment for step tau1
% incr2 = Velocity increment for step tau2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [incr1,incr2] = Increment(tau1,tau2,d)
incr1           = nan;
incr2           = nan;

if  tau1>tau2
    delta=tau1;
    tmp=tau1-tau2;
    incr1   =   (d(tau1+1:length(d))     - d(1:length(d)-delta)).';
    incr2   =   (d(tau2+1:length(d)-tmp) - d(1:length(d)-delta)).';
    
else %     tau1<tau2;
    delta=tau2;
    tmp=tau2-tau1;
    incr1   =   (d(tau1+1:length(d)-tmp) - d(1:length(d)-delta)).';
    incr2   =   (d(tau2+1:length(d))     - d(1:length(d)-delta)).';

% nbins=1001;
% [fuL uL] = hist(incr1,nbins); % compute frequency of u-values at scale r=L
% puL = fuL / ( nansum(fuL) * mean(diff(uL)) ); % relative frequency to approximate p(u,r=L)
% 
% figure
% plot(uL,puL)
% hold on
% [ful ul] = hist(incr2,nbins); % compute frequency of u-values at scale r=lambda
% pul = ful / ( nansum(ful) * mean(diff(ul)) ); % relative frequency to approximate p(u,r=lambda)
% plot(ul,pul)
% set(gca,'yScale','log') 
end