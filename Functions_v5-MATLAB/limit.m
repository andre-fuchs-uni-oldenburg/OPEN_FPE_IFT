%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This function calculates the centers and the width of the bins for incr1 and incr2.
%
% Arguments IN
% incr1 = Velocity increment for step tau1
% incr2 = Velocity increment for step tau2
% increment_bin = number of bins
%
% Arguments OUT
% dx = bin width for incr1
% dy = bin width for incr2
% x = Centers of bins for incr1
% y = Centers of bins for incr2
% dA = dx*dy = area of a single 2D bin this will be used for converting probability to probability 
%              density function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [dx,dy,x,y,dA] = limit(incr1,incr2,increment_bin)  
[~,x]       = hist(incr1,increment_bin);
[~,y]       = hist(incr2,increment_bin);

dx          = mean(diff(x));
dy          = mean(diff(y));

dA          = dx*dy;

%     dx              = range(incr1)/(increment_bin-1);
%     dy              = range(incr2)/(increment_bin-1);
%     x               = min(incr1):dx:max(incr1); 
%     y               = min(incr2):dy:max(incr2);

%     x               = linspace(min(incr1),max(incr1),increment_bin); 
%     y               = linspace(min(incr2),max(incr2),increment_bin);

%     x               = linspace(mean(incr1)-3*std(incr1), mean(incr1)+3*std(incr1),increment_bin);
%     y               = linspace(mean(incr2)-3*std(incr2), mean(incr2)+3*std(incr2),increment_bin);
%     dx              = mean(diff(x));
%     dy              = mean(diff(y));
%     dA              = dx*dy;
end