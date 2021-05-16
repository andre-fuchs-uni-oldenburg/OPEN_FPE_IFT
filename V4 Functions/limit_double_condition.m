%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This function calculates the centers and the width of the bins for incr1, incr2 and incr3.
%
% Arguments IN
% incr1 = Velocity increment for step tau1
% incr2 = Velocity increment for step tau2
% incr3 = Velocity increment for step tau3
% increment_bin = number of bins
%
% Arguments OUT
% dx = bin width for incr1
% dy = bin width for incr2
% dz = bin width for incr3
% x = Centers of bins for incr1
% y = Centers of bins for incr2
% z = Centers of bins for incr3
% dA = dx*dy = area of a single 2D bin this will be used for converting probability to probability 
%              density function
% dV = dx*dy*dz = volume of a single 3D bin this will be used for converting probability to probability 
%                 density function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [dx,dy,dz,x,y,z,dA,dV] = limit_double_condition(incr1,incr2,incr3,increment_bin)
[~,x]           = hist(incr1,increment_bin);
[~,y]           = hist(incr2,increment_bin);
[~,z]           = hist(incr3,increment_bin);

dx              = mean(diff(x));
dy              = mean(diff(y));
dz              = mean(diff(z));

dA              = dx*dy;
dV              = dx*dy*dz;

%     dx              = range(incr1)/(increment_bin-1);
%     dy              = range(incr2)/(increment_bin-1);
%     dz              = range(incr3)/(increment_bin-1);
%     x               = min(incr1):dx:max(incr1); 
%     y               = min(incr2):dy:max(incr2);
%     z               = min(incr3):dz:max(incr3);
%     dA              = dx*dy;
%     dV              = dx*dy*dz;
    
%     x               = linspace(min(incr1),max(incr1),increment_bin); 
%     y               = linspace(min(incr2),max(incr2),increment_bin);
%     z               = linspace(min(incr3),max(incr3),increment_bin);

    
%      %%%%%%%%%%%%%%%%%% NEW NEW NEW NEW
%      tmp=[range(x) range(y) range(z)];
%      ind_tmp=find(tmp==max(tmp));
%      
%     if ind_tmp==1 
%         dx          = mean(diff(x));
%         y           = x;
%         z           = x;
%         dy          = dx;
%         dz          = dx;
%     elseif ind_tmp==2 
%         dy          = mean(diff(y));
%         x           = y;
%         z           = y;
%         dx          = dy;
%         dz          = dy;
%     else
%         dz          = mean(diff(z));
%         x           = z;
%         y           = z;
%         dx          = dz;
%         dy          = dz;
%     end
%     %%%%%%%%%%%%%%%%%% NEW NEW NEW NEW
end   