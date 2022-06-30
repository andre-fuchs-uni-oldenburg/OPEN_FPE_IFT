%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the independent trajectories data.
%
% Arguments IN
% data = filtered data 
% r = scale vector from the start to end of the cascade trajectory
% rind = index-vector from the start to end of the cascade trajectory
% stp = index used to cut v(x) in non-overlapping pieces of length equal to integral scale
% 
% Arguments OUT
% u = independent trajectories
% ind_trajec = index-vector of the trajectories (important vector for the recalculation of the trajectories from data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u,ind_trajec]=create_us(data,r,rind,stp)
i = 0; j = 1;
u           = nan( ceil(numel(data)/stp), numel(r) ); % initialise u-matrix
ind_trajec  = nan( ceil(numel(data)/stp),1 );

while i+rind(1)<numel(data)
	u(j,:)          = data(rind+i) - data(1+i); % velocity increment
    ind_trajec(j)   = i;
	i               = i+stp;
    j               = j+1;
end
% cut out nan's 
u                   = u(1:j-1,:); 
ind_trajec          = ind_trajec(1:j-1,:); 
end