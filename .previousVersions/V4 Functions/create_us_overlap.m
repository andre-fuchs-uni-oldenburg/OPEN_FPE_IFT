%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the overlapping trajectories data.
%
% Arguments IN
% data = filtered data 
% r = scale vector from the start to end of the cascade trajectory
% rind = index-vector from the start to end of the cascade trajectory
% stp = index used to cut v(x) in non-overlapping pieces of length equal to integral scale
% 
% Arguments OUT
% u = overlapping trajectories
% ind_trajec = index-vector of the trajectories (important vector for the recalculation of the trajectories from data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u,ind_trajec]=create_us_overlap(data,r,rind,stp)
end_L           = size(data,1);
u               = nan( end_L-stp+1, numel(r) ); % initialise u-matrix
ind_trajec      = nan( end_L-stp+1,1 );

for i=0:end_L-stp
    u(i+1,:)        = data(rind+i) - data(1+i);           
    ind_trajec(1+i) = i;
end
end