%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arguments IN
% tau1 = first increments in steps
% tau2 = second increments in steps
% incr1 = Velocity increment for step tau1
% incr2 = Velocity increment for step tau2
% num_bin = number of bins
% dx = bin width for incr1
% dy = bin width for incr2
% x = Centers of bins for incr1
% y = Centers of bins for incr2
% dA = dx*dy = area of a single 2D bin this will be used for converting probability to probability 
%              density function
% varargin = Refer to MATLAB documentation
%
% Arguments OUT
% P_AIB = PDF of A i.e incr1 conditioned on B i.e incr2
% P_BIA = PDF of B i.e incr2 conditioned on A i.e incr1
% P_AnB = Joint PDF of A i.e. incr1 and B i.e incr2
% P_A = PDF of A i.e incr1
% P_B = PDF of B i.e incr2
% binP_AIB = 2D binning 
% x_mean_bin = 1D array of centers of bins of incr2
% y_mean_bin = 1D array of centers of bins of incr1
% events = number of elements that fall at the center in each bin
% counter_A = histogram of A
% counter_B = histogram of B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P_AIB,P_BIA,P_AnB,P_A,P_B,binP_AIB,x_mean_bin,y_mean_bin,events,counter_A,counter_B] = distribution(tau1,tau2,incr1,incr2,num_bin,dx,dy,x,y,dA,varargin)
% P_AIB(counter_A>=min_events,counter_B>=min_events)
scatter                        = nan;  

% hist3 creates a bivariate histogram using x y binning
if tau1>tau2
    scatter(1:length(incr1),1)     = incr2(1,1:length(incr1));
    scatter(1:length(incr1),2)     = incr1(1,1:length(incr1));
    counter2                       = hist3(scatter,{y x}); 
else
    scatter(1:length(incr1),1)     = incr1(1,1:length(incr1));
    scatter(1:length(incr1),2)     = incr2(1,1:length(incr1));
    counter2                       = hist3(scatter,{x y}); 
end

for i=1:num_bin
    counter_A(i) = sum(counter2(i,:));
    counter_B(i) = sum(counter2(:,i));
end

% sum(sum(counter2,2))
P_AnB                               = counter2./(length(scatter)*dA); % [p(AnB)]joint probability density
% P_AnB                               = counter2./(size(scatter,1)*size(scatter,2)*dA); % [p(AnB)]joint probability density      
% P_AnB                               = counter2./(length(scatter)*dx); % [p(AnB)]joint probability density     

P_A                                 = nan; 
for j           = 1:num_bin
    if tau1>tau2
        P_A(j)                          = nansum(P_AnB(j,:)).*dx;  
    else
        P_A(j)                          = nansum(P_AnB(j,:)).*dy; 
    end
end

P_B                                 = nan; %[p(B)] PDF of condition  
%P_B = hist(scatter(:,2),y)./(length(scatter)*dy);
for j           = 1:num_bin
    if tau1>tau2
        P_B(j)                          = nansum(P_AnB(:,j)).*dy;  
    else
        P_B(j)                          = nansum(P_AnB(:,j)).*dx; 
    end
end

P_AIB(num_bin,num_bin)                    = nan; %[p(A|B)]
for j           = 1:num_bin
    if P_B(j) > 0 
        P_AIB(:,j)                     = P_AnB(:,j)./P_B(1,j); 
    else 
        P_AIB(:,j)                     = P_AnB(:,j)./(P_B(1,j)+eps);
    end
end


P_BIA(num_bin,num_bin)                    = nan; %[p(B|A)]
for j           = 1:num_bin
    if P_A(j) > 0 
        P_BIA(:,j)                     = P_AnB(:,j)./P_A(1,j); 
    else 
        P_BIA(:,j)                     = P_AnB(:,j)./(P_A(1,j)+eps);
    end
end


if nargin <11
    %%%%%%%%%% Mean of specific bins
    end_L   = length(scatter);
    tmp1    = scatter(:,1);
    tmp2    = scatter(:,2);
    %     x_bin           = min(scatter(:,1)):(range(scatter(:,1))/(aufl)):max(scatter(:,1));
    %     y_bin           = min(scatter(:,2)):(range(scatter(:,2))/(aufl)):max(scatter(:,2));
        x_bin           = linspace(min(scatter(:,1)),max(scatter(:,1)),num_bin+1);
        y_bin           = linspace(min(scatter(:,2)),max(scatter(:,2)),num_bin+1);
    %     if length(x_bin)==aufl+1 & length(y_bin)==aufl+1
            x_mean_bin(1:num_bin) = nan;
            y_mean_bin(1:num_bin) = nan;
            events(1:num_bin)=0;

    %         for q=1:aufl
    %               index_r1=0;
    %               index_r2=0;
    %               
    %               if q<aufl
    %               index_r1 = find(tmp1(1:end_L)>=x_bin(q) & tmp1(1:end_L)<(x_bin(q+1)));
    %               index_r2 = find(tmp2(1:end_L)>=y_bin(q) & tmp2(1:end_L)<(y_bin(q+1)));
    %               else
    %               index_r1 = find(tmp1(1:end_L)>=x_bin(q) & tmp1(1:end_L)<=(x_bin(q+1)));
    %               index_r2 = find(tmp2(1:end_L)>=y_bin(q) & tmp2(1:end_L)<=(y_bin(q+1)));
    %               end
    %               
    %               if length(index_r1)>=1;               
    %               events(q)=length(index_r1);
    %               x_mean_bin(q)=nanmean(tmp1(index_r1)); 
    %               end
    %               if length(index_r2)>=1;
    %               y_mean_bin(q)=nanmean(tmp2(index_r2));
    %               end       
    %         end

             for q=1:num_bin
                if q<num_bin     
                  x_mean_bin(q) = nanmean(tmp1(tmp1>=x_bin(q) & tmp1<x_bin(q+1)));
    %               y_mean_bin(q)=nanmean(tmp2(tmp2>=y_bin(q) & tmp2<y_bin(q+1)));  
    %               events(q)=length(tmp1(tmp1>=x_bin(q) & tmp1<x_bin(q+1)));
                  index         = tmp2>=y_bin(q) & tmp2<y_bin(q+1);
                  y_mean_bin(q) = nanmean(tmp2(index));
                  events(q)     = sum(index);
                else           
                  x_mean_bin(q) = nanmean(tmp1(tmp1>=x_bin(q) & tmp1<=x_bin(q+1)));
    %               y_mean_bin(q)=nanmean(tmp2(tmp2>=y_bin(q) & tmp2<=y_bin(q+1)));
    %               events(q)=length(tmp1(tmp1>=x_bin(q) & tmp1<=x_bin(q+1)));
                  index         = tmp2>=y_bin(q) & tmp2<=y_bin(q+1); 
                  y_mean_bin(q) = nanmean(tmp2(index));
                  events(q)     = sum(index);
                end
            end

    %     end
    % P_AIB(counter_A>=min_events,counter_B>=min_events)

        binP_AIB = NaN(num_bin,num_bin);
    %     binP_AIB_2_4=NaN(aufl,aufl);
    for s=1:num_bin
        for q=1:num_bin
            binP_AIB(q,s)       = x_mean_bin(q)-y_mean_bin(s);
    %         binP_AIB_2_4(q,s)= abs(x_mean_bin(q)-y_mean_bin(s));
        end
    end
    else
    binP_AIB    = nan; 
    x_mean_bin  = nan; 
    y_mean_bin  = nan;
    events      = nan;    
end
end