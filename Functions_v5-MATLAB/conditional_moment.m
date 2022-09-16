%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function estimates the $k$-order conditional moment M^{(k)}\left(u_r,r,\Delta r\right)
% with $k={1,2,3,4}$ for all scales $2\Delta_{EM} <r \leq L$ (specified in the function 
% \textbf{\textit{scale\_steps}}) and for each bin (specified in the function 
% \textbf{\textit{increment\_bin}}) for all values of longitudinal velocity increments $u_r$. For a 
% fixed scale $r$ the conditional moments are calculated for 5 different scales separations 
% $\Delta r=r-r'$ within the range of $\Delta_{EM}\leq \Delta r\leq 2\Delta_{EM}$. 
% The condition $r'<r$ is always fulfilled.
%
% Arguments IN
% low_freq = Frequency at which you would like to use the low-pass filter in Hz
% Fs = Acquisition/Sampling Frequency in Hz
% markov = markov length in number of samples
% m_data = mean of the data
% int_L = Integral length scale in meters
% taylor_L = Taylor length scale in meters
% scale_steps = Number of seperated scales between Integral and Taylor length
% increment_bin = number of bins
% data = filtered data 
% point = Multipoint condition 1=YES or 2=NO
% condition = This input is for multipoint statistics
% tol = This input is for multipoint statistics
%
% Arguments OUT
% evaluated = struct array containing all the information about conditional moments for each scale and for each bin
% step_con_moment = the steps at which the value of conditional moments are calculated
% Enclosed in a "evaluated":
% r is the scale in meters at which moments will be calculated and hence this r will
% be the same at which D1 & D2 will be calculated===>r2
% r_samp is the r in number of samples==>r2
% r_short_sample is r1 ===> (r2>r1) 
% M11 = First order conditional moment
% M21 = Second order conditional moment
% M31 = Third order conditional moment
% M41 = Fourth order conditional moment
% eM1 = Error associated with M11
% eM2 = Error associated with M21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [evaluated, step_con_moment] = conditional_moment(low_freq,Fs,markov,m_data,int_L,taylor_L,scale_steps,increment_bin,data,multi_point,condition,tol)
%steps (in samples) for which we will calculate the conditional moments
step_con_moment     =round(linspace(markov,2*markov,5));

% the low-pass filter frequency or the smallest still resolved frequency in the spectrum corresponds 
% to the smallest scale for which D1 and D2 can be determined
smallest_sample = round((1/low_freq)*Fs)+step_con_moment(end);
smallest_r      = (1/low_freq)*m_data+(step_con_moment(end)*m_data/Fs);
% smallest_r/taylor_L

% log: scales (in meter) which are used for the estimation of Kramers-Moyal coefficients (D1 and D2).
if taylor_L>smallest_r
    r               = unique(round(logspace(log10(taylor_L*10^6),log10(int_L*10^6),scale_steps))./10^6,'stable');
else
    r               = unique(round(logspace(log10(smallest_r*10^6),log10(int_L*10^6),scale_steps))./10^6,'stable');
end

%lin: scales (in meter) which are used for the estimation of Kramers-Moyal coefficients (D1 and D2).
% if taylor_L>smallest_r
%     r               = unique(linspace(taylor_L,int_L,scale_steps),'stable');
% else
%     r               = unique(linspace(smallest_r,int_L,scale_steps),'stable');
% end


%scales (in time and samples) which are used for the estimation of Kramers-Moyal coefficients (D1 and D2).
clear tau step
for l=1:length(r)
    tau(l)  = r(l)/m_data;
    step(l) = round(tau(l)*Fs); 
end
%scales (in samples) which are used for the estimation of Kramers-Moyal coefficients (D1 and D2).
tau2_inter  = unique(step,'stable');
tau2_inter  = tau2_inter(tau2_inter>step_con_moment(end)); 

r           = (tau2_inter./Fs).*m_data;   % scales in meters


    if multi_point==1
        for k=1:length(condition)
            my_field = sprintf('point_eval_%d', k);
            evaluated.(my_field)(1:length(tau2_inter)) = struct('r',nan,'r_samp',nan,'r_short_sample',nan,'P_AIB',NaN(increment_bin,increment_bin),'counter_A',nan,'counter_B',nan,...
                'x_mean_bin',NaN(1,increment_bin),'y_mean_bin',NaN(1,increment_bin),'x',NaN(1,increment_bin),'y',NaN(1,increment_bin),...
                'x_bin_not_nan',NaN(1,increment_bin),'D1',NaN(1,increment_bin),'eD1',NaN(1,increment_bin),'D2',NaN(1,increment_bin),'eD2',NaN(1,increment_bin),...
                'D3',NaN(1,increment_bin),'D4',NaN(1,increment_bin),...
                'D1_opti',NaN(1,increment_bin),'D2_opti',NaN(1,increment_bin),...
                'M11',NaN(increment_bin,length(step_con_moment)),'M21',NaN(increment_bin,length(step_con_moment)),'M31',NaN(increment_bin,length(step_con_moment)),'M41',NaN(increment_bin,length(step_con_moment)),...
                'eM1',NaN(increment_bin,length(step_con_moment)),'eM2',NaN(increment_bin,length(step_con_moment)),'x0',NaN(1,2*increment_bin),'x1',NaN(1,2*increment_bin));
        end
    else
        evaluated(1:length(tau2_inter)) = struct('r',nan,'r_samp',nan,'r_short_sample',nan,'P_AIB',NaN(increment_bin,increment_bin),'counter_A',nan,'counter_B',nan,...
            'x_mean_bin',NaN(1,increment_bin),'y_mean_bin',NaN(1,increment_bin),'x',NaN(1,increment_bin),'y',NaN(1,increment_bin),...
            'x_bin_not_nan',NaN(1,increment_bin),'D1',NaN(1,increment_bin),'eD1',NaN(1,increment_bin),'D2',NaN(1,increment_bin),'eD2',NaN(1,increment_bin),...
            'D3',NaN(1,increment_bin),'D4',NaN(1,increment_bin),...
            'D1_opti',NaN(1,increment_bin),'D2_opti',NaN(1,increment_bin),...
            'M11',NaN(increment_bin,length(step_con_moment)),'M21',NaN(increment_bin,length(step_con_moment)),'M31',NaN(increment_bin,length(step_con_moment)),'M41',NaN(increment_bin,length(step_con_moment)),...
            'eM1',NaN(increment_bin,length(step_con_moment)),'eM2',NaN(increment_bin,length(step_con_moment)),'x0',NaN(1,2*increment_bin),'x1',NaN(1,2*increment_bin)); 
    end


    if multi_point==1
        for k=1:length(condition)
            my_field = sprintf('point_eval_%d', k);
            %evaluated.point_eval_1(1).x_mean_bin
            for scal = 1:length(tau2_inter) 
               tau2         = tau2_inter(scal);
               bin_events   = NaN(increment_bin,length(step_con_moment));
               M11          = nan; 
               M21          = nan;
               M31          = nan;
               M41          = nan;
               eM1          = nan;
               eM2          = nan;

               for jj=length(step_con_moment):-1:1
                   if multi_point==1
                        [incr1,incr2]  = Increment_point((tau2 - step_con_moment(jj)),tau2,data,condition(k),tol);
                   elseif multi_point==0
                        [incr1,incr2]  = Increment((tau2 - step_con_moment(jj)),tau2,data);
                   end

                   [dx,dy,x,y,dA] = limit(incr1,incr2,increment_bin);
                   if length(x)==increment_bin & length(y)==increment_bin
                       %     [M11_direct(:,jj),M21_direct(:,jj),M41_direct(:,jj),x_mean_bin_direct,y_mean_bin_direct,events] = conditional_moments_direct(tau1,tau2,incr1,incr2,aufl,min_events);
                       %     bin_events_direct(:,jj)=events.';
                       [P_AIB,P_BIA,P_AnB,P_A,P_B,binP_AIB,evaluated.(my_field)(scal).x_mean_bin,evaluated.(my_field)(scal).y_mean_bin,events,counter_A,counter_B] = distribution((tau2 - step_con_moment(jj)),tau2,incr1,incr2,increment_bin,dx,dy,x,y,dA);
                       bin_events(:,jj)=events.';
                       %----------------------------- conditional moments --------------------------------
                       for i=1:increment_bin
                           if (tau2 - step_con_moment(jj))>tau2
                               M11(i,jj)    = nansum(binP_AIB(:,i).^1.*P_AIB(:,i).*dy);
                               M21(i,jj)    = nansum(binP_AIB(:,i).^2.*P_AIB(:,i).*dy);
                               M31(i,jj)    = nansum(binP_AIB(:,i).^3.*P_AIB(:,i).*dy);
                               M41(i,jj)    = nansum(binP_AIB(:,i).^4.*P_AIB(:,i).*dy);
                           else
                               M11(i,jj)    = nansum(binP_AIB(:,i).^1.*P_AIB(:,i).*dx);
                               M21(i,jj)    = nansum(binP_AIB(:,i).^2.*P_AIB(:,i).*dx);
                               M31(i,jj)    = nansum(binP_AIB(:,i).^3.*P_AIB(:,i).*dx);
                               M41(i,jj)    = nansum(binP_AIB(:,i).^4.*P_AIB(:,i).*dx);
                           end
                       end
                       %Fehler von den Momenten
                       eM1= real(sqrt(abs((M21 - (M11.^2))./bin_events)));
                       eM2= real(sqrt(abs((M41 - (M21.^2))./bin_events)));
                   end
               end
                 
               evaluated.(my_field)(scal).counter_A=counter_A;
               evaluated.(my_field)(scal).counter_B=counter_B;
               evaluated.(my_field)(scal).x=x;
               evaluated.(my_field)(scal).y=y;

               evaluated.(my_field)(scal).P_AIB=P_AIB;

               evaluated.(my_field)(scal).r=r(scal);
               evaluated.(my_field)(scal).r_samp=tau2;
               evaluated.(my_field)(scal).r_short_sample=(tau2 - step_con_moment(jj));

               evaluated.(my_field)(scal).M11=M11;
               evaluated.(my_field)(scal).M21=M21;
               evaluated.(my_field)(scal).M31=M31;
               evaluated.(my_field)(scal).M41=M41;
               evaluated.(my_field)(scal).eM1=eM1;
               evaluated.(my_field)(scal).eM2=eM2;
               
               disp(['point:' num2str(k) '/'  num2str(length(condition)) ])
               disp(['scal:' num2str(scal) '/'  num2str(length(tau2_inter)) ])
            end   
        end
    else
        k=nan;
        parfor scal = 1:length(tau2_inter) 
               tau2         = tau2_inter(scal);
               bin_events   = NaN(increment_bin,length(step_con_moment));
               M11          = nan; 
               M21          = nan;
               M31          = nan;
               M41          = nan;
               eM1          = nan;
               eM2          = nan;

               for jj=length(step_con_moment):-1:1
                   if multi_point==1
                        [incr1,incr2]  = Increment_point((tau2 - step_con_moment(jj)),tau2,data,condition(k),tol);
                   elseif multi_point==0
                        [incr1,incr2]  = Increment((tau2 - step_con_moment(jj)),tau2,data);
                   end

                   [dx,dy,x,y,dA] = limit(incr1,incr2,increment_bin);
                   if length(x)==increment_bin & length(y)==increment_bin
                       %     [M11_direct(:,jj),M21_direct(:,jj),M41_direct(:,jj),x_mean_bin_direct,y_mean_bin_direct,events] = conditional_moments_direct(tau1,tau2,incr1,incr2,aufl,min_events);
                       %     bin_events_direct(:,jj)=events.';
                       [P_AIB,P_BIA,P_AnB,P_A,P_B,binP_AIB,evaluated(scal).x_mean_bin,evaluated(scal).y_mean_bin,events,counter_A,counter_B] = distribution((tau2 - step_con_moment(jj)),tau2,incr1,incr2,increment_bin,dx,dy,x,y,dA);
                       bin_events(:,jj)=events.';
                       %----------------------------- conditional moments --------------------------------
                       for i=1:increment_bin
                           if (tau2 - step_con_moment(jj))>tau2
                               M11(i,jj)    = nansum(binP_AIB(:,i).^1.*P_AIB(:,i).*dy);
                               M21(i,jj)    = nansum(binP_AIB(:,i).^2.*P_AIB(:,i).*dy);
                               M31(i,jj)    = nansum(binP_AIB(:,i).^3.*P_AIB(:,i).*dy);
                               M41(i,jj)    = nansum(binP_AIB(:,i).^4.*P_AIB(:,i).*dy);
                           else
                               M11(i,jj)    = nansum(binP_AIB(:,i).^1.*P_AIB(:,i).*dx);
                               M21(i,jj)    = nansum(binP_AIB(:,i).^2.*P_AIB(:,i).*dx);
                               M31(i,jj)    = nansum(binP_AIB(:,i).^3.*P_AIB(:,i).*dx);
                               M41(i,jj)    = nansum(binP_AIB(:,i).^4.*P_AIB(:,i).*dx);
                           end
                       end
                       %Fehler von den Momenten
                       eM1= real(sqrt(abs((M21 - (M11.^2))./bin_events)));
                       eM2= real(sqrt(abs((M41 - (M21.^2))./bin_events)));
                   end
               end

               evaluated(scal).counter_A=counter_A;
               evaluated(scal).counter_B=counter_B;
               evaluated(scal).x=x;
               evaluated(scal).y=y;

               evaluated(scal).P_AIB=P_AIB;

               evaluated(scal).r=r(scal);
               evaluated(scal).r_samp=tau2;
               evaluated(scal).r_short_sample=(tau2 - step_con_moment(jj));

               evaluated(scal).M11=M11;
               evaluated(scal).M21=M21;
               evaluated(scal).M31=M31;
               evaluated(scal).M41=M41;
               evaluated(scal).eM1=eM1;
               evaluated(scal).eM2=eM2;
               disp(['scal:' num2str(scal) '/'  num2str(length(tau2_inter)) ])
        end   
    end
end