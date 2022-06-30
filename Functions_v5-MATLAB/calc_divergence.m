%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The aim of this optimization is to minimize a weighted mean square error function in logarithmic 
% space Comment: https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
%
% Arguments IN
% x0 = Initial D1 and D2 before optimization 
% x_mean_bin = 1D array of centers of bins of incr2
% y_mean_bin = 1D array of centers of bins of incr1
% P_B = PDF of B i.e incr2
% P_AnB = Joint PDF of A i.e. incr1 and B i.e incr2
% P_AIB = PDF of A i.e incr1 conditioned on B i.e incr2
% markov = markov length in number of samples 
% x_bin_not_nan = 
% counter_A = histogram of A
% counter_B = histogram of B
% min_events = minimum number of events
% taylor_L = Taylor length scale in meters
% m_data = mean of the data
% Fs = Acquisition/Sampling Frequency in Hz
% tau1 =   is r1 ===> (r2>r1) 
% tau2 =   is the r in number of samples==>is the scale in meters at which moments will be calculated and hence this r will
% be the same at which D1 & D2 will be calculated
%
% Arguments OUT
% d_M = weighted mean square error function in logarithmic space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d_M = calc_divergence(x0,y_mean_bin,x_mean_bin,P_B,P_AnB,P_AIB,markov,x_bin_not_nan,counter_A,counter_B,min_events,taylor_L,m_data,Fs,tau1,tau2,norm_ur,norm_r)
% y_mean_bin(counter_B>=min_events),x_mean_bin(counter_A>=min_events),P_AIB(counter_A>=min_events,counter_B>=min_events)
% tmp=P_AIB(counter_A>=min_events,counter_B>=min_events);
y_input     = y_mean_bin(counter_B>=min_events);
x_input     = x_mean_bin(counter_A>=min_events);
% y_input=y_mean_bin;
% x_input=x_mean_bin;

D1_poly     = x0(1:sum(x_bin_not_nan));
% D2_poly     = x0(length(x_bin_not_nan)+1:end);
D2_poly     = abs(x0(sum(x_bin_not_nan)+1:end));

%% Individual adjustment of the coefficient
% P1_new(1:length(P1)) = P1;
% P2_new(1:length(P2)) = P2;

% P1_new(1)   = x0(1);
% P1_new(2)   = x0(2);
% P1_new(3)   = x0(3);
% P1_new(4)   = x0(4);
% 
% P2_new(1)   = x0(5);
% P2_new(2)   = x0(6);
% P2_new(3)   = x0(7);

% D1_poly     = y.^3.*P1_new(1) +  y.^2.*P1_new(2) + y.*P1_new(3) + P1_new(4);
% D2_poly     = y.^2.*P2_new(1) +  y.*P2_new(2)    + P2_new(3);
% using y
% D1_poly= polyval(coeffvalues(fitresultf1), y);
% D1_poly= polyval(P1_new, y);
% D2_poly= polyval(P2_new, y);


%% short time propagator \cite{Risken} 
% tau=markov;
P_n = ShortTimeProp(y_input,x_input,tau1,tau2,D1_poly,D2_poly,taylor_L,m_data,Fs,norm_ur,norm_r); 


%% weighted mean square error function in logarithmic space; conditional probability densities
P_AIB   = P_AIB(counter_A>=min_events,counter_B>=min_events);
d_M_1   = 0;
d_M_2   = 0;
for j =1:length(y_input)
    for i =1:length(x_input)
         if (P_n(i,j) > 0) && (P_AIB(i,j) > 0)     
            d_M_1 = d_M_1 + ((P_n(i,j) + P_AIB(i,j)).*...
                            (log(P_n(i,j))    - log(P_AIB(i,j))).^2);                        
            d_M_2 = d_M_2 + ((P_n(i,j) + P_AIB(i,j)).*...
                            (log(P_n(i,j)).^2 + log(P_AIB(i,j)).^2));
        end
%         if (P_verbund_n(i,j) < 0)
%             d_M_1 = d_M_1 + 10.^(10);                        
%         end
    end
end
d_M     = d_M_1/d_M_2;


%% Joint probability densities
% P_verbund_n(1:length(x_bin_not_nan),1:length(x_bin_not_nan))=nan;
% for j=1:length(x_bin_not_nan)
% P_verbund_n(:,j) = P_n(:,j)*P_B(1,x_bin_not_nan(j));   %hier wird das P_B bezüglich der experimentellen Daten verwendet
% end
% % % % % % P_verbund_n = P_n.*P_B(x_bin_not_nan);   %hier wird das P_B bezüglich der experimentellen Daten verwendet
% P_n_opti*P_B(evaluated(k).x_bin_not_nan(~isnan(evaluated(k).x_bin_not_nan)))
% P_n_opti_func*P_B
% d_M_1=0;
% d_M_2=0;

% for j =1:(length(x_bin_not_nan))
%     for i =1:(length(x_bin_not_nan))
%          if (P_verbund_n(i,j) > 0) && (P_AnB(x_bin_not_nan(i),x_bin_not_nan(j)) > 0)
%             
%             d_M_1 = d_M_1 + (P_verbund_n(i,j) + P_AnB(x_bin_not_nan(i),x_bin_not_nan(j))).*...
%                             (log(P_verbund_n(i,j))    - log(P_AnB(x_bin_not_nan(i),x_bin_not_nan(j)))).^2;                        
%             d_M_2 = d_M_2 + (P_verbund_n(i,j) + P_AnB(x_bin_not_nan(i),x_bin_not_nan(j))).*...
%                             (log(P_verbund_n(i,j)).^2 + log(P_AnB(x_bin_not_nan(i),x_bin_not_nan(j))).^2);
%         end
% %         if (P_verbund_n(i,j) < 0)
% %             d_M_1 = d_M_1 + 10.^(10);                        
% %         end
%     end
% end
% d_M     = d_M_1/d_M_2;
% % d_M     = d_M_1/(d_M_2+10.^(-10));


%% Hellinger
% d_M_1=0;
% for j =1:length(y_input)
%     for i =1:length(x_input)
%         if (P_n(i,j) > 0) && (P_AIB(i,j) > 0)
%             d_M_1 = d_M_1 + 0.5*(sqrt(P_n(i,j))-sqrt(P_AIB(i,j))).^2;          
%         end
%     end
% end
% d_M     = d_M_1;
% These show no good results when described with the IFT
% probably because too little weight is applied to the extreme edges.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test with conditional probability densities with weighting
% d_M_1=0;
% d_M_2=0;
% for j =1:(length(x_bin_not_nan))
%     for i =1:(length(x_bin_not_nan))
%          if (P_n(i,j) > 0) && (P_AIB(x_bin_not_nan(i),x_bin_not_nan(j)) > 0)
%             
%             d_M_1 = d_M_1 + ((P_n(i,j) + P_AIB(x_bin_not_nan(i),x_bin_not_nan(j))).*...
%                             (log(P_n(i,j))    - log(P_AIB(x_bin_not_nan(i),x_bin_not_nan(j)))).^2).*sqrt(counter2(x_bin_not_nan(i),x_bin_not_nan(j)));                        
%             d_M_2 = d_M_2 + (P_n(i,j) + P_AIB(x_bin_not_nan(i),x_bin_not_nan(j))).*...
%                             (log(P_n(i,j)).^2 + log(P_AIB(x_bin_not_nan(i),x_bin_not_nan(j))).^2);
%         end
% %         if (P_verbund_n(i,j) < 0)
% %             d_M_1 = d_M_1 + 10.^(10);                        
% %         end
%     end
% end
% d_M     = d_M_1/d_M_2;


%% RMS
% d_M_1=0;
% k_tmp=1;
% for j =1:(length(x_bin_not_nan))
%     for i =1:(length(x_bin_not_nan))
%          if (P_n(i,j) > 0) && (P_AIB(x_bin_not_nan(i),x_bin_not_nan(j)) > 0)              
%               d_M_1 = d_M_1 + (P_n(i,j) - P_AIB(x_bin_not_nan(i),x_bin_not_nan(j)))^2;              
%               k_tmp=k_tmp+1;
%          end
%     end
% end
% d_M     = sqrt(d_M_1/k_tmp);


%% simply deviation
% d_M_1=0;
% for j =1:length(y_input)
%     for i =1:length(x_input)
%          if (P_n(i,j) > 0) && (P_AIB(i,j) > 0)         
% %             d_M_1 = d_M_1 + ((P_n(i,j) - P_AIB(x_bin_not_nan(i),x_bin_not_nan(j)))^2/P_AIB(x_bin_not_nan(i),x_bin_not_nan(j)));   
%               d_M_1 = d_M_1 + abs((P_n(i,j) - P_AIB(i,j))); 
%          end
%     end
% end
% d_M     = d_M_1;


%% Kullback-Leibler-Divergenz 
% d_M_1=0;
% for j =1:length(y_input)
%     for i =1:length(x_input)
%          if (P_n(i,j) > 0) && (P_AIB(i,j) > 0)   
%            
%             d_M_1 = d_M_1 + (P_n(i,j).*log(P_n(i,j)./P_AIB(i,j)));
%             
%          end
%     end
% end
% d_M     = d_M_1;


%% Jensen–Shannon-divergence 
% d_M_1=0;
% 
% for j =1:(length(x_bin_not_nan))
%     for i =1:(length(x_bin_not_nan))
%          if (P_n(i,j) > 0) && (P_AIB(x_bin_not_nan(i),x_bin_not_nan(j)) > 0)
%              
%              M=0.5*(P_n(i,j)+P_AIB(x_bin_not_nan(i),x_bin_not_nan(j)));
%              
%              d_M_1 = d_M_1 + (0.5*P_n(i,j)*log(P_n(i,j)/M)+...
%                               0.5*P_AIB(x_bin_not_nan(i),x_bin_not_nan(j))*log(P_AIB(x_bin_not_nan(i),x_bin_not_nan(j))/M));            
%             M=0;
%         end
%     end
% end
% d_M     = d_M_1;


%% Jensen–Shannon-divergence 
% cd('C:\Users\André Fuchs\Dropbox\Matlab_Skripte\mathworks\JSDiv\home\nrazavi\Desktop');
% P=reshape(P_AnB_2,1,[]);
% Q=reshape(P_AnB,1,[]);
% dist=JSDiv(P,Q)
% cd('C:\Users\André Fuchs\Dropbox\Matlab_Skripte\Markov\MIT Optimierung\D12_momente');
% 
% d_M_1=0;
% for j =1:aufl
%     for i =1:aufl
%          if (P_AnB_2(i,j) > 0) && (P_AnB(i,j) > 0)
%             
% %             d_M_1 = d_M_1 +((P_AnB(i,j).*log2(P_AnB(i,j)./(P_AnB(i,j)+P_AnB_2(i,j)))...
% %                           +  P_AnB_2(i,j).*log2(P_AnB_2(i,j)./(P_AnB(i,j)+P_AnB_2(i,j)))).*0.5);
%             d_M_1 = d_M_1 +((P_AnB(i,j).*log2(P_AnB(i,j)./(0.5.*(P_AnB(i,j)+P_AnB_2(i,j))))...
%                           +  P_AnB_2(i,j).*log2(P_AnB_2(i,j)./(0.5.*(P_AnB(i,j)+P_AnB_2(i,j))))).*0.5);
% 
%          end
%     end
% end
% d_M(5,q)    = d_M_1;


%% Jensen-Shannon divergence between two distributions
% cd('C:\Users\André Fuchs\Dropbox\Matlab_Skripte\mathworks\JSDiv\home\nrazavi\Desktop');
% P=reshape(P_verbund_n,1,[]);
% Q=reshape(P_AnB(x_bin_not_nan,x_bin_not_nan),1,[]);
% d_M=JSDiv(P,Q);
% cd('C:\Users\André Fuchs\Dropbox\Matlab_Skripte\Andre\KM_IFT')
% 
% d_M_1=0;
% for j =1:(length(x_bin_not_nan))
%     for i =1:(length(x_bin_not_nan))
%          if (P_verbund_n(i,j) > 0) && (P_AnB(x_bin_not_nan(i),x_bin_not_nan(j)) > 0)
%             
%             d_M_1 = d_M_1 +((P_AnB(x_bin_not_nan(i),x_bin_not_nan(j)).*log(P_AnB(x_bin_not_nan(i),x_bin_not_nan(j))./(P_AnB(x_bin_not_nan(i),x_bin_not_nan(j))+P_verbund_n(i,j)))...
%                           +  P_verbund_n(i,j).*log(P_verbund_n(i,j)./(P_AnB(x_bin_not_nan(i),x_bin_not_nan(j))+P_verbund_n(i,j)))).*0.5);
%          end
%     end
% end
% d_M    = d_M_1;


%% JAN Friedrich 1
% d_M_1=0;
% d_M_2=0;
% for j =1:(length(x_bin_not_nan))
%     for i =1:(length(x_bin_not_nan))
%          if (P_n(i,j) > 0) && (P_AIB(x_bin_not_nan(i),x_bin_not_nan(j)) > 0)
%             
%             d_M_1 = d_M_1 + (P_n(i,j) .* P_AIB(x_bin_not_nan(i),x_bin_not_nan(j)));   
%             
%             d_M_2 = d_M_2 + (P_n(i,j).^2 .* P_AIB(x_bin_not_nan(i),x_bin_not_nan(j)).^2);
%          end
%     end
% end
% d_M     = d_M_1/d_M_2;
end