%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This function is mainly to perform the normalization of the data. Before doing so, it generates
% the pop-up dialog box which asks the user whether to flip the data or not (based on the previous 
% investigation). After that, this function normalizes the entire data with the quantity of 
% $\sigma_{\infty} = \sqrt{2}\sigma$, where $\sigma$ is the standard deviation of the 
% \textbf{data\_filter} (this method is proposed by \cite{renner2001}). This function returns the 
% filtered and normalized data as \textbf{data\_filter}, \mbox{\textbf{siginf} = $\sigma_{\infty}$}
% and \mbox{\textbf{m\_data} = mean of the data} before normalization. In addition scale $r$ is 
% given in units of Taylor length scale $\lambda$. We use this normalization to compare the results 
% of different data sets.
%
% Arguments IN
% data= it is the filtered data
% tmp_flip=Flip hotwire data 1=Yes, 0=No
%
% Arguments OUT
% data= flipped or non-flipped filterd data based on your input as 'tmp_flip'
% siginf=Quantity used to non-dimentionalized the velocity time series
% m_data=mean of the filtered data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [data,siginf,m_data]=normalization(data,tmp_flip,norm_ur)
data    = reshape(data,[],1);
data    = data(~isnan(data));
m_data  = mean(data);
if tmp_flip==1  
    data = flipud(data); 
end

siginf  = sqrt(2).*std(data); % standard variance at infinite scales
if norm_ur==1  
    data    = data./siginf; % dimless flow velocity
end
end