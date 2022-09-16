%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs the pointwise optimization of Kramers-Moyal coefficients 
% $D^{(1,2)}\left(u_r,r\right)$ at each scale and value of velocity increment to minimize possible 
% uncertainties in the absolute values of the Kramers–Moyal coefficients. The object of this 
% optimization is to find the best Fokker-Planck equation to reproduce the conditional PDF's as 
% these are the essential part of the Markov process. This optimization procedure is proposed in 
% \cite{kleinhans2005iterative,Nawroth2007,Reinke2018} and it includes the reconstruction of the 
% conditional probability density functions $p\left(u_{r'} | u_r \right)$ via the short time 
% propagator \cite{Risken} 
% 
% Arguments IN
% evaluated = struct array calculated in the function 'conditional_moment'
% increment_bin = number of bins
% Fs = Acquisition/Sampling Frequency in Hz
% markov = markov length in number of samples
% data = filtered data 
% m_data = mean of the data
% taylor_L = Taylor length scale in meters
% test_opti = weather to plot or not ==> 'Plot? 1=Yes, 0=No'
% multi_point = weather to do multi-point analysis or not 1=Yes, 0=No
% condition = condition for multi-point analysis
% tol = This input is for multipoint statistics
% min_events = minimum number of events
% tol_opti = Optimization: Tolerance of the range of Kramers-Moyal coefficients in % 
%            It is the percentage(For Ex: 0.1 for 10 percent or 0.2 for 20 percent)of D1 or D2 
%            within these limit which you want to optimize these coeffcients D1 & D2
% norm_ur = normalization of the data using $\sigma_\infty$? data 1=Yes, 0=No
% norm_r = normalization of the scale using $\lambda$? data 1=Yes, 0=No
% save_path = path for saving figures and files
% save_name = name for saving files
%
% Arguments OUT
% evaluated = a modified/updated struct 'evaluated' array with optimized value of D1 & D2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [evaluated] = KM_STP_optimization(evaluated,increment_bin,Fs,markov,data,m_data,taylor_L,test_opti,multi_point,condition,tol,min_events,tol_opti,norm_ur,norm_r,save_path,save_name)
if multi_point==1
    if test_opti==1
        k           = askInput({sprintf('Plot: Multi-point number of point condition:')}, {num2str(ceil(length(condition)/2))});
        my_field    = sprintf('point_eval_%d', k);  
        scal        = askInput({sprintf('Which scale number should be optimized?')}, {'1'});
        if sum(evaluated.(my_field)(scal).x_bin_not_nan>0)>0
            [evaluated.(my_field)] = pointwise_KM_STP_optimization(evaluated.(my_field),scal,increment_bin,markov,data,Fs,tol_opti,m_data,taylor_L,test_opti,multi_point,condition,tol,min_events,k,norm_ur,norm_r,save_path,save_name); 
        else
            waitfor(msgbox('There are not enough bins for which the condition events>min_events applies.'));
        end
    else
        for k=1:length(condition)
            my_field = sprintf('point_eval_%d', k);
            for scal = 1:size(evaluated.(my_field),2)
                if sum(evaluated.(my_field)(scal).x_bin_not_nan>0)>0
                    [evaluated.(my_field)]=pointwise_KM_STP_optimization(evaluated.(my_field),scal,increment_bin,markov,data,Fs,tol_opti,m_data,taylor_L,test_opti,multi_point,condition,tol,min_events,k,norm_ur,norm_r);
                end
                disp(['point:' num2str(k) '/'  num2str(length(condition)) ])
                disp(['scal:' num2str(scal) '/'  num2str(size(evaluated.(my_field),2)) ])  
            end
        end      
    end
else
    k=nan;
    %% with x,y_mean bin
    if test_opti==1
        scal        = askInput({sprintf('Which scale number should be optimized?')}, {'1'});  
        [evaluated] = pointwise_KM_STP_optimization(evaluated,scal,increment_bin,markov,data,Fs,tol_opti,m_data,taylor_L,test_opti,multi_point,condition,tol,min_events,k,norm_ur,norm_r,save_path,save_name);   
    else
        f = uifigure;
        d = uiprogressdlg(f,'Title','Please Wait','ShowPercentage','on');
        for scal = 1:size(evaluated,2)
            [evaluated]=pointwise_KM_STP_optimization(evaluated,scal,increment_bin,markov,data,Fs,tol_opti,m_data,taylor_L,test_opti,multi_point,condition,tol,min_events,k,norm_ur,norm_r);
            disp(['scal:' num2str(scal) '/'  num2str(size(evaluated,2)) ])
            d.Value=scal/size(evaluated,2);
        end
        close(d);
        close(f);
    end

    %% x and y same binning 
    % if test_opti==1
    %     scal=str2double(inputdlg({'Which scale number should be optimized?'}))   
    %     [evaluated]=pointwise_KM_STP_optimization_no_mean(evaluated,scal,increment_bin,markov,data,Fs,tol,m_data,taylor_L,test_opti,multi_point,condition,min_events,save_path);   
    % else
    %     for scal = 1:size(evaluated,2)
    %         [evaluated]=pointwise_KM_STP_optimization_no_mean(evaluated,scal,increment_bin,markov,data,Fs,tol,m_data,taylor_L,test_opti,multi_point,condition,min_events);
    %     disp(['scal:' num2str(scal) '/'  num2str(size(evaluated,2)) ])  
    %     end
    % end
end
end   