%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% PART II: Markov analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,ikx,bar,tmp_bar,save_path,save_name,Fs,m_data,kin_vis,increment_bin,f, E_f_no_filter, f_avg_no_filter, E_avg_no_filter, P_avg_no_filter,filter,...
    low_freq,data_filter,f_avg_filter,E_f_avg_filter,k_avg_filter,Ek_avg_filter,r_avg_filter,Dr_avg_filter,Dk_avg_filter,...
    lenght_scales,C2,int_L,taylor_L,int_L_calc, taylor_L_calc,epsi,epsi_calc,diss_scale,Ce, Ce_calc, Re, Re_lambda,...
    flip_data,norm_ur,norm_r,siginf,r_struc,S_exp_2,S_exp_3,S_exp_4,S_exp_5,S_exp_6,S_exp_7,...
    markov,min_events,var_markov,scale_steps,multi_point,condition,tol,evaluated,step_con_moment,plot_moment,...
    test_opti,tol_opti,var_opti,co_KM_opti,co_KM_opti_no_offset,co_KM_non_opti,fitresult_D1_conf,fitresult_D2_conf]=PART_II(...
    data,ikx,bar,tmp_bar,save_path,save_name,Fs,m_data,kin_vis,increment_bin,f, E_f_no_filter, f_avg_no_filter, E_avg_no_filter, P_avg_no_filter,filter,...
    low_freq,data_filter,f_avg_filter,E_f_avg_filter,k_avg_filter,Ek_avg_filter,r_avg_filter,Dr_avg_filter,Dk_avg_filter,...
    lenght_scales,C2,int_L,taylor_L,int_L_calc, taylor_L_calc,epsi,epsi_calc,diss_scale,Ce, Ce_calc, Re, Re_lambda,...
    flip_data,norm_ur,norm_r,siginf,r_struc,S_exp_2,S_exp_3,S_exp_4,S_exp_5,S_exp_6,S_exp_7)

bar = waitbar(10/ikx,'Wilcoxon test','Position', tmp_bar);

%% Markov test
% Following function calculates the Markov length or Einstein-Markov length and generates a plot
markov=wilcoxon_test(data_filter,Fs,m_data,int_L,taylor_L,diss_scale,increment_bin,save_path,save_name,f_avg_filter, E_f_avg_filter);

% One must fix the mininum number of events in a single bin which will make that specific bin valid
min_events=ceil(0.00001*length(data));
if min_events<400
    min_events      = 400;
end
min_events = askInput({sprintf('Minimum number of events in a single bin:')}, {num2str(min_events ,'%1.0f')});

% conditional_PDF_markov(data_filter,Fs,m_data,30,askInput({sprintf('Increment of the second condition:')}, {'0'}),101,min_events,norm_ur,save_path)
conditional_PDF_markov(data_filter,Fs,m_data,markov,askInput({sprintf('Increment of the second condition:')}, {'0'}),increment_bin,min_events,norm_ur,save_path,save_name)
% askInput({sprintf('Number of Bin''s used for the increment:')}, {num2str(increment_bin,'%1.0f')}),min_events,save_path)
var_markov    = askYesno('Do you want to modify the markov lenght and/or the minimum number of events?', 'Yes');
while var_markov==1
    markov = askInput({sprintf('Markov length in samples:')}, {num2str(markov  ,'%1.0f')});
    min_events = askInput({sprintf('Mininum number of events in a single bin:')}, {num2str(min_events ,'%1.0f')});
    conditional_PDF_markov(data_filter,Fs,m_data,markov,askInput({sprintf('Increment of the second condition:')}, {'0'}),increment_bin,min_events,norm_ur,save_path,save_name)
    var_markov    = askYesno('Do you want to modify the markov lenght and/or the minimum number of events?', 'Yes');
end
if ischar(save_path)
    tmp_bar=bar.Position;
    close(bar)
    save(fullfile(save_path,append(save_name,'_','evaluated.mat')),'-v7.3') 
    bar = waitbar(11/ikx,'Calculation of conditional moments','Position', tmp_bar);
end

%% Calculation of conditional moments
scale_steps     = askInput({sprintf(['Number of seperated scales between Integral and Taylor length must be less than:' num2str(ceil((int_L-taylor_L)/(markov/Fs*m_data)),'%1.0f')])}, {num2str(ceil((int_L-taylor_L)/(markov/Fs*m_data)-1),'%1.0f')});
multi_point     = askYesno('Would you like to perform a multi-point analysis?', 'No');
if multi_point==1 
    [condition,tol] = multi_point_condition(data_filter);
    min_events      = 50;
else
    condition = nan;
    tol = nan;
end
% Following function estimates the conditional moments for all scales and for each bin i.e. for all values of velocity increment
[evaluated,step_con_moment] = conditional_moment(low_freq,Fs,markov,m_data,int_L,taylor_L,scale_steps,increment_bin,data_filter,multi_point,condition,tol);
if ischar(save_path)
    tmp_bar=bar.Position;
    close(bar)
    save(fullfile(save_path,append(save_name,'_','evaluated.mat')),'-v7.3') 
    bar = waitbar(12/ikx,'Plot: conditional moments','Position', tmp_bar);
end

%% Plot conditional moments
plot_moment     = askYesno('Do you want to plot conditional moments?', 'Yes');
if plot_moment==1
    plot_conditional_moment(askInput({sprintf('Plot: number of Scale:')}, {num2str(scale_steps)}),askInput({sprintf('Plot: number of a bin (value of the velocity increment u_r):')}, {num2str(ceil(increment_bin/2))}),evaluated,step_con_moment,markov,multi_point,increment_bin,condition,tol,data_filter,save_path,save_name);
end
waitbar(13/ikx,bar,'Calculation of Kramers-Moyal coefficients');

%% Calculation of Kramers-Moyal coefficients (D1 and D2) for each scale and for each bin
[evaluated] = KM_Calculation(increment_bin,min_events,evaluated,step_con_moment,Fs,taylor_L,m_data,multi_point,condition,norm_r);
if ischar(save_path)
    tmp_bar=bar.Position;
    close(bar)
    save(fullfile(save_path,append(save_name,'_','evaluated.mat')),'-v7.3') 
    bar = waitbar(14/ikx,'Plot: Non-optimized Kramers-Moyal coefficients','Position', tmp_bar);
end

%% Plotting non-optimized Kramers-Moyal coefficients (D1 and D2) 
KM_plot_raw(evaluated,Fs,taylor_L,multi_point,condition,norm_ur,norm_r,save_path,save_name)
waitbar(15/ikx,bar,'Optimization of Kramers-Moyal coefficients: conditional PDF''s');

%% Optimization of Kramers-Moyal coefficients (D1 and D2) for each scale and for each bin
test_opti = askYesno('Would you like to consider an example optimization?', 'Yes');
if test_opti==1   
    tol_opti    = askInput({sprintf('Optimization: Tolerance of the range of Kramers-Moyal coefficients in %%')}, {'10'})/100;
    [evaluated] = KM_STP_optimization(evaluated,increment_bin,Fs,markov,data_filter,m_data,taylor_L,test_opti,multi_point,condition,tol,min_events,tol_opti,norm_ur,norm_r,save_path,save_name);
    var_opti    = askYesno('Do you want to modify the tolerance of the range of Kramers-Moyal coefficients?', 'Yes');
    close all
    while var_opti==1
        tol_opti    = askInput({sprintf('Optimization: Tolerance of the range of Kramers-Moyal coefficients in %%')}, {'10'})/100;
        [evaluated] = KM_STP_optimization(evaluated,increment_bin,Fs,markov,data_filter,m_data,taylor_L,test_opti,multi_point,condition,tol,min_events,tol_opti,norm_ur,norm_r,save_path,save_name);
        var_opti    = askYesno('Do you want to modify the tolerance of the range of Kramers-Moyal coefficients?', 'Yes');
    end   
    close all
    [evaluated] = KM_STP_optimization(evaluated,increment_bin,Fs,markov,data_filter,m_data,taylor_L,0,multi_point,condition,tol,min_events,tol_opti,norm_ur,norm_r,save_path,save_name);
else
    tol_opti = askInput({sprintf('Optimization: Tolerance of the range of Kramers-Moyal coefficients in %%')}, {'10'})/100;
    [evaluated] = KM_STP_optimization(evaluated,increment_bin,Fs,markov,data_filter,m_data,taylor_L,0,multi_point,condition,tol,min_events,tol_opti,norm_ur,norm_r,save_path,save_name);
end
if ischar(save_path)
    tmp_bar=bar.Position;
    close(bar)
    save(fullfile(save_path,append(save_name,'_','evaluated.mat')),'-v7.3') 
    bar = waitbar(16/ikx,'Fit: Optimized Kramers-Moyal coefficients','Position', tmp_bar);
end

%% Fitting optimized D1 and D2 and plot
[co_KM_opti,co_KM_opti_no_offset,co_KM_non_opti,fitresult_D1_conf,fitresult_D2_conf] = FIT_KM(evaluated,increment_bin,taylor_L,Fs,m_data,markov,multi_point,condition,norm_ur,norm_r,save_path,save_name);
if ischar(save_path)
    tmp_bar=bar.Position;
    close(bar)
    save(fullfile(save_path,append(save_name,'_','evaluated.mat')),'-v7.3') 
    bar = waitbar(17/ikx,'Plot: Optimized Kramers-Moyal coefficients and Fit','Position', tmp_bar);
end

%% Plotting Kramers-Moyal coefficients (D1 and D2) and Fit
KM_plot(co_KM_opti,evaluated,Fs,taylor_L,int_L,multi_point,condition,norm_ur,norm_r,save_path,save_name)
waitbar(18/ikx,bar,'Plot: Coefficients dij');

%% Plotting different Fits
if multi_point==1
    CO_plot(0,co_KM_opti.(sprintf('point_eval_%d', askInput({sprintf('Plot: Multi-point number of point condition:')}, {num2str(ceil(length(condition)*1/4))}))), ...
            co_KM_opti.(sprintf('point_eval_%d', askInput({sprintf('Plot: Multi-point number of point condition:')}, {num2str(ceil(length(condition)*3/4))}))),...
            evaluated.point_eval_1,taylor_L,int_L,norm_ur,norm_r,siginf);
else
    CO_plot(0,co_KM_opti,co_KM_non_opti,evaluated,taylor_L,int_L,norm_ur,norm_r,siginf);
end
end