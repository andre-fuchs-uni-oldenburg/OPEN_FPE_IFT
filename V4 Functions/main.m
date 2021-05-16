%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    %%%%%%%% V 4.0 %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This command lets you navigate through your local directories in order to select and load the data file(.mat file extension containing hot wire data in unites of m/sec. 
waitfor(msgbox('Select the data file. The name of the variable/time series to be analyzed must be ''data'' (1D column vector).'));
% load_data     = askYesno('Do you want to load .mat file (Yes) other data type (No=Import Tool)', 'Yes');
% if load_data==1
%     uiopen;
% else
    import      = uiimport; 
    data        = import.data;
    clear import
% end
if size(data,2)>1
    data = data.';
end
data_length     = askInput({sprintf('How much %% of data do you want to use for the analysis?')}, {'100'})/100;
data            = data(1:ceil(length(data)*data_length));
ikx             = 21;     
% bar = waitbar(1/ikx,'Loading data','Units','normalized','Position', [0 0.1 0.3 0.1]); pause(.5)
bar             = waitbar(1/ikx,'Loading data');

%% Open dialog box for saving files
waitfor(msgbox('Select a folder for saving figures and files.'));
save_path = uigetdir;

%% Acquisition/Sampling Frequency
Fs      = str2double(inputdlg({'Sampling frequency in Hz:'})); 
m_data  = nanmean(data);

%% Enter the kinematic viscosity of fluid based on your experimental conditions
kin_vis = askInput({sprintf('Enter the kinematic viscosity in m^2/s')}, {'15.32e-6'});

%% Estimation of the number of bin's
% increment_bin = round(sqrt(size(data,1))); % take square root of number of data-points as choice for number of bins
increment_bin = ceil((range(data)/nanstd(data)*10));
if mod(increment_bin,2) == 0
    increment_bin = increment_bin +1;
end
increment_bin = askInput({sprintf('Number of Bin''s used:')}, {num2str(increment_bin,'%1.0f')});

%% plot_stationarity is the function to plot the mean, std, skewness and kurtosis of sections of a length of 5% of the data  
waitbar(2/ikx,bar,'Check: Stationarity');
[data] = plot_stationarity(data,save_path);

%% plot_pdf is the function to plot probability density function(PDF) of the data - For more details open this function script
waitbar(3/ikx,bar,'Plot: PDF');
plot_pdf(data,increment_bin,save_path) 

%% Spectrum of raw data - For more details open this function script
waitbar(4/ikx,bar,'Plot: Spectrum of raw data');
[f, E_f_no_filter, f_avg_no_filter, E_avg_no_filter, P_avg_no_filter] = spectrum(data,Fs,increment_bin);

%% lowpass filter of data
filter     = askYesno('Do you want to filter (Low-pass filter) the data?', 'Yes');
if filter==1
    waitfor(msgbox('Use interactive data cursor to select a cut-off Frequency.'));
    uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
    uiwait(gcf);
    uicontrol('Position',[10 10 100 40],'String','Processing','Callback','uiresume(gcbf)');
    %Choose the frequency from figure(1): ESD at which you want to use low pass filter
    low_freq = str2double(inputdlg({'Cut-off Frequency (Low-pass filter) in Hz (less than Fs/2):'}));
    close
    waitbar(5/ikx,bar,'Lowpass filter of data');
    % Filter the actual data in order to elliminate unwanted experimental noise at higher frequency
    [data_filter,f_avg_filter,E_f_avg_filter,k_avg_filter,Ek_avg_filter,r_avg_filter,Dr_avg_filter,Dk_avg_filter] = frequency_filter(data,Fs,low_freq,kin_vis,filter,save_path); 
else
	data_filter = data;
    low_freq    = Fs/2;
    close
    [data_filter,f_avg_filter,E_f_avg_filter,k_avg_filter,Ek_avg_filter,r_avg_filter,Dr_avg_filter,Dk_avg_filter] = frequency_filter(data,Fs,low_freq,kin_vis,filter,save_path); 
end
if ischar(save_path)
    tmp_bar=bar.Position;
    close(bar)
    save(fullfile(save_path,'evaluated.mat'),'-v7.3') 
    bar = waitbar(6/ikx,'Estimation of fundamental length scales','Position', tmp_bar);
end

%% Estimation of fundamental length scales
lenght_scales     = askYesno('Do you want to estimate the fundamental length scales?', 'Yes');
if lenght_scales==1
    C2 = askInput({sprintf('Enter the value of Kolmogorov constant, C_2 associated with second order structure function(in between 2.0 to 2.2)')}, {'2.1'});
    [int_L,taylor_L,int_L_calc, taylor_L_calc,epsi,epsi_calc,diss_scale,Ce, Ce_calc, Re, Re_lambda ] = length_scales(data_filter,Fs,low_freq,kin_vis,C2,increment_bin,m_data,save_path);
else
    int_L       = str2double(inputdlg({'Enter the integral length scale in meter:'}));
    taylor_L    = str2double(inputdlg({'Enter the Taylor length scale in meter:'})); 
    diss_scale  = str2double(inputdlg({'Enter the Kolmogorov length scale in meter:'})); 
end

%round to the next larger value in samples
diss_scale  = ceil(Fs*diss_scale/m_data)*m_data/Fs;
taylor_L    = ceil(Fs*taylor_L/m_data)*m_data/Fs;
int_L       = ceil(Fs*int_L/m_data)*m_data/Fs;

if ischar(save_path)
    tmp_bar=bar.Position;
    close(bar)
    save(fullfile(save_path,'evaluated.mat'),'-v7.3')  
    bar = waitbar(7/ikx,'Normalization of the data','Position', tmp_bar);
end

%% Weather to flip the Hot-wire data or not?
flip_data       = askYesno('Do you want to check if you have to flip the data?', 'Yes');
norm_ur         = askYesno('Do you want to perform the normalization of the data using $\sigma_\infty$?', 'Yes');
norm_r          = askYesno('Do you want to perform the normalization of the scale using $\lambda$?', 'Yes');
if flip_data==1
    Struc_flip_test(Fs,data_filter,int_L,taylor_L,save_path)
    % Normalize the timeseries by sigma_infinity
    [data_filter,siginf,m_data] = normalization(data_filter,askYesno('Do you want to flip hotwire data?', 'Yes'),norm_ur);
    close
else
    [data_filter,siginf,m_data] = normalization(data_filter,0,norm_ur);
end
if ischar(save_path)
    tmp_bar=bar.Position;
    close(bar)
    save(fullfile(save_path,'evaluated.mat'),'-v7.3')  
    bar = waitbar(8/ikx,'Plot: Increment PDF','Position', tmp_bar);
end

%% Is the function to plot probability density function(PDF) for the velocity increment at integral and Taylor length scale  
plot_increment_pdf(data_filter,301,save_path,int_L,taylor_L,m_data,Fs,diss_scale) 

%% Calculation of structure functions
[r_struc,S_exp_2,S_exp_3,S_exp_4,S_exp_5,S_exp_6,S_exp_7]=plot_struc_function(askInput({sprintf('Number of seperated scales between Integral and Taylor (calculation of structure functions):')}, {num2str(30,'%1.0f')}),Fs,data_filter,m_data,int_L,taylor_L,askInput({sprintf('intermittency coefficient mu:')}, {num2str(0.26,'%1.2f')}),askInput({sprintf('Beta Model Coefficient D :')}, {num2str(2.8,'%1.2f')}),norm_r,save_path);

%% Markov test
% Following function calculates the Markov length or Einstein-Markov length and generates a plot
waitbar(9/ikx,bar,'Wilcoxon test');
markov=wilcoxon_test(data_filter,Fs,m_data,int_L,taylor_L,diss_scale,increment_bin,save_path,f_avg_filter, E_f_avg_filter);

% One must fix the mininum number of events in a single bin which will make that specific bin valid
min_events=ceil(0.00001*length(data));
if min_events<400
    min_events      = 400;
end
min_events = askInput({sprintf('Minimum number of events in a single bin:')}, {num2str(min_events ,'%1.0f')});

conditional_PDF_markov(data_filter,Fs,m_data,markov,askInput({sprintf('Increment of the second condition:')}, {'0'}),increment_bin,min_events,norm_ur,save_path)
% askInput({sprintf('Number of Bin''s used for the increment:')}, {num2str(increment_bin,'%1.0f')}),min_events,save_path)
var_markov    = askYesno('Do you want to modify the markov lenght and/or the minimum number of events?', 'Yes');
while var_markov==1
    markov = askInput({sprintf('Markov length in samples:')}, {num2str(markov  ,'%1.0f')});
    min_events = askInput({sprintf('Mininum number of events in a single bin:')}, {num2str(min_events ,'%1.0f')});
    conditional_PDF_markov(data_filter,Fs,m_data,markov,askInput({sprintf('Increment of the second condition:')}, {'0'}),increment_bin,min_events,norm_ur,save_path)
    var_markov    = askYesno('Do you want to modify the markov lenght and/or the minimum number of events?', 'Yes');
end
if ischar(save_path)
    tmp_bar=bar.Position;
    close(bar)
    save(fullfile(save_path,'evaluated.mat'),'-v7.3') 
    bar = waitbar(10/ikx,'Calculation of conditional moments','Position', tmp_bar);
end

%% Calculation of conditional moments
scale_steps     = askInput({sprintf(['Number of seperated scales between Integral and Taylor length must be less than:' num2str(ceil((int_L-taylor_L)/(markov/Fs*m_data)),'%1.0f')])}, {num2str(ceil((int_L-taylor_L)/(markov/Fs*m_data)-1),'%1.0f')});
multi_point     = askYesno('Would you like to perform a multi-point analysis?', 'Yes');
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
    save(fullfile(save_path,'evaluated.mat'),'-v7.3') 
    bar = waitbar(11/ikx,'Plot: conditional moments','Position', tmp_bar);
end

%% Plot conditional moments
plot_moment     = askYesno('Do you want to plot conditional moments?', 'Yes');
if plot_moment==1
    plot_conditional_moment(askInput({sprintf('Plot: number of Scale:')}, {num2str(scale_steps)}),askInput({sprintf('Plot: number of a bin (value of the velocity increment u_r):')}, {num2str(ceil(increment_bin/2))}),evaluated,step_con_moment,markov,multi_point,increment_bin,condition,tol,data_filter,save_path)
end
waitbar(12/ikx,bar,'Calculation of Kramers-Moyal coefficients');

%% Calculation of Kramers-Moyal coefficients (D1 and D2) for each scale and for each bin
[evaluated] = KM_Calculation(increment_bin,min_events,evaluated,step_con_moment,Fs,taylor_L,m_data,multi_point,condition,norm_r);
if ischar(save_path)
    tmp_bar=bar.Position;
    close(bar)
    save(fullfile(save_path,'evaluated.mat'),'-v7.3') 
    bar = waitbar(13/ikx,'Plot: Non-optimized Kramers-Moyal coefficients','Position', tmp_bar);
end

%% Plotting non-optimized Kramers-Moyal coefficients (D1 and D2) 
KM_plot_raw(evaluated,Fs,taylor_L,multi_point,condition,norm_ur,norm_r,save_path)
waitbar(14/ikx,bar,'Optimization of Kramers-Moyal coefficients: conditional PDF''s');

%% Optimization of Kramers-Moyal coefficients (D1 and D2) for each scale and for each bin
test_opti = askYesno('Would you like to consider an example optimization?', 'Yes');
if test_opti==1   
    tol_opti    = askInput({sprintf('Optimization: Tolerance of the range of Kramers-Moyal coefficients in %%')}, {'10'})/100;
    [evaluated] = KM_STP_optimization(evaluated,increment_bin,Fs,markov,data_filter,m_data,taylor_L,test_opti,multi_point,condition,tol,min_events,tol_opti,norm_ur,norm_r,save_path);
    var_opti    = askYesno('Do you want to modify the tolerance of the range of Kramers-Moyal coefficients?', 'Yes');
    close all
    while var_opti==1
        tol_opti    = askInput({sprintf('Optimization: Tolerance of the range of Kramers-Moyal coefficients in %%')}, {'10'})/100;
        [evaluated] = KM_STP_optimization(evaluated,increment_bin,Fs,markov,data_filter,m_data,taylor_L,test_opti,multi_point,condition,tol,min_events,tol_opti,norm_ur,norm_r,save_path);
        var_opti    = askYesno('Do you want to modify the tolerance of the range of Kramers-Moyal coefficients?', 'Yes');
    end   
    close all
    [evaluated] = KM_STP_optimization(evaluated,increment_bin,Fs,markov,data_filter,m_data,taylor_L,0,multi_point,condition,tol,min_events,tol_opti,norm_ur,norm_r,save_path);
else
    tol_opti = askInput({sprintf('Optimization: Tolerance of the range of Kramers-Moyal coefficients in %%')}, {'10'})/100;
    [evaluated] = KM_STP_optimization(evaluated,increment_bin,Fs,markov,data_filter,m_data,taylor_L,0,multi_point,condition,tol,min_events,tol_opti,norm_ur,norm_r,save_path);
end
if ischar(save_path)
    tmp_bar=bar.Position;
    close(bar)
    save(fullfile(save_path,'evaluated.mat'),'-v7.3') 
    bar = waitbar(15/ikx,'Fit: Optimized Kramers-Moyal coefficients','Position', tmp_bar);
end

%% Fitting optimized D1 and D2 and plot
[co_KM_opti,co_KM_opti_no_offset,co_KM_non_opti,fitresult_D1_conf,fitresult_D2_conf] = FIT_KM(evaluated,increment_bin,taylor_L,Fs,m_data,markov,multi_point,condition,norm_ur,norm_r,save_path);
if ischar(save_path)
    tmp_bar=bar.Position;
    close(bar)
    save(fullfile(save_path,'evaluated.mat'),'-v7.3') 
    bar = waitbar(16/ikx,'Plot: Optimized Kramers-Moyal coefficients and Fit','Position', tmp_bar);
end

%% Plotting Kramers-Moyal coefficients (D1 and D2) and Fit
KM_plot(co_KM_opti,evaluated,Fs,taylor_L,int_L,multi_point,condition,norm_ur,norm_r,save_path)
waitbar(17/ikx,bar,'Plot: Coefficients dij');

%% Plotting different Fits
if multi_point==1
    CO_plot(co_KM_opti.(sprintf('point_eval_%d', askInput({sprintf('Plot: Multi-point number of point condition:')}, {num2str(ceil(length(condition)*1/4))}))), ...
            co_KM_opti.(sprintf('point_eval_%d', askInput({sprintf('Plot: Multi-point number of point condition:')}, {num2str(ceil(length(condition)*3/4))}))),...
            evaluated.point_eval_1,taylor_L,int_L,norm_ur,norm_r);
else
    CO_plot(co_KM_opti,co_KM_non_opti,evaluated,taylor_L,int_L,norm_ur,norm_r);
end
waitbar(18/ikx,bar,'Calculation of Entropy: Conditional PDF''s optimized Kramers-Moyal coefficients');

%% Verification:  INTEGRAL FLUCTUATION THEOREM - Calculation of Entropy
waitfor(msgbox('The calculation leading towards the integral fluctuation theorem will be performed.'));
close all
trajec          = askYesno('Would you like to adjust the start/end of the cascade trajectory? ?', 'Yes'); 
z               = askInput({sprintf('The total entropy variation should be calculated for overlapping (z=1) or independent cascade trajectories (z=3):')},{num2str(3)});
% dr_ind        = min([round((1/low_freq)*Fs)*2+1 markov round(Fs*taylor_L/m_data)])
dr_ind          = askInput({sprintf('Enter the separation of scales/step increment (in samples) referred to the sequence from large to small scales in the cascade trajectory:')},{num2str(markov)});
data_length     = askInput({sprintf('How much %% of data do you want do use for the analysis?')}, {'100'})/100;

[Sm,Ds,DS,r,ind_trajec,rind,dr,r_s,~,~,A,~,~,~,~,~,~] = checkFT(evaluated,int_L,taylor_L,Fs,data_filter,m_data,z,co_KM_opti,data_length,markov,trajec,dr_ind,norm_ur,norm_r);

% Plot of Entropy
plot_entropy(Sm,Ds,DS,increment_bin,save_path)

if ischar(save_path)
    tmp_bar=bar.Position;
    close(bar)
    save(fullfile(save_path,'evaluated.mat'),'-v7.3') 
end


%% Optimization of Kramers-Moyal coefficients Fit: IFT opti
close all
opti_IFT    = askYesno('Do you want to perform the pointwise optimization of Kramers-Moyal coefficients towards the integral fluctuation theorem?', 'Yes');
if opti_IFT==1
    bar = waitbar(19/ikx,'Optimization of Kramers-Moyal coefficients: IFT','Position', tmp_bar);
    z               = 3;
    dr_ind          = 1;
    iter            = askInput({sprintf('Enter the maximum number of iteration used for the optimization algorithms:')}, {'5'});
    tol_D1          = askInput({sprintf('Optimization: Tolerance of the range of Kramers-Moyal coefficients in %%')}, {'10'})/100;
    tol_D2          = tol_D1;

    [co_IFT_opti,history]=OPTI_IFT_dij(tol_D1,tol_D2,evaluated,co_KM_opti,fitresult_D1_conf,fitresult_D2_conf,int_L,taylor_L,m_data,Fs,...
        data_filter,z,iter,markov,trajec,dr_ind,norm_ur,norm_r,save_path);
    % [history]       = DS_iter(tol_D1,tol_D2,history,evaluated,Fs,int_L,taylor_L,scale_steps,co_KM_opti,m_data,data_filter,z,markov,trajec,dr_ind,askYesno('co-calculation?', 'Yes'),askYesno('DS-calculation?', 'Yes'));
    % find(history.fval==min(history.fval),1,'first')
%   KM_plot(co_IFT_opti,evaluated,Fs,taylor_L,int_L,multi_point,condition,norm_ur,norm_r,save_path)
%   CO_plot(co_IFT_opti,co_KM_opti,evaluated,taylor_L,int_L,norm_ur,norm_r);
    if ischar(save_path)
        close all
        tmp_bar=bar.Position;
        close(bar)
        save(fullfile(save_path,'evaluated.mat'),'-v7.3') 
        bar = waitbar(20/ikx,'Calculation of Entropy: IFT optimized Kramers-Moyal coefficients','Position', tmp_bar);
    end

    z               = 1;
    data_length     = 1;
    [Sm,Ds,DS,r,ind_trajec,rind,dr,r_s,~,~,A,~,~,~,~,~,~] = checkFT(evaluated,int_L,taylor_L,Fs,data_filter,m_data,z,co_IFT_opti,data_length,markov,trajec,dr_ind,norm_ur,norm_r);
    plot_entropy(Sm,Ds,DS,increment_bin,save_path)
    if ischar(save_path)
        tmp_bar=bar.Position;
        close(bar)
        save(fullfile(save_path,'evaluated.mat'),'-v7.3') 
%         bar = waitbar(21/ikx,'Calculation of structure functions','Position', tmp_bar);
    end
end

%% Verification: Calculation of structure functions (2,3,6)
% [r_struc,S_exp_2,S_exp_3,S_exp_6,S_KM_ST_opti_func_2,S_KM_ST_opti_func_3,S_KM_ST_opti_func_6]=Struc_recon(60,increment_bin,co_IFT_opti,m_data,Fs,data_filter,evaluated,int_L,taylor_L,dr_ind);
% if ischar(save_path)
%     tmp_bar=bar.Position;
%     close(bar)
%     save(fullfile(save_path,'evaluated.mat'),'-v7.3') 
%     bar = waitbar(22/ikx,'Plot: Instantaneous stationary distribution','Position', tmp_bar);
% end

%% Plot the instantaneous stationary distribution of the FPE for a fixed scale r
% [F,D,Z] = pdf_in_st_dist(evaluated,int_L,taylor_L,m_data,Fs,data_filter,ind_trajec,rind,co_KM_opti);
% [F,D,Z] = pdf_in_st_dist(evaluated,int_L,taylor_L,m_data,Fs,data_filter,ind_trajec,rind,co_IFT_opti,increment_bin);

%% The analysis is complete.
close all
delete(bar)
% Save results
if ischar(save_path)  
   save(fullfile(save_path,'evaluated.mat'),'-v7.3') 
end
waitfor(msgbox('Analysis is completed.','Success','help'));