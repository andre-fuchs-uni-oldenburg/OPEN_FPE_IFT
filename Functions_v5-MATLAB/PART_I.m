%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% PART I: Standard turbulence analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,ikx,bar,tmp_bar,save_path,save_name,Fs,m_data,kin_vis,increment_bin,f, E_f_no_filter, f_avg_no_filter, E_avg_no_filter, P_avg_no_filter,filter,...
    low_freq,data_filter,f_avg_filter,E_f_avg_filter,k_avg_filter,Ek_avg_filter,r_avg_filter,Dr_avg_filter,Dk_avg_filter,...
    lenght_scales,C2,int_L,taylor_L,int_L_calc, taylor_L_calc,epsi,epsi_calc,diss_scale,Ce, Ce_calc, Re, Re_lambda,...
    flip_data,norm_ur,norm_r,siginf,r_struc,S_exp_2,S_exp_3,S_exp_4,S_exp_5,S_exp_6,S_exp_7]=PART_I
%% This command lets you navigate through your local directories in order to select and load the data file(.mat file extension containing hot wire data in unites of m/sec. 
waitfor(msgbox('Select the data file. The name of the variable/time series to be analyzed must be ''data'' (1D column vector).'));
% load_data     = askYesno('Do you want to load .mat file (Yes) other data type (No=Import Tool)', 'Yes');
% if load_data==1
%     uiopen;
% else
import  = uiimport; 
import  = cell2struct(struct2cell(import), {'data'});
data    = import.data;
clear import
% end
if size(data,2)>1
    data = data.';
end
data_length     = askInput({sprintf('How much %% of data do you want to use for the analysis?')}, {'100'})/100;
data            = data(1:ceil(length(data)*data_length));
ikx             = 22;     
% bar = waitbar(1/ikx,'Loading data','Units','normalized','Position', [0 0.1 0.3 0.1]); pause(.5)
bar             = waitbar(1/ikx,'Loading data');

%% Open dialog box for saving files
waitfor(msgbox('Select a folder for saving figures and files.'));
save_path   = uigetdir;
if save_path==0
	save_name = save_path;
else
	save_name =  char(inputdlg({'Enter the name for saving files:'}));
end

%% Acquisition/Sampling Frequency
Fs      = str2double(inputdlg({'Sampling frequency in Hz:'})); 
m_data  = nanmean(data);
if ~m_data>0
	m_data = askInput({sprintf('The mean value is smaller than 0 m/s. Should the mean value be set equal to 1 m/s (see readme)?')}, {'1'});
end

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
[data] = plot_stationarity(data,save_path,save_name);

%% plot_pdf is the function to plot probability density function(PDF) of the data - For more details open this function script
waitbar(3/ikx,bar,'Plot: PDF');
plot_pdf(data,increment_bin,save_path,save_name) 

%% Spectrum of raw data - For more details open this function script
waitbar(4/ikx,bar,'Plot: Spectrum of raw data');
[f, E_f_no_filter, f_avg_no_filter, E_avg_no_filter, P_avg_no_filter,K41_index] = spectrum(data,Fs,increment_bin);

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
    [data_filter,f_avg_filter,E_f_avg_filter,k_avg_filter,Ek_avg_filter,r_avg_filter,Dr_avg_filter,Dk_avg_filter] = frequency_filter(data,Fs,low_freq,K41_index,kin_vis,filter,save_path,save_name); 
else
	data_filter = data;
    low_freq    = Fs/2;
    close
    [data_filter,f_avg_filter,E_f_avg_filter,k_avg_filter,Ek_avg_filter,r_avg_filter,Dr_avg_filter,Dk_avg_filter] = frequency_filter(data,Fs,low_freq,K41_index,kin_vis,filter,save_path,save_name); 
end
if ischar(save_path)
    tmp_bar=bar.Position;
    close(bar)
    save(fullfile(save_path,append(save_name,'_','evaluated.mat')),'-v7.3') 
    bar = waitbar(6/ikx,'Estimation of fundamental length scales','Position', tmp_bar);
end

%% Estimation of fundamental length scales
lenght_scales     = askYesno('Do you want to estimate the fundamental length scales?', 'Yes');
if lenght_scales==1
    C2 = askInput({sprintf('Enter the value of constant, C_2 associated with second order structure function:')}, {'2.0'});
    [int_L,taylor_L,int_L_calc, taylor_L_calc,epsi,epsi_calc,diss_scale,Ce, Ce_calc, Re, Re_lambda ] = length_scales(data_filter,Fs,low_freq,kin_vis,C2,increment_bin,m_data,save_path,save_name);
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
    save(fullfile(save_path,append(save_name,'_','evaluated.mat')),'-v7.3')  
    bar = waitbar(7/ikx,'Normalization of the data','Position', tmp_bar);
end

%% Wether to flip the Hot-wire data or not?
flip_data       = askYesno('Do you want to check if you have to flip the data?', 'Yes');
norm_ur        = askYesno('Do you want to perform the normalization of the data using $\sigma_\infty$?', 'Yes');
norm_r          = askYesno('Do you want to perform the normalization of the scale using $\lambda$?', 'Yes');
if flip_data==1
    Struc_flip_test(Fs,data_filter,int_L,taylor_L,save_path,save_name)
    % Normalize the timeseries by sigma_infinity
    [data_filter,siginf,m_data] = normalization(data_filter,askYesno('Do you want to flip hotwire data?', 'Yes'),norm_ur);
    close
else
    [data_filter,siginf,m_data] = normalization(data_filter,0,norm_ur);
end
if ischar(save_path)
    tmp_bar=bar.Position;
    close(bar)
    save(fullfile(save_path,append(save_name,'_','evaluated.mat')),'-v7.3')  
    bar = waitbar(8/ikx,'Plot: Increment PDF','Position', tmp_bar);
end

%% Is the function to plot probability density function(PDF) for the velocity increment at integral and Taylor length scale  
plot_increment_pdf(data_filter,301,save_path,save_name,int_L,taylor_L,m_data,Fs,norm_ur,diss_scale) 

%% Calculation of structure functions
waitbar(9/ikx,bar,'Calculation of structure functions');
[r_struc,S_exp_2,S_exp_3,S_exp_4,S_exp_5,S_exp_6,S_exp_7]=plot_struc_function(askInput({sprintf('Number of seperated scales between Integral and Taylor (calculation of structure functions):')}, {num2str(30,'%1.0f')}),Fs,data_filter,m_data,int_L,taylor_L,askInput({sprintf('intermittency coefficient mu:')}, {num2str(0.26,'%1.2f')}),askInput({sprintf('Beta Model Coefficient D :')}, {num2str(2.8,'%1.2f')}),norm_r,save_path,save_name);
if ischar(save_path)
    tmp_bar=bar.Position;
    close(bar)
    save(fullfile(save_path,append(save_name,'_','evaluated.mat')),'-v7.3')
end
end