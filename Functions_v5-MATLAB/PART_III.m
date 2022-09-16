%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% PART III: Entropy analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,ikx,bar,tmp_bar,save_path,save_name,Fs,m_data,kin_vis,increment_bin,f, E_f_no_filter, f_avg_no_filter, E_avg_no_filter, P_avg_no_filter,filter,...
    low_freq,data_filter,f_avg_filter,E_f_avg_filter,k_avg_filter,Ek_avg_filter,r_avg_filter,Dr_avg_filter,Dk_avg_filter,...
    lenght_scales,C2,int_L,taylor_L,int_L_calc, taylor_L_calc,epsi,epsi_calc,diss_scale,Ce, Ce_calc, Re, Re_lambda,...
    flip_data,norm_ur,norm_r,siginf,r_struc,S_exp_2,S_exp_3,S_exp_4,S_exp_5,S_exp_6,S_exp_7,...
    markov,min_events,var_markov,scale_steps,multi_point,condition,tol,evaluated,step_con_moment,plot_moment,...
    test_opti,tol_opti,var_opti,co_KM_opti,co_KM_opti_no_offset,co_KM_non_opti,fitresult_D1_conf,fitresult_D2_conf,...
    trajec,z,dr_ind,data_length,Sm,Ds,DS,r,ind_trajec,rind,dr,r_s,A,iter,tol_D1,tol_D2,co_IFT_opti,history]=PART_III(...
    data,ikx,bar,tmp_bar,save_path,save_name,Fs,m_data,kin_vis,increment_bin,f, E_f_no_filter, f_avg_no_filter, E_avg_no_filter, P_avg_no_filter,filter,...
    low_freq,data_filter,f_avg_filter,E_f_avg_filter,k_avg_filter,Ek_avg_filter,r_avg_filter,Dr_avg_filter,Dk_avg_filter,...
    lenght_scales,C2,int_L,taylor_L,int_L_calc, taylor_L_calc,epsi,epsi_calc,diss_scale,Ce, Ce_calc, Re, Re_lambda,...
    flip_data,norm_ur,norm_r,siginf,r_struc,S_exp_2,S_exp_3,S_exp_4,S_exp_5,S_exp_6,S_exp_7,...
    markov,min_events,var_markov,scale_steps,multi_point,condition,tol,evaluated,step_con_moment,plot_moment,...
    test_opti,tol_opti,var_opti,co_KM_opti,co_KM_opti_no_offset,co_KM_non_opti,fitresult_D1_conf,fitresult_D2_conf)

bar = waitbar(19/ikx,'Calculation of Entropy: Conditional PDF''s optimized Kramers-Moyal coefficients','Position', tmp_bar);

%% Verification:  INTEGRAL FLUCTUATION THEOREM - Calculation of Entropy
waitfor(msgbox('The calculation leading towards the integral fluctuation theorem will be performed.'));
close all
trajec          = askYesno('Would you like to adjust the start/end of the cascade trajectory? ?', 'Yes');
answer = questdlg('The total entropy variation should be calculated for overlapping or independent cascade trajectories?','Cascade Trajectories', ...
     'overlapping cascade trajectories', ...
     'independent cascade trajectories', ...
     'independent cascade trajectories');
switch answer
    case 'overlapping cascade trajectories'
        z=1;
    case 'independent cascade trajectories'
        z=3;
end
% z               = askInput({sprintf('The total entropy variation should be calculated for overlapping (z=1) or independent cascade trajectories (z=3):')},{num2str(3)});
% dr_ind        = min([round((1/low_freq)*Fs)*2+1 markov round(Fs*taylor_L/m_data)])
dr_ind          = askInput({sprintf('Enter the separation of scales/step increment (in samples) referred to the sequence from large to small scales in the cascade trajectory:')},{num2str(markov)});
data_length     = askInput({sprintf('How much %% of data do you want do use for the analysis?')}, {'100'})/100;

[Sm,Ds,DS,r,ind_trajec,rind,dr,r_s,~,~,A,~,~,~,~,~,~] = checkFT(evaluated,int_L,taylor_L,Fs,data_filter,m_data,z,co_KM_opti,data_length,markov,trajec,dr_ind,norm_ur,norm_r);

% Plot of Entropy
plot_entropy(Sm,Ds,DS,increment_bin,save_path,save_name)

if ischar(save_path)
    tmp_bar=bar.Position;
    close(bar)
    save(fullfile(save_path,append(save_name,'_','evaluated.mat')),'-v7.3') 
end

%% Optimization of Kramers-Moyal coefficients Fit: IFT opti
close all
opti_IFT    = askYesno('Do you want to perform the pointwise optimization of Kramers-Moyal coefficients relatively to the integral fluctuation theorem?', 'Yes');
if opti_IFT==1
    bar = waitbar(20/ikx,'Optimization of Kramers-Moyal coefficients: IFT','Position', tmp_bar);
%     z               = 3;
%     dr_ind          = 1;
	dr_ind          = askInput({sprintf('Enter the separation of scales/step increment (in samples) referred to the sequence from large to small scales in the cascade trajectory:')},{num2str(1)});
    iter            = askInput({sprintf('Enter the maximum number of iteration used for the optimization algorithms:')}, {'5'});
    tol_D1          = askInput({sprintf('Optimization: Tolerance of the range of Kramers-Moyal coefficients in %%')}, {'10'})/100;
    tol_D2          = tol_D1;

    [co_IFT_opti,history]=OPTI_IFT_dij(tol_D1,tol_D2,evaluated,co_KM_opti,fitresult_D1_conf,fitresult_D2_conf,int_L,taylor_L,m_data,Fs,...
        data_filter,3,iter,markov,trajec,dr_ind,norm_ur,norm_r,save_path,save_name);
  % [history]       = DS_iter(tol_D1,tol_D2,history,evaluated,Fs,int_L,taylor_L,scale_steps,co_KM_opti,m_data,data_filter,z,markov,trajec,dr_ind,askYesno('co-calculation?', 'Yes'),askYesno('DS-calculation?', 'Yes'));
    % find(history.fval==min(history.fval),1,'first')
    KM_plot(co_IFT_opti,evaluated,Fs,taylor_L,int_L,multi_point,condition,norm_ur,norm_r,save_path,save_name)
    CO_plot(0,co_IFT_opti,co_KM_opti,evaluated,taylor_L,int_L,norm_ur,norm_r,siginf);
    if ischar(save_path)
        close all
        tmp_bar=bar.Position;
        close(bar)
        save(fullfile(save_path,append(save_name,'_','evaluated.mat')),'-v7.3') 
        bar = waitbar(21/ikx,'Calculation of Entropy: IFT optimized Kramers-Moyal coefficients','Position', tmp_bar);
    end

       
    answer = questdlg('The total entropy variation should be calculated for overlapping or independent cascade trajectories (using the optimized (IFT) Kramers-Moyal coefficients)?','Cascade Trajectories', ...
     'overlapping cascade trajectories', ...
     'independent cascade trajectories', ...
     'independent cascade trajectories');
    switch answer
        case 'overlapping cascade trajectories'
            z=1;
        case 'independent cascade trajectories'
            z=3;
    end  
%     data_length     = askInput({sprintf('How much %% of data do you want do use for the analysis?')}, {'100'})/100;
    [Sm,Ds,DS,r,ind_trajec,rind,dr,r_s,~,~,A,~,~,~,~,~,~] = checkFT(evaluated,int_L,taylor_L,Fs,data_filter,m_data,z,co_IFT_opti,data_length,markov,trajec,dr_ind,norm_ur,norm_r);
    plot_entropy(Sm,Ds,DS,increment_bin,save_path,save_name)
    if ischar(save_path)
        tmp_bar=bar.Position;
        close(bar)
        save(fullfile(save_path,append(save_name,'_','evaluated.mat')),'-v7.3') 
    end
end
end