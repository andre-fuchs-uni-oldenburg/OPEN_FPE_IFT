%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The calculation leading towards the integral fluctuation theorem will be done. In the spirit of 
% non-equilibrium stochastic thermodynamics \cite{seifert2012stochastic} it is possible to associate 
% with every individual cascade trajectory $\left[u(\cdot) \right]$ a total entropy variation 
% $\Delta S_{tot}$. In this investigation it is assumed that a single cascade trajectory represents 
% one realization of the turbulent cascade process and a large number of these trajectories reflect 
% the statistics caused by the process.
% The set of measured cascade trajectories results in a set of total entropy variation values 
% $\Delta S_{tot}$. This function calculates the system entropy, medium entropy and the total 
% entropy variation for all the independent cascade trajectories.
% 
% Arguments IN
% evaluated = a modified/updated struct 'evaluated' array with in function 'KM_STP_optimization'
% int_L = Integral length scale in meters
% taylor_L = Taylor length scale in meters
% Fs = Acquisition/Sampling Frequency in Hz
% data = filtered data 
% m_data = mean of the data
% z = This is the parameter which decides how to compute entropy using different methods using either 
% overlapping or independent trajectories
% z = 1 ==> Overlapping trajectories 
% z = 3 ==> Independent trajectories
% co = Coefficients $d_{ij}(r)$ of the optimized Kramers-Moyal coefficients, where first entry is 
%      zeroth power and power laws for r-dependency with exponents co.ae and co.be, for each power 
%      in D1 and D2 respectively
% data_length = A value between 0 & 1 ==> how much of data you want to consider for this calculation
% data_length = 1 for all the data i.e. data(1:end)
% data_length = 0.5 for the data i.e data(1:end/2)   <== For 50% of the data 
% markov = markov length in number of samples
% trajec = 1 ==> The start/end of the cascade trajectory will be adjusted.
% dr_ind = Separation of scales/step increment (in samples) referred to the sequence from large to small scales in the cascade trajectory
% norm_ur = normalization of the data using $\sigma_\infty$? data 1=Yes, 0=No
% norm_r = normalization of the scale using $\lambda$? data 1=Yes, 0=No
% 
% Arguments OUT
% Sm = Entropy of the medium
% Ds = Entropy of the system or the Shanon entropy
% DS = Total entropy production=Sm+Ds
% r =  Normalized scale vector from the start to end of the cascade trajectory
% ind_trajec = index-vector of the trajectories (important vector for the recalculation of the trajectories from data)
% rind = index-vector from the start to end of the cascade trajectory
% dr = separation of scales/step increment (in normalized scale) referred to the sequence from 
%      large to small scales in the cascade trajectory
% r_s = mid-point scale vector (Stratonovich convention)
% u = used trajectoires (only z=3 ==> Independent trajectories) 
% A = action functional, pathintegral of Lagrangian
% Lag,p,H,tmpvar,H1,H2,ur_s_tmp (are not yet included/used in the current version)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Sm,Ds,DS,r,ind_trajec,rind,dr,r_s,u,Lag,A,p,H,tmpvar,H1,H2,ur_s_tmp]=checkFT(evaluated,int_L,taylor_L,Fs,data,m_data,z,co,data_length,markov,trajec,dr_ind,norm_ur,norm_r)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');

%% Fix length of trajectories
if trajec==1
        L       = 0;
        % Generally each Overlapping/Non-Overlapping trajectory starts at L and
        % ends at Lambda for the calculation of the entropy; but its a choice
        % So we ask you at which biggest scale you want to start the trajectory
        tmp     = str2double(inputdlg({'Adjustment: Start of trajectory 1=Yes, 0=No'}));
        if tmp == 1
            % Here, if you want to start your trajectories at scale smaller than L then
            % you should enter a value between 0 & 1 
            % For ex: if your L=1 meters and you enter 0.9 your trajectory will start at 0.9 meters
            L   = str2double(inputdlg({'Start of trajectory in multiples of Integral length scale:'}))*int_L/taylor_L;
        else
            L   = int_L/taylor_L; % dimless integral scale
        end
        % So we ask you at which smallest scale you want to end the trajectory
        tmp     = str2double(inputdlg({'Adjustment: End of trajectory 1=Yes, 0=No'}));
        if tmp == 1
            % Here, if you want to end your trajectories at scale bigger than lambda then
            % you should enter a value between 1 and (L/lambda) 
            % For ex: if your lambda=1 cm and if you enter 5 your
            % trajectory will end at 5 cms
            l   = str2double(inputdlg({'End of trajectory in multiples of Taylor length scale:'}))*taylor_L/taylor_L;
        else
            l   = evaluated(1).r/taylor_L;% dimless smallest scale by D1 and D2
        end
else % start=L, end=lambda
    L   = int_L/taylor_L;
    l   = 1;
end

x       = m_data.*(1:ceil(1.1*L*Fs*taylor_L/m_data)).'./Fs;
x       = x/taylor_L; % dimless position
Lind    = find(abs(x(:,1)-L  )==min(abs(x(:,1)-L  )), 1, 'last' ); % index for x showing the LAST  entry where x is within range lambda..L
lind    = find(abs(x(:,1)-l  )==min(abs(x(:,1)-l  )), 1, 'first'); % index for x showing the FIRST entry where x is within range lambda..L

% Index for calculating the increment using data
Lind    = Lind+1;   
lind    = lind+1;

if lind<2
    lind = 2;
end

rindful = 1:Lind;                       % index-vector for all scales
rind    = fliplr(rindful(lind:Lind));   % index-vector only for scales from start to end of the cascade trajectory
% rind = round(linspace(Lind,lind,floor(L/l)+1)); % index-vector for scales in inertial range in steps of lambda
stp     = Lind;                         % index used to cut v(x) in non-overlapping pieces of length equal to integral scale
% R(k).r   = x(rind); % scales from index-vector

if norm_r==1
    r       = x(rind); % scales from index-vector 
else
    r       = x(rind)*taylor_L; % scales from index-vector 
end

clear x
clear Sm Ds DS tmp_Sm tmp_Ds tmp_DS not_nan_index u x



%% calc Entropy - overlapping trajectories
if z==1
    %% PDF of start/end of the cascade trajectory ==> Entropy of the system or the Shanon entropy
    rind_tmp    = rind(1:dr_ind:end);
    end_L       = size(data,1);
    tmp_size    = size(data,1)-stp+1;
    u_tmp       = nan( tmp_size, 2 );

    for i=0:end_L-stp
        %% u
        % u_tmp(i+1,:) = (data([rind_tmp(1)+i,rind_tmp(end)+i])- data(1+i)).';
        %% u_s
        u_tmp(i+1,:) = (data([rind_tmp(2)+i,rind_tmp(end-1)+i])- data(1+i)).';
    end   

    nbins=301;
    [fuL, uL]       = hist(u_tmp(:,1),nbins); 
    puL             = fuL / ( nansum(fuL) * mean(diff(uL)) ); % relative frequency to approximate p(u,r=L)

    [ful, ul]       = hist(u_tmp(:,2),uL); 
    pul             = ful / ( nansum(ful) * mean(diff(ul)) ); % relative frequency to approximate p(u,r=lambda)

    figure
    plot(uL,puL,'LineWidth',2)
    hold on
    plot(ul,pul,'LineWidth',2)
    norm1 = pdf(fitdist(u_tmp(:,1),'Normal'),uL);
    plot(uL,norm1,'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176])
    norm1 = pdf(fitdist(u_tmp(:,2),'Normal'),uL);
    plot(ul,norm1,'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176])
    set(gca,'yScale','log') 
    set(gcf, 'Color', 'w')
    set(gca, 'FontSize',18)
    xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
    ylabel('$p(u_r / \sigma_\infty)$','interpreter','latex')
    axis square
    legend('u_L','u_\lambda');
    xlim([min(uL) max(uL)])
    ylim([0.0000001*max(puL) 2*max(pul)])
    tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
    uiwait(gcf);
    delete(tmp_ui);
    close
    clear u_tmp

    
    %% Preallocating (considering the available RAM)
    data        = data(1:ceil(length(data)*data_length));
    giga        = askInput({sprintf('Enter available RAM in GB:')},{num2str(8)}) * 1.07*10^9; %1 GB in bytes   
    k_split     = ceil((12*4*length(data)*length(r))/giga)
    data        = reshape(data((1:floor(numel(data)/k_split)*k_split)),[],k_split);

    A           = nan( size(data,1)-stp+1,k_split );
    Sm          = nan( size(data,1)-stp+1,k_split );
    Ds          = nan( size(data,1)-stp+1,k_split );
    DS          = nan( size(data,1)-stp+1,k_split );
    ind_trajec  = nan( size(data,1)-stp+1,k_split );
    %  Q_hk     = nan( size(data,1)-stp+1,k_split );
    %  Q_med    = nan( size(data,1)-stp+1,k_split );
    %  Q_ex     = nan( size(data,1)-stp+1,k_split );
    %  Q_med_xi = nan; 

    clear tmp
    tmp(1:k_split) = struct('Sm_tmp',nan(size(data,1)-stp+1,1),'Ds_tmp',nan(size(data,1)-stp+1,1),'DS_tmp',nan(size(data,1)-stp+1,1), 'ind_trajec_tmp',nan(size(data,1)-stp+1,1), 'A_tmp',nan(size(data,1)-stp+1,1));
    % tmp(1:k_split) = struct('Sm_tmp',nan(size(data,1)-stp+1,1),'Ds_tmp',nan(size(data,1)-stp+1,1),'DS_tmp',nan(size(data,1)-stp+1,1), 'ind_trajec_tmp',nan(size(data,1)-stp+1,1),...
    %'Q_hk_tmp',nan(size(data,1)-stp+1,1),'Q_med_tmp',nan(size(data,1)-stp+1,1),'Q_ex',nan(size(data,1)-stp+1,1)); 
    f = uifigure;
    d = uiprogressdlg(f,'Title','Please Wait','ShowPercentage','on');
    
    
    %% calc Entropy
    for k=1:k_split
        [u,ind_trajec_tmp]                                      = create_us_overlap(data(:,k),r,rind,stp);     
        [Sm_tmp,Ds_tmp,DS_tmp,dr,r_s,~,~,~,A_tmp,~,~,~,~,~,~]   = calcDS(z,u,r.', co,markov,taylor_L,Fs,m_data,dr_ind, uL, puL, ul, pul);
        tmp(k).A_tmp(1:length(A_tmp))                           = A_tmp;
        tmp(k).Sm_tmp(1:length(Sm_tmp))                         = Sm_tmp;
        tmp(k).Ds_tmp(1:length(Ds_tmp))                         = Ds_tmp;
        tmp(k).DS_tmp(1:length(DS_tmp))                         = DS_tmp;
        % [Q_hk_tmp, Q_med_tmp, Q_ex_tmp,~]                     = heat(u,r.', co);
        % tmp(k).Q_hk_tmp(1:length(Q_hk_tmp))                   = Q_hk_tmp;
        % tmp(k).Q_med_tmp(1:length(Q_med_tmp))                 = Q_med_tmp;
        % tmp(k).Q_ex_tmp(1:length(Q_ex_tmp))                   = Q_ex_tmp;

        if k==1
            tmp(k).ind_trajec_tmp(1:length(ind_trajec_tmp))     = ind_trajec_tmp;
        else
            % tmp(k).ind_trajec_tmp(1:length(ind_trajec_tmp))   = (size(data,1)-stp)*(k-1)+ind_trajec_tmp+1;
            % tmp(k).ind_trajec_tmp(1:length(ind_trajec_tmp))   = tmp(k-1).ind_trajec_tmp(length(ind_trajec_tmp))+ind_trajec_tmp+1;
            % tmp(k).ind_trajec_tmp(1:length(ind_trajec_tmp))   = tmp(k-1).ind_trajec_tmp(length(ind_trajec_tmp))+stp+ind_trajec_tmp+1;
            tmp(k).ind_trajec_tmp(1:length(ind_trajec_tmp))     = tmp(k-1).ind_trajec_tmp(length(ind_trajec_tmp))+stp+ind_trajec_tmp;
        end   
        d.Value=k/k_split;
		k/k_split
    end
    
    close(d);
    close(f);

    for k=1:k_split
        A(1:length(tmp(k).A_tmp),k)                             = tmp(k).A_tmp;
        Sm(1:length(tmp(k).Sm_tmp),k)                           = tmp(k).Sm_tmp;
        Ds(1:length(tmp(k).Ds_tmp),k)                           = tmp(k).Ds_tmp;
        DS(1:length(tmp(k).DS_tmp),k)                           = tmp(k).DS_tmp;
        ind_trajec(1:length(tmp(k).ind_trajec_tmp),k)           = tmp(k).ind_trajec_tmp;
    %     Q_hk(1:length(tmp(k).Q_hk_tmp),k)                     = tmp(k).Q_hk_tmp;
    %     Q_med(1:length(tmp(k).Q_med_tmp),k)                   = tmp(k).Q_med_tmp;
    %     Q_ex(1:length(tmp(k).Q_ex_tmp),k)                     = tmp(k).Q_ex_tmp; 
    end
    clear tmp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calc Entropy - non-overlapping trajectories
elseif z==3
    %% Preallocating (considering the available RAM)
    tmp_size    = ceil(numel(data)/stp);
    data        = data(1:ceil(length(data)*data_length));
%     giga        = askInput({sprintf('Enter available RAM in GB:')},{num2str(8)}) * 1.07*10^9; %1 GB in bytes 
    giga        = 8* 1.07*10^9;
    k_split     = ceil((20*4*tmp_size*length(r(1:dr_ind:end)))/giga);
    data        = reshape(data((1:floor(numel(data)/k_split)*k_split)),[],k_split);
   
    Sm          = nan( tmp_size, k_split );
    Ds          = nan( tmp_size, k_split );
    DS          = nan( tmp_size, k_split );
    ind_trajec  = nan( tmp_size, k_split );
    A           = nan( tmp_size, k_split );    
%     Lag         = nan( tmp_size, k_split );
%     p           = nan( tmp_size, k_split );
%     H           = nan( tmp_size, k_split );
%     tmpvar      = nan( tmp_size, k_split );
%     H1          = nan( tmp_size, k_split );
%     H2          = nan( tmp_size, k_split );
%     u           = nan( tmp_size, k_split );
%     Q_hk      = nan( tmp_size,1 );
%     Q_med     = nan( tmp_size,1 );
%     Q_ex      = nan( tmp_size,1 );
%     not_nan_index = nan( tmp_size, k_split );
 
    
    %% calc Entropy
    f = uifigure;
    d = uiprogressdlg(f,'Title','Please Wait','ShowPercentage','on');
    for k=1:k_split
        [u_tmp,ind_trajec_tmp]                  = create_us(data(:,k),r,rind,stp);
        [Sm_tmp,Ds_tmp,DS_tmp,dr,r_s,~,~,L_tmp,A_tmp,p_tmp,H_tmp,tmpvar_tmp,H1_tmp,H2_tmp,ur_s_tmp]           = calcDS(z,u_tmp,r.',co,markov,taylor_L,Fs,m_data,dr_ind);
    %   [Q_hk_tmp, Q_med_tmp, Q_ex_tmp,Q_med_xi] = heat(u_tmp,r.', co, r(end)); 

        A(1:sum(~isnan(DS_tmp),1),k)            = A_tmp(~isnan(DS_tmp));
        ind_trajec(1:sum(~isnan(DS_tmp),1),k)   = ind_trajec_tmp(~isnan(DS_tmp));
        Sm(1:sum(~isnan(DS_tmp),1),k)           = Sm_tmp(~isnan(DS_tmp));
        Ds(1:sum(~isnan(DS_tmp),1),k)           = Ds_tmp(~isnan(DS_tmp));
        DS(1:sum(~isnan(DS_tmp),1),k)           = DS_tmp(~isnan(DS_tmp));
        
        if k_split==1
%             Lag(1:sum(~isnan(DS_tmp),1),:)         = L_tmp(~isnan(DS_tmp),:);
%             p(1:sum(~isnan(DS_tmp),1),:)           = p_tmp(~isnan(DS_tmp),:); 
%             H(1:sum(~isnan(DS_tmp),1),:)           = H_tmp(~isnan(DS_tmp),:); 
%             tmpvar(1:sum(~isnan(DS_tmp),1),:)      = tmpvar_tmp(~isnan(DS_tmp),:);
            ur_s_tmp(1:sum(~isnan(DS_tmp),1),:)    = ur_s_tmp(~isnan(DS_tmp),:);
%             H1(1:sum(~isnan(DS_tmp),1),:)          = H1_tmp(~isnan(DS_tmp),:);   
%             H2(1:sum(~isnan(DS_tmp),1),:)          = H2_tmp(~isnan(DS_tmp),:);
            u(1:sum(~isnan(DS_tmp),1),:)           = u_tmp(~isnan(DS_tmp),:);
        end
%             Q_hk(1:sum(~isnan(DS_tmp),1),k)           = Q_hk_tmp(~isnan(DS_tmp));
%             Q_med(1:sum(~isnan(DS_tmp),1),k)          = Q_med_tmp(~isnan(DS_tmp));
%             Q_ex(1:sum(~isnan(DS_tmp),1),k)           = Q_ex_tmp(~isnan(DS_tmp));
%             not_nan                                   = (0:tmp_size).';
%             not_nan_index(1:sum(~isnan(DS_tmp),1),k)  = not_nan(~isnan(DS_tmp));
        clear ind_trajec_tmp Sm_tmp Ds_tmp DS_tmp u_tmp

     d.Value=k/k_split;
    end
    close(d);
    close(f);
end

A       = reshape(A,[],1);
Sm      = reshape(Sm,[],1);
Ds      = reshape(Ds,[],1);
DS      = reshape(DS,[],1);
% Q_hk  = reshape(Q_hk,[],1);
% Q_med = reshape(Q_med,[],1);
% Q_ex  = reshape(Q_ex,[],1);

% cut-out NaN's due to vanishing or negative D2
DS(isinf(DS))   = nan;
A               = A(~isnan(DS));

% if z==1
    Lag             = nan;
    p               = nan;
    H               = nan;
    tmpvar          = nan;
    H1              = nan;
    H2              = nan;
    if k_split>1
        ur_s_tmp        = nan;
        u               = nan;
    end
% end

% if z==3
    % L      = reshape(L,[],1);
    % p      = reshape(p,[],1);
    % H      = reshape(H,[],1);
    % tmpvar = reshape(tmpvar,[],1);
    % H1     = reshape(H1,[],1);
    % H2     = reshape(H2,[],1);
% end
% if z==3
%     Lag             = Lag(~isnan(DS),:);
%     p               = p(~isnan(DS),:);
%     H               = H(~isnan(DS),:);
%     tmpvar          = tmpvar(~isnan(DS),:);
%     ur_s_tmp        = ur_s_tmp(~isnan(DS),:);
%     H1              = H1(~isnan(DS),:);
%     H2              = H2(~isnan(DS),:);
 if k_split==1
     ur_s_tmp        = ur_s_tmp(~isnan(DS),:);
     u               = u(~isnan(DS),:);
 end
% end

Sm              = Sm(~isnan(DS));
Ds              = Ds(~isnan(DS));
% Q_hk  = Q_hk(~isnan(DS));
% Q_med = Q_med(~isnan(DS));
% Q_ex  = Q_ex(~isnan(DS));

ind_trajec      = reshape(ind_trajec,[],1);
ind_trajec      = ind_trajec(~isnan(DS));
% not_nan_index = ~isnan(DS);

DS              = DS(~isnan(DS));

if z==1
    Sm          = Sm(1:end-1);
    Ds          = Ds(1:end-1);
    DS          = DS(1:end-1); 
    ind_trajec  = ind_trajec(1:end-1);
end
r_s = r_s(1,:);
end