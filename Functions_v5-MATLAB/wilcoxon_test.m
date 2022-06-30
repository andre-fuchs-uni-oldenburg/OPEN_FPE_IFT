%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function determines the Einstein-Markov length $\Delta_{EM}$\cite{renner2001}. Above this 
% length scale, the Markov properties hold and below this length scale, the Markov properties cease 
% to hold. The Wilcoxon test is a parameter-free procedure to compare two empirically determined 
% probability distributions (two data sets of velocity increments) (see \cite{Lueck2006markov} for 
% details). It is a quantitative test that determines the $\Delta_{EM}$. A sufficient resolution 
% in measurement below Taylor's length scale is expected to perform this test. Also, a vertical 
% dashed line at the Taylor length scale $\lambda$ will be added to the plot.
% COMMENT: for mulitpoint include Increment_point(tau1,tau2,d,condition,tol) 
%
% Arguments IN
% data = it is the filtered data 
% Fs =Acquisition/Sampling Frequency in Hz
% m_data = mean of the filtered data
% int_L = Integral length scale in meters
% taylor_L = Taylor length scale in meters
% diss_scale = Kolmogorv lenght scale in meters
% increment_bin = The number of bins in which you would like to divide your data
% f_avg = Frequency with smoothing;
% E_avg = Energy spectral density(ESD) with smoothing
%
% Arguments OUT
% markov = Einstein-Markov length scale in samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [markov,dt,U]=wilcoxon_test(data,Fs,m_data,int_L,taylor_L,diss_scale,increment_bin,save_path,save_name,f_avg, E_avg)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');
                       
NumBin          = increment_bin;
Numincr         = increment_bin*10;
% Numincr         = 60;

% m_data*length(data)/Fs
taylor_L_tmp    = taylor_L;
int_L_tmp       = int_L;

taylor_L        = 1*m_data/Fs;
% int_L = ceil(Fs*10*int_L/m_data)*m_data/Fs;
% floor(int_Ls/3)

% U=1.1;
L_ind=1;
% while min(U) > 1
    int_L=L_ind*int_L;
    
    int_Ls  = round(Fs*int_L/m_data);       % Get the integral length scale in unit of samples
    Nint    = floor(length(data)/int_Ls);   % Get the number of statistically independant intervals
    data    = data(1:Nint*int_Ls);          % Set data to a round number of independant intervals by removing the remainders at the end.

    % if taylor_L*10<int_L
    % % Incr=linspace(m_data/Fs,taylor_L*10,Numincr);               % Set the increments in samples "Delta r". We always use the whole range available for an independant interval 
    % % Incr=linspace(m_data/Fs,(int_L-taylor_L)/3,Numincr);        
    % Incr=logspace(log10(taylor_L/50),log10(taylor_L*10),Numincr); 
    % else
    % Incr=logspace(log10(taylor_L/50),log10(int_L-taylor_L),Numincr);
    % end    
%     Incr    = unique(round(logspace(log10(taylor_L*10^6),log10((int_L-taylor_L)*10^6),Numincr))./10^6,'stable');  % Give the number of increments  
    Incr    = unique(round(logspace(log10(taylor_L*10^6),log10((int_L-taylor_L)/2.5*10^6),Numincr))./10^6,'stable');  % Give the number of increments
    
    dt      = unique(floor(Incr*Fs/m_data),'stable');              % Calculate the scale increment in unit of samples
    %%
%     NumBin          = 201
%     Numincr         = 60
%     dt              = unique(round(linspace(1,150,Numincr)),'stable');  % Give the number of increments
    %%
    
    dt      = dt(dt>0);
    

    Numincr=length(dt);

    dr=nan(Numincr,1);                                      % Preallocate a dr matrix to store the values during the loop and use it to plot at the end
    U=nan(Numincr,1);                                       % Preallocate a vector that will contain the mean values of T=abs(Q-mean(Q))/std(Q). 

    for ii = 1:Numincr
        ii/Numincr
                                                            % We need to erase the previous r scale for each iteration.
        t=round(Fs.*taylor_L*ones(1,4)./m_data);            % Calculate the scale in unit of samples % dt(ii)=floor(Incr(ii)*Fs/m_data);                                                    
        dr(ii)=dt(ii)*m_data/Fs;                            % We calculate dr from dt because dt MUST be a constant round number.
        t=t+(0:3)*dt(ii);                                   % Then we build a vector that contains the incremented scale
        r=t*m_data/Fs;                                      % And we calculate r and its incremented values from t because t=constant round number.

        AB=nan(Nint,4);                                     % Preallocate a matrix to store the values at each t, in each interval, to calculate the V increments.

        if max(t)<=int_Ls
            for i=1:Nint                                        
                interv=data(1+(i-1)*int_Ls:i*int_Ls);           % 1st interval from 1 to int_Ls, 2nd from int_Ls+1 to 2*int_Ls......
                AB(i,:)=[interv(t(1)), interv(t(2)), interv(t(3)), interv(t(4))]; % each line corresponds to 1 interval. From each interval, we extract 3 increments (1,2,3*Delta t), and the 1st one is the reference (t).
            end
        Vr=nan(Nint,length(t)-1);                           % Build a matrix with the right size to contain 3 vectors of Nint length.
        Vr(:,1)=AB(:,2)-AB(:,1);                            % 1st column = 1st increment = 1st delta t - reference t
        Vr(:,2)=AB(:,3)-AB(:,1);                            % 1st column = 2nd increment = 2nd delta t - reference t
        Vr(:,3)=AB(:,4)-AB(:,1);                            % 1st column = 3rd increment = 3rd delta t - reference t
        clear interv i AB

        [~,Binctr]=hist(Vr(:,2),NumBin);                    % Create a binning, that only depends on the 2nd increment (1st condition)
        Binwidth=Binctr(2)-Binctr(1);                       % We also need the bin size to create their boundaries.

        C3=abs(Vr(:,3))<=Binwidth/2;                        % The second condition is only observed around 0 (highest data concentration). We take the same binsize as 1st condition
        V13=Vr(C3==1,1);                                    % Then we find the values of 1st and 2nd increments considering the 2nd condition.
        V23=Vr(C3==1,2);                                    % p(V1|V3) line 44(above) and p(V2|V3) line 45 (this line)

        T=nan(NumBin,1);                                    % preallocate an empty vector for T=abs(Q-mean(Q))/std(Q)

        for i = 1: NumBin                                   % We will now calculate Q and T for each bin
            BP=Vr(:,2)<(Binctr(i)+Binwidth/2) & Vr(:,2)>=(Binctr(i)-Binwidth/2); % Vector that contains the positions of the values in the bin i, for the single conditionned V1
            m=sum(BP);                                      % number of values in the bin i, with a single condition
            V11=Vr(BP==1,1);                                % Single conditionned V1 (1st increment), for the bin i
            BP2=V23<(Binctr(i)+Binwidth/2) & V23>=(Binctr(i)-Binwidth/2); % Vector that contains the positions of the values in the bin i, for the double conditionned V1
            n=sum(BP2);                                     % number of values in the bin i, with a double condition
            V12=V13(BP2==1);                                % Double conditionned V1 (1st increment), for the bin i

            if n>=25 && m>=25                               % We need this condition, to be sure that we have enough data to calculate mean and have a Gaussian behaviour.
    %       if n>=1 && m>=1 
                sortedx = sort(V11);                        % This line is unnecessary, but it follows what is explained in M. Waechter's article in 2004 (rough surfaces, eq (18))
                sortedy = sort(V12);                        % This is however very important. We need to count how many values of V11 are below each single values of V12
                Q = 0;
                for idx = 1:length(sortedy)                 % So we loop on the length of V12 to take all the values 1 by 1
                    Q = Q + sum(sortedy(idx) > sortedx);    % We count how many values of V11 are under V12(idx), and sum them to get "the total number of inversions"
                end
                Qm=n*m/2;                                   % Theoretical mean(Q), according to M. Waechter's article in 2004 (rough surfaces, eq (19))
                SigQ=sqrt(n*m*(n+m+1)/12);                  % Theoretical std(Q), according to M. Waechter's article in 2004 (rough surfaces, eq (19))
                T(i)=abs(Q-Qm)/SigQ; 
            end
        end
            U(ii)=nanmean(T)/sqrt(2/pi);                    % from T we calculate the mean. This mean is supposed to be close to sqrt(2/pi) if the test is fulfilled. 
                                                            % Thus, we normalize mean(T) by this value, and we should have something around 1.
        end
    end
    
%     L_ind=L_ind+0.5;
% end
% max(dt)
% largest dt where U ist not nan aprox int_Ls/3
int_Ls/3
int_Ls

% Numincr/sum(~isnan(U))
% dt=dt(~isnan(U));
% U=U(~isnan(U));

h(1) = figure;
subplot(1,2,1)
% Once we have all the values of mean(T), for all the increment values, we can plot them.
loglog(dt,U,'MarkerEdgeColor','k','MarkerSize',6,'Marker','o', 'LineStyle','none','LineWidth',1.5)
hold on
vline(floor(taylor_L_tmp/(m_data/Fs)),'k') 
% vline(ceil(0.9*(Fs*taylor_L/m_data)),'k')
ylabel('$\frac{\left<t(r,\Delta r)\right>}{\sqrt{2/\pi}}$','Interpreter','latex')% mean(T) is normalized by sqrt(2/pi)
xlabel('$\Delta r/sample$','Interpreter','latex')            % ('$\frac{dr}{\lambda}$','Interpreter','latex')% dr is normalized by the taylor length
% grid on
% hline=refline([0 1]);                                   % The reference line is taken at 1 as mean(T) has been normalized
% hline.Color = 'r';
plot(linspace(0,1000,1000),ones(1000,1),'r','LineWidth',2)
% title(['$\lambda=' num2str(num2str(floor(taylor_L/(m_data/Fs)))) 'sample$'],'interpreter','latex')
set(gca,'FontSize',16)
axis square
set(gcf, 'Color', 'w')
xlim([min(dt) max(dt)])
ylim([0.5 max(U)])
fig_setup
set(gca,'XTick',[10^0  10^1  10^2 10^3]);
set(gca,'yTick',[10^0  10^1  10^2 10^3]);
legend off
txt = {'(a)'};
a=text(4,0.5,txt);
a.FontSize = 22;
a.Units='normalized';

% h(2) = figure;
subplot(1,2,2)
% Once we have all the values of mean(T), for all the increment values, we can plot them.
plot(dt,U,'MarkerEdgeColor','k','MarkerSize',6,'Marker','o', 'LineStyle','none','LineWidth',1.5)                                        
hold on
vline(floor(taylor_L_tmp/(m_data/Fs)),'k') 
ylabel('$\frac{\left<t(r,\Delta r)\right>}{\sqrt{2/\pi}}$','Interpreter','latex')% mean(T) is normalized by sqrt(2/pi)
xlabel('$\Delta r/sample$','Interpreter','latex')            % ('$\frac{dr}{\lambda}$','Interpreter','latex')% dr is normalized by the taylor length
% grid on
% title(['$\lambda=' num2str(num2str(floor(taylor_L/(m_data/Fs)))) 'sample$'],'interpreter','latex')
set(gca,'FontSize',16)
set(gca,'xScale','log')
% hline=refline([0 1]);                                   % The reference line is taken at 1 as mean(T) has been normalized
% hline.Color = 'r';
plot(linspace(0,1000,1000),ones(1000,1),'r','LineWidth',2)
axis square
set(gcf, 'Color', 'w')
xlim([min(dt) max(dt)])
ylim([0 max(U)])
datacursormode on
fig_setup
set(gca,'XTick',[10^0  10^1  10^2 10^3]);
legend off
txt = {'(b)'};
b=text(4,0.5,txt);
b.FontSize = 22;
b.Units='normalized';
pos_txt=[-0.19   0.9 0];
a.Position=pos_txt;
b.Position=pos_txt;

set(gcf,'Position',[238 898 921 420]) 


waitfor(msgbox('Use interactive data cursor to select the Markov length in samples (see readme).'));
tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf);
delete(tmp_ui);

markov = askInput({sprintf('Markov length in samples:')}, {num2str(0.9*(Fs*taylor_L_tmp/m_data)  ,'%1.0f')});

h(2) = figure;
set(gcf, 'Color', 'w')
yyaxis left
loglog(f_avg,E_avg,'-','LineWidth',2)%% Plot of Frequency Vs. Power Spectral Density(m^2/sec) with smoothing
xlabel('f / Hz','interpreter','latex');
ylabel('E(f) / $\frac{m^2}{s}$', 'interpreter','latex')
set(gca,'FontSize',10)
set(gcf, 'Color', 'w')
xlim([min(f_avg) max(f_avg)]*2)
% xlim([min(f_avg) 10^2])
grid on
axis square
datacursormode on
fig_setup
grid off
legend off
set(gca,'XTick',[10^-2 10^0  10^2  10^4]);
vline((1/(int_L_tmp/m_data))/(2*pi),'k')
vline((1/(taylor_L_tmp/m_data))/(2*pi),'k')
vline((1/(diss_scale/m_data))/(2*pi),'k')
vline((1/((markov*m_data/Fs)/m_data))/(2*pi),'r')

yyaxis right
loglog(f_avg,E_avg.*f_avg.^(5/3),'-','LineWidth',2)%% Plot of Frequency Vs. Power Spectral Density(m^2/sec) with smoothing
ylabel('$E(f)f^{5/3}$ / $\frac{m^2}{s}$', 'interpreter','latex')
xlim([min(f_avg) max(f_avg)]*2)
ylim([min(E_avg.*f_avg.^(5/3)) max(E_avg.*f_avg.^(5/3))]*2)


if ischar(save_path)     
    savefig(h,fullfile(save_path,append(save_name,'_','wilcoxon_test.fig')),'compact')
    for a = 1:length(h)
        exportgraphics(h(a),fullfile(save_path,append(save_name,'_',sprintf('wilcoxon_test_%d.png', a))))
    end
end
tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf);
delete(tmp_ui);
close all
end