%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% In this function the integral length scale $L$, Taylor length scale $\lambda$ and Kolmogorov 
% length scale are estimated using different methods of calculation. Within this function, pop-up 
% dialog boxes will be generated to enter the values of the integral, Taylor length scale in $m$ 
% and energy dissipation rate in $m^2/s^3$ on which the further processing of data (solving the 
% FPE and extracting cascade trajectories) should be referred. The entered length scales will be 
% round towards the nearest integer sample. The proposed value in the pop-up dialog box is the 
% median length scale for all methods.
% COMMENT: All velocity increments in this script are taken as right hand increment i.e. 
% dV=V(t2)-V(t1) where t2>t1; t1,t2 in seconds
% 
% Arguments IN
% data = it is the filtered data 
% Fs = Acquisition/Sampling Frequency in Hz
% low_freq = Frequency at which you would like to use the low-pass filter in Hz
% kin_vis = Kinematic viscosity of fluid in m^2/sec
% C2 = constant, C2 associated with second order structure function(in between 2.0 to 2.2)
%
% Arguments OUT
% int_L = Integral length scale in meters
% taylor_L = Taylors length scale in meters
% int_L_calc = 1D array of integral length scale calculated by different methods/formula
% taylor_L_calc = 1D array of Taylor's length scale calculated by different methods/formula
% epsi = mean energy dissipation rate
% epsi_calc = 1D array of mean energy dissipation rate calculated by different methods/formula
% diss_scale = Kolmogorv lenght scale in meters
% Ce = is the non-dimentional / normalized energy dissipation rate
% Ce_calc = 1D array of Ce calculated by different methods/formula
% Re = Reynolds Number  based on integral length scale
% Re_lambda= Reynolds Number based on Taylor's length scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [int_L,taylor_L,int_L_calc, taylor_L_calc,epsi, epsi_calc,diss_scale,Ce, Ce_calc, Re, Re_lambda]=length_scales(data,Fs,low_freq,kin_vis,C2,increment_bin,m_data,save_path,save_name)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');

maxnumattempts      = 2;


% Pressure calculation: 
% http://www.peacesoftware.de/einigewerte/luft_e.html
% https://www.einheiten-umrechnen.de/einheiten-rechner.php?typ=druck
% The higher the pressure the lower the kinematic viscosity, 
% the higher the temperature the higher the kinematic viscosity.
% High Pressure: 1040hPa = 1,04 Bar, at 0 Celsius: 13.00*10^-6 at 20 Celsius: 14.73*10^-6
%                1000hPa = 1 Bar,    at 0 Celsius: 13.52*10^-6 at 20 Celsius: 15.32*10^-6
% Low pressure:  960hPa  = 0,96 Bar, at 0 Celsius: 17.24*10^-6 at 20 Celsius: 14.73*10^-6

% SHREK: Helium 
% http://hbcponline.com/faces/documents/06_45/06_45_0001.xhtml?search=true
% 10^-6 *10*10^-4*3.259/0.1448402
% kin_vis = 3.259*10^-6 / (0.1448402*10^3) % at 2 Kelvin
% kin_vis = 1.468*10^-6 / (0.1456217*10^3) % at 2.5 Kelvin
% kin_vis = 2.1*10^-8

if mod(length(data),2)==1
    data=data(1:end-1);% This is to make length of data array even
end  

% mean of the filtered data which is nothing but the mean of unfiltered data==>Filtering does not change the mean value
% m_data          = mean(data); 
% standard deviation of the filtered data which will always be lesser than the standard deviation of the unfiltered data
uprime          = std(data); 
% variance of the data
uprimesquared   = uprime^2;

% PSD of the fluctuations after filtering; mean value does not play a role in spectrum
daten           = data-mean(data);
L               = length(daten); %Length of the data

% spek              = abs(fft(daten,L)).^2/L;           % Power spectral density(PSD)
spek            = abs(fft(daten,L)).^2/(L*Fs);          % Energy spectral density(ESD) using a fft
spek            = 2.*spek(2:L/2+1);                     % FFt will yield half number of unique points

f = Fs/2*linspace(0,1,L/2+1).';                         % Nyquist frequency of the signal==Fs/2
f = f(2:end);                                           % Remove zero Hz component


% moving average with equally spaced frequency interval in log-space 
plot_length         = increment_bin*10;
if mod(plot_length,2) == 0
    plot_length = plot_length +1;
end

intervall       = unique(round(logspace(0,log10(L/2),plot_length)),'stable');
plot_length     = length(intervall);

% Initializing the arrays/preallocation
spek_smooth         = zeros(plot_length-2,1);
x_achse_f_spektrum  = zeros(plot_length-2,1);

% Averaging of spectrum and hence smoothing
for i=2:(plot_length-1)
    x_achse_f_spektrum(i-1,1) = mean(f(intervall(i-1):intervall(i+1)));
    spek_smooth(i-1,1)        = mean(spek(intervall(i-1):intervall(i+1)));
end
% ------------------------Windowing / averaging end ------------------------  

E_k                     = ((spek_smooth.*m_data)./(2*pi));
% 2-pi - k Spektrum; Wavenumbers(1/m)
x_achse_k_spektrum_2pi  = (2.*pi.*x_achse_f_spektrum)./m_data;



% Wavenumber Spectrum (m^3/sec^2)
k           = (2.*pi.*f(1:end))./m_data;
Ek          = ((spek(1:end).*m_data)./(2.*pi));
low_freq_k  = (2*pi*low_freq)/m_data;


h(1) = figure;
% Plot of Frequency Vs. Power Spectral Density(m^2/sec) without smoothing
loglog(f,spek) 
hold on
% Plot of Frequency Vs. Power Spectral Density(m^2/sec) with smoothing
loglog(x_achse_f_spektrum,spek_smooth,'-','LineWidth',2)
xlabel('f / Hz','interpreter','latex');
ylabel('E(f) / $\frac{m^2}{s}$', 'interpreter','latex')
xlabel('$f\ (Hz)$','interpreter','latex');
ylabel('$E\ (m^2/s)$', 'interpreter','latex')

% title('energy density spectrum','interpreter','latex')
set(gca,'FontSize',10)
set(gcf, 'Color', 'w')
% xlim([min(x_achse_f_spektrum) max(x_achse_f_spektrum)]*2)
xlim([min(x_achse_f_spektrum) 10^2])
grid on
axis square
% datacursormode on
fig_setup
legend off
grid off
set(gca,'XTick',[10^-2 10^0  10^2  10^4]);
waitfor(msgbox('Use interactive data cursor to select the range of frequency which will be used to linearly extrapolate the value of ESD at a frequency of 0 Hz (see readme).'));
% tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
% uiwait(gcf); 
% delete(tmp_ui);



%% First Method==> Calculating Ef0 with using Robust property
% Calculate the integral length scale (macro scale) Roach Eqn 15
% In this method we will ask you to choose two values of frequencies in Hz such as 
% f_start & f_end ===> [f_start<f_end ]<== 
% Must condition and so we will extrapolat a line designed bz f_start & f_end to find
% energy content at the limit f-->0
trying              = true;
attemptcount        = 1;
while trying
    try
        [I,yy]  = ginput(2);
        if I(1)>I(2)
            I   = flip(I);
            yy  = flip(yy);
        end     
        tmp_ui = uicontrol('Position',[10 10 100 40],'String','Processing','Callback','uiresume(gcbf)');

        f_start_plot     = askInput({sprintf('Lowest frequency of constant range in Hz:')}, {num2str(I(1),'%.8f')});
        f_end_plot       = askInput({sprintf('Highest frequency of constant range in Hz:')}, {num2str(I(2),'%.8f')});
        f_start     = find(f>f_start_plot, 1,'first');
        f_end       = find(f<=f_end_plot, 1,'last');

        % Ef0       = mean(spek(f_start:f_end));
        % Ef0       = median(spek(f_start:f_end));

        % IntegralLengthScale_1   = (Ef0*m_data)/(4*uprimesquared);           % Integral Scale in meters

        %Second Method==> Calculating Ef0 from fit without using Robust property
        [xData, yData] = prepareCurveData( f(f_start:f_end),(spek(f_start:f_end)));
        
        ft = fittype( 'poly1' );
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        % opts.Algorithm = 'Levenberg-Marquardt';
        % opts.Display = 'Off';
        opts.Lower = [0 0 ];
        % opts.StartPoint = [0.5 0.1];
        [fitresult, gof]        = fit( xData, yData, ft, opts );
        IntegralLengthScale_2   = (fitresult.p2*m_data)/(4*uprimesquared);   % Integral Scale in meters
        fitresult.p1=0;
        plot(f,feval(fitresult,f),'k','LineWidth',2,'LineStyle','-')
        vline(f_start_plot,'k')              % Construct a verticle line at f_start
        vline(f_end_plot,'k')                % Construct a verticle line at f_end
        legend('raw','averaged','Location','southwest')
        
        trying = false;
    catch
        waitfor(msgbox('Please try again: f_start should be smaller than f_end and the decimal separator should be "."(period).'));   
        if attemptcount>maxnumattempts
            waitfor(msgbox('Function failed.', 'Error','error'));   
            error('Function failed after %d times',maxnumattempts)
        end
        attemptcount=attemptcount+1;
    end  
end       
      
        
% %Third Method==> Calculating Ef0 from fit with using Robust property='LAR' specifies the least absolute residual method.
% [xData, yData] = prepareCurveData( f(f_start:f_end),spek(f_start:f_end));
% ft = fittype( 'a*x*0+b', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Algorithm = 'Levenberg-Marquardt';
% opts.Display = 'Off';
% opts.Robust = 'LAR';
% opts.StartPoint         = [0.0 0.01];
% [fitresult, gof]        = fit( xData, yData, ft, opts );
% IntegralLengthScale_3   = (fitresult.b*m_data)/(4*uprimesquared);
% 
% %%Fourth Method==> Calculating Ef0 from fit with using Robust property='Bisquare' specifies the bisquare weights method.
% [xData, yData] = prepareCurveData( f(f_start:f_end),spek(f_start:f_end));
% ft = fittype( 'a*x*0+b', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Algorithm = 'Levenberg-Marquardt';
% opts.Display = 'Off';
% opts.Robust = 'Bisquare';
% opts.StartPoint         = [0.0 0.01];
% [fitresult, gof]        = fit( xData, yData, ft, opts );
% IntegralLengthScale_4   = (fitresult.b*m_data)/(4*uprimesquared);



%% Dissipation from different methods
trying              = true;
attemptcount        = 1;
while trying
    try
        clear r tau step stepp 
        r                   = unique(round(logspace(log10((1./Fs).*m_data*10^6),log10(20*IntegralLengthScale_2*10^6),increment_bin))./10^6,'stable');

        tau(1:length(r))    = nan;
        step(1:length(r))   = nan;

        for z=1:length(r)
            tau(z)          = r(z)/m_data;
            step(z)         = round(tau(z)*Fs);  
        end
        stepp              = unique(step);
        stepp              = stepp(stepp>0);

        clear r R S_exp_3 S_exp_2
        r(1:length(stepp)) = nan;
        R(1:length(stepp)) = nan;
        S_exp_3(1:length(stepp))=nan; % Third order structure function
        S_exp_2(1:length(stepp))=nan; % Second order structure function

        parfor z=1:length(stepp)
        %     z/length(stepp)   
            tau         = stepp(z)              % In steps, just a number
            r(z)        = (tau./Fs).*m_data;    % in meters
            incr        = (data(tau+1:L) - data(1:L-tau)).';   % velocity Increment in steps of tau
            S_exp_2(z)  = mean(incr.^2);        % Second order structure function
            S_exp_3(z)  = mean(incr.^3);        % Third order structure function   
            % Autocorr André
            R(1,z)      = mean(daten(tau+1:L).*daten(1:L-tau))/uprimesquared;
        end

        % HIT = Homogeneous Isotropic Turbulence
        eps_2   = (1./(r)).*((S_exp_2./C2).^(3/2)); % Estimation of dissipation from S_exp_2 using HIT
        eps_3   = -(5/4).*(S_exp_3./r);             % Estimation of dissipation from S_exp_3 using HIT 

        I0      = find( eps_2==max(eps_2),1);       % Finding the maximum value of dissipation in inertial range
        I1      = find(-eps_3==max(-eps_3),1);      % Finding the maximum value of dissipation in inertial range

        % dissipation rate m^2.s^-3 from S2
        epsi_S2 = mean(eps_2(1,I0-2:I0+2));         % finding the mean amongst 5 points closest to peak value;
        % dissipation rate m^2.s^-3 from S3
        epsi_S3 = mean(-eps_3(1,I1-2:I1+2));        % finding the mean amongst 5 points closest to peak value;
        
        trying = false;
    catch
%         waitfor(msgbox(' f_start should be smaller than f_end should be. ')); 
        increment_bin = askInput({sprintf('Please try again by increasing the number of Bin''s used (must be larger than:):')}, {num2str(increment_bin,'%1.0f')});      
        if attemptcount>maxnumattempts
            waitfor(msgbox('Function failed.', 'Error','error'));   
            error('Function failed after %d times',maxnumattempts)
        end
        attemptcount=attemptcount+1;
    end  
end       
        
delete(tmp_ui);


% figure
% loglog(r,S_exp_2)


%% Integral length scale from Auto-corelation or structure function 
%% Autocorr André
low_freq_r  =  m_data/low_freq;  

% Full integration 
int_L_full  = trapz(r,R);
% trapz(r,abs(R))

% integrate up to the first zero crossing
tmp_zero=(find(R<0,1,'first'))-1;
if isempty(tmp_zero)
   int_L_0     = nan;
else
    % r(tmp_zero)
    int_L_0     = trapz(r(1:tmp_zero),R(1:tmp_zero));
end

% integrate up to the first 1/e crossing
tmp_exp=(find(R<1/exp(1),1,'first'));
if isempty(tmp_exp)
    int_L_1_exp=nan;
    theo_exp_Funktion = nan;
else
    int_L_1_exp = trapz(r(1:tmp_exp),R(1:tmp_exp));
    % local minimum 
    % index1 = find(diff(sign(diff(R))) ~= 0,2,'first')+1;
    % r(index1)
    % int_L_min=trapz(r(1:index1),abs(R(1:index1)));

    %exp Fit to R; Fitregion r=1/exp(1):2*low_freq_r
    tmp_up      = (find(r>4*low_freq_r,1,'first'));
    tmp_low     = (find(R<1/exp(1),1,'first'));
    
    if tmp_up<tmp_low
        [xData, yData] = prepareCurveData( r(tmp_up:tmp_low), R (tmp_up:tmp_low));
        ft = fittype( 'exp1' );
        % ft = fittype( 'exp(-x/b)', 'independent', 'x', 'dependent', 'y' )
        % exp(-(x+a)/b)
        opts                = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Algorithm      = 'Levenberg-Marquardt';
        opts.Display        = 'Off';
        opts.Robust         = 'Bisquare';
        [fitresult, gof]    = fit( xData, yData, ft, opts );
        theo_exp_Funktion   = ((fitresult.a/-fitresult.b));
        % theo_exp_Funktion=sum(R(1:tmp_low).*diff(r(1:tmp_low+1)));
        % theo_exp_Funktion=(fitresult.a/-fitresult.b * exp(-fitresult.b*r(tmp_up)));
        % theo_exp_Funktion=sum(R(1:tmp_up).*diff(r(1:tmp_up+1)))+(fitresult.a/-fitresult.b * exp(fitresult.b*r(tmp_up)));
        % theo_exp_Funktion=sum(R(2:tmp_up).*diff(r(1:tmp_up)))+(fitresult.a/-fitresult.b * exp(fitresult.b*r(tmp_up)));
        
        % Plot the Auto-correlation coefficient Vs. Scale, r and exp fit
        h(2) = figure;
        if ~isempty(tmp_exp)
            plot(r,fitresult.a*exp(fitresult.b*r),'-r','LineWidth',2)
            % plot( fitresult, r, R )
            hold on
        end
        plot(r,R,'MarkerEdgeColor','k','MarkerSize',6,'Marker','o', 'LineStyle','none','LineWidth',1.5)
        ylabel('$R_{\widetilde{u}\widetilde{u}}$', 'interpreter','latex')
        xlabel('r/m', 'interpreter','latex')
        xlabel('$r\ (m)$','interpreter','latex');
        % title(['integral length scale L = ' num2str(theo_exp_Funktion) 'm'])
        % title(['first zero crossing = ' num2str(min(r(find(R<0)-1))) 'm'])
        legend('exp fit', 'Location', 'NorthEast' );
        axis square
        set(gcf, 'Color', 'w')
        set(gca, 'FontSize',18)
        % xlim([0 150*low_freq_r])
        if ~isempty(tmp_exp)
            xlim([0 4*r(tmp_low)])
            vline(r(tmp_low),'k')
            vline(r(tmp_up),'k')
        % vline(r(find(R<1/exp(1),1,'first')),'r')
        end
        fig_setup
    else
        theo_exp_Funktion=nan;  
    end
end
tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf);
delete(tmp_ui);
tmp_ui=uicontrol('Position',[10 10 100 40],'String','Processing','Callback','uiresume(gcbf)');
pause(1)



%% Following calculation ensures that auto-correlation fluctuates around zero. For that you need to take the value of standard deviation of velocity calculated from the 2nd order structure function
Ruu     = 1-(S_exp_2./(mean(S_exp_2((r>10*IntegralLengthScale_2)))));   % Auto-correlation coefficient Using std from S_exp_2 ==>This is a dimentionless quantity
Ruu     = [1,Ruu];                                                      % Auto-correlation is 1 at scale=0
r_Ruu   = [0,r];                                                        % Auto-correlation is 1 at scale=0
LL      = (cumsum(Ruu(1,1:end-1).*diff(r_Ruu)))';                       % Cumulative integration of Ruu with respect to r
pks     = findpeaks(LL,3);                                              % Finding the peaks in LL
Npks    = size(pks,1);

if Npks == 0 
IntegralLengthScale_5 = trapz(r_Ruu,Ruu);       % If no peaks are found then simply take the area under the curve which will be in meters
else  
IntegralLengthScale_5 = pks(1);                 % in meters At this length we should reach an asymptotic value on the plot of r Vs. LL
end

lambda_S2   = (sqrt(15*kin_vis*uprimesquared./epsi_S2));    % Taylors length scale in meters from epsi_S2 and hence from S_exp_2
lambda_S3   = (sqrt(15*kin_vis*uprimesquared./epsi_S3));    % Taylors length scale in meters from epsi_S3 and hence from S_exp_3



%% Integral length scale from auto-corelation function 'xcorr' of MATLAB 
% tau                     = 1/Fs;             % time in seconds/Time per sample in secs
% [Ruu_xcorr,tau_xcorr]   = xcorr(data - mean(data),size(data,1)/2,'unbiased');
% Ruu_xcorr               = Ruu_xcorr(tau_xcorr>=0);
% tau_xcorr               = tau_xcorr(tau_xcorr>=0);
% Ruu_xcorr               = Ruu_xcorr./max(Ruu_xcorr);
% dist_xcorr              = tau_xcorr.*tau.*m_data;
% Ruu_xcorr               = Ruu_xcorr';
% LL_xcorr                = (cumsum(Ruu_xcorr(1,1:end-1).*diff(dist_xcorr(1,1:end))));
% pks_xcorr               = findpeaks(LL_xcorr);
% Npks_xcorr              = size(pks_xcorr,1);
% 
% if Npks_xcorr == 0 
% IntegralLengthScale_6   = trapz(dist_xcorr,Ruu_xcorr);  % If no peaks are found then simply take the area under the curve
% else  
% IntegralLengthScale_6   = pks_xcorr(1);                 % in meters at this length we should reach an asymptotic value on plot of r Vs. LL
% end
% 
% figure;
% plot(r,R,'-k','LineWidth',2)
% hold on
% plot(r_Ruu,Ruu,'LineWidth',2)
% plot(dist_xcorr,Ruu_xcorr,'LineWidth',2)
% xlim([0 1500*low_freq_r])


% figure;
% plot(r_Ruu(1:end-1),LL,'LineWidth',2)
% hold on
% plot(dist_xcorr(1:end-1),LL_xcorr,'LineWidth',2)



%% Zero-crossing method: Martin Obligado
delete(tmp_ui);
[int_L_zero_cros,lambda_zero_cros,lambda1_zero_cros,Eps_zero_cros,h]=ZeroCrossings(data,Fs,m_data,kin_vis,low_freq,x_achse_k_spektrum_2pi,E_k,h);

%% 
% int_L_calc  = [IntegralLengthScale_2,int_L_full,int_L_0,int_L_1_exp,theo_exp_Funktion,IntegralLengthScale_5,IntegralLengthScale_6,int_L_zero_cros];%% meter
int_L_calc  = [IntegralLengthScale_2,int_L_full,int_L_0,int_L_1_exp,theo_exp_Funktion,IntegralLengthScale_5,int_L_zero_cros];%% meter

% Choose the best estimate of integral and Taylor length scale on which you would like to rely on
h(3) = figure;
plot(int_L_calc,'MarkerEdgeColor','k','MarkerSize',6,'Marker','o', 'LineStyle','none','LineWidth',1.5)
xlabel('method', 'interpreter','latex')
ylabel('L/m', 'interpreter','latex')
ylabel('$L\ (m)$','interpreter','latex');
axis square
set(gcf, 'Color', 'w')
set(gca, 'FontSize',18)
fig_setup
legend off
set(gca,'XTick',[1 2 3 4 5 6 7 8]);
datacursormode on
tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf);
delete(tmp_ui);
int_L    = askInput({sprintf('Enter the integral length scale in meter. Median integral length scale in meter:')}, {num2str(nanmedian(int_L_calc),'%.5f')});






%% Taylor length by different methods
%% using Autocorrelation function
trying              = true;
attemptcount        = 1;
while trying
    try
        figure;
        plot(r_Ruu,Ruu,'MarkerEdgeColor','k','MarkerSize',6,'Marker','o', 'LineStyle','none','LineWidth',1.5)
        ylabel('$R_{\widetilde{u}\widetilde{u}}$', 'interpreter','latex')
        xlabel('r/m', 'interpreter','latex')
        xlabel('$r\ (m)$','interpreter','latex');
        axis square
        set(gcf, 'Color', 'w')
        set(gca, 'FontSize',18)
        ylim([0.7 1.05])
        xlim([0 r_Ruu(find(Ruu<0.7,1,'first'))])
        datacursormode on
        waitfor(msgbox('Use interactive data cursor to select the range of scales which will be used to estimate Taylor lenght scalee (see readme).'));
        [I,yy]  = ginput(1);
        i_tmp = askInput({sprintf(['Fit Region Auto-correlation Coefficient e.g. = ' ,num2str(0.98,'%1.2f')])}, {num2str(yy(1),'%.2f')});
        close

     
        [xData, yData] = prepareCurveData( r_Ruu(1:(find(R<i_tmp,1,'first'))), R (1:(find(R<i_tmp,1,'first'))));
        ft = fittype( '1-(x^2/a^2)', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Lower = 0;
        opts.StartPoint = 0.002;
        [fitresult, gof] = fit( xData, yData, ft, opts );
           
        h(4) = figure;
        plot(r_Ruu,feval(fitresult,r_Ruu),'r','LineWidth',2)
        hold on
        plot(r_Ruu,Ruu,'MarkerEdgeColor','k','MarkerSize',6,'Marker','o', 'LineStyle','none','LineWidth',1.5)
        ylabel('$R_{\widetilde{u}\widetilde{u}}$', 'interpreter','latex')
        xlabel('r/m', 'interpreter','latex')
        xlabel('$r\ (m)$','interpreter','latex');
        axis square
        set(gcf, 'Color', 'w')
        set(gca, 'FontSize',18)
        ylim([0.7 1.05])
        xlim([0 r_Ruu(find(Ruu<0.7,1,'first'))])
        vline(r_Ruu((find(R<i_tmp,1,'first'))),'k')
        lambda_autocorr=fitresult.a;
        % legend('','$R_{\widetilde{u}\widetilde{u}}(r)=1-\frac{r^2}{\lambda^2}$', 'Location', 'NorthEast' );
        legend('$R_{\widetilde{u}\widetilde{u}}(r)=1-\frac{r^2}{\lambda^2}$', 'Location', 'NorthEast' );
        title(['$\lambda$ = ' num2str(lambda_autocorr,'%1.4f') ' m'],'interpreter','latex')
        fig_setup
               
        trying = false;
    catch
        waitfor(msgbox('Please try again'));   
%         askInput({sprintf('Please try again by increasing the number of Bin''s used:')}, {num2str(increment_bin,'%1.0f')});
        if attemptcount>maxnumattempts
            waitfor(msgbox('Function failed.', 'Error','error'));   
            error('Function failed after %d times',maxnumattempts)
        end
        attemptcount=attemptcount+1;
    end  
end       

tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf);
delete(tmp_ui); 




%% numerical differentiation
%<(du/dx)^2>
du_dx           = mean((diff(daten)./(m_data*(1/Fs))).^2);
lambda_direct   = sqrt(uprimesquared/du_dx);

% Aronson and Löfdahl
tmp = sqrt((uprimesquared.*r.^2)./S_exp_2);
h(5) = figure;
plot(r,tmp,'MarkerEdgeColor','k','MarkerSize',6,'Marker','o', 'LineStyle','none','LineWidth',1.5)
hold on
xlabel('r/m', 'interpreter','latex')
xlabel('$r\ (m)$','interpreter','latex');
ylabel('$\sqrt{\frac{u''^2 r^2}{S^{2}(r)}}/m$', 'interpreter','latex')
ylabel('$\sqrt{\frac{u''^2 r^2}{S^{2}}}\ (m)$', 'interpreter','latex')
axis square
set(gcf, 'Color', 'w')
set(gca, 'FontSize',18)
xlim([0 25*low_freq_r])
datacursormode on
waitfor(msgbox('Use interactive data cursor to select the range of scales which will be used to estimate Taylor lenght scalee (see readme).'));

tmp_up      = (find(r>4*low_freq_r,1,'first'));
vline(r(tmp_up),'k')

trying              = true;
attemptcount        = 1;
while trying
    try
        [I,yy]  = ginput(1);
        % i_tmp       = str2double(inputdlg({'Fit Region rmin ='}));
        % tmp_up      = (find(r>i_tmp,1,'first'));
        i_tmp       = askInput({sprintf(['Fit Region r_max e.g. = ' ,num2str(20*low_freq_r,'%0.2f')])}, {num2str(I(1),'%.2f')});
        tmp_low     = (find(r<i_tmp,1,'last'));
        [xData, yData] = prepareCurveData(r(tmp_up:tmp_low),tmp(tmp_up:tmp_low));
        ft = fittype( 'poly1' );
        [fitresult, gof] = fit( xData, yData, ft );
        plot(linspace(0,1,1001),feval(fitresult,linspace(0,1,1001)),'r','LineWidth',2)
        fitresult.p1=0;
        plot(linspace(0,1,1001),feval(fitresult,linspace(0,1,1001)),'k','LineWidth',2,'LineStyle','-')
        vline(r(tmp_low),'k')
        A_L_lambda =   fitresult.p2;
        fig_setup
        legend off
        
        trying = false;
    catch
        waitfor(msgbox('Please try again'));   
        if attemptcount>maxnumattempts
            waitfor(msgbox('Function failed.', 'Error','error'));   
            error('Function failed after %d times',maxnumattempts)
        end
        attemptcount=attemptcount+1;
    end  
end               
        
tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf);
delete(tmp_ui);


% Dissipationspectrum
k_end                   = find(k<=low_freq_k, 1,'last');
lambda_spec_low_freq    = sqrt(trapz(k(1:k_end),Ek(1:k_end))/trapz(k(1:k_end),Ek(1:k_end).*k(1:k_end).^2));
lambda_spec_full        = sqrt(uprimesquared/trapz(k,Ek.*k.^2));
% trapz(k,2.*Ek.*k.^2) anstatt trapz(k,Ek.*k.^2)
% lambda_1=sqrt(uprimesquared/trapz(k,2.*Ek.*k.^2))
% lambda_2=sqrt(trapz(k,2.*Ek)/trapz(k,2.*Ek.*k.^2))

% take everthing of the specturm
% lambda_1=sqrt(uprimesquared/trapz(k,Ek.*k.^2));
% lambda_1=sqrt(uprimesquared/trapz(x_achse_k_spektrum_2pi(2:end),E_k.*x_achse_k_spektrum_2pi(2:end).^2))
% lambda_2=sqrt(trapz(k,Ek)/trapz(x_achse_k_spektrum_2pi(2:end),E_k.*x_achse_k_spektrum_2pi(2:end).^2))

% find Minimum at high frequencies
% k_start           = find(k>0.2*10^4, 1,'first');
% k_end             = find(k<=5*10^4, 1,'last');
% lambda_1_kurz     = sqrt(uprimesquared/trapz(k(1:find(spek==min(spek(k_start:k_end)))),Ek(1:find(spek==min(spek(k_start:k_end)))).*k(1:find(spek==min(spek(k_start:k_end)))).^2));

% lambda_4  = sqrt(uprimesquared/du_dx_4);
% lambda_6    = sqrt((10*kin_vis/m_data) * (uprimesquared/(uprimesquared/(m_data*(1/Fs)))));

% taylor_L_calc = [lambda_autocorr,lambda_direct,A_L_lambda,lambda_spec_low_freq,lambda_spec_full, lambda_S2, lambda_S3,lambda_zero_cros,lambda1_zero_cros];%% meter
taylor_L_calc = [lambda_autocorr,lambda_direct,A_L_lambda,lambda_spec_low_freq,lambda_spec_full,lambda_zero_cros,lambda1_zero_cros];%% meter


h(6) = figure;
plot(taylor_L_calc,'MarkerEdgeColor','k','MarkerSize',6,'Marker','o', 'LineStyle','none','LineWidth',1.5)
xlabel('method', 'interpreter','latex')
ylabel('$\lambda/m$', 'interpreter','latex')
ylabel('$\lambda\ (m)$', 'interpreter','latex')
axis square
set(gcf, 'Color', 'w')
set(gca, 'FontSize',18)
fig_setup
legend off
set(gca,'XTick',[1 2 3 4 5 6 7]);
datacursormode on
tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf);
delete(tmp_ui);
taylor_L = askInput({sprintf('Enter the Taylor length scale in meter. Median Taylor length scale in meter:')}, {num2str(nanmedian(taylor_L_calc),'%.5f')});



%% Dissipation from different methods
% epsi_1  = 15*kin_vis*trapz(k(1:k_end),Ek(1:k_end).*k(1:k_end).^2); %% From dissipation Spectrum untill cutoff frequency
% epsi_2  = 15*kin_vis*trapz(k,Ek.*k.^2); %%From entire dissipation Spectrum

% epsi_2  = trapz(k,2.*kin_vis.*Ek.*k.^2);
% epsi_3  = 15*kin_vis*du_dx_3;
% epsi_4  = 15*kin_vis*du_dx_4;
% epsi_5  = 15*kin_vis*du_dx_5;

% epsi_5 = (m_data/2) * (mean(diff(daten.^2)/(m_data*(1/Fs))));
% epsi_6 = (3*m_data/2) * (var(diff(daten)./(m_data*(1/Fs))));
% epsi_7 = 15*kin_vis*(uprimesquared/lambda_1^2); 
% epsi_7a= 30*kin_vis*(uprimesquared/lambda_1^2) %longitudinale microscale
% epsi_7b= 15*kin_vis*(uprimesquared/(lambda_1/sqrt(2))^2) %trasnversale microscale
% epsi_8 = m_data^3/int_L

% epsi_1/m_data
% var(daten)/(m_data*(1/Fs))
% std(daten)/(m_data*(1/Fs))
% epsi_calc   = [epsi_1, epsi_2, epsi_3, epsi_5, epsi_S2,epsi_S3];


h(7) = figure;
semilogx(r,eps_2,'LineWidth',2)
hold on
semilogx(r,-eps_3,'LineWidth',2)
xlim([min(r)*0.5 max(r)*2])
xlabel('r/m', 'interpreter','latex')
xlabel('$r\ (m)$','interpreter','latex');
ylabel('$\epsilon(r)/\frac{m^2}{s^3}$', 'interpreter','latex')
ylabel('$\epsilon\ (m^2/s^3)$','interpreter','latex');
set(gca,'FontSize',16)
set(gcf, 'Color', 'w')
fig_setup
axis square
% vline(r(I0),'r')                          
% vline(r(I1),'k')                          

plot(linspace(min(r)*0.5,max(r)*2,1000),eps_2(I0).*ones(1000,1),'k','LineWidth',2)
plot(linspace(min(r)*0.5,max(r)*2,1000),-eps_3(I1).*ones(1000,1),'k','LineWidth',2)
legend('$\frac{1}{r}\left[\frac{S^{2}(r)}{C_{2}}\right]^{3/2}$','$-\frac{5}{4}\left[\frac{S^{3}(r)}{r}\right]$','Location','south')

tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf);
delete(tmp_ui);

etaflow_S2  = ((kin_vis^3/epsi_S2)^(1/4));                  % Kolmogorov length scale in meters from epsi_S2 and hence from S_exp_2
etaflow_S3  = ((kin_vis^3/epsi_S3)^(1/4));                  % Kolmogorov length scale in meters from epsi_S3 and hence from S_exp_3



%% Spectrum-Fit
epsi_calc   = [epsi_S2,epsi_S3,(15*kin_vis*(uprimesquared/taylor_L^2))];

h(8) = figure;
loglog(x_achse_k_spektrum_2pi,E_k,'-','LineWidth',2)
hold on
a=find(x_achse_k_spektrum_2pi<=(2*pi*((1/(int_L/m_data))/(2*pi))/m_data), 1,'last');
b=find(x_achse_k_spektrum_2pi>=(2*pi*((1/(taylor_L/m_data))/(2*pi))/m_data), 1,'first');

% b=find(x_achse_k_spektrum_2pi>=(1)/diss_scale, 1,'first'); 
% a=1;
% b=length(x_achse_k_spektrum_2pi);
[xData, yData] = prepareCurveData( x_achse_k_spektrum_2pi(a:b), E_k(a:b) );
% [xData, yData] = prepareCurveData( x_achse_f_spektrum,spek_smooth);
% [xData, yData] = prepareCurveData(k, Ek );

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.DiffMinChange = 1e-12;
opts.Display = 'Off';
opts.MaxFunEvals = 1200;
opts.MaxIter = 1200;

ft = fittype( 'a*b^(2/3)*(x-d)^(-5/3)', 'independent', 'x', 'dependent', 'y' );
% opts.Lower = [0 median(epsi_calc) -Inf];
% opts.Upper = [Inf median(epsi_calc) Inf];
opts.Lower = [0 0 -Inf];
opts.Upper = [Inf Inf Inf];
opts.StartPoint = [0.5 nanmedian(epsi_calc) 0];
try 
    [fitresult, gof] = fit( xData, yData, ft, opts );
    ft2 = fittype( 'a*b^(2/3)*(x)^(-5/3)', 'independent', 'x', 'dependent', 'y' );
    % opts.Lower = [0 nanmedian(epsi_calc)];
    % opts.Upper = [Inf nanmedian(epsi_calc)];
    opts.Lower = [0 0];
    opts.Upper = [Inf Inf];
    opts.StartPoint = [0.5 nanmedian(epsi_calc)];
    [fitresult2, gof] = fit( xData, yData, ft2, opts );
    plot(x_achse_k_spektrum_2pi, feval(fitresult,x_achse_k_spektrum_2pi),'-','LineWidth',2)
    plot(x_achse_k_spektrum_2pi, feval(fitresult2,x_achse_k_spektrum_2pi),'-','LineWidth',2)
catch 
end
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('$k= / \frac{1}{m}$', 'interpreter','latex');
ylabel('$E(k)$ / $\frac{m^3}{s^2}$', 'interpreter','latex')
xlabel('$k\ (1/m)$','interpreter','latex');
ylabel('$\frac{E <u>}{2 \pi}\ (m^3/s^2)$', 'interpreter','latex')

xlim([min(x_achse_k_spektrum_2pi) 2*max(x_achse_k_spektrum_2pi)])
ylim([min(E_k) 2*max(E_k)])
axis square
% vline(x_achse_k_spektrum_2pi(a),'k')
% vline(x_achse_k_spektrum_2pi(b),'k')
% vline((1)/diss_scale,'k')
set(gcf, 'Color', 'w')
legend('$\frac{E <u>}{2 \pi}$','$C_k\langle \epsilon \rangle^{2/3}(k-k_0)^{-5/3}$','$C_k\langle \epsilon \rangle^{2/3}k^{-5/3}$','Location','south')
fig_setup
set(gca,'XTick',[10^-2 10^0  10^2  10^4]);
% fitresult;
% fitresult2;
tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf);
delete(tmp_ui);

% figure
% plot( fitresult, xData, yData );
% hold on
% plot( fitresult2, xData, yData );
% plot( fitresult, xData, yData );
% loglog(x_achse_k_spektrum_2pi,E_k)
% set(gca,'XScale','log');
% set(gca,'YScale','log');


epsi_calc   = [epsi_S2,epsi_S3, (15*kin_vis*(uprimesquared/taylor_L^2)), fitresult2.b,fitresult.b, Eps_zero_cros];

h(9) = figure;
plot(epsi_calc,'MarkerEdgeColor','k','MarkerSize',6,'Marker','o', 'LineStyle','none','LineWidth',1.5)
xlabel('method', 'interpreter','latex')
ylabel('$\langle \epsilon \rangle/\frac{m^2}{s^3}$', 'interpreter','latex')
ylabel('$\langle \epsilon \rangle\ (m^2/s^3)$', 'interpreter','latex')
axis square
set(gcf, 'Color', 'w')
set(gca, 'FontSize',18)
fig_setup
legend off
set(gca,'XTick',[1 2 3 4 5 6 7 ]);
datacursormode on
tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf);
delete(tmp_ui);
epsi = askInput({sprintf('Enter the energy dissipation rate in $m^2/s^3$. Median energy dissipation rate in $m^2/s^3$:')}, {num2str(nanmedian(epsi_calc),'%.5f')});




Re          = m_data*int_L/kin_vis;
% Re_lambda   = uprime*lambda_1_kurz_low_freq/kin_vis;
Re_lambda   = uprime*taylor_L/kin_vis; 
% Re_lambda = uprime*lambda_1/kin_vis;

% Calculate the dissipation length scale (micro scale) Roach Eqn 8
% Z=((2*pi^2)/(m_data^2*uprimesquared))*trapz((f(2:end).^2),spek(2:L/2+1)); % Use trapezoid rule for integration
% DissipationScale=sqrt(1/Z)

DissipationScale_1      = ((kin_vis^3)/epsi)^0.25;
% DissipationScale_2    = sqrt(1/(trapz(f,spek.*f.^2)*((2*pi^2)/(m_data^2*uprimesquared))))
% DissipationScale_3      = int_L/(Re^0.75);

% kolmogorov_calc         = [DissipationScale_1, etaflow_S2,etaflow_S3];%% meter
diss_scale         = DissipationScale_1;

%%Non dimentional dissipation constant
Ce      = 15*(int_L/taylor_L)/Re_lambda;
Ce_b    = epsi*int_L/uprime^3;

% Ce      = 15*(int_L/lambda_1_kurz_low_freq)/Re_lambda;
% Ce_1    = epsi_1*int_L/uprime^3;
% Ce_2    = epsi_2*int_L/uprime^3;
% Ce_3    = epsi_3*int_L/uprime^3;
% Ce_5    = epsi_5*int_L/uprime^3;

% Ce_calc =[Ce, Ce_1, Ce_2, Ce_3, Ce_5];
Ce_calc = [Ce, Ce_b];
Ce      = mean(Ce_calc);


%% Plot Spectrum
h(10) = figure;
set(gcf, 'Color', 'w')
yyaxis left
loglog(x_achse_k_spektrum_2pi,E_k,'-','LineWidth',2)
% hold on
% loglog(x_achse_k_spektrum_2pi,(x_achse_k_spektrum_2pi.^(-5/3).*ones(length(E_k)).*10),'-','LineWidth',1,'color','k')
% loglog(x_achse_k_spektrum_2pi,(x_achse_k_spektrum_2pi.^(-2).*ones(length(E_k)).*0.7.*10^2),'-','LineWidth',1,'color','k')
% loglog(x_achse_k_spektrum_2pi,(x_achse_k_spektrum_2pi.^(-3).*ones(length(E_k)).*0.6.*10^5),'-','LineWidth',1,'color','k')
% xlabel('$k=\frac{2 \pi f}{<u>} / m^{-1}$', 'interpreter','latex');
xlabel('$k= / \frac{1}{m}$', 'interpreter','latex');
% ylabel('$E(k)=\frac{E(f) <u>}{2 \pi}$ / $\frac{m^3}{s^2}$', 'interpreter','latex')
ylabel('$E(k)$ / $\frac{m^3}{s^2}$', 'interpreter','latex')
%xlabel('k');
%ylabel('E(k)')
xlabel('$k\ (1/m)$','interpreter','latex');
ylabel('$\frac{E <u>}{2 \pi}\ (m^3/s^2)$', 'interpreter','latex')
xlim([min(x_achse_k_spektrum_2pi) 2*max(x_achse_k_spektrum_2pi)])
% ylim([10^-10 1])
axis square
hold on
yyaxis right
% loglog(x_achse_k_spektrum_2pi,15.*kin_vis.*E_k.*(x_achse_k_spektrum_2pi).^2,'-','LineWidth',2)
loglog(x_achse_k_spektrum_2pi,E_k.*(x_achse_k_spektrum_2pi).^2,'-','LineWidth',2)
ylabel('E(k)$\cdot k^2$ / $\frac{m}{s^2}$', 'interpreter','latex')
ylabel('$\frac{E <u>}{2 \pi} k^2\ (m/s^2)$', 'interpreter','latex')

set(gca,'YScale','log');
set(gca,'FontSize',24)
% vline(low_freq_k,'r')
% vline((1/(int_L/m_data))/m_data,'k')
vline((2*pi*((1/(int_L/m_data))/(2*pi))/m_data),'k')
vline((2*pi*((1/(taylor_L/m_data))/(2*pi))/m_data),'k')
vline((2*pi*((1/(diss_scale/m_data))/(2*pi))/m_data),'k')
fig_setup
set(gca,'XTick',[10^-2 10^0  10^2  10^4]);
legend off



h(11) = figure;
set(gcf, 'Color', 'w')
yyaxis left
% Plot of Frequency Vs. Power Spectral Density(m^2/sec) with smoothing
loglog(x_achse_f_spektrum,spek_smooth,'-','LineWidth',2)
xlabel('f / Hz','interpreter','latex');
ylabel('E(f) / $\frac{m^2}{s}$', 'interpreter','latex')
xlabel('$f\ (Hz)$','interpreter','latex');
ylabel('$E\ (m^2/s)$', 'interpreter','latex')
set(gca,'FontSize',10)
set(gcf, 'Color', 'w')
xlim([min(x_achse_f_spektrum) max(x_achse_f_spektrum)]*2)
% xlim([min(x_achse_f_spektrum) 10^2])
grid on
axis square
datacursormode on
fig_setup
grid off
legend off
set(gca,'XTick',[10^-2 10^0  10^2  10^4]);
vline((1/(int_L/m_data))/(2*pi),'k')
vline((1/(taylor_L/m_data))/(2*pi),'k')
vline((1/(diss_scale/m_data))/(2*pi),'k')

yyaxis right
% Plot of Frequency Vs. Power Spectral Density(m^2/sec) with smoothing
loglog(x_achse_f_spektrum,spek_smooth.*x_achse_f_spektrum.^(5/3),'-','LineWidth',2)
ylabel('$E(f)f^{5/3}$ / $\frac{m^2}{s}$', 'interpreter','latex')
ylabel('$E\ f^{5/3}\ (m^2/s)$', 'interpreter','latex')
xlim([min(x_achse_f_spektrum) max(x_achse_f_spektrum)]*2)
ylim([min(spek_smooth.*x_achse_f_spektrum.^(5/3)) max(spek_smooth.*x_achse_f_spektrum.^(5/3))]*2)




% k Spektrum; Scales,r (m)
x_achse_spektrum_meter  = m_data./(2*pi*x_achse_f_spektrum);
Dr_avg_filter=(spek_smooth.*x_achse_f_spektrum.^2)./m_data;
Dk_avg_filter=15.*kin_vis.*E_k.*(x_achse_k_spektrum_2pi).^2;



figure;
set(gcf, 'Color', 'w')
subplot(2,2,1)
loglog(x_achse_f_spektrum,spek_smooth,'-','LineWidth',1)%% Plot of Frequency Vs. Power Spectral Density(m^2/sec) with smoothing & Filtering
xlabel('f / Hz','interpreter','latex');
ylabel('E(f) / $\frac{m^2}{s}$', 'interpreter','latex')
xlabel('$f\ (Hz)$','interpreter','latex');
ylabel('$E\ (m^2/s)$', 'interpreter','latex')
xlim([min(x_achse_f_spektrum) 2*max(x_achse_f_spektrum)])
ylim([min(spek_smooth) max(spek_smooth)]*2)
title('frequency domain','interpreter','latex')
set(gca,'FontSize',12)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
% plot([10^3 10^3],[10^-6 10^-2],'y = x ','-','LineWidth',1) 
% h = vline(10^3,'k','The Answer')
% text(245,spek_smooth(find(x_achse_f_spektrum>245,1,'first')),'\leftarrow sin(\pi) = 0','FontSize',16)
set(gca,'XTick',[10^-2 10^0  10^2  10^4]);
hold on
vline((1/(int_L/m_data))/(2*pi),'k')
vline((1/(taylor_L/m_data))/(2*pi),'k')
vline((1/(diss_scale/m_data))/(2*pi),'k')
txt = {'(a)'};
a=text(4,0.5,txt);
a.Units='normalized';

subplot(2,2,2)
loglog(x_achse_k_spektrum_2pi,E_k,'-','LineWidth',1)%% Plot of Wavenumber(k) Vs. Energy(k)(m^3/sec^2) with smoothing & Filtering
% loglog(x_achse_k_spektrum_2pi,x_achse_k_spektrum_2pi.^(-5/3).*(max(spek_smooth).*200),'--','color','k','LineWidth',2)
xlabel('$k=\frac{2 \pi f}{<u>} / m^{-1}$', 'interpreter','latex');
ylabel('$E(k)=\frac{E(f) <u>}{2 \pi}$ / $\frac{m^3}{s^2}$', 'interpreter','latex')
xlabel('$k=\frac{2 \pi f}{<u>}\ (1/m)$','interpreter','latex');
ylabel('$\frac{E <u>}{2 \pi}\ (m^3/s^2)$', 'interpreter','latex')
xlim([min(x_achse_k_spektrum_2pi) 2*max(x_achse_k_spektrum_2pi)])
ylim([min(E_k) max(E_k)]*2)
title('wave number domain','interpreter','latex')
set(gca,'FontSize',12)
set(gca,'XTick',[10^-2 10^0  10^2  10^4]);
hold on
vline((2*pi*((1/(int_L/m_data))/(2*pi))/m_data),'k')
vline((2*pi*((1/(taylor_L/m_data))/(2*pi))/m_data),'k')
vline((2*pi*((1/(diss_scale/m_data))/(2*pi))/m_data),'k')
txt = {'(b)'};
b=text(4,0.5,txt);
b.Units='normalized';

subplot4=subplot(2,2,3);
loglog(x_achse_spektrum_meter,Dr_avg_filter,'-','LineWidth',1)%% Plot of Scales,r(m) Vs. Energy(r)(m/sec^2) with smoothing & Filtering
xlabel('$r=\frac{<u>}{2 \pi f} / m$', 'interpreter','latex');
ylabel('$E(r)=\frac{E(f) f^2}{<u>}$ / $\frac{m}{s^2}$', 'interpreter','latex')
xlabel('$r=\frac{<u>}{2 \pi f}\ (m)$','interpreter','latex');
ylabel('$\frac{E\ f^2}{<u>}\ (m/s^2)$', 'interpreter','latex')
set(subplot4, 'Xdir', 'reverse')
xlim([0.5* min(x_achse_spektrum_meter) max(x_achse_spektrum_meter)])
ylim([min((spek_smooth.*x_achse_f_spektrum.^2)./m_data) max((spek_smooth.*x_achse_f_spektrum.^2)./m_data)]*2)
set(gca,'FontSize',12)
set(gca,'XTick',[10^-2 10^0  10^2  10^4]);
title('dissipation spectrum','interpreter','latex')
vline(int_L,'k')
vline(taylor_L,'k')
vline(diss_scale,'k')
txt = {'(c)'};
c=text(4,0.5,txt);
c.Units='normalized';



% Dissipation Spectrum
% kin_vis = 15.328737178409*10^-6;
subplot(2,2,4)
loglog(x_achse_k_spektrum_2pi,Dk_avg_filter,'-','LineWidth',1)%% Plot of Wavenumber(k) Vs. Dissipation(k)(m/sec^2) with smoothing & Filtering
xlabel('$k=\frac{2 \pi f}{<u>} / m^{-1}$', 'interpreter','latex');
ylabel('E(k)$\cdot k^2$ / $\frac{m}{s^2}$', 'interpreter','latex')
xlabel('$k=\frac{2 \pi f}{<u>}\ (1/m)$','interpreter','latex');
ylabel('$\frac{E <u>}{2 \pi}\ k^2 \ (m/s^2)$', 'interpreter','latex')
xlim([min(x_achse_k_spektrum_2pi) 2*max(x_achse_k_spektrum_2pi)])
ylim([min(15.*kin_vis.*E_k.*(x_achse_k_spektrum_2pi).^2) max(15.*kin_vis.*E_k.*(x_achse_k_spektrum_2pi).^2)]*2)
title('dissipation spectrum','interpreter','latex')
set(gca,'FontSize',12)
set(gca,'XTick',[10^-2 10^0  10^2  10^4]);
datacursormode on
hold on
vline((2*pi*((1/(int_L/m_data))/(2*pi))/m_data),'k')
vline((2*pi*((1/(taylor_L/m_data))/(2*pi))/m_data),'k')
vline((2*pi*((1/(diss_scale/m_data))/(2*pi))/m_data),'k')
txt = {'(d)'};
d=text(4,0.5,txt);
d.Units='normalized';

% font=22;
% a.FontSize = font;
% b.FontSize = font;
% c.FontSize = font;
% d.FontSize = font;

pos_txt=[-0.25   0.9];
a.Position=pos_txt;
b.Position=pos_txt;
c.Position=pos_txt;
d.Position=pos_txt;








if ischar(save_path)
    savefig(h,fullfile(save_path,append(save_name,'_','length_scales.fig')),'compact')
    for a = 1:length(h)
       exportgraphics(h(a),fullfile(save_path,append(save_name,'_',sprintf('length_scales_%d.png', a))))
    end
end
close all
end