%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This function returns the filtered data using the low pass filter at the previously set frequency
% (using $\textbf{\textit{butter}}$ and $\textbf{\textit{filtfilt}}$ function of MATLAB). The 
% filtered data named as $\textbf{data\_filter}$ for all the further data post-processing. If the 
% filtering was negated in the previous step, $\textbf{data\_filter}$ and $\textbf{data}$ are equal.
% In addition, different representation/normalization of the energy spectrum density with respect to
% frequency $f$, scale $r$, wave number $k$ will also be plotted. 
%
% Arguments IN
% data = 1D array of data to which you would like to apply the low pass filter
% Fs = Acquisition/Sampling Frequency in Hz
% low_freq = Frequency at which you would like to use the low-pass filter in Hz
% kin_vis = Kinematic viscosity of fluid in m^2/sec
%
% Arguments OUT
% tmps = this is the filtered data which will be returned by this function
% x_achse_k_spektrum_2pi_filter = 2 pi k waven umbers(1/m) after filtering
% spek_smooth_filter = Frequency Spectrum after filtering
% E_k = Wavenumber Spectrum (m^3/sec^2) after filtering
% x_achse_spektrum_meter_filter = r Spektrum; Scales,r (m)  
% Dr_avg_filter = Dissipation Spectrum; scale domain; after filtering
% Dk_avg_filter = Dissipation Spectrum; wave number domain; after filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [tmps,x_achse_f_spektrum_filter,spek_smooth_filter,x_achse_k_spektrum_2pi_filter,E_k,x_achse_spektrum_meter_filter,Dr_avg_filter,Dk_avg_filter] = frequency_filter(data, Fs, low_freq,I,kin_vis,filter,save_path,save_name)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');

close
if mod(length(data),2)==1
    data=data(1:end-1); % This is to make length of data array even
end  

m_data=nanmean(data); % Mean of the data
if filter==1
    % butter is MATLAB function to design a filter of desired order and at desired cut-off frequency
    order=20;
    [b,a] = butter(order,low_freq/(Fs/2),'low');
    % tmps will be the array returned by this function
    tmps = filtfilt (b, a, data); % For filtfilt  see the MATLAB documentation

    % For a fixed cutoff frequency used in butter filter there is a threashold of order of butter 
    % filter at which we can apply the filter successfully.
    % Next loop tries to design a highest order of butter filter for a fixed cut-off low frequency
    % We repeat line 16,17 & 18 in the following while loop
    while (sum(isnan(tmps))>0 || (length(tmps)-nnz(tmps))>0 || var(tmps)>var(data))
        % sum(isnan(tmps))>0 ==> To avoid nan entries in filtered data
        % (length(tmps)-nnz(tmps))>0 ==> To avoid zeros in filterd data 
        % var(tmps)>var(data) ==>variance of the filtered data will always be
        % lesser than the variance of the unfilterd data
        clear tmps;
        order=order-1;
        clear b a;
        [b,a] = butter(order,low_freq/(Fs/2),'low');
        tmps = filtfilt (b, a, data); %tmps will be the array returned by this function
    end 
    % PSD of the fluctuations after filtering; mean value does not play a role in spectrum
    data_filter     = tmps-m_data;
else
    tmps            = data;
end

% PSD of the fluctuations before filtering; mean value does not play a role in spectrum
data          = data-m_data; 
L               = length(data); % length of the data

% spek              = abs(fft(data,L)).^2/L; % Power spectral density(PSD)
spek    = abs(fft(data,L)).^2/(L*Fs);        % Energy spectral density(ESD) using a fft
spek    = 2.*spek(2:L/2+1);                  % FFt will yield half number of unique points 

f       = Fs/2*linspace(0,1,L/2+1).';        % Nyquist frequency of the signal==Fs/2
f       = f(2:end);                          % Remove zero Hz component

% moving average with equally spaced frequency interval in log-space 
intervall = unique(round(logspace(0,log10(L/2),3001)),'stable');

% moving average with equally spaced frequency interval in lin-space 
% intervall=unique(round(linspace(1,L/2,plot_length)),'stable');

plot_length = length(intervall);             % Number of points for plotting averaged spectrum

% Initializing the arrays/preallocation
spek_smooth         = zeros(plot_length-2,1);
x_achse_f_spektrum  = zeros(plot_length-2,1);

% Averaging of spectrum and hence smoothing
for i=2:(plot_length-1)
    x_achse_f_spektrum(i-1,1) = mean(f(intervall(i-1):intervall(i+1)));
    spek_smooth(i-1,1)        = mean(spek(intervall(i-1):intervall(i+1)));
end
% ------------------------Windowing / averaging end ------------------------ 



%%%%%%%%% Filter data
%% PSD of the fluctuations after filtering
if filter==0
    data_filter     = data;    
end  
    
spek_filter              = abs(fft(data_filter,L)).^2/(L*Fs);   % Energy spectral density(ESD)
spek_filter              = 2.*spek_filter(1:L/2+1);             % FFt will yield half number of unique points
                                                                
f_filter = Fs/2*linspace(0,1,L/2+1).';                          % Nyquist frequency of the signal

% %moving average with equally spaced frequency interval in log-space
% intervall=unique(round(logspace(0,log10(L/2),plot_length)),'stable');
% plot_length=length(intervall);

spek_smooth_filter          = zeros(plot_length-2,1);
x_achse_f_spektrum_filter   = zeros(plot_length-2,1);

% Averaging of spectrum and hence smoothing
for i=2:(plot_length-1)
    spek_smooth_filter(i-1,1)           = mean(spek_filter(intervall(i-1):intervall(i+1)));
    x_achse_f_spektrum_filter(i-1,1)    = mean(f_filter(intervall(i-1):intervall(i+1)));
end
% ------------------------Windowing / averaging end ------------------------  


id          = round(logspace(0,log10(length(spek)),10^5));
id(end)     = length(spek);
f_plot      = f(id);
spek_plot   = spek(id);

%%%PLOTTING
% Plot of Frequency Vs. Energy spectral density(ESD)
h(1) = figure;
% Plot of Frequency Vs. Energy spectral density(ESD) without smoothing
loglog(f_plot,spek_plot)
hold on
%Plot of Frequency Vs. Power Spectral Density(m^2/sec) with smoothing
loglog(x_achse_f_spektrum,spek_smooth,'LineWidth',2)
% Plot of Frequency Vs. Power Spectral Density(m^2/sec) with smoothing & Filtering
loglog(x_achse_f_spektrum_filter,spek_smooth_filter,'LineWidth',2)
set(groot, 'defaultLegendInterpreter','latex');
% xlabel('f / Hz','interpreter','latex');
% ylabel('E(f) / $\frac{m^2}{s}$', 'interpreter','latex')
xlabel('$f\ (Hz)$','interpreter','latex');
ylabel('$E\ (m^2/s)$', 'interpreter','latex')
title(['Low-pass filter frequency =',num2str(low_freq,'%1.0f'), ' Hz'],'interpreter','latex')
set(gca,'FontSize',18)
set(gcf, 'Color', 'w')
xlim([min(x_achse_f_spektrum_filter) max(x_achse_f_spektrum_filter)]*2)
ylim([min(spek) max(spek)]*2)
axis square 
% vline(low_freq,'k')

% vline(m_data/int_L,'k')
% vline(m_data/taylor_L,'k')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
% title(['$ti=$',num2str(std(data)/m_data,'%1.3f')],'interpreter','latex') %% Turbulence Intensity of the flow
legend('raw','averaged','averaged \& low-pass filtered','Location','southwest')
set(gca,'XTick',[10^-4 10^-2 10^0  10^2  10^4]);


% savefig(fullfile(save_path,append(save_name,'_','spectrum_filter.fig')),'compact')
fig_setup

f_start            = find(f_plot>I(1), 1,'first');
f_end              = find(f_plot<=I(2), 1,'last');

[xData, yData] = prepareCurveData( f_plot(f_start:f_end), spek_plot(f_start:f_end));
ft = fittype( 'a*x^(-5/3)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
[fitresult, gof] = fit( xData, yData, ft, opts );

plot(f_plot,feval(fitresult,f_plot),'Color',[0 0 0],'LineWidth',2,'LineStyle','--')
legend('raw','averaged','averaged \& low-pass filtered','$f^{(-5/3)}$','Location','southwest')
ylim([min(spek_plot) max(spek_plot)]*4)

tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf); 
delete(tmp_ui);

if ischar(save_path)
   % savefig(fullfile(save_path,append(save_name,'_','spectrum_var_domain.fig')))
%    savefig(h,fullfile(save_path,append(save_name,'_','spectrum_filter.fig')),'compact')
%    for a = 1:length(h)
       exportgraphics(h(1),fullfile(save_path,append(save_name,'_',sprintf('spectrum_filter_%d.png', 1))))
%        export_fig(h(a), fullfile(save_path,append(save_name,'_',sprintf('spectrum_filter_%d.png', a))));
%    end
end


% Wavenumber Spectrum (m^3/sec^2)
E_k = ((spek_smooth_filter.*m_data)./(2*pi)); 
% 2-pi - k Spektrum; Wavenumbers(1/m)
x_achse_k_spektrum_2pi_filter = (2.*pi.*x_achse_f_spektrum_filter)./m_data;
% k Spektrum; Scales,r (m)
x_achse_spektrum_meter_filter = m_data./(2*pi*x_achse_f_spektrum_filter); 

Dr_avg_filter=(spek_smooth_filter.*x_achse_f_spektrum_filter.^2)./m_data;
Dk_avg_filter=15.*kin_vis.*E_k.*(x_achse_k_spektrum_2pi_filter).^2;


h(2) = figure;
set(gcf, 'Color', 'w')
subplot(2,2,1)
loglog(x_achse_f_spektrum_filter,spek_smooth_filter,'-','LineWidth',1)%% Plot of Frequency Vs. Power Spectral Density(m^2/sec) with smoothing & Filtering
% xlabel('f / Hz','interpreter','latex');
% ylabel('E(f) / $\frac{m^2}{s}$', 'interpreter','latex')
xlabel('$f\ (Hz)$','interpreter','latex');
ylabel('$E\ (m^2/s)$', 'interpreter','latex')
xlim([min(x_achse_f_spektrum_filter) 2*max(x_achse_f_spektrum_filter)])
ylim([min(spek_smooth_filter) max(spek_smooth_filter)]*2)
title('frequency domain','interpreter','latex')
set(gca,'FontSize',12)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
% plot([10^3 10^3],[10^-6 10^-2],'y = x ','-','LineWidth',1) 
% h = vline(10^3,'k','The Answer')
% text(245,spek_smooth(find(x_achse_f_spektrum>245,1,'first')),'\leftarrow sin(\pi) = 0','FontSize',16)
set(gca,'XTick',sort([10^-2 10^0  10^2 10^4 low_freq]));
hold on
vline(low_freq,'k')

txt = {'(a)'};
a=text(4,0.5,txt);
a.Units='normalized';

subplot(2,2,2)
loglog(x_achse_k_spektrum_2pi_filter,E_k,'-','LineWidth',1)%% Plot of Wavenumber(k) Vs. Energy(k)(m^3/sec^2) with smoothing & Filtering
% loglog(x_achse_k_spektrum_2pi_filter,x_achse_k_spektrum_2pi_filter.^(-5/3).*(max(spek_smooth).*200),'--','color','k','LineWidth',2)
% xlabel('$k=\frac{2 \pi f}{<u>} / m^{-1}$', 'interpreter','latex');
% ylabel('$E(k)=\frac{E(f) <u>}{2 \pi}$ / $\frac{m^3}{s^2}$', 'interpreter','latex')
xlabel('$k=\frac{2 \pi f}{<u>}\ (1/m)$','interpreter','latex');
ylabel('$\frac{E <u>}{2 \pi}\ (m^3/s^2)$', 'interpreter','latex')
xlim([min(x_achse_k_spektrum_2pi_filter) 2*max(x_achse_k_spektrum_2pi_filter)])
ylim([min(E_k) max(E_k)]*2)
title('wave number domain','interpreter','latex')
set(gca,'FontSize',12)
set(gca,'XTick',sort([10^-2 10^0  10^2 10^4 round(2*pi*low_freq/m_data)]));
set(gca,'XTick',sort([10^-2 10^0  10^2 round(2*pi*low_freq/m_data)]));
hold on
vline(2*pi*low_freq/m_data,'k')
txt = {'(b)'};
b=text(4,0.5,txt);
b.Units='normalized';

subplot4=subplot(2,2,3);
loglog(x_achse_spektrum_meter_filter,Dr_avg_filter,'-','LineWidth',1)%% Plot of Scales,r(m) Vs. Energy(r)(m/sec^2) with smoothing & Filtering
% xlabel('$r=\frac{<u>}{2 \pi f} / m$', 'interpreter','latex');
% ylabel('$E(r)=\frac{E(f) f^2}{<u>}$ / $\frac{m}{s^2}$', 'interpreter','latex')
xlabel('$r=\frac{<u>}{2 \pi f}\ (m)$','interpreter','latex');
ylabel('$\frac{E\ f^2}{<u>}\ (m/s^2)$', 'interpreter','latex')
set(subplot4, 'Xdir', 'reverse')
xlim([0.5* min(x_achse_spektrum_meter_filter) max(x_achse_spektrum_meter_filter)])
ylim([min((spek_smooth_filter.*x_achse_f_spektrum_filter.^2)./m_data) max((spek_smooth_filter.*x_achse_f_spektrum_filter.^2)./m_data)]*2)
set(gca,'FontSize',12)
set(gca,'XTick',sort([10^-2 10^0  10^2 10^4 m_data/(2*pi*low_freq)]));
title('dissipation spectrum','interpreter','latex')
hold on
vline(m_data/(2*pi*low_freq),'k')
txt = {'(c)'};
c=text(4,0.5,txt);
c.Units='normalized';


% Dissipation Spectrum
% kin_vis = 15.328737178409*10^-6;
subplot(2,2,4)
loglog(x_achse_k_spektrum_2pi_filter,Dk_avg_filter,'-','LineWidth',1)%% Plot of Wavenumber(k) Vs. Dissipation(k)(m/sec^2) with smoothing & Filtering
% xlabel('$k=\frac{2 \pi f}{<u>} / m^{-1}$', 'interpreter','latex');
% ylabel('E(k)$\cdot k^2$ / $\frac{m}{s^2}$', 'interpreter','latex')
xlabel('$k=\frac{2 \pi f}{<u>}\ (1/m)$','interpreter','latex');
ylabel('$\frac{E <u>}{2 \pi}\ k^2 \ (m/s^2)$', 'interpreter','latex')

xlim([min(x_achse_k_spektrum_2pi_filter) 2*max(x_achse_k_spektrum_2pi_filter)])
ylim([min(15.*kin_vis.*E_k.*(x_achse_k_spektrum_2pi_filter).^2) max(15.*kin_vis.*E_k.*(x_achse_k_spektrum_2pi_filter).^2)]*2)
title('dissipation spectrum','interpreter','latex')
set(gca,'FontSize',12)
set(gca,'XTick',sort([10^-2 10^0  10^2 10^4 round(2*pi*low_freq/m_data)]));
set(gca,'XTick',sort([10^-2 10^0  10^2 round(2*pi*low_freq/m_data)]));
hold on
vline(2*pi*low_freq/m_data,'k')
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

datacursormode on


tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf); 
delete(tmp_ui);

if ischar(save_path)
   % savefig(fullfile(save_path,append(save_name,'_','spectrum_var_domain.fig')))
   savefig(h,fullfile(save_path,append(save_name,'_','spectrum_filter.fig')),'compact')
%    for a = 1:length(h)
       exportgraphics(h(2),fullfile(save_path,append(save_name,'_',sprintf('spectrum_filter_%d.png', 2))))
%        export_fig(h(a), fullfile(save_path,append(save_name,'_',sprintf('spectrum_filter_%d.png', a))));
%    end
end
close all
end