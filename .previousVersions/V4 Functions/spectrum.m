%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This function calculates the energy spectral density (ESD) of the hot wire signal using 
% $\textbf{\textit{fft}}$ function of MATLAB. Also, the ESD with and without averaging 
% (moving average with equally spaced frequency interval in log-space) as a function of frequency 
% will be plotted.
%
% Arguments IN
% data = 1D array for which you would like to plot power spectram density(PSD)
% 
% Arguments OUT
% f = Frequency without smoothing
% spek = Energy spectral density(ESD) without smoothing
% x_achse_f_spektrum = Frequency with smoothing
% spek_smooth = Energy spectral density(ESD) with smoothing
% spek_power = Power spectral density(PSD) with smoothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function  [f, spek, x_achse_f_spektrum, spek_smooth,spek_power]=spectrum(data,Fs,increment_bin)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');

% PSD of the fluctuations; mean value does not play a role in spectrum
data           = data-mean(data); 

if mod(length(data),2)==1
    data=data(1:end-1); % This is to make length of data array even
end    

L                 = length(data); 

% please have a look at 'fft' function of matlab in the MATLAB documentation
% spek              = abs(fft(data,L)).^2/L;        % Power spectral density(PSD)
spek              = abs(fft(data,L)).^2/(L*Fs);     % Energy spectral density(ESD) using a fft

spek              = 2.*spek(2:L/2+1);               % FFt will yield half number of unique points

f                 = Fs/2*linspace(0,1,L/2+1).';     % Nyquist frequency of the signal==Fs/2
f                 = f(2:end);                       % Remove zero Hz component 


% moving average with equally spaced frequency interval in log-space 
plot_length         = increment_bin*10;
if mod(plot_length,2) == 0
    plot_length = plot_length +1;
end

intervall=unique(round(logspace(0,log10(L/2),plot_length)),'stable');

% moving average with equally spaced frequency interval in lin-space 
% intervall=unique(round(linspace(1,L/2,plot_length)),'stable');

plot_length=length(intervall); % Number of points for plotting averaged spectrum

% Initializing the arrays/preallocation
spek_smooth=zeros(plot_length-2,1);
x_achse_f_spektrum=zeros(plot_length-2,1);

% Averaging of spectrum and hence smoothing
for i=2:(plot_length-1)
    x_achse_f_spektrum(i-1,1) = mean(f(intervall(i-1):intervall(i+1)));
    spek_smooth(i-1,1)        = mean(spek(intervall(i-1):intervall(i+1)));
end
spek_power = spek_smooth.*Fs;

% ------------------------Windowing / averaging end ------------------------  
figure
% Plot of Frequency Vs. Energy spectral density(ESD) without smoothing
loglog(f,spek)
hold on
% Plot of Frequency Vs. Energy spectral density(ESD) with smoothing
loglog(x_achse_f_spektrum,spek_smooth,'-','LineWidth',2)
xlabel('f / Hz','interpreter','latex');
ylabel('E(f) / $\frac{m^2}{s}$', 'interpreter','latex')
title('Energy density spectrum','Interpreter','latex')
set(gca,'FontSize',18)
set(gcf, 'Color', 'w')
xlim([min(x_achse_f_spektrum) max(x_achse_f_spektrum)]*2)
legend('raw','averaged')
axis square
datacursormode on
set(gca,'XTick',[10^-2 10^0  10^2  10^4]);
fig_setup
end