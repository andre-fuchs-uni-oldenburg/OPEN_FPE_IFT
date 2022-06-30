%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This function plots the mean, standard deviation, skewness and kurtosis of sections of a length of
% 5\% of the data to check the stationarity of the data. In the title of the figure the number of 
% nan's in the dataset and the turbulence intensity is printed.
% 
% Arguments IN
% data = 1D array of data
% 
% Arguments OUT
% data = 1D array of data without nan's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data]=plot_stationarity(data,save_path,save_name)
%% Include a latex function for plotting all figure properties in latex format
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');

percent    = 5;
k_split    = 100/percent;
if numel(data)<k_split
    k_split=numel(data);
end
tmp        = reshape(data((1: floor(numel(data)/k_split)*k_split)),[],k_split);


x          = (1:k_split)*percent;

h(1) = figure;
% tiledlayout(1,2);
% nexttile;
subplot(1,2,1)
plot(x,nanmean(tmp),'LineWidth',2,'MarkerSize',8,'Marker','o','LineStyle','none')
hold on
plot(x,nanstd(tmp),'LineWidth',2,'MarkerSize',8,'Marker','o','LineStyle','none')
plot(x,skewness(tmp),'LineWidth',2,'MarkerSize',8,'Marker','o','LineStyle','none')
plot(x,kurtosis(tmp),'LineWidth',2,'MarkerSize',8,'Marker','o','LineStyle','none')
% plot(x,kurtosis(data).*ones(numel(x),1),'k','LineStyle','--','LineWidth',2);
xlim([0 105])
axis square 
% xlabel('\% of data', 'interpreter','latex')
xlabel('percentage of data', 'interpreter','latex')
% ylabel('m/s', 'interpreter','latex')
set(gca,'FontSize',18)
set(gcf, 'Color', 'w')
% legend('$\left<u\right>\ (m/s)$','std/$\frac{m}{s}$','skewnss','kurtosis','Location','east')
legend('mean $(m/s)$','std $(m/s)$','skewnss','kurtosis','Location','east')
% legend('mean/$\frac{m}{s}$','std/$\frac{m}{s}$','skewnss','kurtosis','Location','east')
title(['Number of nans = ',num2str(sum(isnan(data)),'%i')],'interpreter','latex')
fig_setup
% if 0.8*min([min(nanmean(tmp)) min(nanstd(tmp)) min(skewness(tmp)) min(kurtosis(tmp))])<0
%     ylim([ 1.5*min([min(nanmean(tmp)) min(nanstd(tmp)) min(skewness(tmp)) min(kurtosis(tmp))])...
%            1.2*max([max(nanmean(tmp)) max(nanstd(tmp)) max(skewness(tmp)) max(kurtosis(tmp))]) ]);
% else
%     ylim([ 0.5*min([min(nanmean(tmp)) min(nanstd(tmp)) min(skewness(tmp)) min(kurtosis(tmp))])...
%            1.2*max([max(nanmean(tmp)) max(nanstd(tmp)) max(skewness(tmp)) max(kurtosis(tmp))]) ]);
% end 
% ylim([-0.5 3.5])
txt = {'(a)'};
a=text(4,0.5,txt);
a.FontSize = 22;
a.Units='normalized';


% h(2) = figure; 
% nexttile;
percent    = 0.001;
k_split    = 100/percent;

if numel(data)<k_split
    k_split=numel(data);
end
    
tmp        = reshape(data((1: floor(numel(data)/k_split)*k_split)),[],k_split);
x          = (1:k_split)*percent;

subplot(1,2,2)
% plot((1:length(data))*100/length(data),data,'k')
% hold on
plot(x,tmp(1,:),'k')
axis square 
% xlabel('\% of data', 'interpreter','latex')
xlabel('percentage of data', 'interpreter','latex')
set(gca,'FontSize',18)
set(gcf, 'Color', 'w')
ylabel('$u\ (m/s)$', 'interpreter','latex')
% ylabel('$u / \frac{m}{s}$', 'interpreter','latex')
title(['Ti = ',num2str(nanstd(data)/nanmean(data)*100,'%1.2f'), '\,\%'],'interpreter','latex')
fig_setup
legend off
txt = {'(b)'};
b=text(4,0.5,txt);
b.FontSize = 22;
b.Units='normalized';
pos_txt=[-0.25   0.9 0];
a.Position=pos_txt;
b.Position=pos_txt;
% ylim([0 4])

% set(gcf,'Units','normalized','Position',[0.0926 0.6229 0.3598 0.2917])
uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf);

if ischar(save_path)
   savefig(h,fullfile(save_path,append(save_name,'_','stationarity.fig')),'compact')
   for a = 1:length(h)
       exportgraphics(h(a),fullfile(save_path,append(save_name,'_',sprintf('stationarity_%d.png', a))))
   end
end


if sum(isnan(data))>0
waitfor(msgbox('The time series contains missing entries. Select a method to fill these entries (see readme).','Success','help'))
List = { 'Nearest non-missing value','Linear interpolation of neighboring, non-missing values','Piecewise cubic spline interpolation','moving average','moving median'}; 
[idx, tf] = listdlg('ListString', List,'SelectionMode', 'Single', 'PromptString', 'Fill missing entries using method:', 'Initialvalue', 3);

if idx==1
	data = fillmissing(data,'nearest');
elseif idx==2
	data= fillmissing(data,'linear');
elseif idx==3 
	data= fillmissing(data,'spline');
elseif idx==4
	data= fillmissing(data,'movmean',str2double(inputdlg({'Window of length in sample:'}))); 
elseif idx==5
	data= fillmissing(data,'movmedian',str2double(inputdlg({'Window of length in sample:'}))); 
end
%     delete_nan     = askYesno('Do you want to delete nans?', 'Yes');
%     if delete_nan==1
%         data=data(~isnan(data));
%     end
end
close
end