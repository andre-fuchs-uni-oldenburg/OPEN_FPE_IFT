%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots the probability density function (PDF) of the data with the specified number 
% of bins. It also plots the Gaussian distribution which has the same standard deviation and mean as
% of the data. In the title of the figure the range of the data (difference between the maximum and
% minimum values of sample data), the skewness and flatness of the data is printed.
% COMMENT: hist is a standard matlab function==>Look in to standard MATLAB documentation
%
% Arguments IN
% data  = 1D array for which you would like to plot the probability density function(PDF)
% nbins = The number of bins in which you would like to divide your data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function plot_pdf(data,nbins,save_path,save_name)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');

[fuL, uL] = hist(data,nbins); % fuL is the counters & uL is the center of the bin
% Converting fuL to puL in order to change counters to probability density fuction
% Probability = fuL./nansum(fuL)==> Sum of Probability array should be 1.0<== Definition of Probability 
% puL = [Probability./mean(diff(uL)] ==> Probability density function
% ==>trapz(uL(1:end),puL) must be 1.0 <== Definition of Probability density function
puL = fuL / ( nansum(fuL) * mean(diff(uL)) );

% Plotting the PDF
h(1)=figure;
% Actual PDF of the data
plot(uL,puL,'MarkerEdgeColor','k','MarkerSize',6,'Marker','o', 'LineStyle','none')
hold on
axis square   
set(gca,'yScale','log') 
set(gca,'FontSize',22)
set(gcf, 'Color', 'w')
xlim([min(data) max(data)])
ylim([0.0000001*max(puL) 5*max(puL)])
xPos = mean(data);
% Plot a vertical line at x=mean(data) to visualize a flatness or skewness
plot([xPos xPos], get(gca,'ylim'),'LineWidth',2,'LineStyle','--','Color','k'); 
% Adapts to y limits of current axes
% Plot a gaussian distribution PDF with same mean and same standard deviation as of our original
% data for comparison.
gaus_bin    = linspace(min(data),max(data),nbins);
norm1       = pdf(fitdist(data,'Normal'),gaus_bin);
plot(gaus_bin,norm1,'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176])
% ylabel('$p(u)$','interpreter','latex')
% xlabel('$u / \frac{m}{s}$', 'interpreter','latex')
xlabel('$u\ (m/s)$', 'interpreter','latex')
ylabel('PDF','interpreter','latex')

% title(['range = ',num2str(range(data),'%1.2f'), ';',...
%     ' skewness = ',num2str(skewness(data),'%1.2f'), ';',...
%     ' kurtosis = ',num2str(kurtosis(data),'%1.2f'),...
%     ],'interpreter','latex')
title(['range = ',num2str(range(data),'%1.2f'), ';',...
          ' S = ',num2str(skewness(data),'%1.2f'), ';',...
          ' K = ',num2str(kurtosis(data),'%1.2f'),...
    ],'interpreter','latex')
legend1 = legend('Original data','Mean of data','Gaussian Distribution','Location','south');
set(legend1,'FontSize',18);
fig_setup
uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf); 
if ischar(save_path)
    savefig(h,fullfile(save_path,append(save_name,'_','pdf.fig')),'compact')
    for a = 1:length(h)
        exportgraphics(h(a),fullfile(save_path,append(save_name,'_',sprintf('pdf_%d.png', a))))
    end
end
close
end