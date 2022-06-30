%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arguments IN
% data_filter = filtered data 
%
% Arguments OUT
% condition = This input is for multipoint statistics
% tol = This input is for multipoint statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [condition,tol]=multi_point_condition(data_filter)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');

bin             = askInput({sprintf('Number of Bin''s used to perform a multi-point analysis?:')}, {num2str(6,'%1.0f')});
tmp             = linspace(min(data_filter),max(data_filter),bin+1);
[~,condition]   = hist(data_filter,bin);
tol             = mean(diff(condition))/2;

% percent    = 0.001;
% k_split    = 100/percent;
% tmp_data   = reshape(data_filter((1:floor(numel(data_filter)/k_split)*k_split)),[],k_split);
% x          = (1:k_split)*percent;


figure
% plot(x,tmp_data(1,:),'k')
plot((1:length(data_filter))/length(data_filter),data_filter,'k')
hold on
for i=1:length(tmp)
    yline(tmp(i),'--r');
end
for i=1:length(condition)
    yline(condition(i),'--b');
end
axis square 
xlabel('\% of data', 'interpreter','latex')
set(gca,'FontSize',18)
set(gcf, 'Color', 'w')
ylabel('$u / \frac{m}{s}$', 'interpreter','latex')
fig_setup
legend off

%     increment_bin = ceil((range(data)/std(data)*10));
%     if mod(increment_bin,2) == 0
%         increment_bin = increment_bin +1;
%     end
%     increment_bin = askInput({sprintf('Number of Bin''s used:')}, {num2str(increment_bin,'%1.0f')});

% min_events=ceil(0.00001*length(data));
% if min_events<400
%     min_events      = 400;
% end
% min_events = askInput({sprintf('Minimum number of events in a single bin:')}, {num2str(min_events ,'%1.0f')});
end