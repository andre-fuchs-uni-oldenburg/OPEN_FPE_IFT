%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot: comparision co_KM_opti and co_IFT_opti dij coefficients
% Arguments IN
% data = filtered data 
% stp = index used to cut v(x) in non-overlapping pieces of length equal to integral scale
% r =  scale vector from the start to end of the cascade trajectory
% rind = index-vector from the start to end of the cascade trajectory
% co_KM_opti = Coefficients $d_{ij}(r)$ of the optimized Kramers-Moyal coefficients towards the best 
%      Fokker-Planck equation to reproduce the conditional PDF's using the surface fits
% co_IFT_opti = Coefficients $d_{ij}(r)$ of the optimized Kramers-Moyal coefficients towards the 
%               integral fluctuation theorem using the surface fits
% evaluated = a modified/updated struct 'evaluated' array with in function 'KM_STP_optimization'
% Fs = Acquisition/Sampling Frequency in Hz
% taylor_L = Taylor length scale in meters
% lb = Lower bound
% ub = Upper bound
% t = length(r)
% history = details of the optimization for every single iteration 
% index = index of minimized error function ==> history
% d11,d20,d21,d22 = Coefficients $d_{ij}(r)$ of co_IFT_opti
% m_data = mean of the data
% markov = markov length in number of samples
% trajec = 1 ==> The start/end of the cascade trajectory will be adjusted.
% dr_ind = Separation of scales/step increment (in samples) referred to the sequence from large to small scales in the cascade trajectory
% z = 3 ==> Independent trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OPTI_IFT_plot(data,stp,r,rind,co_KM_opti,co_IFT_opti,evaluated,Fs,taylor_L,lb,ub,t,history,index,d11,d20,d21,d22,m_data,markov,trajec,dr_ind,z,norm_ur,norm_r,save_path,save_name)
tmp_size    = ceil(size(data,1)/stp)-1;
[u]         = create_us(data,r,rind,stp);

%% co_KM_opti
Sm = nan( tmp_size, 1 );
Ds = nan( tmp_size, 1 );
DS = nan( tmp_size, 1 );

[Sm_tmp,Ds_tmp,DS_tmp,~,~,~]    = calcDS(z,u,r.',co_KM_opti,markov,taylor_L,Fs,m_data,dr_ind);
Sm(1:sum(~isnan(DS_tmp),1),1)   = Sm_tmp(~isnan(DS_tmp));
Ds(1:sum(~isnan(DS_tmp),1),1)   = Ds_tmp(~isnan(DS_tmp));
DS(1:sum(~isnan(DS_tmp),1),1)   = DS_tmp(~isnan(DS_tmp));

DS(isinf(DS))           = nan;
Sm                      = reshape(Sm,[],1);
Sm                      = Sm(~isnan(DS));
Ds                      = reshape(Ds,[],1);
Ds                      = Ds(~isnan(DS));
DS                      = reshape(DS,[],1);
DS                      = DS(~isnan(DS));

N = size(DS,1); % number of trajectories u(r)
NN = 20;
Nvec_alt = round(logspace(1,log10(N),NN));
FTvec_alt = nan(1,NN);
FTvec_err_alt = nan(1,NN);
for n=1:NN
    FTvec_alt(n) = nanmean( exp( -DS(1:Nvec_alt(n)) ) );
%     FTvec(n) = median( exp( -DS(1:Nvec(n)) ) );
    FTvec_err_alt(n)=nanstd( exp( -DS(1:Nvec_alt(n)) ) )/sqrt(Nvec_alt(n));
end




%% co_IFT_opti
Sm = nan( tmp_size, 1 );
Ds = nan( tmp_size, 1 );
DS = nan( tmp_size, 1 );
[Sm_tmp,Ds_tmp,DS_tmp,~,~,~]    = calcDS(z,u,r.',co_IFT_opti,markov,taylor_L,Fs,m_data,dr_ind);
Sm(1:sum(~isnan(DS_tmp),1),1)   = Sm_tmp(~isnan(DS_tmp));
Ds(1:sum(~isnan(DS_tmp),1),1)   = Ds_tmp(~isnan(DS_tmp));
DS(1:sum(~isnan(DS_tmp),1),1)   = DS_tmp(~isnan(DS_tmp));

DS(isinf(DS))           = nan;
Sm                      = reshape(Sm,[],1);
Sm                      = Sm(~isnan(DS));
Ds                      = reshape(Ds,[],1);
Ds                      = Ds(~isnan(DS));
DS                      = reshape(DS,[],1);
DS                      = DS(~isnan(DS));


N = size(DS,1); % number of trajectories u(r)
NN = 20;
Nvec = round(logspace(1,log10(N),NN));
FTvec = nan(1,NN);
FTvec_err = nan(1,NN);
for n=1:NN
    FTvec(n) = nanmean( exp( -DS(1:Nvec(n)) ) );
%     FTvec(n) = median( exp( -DS(1:Nvec(n)) ) );
    FTvec_err(n)=nanstd( exp( -DS(1:Nvec(n)) ) )/sqrt(Nvec(n));
end


h(1) =  figure;
errorbar(Nvec_alt,FTvec_alt,FTvec_err_alt,'MarkerSize',8,'Marker','*','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',3,'LineStyle','none')
hold on
errorbar(Nvec,FTvec,FTvec_err,'MarkerSize',8,'Marker','*','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',3,'LineStyle','none')
plot(Nvec,ones(numel(Nvec),1),'k','LineStyle','--','LineWidth',2);
axis square      
set(gca, 'FontSize',18)
set(gcf, 'Color', 'w')
xlabel('$N$','interpreter','latex')
ylabel('$\langle e^{\mathrm{-}\Delta S_{tot}} \rangle_N$','interpreter','latex')
set(gca,'xScale','log') 
% set(gca,'yScale','log') 
xlim([min(Nvec) max(Nvec)])
% set(gca,'XTick',[10 1000 100000])
% ylim([0.1 1.1])
title(['$\langle e^{\mathrm{-}\Delta S_{tot}} \rangle_{max(N)}=$',num2str(nanmean( exp( -DS )),'%1.2f')],'interpreter','latex')
fig_setup
legend off


h(2) = figure;
subplot(1,2,1)
errorbar(Nvec,abs(1-FTvec),FTvec_err,'MarkerSize',8,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',3,'LineStyle','none','Color',[0 0 0])
% plot(Nvec,abs(1-FTvec),'MarkerSize',8,'Marker','*','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',4,'LineStyle','none')
hold on 
[xData, yData, weights] = prepareCurveData( Nvec, abs(1-FTvec), FTvec_err );
% Set up fittype and options.
ft = fittype( 'power1' );
ft_2 = fittype( 'a*x^-0.5', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Weights = weights;
% [fit_IFT, gof] = fit( xData, yData, ft, opts );
% plot(Nvec,fit_IFT.a*Nvec.^fit_IFT.b,'k','LineWidth',4)

[fit_IFT_2, gof] = fit( xData, yData, ft_2, opts );
plot(Nvec,fit_IFT_2.a*Nvec.^-0.5,'r','LineWidth',4)
% plot(Nvec,abs(1-FTvec),'MarkerSize',8,'Marker','*','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',4,'LineStyle','none')
axis square      
set(gca, 'FontSize',18)
set(gcf, 'Color', 'w')
xlabel('$N$','interpreter','latex')
ylabel('$|1-\langle e^{\mathrm{-}\Delta S_{tot}} \rangle_N|$','interpreter','latex')
set(gca,'xScale','log') 
% set(gca,'yScale','log') 
xlim([min(Nvec) max(Nvec)])
% set(gca,'XTick',[10 1000 100000])
% ylim([-0.3 1.2])
% title(['$\langle e^{\mathrm{-}\Delta S_{tot}} \rangle_{max(N)}=$',num2str(nanmean( exp( -DS )),'%1.4f')],'interpreter','latex')
plot(Nvec,0.*ones(numel(Nvec),1),'k','LineStyle','--','LineWidth',2);
legend({'exp','$aN^{-0.5}$'},'Interpreter','latex','Location','northeast')


% figure
subplot(1,2,2)
% errorbar(Nvec,abs(1-FTvec),FTvec_err,'MarkerSize',8,'Marker','*','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',4,'LineStyle','none')
plot(Nvec,abs(1-FTvec),'MarkerSize',8,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',3,'LineStyle','none')
hold on 
[xData, yData, weights] = prepareCurveData( Nvec, abs(1-FTvec), FTvec_err );
% Set up fittype and options.
ft = fittype( 'power1' );
ft_2 = fittype( 'a*x^-0.5', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Weights = weights;
% [fit_IFT, gof] = fit( xData, yData, ft, opts );
% plot(Nvec,fit_IFT.a*Nvec.^fit_IFT.b,'k','LineWidth',4)

[fit_IFT_2, gof] = fit( xData, yData, ft_2, opts );
plot(Nvec,fit_IFT_2.a*Nvec.^-0.5,'r','LineWidth',4)
% plot(Nvec,abs(1-FTvec),'MarkerSize',8,'Marker','*','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',4,'LineStyle','none')
axis square      
set(gca, 'FontSize',18)
set(gcf, 'Color', 'w')
xlabel('$N$','interpreter','latex')
ylabel('$|1-\langle e^{\mathrm{-}\Delta S_{tot}} \rangle_N|$','interpreter','latex')
set(gca,'xScale','log') 
set(gca,'yScale','log') 
xlim([min(Nvec) max(Nvec)])
% set(gca,'XTick',[10 1000 100000])
% ylim([-0.3 1.2])
% title(['$\langle e^{\mathrm{-}\Delta S_{tot}} \rangle_{max(N)}=$',num2str(nanmean( exp( -DS )),'%1.4f')],'interpreter','latex')
legend({'exp','$aN^{-0.5}$'},'Interpreter','latex','Location','northeast')

   

aufl=size(evaluated(1).D1,2);
xx=nan(length(evaluated),aufl);
y=nan(length(evaluated),aufl);
z1=nan(length(evaluated),aufl);
z1e=nan(length(evaluated),aufl);
z2=nan(length(evaluated),aufl);
z2e=nan(length(evaluated),aufl);
z1_opti=nan(length(evaluated),aufl);
z2_opti=nan(length(evaluated),aufl);
z4=nan(length(evaluated),aufl);

% for i=find((evaluated.r)/taylor_L>4):length(evaluated) 
for i=1:length(evaluated)
    xx(i,1:sum(evaluated(i).x_bin_not_nan))=evaluated(i).y_mean_bin(evaluated(i).x_bin_not_nan);
    if norm_r==1
    y(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).r./taylor_L);
    else
    y(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).r);   
    end 
    z1(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).D1(evaluated(i).x_bin_not_nan));
    z1e(i,1:sum(evaluated(i).x_bin_not_nan))=abs(1./(evaluated(i).eD1(evaluated(i).x_bin_not_nan)));

    z2(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).D2(evaluated(i).x_bin_not_nan));
    z2e(i,1:sum(evaluated(i).x_bin_not_nan))=abs(1./(evaluated(i).eD2(evaluated(i).x_bin_not_nan)));

    z1_opti(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).D1_opti);
    z2_opti(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).D2_opti);

    z4(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).D4(evaluated(i).x_bin_not_nan));
end    


% yyy=logspace(log10(min(y(:,1))),log10(max(y(:,1))),t);
yyy=y(:,1).';

h(3) =  figure;
subplot(2,2,1);
plot(yyy,...
    d11,'k','LineWidth',2,'Marker','o','LineStyle','none')
hold on
plot(yyy,...
    history.x1_iter(index,1:t),'LineWidth',2,'Marker','+','LineStyle','none')
plot(yyy,...
    lb(1:t),'LineWidth',2,'LineStyle','--')
plot(yyy,...
    ub(1:t),'LineWidth',2,'LineStyle','--')
axis square
% ylabel('$d_{11}(r)$', 'interpreter','latex')
ylabel('$d_{11}$', 'interpreter','latex')
    if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
    else
    xlabel('$r/m $', 'interpreter','latex')
    end
xlim([min(y(:,1)) max(y(:,1))])
set(gca, 'FontSize',24)
% set(gca,'xScale','log') 
%
subplot(2,2,2);
plot(yyy,...
    d20,'k','LineWidth',2,'Marker','o','LineStyle','none')
hold on
plot(yyy,...
    history.x1_iter(index,t*1+1:(t*1+1)+t-1),'LineWidth',2,'Marker','+','LineStyle','none')
plot(yyy,...
    lb(t*1+1:(t*1+1)+t-1),'LineWidth',2,'LineStyle','--')
plot(yyy,...
    ub(t*1+1:(t*1+1)+t-1),'LineWidth',2,'LineStyle','--')
axis square
% ylabel('$d_{d20}(r)$', 'interpreter','latex')
ylabel('$d_{d20}$', 'interpreter','latex')
    if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
    else
    xlabel('$r/m $', 'interpreter','latex')
    end
xlim([min(y(:,1)) max(y(:,1))])
set(gca, 'FontSize',24)
% set(gca,'xScale','log') 
%
subplot(2,2,3);
plot(yyy,...
    d21,'k','LineWidth',2,'Marker','o','LineStyle','none')
hold on
plot(yyy,...
    history.x1_iter(index,t*2+1:(t*2+1)+t-1),'LineWidth',2,'Marker','+','LineStyle','none')
plot(yyy,...
    lb(t*2+1:(t*2+1)+t-1),'LineWidth',2,'LineStyle','--')
plot(yyy,...
    ub(t*2+1:(t*2+1)+t-1),'LineWidth',2,'LineStyle','--')
axis square
% ylabel('$d_{d21}(r)$', 'interpreter','latex')
ylabel('$d_{d21}$', 'interpreter','latex')
    if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
    else
    xlabel('$r/m $', 'interpreter','latex')
    end
xlim([min(y(:,1)) max(y(:,1))])
set(gca, 'FontSize',24)
% set(gca,'xScale','log') 
%
subplot(2,2,4);
plot(yyy,...
    d22,'k','LineWidth',2,'Marker','o','LineStyle','none')
hold on
plot(yyy,...
    history.x1_iter(index,t*3+1:(t*3+1)+t-1),'LineWidth',2,'Marker','+','LineStyle','none')
plot(yyy,...
    lb(t*3+1:(t*3+1)+t-1),'LineWidth',2,'LineStyle','--')
plot(yyy,...
    ub(t*3+1:(t*3+1)+t-1),'LineWidth',2,'LineStyle','--')
axis square
% ylabel('$d_{d22}(r)$', 'interpreter','latex')
ylabel('$d_{d22}$', 'interpreter','latex')
    if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
    else
    xlabel('$r/m $', 'interpreter','latex')
    end
xlim([min(y(:,1)) max(y(:,1))])
set(gca, 'FontSize',24)
set(gcf, 'Color', 'w')
% set(gca,'xScale','log') 



h(4) = figure;
subplot(2,2,1);
plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
    co_KM_opti.a(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_KM_opti.ea(1)+co_KM_opti.a(2),'k','LineWidth',2)
hold on
plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
    co_IFT_opti.a(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_IFT_opti.ea(1)+co_IFT_opti.a(2),'LineWidth',2)
plot(yyy,...
    lb(1:t),'LineWidth',2,'LineStyle','--')
plot(yyy,...
    ub(1:t),'LineWidth',2,'LineStyle','--')
axis square
ylabel('$d_{11}$', 'interpreter','latex')
    if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
    else
    xlabel('$r/m $', 'interpreter','latex')
    end
xlim([min(y(:,1)) max(y(:,1))])
set(gca, 'FontSize',24)
% set(gca,'xScale','log') 
%
subplot(2,2,2);
plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
    co_KM_opti.b(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_KM_opti.eb(1)+co_KM_opti.b(2),'k','LineWidth',2)
hold on
plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
    co_IFT_opti.b(1).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_IFT_opti.eb(1)+co_IFT_opti.b(2),'LineWidth',2)
plot(yyy,...
    lb(t*1+1:(t*1+1)+t-1),'LineWidth',2,'LineStyle','--')
plot(yyy,...
    ub(t*1+1:(t*1+1)+t-1),'LineWidth',2,'LineStyle','--')
axis square
ylabel('$d_{d20}$', 'interpreter','latex')
    if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
    else
    xlabel('$r/m $', 'interpreter','latex')
    end
xlim([min(y(:,1)) max(y(:,1))])
set(gca, 'FontSize',24)
% set(gca,'xScale','log') 
%
subplot(2,2,3);
plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
    co_KM_opti.b(3).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_KM_opti.eb(2)+co_KM_opti.b(4),'k','LineWidth',2)
hold on
plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
    co_IFT_opti.b(3).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_IFT_opti.eb(2)+co_IFT_opti.b(4),'LineWidth',2)
plot(yyy,...
    lb(t*2+1:(t*2+1)+t-1),'LineWidth',2,'LineStyle','--')
plot(yyy,...
    ub(t*2+1:(t*2+1)+t-1),'LineWidth',2,'LineStyle','--')
axis square
ylabel('$d_{d21}$', 'interpreter','latex')
    if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
    else
    xlabel('$r/m $', 'interpreter','latex')
    end
xlim([min(y(:,1)) max(y(:,1))])
set(gca, 'FontSize',24)
% set(gca,'xScale','log') 
%
subplot(2,2,4);
plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
    co_KM_opti.b(5).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_KM_opti.eb(3)+co_KM_opti.b(6),'k','LineWidth',2)
hold on
plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
    co_IFT_opti.b(5).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^co_IFT_opti.eb(3)+co_IFT_opti.b(6),'LineWidth',2)
plot(yyy,...
    lb(t*3+1:(t*3+1)+t-1),'LineWidth',2,'LineStyle','--')
plot(yyy,...
    ub(t*3+1:(t*3+1)+t-1),'LineWidth',2,'LineStyle','--')
axis square
ylabel('$d_{d22}$', 'interpreter','latex')
    if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
    else
    xlabel('$r/m $', 'interpreter','latex')
    end
xlim([min(y(:,1)) max(y(:,1))])
set(gca, 'FontSize',24)
set(gcf, 'Color', 'w')
% set(gca,'xScale','log') 



    

[xData, yData, zData, weights] = prepareSurfaceData( xx, y, z1_opti, z1e );
tmp_u=repmat(linspace(min(xData),max(xData),100),100,1);
tmp_r=repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100);

h(5) = figure;
scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
hold on
surf(tmp_u,tmp_r,(co_IFT_opti.a(1).*tmp_r.^co_IFT_opti.ea(1)+co_IFT_opti.a(2)).*tmp_u,'EdgeColor','none')
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca, 'FontSize',18)

if norm_ur==1
xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
else
xlabel('$u_r$', 'interpreter','latex')   
end
if norm_r==1
ylabel('$r / \lambda$', 'interpreter','latex')
else
ylabel('$r/m $', 'interpreter','latex')
end
zlabel('$D1 (u_r,r)$', 'interpreter','latex')
% zlabel('$D1 (x,r) / [\frac{m}{s^2}]$', 'interpreter','latex')
colormap (parula(200))
% colorbar
grid off
xlim([min(xData) max(xData)])
ylim([min(y(:,1)) max(y(:,1))])
% zlim([-0.005 0.005])
axis square
set(gca,'YDir','reverse');


[xData, yData, zData, weights] = prepareSurfaceData( xx, y, z2_opti, z2e );
tmp_u=repmat(linspace(min(xData),max(xData),100),100,1);
tmp_r=repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100);
h(6) = figure;
scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
hold on
surf(tmp_u,tmp_r,...
    ((co_IFT_opti.b(1).*tmp_r.^co_IFT_opti.eb(1)+co_IFT_opti.b(2)).*tmp_u.^0)+...
    ((co_IFT_opti.b(3).*tmp_r.^co_IFT_opti.eb(2)+co_IFT_opti.b(4)).*tmp_u.^1)+...
    ((co_IFT_opti.b(5).*tmp_r.^co_IFT_opti.eb(3)+co_IFT_opti.b(6)).*tmp_u.^2),'EdgeColor','none')
set(gca,'yScale','log') 
% set(gca,'zScale','log')  
set(gcf, 'Color', 'w')
set(gca, 'FontSize',18)
if norm_ur==1
xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
else
xlabel('$u_r$', 'interpreter','latex')   
end
if norm_r==1
ylabel('$r / \lambda$', 'interpreter','latex')
else
ylabel('$r/m $', 'interpreter','latex')
end
% zlabel('$D2 (x,r) / [\frac{m^2}{s^3}]$', 'interpreter','latex')
zlabel('$D2 (u_r,r)$', 'interpreter','latex')
colormap (parula(200))
% colorbar
grid off
xlim([min(xData) max(xData)])
ylim([min(y(:,1)) max(y(:,1))])
% zlim([0 2*max(zData)])
% caxis([0 max(zData)])
axis square
% close
set(gca,'YDir','reverse');

if ischar(save_path)
    savefig(h,fullfile(save_path,append(save_name,'_','opti_IFT.fig')),'compact')
	    for a = 1:length(h)
        exportgraphics(h(a),fullfile(save_path,append(save_name,'_',sprintf('opti_IFT_%d.png', a))))
    end
end
for i=1:6
    tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
    uiwait(gcf);
    delete(tmp_ui);
    close
end
end