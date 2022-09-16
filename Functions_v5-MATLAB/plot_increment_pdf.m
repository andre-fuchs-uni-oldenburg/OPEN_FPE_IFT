%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots the probability density function (PDF) of the velocity increments 
% at the scale $r=L$, $r=\lambda$ and  $r=\eta$. The colored dashed line correspond to Castaing fits 
% (form factor $\Lambda_r$ \cite{Castaing_1990}) and grey dashed line to Gaussian fits.
%
% Arguments IN
% data = 1D array of data to which you would like to apply the low pass filter
% m_data = mean of the data
% nbins = The number of bins in which you would like to divide your data
% Fs = Acquisition/Sampling Frequency in Hz
% int_L = Integral length scale in meters
% taylor_L = Taylors length scale in meters
% diss_scale = Kolmogorv lenght scale in meters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_increment_pdf(data,nbins,save_path,save_name,int_L,taylor_L,m_data,Fs,norm_ur,varargin)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaulttextinterpreter','latex');   
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTitle','latex');

tau1            = round((int_L/m_data)*Fs);
tau2            = round((taylor_L/m_data)*Fs);

[incr1,incr2]   = Increment(tau1,tau2,data);

incr1           = incr1.';
incr2           = incr2.';

if nargin >= 9 
    tau3                = round((varargin{1}/m_data)*Fs);
    [incr1,incr2,incr3] = Increment_double_condition(tau1,tau2,tau3,data);
    incr1               = incr1.';
    incr2               = incr2.';
    incr3               = incr3.';
end

%     incr1               = incr1.'/std(incr1);
%     incr2               = incr2.'/std(incr2);
%     incr3               = incr3.'/std(incr3);

% x = linspace(-6,6,nbins);
% [fuL, uL]       = hist(incr1,x); 

[fuL, uL]       = hist(incr1,nbins); 
puL             = fuL / ( nansum(fuL) * mean(diff(uL)) );

[ful, ul]       = hist(incr2,uL); 
pul             = ful / ( nansum(ful) * mean(diff(ul)) );

h(1) = figure;
plot(uL,puL,'LineWidth',2)
hold on
plot(ul,pul,'LineWidth',2)
ylim([0.000001*max(puL) 2*max(pul)])

if nargin >= 8   
    [fue, ue]       = hist(incr3,uL); 
    pue             = fue / ( nansum(fue) * mean(diff(ue)) );
    plot(ue,pue,'k','LineWidth',2)

%%    castaing
    Lam_sq_L   = log(kurtosis(incr1)/3)/4;
    sig_sq          = mean(incr1.^2)*exp(-2*Lam_sq_L);
    dX                = 0.01;
    lim                = 20;
%     lim                = range(data);
    intbound       = lim/dX;
    b                   = exp(Lam_sq_L);
   %asymmetry   = skewness(incr1);
    asymmetry    = 0;
    pdf_Cast       = nan(1,length(uL));
    X                   = ((0:intbound) + 1) * dX; 
    
    for k=1:length(uL)
        P           = trapz(X,((1./(X.^2)).* ...
        exp(-(0.5 .* (uL(k)./X).^2.* (1 + asymmetry .* uL(k)./(X .* sqrt(1 + (uL(k)./X).^2))))) .* ...
        exp((-1 .*log(X./sqrt(sig_sq)).^2)./(2 .* Lam_sq_L))...
        ));

        pdf_Cast(k) = P/(2 * pi* sqrt(Lam_sq_L));     
    end
    pdf_Cast_L      = pdf_Cast/trapz(uL,pdf_Cast);
    plot(uL,pdf_Cast_L,'LineWidth',2,'LineStyle','--','Color',[0    0.4470    0.7410])

    Lam_sq_l         = log(kurtosis(incr2)/3)/4;
    sig_sq              = mean(incr2.^2)*exp(-2*Lam_sq_l);
    b                       = exp(Lam_sq_l);
    
    for k=1:length(ul)    
        P           = trapz(X,((1./(X.^2)).* ...
        exp(-(0.5 .* (ul(k)./X).^2.* (1 + asymmetry .* ul(k)./(X .* sqrt(1 + (ul(k)./X).^2))))) .* ...
        exp((-1 .*log(X./sqrt(sig_sq)).^2)./(2 .* Lam_sq_l))...
        ));

        pdf_Cast(k) = P/(2 * pi* sqrt(Lam_sq_l));     
    end
    pdf_Cast_l      = pdf_Cast/trapz(ul,pdf_Cast);
    plot(ul,pdf_Cast_l,'LineWidth',2,'LineStyle','--','Color',[0.8500    0.3250    0.0980])

    Lam_sq_e        = log(kurtosis(incr3)/3)/4;
    sig_sq          = mean(incr3.^2)*exp(-2*Lam_sq_e);
    b               = exp(Lam_sq_e);
    
    for k=1:length(ue)    
        P           = trapz(X,((1./(X.^2)).* ...
        exp(-(0.5 .* (ue(k)./X).^2.* (1 + asymmetry .* ue(k)./(X .* sqrt(1 + (ue(k)./X).^2))))) .* ...
        exp((-1 .*log(X./sqrt(sig_sq)).^2)./(2 .* Lam_sq_e))...
        ));

        pdf_Cast(k) = P/(2 * pi* sqrt(Lam_sq_e));     
    end
    pdf_Cast_e      = pdf_Cast/trapz(ue,pdf_Cast);
    plot(ue,pdf_Cast_e,'LineWidth',2,'LineStyle','--','Color','k')


    norm1 = pdf(fitdist(incr1,'Normal'),uL);
    plot(uL,norm1,'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176])
    norm1 = pdf(fitdist(incr2,'Normal'),uL);
    plot(uL,norm1,'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176])
    norm1 = pdf(fitdist(incr3,'Normal'),uL);
    plot(uL,norm1,'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176])
    legend({'$u_L$','$u_\lambda$','$u_\eta$'},'Interpreter','latex','Location','northeast','FontSize',18)  
    legend('$u_L$','$u_\lambda$','$u_\eta$',['Cast. $\Lambda_r=$',num2str(Lam_sq_L,'%1.2f')],['Cast. $\Lambda_r=$',num2str(Lam_sq_l,'%1.2f')],['Cast. $\Lambda_r=$',num2str(Lam_sq_e,'%1.2f')],'Gaussian fit','FontSize',18);
    ylim([0.000001*max(puL) 2*max(pue)])
    title(['$S(u_L)$ = ',num2str(skewness(incr1),'%1.2f'), ';',...
          ' $S(u_\lambda)$ = ',num2str(skewness(incr2),'%1.2f'), ';',...
          ' $S(u_\eta)$ = ',num2str(skewness(incr3),'%1.2f'), ';',...
    ],'interpreter','latex')
else
    %%    castaing
    Lam_sq_L    = log(kurtosis(incr1)/3)/4;
    sig_sq      = mean(incr1.^2)*exp(-2*Lam_sq_L);
    dX          = 0.01;
    lim         = 20;
    intbound    = lim/dX;
    b           = exp(Lam_sq_L);
   %asymmetry   = skewness(incr1);
    asymmetry   = 0;
    pdf_Cast    = nan(1,length(uL));
    X = ((0:intbound) + 1) * dX;
    
    for k=1:length(uL)
        P           = trapz(X,((1./(X.^2)).* ...
        exp(-(0.5 .* (uL(k)./X).^2.* (1 + asymmetry .* uL(k)./(X .* sqrt(1 + (uL(k)./X).^2))))) .* ...
        exp((-1 .*log(X./sqrt(sig_sq)).^2)./(2 .* Lam_sq_L))...
        ));

        pdf_Cast(k) = P/(2 * pi* sqrt(Lam_sq_L));     
    end
    pdf_Cast_L      = pdf_Cast/trapz(uL,pdf_Cast);
    plot(uL,pdf_Cast_L,'LineWidth',2,'LineStyle','--','Color',[0    0.4470    0.7410])

    Lam_sq_l        = log(kurtosis(incr2)/3)/4;
    sig_sq          = mean(incr2.^2)*exp(-2*Lam_sq_l);
    b               = exp(Lam_sq_l);
    
    for k=1:length(ul)    
        P           = trapz(X,((1./(X.^2)).* ...
        exp(-(0.5 .* (ul(k)./X).^2.* (1 + asymmetry .* ul(k)./(X .* sqrt(1 + (ul(k)./X).^2))))) .* ...
        exp((-1 .*log(X./sqrt(sig_sq)).^2)./(2 .* Lam_sq_l))...
        ));

        pdf_Cast(k) = P/(2 * pi* sqrt(Lam_sq_l));     
    end
    pdf_Cast_l      = pdf_Cast/trapz(ul,pdf_Cast);
    plot(ul,pdf_Cast_l,'LineWidth',2,'LineStyle','--','Color',[0.8500    0.3250    0.0980])

    norm1 = pdf(fitdist(incr1,'Normal'),uL);
    plot(uL,norm1,'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176])
    norm1 = pdf(fitdist(incr2,'Normal'),uL);
    plot(uL,norm1,'LineWidth',2,'LineStyle','--','Color',[0.501960813999176 0.501960813999176 0.501960813999176])
%     legend({'$u_L$','$u_\lambda$'},'Interpreter','latex','Location','northeast','FontSize',18)
    legend('$u_L$','$u_\lambda$',['Cast. $\Lambda_r=$',num2str(Lam_sq_L,'%1.2f')],['Cast. $\Lambda_r=$',num2str(Lam_sq_l,'%1.2f')]);

    title(['$S(u_L)$ = ',num2str(skewness(incr1),'%1.2f'), ';',...
          ' $S(u_\lambda)$ = ',num2str(skewness(incr2),'%1.2f'), ';',...
    ],'interpreter','latex')
end
set(gca,'yScale','log') 
set(gcf, 'Color', 'w')
set(gca,'FontSize',22)
% ylabel('$p(u_r / \sigma_\infty)$','interpreter','latex')

ylabel('PDF','interpreter','latex')

if norm_ur==1  
    xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
else
    xlabel('$u_r\ (m/s)$','interpreter','latex');
end



axis square
xlim([min(uL) max(uL)])
fig_setup
%     set(gca,'XTick',[10^-3 10^-2  10^-1  10^0]);

tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf);
delete(tmp_ui);

if ischar(save_path)
savefig(h,fullfile(save_path,append(save_name,'_','increment_pdf.fig')),'compact')
    for a = 1:length(h)
        exportgraphics(h(a),fullfile(save_path,append(save_name,'_',sprintf('increment_pdf_%d.png', a))))
    end
end
close
end