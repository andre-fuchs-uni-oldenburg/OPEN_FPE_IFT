function [Lint,lambda,lambda1,EpsNew,h]=ZeroCrossings(data,Fs,m_data,nu,low_freq,x_achse_k_spektrum_2pi,E_k,h)
if low_freq== Fs/2
    fc  = low_freq-1;
else
    fc  = low_freq;
end

% f0  = 0.02;
f0  = 50/(length(data)/Fs);
% I did it in a way where we only get lambda and L
FilterSizes         = logspace(log10(f0),log10(fc),60);

etac                = m_data./FilterSizes;

[ns,l,sigmavoro]    = computeZeros(data,Fs,FilterSizes,m_data);

C       = 1;
lambda  = l/C/pi;
II      = ~isnan(sigmavoro);
Lint    = interp1(sigmavoro(II)./sqrt(1/2),etac(II),1,'linear','extrap');


% 
% h(13) = figure;
% loglog(1./etac,ns.*l,'MarkerEdgeColor','k','MarkerSize',6,'Marker','o', 'LineStyle','none','LineWidth',1.5);
% % loglog(Lint./etac,ns.*l,'MarkerEdgeColor','k','MarkerSize',6,'Marker','o', 'LineStyle','none','LineWidth',1.5);
% ylabel('$n_s l$','FontName','times');
% xlabel('$\eta_c^{-1}~(m^{-1})$','FontName','times');
% axis square
% set(gcf, 'Color', 'w')
% set(gca, 'FontSize',18)
% fig_setup
% legend off
% % xlim([0.8*min(1./etac) 2*max(1./etac)])
% % ylim([0.05 1.5])
% 
% 
h(12) = figure;
loglog(1./etac(II),sigmavoro(II)./sqrt(1/2),'MarkerEdgeColor','k','MarkerSize',6,'Marker','o', 'LineStyle','none','LineWidth',1.5);
% loglog(Lint./etac(II),sigmavoro(II)./sqrt(1/2),'MarkerEdgeColor','k','MarkerSize',6,'Marker','o', 'LineStyle','none','LineWidth',1.5);
hold on;
plot(1./etac(II),ones(numel(etac(II)),1),'k','LineStyle','--','LineWidth',2);
% plot(Lint./etac(II),ones(numel(etac(II)),1),'k');
% xlim([0.8*min(1./etac) 2*max(1./etac)])
% ylim([0.7 2])
tmp=1./etac(II);
vline(tmp((abs(sigmavoro(II)./sqrt(1/2)-1))==min(abs(sigmavoro(II)./sqrt(1/2)-1))),'k')
% plot(ones(numel(etac(II)),1),sigmavoro(II)./sqrt(1/2),'k')
ylabel('$\sigma_{voro}/\sigma_{RPP}$','FontName','times');
xlabel('$\eta_c^{-1}~(m^{-1})$','FontName','times');
xlabel('$\eta_c^{-1}\ (1/m)$', 'interpreter','latex');
axis square
set(gcf, 'Color', 'w')
set(gca, 'FontSize',18)
fig_setup
legend off
hold off
set(gca,'XTick',[10^-1 10^0  10^1  10^2  10^3 10^4]);

tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf);
delete(tmp_ui);

% figure
% plot(ns.*l)

%Martin: The figure shows the plot of ns*l vs 1/etac
h(13) = figure;
loglog(1./etac,ns.*l,'MarkerEdgeColor','k','MarkerSize',6,'Marker','o', 'LineStyle','none','LineWidth',1.5);
% loglog(Lint./etac,ns.*l,'MarkerEdgeColor','k','MarkerSize',6,'Marker','o', 'LineStyle','none','LineWidth',1.5);
ylabel('$n_s l$','FontName','times');
xlabel('$\eta_c^{-1}~(m^{-1})$','FontName','times');
xlabel('$\eta_c^{-1}\ (1/m)$', 'interpreter','latex');
axis square
set(gcf, 'Color', 'w')
set(gca, 'FontSize',18)
fig_setup
legend off
% xlim([0.8*min(1./etac) 2*max(1./etac)])
% ylim([0.05 1.5])
set(gca,'XTick',[10^-1 10^0  10^1  10^2  10^3 10^4]);
tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf);
delete(tmp_ui);








% If the velocity derivatives are resolved we can get a better approximation of lambda;
[b,a]   = butter(6,FilterSizes(end)/(Fs/2));
UU      = filter(b,a,data);
UUp1    = gradient(UU,1/Fs)./m_data;
% flat1   = kurtosis(data);

mdu     = mean(abs(UUp1));
sdu     = std(UUp1);
C1      = sdu./mdu.*sqrt(2/pi);
lambda1 = l/C1/pi;

[EpsNew,h]        = PostProcEps(data,m_data,Fs,nu,x_achse_k_spektrum_2pi,E_k,h);
end


function [ns,l,sigmavoro]=computeZeros(data,Fs,FilterSizes,m_data)
for j=1:numel(FilterSizes)
    [b,a]       = butter(6,FilterSizes(j)/(Fs/2));
    UU          = filter(b,a,data);
    sU          = UU-mean(UU);
    S           = sign(sU);
    I           = find((abs(diff(S)))>1);
    ns(j)       = numel(I)./numel(data).*Fs./m_data;
    voro        = (diff(I(2:end))+diff(I(1:end-1)))./2./Fs.*m_data;
    sigmavoro(j)= (std(voro)./mean(voro));
end
l               = mean(diff(I)./Fs.*m_data);
end


function [EpsNew,h]=PostProcEps(data,m_data,Fs,nu,x_achse_k_spektrum_2pi,E_k,h)
% UU              = m_data;
% UUstd           = std(data);
% U_tmp1          = m_data; 
% Fs              = Fs.*2.*pi./m_data;

% Nfft            = 4096*4;
% [P,f]           = pwelch(data-m_data,hanning(Nfft),Nfft/4,Nfft,(Fs.*2.*pi./m_data));
    
% a               = 0;
% while a==0
%     FIG=figure;
%     loglog(P.*(f).^2.*15.*nu,'LineWidth',2);
%     set(FIG, 'Position', [1, 1, 1000, 800]);
%     xlim([1 numel(P)]

%     b=1;
%     while b==1
        %Andre
        y_tmp = E_k.*(x_achse_k_spektrum_2pi).^2.*15.*nu;
        h(14) = figure;
        loglog(x_achse_k_spektrum_2pi,y_tmp,'-','LineWidth',2)
        hold on
        xlabel('$k= / \frac{1}{m}$', 'interpreter','latex');
        ylabel('E(k)$\cdot k^2 \cdot 15\nu $ / $\frac{m^3}{s^3}$', 'interpreter','latex')
        xlabel('$k\ (1/m)$','interpreter','latex');
%         xlabel('$k=\frac{2 \pi f}{<u>} / m^{-1}$', 'interpreter','latex');
        ylabel('$E k^2 15\nu\ (m^3/s^3)$', 'interpreter','latex')
        
        set(gca,'YScale','log');
        axis square
        set(gcf, 'Color', 'w')
        set(gca, 'FontSize',18)
        fig_setup
        legend off
        set(gca,'XTick',[10^-2 10^0  10^2  10^4]);
        xlim([0.0001*max(x_achse_k_spektrum_2pi) 1.5*max(x_achse_k_spektrum_2pi) ])
        ylim([0.00001*max(y_tmp) 2*max(y_tmp)])
%         legend off
%         txt = {'(b)'};
%         b=text(4,0.5,txt);
%         b.FontSize = 22;
%         b.Units='normalized';
%         pos_txt=[-0.19   0.9 0];
%         a.Position=pos_txt;
%         b.Position=pos_txt;
%         set(gca,'XTick',[10^-1 10^0  10^1  10^2  10^3 10^4]);
        
        
        waitfor(msgbox('Use interactive data cursor to select the fit range of wave numbers (see readme) is used to model the the dissipation spectrum.'));
        
        [I,yy]  = ginput(2);
        if I(1)>I(2)
            I   = flip(I);
            yy  = flip(yy);
        end     
%         I1      = floor(I(1));
%         I2      = floor(I(2));

%         I(2,1)  = (2.*pi.*fc)./m_data;
%         yy(2,1) = y_tmp((abs(x_achse_k_spektrum_2pi-I(2,1)))==min(abs(x_achse_k_spektrum_2pi-I(2,1))))
        ind_1   = find(abs(x_achse_k_spektrum_2pi-I(1,1))==min(abs(x_achse_k_spektrum_2pi-I(1,1))));
        ind_2   = find(abs(x_achse_k_spektrum_2pi-I(2,1))==min(abs(x_achse_k_spektrum_2pi-I(2,1))));
        
        I1      = x_achse_k_spektrum_2pi(ind_1);
        I2      = x_achse_k_spektrum_2pi(ind_2);
        
        yy(1,1) = y_tmp(ind_1);
        yy(2,1) = y_tmp(ind_2);
        
        plot(I,yy,'MarkerEdgeColor','k','MarkerSize',6,'Marker','o', 'LineStyle','none','LineWidth',1.5);
%         xlim([1 max(f.*2.*pi/m_data)])

        % Andre
        [xData, yData] = prepareCurveData( x_achse_k_spektrum_2pi(ind_1:ind_2), y_tmp(ind_1:ind_2));
        % Set up fittype and options.
        ft      = fittype( 'power1' );
        opts    = fitoptions( 'Method', 'NonlinearLeastSquares' );
        [fitresult, gof] = fit( xData, yData, ft, opts );  
        plot(x_achse_k_spektrum_2pi(ind_1:end), feval(fitresult,x_achse_k_spektrum_2pi(ind_1:end)),'--','LineWidth',2,'Color',[0 0 0])
%         fitresult
        b       = askYesno('Would you like to repeat the selection?', 'Yes');
        if b==1
            close
        end 
%     end


%         tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
%         uiwait(gcf);
%         delete(tmp_ui);
            
    
    if I>max(x_achse_k_spektrum_2pi)
        EpsNew  = trapz(x_achse_k_spektrum_2pi,y_tmp);
%         EpsNew  = trapz(f,P.*(f).^2.*15.*nu);
        a=1;
    else
%         pow1    = @(p,x)p(1).*abs(x).^p(2);
%         x1      = [1e6 -4];        % If the fit is no good, you can change here the first element
%         [p,resnorm,~,exitflag,output] = lsqcurvefit(pow1,x1,x_achse_k_spektrum_2pi(ind_1:ind_2),E_k(ind_1:ind_2).*(x_achse_k_spektrum_2pi(ind_1:ind_2)).^2.*15.*nu);
        
%         [p, sig, ff, covp,corp,r2,rv]=mf_flsqr(f(I1:I2),P(I1:I2).*(f(I1:I2)).^2.*15.*nu,[],[1e6 1e-20 -4 1e-20],[1 0 1 0],@pow1,[0.01 5000 1e-6]);        
        
%         if p(2)>=-1
%             EpsNew  = trapz(x_achse_k_spektrum_2pi(1:ind_2),y_tmp(1:ind_2));
%         else   
%             [xData, yData] = prepareCurveData( x_achse_k_spektrum_2pi(ind_1:ind_2), y_tmp(ind_1:ind_2));
%             EpsNew = trapz(x_achse_k_spektrum_2pi(1:ind_1),y_tmp(1:ind_1))+trapz(xData,feval(fitresult,xData));
			EpsNew = trapz(x_achse_k_spektrum_2pi(1:ind_1),y_tmp(1:ind_1))+trapz(x_achse_k_spektrum_2pi(ind_1:end),feval(fitresult,x_achse_k_spektrum_2pi(ind_1:end)));
            % not more than 5% should be modeld  
%             p(1).*ind_1.^(p(3)+1)/(p(3)+1);
%         end
        
%         fcut = floor(x_achse_k_spektrum_2pi(ind_2));
%         figure;
%         loglog(x_achse_k_spektrum_2pi,E_k.*(f).^2.*15.*nu,'LineWidth',2);
%         hold on;
%         loglog(x_achse_k_spektrum_2pi(ind_1:end),pow1(x_achse_k_spektrum_2pi(ind_1:end),p),'r','LineWidth',2)
%         a   = input('If Ok press 1, If wrong press 0');
%         tmp =(-p(1).*ind_2.^(p(3)+1)/(p(3)+1)./EpsNew);
%         display(tmp)
    end
    

    title('dissipation spectrum','interpreter','latex')
tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
uiwait(gcf);
delete(tmp_ui);
% close all
end