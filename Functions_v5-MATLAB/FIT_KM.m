%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs the surface fit with a linear function for $D^{(1)}\left(u_r,r\right)$ 
% and a parabolic function for $D^{(2)}\left(u_r,r\right)$ to the optimized and non-optimized KMC's.
% Coefficients $d_{ij}(r)$ in the fits are functions of scale~$r$ of the form 
% $\alpha (r/\lambda)^{\beta}+\gamma$. After fitting, this function plots the optimized $D^{(1,2)}$ 
% and its surface fits. This function also plots the parameters $d_{11}$, $d_{20}$, $d_{21}$ and 
% $d_{22}$ as a function of $\frac{r}{\lambda}$ for optimized and non-optimized 
% $D^{(1,2)}\left(u_r,r\right)$.
%
% Arguments IN
% evaluated_point = struct array calculated in the function 'conditional_moment'
% increment_bin = number of bins
% taylor_L = Taylor length scale in meters
% Fs = Acquisition/Sampling Frequency in Hz
% m_data = mean of the data
% markov = markov length in number of samples
% multi_point = weather to do multi-point analysis or not 1=Yes, 0=No
% condition = condition for multi-point analysis
% norm_ur = normalization of the data using $\sigma_\infty$? data 1=Yes, 0=No
% norm_r = normalization of the scale using $\lambda$? data 1=Yes, 0=No
% save_path = path for saving figures and files
% save_name = name for saving files
%
% Arguments OUT
% co_KM_opti = Coefficients $d_{ij}(r)$ of the optimized Kramers-Moyal coefficients using the surface fits
% co_KM_opti_old = Coefficients $d_{ij}(r)$ of the optimized Kramers-Moyal coefficients using the
% surface fits without an offset
% co_no_opti = Coefficients $d_{ij}(r)$ of the non-optimized Kramers-Moyal coefficients using the surface fits
% fitresult_D1_conf = Confidence intervals for fit coefficients of D1
% fitresult_D2_conf = Confidence intervals for fit coefficients of D2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[co_KM_opti,co_KM_opti_old,co_no_opti,fitresult_D1_conf,fitresult_D2_conf]=FIT_KM(evaluated_point,increment_bin,taylor_L,Fs,m_data,markov,multi_point,condition,norm_ur,norm_r,save_path,save_name)

% taylor_L=1;

if multi_point==1
    for k=1:length(condition)
            my_field = sprintf('point_eval_%d', k);              
            co_KM_opti.(my_field)               = struct('a',NaN(1,2),'ea',nan,'b',NaN(1,6),'eb',NaN(1,3));
            co_KM_opti_old.(my_field)           = struct('a',NaN(1,2),'ea',nan,'b',NaN(1,6),'eb',NaN(1,3));
            co_no_opti.(my_field)               = struct('a',NaN(1,2),'ea',nan,'b',NaN(1,6),'eb',NaN(1,3));
    end
    
    for k=1:length(condition)
            my_field    = sprintf('point_eval_%d', k);    
            clear evaluated
            evaluated   = evaluated_point.(my_field);
            
            tmp=0;
            for scal    = 1:length(evaluated)
                tmp     = tmp+sum(evaluated(scal).x_bin_not_nan>0);
            end
            
            if tmp>0
                xx      = nan(length(evaluated),increment_bin);
                y       = nan(length(evaluated),increment_bin);
                z1      = nan(length(evaluated),increment_bin);
                z1e     = nan(length(evaluated),increment_bin);
                z2      = nan(length(evaluated),increment_bin);
                z2e     = nan(length(evaluated),increment_bin);
                z1_opti = nan(length(evaluated),increment_bin);
                z2_opti = nan(length(evaluated),increment_bin);
                z4      = nan(length(evaluated),increment_bin);

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

                %TEST Beginn
                % y=y-(((markov./Fs).*m_data)/taylor_L);
                %TEst Ende


                %% D1-Fit
                [fitresult_D1, xData, yData, zData, weights]    = FIT_D1(xx, y, z1_opti, z1e);
                fitresult_D1_conf                               = confint(fitresult_D1);

                h(1) = figure;
                % surf(xx, y, z1_opti,'EdgeColor','none')
                scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
                % scatter3(xData,yData,yData.*zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
                hold on
                surf(repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100),feval(fitresult_D1,repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100)),'EdgeColor','none')
%                 surf(repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100),...
%                     (-0.36.*repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100).^(-1).*repmat(linspace(min(xData),max(xData),100),100,1)),'FaceColor',[0.501960813999176 0.501960813999176 0.501960813999176],'EdgeColor',[0.501960813999176 0.501960813999176 0.501960813999176])

                % surf(repmat(linspace(-4,4,100),100,1),repmat(linspace(min(y(:,1)),max(y(:,1)),100).',1,100),feval(fitresult_D1,repmat(linspace(-4,4,100),100,1),repmat(linspace(min(y(:,1)),max(y(:,1)),100).',1,100)))
                set(gca,'yScale','log') 
                set(gcf, 'Color', 'w')
                set(gca, 'FontSize',18)
                if norm_ur==1
                    xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
                else
                    xlabel('$u_r\ (m/s)$', 'interpreter','latex')   
                end
                if norm_r==1
                    ylabel('$r / \lambda$', 'interpreter','latex')
                else
                    ylabel('$r\ (m)$', 'interpreter','latex')
                end
                zlabel('$D^{(1)} (u_r,r)$', 'interpreter','latex')
                % zlabel('$D1 (x,r) / [\frac{m}{s^2}]$', 'interpreter','latex')
                colormap (parula(200))
                % colorbar
                grid off
                xlim([min(xData) max(xData)])
                ylim([min(y(:,1)) max(y(:,1))])
                % zlim([-0.005 0.005])
                axis square
                set(gca,'YDir','reverse');
                title('optimized KM, optimized fit','interpreter','latex');
                % title(['$(', num2str(fitresult_D1.a1,'%1.1E'), ' * r / \lambda^{' num2str(fitresult_D1.ra1,'%1.1E') '} ) * u_r $'],'interpreter','latex')   
                % title(['$\left(',num2str(fitresult_D1.a1),'\cdot r/\lambda^{',num2str(fitresult_D1.ra1),'} +',num2str(fitresult_D1.a11),' \right)\cdot u_r $'],'interpreter','latex')   
                rotate3d on
                colorbar('north','TickLabelInterpreter','latex','AxisLocation','out');
                fig_setup
                legend off
                txt = {'(a)'};
                a=text(4,0.5,txt);
                a.Units='normalized';

                %%%%%  2019 %%%%%%%%
                % optimized coefficients 
                % D1(x,y) = (a1*y^ra1+a11)*x
                co_KM_opti.(my_field).a  = [fitresult_D1.a1 fitresult_D1.a11];
                co_KM_opti.(my_field).ea = [fitresult_D1.ra1];

                % non-optimized coefficients 
                [fitresult_D1_o, xData_o1, yData_o1, zData_o1, weights] = FIT_D1(xx, y, z1, z1e);
                %%%%%  2019 %%%%%%%%
                % D1(x,y) = (a1*y^ra1+a11)*x
                co_no_opti.(my_field).a  = [fitresult_D1_o.a1 fitresult_D1_o.a11];
                co_no_opti.(my_field).ea = [fitresult_D1_o.ra1];

                % optimized coefficients no offset
                % D1(x,y) = (a1*y^ra1)*x
                ft = fittype( '(a1*y^ra1)*x', 'independent', {'x', 'y'}, 'dependent', 'z' );
                opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                % opts.DiffMinChange = 1e-12;
                % opts.Display = 'Off';
                % opts.MaxFunEvals = 1000;
                % opts.MaxIter = 1000;
                % opts.TolFun = 1e-12;
                % opts.TolX = 1e-12;
                opts.Weights = weights;

                % opts.StartPoint = [-0.1  -1];
                % opts.Lower = [-2 -2];
                % opts.Upper = [0 0];

                opts.StartPoint = [0  0];
                opts.Lower = [-1 -Inf];
                opts.Upper = [1 0];

                [fitresult_D1_old, gof] = fit( [xData, yData], zData, ft, opts );
                co_KM_opti_old.(my_field).a  = [fitresult_D1_old.a1 0];
                co_KM_opti_old.(my_field).ea = [fitresult_D1_old.ra1];

                % figure
                % plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),fitresult_D1.a0.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D1.ra0)
                % hold on
                % plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),fitresult_D1.a1.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D1.ra1)
                % plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),fitresult_D1.a2.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^0)
                % plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),fitresult_D1.a3.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D1.ra3)
                % set(gca,'xScale','log') 
                % axis square
                % set(gcf, 'Color', 'w')
                % % grid on
                % xlabel('$r / \lambda$', 'interpreter','latex')
                % legend('d_{10}(r)','d_{11}(r)','d_{12}(r)','d_{13}(r)')
                % xlim([min(y(:,1)) max(y(:,1))])



                %% D2-Fit
                [fitresult_D2, xData, yData, zData, weights] = FIT_D2(xx, y, z2_opti, z2e);
                fitresult_D2_conf                            = confint(fitresult_D2);
                h(2) = figure;
                % surf(xx, y, z2_opti,'EdgeColor','none')
                scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
                % scatter3(xData,yData,yData.*zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
                hold on
                surf(repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100),feval(fitresult_D2,repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100)),'EdgeColor','none')
%                 surf(repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100),...
%                     (0.0144.*repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100).^(-1).*repmat(linspace(min(xData),max(xData),100),100,1).^(2)),'FaceColor',[0.501960813999176 0.501960813999176 0.501960813999176],'EdgeColor',[0.501960813999176 0.501960813999176 0.501960813999176])
                % surf(repmat(linspace(-4,4,100),100,1),repmat(linspace(min(y(:,1)),max(y(:,1)),100).',1,100),feval(fitresult_D1,repmat(linspace(-4,4,100),100,1),repmat(linspace(min(y(:,1)),max(y(:,1)),100).',1,100)))
                set(gca,'yScale','log') 
                % set(gca,'zScale','log')  
                set(gcf, 'Color', 'w')
                set(gca, 'FontSize',18)
                if norm_ur==1
                    xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
                else
                    xlabel('$u_r\ (m/s)$', 'interpreter','latex')   
                end
                if norm_r==1
                    ylabel('$r / \lambda$', 'interpreter','latex')
                else
                    ylabel('$r\ (m)$', 'interpreter','latex')
                end
                % zlabel('$D2 (x,r) / [\frac{m^2}{s^3}]$', 'interpreter','latex')
                zlabel('$D^{(2)} (u_r,r)$', 'interpreter','latex')
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
                set(gca,'TickLabelInterpreter','latex');
                title('optimized KM, optimized fit','interpreter','latex');
                % title(['$(', num2str(fitresult_D1.a1,'%1.1E'), ' * r / \lambda^{' num2str(fitresult_D1.ra1,'%1.1E') '} ) * u_r $'],'interpreter','latex')   
                % title(['$\left(',num2str(fitresult_D1.a1),'\cdot r/\lambda^{',num2str(fitresult_D1.ra1),'} +',num2str(fitresult_D1.a11),' \right)\cdot u_r $'],'interpreter','latex')   
                rotate3d on
                colorbar('north','TickLabelInterpreter','latex','AxisLocation','out');
                fig_setup
                legend off
                txt = {'(b)'};
                b=text(4,0.5,txt);
                b.Units='normalized';

                % optimized coefficients 
                %%%%%  2019 %%%%%%%%
                % D2(x,y) = (b0*y^rb0+b00) + (b1*y^rb1+b11)*x + (b2*y^rb2+b22)*x^2
                co_KM_opti.(my_field).b  = [fitresult_D2.b0 fitresult_D2.b00 fitresult_D2.b1 fitresult_D2.b11 fitresult_D2.b2 fitresult_D2.b22];
                co_KM_opti.(my_field).eb = [fitresult_D2.rb0 fitresult_D2.rb1 fitresult_D2.rb2 ];     % for D2

                % non-optimized coefficients 
                [fitresult_D2_o, xData_o2, yData_o2, zData_o2, weights] = FIT_D2(xx, y, z2, z2e);
                %%%%%  2019 %%%%%%%%
                % D2(x,y) = (b0*y^rb0+b00) + (b1*y^rb1+b11)*x + (b2*y^rb2+b22)*x^2
                co_no_opti.(my_field).b  = [fitresult_D2_o.b0 fitresult_D2_o.b00 fitresult_D2_o.b1 fitresult_D2_o.b11 fitresult_D2_o.b2 fitresult_D2_o.b22];
                co_no_opti.(my_field).eb = [ fitresult_D2_o.rb0 fitresult_D2_o.rb1 fitresult_D2_o.rb2 ];     % for D2


                % optimized coefficients no offset
                ft = fittype( '(b0*y^rb0) + (b1*y^rb1)*x + (b2*y^rb2)*x^2', 'independent', {'x', 'y'}, 'dependent', 'z' );
                opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                % opts.DiffMinChange = 1e-12;
                % opts.Display = 'Off';
                % opts.MaxFunEvals = 1000;
                % opts.MaxIter = 1000;
                % opts.TolFun = 1e-12;
                % opts.TolX = 1e-12;
                opts.Weights = weights;

                % opts.StartPoint = [0.2 -0.02  0.007 0.001 -1 -0.3];
                % opts.Lower = [0 -1  0 0 -2 -2];
                % opts.Upper = [1  1  1 2  2  0];

                opts.StartPoint =   [0 0 0     0 0 0];
                opts.Lower =        [-1 -1  -1   0 -Inf -Inf];
                opts.Upper =        [1  1  1   Inf Inf 0];

                [fitresult_D2_old, gof] = fit( [xData, yData], zData, ft, opts );
                co_KM_opti_old.(my_field).b  = [fitresult_D2_old.b0 0 fitresult_D2_old.b1 0 fitresult_D2_old.b2 0];
                co_KM_opti_old.(my_field).eb = [fitresult_D2_old.rb0 fitresult_D2_old.rb1 fitresult_D2_old.rb2 ];     % for D2


                %% d_ij(r) plots
                h(3) = figure;
                subplot(2,2,1);
                plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
                    fitresult_D1.a1.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D1.ra1+fitresult_D1.a11,'LineWidth',2)
                hold on
                plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
                    fitresult_D1_o.a1.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D1_o.ra1+fitresult_D1_o.a11,'LineWidth',2)
%                 plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%                     fitresult_D1_old.a1.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D1_old.ra1,'LineWidth',2,'Color','k')
%                 plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%                     (-0.36).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^(-1),'LineWidth',2,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
                axis square
                ylabel('$d_{11}$', 'interpreter','latex')          
                if norm_r==1
                xlabel('$r / \lambda$', 'interpreter','latex')
                else
                xlabel('$r/m $', 'interpreter','latex')
                end 
                xlim([min(y(:,1)) max(y(:,1))])
                set(gca, 'FontSize',12)
                %
                subplot(2,2,2);
                plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
                    fitresult_D2.b0.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D2.rb0+fitresult_D2.b00,'LineWidth',2)
                hold on
                plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
                    fitresult_D2_o.b0.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D2_o.rb0+fitresult_D2_o.b00,'LineWidth',2)

%                 plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%                     fitresult_D2_old.b0.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D2_old.rb0,'LineWidth',2,'Color','k')
                axis square
                ylabel('$d_{20}$', 'interpreter','latex')
                if norm_r==1
                xlabel('$r / \lambda$', 'interpreter','latex')
                else
                xlabel('$r/m $', 'interpreter','latex')
                end 
                xlim([min(y(:,1)) max(y(:,1))])
                set(gca, 'FontSize',12)
                %
                subplot(2,2,3);
                plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
                    fitresult_D2.b1.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D2.rb1+fitresult_D2.b11,'LineWidth',2)
                hold on
                plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
                    fitresult_D2_o.b1.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D2_o.rb1+fitresult_D2_o.b11,'LineWidth',2)
%                 plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%                     fitresult_D2_old.b1.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D2_old.rb1,'LineWidth',2,'Color','k')
                axis square
                ylabel('$d_{21}$', 'interpreter','latex')
                if norm_r==1
                xlabel('$r / \lambda$', 'interpreter','latex')
                else
                xlabel('$r/m $', 'interpreter','latex')
                end 
                xlim([min(y(:,1)) max(y(:,1))])
                set(gca, 'FontSize',12)
                %
                subplot(2,2,4);
                plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
                    fitresult_D2.b2.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D2.rb2+fitresult_D2.b22,'LineWidth',2)
                hold on
                plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
                    fitresult_D2_o.b2.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D2_o.rb2+fitresult_D2_o.b22,'LineWidth',2)
%                 plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%                     fitresult_D2_old.b2.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D2_old.rb2,'LineWidth',2,'Color','k')
%                 plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%                     0.0144.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^(-1),'LineWidth',2,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
                axis square
                ylabel('$d_{22}$', 'interpreter','latex')
                if norm_r==1
                xlabel('$r / \lambda$', 'interpreter','latex')
                else
                xlabel('$r/m $', 'interpreter','latex')
                end 
                xlim([min(y(:,1)) max(y(:,1))])
                set(gca, 'FontSize',14)
                set(gca,'TickLabelInterpreter','latex');
                set(gcf, 'Color', 'w')
                %legend('optimized','non-optimized','Location','NorthWest')
%                 legend({'Optimized','Non-optimized','Optimized no offset','K62'},'Position',[0.45 0.5 0.1 0.05])
                legend({'Optimized','Non-optimized','K62'},'Position',[0.45 0.5 0.1 0.05])

                close all
                disp(['point:' num2str(k) '/'  num2str(length(condition)) ])
            end
    end
else
    evaluated   = evaluated_point;
    xx          = nan(length(evaluated),increment_bin);
    y           = nan(length(evaluated),increment_bin);
    z1          = nan(length(evaluated),increment_bin);
    z1e         = nan(length(evaluated),increment_bin);
    z2          = nan(length(evaluated),increment_bin);
    z2e         = nan(length(evaluated),increment_bin);
    z1_opti     = nan(length(evaluated),increment_bin);
    z2_opti     = nan(length(evaluated),increment_bin);
    z4          = nan(length(evaluated),increment_bin);
    
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
    
%     z1_opti=z1_opti.*0.5;
%     z2_opti=z2_opti.*0.5;

    %TEST Beginn
    % y=y-(((markov./Fs).*m_data)/taylor_L);
    %TEst Ende


    %% D1-Fit
    [fitresult_D1, xData, yData, zData, weights]    = FIT_D1(xx, y, z1_opti, z1e);
    fitresult_D1_conf                               = confint(fitresult_D1);

    h(1) = figure;
    % surf(xx, y, z1_opti,'EdgeColor','none')
    scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
    % scatter3(xData,yData,yData.*zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
    hold on
    surf(repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100),feval(fitresult_D1,repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100)),'EdgeColor','none')
%     surf(repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100),...
%         (-0.36.*repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100).^(-1).*repmat(linspace(min(xData),max(xData),100),100,1)),'FaceColor',[0.501960813999176 0.501960813999176 0.501960813999176],'EdgeColor',[0.501960813999176 0.501960813999176 0.501960813999176])

    % surf(repmat(linspace(-4,4,100),100,1),repmat(linspace(min(y(:,1)),max(y(:,1)),100).',1,100),feval(fitresult_D1,repmat(linspace(-4,4,100),100,1),repmat(linspace(min(y(:,1)),max(y(:,1)),100).',1,100)))
    set(gca,'yScale','log') 
    set(gcf, 'Color', 'w')
    set(gca, 'FontSize',18)
    if norm_ur==1
        xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
    else
        xlabel('$u_r\ (m/s)$', 'interpreter','latex')   
    end
    if norm_r==1
        ylabel('$r / \lambda$', 'interpreter','latex')
    else
        ylabel('$r\ (m)$', 'interpreter','latex')
    end
    zlabel('$D^{(1)} (u_r,r)$', 'interpreter','latex')
    % zlabel('$D1 (x,r) / [\frac{m}{s^2}]$', 'interpreter','latex')
    colormap (parula(200))
    % colorbar
    grid off
    xlim([min(xData) max(xData)])
    ylim([min(y(:,1)) max(y(:,1))])
    % zlim([-0.005 0.005])
    axis square
    set(gca,'YDir','reverse');
%     title('optimized KM, optimized fit','interpreter','latex');
    % title(['$(', num2str(fitresult_D1.a1,'%1.1E'), ' * r / \lambda^{' num2str(fitresult_D1.ra1,'%1.1E') '} ) * u_r $'],'interpreter','latex')   
    % title(['$\left(',num2str(fitresult_D1.a1),'\cdot r/\lambda^{',num2str(fitresult_D1.ra1),'} +',num2str(fitresult_D1.a11),' \right)\cdot u_r $'],'interpreter','latex')   
    rotate3d on
    colorbar('north','TickLabelInterpreter','latex','AxisLocation','out');
    fig_setup
    legend off
    txt = {'(a)'};
    a=text(4,0.5,txt);
    a.Units='normalized';


    %%%%%  2019 %%%%%%%%
    % optimized coefficients 
    % D1(x,y) = (a1*y^ra1+a11)*x
    co_KM_opti.a  = [fitresult_D1.a1 fitresult_D1.a11];
    co_KM_opti.ea = [fitresult_D1.ra1];

    % non-optimized coefficients 
    [fitresult_D1_o, xData_o1, yData_o1, zData_o1, weights] = FIT_D1(xx, y, z1, z1e);
    %%%%%  2019 %%%%%%%%
    % D1(x,y) = (a1*y^ra1+a11)*x
    co_no_opti.a  = [fitresult_D1_o.a1 fitresult_D1_o.a11];
    co_no_opti.ea = [fitresult_D1_o.ra1];

    % optimized coefficients no offset
    % D1(x,y) = (a1*y^ra1)*x
    ft = fittype( '(a1*y^ra1)*x', 'independent', {'x', 'y'}, 'dependent', 'z' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    % opts.DiffMinChange = 1e-12;
    % opts.Display = 'Off';
    % opts.MaxFunEvals = 1000;
    % opts.MaxIter = 1000;
    % opts.TolFun = 1e-12;
    % opts.TolX = 1e-12;
    opts.Weights = weights;

    % opts.StartPoint = [-0.1  -1];
    % opts.Lower = [-2 -2];
    % opts.Upper = [0 0];

    opts.StartPoint = [0  0];
    opts.Lower = [-1 -Inf];
    opts.Upper = [1 0];

    [fitresult_D1_old, gof] = fit( [xData, yData], zData, ft, opts );
    co_KM_opti_old.a  = [fitresult_D1_old.a1 0];
    co_KM_opti_old.ea = [fitresult_D1_old.ra1];

    % figure
    % plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),fitresult_D1.a0.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D1.ra0)
    % hold on
    % plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),fitresult_D1.a1.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D1.ra1)
    % plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),fitresult_D1.a2.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^0)
    % plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),fitresult_D1.a3.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D1.ra3)
    % set(gca,'xScale','log') 
    % axis square
    % set(gcf, 'Color', 'w')
    % % grid on
    % xlabel('$r / \lambda$', 'interpreter','latex')
    % legend('d_{10}(r)','d_{11}(r)','d_{12}(r)','d_{13}(r)')
    % xlim([min(y(:,1)) max(y(:,1))])



    %% D2-Fit
    [fitresult_D2, xData, yData, zData, weights] = FIT_D2(xx, y, z2_opti, z2e);
    fitresult_D2_conf                            = confint(fitresult_D2);
    h(2) = figure;
    % surf(xx, y, z2_opti,'EdgeColor','none')
    scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
    % scatter3(xData,yData,yData.*zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
    hold on
    surf(repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100),feval(fitresult_D2,repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100)),'EdgeColor','none')
%     surf(repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100),...
%         (0.0144.*repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100).^(-1).*repmat(linspace(min(xData),max(xData),100),100,1).^(2)),'FaceColor',[0.501960813999176 0.501960813999176 0.501960813999176],'EdgeColor',[0.501960813999176 0.501960813999176 0.501960813999176])
    % surf(repmat(linspace(-4,4,100),100,1),repmat(linspace(min(y(:,1)),max(y(:,1)),100).',1,100),feval(fitresult_D1,repmat(linspace(-4,4,100),100,1),repmat(linspace(min(y(:,1)),max(y(:,1)),100).',1,100)))
    set(gca,'yScale','log') 
    % set(gca,'zScale','log')  
    set(gcf, 'Color', 'w')
    set(gca, 'FontSize',18)
    if norm_ur==1
        xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
    else
        xlabel('$u_r\ (m/s)$', 'interpreter','latex')   
    end
    if norm_r==1
        ylabel('$r / \lambda$', 'interpreter','latex')
    else
        ylabel('$r\ (m)$', 'interpreter','latex')
    end
    % zlabel('$D2 (x,r) / [\frac{m^2}{s^3}]$', 'interpreter','latex')
    zlabel('$D^{(2)} (u_r,r)$', 'interpreter','latex')
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
    set(gca,'TickLabelInterpreter','latex');
%     title('optimized KM, optimized fit','interpreter','latex');
    % title(['$(', num2str(fitresult_D1.a1,'%1.1E'), ' * r / \lambda^{' num2str(fitresult_D1.ra1,'%1.1E') '} ) * u_r $'],'interpreter','latex')   
    % title(['$\left(',num2str(fitresult_D1.a1),'\cdot r/\lambda^{',num2str(fitresult_D1.ra1),'} +',num2str(fitresult_D1.a11),' \right)\cdot u_r $'],'interpreter','latex')   
    rotate3d on
    colorbar('north','TickLabelInterpreter','latex','AxisLocation','out');
    fig_setup
    legend off
    txt = {'(b)'};
    b=text(4,0.5,txt);
    b.Units='normalized';

    font=22;
    a.FontSize = font;
    b.FontSize = font;
    pos_txt=[-0.22   0.8];
    a.Position=pos_txt;
    b.Position=[-0.25  0.8];

    % optimized coefficients 
    %%%%%  2019 %%%%%%%%
    % D2(x,y) = (b0*y^rb0+b00) + (b1*y^rb1+b11)*x + (b2*y^rb2+b22)*x^2
    co_KM_opti.b  = [fitresult_D2.b0 fitresult_D2.b00 fitresult_D2.b1 fitresult_D2.b11 fitresult_D2.b2 fitresult_D2.b22];
    co_KM_opti.eb = [fitresult_D2.rb0 fitresult_D2.rb1 fitresult_D2.rb2 ];     % for D2

    % non-optimized coefficients 
    [fitresult_D2_o, xData_o2, yData_o2, zData_o2, weights] = FIT_D2(xx, y, z2, z2e);
    %%%%%  2019 %%%%%%%%
    % D2(x,y) = (b0*y^rb0+b00) + (b1*y^rb1+b11)*x + (b2*y^rb2+b22)*x^2
    co_no_opti.b  = [fitresult_D2_o.b0 fitresult_D2_o.b00 fitresult_D2_o.b1 fitresult_D2_o.b11 fitresult_D2_o.b2 fitresult_D2_o.b22];
    co_no_opti.eb = [ fitresult_D2_o.rb0 fitresult_D2_o.rb1 fitresult_D2_o.rb2 ];     % for D2


    % optimized coefficients no offset
    ft = fittype( '(b0*y^rb0) + (b1*y^rb1)*x + (b2*y^rb2)*x^2', 'independent', {'x', 'y'}, 'dependent', 'z' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    % opts.DiffMinChange = 1e-12;
    % opts.Display = 'Off';
    % opts.MaxFunEvals = 1000;
    % opts.MaxIter = 1000;
    % opts.TolFun = 1e-12;
    % opts.TolX = 1e-12;
    opts.Weights = weights;

    % opts.StartPoint = [0.2 -0.02  0.007 0.001 -1 -0.3];
    % opts.Lower = [0 -1  0 0 -2 -2];
    % opts.Upper = [1  1  1 2  2  0];

    opts.StartPoint =   [0 0 0     0 0 0];
    opts.Lower =        [-1 -1  -1   0 -Inf -Inf];
    opts.Upper =        [1  1  1   Inf Inf 0];

    [fitresult_D2_old, gof] = fit( [xData, yData], zData, ft, opts );
    co_KM_opti_old.b  = [fitresult_D2_old.b0 0 fitresult_D2_old.b1 0 fitresult_D2_old.b2 0];
    co_KM_opti_old.eb = [fitresult_D2_old.rb0 fitresult_D2_old.rb1 fitresult_D2_old.rb2 ];     % for D2


    %% d_ij(r) plots
    h(3) = figure;
    subplot(2,2,1);
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        fitresult_D1.a1.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D1.ra1+fitresult_D1.a11,'LineWidth',2)
    hold on
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        fitresult_D1_o.a1.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D1_o.ra1+fitresult_D1_o.a11,'LineWidth',2)
%     plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%         fitresult_D1_old.a1.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D1_old.ra1,'LineWidth',2,'Color','k')
%     plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%         (-0.36).*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^(-1),'LineWidth',2,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
    axis square
    ylabel('$d_{11}$', 'interpreter','latex')
    if norm_r==1
    xlabel('$r / \lambda$', 'interpreter','latex')
    else
    xlabel('$r/m $', 'interpreter','latex')
    end 
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',14)
    fig_setup
    legend off
    txt = {'(a)'};
    a=text(4,0.5,txt);
    a.Units='normalized';
    %
    subplot(2,2,2);
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        fitresult_D2.b0.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D2.rb0+fitresult_D2.b00,'LineWidth',2)
    hold on
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        fitresult_D2_o.b0.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D2_o.rb0+fitresult_D2_o.b00,'LineWidth',2)

%     plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%         fitresult_D2_old.b0.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D2_old.rb0,'LineWidth',2,'Color','k')
    axis square
    ylabel('$d_{20}$', 'interpreter','latex')
    if norm_r==1
        xlabel('$r / \lambda$', 'interpreter','latex')
    else
        xlabel('$r\ (m)$', 'interpreter','latex')
    end
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',14)
    fig_setup
    legend off
    txt = {'(b)'};
    b=text(4,0.5,txt);
    b.Units='normalized';

    %
    subplot(2,2,3);
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        fitresult_D2.b1.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D2.rb1+fitresult_D2.b11,'LineWidth',2)
    hold on
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        fitresult_D2_o.b1.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D2_o.rb1+fitresult_D2_o.b11,'LineWidth',2)
%     plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%         fitresult_D2_old.b1.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D2_old.rb1,'LineWidth',2,'Color','k')
    axis square
    ylabel('$d_{21}$', 'interpreter','latex')
    if norm_r==1
        xlabel('$r / \lambda$', 'interpreter','latex')
    else
        xlabel('$r\ (m)$', 'interpreter','latex')
    end
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',14)
    fig_setup
    legend off
    txt = {'(c)'};
    c=text(4,0.5,txt);
    c.Units='normalized';


    %
    subplot(2,2,4);
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        fitresult_D2.b2.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D2.rb2+fitresult_D2.b22,'LineWidth',2)
    hold on
    plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
        fitresult_D2_o.b2.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D2_o.rb2+fitresult_D2_o.b22,'LineWidth',2)
%     plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%         fitresult_D2_old.b2.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^fitresult_D2_old.rb2,'LineWidth',2,'Color','k')
%     plot(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100),...
%         0.0144.*logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).^(-1),'LineWidth',2,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
    axis square
    ylabel('$d_{22}$', 'interpreter','latex')
    if norm_r==1
        xlabel('$r / \lambda$', 'interpreter','latex')
    else
        xlabel('$r\ (m)$', 'interpreter','latex')
    end
    xlim([min(y(:,1)) max(y(:,1))])
    set(gca, 'FontSize',14)
    set(gca,'TickLabelInterpreter','latex');
    set(gcf, 'Color', 'w')
    %legend('optimized','non-optimized','Location','NorthWest')
%     legend({'Optimized','Non-optimized','Optimized no offset','K62'},'Position',[0.45 0.5 0.1 0.05])
    legend({'Optimized','Non-optimized','K62'},'Position',[0.45 0.5 0.1 0.05])
    fig_setup
% legend off
    txt = {'(d)'};
    d=text(4,0.5,txt);
    d.Units='normalized';

    font=22;
    a.FontSize = font;
    b.FontSize = font;
    c.FontSize = font;
    d.FontSize = font;

    pos_txt=[-0.68   0.9];
    a.Position=pos_txt;
    b.Position=pos_txt;
    c.Position=pos_txt;
    d.Position=pos_txt;

    %% 3D plots
    figure
    % surf(xx, y, z1_opti,'EdgeColor','none')
    scatter3(xData_o1, yData_o1, zData_o1,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
    % scatter3(xData,yData,yData.*zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
    hold on
    surf(repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100),feval(fitresult_D1,repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100)),'EdgeColor','none')
    % surf(repmat(linspace(-4,4,100),100,1),repmat(linspace(min(y(:,1)),max(y(:,1)),100).',1,100),feval(fitresult_D1,repmat(linspace(-4,4,100),100,1),repmat(linspace(min(y(:,1)),max(y(:,1)),100).',1,100)))
    set(gca,'yScale','log') 
    set(gcf, 'Color', 'w')
    set(gca, 'FontSize',18)
    if norm_ur==1
        xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
    else
        xlabel('$u_r\ (m/s)$', 'interpreter','latex')   
    end
    if norm_r==1
        ylabel('$r / \lambda$', 'interpreter','latex')
    else
        ylabel('$r\ (m)$', 'interpreter','latex')
    end
    zlabel('$D^{(1)} (u_r,r)$', 'interpreter','latex')
    % zlabel('$D1 (x,r) / [\frac{m}{s^2}]$', 'interpreter','latex')
    colormap (parula(200))
    % colorbar
    grid off
    xlim([min(xData) max(xData)])
    ylim([min(y(:,1)) max(y(:,1))])
    % zlim([-0.005 0.005])
    axis square
    set(gca,'YDir','reverse');
    set(gca,'TickLabelInterpreter','latex');
    title('non-optimized KM, optimized fit','interpreter','latex');
    close


    figure
    % surf(xx, y, z2_opti,'EdgeColor','none')
    scatter3(xData_o2,yData_o2,zData_o2,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
    % scatter3(xData,yData,yData.*zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
    hold on
    surf(repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100),feval(fitresult_D2,repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100)),'EdgeColor','none')
    % surf(repmat(linspace(-4,4,100),100,1),repmat(linspace(min(y(:,1)),max(y(:,1)),100).',1,100),feval(fitresult_D1,repmat(linspace(-4,4,100),100,1),repmat(linspace(min(y(:,1)),max(y(:,1)),100).',1,100)))
    set(gca,'yScale','log') 
    % set(gca,'zScale','log')  
    set(gcf, 'Color', 'w')
    set(gca, 'FontSize',18)
    if norm_ur==1
        xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
    else
        xlabel('$u_r\ (m/s)$', 'interpreter','latex')   
    end
    if norm_r==1
        ylabel('$r / \lambda$', 'interpreter','latex')
    else
        ylabel('$r\ (m)$', 'interpreter','latex')
    end
    % zlabel('$D2 (x,r) / [\frac{m^2}{s^3}]$', 'interpreter','latex')
    zlabel('$D^{(2)} (u_r,r)$', 'interpreter','latex')
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
    set(gca,'TickLabelInterpreter','latex');
    % title(['$(', num2str(fitresult_D1.a1,'%1.1E'), ' * r / \lambda^{' num2str(fitresult_D1.ra1,'%1.1E') '} ) * u_r $'],'interpreter','latex')   
    % title(['$\left(',num2str(fitresult_D1.a1),'\cdot r/\lambda^{',num2str(fitresult_D1.ra1),'} +',num2str(fitresult_D1.a11),' \right)\cdot u_r $'],'interpreter','latex')   
    title('non-optimized KM, optimized fit','interpreter','latex');
    close

    figure
    scatter3(xData_o1, yData_o1, zData_o1,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
    hold on
    surf(repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100),feval(fitresult_D1_o,repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100)),'EdgeColor','none')
    set(gca,'yScale','log') 
    set(gcf, 'Color', 'w')
    set(gca, 'FontSize',18)
    if norm_ur==1
        xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
    else
        xlabel('$u_r\ (m/s)$', 'interpreter','latex')   
    end
    if norm_r==1
        ylabel('$r / \lambda$', 'interpreter','latex')
    else
        ylabel('$r\ (m)$', 'interpreter','latex')
    end
    zlabel('$D^{(1)} (u_r,r)$', 'interpreter','latex')
    colormap (parula(200))
    grid off
    xlim([min(xData) max(xData)])
    ylim([min(y(:,1)) max(y(:,1))])
    axis square
    set(gca,'YDir','reverse');
    set(gca,'TickLabelInterpreter','latex');
    title('non-optimized KM, non-optimized fit','interpreter','latex');
    close


    figure
    scatter3(xData_o2,yData_o2,zData_o2,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
    hold on
    surf(repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100),feval(fitresult_D2_o,repmat(linspace(min(xData),max(xData),100),100,1),repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100)),'EdgeColor','none')
    set(gca,'yScale','log') 
    set(gcf, 'Color', 'w')
    set(gca, 'FontSize',18)
    if norm_ur==1
        xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
    else
        xlabel('$u_r\ (m/s)$', 'interpreter','latex')   
    end
    if norm_r==1
        ylabel('$r / \lambda$', 'interpreter','latex')
    else
        ylabel('$r\ (m)$', 'interpreter','latex')
    end
    zlabel('$D^{(2)} (u_r,r)$', 'interpreter','latex')
    colormap (parula(200))
    grid off
    xlim([min(xData) max(xData)])
    ylim([min(y(:,1)) max(y(:,1))])
    axis square
    set(gca,'YDir','reverse');
    set(gca,'TickLabelInterpreter','latex');
    title('non-optimized KM, non-optimized fit','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    close
    
    
    

    z1=nan(length(evaluated),increment_bin);
    z2=nan(length(evaluated),increment_bin);
    xx_no_opti  = nan(length(evaluated),increment_bin);
    
    for i=1:length(evaluated)
%         xx_no_opti(i,:)=evaluated(i).y_mean_bin;
%         z1(i,:)=(evaluated(i).D1);
%         z2(i,:)=(evaluated(i).D2);
        xx_no_opti=xx;
        z1(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).D1(evaluated(i).x_bin_not_nan));
        z2(i,1:sum(evaluated(i).x_bin_not_nan))=(evaluated(i).D2(evaluated(i).x_bin_not_nan));
    end
    
    [xData, yData, zData] = prepareSurfaceData( xx, y, z1_opti);
    tmp_u=repmat(linspace(min(xData),max(xData),100),100,1);
    tmp_r=repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100);

    h(4) = figure;
%     scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
    scatter3(xData,yData,zData,'MarkerFaceColor',[0    0.4470    0.7410],'MarkerEdgeColor',[1 1 1])
    hold on
    [xData, yData, zData] = prepareSurfaceData( xx_no_opti, y, z1);
    scatter3(xData,yData,zData,'MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerEdgeColor',[1 1 1])
    set(gca,'yScale','log') 
    set(gcf, 'Color', 'w')
    set(gca, 'FontSize',18)
    if norm_ur==1
        xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
    else
        xlabel('$u_r\ (m/s)$', 'interpreter','latex')   
    end
    if norm_r==1
        ylabel('$r / \lambda$', 'interpreter','latex')
    else
        ylabel('$r\ (m)$', 'interpreter','latex')
    end
    zlabel('$D^{(1)} (u_r,r)$', 'interpreter','latex')
    % zlabel('$D1 (x,r) / [\frac{m}{s^2}]$', 'interpreter','latex')
    colormap (parula(200))
    % colorbar
    grid off
    xlim([min(xData) max(xData)])
    ylim([min(y(:,1)) max(y(:,1))])
    % zlim([-0.005 0.005])
    axis square
    set(gca,'YDir','reverse');
    rotate3d on
    legend({'Optimized','Non-optimized'},'Position',[0.5 0.7 0.3 0.1])



    [xData, yData, zData] = prepareSurfaceData( xx, y, z2_opti);
    tmp_u=repmat(linspace(min(xData),max(xData),100),100,1);
    tmp_r=repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100);
    h(5) = figure;
%     scatter3(xData,yData,zData,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1])
    scatter3(xData,yData,zData,'MarkerFaceColor',[0    0.4470    0.7410],'MarkerEdgeColor',[1 1 1])
    hold on
    [xData, yData, zData] = prepareSurfaceData( xx_no_opti, y, z2);
    tmp_u=repmat(linspace(min(xData),max(xData),100),100,1);
    tmp_r=repmat(logspace(log10(min(y(:,1))),log10(max(y(:,1))),100).',1,100);
    scatter3(xData,yData,zData,'MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerEdgeColor',[1 1 1])
    % surf(tmp_u,tmp_r,...
    %     (d20_R+d21_R.*tmp_u+d22_R.*tmp_u.^2)...
    %     ,'FaceColor',[0.831372559070587 0.815686285495758 0.7843137383461])
    set(gca,'yScale','log') 
    % set(gca,'zScale','log')  
    set(gcf, 'Color', 'w')
    set(gca, 'FontSize',18)
    if norm_ur==1
        xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
    else
        xlabel('$u_r\ (m/s)$', 'interpreter','latex')   
    end
    if norm_r==1
        ylabel('$r / \lambda$', 'interpreter','latex')
    else
        ylabel('$r\ (m)$', 'interpreter','latex')
    end
    % zlabel('$D2 (x,r) / [\frac{m^2}{s^3}]$', 'interpreter','latex')
    zlabel('$D^{(2)} (u_r,r)$', 'interpreter','latex')
    colormap (parula(200))
    % colorbar
    grid off
    xlim([min(xData) max(xData)])
    ylim([min(y(:,1)) max(y(:,1))])
    % zlim([0 2*max(zData)])
    caxis([0 max(zData)])
    axis square
    % close
    set(gca,'YDir','reverse');
    rotate3d on
    legend({'Optimized','Non-optimized'},'Position',[0.5 0.7 0.3 0.1])
    
    if ischar(save_path)
        savefig(h,fullfile(save_path,append(save_name,'_','KM_Fit.fig')),'compact')
        for a = 1:length(h)
            exportgraphics(h(a),fullfile(save_path,append(save_name,'_',sprintf('KM_Fit_%d.png', a))))
        end
    end 
    
end
for i=1:5
    tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
    uiwait(gcf);
    delete(tmp_ui);
    close
end
end