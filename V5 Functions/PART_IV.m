%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% PART IV: Consistency check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,ikx,bar,tmp_bar,save_path,save_name,Fs,m_data,kin_vis,increment_bin,f, E_f_no_filter, f_avg_no_filter, E_avg_no_filter, P_avg_no_filter,filter,...
    low_freq,data_filter,f_avg_filter,E_f_avg_filter,k_avg_filter,Ek_avg_filter,r_avg_filter,Dr_avg_filter,Dk_avg_filter,...
    lenght_scales,C2,int_L,taylor_L,int_L_calc, taylor_L_calc,epsi,epsi_calc,diss_scale,Ce, Ce_calc, Re, Re_lambda,...
    flip_data,norm_ur,norm_r,siginf,r_struc,S_exp_2,S_exp_3,S_exp_4,S_exp_5,S_exp_6,S_exp_7,...
    markov,min_events,var_markov,scale_steps,multi_point,condition,tol,evaluated,step_con_moment,plot_moment,...
    test_opti,tol_opti,var_opti,co_KM_opti,co_KM_opti_no_offset,co_KM_non_opti,fitresult_D1_conf,fitresult_D2_conf,...
    trajec,z,dr_ind,data_length,Sm,Ds,DS,r,ind_trajec,rind,dr,r_s,A,iter,tol_D1,tol_D2,co_IFT_opti,history]=PART_IV(...
    data,ikx,bar,tmp_bar,save_path,save_name,Fs,m_data,kin_vis,increment_bin,f, E_f_no_filter, f_avg_no_filter, E_avg_no_filter, P_avg_no_filter,filter,...
    low_freq,data_filter,f_avg_filter,E_f_avg_filter,k_avg_filter,Ek_avg_filter,r_avg_filter,Dr_avg_filter,Dk_avg_filter,...
    lenght_scales,C2,int_L,taylor_L,int_L_calc, taylor_L_calc,epsi,epsi_calc,diss_scale,Ce, Ce_calc, Re, Re_lambda,...
    flip_data,norm_ur,norm_r,siginf,r_struc,S_exp_2,S_exp_3,S_exp_4,S_exp_5,S_exp_6,S_exp_7,...
    markov,min_events,var_markov,scale_steps,multi_point,condition,tol,evaluated,step_con_moment,plot_moment,...
    test_opti,tol_opti,var_opti,co_KM_opti,co_KM_opti_no_offset,co_KM_non_opti,fitresult_D1_conf,fitresult_D2_conf,...
	trajec,z,dr_ind,data_length,Sm,Ds,DS,r,ind_trajec,rind,dr,r_s,A,iter,tol_D1,tol_D2,co_IFT_opti,history)

	bar=waitbar(22/ikx,'Consistency check: Reconstruction of pdfs & structure functions','Position', tmp_bar);

	answer = questdlg('Select the Kramers-Moyal coefficients to be used for the consistency check:','Consistency check', ...
     	'non-optimized', ...
     	'optimized (conditional PDF)', ...
     	'optimized (IFT)', ...
     	'optimized (IFT)');
	switch answer
    	case 'non-optimized'
       	co = co_KM_non_opti;
       	co_recon = 1;
    	case 'optimized (conditional PDF)'
       	co = co_KM_opti; 
       	co_recon = 2;
    	case 'optimized (IFT)'
       	co = co_IFT_opti;   
       	co_recon = 3;
	end

	iter_recon  = askYesno('Do you want the reconstruction to be iterative (see readme)?', 'Yes');
		
	taylor_L_sample    = ceil(Fs*taylor_L/m_data);
	
	tau2_total    = flip(ceil(Fs*int_L/m_data) : -markov : ceil(Fs*taylor_L/m_data));
	
	if tau2_total(1)-markov < evaluated(1).r_short_sample
    	tau2_total = tau2_total(2:end);
	end
	scale_steps_recon   = length(tau2_total);
	
	for i = scale_steps_recon:-1:1
    	i
    	tau1            = tau2_total(i)-markov;                   %min(tau1) is Taylor length or Markov length
    	tau2            = tau2_total(i);                 %max(tau2) is int_L% put one markov more to start with at L
	
    	[incr1,incr2]   = Increment(tau1,tau2,data_filter);
	%     [fuL, uL]       = hist(incr1,increment_bin);
	%     PA_Gauss(i,:) = pdf(fitdist(incr1.','Normal'),uL);
	%     clear uL fuL
	
    	[dx,dy,x,y,dA]  = limit(incr1,incr2,increment_bin); 
    	r__recon(i)            = (tau1./Fs).*m_data;
    	samples(i)      = tau1;
	
    	[P_AIB,P_BIA,P_AnB,P_A,P_B,binP_AIB,x_mean_bin,y_mean_bin,events,counter_A,counter_B] = distribution(tau1,tau2,incr1,incr2,increment_bin,dx,dy,x,y,dA,1);
	%     [P_AIB,P_BIA,P_AnB,P_A,P_B,binP_AIB,x_mean_bin,y_mean_bin,events,counter_A,counter_B] = distribution(tau1,tau2,incr1,incr2,increment_bin,dx,dy,x,y,dA);
	%     y = y_mean_bin;
	%     x = x_mean_bin;
    	
    	%P_A_STP            = (P_AIB*P_B.'.*dy).';%From experiment
    	%% this is to use the PDF from STP for reconstruction at next smaller scale
    	if iter_recon==1
        	if i==scale_steps_recon
            	P_B = P_B;
        	else
            	P_B = PASTP(i+1,:);
            	clear P_A_STP
        	end
    	end
    	
    	%%%%%%%%%%%%%%%%%%%%%%%
    	r_2     =   ((tau2./Fs).*m_data)/taylor_L;
    	D1_poly =    (co.a(1).*r_2.^co.ea(1)+co.a(2)).*y;
    	D2_poly =   ((co.b(1).*r_2.^co.eb(1)+co.b(2)).*y.^0)+...
                	((co.b(3).*r_2.^co.eb(2)+co.b(4)).*y.^1)+...
                	((co.b(5).*r_2.^co.eb(3)+co.b(6)).*y.^2);
            	
    	P_n_opti_func      = ShortTimeProp(y,x,tau1,tau2,D1_poly,D2_poly,taylor_L,m_data,Fs,norm_ur,norm_r);     
	
	
	%     sum(sum(isnan(P_n_opti_func)))
    	P_n_opti_func(isnan(P_n_opti_func))=0;
    	
	%     size(y)
	%     size(x)
	%     size(P_n_opti_func)
	%     size(P_B.')
    	
	%     x_A_STP            = x;
    	P_A_STP            = (P_n_opti_func*P_B.'.*dy).';
	%     size(P_A_STP)
    	%%%%%%%%%%%%%%%%%%%%%%
	
    	% % % % figure(1000+i)
    	% % % % %plot(y,P_B.','LineWidth',2)%from experiment
    	% % % % hold on
    	% % % % plot(x,P_A.','LineWidth',2)%from experiment
    	% % % % %     P_B_STP(k,:)            = (P_BIA*P_A.'.*dx).';
    	% % % % %     plot(x,P_B_STP(k,:),'k','LineWidth',2)
    	% % % % plot(x_A_STP,P_A_STP,'r--','LineWidth',2)%from STP
    	% % % % set(gca,'yScale','log')
    	% % % % set(gcf, 'Color', 'w')
    	% % % % set(gca,'FontSize',22)
    	% % % % xlabel('$u_{r} / \sigma_\infty$', 'interpreter','latex')
    	% % % % ylabel('$p(u_{r} / \sigma_\infty)$','interpreter','latex')
    	% % % % axis square
    	
	
    	
	
    	X(i,:)          = x;
    	PA(i,:)         = P_A;
    	PASTP(i,:)      = P_A_STP;
    	rms_incr1(i)    = rms(incr1.*siginf);
    	
	%     scale_r(i)      = (samples(i).*m_data)./Fs;
    	
	%     S2_exp(i)   = nansum(((X(i,:)).^2).*PA(i,:)     .*mean(diff(X(i,:))));
	%     S3_exp(i)   = nansum(((X(i,:)).^3).*PA(i,:)     .*mean(diff(X(i,:))));
	%     S2_recon(i) = nansum(((X(i,:)).^2).*PASTP(i,:)  .*mean(diff(X(i,:))));
	%     S3_recon(i) = nansum(((X(i,:)).^3).*PASTP(i,:)  .*mean(diff(X(i,:))));
	
	
	
    	P_A_STP                 = P_A_STP(~isnan(x));
    	x                       = x(~isnan(x));
	
    	S_KM_ST_opti_func_2(i)  = trapz(x,(x.^2.*P_A_STP));
    	S_KM_ST_opti_func_3(i)  = trapz(x,(x.^3.*P_A_STP));
    	S_KM_ST_opti_func_4(i)  = trapz(x,(x.^4.*P_A_STP));
    	S_KM_ST_opti_func_5(i)  = trapz(x,(x.^5.*P_A_STP));
    	S_KM_ST_opti_func_6(i)  = trapz(x,(x.^6.*P_A_STP));  
    	S_KM_ST_opti_func_7(i)  = trapz(x,(x.^7.*P_A_STP));  
    	
    	clear x P_A P_B P_n_opti_func
	end
	
	
	% d11 = -0.362;
	% d22 = 0.0144;
	d11 = (-1/3);
	d22 = 0;
        
%     D1_poly =    (co.a(1).*r_2.^co.ea(1)+co.a(2)).*y;
%     D2_poly =   ((co.b(1).*r_2.^co.eb(1)+co.b(2)).*y.^0)+...
%                 	((co.b(3).*r_2.^co.eb(2)+co.b(4)).*y.^1)+...
%                 	((co.b(5).*r_2.^co.eb(3)+co.b(6)).*y.^2);
	

	for i_stru=length(r_struc):-1:2
    	
    	if i_stru==length(S_exp_3)
        	S_exp_3_recon(i_stru)=S_exp_3(i_stru);
    	end
    	
    	r_tmp=r_struc(i_stru-1)*taylor_L;
	%     r_tmp=1;
    	
     	if iter_recon==1
         	S_exp_3_recon(i_stru-1) = (...
                                        	(-3*(d11/r_tmp + 2*d22/r_tmp)*S_exp_3_recon(i_stru))*...
                                        	((r_struc(i_stru-1)-r_struc(i_stru))*taylor_L)...
                                        	)...
                                        	+S_exp_3_recon(i_stru);
     	else
         	S_exp_3_recon(i_stru-1) = (...
                                        	(-3*(d11/r_tmp + 2*d22/r_tmp)*S_exp_3(i_stru))*...
                                        	((r_struc(i_stru-1)-r_struc(i_stru))*taylor_L)...
                                        	)...
                                        	+S_exp_3(i_stru);
     	end
	end
	S_exp_3_recon(end)=nan;
	% 
	% 
	% fontsize=22;
	% line=2;
	% 
	% h(1) = figure;
	% set(gcf, 'Color', 'w')
	% subplot(3,2,1)
	% plot(r_struc,S_exp_2,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
	% hold on
	% plot(r./taylor_L,S_KM_ST_opti_func_2,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 
	% set(gca,'xScale','log') 
	% set(gca,'yScale','log') 
	% set(gcf, 'Color', 'w')
	% set(gca,'FontSize',fontsize)
	% xlabel('$r/\lambda$', 'interpreter','latex')
	% % xlabel('$r / \lambda$', 'interpreter','latex')
	% ylabel('$S^2(r)$', 'interpreter','latex')
	% xlim([min(r_struc)*0.5 max(r_struc)*2])
	% % ylim([min(S_exp_2) max(S_exp_2)])
	% vline(evaluated(1).r/taylor_L,'r')
	% vline(taylor_L/taylor_L,'k')
	% vline(int_L/taylor_L,'k')
	% set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0]);
	% axis square
	% 
	% 
	% subplot(3,2,2)
	% plot(r_struc,-S_exp_3,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
	% hold on
	% plot(r./taylor_L,-S_KM_ST_opti_func_3,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 
	% set(gca,'xScale','log')
	% set(gca,'yScale','log') 
	% set(gcf, 'Color', 'w')
	% set(gca,'FontSize',fontsize)
	% % xlabel('$r$', 'interpreter','latex')
	% xlabel('$r / \lambda$', 'interpreter','latex')
	% ylabel('$-S^3(r)$', 'interpreter','latex')
	% xlim([min(r_struc)*0.5 max(r_struc)*2])
	% % ylim([min(-S_exp_3) max(-S_exp_3)])
	% vline(evaluated(1).r/taylor_L,'r')
	% vline(taylor_L/taylor_L,'k')
	% vline(int_L/taylor_L,'k')
	% set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0]);
	% axis square
	% 
	% 
	% subplot(3,2,3)
	% plot(r_struc,S_exp_4,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
	% hold on
	% plot(r./taylor_L,S_KM_ST_opti_func_4,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
	% set(gca,'xScale','log')
	% set(gca,'yScale','log') 
	% set(gcf, 'Color', 'w')
	% set(gca,'FontSize',fontsize)
	% % xlabel('$r$', 'interpreter','latex')
	% xlabel('$r / \lambda$', 'interpreter','latex')
	% ylabel('$S^4(r)$', 'interpreter','latex')
	% xlim([min(r_struc)*0.5 max(r_struc)*2])
	% % ylim([min(S_exp_4) max(S_exp_4)])
	% vline(evaluated(1).r/taylor_L,'r')
	% vline(taylor_L/taylor_L,'k')
	% vline(int_L/taylor_L,'k')
	% set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0]);
	% axis square
	% 
	% subplot(3,2,4)
	% plot(r_struc,-S_exp_5,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
	% hold on
	% plot(r./taylor_L,-S_KM_ST_opti_func_5,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 
	% set(gca,'xScale','log')
	% set(gca,'yScale','log') 
	% set(gcf, 'Color', 'w')
	% set(gca,'FontSize',fontsize)
	% % xlabel('$r$', 'interpreter','latex')
	% xlabel('$r / \lambda$', 'interpreter','latex')
	% ylabel('$-S^5(r)$', 'interpreter','latex')
	% xlim([min(r_struc)*0.5 max(r_struc)*2])
	% % ylim([min(-S_exp_5) max(-S_exp_5)])
	% vline(evaluated(1).r/taylor_L,'r')
	% vline(taylor_L/taylor_L,'k')
	% vline(int_L/taylor_L,'k')
	% set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0]);
	% axis square
	% 
	% subplot(3,2,5)
	% plot(r_struc,S_exp_6,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
	% hold on
	% plot(r./taylor_L,S_KM_ST_opti_func_6,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
	% set(gca,'xScale','log')
	% set(gca,'yScale','log') 
	% set(gcf, 'Color', 'w')
	% set(gca,'FontSize',fontsize)
	% % xlabel('$r$', 'interpreter','latex')
	% xlabel('$r / \lambda$', 'interpreter','latex')
	% ylabel('$S^6(r)$', 'interpreter','latex')
	% xlim([min(r_struc)*0.5 max(r_struc)*2])
	% % ylim([min(S_exp_6) max(S_exp_6)])
	% vline(evaluated(1).r/taylor_L,'r')
	% vline(taylor_L/taylor_L,'k')
	% vline(int_L/taylor_L,'k')
	% set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0]);
	% axis square 
	% 
	% subplot(3,2,6)
	% plot(r_struc,-S_exp_7,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
	% hold on
	% plot(r./taylor_L,-S_KM_ST_opti_func_7,'LineWidth',line,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 
	% set(gca,'xScale','log')
	% set(gca,'yScale','log') 
	% set(gcf, 'Color', 'w')
	% set(gca,'FontSize',fontsize)
	% % xlabel('$r$', 'interpreter','latex')
	% xlabel('$r / \lambda$', 'interpreter','latex')
	% ylabel('$-S^7(r)$', 'interpreter','latex')
	% xlim([min(r_struc)*0.5 max(r_struc)*2])
	% % ylim([min(-S_exp_7) max(-S_exp_7)])
	% vline(evaluated(1).r/taylor_L,'r')
	% vline(taylor_L/taylor_L,'k')
	% vline(int_L/taylor_L,'k')
	% set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0]);
	% axis square
	% % set(h(1),'Position',[1 46 835 1299]);
	
	
	
	
	
	%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	% [0.501960813999176 0.501960813999176 0.501960813999176]
	fontsize=22;
	line=2;
	h(1) = figure;
	subplot(1,2,1)
	p1=plot(r_struc,S_exp_2,'MarkerSize',8,'Marker','+','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
	hold on
	p2=plot(r_struc,S_exp_4,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
	p3=plot(r_struc,S_exp_6,'MarkerSize',8,'Marker','*','LineWidth',2,'LineStyle','none','Color',[0 0 0]); 
	plot(r__recon./taylor_L,S_KM_ST_opti_func_2,'r','LineWidth',line) 
	plot(r__recon./taylor_L,S_KM_ST_opti_func_4,'r','LineWidth',line) 
	plot(r__recon./taylor_L,S_KM_ST_opti_func_6,'r','LineWidth',line) 
	set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
	set(gca,'xScale','log')
	set(gca,'yScale','log') 
	axis square
	set(gcf, 'Color', 'w')
	set(gca,'FontSize',fontsize)
	xlabel('$r/m$', 'interpreter','latex')
	if norm_r==1
    	xlabel('$r / \lambda$', 'interpreter','latex')
	else
    	xlabel('$r\ (m)$', 'interpreter','latex')
	end
	xlim([min(r_struc)*0.5 max(r_struc)*2])
	% ylim([0.5*min([S_exp_2 S_exp_4 S_exp_6]) 2*max([S_exp_2 S_exp_4 S_exp_6])])
	legend([p1,p2,p3],{'$S^2(r)$','$S^4(r)$','$S^6(r)$'},'Interpreter','latex','Location','northwest','FontSize',18)  
	legend([p1,p2,p3],{'$S^2$','$S^4$','$S^6$'},'Interpreter','latex','Location','northwest','FontSize',18)  
	fig_setup
	txt = {'(a)'};
	a=text(4,0.5,txt);
	a.Units='normalized';
	
	subplot(1,2,2)
	p4=plot(r_struc,-S_exp_3,'MarkerSize',8,'Marker','+','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
	hold on
	p5=plot(r_struc,-S_exp_5,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0]);
	p6=plot(r_struc,-S_exp_7,'MarkerSize',8,'Marker','*','LineWidth',2,'LineStyle','none','Color',[0 0 0]); 
	plot(r__recon./taylor_L,-S_KM_ST_opti_func_3,'r','LineWidth',line)
	plot(r__recon./taylor_L,-S_KM_ST_opti_func_5,'r','LineWidth',line)
	plot(r__recon./taylor_L,-S_KM_ST_opti_func_7,'r','LineWidth',line)
	% plot(r__recon,feval(fitresult,r),'LineWidth',2,'LineStyle','--','Color','r')
% 	plot(r_struc(1:end-1),-S_exp_3_recon(1:end-1),'LineWidth',2)
	set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
	set(gca,'xScale','log')
	set(gca,'yScale','log') 
	axis square
	set(gcf, 'Color', 'w')
	set(gca,'FontSize',fontsize)
	xlabel('$r/m$', 'interpreter','latex')
	if norm_r==1
    	xlabel('$r / \lambda$', 'interpreter','latex')
	else
    	xlabel('$r\ (m)$', 'interpreter','latex')
	end
	xlim([min(r_struc)*0.5 max(r_struc)*2])
	% ylim([0.5*min([-S_exp_3 -S_exp_5 -S_exp_7]) 2*max([-S_exp_3 -S_exp_5 -S_exp_7])])
	legend([p4,p5,p6],{'$-S^3(r)$','$-S^5(r)$','$-S^7(r)$'},'Interpreter','latex','Location','northwest','FontSize',18)  
	legend([p4,p5,p6],{'$-S^3$','$-S^5$','$-S^7$'},'Interpreter','latex','Location','northwest','FontSize',18)  
	% vline(taylor_L,'k')
	% vline(int_L,'k')
	% vline(taylor_L/taylor_L,'k')
	% vline(int_L/taylor_L,'k')
	fig_setup
	txt = {'(b)'};
	b=text(4,0.5,txt);
	b.Units='normalized';
	
	font=22;
	a.FontSize = font;
	b.FontSize = font;
	
	pos_txt=[-0.22   0.9];
	a.Position=pos_txt;
	b.Position=pos_txt;
	
	
	
	
	%% Intermittency plot of velocity increment at least at 6 different scales
	h(2) = figure;
	i_plot(:) = unique(round(linspace(1,size(X,1),6)));
	c = flip(repmat(linspace(0,0.6,length(i_plot)).',1,3));
	
	set(gcf, 'Color', 'w')
	subplot(1,2,1)
	for i=1:length(i_plot)
    	if i==1
        	p1  = semilogy(X(i_plot(i),:),(10^(i-1))*PA(i_plot(i),:).','MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',c(i,:));%from experiment
    	elseif i==length(i_plot)
        	p2  = semilogy(X(i_plot(i),:),(10^(i-1))*PA(i_plot(i),:).','MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',c(i,:));%from experiment 
    	else
        	semilogy(X(i_plot(i),:),(10^(i-1))*PA(i_plot(i),:).','MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',c(i,:))%from experiment 
    	end
    	hold on;
    	set(gcf, 'Color', 'w')
    	axis square
    	if norm_ur==1  
        	xlabel('$u_r / \sigma_\infty$', 'interpreter','latex')
    	else
        	xlabel('$u_r\ (m/s)$','interpreter','latex');
    	end
    	ylabel('$p(u_r / \sigma_\infty)$','interpreter','latex')
    	ylabel('PDF','interpreter','latex')
	%     xlim([-5 5])
    	fig_setup
	end
	% tmp=[0 0 1;0.0199999995529652 0.0199999995529652 1;0.0399999991059303 0.0399999991059303 1;0.0599999986588955 0.0599999986588955 1;0.0799999982118607 0.0799999982118607 1;0.100000001490116 0.100000001490116 1;0.119999997317791 0.119999997317791 1;0.140000000596046 0.140000000596046 1;0.159999996423721 0.159999996423721 1;0.180000007152557 0.180000007152557 1;0.200000002980232 0.200000002980232 1;0.219999998807907 0.219999998807907 1;0.239999994635582 0.239999994635582 1;0.259999990463257 0.259999990463257 1;0.280000001192093 0.280000001192093 1;0.300000011920929 0.300000011920929 1;0.319999992847443 0.319999992847443 1;0.340000003576279 0.340000003576279 1;0.360000014305115 0.360000014305115 1;0.379999995231628 0.379999995231628 1;0.400000005960464 0.400000005960464 1;0.419999986886978 0.419999986886978 1;0.439999997615814 0.439999997615814 1;0.46000000834465 0.46000000834465 1;0.479999989271164 0.479999989271164 1;0.5 0.5 1;0.519999980926514 0.519999980926514 1;0.540000021457672 0.540000021457672 1;0.560000002384186 0.560000002384186 1;0.579999983310699 0.579999983310699 1;0.600000023841858 0.600000023841858 1;0.620000004768372 0.620000004768372 1;0.639999985694885 0.639999985694885 1;0.660000026226044 0.660000026226044 1;0.680000007152557 0.680000007152557 1;0.699999988079071 0.699999988079071 1;0.720000028610229 0.720000028610229 1;0.740000009536743 0.740000009536743 1;0.759999990463257 0.759999990463257 1;0.779999971389771 0.779999971389771 1;0.800000011920929 0.800000011920929 1;0.819999992847443 0.819999992847443 1;0.839999973773956 0.839999973773956 1;0.860000014305115 0.860000014305115 1;0.879999995231628 0.879999995231628 1;0.899999976158142 0.899999976158142 1;0.920000016689301 0.920000016689301 1;0.939999997615814 0.939999997615814 1;0.959999978542328 0.959999978542328 1;0.980000019073486 0.980000019073486 1;1 1 1;1 0.980392158031464 0.980392158031464;1 0.960784316062927 0.960784316062927;1 0.941176474094391 0.941176474094391;1 0.921568632125854 0.921568632125854;1 0.901960790157318 0.901960790157318;1 0.882352948188782 0.882352948188782;1 0.862745106220245 0.862745106220245;1 0.843137264251709 0.843137264251709;1 0.823529422283173 0.823529422283173;1 0.803921580314636 0.803921580314636;1 0.7843137383461 0.7843137383461;1 0.764705896377563 0.764705896377563;1 0.745098054409027 0.745098054409027;1 0.725490212440491 0.725490212440491;1 0.705882370471954 0.705882370471954;1 0.686274528503418 0.686274528503418;1 0.666666686534882 0.666666686534882;1 0.647058844566345 0.647058844566345;1 0.627451002597809 0.627451002597809;1 0.607843160629272 0.607843160629272;1 0.588235318660736 0.588235318660736;1 0.5686274766922 0.5686274766922;1 0.549019634723663 0.549019634723663;1 0.529411792755127 0.529411792755127;1 0.509803950786591 0.509803950786591;1 0.490196079015732 0.490196079015732;1 0.470588237047195 0.470588237047195;1 0.450980395078659 0.450980395078659;1 0.431372553110123 0.431372553110123;1 0.411764711141586 0.411764711141586;1 0.39215686917305 0.39215686917305;1 0.372549027204514 0.372549027204514;1 0.352941185235977 0.352941185235977;1 0.333333343267441 0.333333343267441;1 0.313725501298904 0.313725501298904;1 0.294117659330368 0.294117659330368;1 0.274509817361832 0.274509817361832;1 0.254901975393295 0.254901975393295;1 0.235294118523598 0.235294118523598;1 0.215686276555061 0.215686276555061;1 0.196078434586525 0.196078434586525;1 0.176470592617989 0.176470592617989;1 0.156862750649452 0.156862750649452;1 0.137254908680916 0.137254908680916;1 0.117647059261799 0.117647059261799;1 0.0980392172932625 0.0980392172932625;1 0.0784313753247261 0.0784313753247261;1 0.0588235296308994 0.0588235296308994;1 0.0392156876623631 0.0392156876623631;1 0.0196078438311815 0.0196078438311815;1 0 0];
	%  
	% set(gca,'CLim',[neg pos],'Colormap',...
	%     tmp,'FontSize',18);
	% colorbar;
	% grid off
	%% Reconstruction
	for i=1:length(i_plot)
    	semilogy(X(i_plot(i),:),(10^(i-1))*PASTP(i_plot(i),:).','r','LineWidth',2)
    	%semilogy(X(i_plot(i),:),(10^(i-1))*PA_Gauss(i_plot(i),:).','k--','LineWidth',2)
	end
	legend([p1,p2],{['$r/\lambda=$',num2str(samples(1),'%1.2f')],['$r/\lambda=$',num2str(samples(i_plot(length(i_plot))),'%1.2f')]},'Interpreter','latex','Location','south','FontSize',18)  
	legend([p1,p2],{['$r/\lambda=$',num2str(samples(1)/taylor_L_sample,'%1.2f')],['$r/\lambda=$',num2str(samples(i_plot(length(i_plot)))/taylor_L_sample,'%1.2f')]},'Interpreter','latex','Location','south','FontSize',18)  
	txt = {'(a)'};
	a=text(4,0.5,txt);
	a.Units='normalized';
	
	% semilogy(X(i,:)./rms_incr1(i),(10^-(i-1))*PA(i,:).','ko','LineWidth',2)
	
	subplot(1,2,2)
	for i=1:length(i_plot)
    	if i==1
        	p1  = semilogy(X(i_plot(i),:).*siginf./rms_incr1(i_plot(i)),(10^-(i-1))*PA(i_plot(i),:).','MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',c(i,:));%from experiment
    	elseif i==length(i_plot)
        	p2  = semilogy(X(i_plot(i),:).*siginf./rms_incr1(i_plot(i)),(10^-(i-1))*PA(i_plot(i),:).','MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',c(i,:));%from experiment 
    	else
        	semilogy(X(i_plot(i),:).*siginf./rms_incr1(i_plot(i)),(10^-(i-1))*PA(i_plot(i),:).','MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',c(i,:))%from experiment 
    	end
    	hold on;
    	set(gcf, 'Color', 'w')
    	axis square
    	if norm_ur==1  
    	xlabel('$u_r / rms(u_r)$', 'interpreter','latex')
    	else
        	xlabel('$u_r\ (m/s)$','interpreter','latex');
    	end
    	ylabel('$p(u_r / rms(u_r))$','interpreter','latex')
    	ylabel('PDF','interpreter','latex')
	%     xlim([-5 5])
    	fig_setup
	end
	%% Reconstruction
	for i=1:length(i_plot)
    	semilogy(X(i_plot(i),:).*siginf./rms_incr1(i_plot(i)),(10^-(i-1))*PASTP(i_plot(i),:).','r','LineWidth',2)
    	%semilogy(X(i_plot(i),:),(10^(i-1))*PA_Gauss(i_plot(i),:).','k--','LineWidth',2)
	end
	legend([p1,p2],{['$r/\lambda=$',num2str(samples(1),'%1.2f')],['$r/\lambda=$',num2str(samples(i_plot(length(i_plot))),'%1.2f')]},'Interpreter','latex','Location','south','FontSize',18)  
	legend([p1,p2],{['$r/\lambda=$',num2str(samples(1)/taylor_L_sample,'%1.2f')],['$r/\lambda=$',num2str(samples(i_plot(length(i_plot)))/taylor_L_sample,'%1.2f')]},'Interpreter','latex','Location','south','FontSize',18)  
	txt = {'(b)'};
	b=text(4,0.5,txt);
	b.Units='normalized';
	
	font=22;
	a.FontSize = font;
	b.FontSize = font;
	
	pos_txt=[-0.55   0.9];
	a.Position=pos_txt;
	b.Position=pos_txt;
    	
	
	% if ischar(save_path)
	%     savefig(h,fullfile(save_path,'Struc_recon.fig'),'compact')
	%     for a = 1:length(h)
	%         exportgraphics(h(a),fullfile(save_path,sprintf('Struc_recon_%d.png', a)))
	%     end
	% end
	% % set(gca,'ColorScale','log')
	% 
	% for i=1:2
	%     tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
	%     uiwait(gcf);
	%     delete(tmp_ui);
	%     close
	% end
	
	
	
	% figure
	% plot(samples)
	% 
	% figure
	% set(gcf, 'Color', 'w')
	% subplot(1,2,1)
	% plot(r__recon./taylor_L,S2_exp,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
	% hold on
	% plot(r__recon./taylor_L,S2_recon,'LineWidth',2,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 
	% set(gca,'xScale','log') 
	% set(gca,'yScale','log') 
	% set(gcf, 'Color', 'w')
	% % set(gca,'FontSize',fontsize)
	% xlabel('$r/\lambda$', 'interpreter','latex')
	% % xlabel('$r / \lambda$', 'interpreter','latex')
	% ylabel('$S^2(r)$', 'interpreter','latex')
	% % xlim([min(r./taylor_L)*0.5 max(r./taylor_L)*2])
	% % ylim([min(S_exp_2) max(S_exp_2)])
	% % vline(evaluated(1).r__recon/taylor_L,'r')
	% vline(taylor_L/taylor_L,'k')
	% vline(int_L/taylor_L,'k')
	% % set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0]);
	% axis square
	% fig_setup
	% legend off
	% 
	% subplot(1,2,2)
	% plot(r__recon./taylor_L,-S3_exp,'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
	% hold on
	% plot(r__recon./taylor_L,-S3_recon,'LineWidth',2,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]) 
	% set(gca,'xScale','log')
	% set(gca,'yScale','log') 
	% set(gcf, 'Color', 'w')
	% % set(gca,'FontSize',fontsize)
	% % xlabel('$r$', 'interpreter','latex')
	% xlabel('$r / \lambda$', 'interpreter','latex')
	% ylabel('$-S^3(r)$', 'interpreter','latex')
	% % xlim([min(r./taylor_L)*0.5 max(r./taylor_L)*2])
	% % ylim([min(-S_exp_3) max(-S_exp_3)])
	% % vline(evaluated(1).r__recon/taylor_L,'r')
	% vline(taylor_L/taylor_L,'k')
	% vline(int_L/taylor_L,'k')
	% % set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0]);
	% axis square
	% fig_setup
	% legend off
	
	if ischar(save_path)
    	savefig(h,fullfile(save_path,'KM_raw.fig'),'compact')
    	for a = 1:length(h)
        	exportgraphics(h(a),fullfile(save_path,sprintf('KM_raw_%d.png', a)))
    	end
	end
	
	for i=1:2
    	tmp_ui=uicontrol('Position',[10 10 100 40],'String','Continue','Callback','uiresume(gcbf)');
    	uiwait(gcf);
    	delete(tmp_ui);
    	close
	end

	close(bar)
end