%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the Kramers-Moyal coefficients $D^{(k)}\left(u_r,r\right)$ with $k={1,2,3,4}$
% for all scales (specified in \textbf{scale\_steps}) and for each bin (specified in \textbf{increment\_bin}) for 
% all values of longitudinal velocity increments by a linear extrapolation in $\Delta r$ of the 
% $k$-order conditional moments $M^{k}\left(u_r,r\right)$ (see Fig. \ref{fig:con_mom}) and the 
% function \textbf{KM\_plot\_raw} plots them accordingly. With $r'<r$:
%
% This limit approximation leads to uncertainties in the absolute values of the Kramers–Moyal 
% coefficients, whereas the functional forms of $D^{(k)}\left(u_r,r\right)$ are commonly well estimated. 
% In order to overcome this problem, the optimization algorithm described below is performed.
%
% Estimation of Kramers-Moyal coefficients (D1 and D2) for each scale and for each bin are given 
% by derivatives of the corresponding conditional moments % For each bin, a linear fit is computed 
% for all steps in steps.
%
% Arguments IN
% increment_bin = number of bins
% min_events = minimum number of events
% evaluated = struct array calculated in the function 'conditional_moment'
% step_con_moment = the steps at which the value of conditional moments are calculated
% Fs = Acquisition/Sampling Frequency in Hz
% taylor_L = Taylor length scale in meters
% m_data = mean of the data
% multi_point = Multipoint condition 1=YES or 2=NO
% condition = This input is for multipoint statistics
% norm_r = normalization of the scale using $\lambda$? data 1=Yes, 0=No
%
% Arguments OUT
% evaluated = a modified/updated struct 'evaluated' array
% Enclosed in a "evaluated":
% r is the scale in meters at which moments will be calculated and hence this r will
% be the same at which D1 & D2 will be calculated===>r2
% r_samp is nothing but the r in number of samples==>r2
% r_short_sample is nothing but r1 ===> (r2>r1) 
% D1 = Drift coefficient
% eD1 = error associated with drift coefficient
% D2 = Diffusion coefficient
% eD2 = error associated with diffusion coefficient
% D3 = Third order Kramers-Moyal coefficient
% D4 = Fourth order Kramers-Moyal coefficient
% D1_opti = Optimised D1
% D2_opti = Optimised D2
% M11 = First order conditional moment
% M21 = Second order conditional moment
% M31 = Third order conditional moment
% M41 = Fourth order conditional moment
% eM1 = Error associated with M11
% eM2 = Error associated with M21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [evaluated]=KM_Calculation(increment_bin,min_events,evaluated,step_con_moment,Fs,taylor_L,m_data,multi_point,condition,norm_r)
%------------------------------- D1 und D2 --------------------------------
%--------------------------------------------------------------------------
% samples of taylor length for normalization r/lambda, or steps/step_lambda
if norm_r==1 
step_taylor_L   = (taylor_L/(m_data/Fs));
% step_taylor_L = ceil(Fs*taylor_L/m_data);  
else
step_taylor_L   = 1;
end

if multi_point==1
        for k=1:length(condition)
            my_field = sprintf('point_eval_%d', k);

            % This loop is for each scale at which KMC's will be calculated
            for scal = 1:size(evaluated.(my_field),2) 
                % Initialising all the arrays/variables
                D1(1:increment_bin)  = nan;
                eD1(1:increment_bin) = nan;
                D2(1:increment_bin)  = nan;
                D2b(1:increment_bin) = nan; 
                eD2(1:increment_bin) = nan; 
                eD2b(1:increment_bin)= nan;
                D3(1:increment_bin)  = nan; 
                D4(1:increment_bin)  = nan;

            %     if length(evaluated.(my_field)(scal).x_mean_bin)==aufl & length(evaluated.(my_field)(scal).y_mean_bin)==aufl
                    for i=1:increment_bin %%This loop is for each bin of a specific scale at which KMC's will be calculated
                        xData1=nan;
                        yData1=nan;
                        xData2=nan;
                        yData2=nan;
                        xData3=nan;
                        yData3=nan;
                        xData4=nan;
                        yData4=nan;
            %             if  min(evaluated.(my_field)(scal).counter_B)>=min_events
                            %             if  ~any(imag(evaluated.(my_field)(scal).eM1(i,:))) && ~any(imag(evaluated.(my_field)(scal).eM2(i,:)))
                            if  evaluated.(my_field)(scal).eM1(i,:)>0 & evaluated.(my_field)(scal).eM2(i,:)>0

                                % D1
            %% SLOPE METHOD
                                % [xData1, yData1, weights1] = prepareCurveData(step_con_moment,evaluated.(my_field)(scal).M11(i,:),(1./evaluated.(my_field)(scal).eM1(i,:)));
                                % ft1 = fittype( 'poly1' );
                                % opts = fitoptions( 'Method', 'LinearLeastSquares' );
                                % opts.Weights = weights1;
                                % [fitresult1, gof1] = fit( xData1, yData1, ft1, opts );
                                % plot(fitresult1, xData1, yData1)
                                % D1(i)=fitresult1.p1.*Fs;
            %% INTERCEPT METHOD
                                % Prepare data to plot the curve of (r/lambda) Vs. M11
                                % So we will look at the value of M11 at which (r/lambda)
                                % tends to zero & that would be our D1
                                % Here M11 and (r/lambda) are dimentionless and hence
                                % D1 will also be dimentionless

                                [xData1, yData1, weights1] = prepareCurveData((step_con_moment./step_taylor_L),(evaluated.(my_field)(scal).M11(i,:)./(step_con_moment./step_taylor_L)),((1./evaluated.(my_field)(scal).eM1(i,:))./(step_con_moment./step_taylor_L)));
                                ft1 = fittype( 'poly1' );
                                opts = fitoptions( 'Method', 'LinearLeastSquares' );
                                opts.Weights = weights1;
                                [fitresult1, gof1] = fit( xData1, yData1, ft1, opts );
                                %                 plot(fitresult1, xData1, yData1)
                                D1(i)=fitresult1.p2; %%This is the intercept/offset on y-axis i.e. M11

                                %--------------------------------------------------------------------------- 
                                %--------------------------------------------------------------------------- 
                                % D2 with correction METHOD
                                [xData2, yData2, weights2] = prepareCurveData((step_con_moment./step_taylor_L),...
                                    (((evaluated.(my_field)(scal).M21(i,:)-((step_con_moment./step_taylor_L).*fitresult1.p2).^2)./2)./(step_con_moment./step_taylor_L)),...
                                    ((1./evaluated.(my_field)(scal).eM2(i,:))./(step_con_moment./step_taylor_L)));
                                %[xData2, yData2, weights2] = prepareCurveData(step_con_moment,(evaluated.(my_field)(scal).M21(i,:)-((fitresult1.p1.*step_con_moment).^2)),(1./evaluated.(my_field)(scal).eM2(i,:)));
                                ft2 = fittype( 'poly1' );
                                opts = fitoptions( 'Method', 'LinearLeastSquares' );
                                opts.Weights = weights2;
                                [fitresult2, gof2] = fit( xData2, yData2, ft2, opts );
                                %                 D2(i)=fitresult2.p1.*Fs/2;
                                %Slope
            %                     D2(i)=fitresult2.p1.*Fs;
                                %Intercept
                                D2(i)=fitresult2.p2;%%This is the intercept
                                %--------------------------------------------------------------------------- 
                                %---------------------------------------------------------------------------
                                % D2 without correction
            %% Slope
                                % [xData2b, yData2b, weights2b] = prepareCurveData(step_con_moment,(evaluated.(my_field)(scal).M21(i,:)./2),(1./evaluated.(my_field)(scal).eM2(i,:)));
                                % ft2b = fittype( 'poly1' );
                                % opts = fitoptions( 'Method', 'LinearLeastSquares' );
                                % opts.Weights = weights2b;
                                % [fitresult2b, gof2b] = fit( xData2b, yData2b, ft2b, opts );
                                % D2b(i)=fitresult2b.p1.*Fs;
            %% INTERCEPT METHOD WITHOUT CORRECTION TO D2
                                % Prepare data to plot the curve of (r/lambda) Vs. M21
                                % So we will look at the value of M21 at which (r/lambda)
                                % tends to zero & that would be our D2
                                % Here M21 and (r/lambda) are dimentionless and hence
                                % D2 will also be dimentionless
                                [xData2b, yData2b, weights2b] = prepareCurveData((step_con_moment./step_taylor_L),((evaluated.(my_field)(scal).M21(i,:)./2)./(step_con_moment./step_taylor_L)),((1./evaluated.(my_field)(scal).eM2(i,:))./(step_con_moment./step_taylor_L)));
                                ft2b = fittype( 'poly1' );
                                opts = fitoptions( 'Method', 'LinearLeastSquares' );
                                opts.Weights = weights2b;
                                [fitresult2b, gof2b] = fit( xData2b, yData2b, ft2b, opts );
                                D2b(i)=fitresult2b.p2;


                                %--------------------------------------------------------------------------- 
                                %---------------------------------------------------------------------------
                                % D3 with correction METHOD
                                [xData3, yData3] = prepareCurveData((step_con_moment./step_taylor_L),(evaluated.(my_field)(scal).M31(i,:)-(6.*(step_con_moment./step_taylor_L).^2.*fitresult1.p1.*fitresult2.p1))./6);
                                ft3 = fittype( 'poly1' );
                                opts = fitoptions( 'Method', 'LinearLeastSquares' );
                                [fitresult3, gof3] = fit( xData3, yData3, ft3, opts );
            %% Slope
                                % D3(i)=fitresult3.p1.*Fs; 
                                %Intercept
                                D3(i)=fitresult3.p2;

            %% INTERCEPT METHOD WITHOUT CORRECTION TO D3
                                [xData3b, yData3b] = prepareCurveData((step_con_moment./step_taylor_L),((evaluated.(my_field)(scal).M31(i,:)./6)./(step_con_moment./step_taylor_L)));
                                ft3b = fittype( 'poly1' );
                                opts = fitoptions( 'Method', 'LinearLeastSquares' );
                                [fitresult3b, gof3b] = fit( xData3b, yData3b, ft3b, opts );
                                D3b(i)=fitresult3b.p2;


                                %--------------------------------------------------------------------------- 
                                %---------------------------------------------------------------------------
                                % D4 with correction METHOD
                                [xData4, yData4] = prepareCurveData((step_con_moment./step_taylor_L),(evaluated.(my_field)(scal).M41(i,:)-(12.*(step_con_moment./step_taylor_L).^2.*(2.*fitresult1.p1.*fitresult3.p1+fitresult2.p1.^2)))./24);  
                                ft4 = fittype( 'poly1' );
                                opts = fitoptions( 'Method', 'LinearLeastSquares' );
                                [fitresult4, gof4] = fit( xData4, yData4, ft4, opts );
            %% Slope
                                % D4(i)=fitresult4.p1.*Fs;
            %% Intercept
                                D4(i)=fitresult4.p2;

            %% INTERCEPT METHOD WITHOUT CORRECTION TO D4
                                [xData4b, yData4b] = prepareCurveData((step_con_moment./step_taylor_L),((evaluated.(my_field)(scal).M41(i,:)./24)./(step_con_moment./step_taylor_L)));
                                ft4b = fittype( 'poly1' );
                                opts = fitoptions( 'Method', 'LinearLeastSquares' );
                                [fitresult4b, gof4b] = fit( xData4b, yData4b, ft4b, opts );
                                D4b(i)=fitresult4b.p2;


                                %-----------------------------------------------
                                %Error calc for D1 and D2b
                                eD1(i)  = real(sqrt(abs((2.*D2b(i)-D1(i).^2)/max(evaluated.(my_field)(scal).counter_B))));
                                %eD2(i)  = real(sqrt(abs((2.*D4(i)-D2(i).^2)/max(evaluated.(my_field)(scal).counter_B))));
                                eD2b(i) = real(sqrt(abs((2.*D4b(i)-D2b(i).^2)/max(evaluated.(my_field)(scal).counter_B))));
                            end
                    end

                D2          = D2+abs(min(D2));
                D2(D2<0)    = 0;
                D2          = abs(D2);   
                D2b(D2b<0)  = 0;

                evaluated.(my_field)(scal).x_bin_not_nan   = (evaluated.(my_field)(scal).counter_B>=min_events & ~isnan(D1));

                evaluated.(my_field)(scal).D1(1,:)         = D1;
                evaluated.(my_field)(scal).eD1(1,:)        = eD1;

                evaluated.(my_field)(scal).D2(1,:)         = abs(D2b);
                evaluated.(my_field)(scal).eD2(1,:)        = eD2b;

                evaluated.(my_field)(scal).D3(1,:)         = abs(D3);
                evaluated.(my_field)(scal).D4(1,:)         = abs(D4);
                
                disp(['point:' num2str(k) '/'  num2str(length(condition)) ])
                disp(['scal:' num2str(scal) '/'  num2str(size(evaluated.(my_field),2)) ])
            end
            
            
        end
else
    % This loop is for each scale at which KMC's will be calculated
    f = uifigure;
    d = uiprogressdlg(f,'Title','Please Wait','ShowPercentage','on');
    for scal = 1:size(evaluated,2) 
        % Initialising all the arrays/variables
        D1(1:increment_bin)  = nan;
        eD1(1:increment_bin) = nan;
        D2(1:increment_bin)  = nan;
        D2b(1:increment_bin) = nan; 
        eD2(1:increment_bin) = nan; 
        eD2b(1:increment_bin)= nan;
        D3(1:increment_bin)  = nan; 
        D4(1:increment_bin)  = nan;

    %     if length(evaluated(scal).x_mean_bin)==aufl & length(evaluated(scal).y_mean_bin)==aufl
            for i=1:increment_bin %%This loop is for each bin of a specific scale at which KMC's will be calculated
                xData1=nan;
                yData1=nan;
                xData2=nan;
                yData2=nan;
                xData3=nan;
                yData3=nan;
                xData4=nan;
                yData4=nan;
    %             if  min(evaluated(scal).counter_B)>=min_events
                    %             if  ~any(imag(evaluated(scal).eM1(i,:))) && ~any(imag(evaluated(scal).eM2(i,:)))
                    if  evaluated(scal).eM1(i,:)>0 & evaluated(scal).eM2(i,:)>0

                        % D1
        %% SLOPE METHOD
                        % [xData1, yData1, weights1] = prepareCurveData(step_con_moment,evaluated(scal).M11(i,:),(1./evaluated(scal).eM1(i,:)));
                        % ft1 = fittype( 'poly1' );
                        % opts = fitoptions( 'Method', 'LinearLeastSquares' );
                        % opts.Weights = weights1;
                        % [fitresult1, gof1] = fit( xData1, yData1, ft1, opts );
                        % plot(fitresult1, xData1, yData1)
                        % D1(i)=fitresult1.p1.*Fs;
        %% INTERCEPT METHOD
                        % Prepare data to plot the curve of (r/lambda) Vs. M11
                        % So we will look at the value of M11 at which (r/lambda)
                        % tends to zero & that would be our D1
                        % Here M11 and (r/lambda) are dimentionless and hence
                        % D1 will also be dimentionless
                        [xData1, yData1, weights1] = prepareCurveData((step_con_moment./step_taylor_L),(evaluated(scal).M11(i,:)./(step_con_moment./step_taylor_L)),((1./evaluated(scal).eM1(i,:))./(step_con_moment./step_taylor_L)));
                        ft1 = fittype( 'poly1' );
                        opts = fitoptions( 'Method', 'LinearLeastSquares' );
                        opts.Weights = weights1;
                        [fitresult1, gof1] = fit( xData1, yData1, ft1, opts );
                        %                 plot(fitresult1, xData1, yData1)
                        D1(i)=fitresult1.p2; %%This is the intercept/offset on y-axis i.e. M11

                        %--------------------------------------------------------------------------- 
                        %--------------------------------------------------------------------------- 
                        % D2 with correction METHOD
                        [xData2, yData2, weights2] = prepareCurveData((step_con_moment./step_taylor_L),...
                            (((evaluated(scal).M21(i,:)-((step_con_moment./step_taylor_L).*fitresult1.p2).^2)./2)./(step_con_moment./step_taylor_L)),...
                            ((1./evaluated(scal).eM2(i,:))./(step_con_moment./step_taylor_L)));
                        %[xData2, yData2, weights2] = prepareCurveData(step_con_moment,(evaluated(scal).M21(i,:)-((fitresult1.p1.*step_con_moment).^2)),(1./evaluated(scal).eM2(i,:)));
                        ft2 = fittype( 'poly1' );
                        opts = fitoptions( 'Method', 'LinearLeastSquares' );
                        opts.Weights = weights2;
                        [fitresult2, gof2] = fit( xData2, yData2, ft2, opts );
                        %                 D2(i)=fitresult2.p1.*Fs/2;
        %% Slope
                        % D2(i)=fitresult2.p1.*Fs;
        %% Intercept
                        D2(i)=fitresult2.p2;
                        %--------------------------------------------------------------------------- 
                        %---------------------------------------------------------------------------
                        % D2 without correction
        %% Slope
                        % [xData2b, yData2b, weights2b] = prepareCurveData(step_con_moment,(evaluated(scal).M21(i,:)./2),(1./evaluated(scal).eM2(i,:)));
                        % ft2b = fittype( 'poly1' );
                        % opts = fitoptions( 'Method', 'LinearLeastSquares' );
                        % opts.Weights = weights2b;
                        % [fitresult2b, gof2b] = fit( xData2b, yData2b, ft2b, opts );
                        % D2b(i)=fitresult2b.p1.*Fs;
        %% INTERCEPT METHOD WITHOUT CORRECTION TO D2
                        % Prepare data to plot the curve of (r/lambda) Vs. M21
                        % So we will look at the value of M21 at which (r/lambda)
                        % tends to zero & that would be our D2
                        % Here M21 and (r/lambda) are dimentionless and hence
                        % D2 will also be dimentionless
                        [xData2b, yData2b, weights2b] = prepareCurveData((step_con_moment./step_taylor_L),((evaluated(scal).M21(i,:)./2)./(step_con_moment./step_taylor_L)),((1./evaluated(scal).eM2(i,:))./(step_con_moment./step_taylor_L)));
                        ft2b = fittype( 'poly1' );
                        opts = fitoptions( 'Method', 'LinearLeastSquares' );
                        opts.Weights = weights2b;
                        [fitresult2b, gof2b] = fit( xData2b, yData2b, ft2b, opts );
                        D2b(i)=fitresult2b.p2;


                        %--------------------------------------------------------------------------- 
                        %---------------------------------------------------------------------------
                        % D3 with correction METHOD
                        [xData3, yData3] = prepareCurveData((step_con_moment./step_taylor_L),(evaluated(scal).M31(i,:)-(6.*(step_con_moment./step_taylor_L).^2.*fitresult1.p1.*fitresult2.p1))./6);
                        ft3 = fittype( 'poly1' );
                        opts = fitoptions( 'Method', 'LinearLeastSquares' );
                        [fitresult3, gof3] = fit( xData3, yData3, ft3, opts );
        %% Slope
                        % D3(i)=fitresult3.p1.*Fs; 
        %% Intercept
                        D3(i)=fitresult3.p2;

        %% INTERCEPT METHOD WITHOUT CORRECTION TO D3
                        [xData3b, yData3b] = prepareCurveData((step_con_moment./step_taylor_L),((evaluated(scal).M31(i,:)./6)./(step_con_moment./step_taylor_L)));
                        ft3b = fittype( 'poly1' );
                        opts = fitoptions( 'Method', 'LinearLeastSquares' );
                        [fitresult3b, gof3b] = fit( xData3b, yData3b, ft3b, opts );
                        D3b(i)=fitresult3b.p2;


                        %--------------------------------------------------------------------------- 
                        %---------------------------------------------------------------------------
                        % D4 with correction METHOD
                        [xData4, yData4] = prepareCurveData((step_con_moment./step_taylor_L),(evaluated(scal).M41(i,:)-(12.*(step_con_moment./step_taylor_L).^2.*(2.*fitresult1.p1.*fitresult3.p1+fitresult2.p1.^2)))./24);  
                        ft4 = fittype( 'poly1' );
                        opts = fitoptions( 'Method', 'LinearLeastSquares' );
                        [fitresult4, gof4] = fit( xData4, yData4, ft4, opts );
        %% Slope
                        % D4(i)=fitresult4.p1.*Fs;
        %% Intercept
                        D4(i)=fitresult4.p2;

        %% INTERCEPT METHOD WITHOUT CORRECTION TO D4
                        [xData4b, yData4b] = prepareCurveData((step_con_moment./step_taylor_L),((evaluated(scal).M41(i,:)./24)./(step_con_moment./step_taylor_L)));
                        ft4b = fittype( 'poly1' );
                        opts = fitoptions( 'Method', 'LinearLeastSquares' );
                        [fitresult4b, gof4b] = fit( xData4b, yData4b, ft4b, opts );
                        D4b(i)=fitresult4b.p2;


                        %-----------------------------------------------
                        %Error calc for D1 and D2b
                        eD1(i)  = real(sqrt(abs((2.*D2b(i)-D1(i).^2)/max(evaluated(scal).counter_B))));
                        % eD2(i)  = real(sqrt(abs((2.*D4(i)-D2(i).^2)/max(evaluated(scal).counter_B))));
                        eD2b(i) = real(sqrt(abs((2.*D4b(i)-D2b(i).^2)/max(evaluated(scal).counter_B))));
                    end
            end

        D2          = D2+abs(min(D2));
        D2(D2<0)    = 0;
        D2          = abs(D2);   
        D2b(D2b<0)  = 0;

        evaluated(scal).x_bin_not_nan   = (evaluated(scal).counter_B>=min_events & ~isnan(D1));

        evaluated(scal).D1(1,:)         = D1;
        evaluated(scal).eD1(1,:)        = eD1;

        evaluated(scal).D2(1,:)         = abs(D2b);
        evaluated(scal).eD2(1,:)        = eD2b;

        evaluated(scal).D3(1,:)         = abs(D3);
        evaluated(scal).D4(1,:)         = abs(D4);

        disp(['scal:' num2str(scal) '/'  num2str(size(evaluated,2)) ])
        d.Value=scal/size(evaluated,2);
    end
close(d);
close(f);
end
end