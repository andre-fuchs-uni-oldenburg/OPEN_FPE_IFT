
import numpy as np
from scipy.stats import kurtosis, skew
import matplotlib.pyplot as plt

def chunks(arrl, n):
    sarr = []
    for i in range(0, len(arrl), n):
        sarr.append(  arrl[i:i+n] )
    return sarr;

def matrange(arr):
    return np.max(arr)-np.min(arr) ;


####################################################################################################
####################################################################################################
# This function plots the mean, standard deviation, skewness and kurtosis of sections of a length of
# 5\# of the data to check the stationarity of the data. In the title of the figure the number of
# nan's in the dataset and the turbulence intensity is printed.
#
# Arguments IN
# data = 1D array of data
#
# Arguments OUT
# data = 1D array of data without nan's
####################################################################################################
####################################################################################################

def plot_stationarity(data,outpath,percent=5):

    #percent = 5
    k_split = int(100 / percent)
    
    if data.size < k_split:
        k_split = data.size

    #tmp = np.array_split(data, k_split)   
    tmp = chunks(data, k_split)   

    if tmp[0].shape != tmp[-1].shape :
        tmp = tmp[:-1]
    tmp = np.array(tmp)

    x = np.linspace(0,100,tmp.shape[1])

    plt.clf()
    kwargs = { 'linewidth': 2,'markersize': 8,'marker':'o','linestyle':'none' , 'fillstyle' : 'none' } 
    plt.plot(x,np.nanmean(tmp,axis=0), label=r'mean $(m/s)$', **kwargs )
    plt.plot(x,np.nanstd(tmp,axis=0), label=r'std $(m/s)$', **kwargs )
    plt.plot(x,skew(tmp,axis=0),  label='skewnss', **kwargs )
    plt.plot(x,kurtosis(tmp,axis=0),  label='kurtosis', **kwargs )

    plt.xlim(np.array([-5,105]))
    plt.xlabel(r'Percentage of data')
    #plt.ylabel('m/s')
    plt.legend()#location='east')
    plt.title('Number of nans = '+str( np.count_nonzero(np.isnan(data)) ))
    plt.savefig(outpath+'plot_precent_vs_datastats.png')
    plt.show()
    print('>> Plot saved to '+outpath+'plot_precent_vs_datastats.png')
    


    plt.clf()
    plt.plot(np.linspace(0,100,data.size),data,'k')
    plt.xlim(np.array([-5,105]))
    plt.xlabel(r'Percentage of data')
    plt.ylabel( r'$u\ (m/s)$')
    plt.title('Ti = '+str( np.round(np.nanstd(data)/np.nanmean(data)*100, 2) )+'%' )
    plt.savefig(outpath+'plot_precent_vs_data.png')
    plt.show()
    print('>> Plot saved to '+outpath+'plot_precent_vs_data.png')
    
    return ;



####################################################################################################
####################################################################################################
# This function plots the probability density function (PDF) of the data with the specified number
# of bins. It also plots the Gaussian distribution which has the same standard deviation and mean as
# of the data. In the title of the figure the range of the data (difference between the maximum and
# minimum values of sample data), the skewness and flatness of the data is printed.
# COMMENT: hist is a standard matlab function==>Look in to standard MATLAB documentation
#
# Arguments IN
# data  = 1D array for which you would like to plot the probability density function(PDF)
# nbins = The number of bins in which you would like to divide your data
####################################################################################################
####################################################################################################

from scipy.stats import norm

def plot_pdf(data,nbins,outpath):
    fuL, uL = np.histogram(data,nbins, density=True)  # fuL is the counters & uL is the center of the bin

    #Converting fuL to puL in order to change counters to probability density fuction
    #Probability = fuL./nansum(fuL)==> Sum of Probability array should be 1.0<== Definition of Probability 
    #puL = [Probability./mean(diff(uL)] ==> Probability density function
    #==>trapz(uL(1:end),puL) must be 1.0 <== Definition of Probability density function
    #puL = fuL / ( np.nansum(fuL) * np.mean(np.diff(uL)) ) #density=True sets Probability density function

    uL = (uL[1:]+uL[:-1])*0.5 #bin_centers


    #'MarkerEdgeColor','k','MarkerSize',6,'Marker','o','LineStyle','none')
    plt.clf()
    # Actual PDF of the data
    plt.semilogy(uL,fuL, label='Original data', marker='o', markersize=6, markerfacecolor=None, markeredgecolor='black', linestyle='none' , fillstyle='none' )

    xPos = np.mean(data)
    yPos = np.max(fuL)
    # Plot a vertical line at x=mean(data) to visualize a flatness or skewness
    plt.plot([xPos,xPos], [0,yPos] , label='Mean of data' , linestyle='--',linewidth=2, color ='black')

    # Plot a gaussian distribution PDF with same mean and same standard deviation as of our original data for comparison.
    mu, sigma = norm.fit(data)
    gaus_bin = np.linspace(np.amin(data),np.amax(data),nbins)  
    pdf = norm.pdf(gaus_bin,loc=mu,scale=sigma)
    plt.plot(gaus_bin, pdf, label='Gaussian Distribution' , linestyle='--',linewidth=2, color ='grey') # [( 0.501960813999176,0.501960813999176,0.501960813999176 )])

    plt.legend()
    plt.xlabel(r'$u\ (m/s)$')
    plt.ylabel(r'PDF')

    plttitle = 'range = '+str( np.round(matrange(data),2)) + ' S = '+str( np.round(skew(data),2))+' K = '+str( np.round(kurtosis(data),2)) 
    plt.title(plttitle )

    plt.savefig(outpath+'plot_data_pdf.png')
    plt.show()
    print('>> Plot saved to '+outpath+'plot_data_pdf.png')

    return ;




def plotSpectrum(data,outpath,fsamp,increment_bin):
    # PSD of the fluctuations; mean value does not play a role in spectrum
    data    = data - np.nanmean(data)

    # This is to make length of data array even
    if data.size%2 != 0 : 
        data = data[:-1]

    """
    L = data.size
    #% spek              = abs(fft(data,L)).^2/L;        % Power spectral density(PSD)
    spek              = np.square( np.abs(np.fft(data,L)) )/ (L*fsamp)    # % Energy spectral density(ESD) using a fft
    spek              = 2.*spek[2:L/2+1]               #% FFt will yield half number of unique points
    f                 = fsamp/2*np.linspace(0,1,L/2+1) #     % Nyquist frequency of the signal==Fs/2
    f                 = f[2:]                       #% Remove zero Hz component 


    plot_length         = increment_bin*10
    if plot_length%2 == 0
        plot_length = plot_length +1 

    intervall=unique(round(logspace(0,log10(L/2),plot_length)),'stable');

    # moving average with equally spaced frequency interval in lin-space 
    # intervall=unique(round(linspace(1,L/2,plot_length)),'stable');

    plot_length= len(intervall) #; % Number of points for plotting averaged spectrum

    #% Initializing the arrays/preallocation
    spek_smooth = np.zeros(plot_length-2,1);
    x_achse_f_spektrum = np.zeros(plot_length-2,1);


    # Averaging of spectrum and hence smoothing
    for i=2:(plot_length-1)
        x_achse_f_spektrum(i-1,1) = mean(f(intervall(i-1):intervall(i+1)));
        spek_smooth(i-1,1)        = mean(spek(intervall(i-1):intervall(i+1)));
    end
    spek_power = spek_smooth.*Fs;
    """

    import scipy.signal

    # f contains the frequency components
    # S is the PSD
    # between computing the power spectral density (‘density’) where Pxx has units of V**2/Hz and computing the power spectrum (‘spectrum’) where Pxx has units of V**2, if x is measured in V and fs is measured in Hz. Defaults to ‘density’

    #(f, S) = scipy.signal.periodogram(data, fsamp, scaling='density') #PSD
    #(f, S) = scipy.signal.periodogram(data, fsamp, scaling='spectrum') #ESD


    #https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.welch.html
    def getSpect(data, fsamp,segmentLength=None):
        (f, S) = scipy.signal.welch(data, fsamp, window='hann', nperseg=segmentLength, noverlap=None, nfft=None, detrend='constant', return_onesided=True, scaling='spectrum', average='mean') 
        return f,S

    segmentLength = 1024*4#data.size//40 ### CHECK

    plt.clf()
    #Plot of Frequency Vs. Energy spectral density(ESD) without smoothing
    #f_plot,spek_plot = getSpect(data, fsamp,segmentLength)

    f_plot,spek_plot = scipy.signal.periodogram(data, fsamp, window='hann', nfft=None, detrend='constant', return_onesided=True, scaling='spectrum')
    plt.loglog(f_plot,spek_plot, 'r-', label='raw') #semilogy

    #Plot of Frequency Vs. Energy spectral density(ESD) with smoothing
    data_chunks = np.array_split(data,40) #get n=40 chunks
    f_plot_smooth, spek_smooth = [] ,[]
    print("Calculating averaged spectrum (may take a while)..")
    for chnk in data_chunks:
        fp, sp = getSpect(chnk, fsamp,segmentLength)
        f_plot_smooth.append(fp) 
        spek_smooth.append(sp) 
    f_plot_smooth, spek_smooth = np.array( f_plot_smooth).mean(axis=0) , np.array( spek_smooth).mean(axis=0)
    plt.loglog(f_plot_smooth,spek_smooth,'b-', label='averaged')

    plt.xlabel(r'$f\ (Hz)$' )
    plt.ylabel(r'$E\ (m^2/s)$')
    plt.title('Energy density spectrum' )
    #plt.xlim(min(f_plot_smooth)*2, max(f_plot_smooth)*2)
    plt.title('Energy density spectrum')

    plt.legend()
    plt.savefig(outpath+'plot_data_spectrum.png')
    print('>> Plot saved to '+outpath+'plot_data_spectrum.png')

    #plt.show()

    return ;






