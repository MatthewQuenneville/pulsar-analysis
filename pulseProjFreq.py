#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pylab as plt
import pulseFinder as pf
from math import factorial
import warnings
import pulseSpec as ps
from scipy.special import erfc
from scipy.optimize import curve_fit

# Resolution to use for searching in seconds. Must be larger than or
# equal to phase bin size.
searchRes=1.0/10000

pulseWidth=0.0001 # 100 microseconds

leadWidth=0.0001
trailWidth=0.0003

def expModGauss(x,sigma):
    mu=0.0
    gamma=1.0
    return (gamma/2.)*np.exp((gamma/2)*(2*mu+gamma*sigma*sigma-2*x))*erfc((mu+gamma*sigma*sigma-x)/(np.sqrt(2)*sigma))

if __name__ == "__main__":
    # Load files
    w,binWidth,_=pf.loadFiles(sys.argv[1:])

    # Get run information
    deltat=pf.getDeltaT(sys.argv[1])
    telescope=pf.getTelescope(sys.argv[1])
    startTime=pf.getStartTime(sys.argv[1])
    freqBand=pf.getFrequencyBand(telescope)

    # Rebin to find giant pulses, then resolve pulses with finer binning
    nSearchBins=min(w.shape[1],int(round(deltat/searchRes)))
    
    w_rebin=pf.rebin(w,nSearchBins)
    timeSeries_rebin=pf.getTimeSeries(w_rebin)
    timeSeries=pf.getTimeSeries(w)
    pulseList=pf.getPulses(timeSeries_rebin,binWidth=searchRes)
    if nSearchBins<w.shape[1]:
        pulseList=[(pf.resolvePulse(
                    timeSeries,int(pos*w.shape[1]/nSearchBins),
                    binWidth=binWidth,searchRadius=1.0/10000),height) 
                   for (pos,height) in pulseList]
    try:
        largestPulse=pulseList[0][0]
    except IndexError:
        print "Error, no giant pulse found in "+telescope+" for start time:"
        print startTime.iso
        sys.exit()
    
    # Find range of pulse to plot
    trailBins=int(np.ceil(trailWidth/binWidth))       
    leadBins=int(np.ceil(leadWidth/binWidth))
    
    nPulseBins=int(np.ceil(pulseWidth/binWidth))
    pulseRange=range(largestPulse-leadBins,largestPulse+trailBins)
    offRange=range(largestPulse-leadBins-(trailBins+leadBins),
                   largestPulse-leadBins)

    # Add entries to dynamic spectra and frequency band dictionaries
    dynamicSpec=ps.dynSpec(w,indices=pulseRange,normChan=False)
    dynamicSpec_BG=ps.dynSpec(w,indices=offRange,normChan=False)

    if w.shape[-1]==4:
        profile=dynamicSpec[:,:,(0,3)].sum(0).sum(-1)
    else:
        profile=dynamicSpec.sum(0)
    pulseBins=sorted(range(len(pulseRange)),key=lambda x: profile[x],
                      reverse=True)[:nPulseBins]

    # Plot spectra    
    spec=dynamicSpec[:,pulseBins,:].sum(1)-dynamicSpec_BG.mean(1)*nPulseBins

    freqList=[freqBand[0]+i*(freqBand[1]-freqBand[0])/spec.shape[0] 
              for i in range(spec.shape[0])]
    
    # Plot spectra of each polarization
    if spec.shape[-1]==4:
        f,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex='col',sharey='row')
        ax1.plot(freqList,spec[:,0],'.')
        ax1.set_ylabel('Intensity (Tsys)')
        ax1.set_title('Polarization 0')

        ax2.plot(freqList,spec[:,3],'.')
        ax2.set_title('Polarization 3')

        ax3.plot(freqList,spec[:,1],'.')
        ax3.set_xlabel('Frequency (MHz)')
        ax3.set_ylabel('Intensity (Tsys)')
        ax3.set_title('Polarization 1') 

        ax4.plot(freqList,spec[:,2],'.')
        ax4.set_xlabel('Frequency (MHz)')
        ax4.set_title('Polarization 2')

        plt.suptitle('Spectra',size=16)
        plt.show()
        
    else:
        plt.plot(freqList,spec,'.')
        plt.title('Spectrum')
        plt.ylabel('Intensity (Tsys)')
        plt.xlabel('Frequency Channel')
        plt.show()
        
    # Plot histogram of spectral noise
    if spec.shape[-1]==4:
        # Normalize intensity
        spec1=spec[:,0]/np.mean(spec[:,0])
        spec2=spec[:,3]/np.mean(spec[:,3])
                               
        # Plot histograms and plot if lmfit is found
        xmin=min(min(spec2),min(spec1))
        xmax=max(max(spec2),max(spec1))
        
        bins=np.linspace(np.floor(xmin),np.ceil(xmax),50)
        x_fine=np.linspace(np.floor(xmin),np.ceil(xmax),1000)

        f,(ax1,ax2) = plt.subplots(1,2,sharex='col',sharey='row')
        specHist1=ax1.hist(spec1,normed=True,bins=bins,label='Data')
        binCenters=np.array([(bins[i+1]+bins[i])/2. for i in range(len(bins)-1)
                             if not(specHist1[0][i]==0)])
        binEntries=np.array([i for i in specHist1[0] if not i==0])

        # Plot exponentially modified gaussian with parameters
        # estimated via fitting, and direct parameter estimation
        popt,pcov=curve_fit(expModGauss,binCenters,binEntries,
                            sigma=np.power(binEntries,0.5))
        if np.std(spec1)>1.0:
            ax1.plot(x_fine,expModGauss(x_fine,np.sqrt(np.var(spec1)-1.0)),label='Estimated')
        if popt[0]>0.0:
            ax1.plot(x_fine,expModGauss(x_fine,popt[0]),label='Fitted')
        
        ax1.legend()
        ax1.set_xlim(min(bins),max(bins))
        ax1.set_xlabel("Intensity")
        ax1.set_title('Polarization 0')
        specHist2=ax2.hist(spec2,normed=True,bins=bins,label='Data')
        binCenters=np.array([(bins[i+1]+bins[i])/2. for i in range(len(bins)-1)
                             if not(specHist2[0][i]==0)])
        binEntries=np.array([i for i in specHist2[0] if not i==0])
        # Plot exponentially modified gaussian with parameters
        # estimated via fitting, and direct parameter estimation
        popt,pcov=curve_fit(expModGauss,binCenters,binEntries,
                            sigma=np.power(binEntries,0.5))
        if np.std(spec2)>1.0:
            ax2.plot(x_fine,expModGauss(x_fine,np.sqrt(np.var(spec2)-1.0)),label='Estimated')
        if popt[0]>0.0:
            ax2.plot(x_fine,expModGauss(x_fine,popt[0]),label='Fitted')
        
        ax2.legend()
        ax2.set_title('Polarization 3')
        ax2.set_xlim(min(bins),max(bins))
        ax2.set_xlabel("Intensity")
        plt.show()
    else:
        # Normalize intensity
        specNorm=spec/np.mean(spec)

        # Plot histograms and fit if lmfit is found
        xmin=min(specNorm)
        xmax=max(specNorm)
        
        bins=np.linspace(np.floor(xmin),np.ceil(xmax),50)
        x_fine=np.linspace(np.floor(xmin),np.ceil(xmax),1000)
        specHist=plt.hist(specNorm,normed=True,bins=bins,label='Data')
        binCenters=np.array([(bins[i+1]+bins[i])/2. for i in range(len(bins)-1) if not(specHist[0][i]==0)])
        binEntries=np.array([i for i in specHist[0] if not i==0])
        # Plot exponentially modified gaussian with parameters
        # estimated via fitting, and direct parameter estimation
        popt,pcov=curve_fit(expModGauss,binCenters,binEntries,
                            sigma=np.power(binEntries,0.5))
        if np.std(spec)>1.0:
            ax1.plot(x_fine,expModGauss(x_fine,np.sqrt(np.var(spec)-1.0)),label='Estimated')
        if popt[0]>0.0:
            plt.plot(x_fine,expModGauss(x_fine,popt[0]),label='Fitted')
        plt.legend()
        plt.yscale('log')
        plt.xlim(min(bins),max(bins))
        plt.xlabel("Intensity")
        plt.show()


    # Plot fourier transforms of spectra
    if spec.shape[-1]==4:
        f,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex='col',sharey='row')
        sample=np.fft.fftshift(np.fft.fftfreq(
                spec.shape[0],(freqBand[1]-freqBand[0])/spec.shape[0]))
        ax1.plot(sample,np.abs(np.fft.fftshift(np.fft.fft(spec[:,0]))),'.')
        ax1.set_ylabel('abs(Fourier Amplitude)')
        ax1.set_title('Polarization 0')
        ax2.plot(sample,np.abs(np.fft.fftshift(np.fft.fft(spec[:,3]))),'.')
        ax2.set_title('Polarization 3')
        ax3.plot(sample,np.abs(np.fft.fftshift(np.fft.fft(spec[:,1]))),'.')
        ax3.set_ylabel('abs(Fourier Amplitude)')
        ax3.set_xlabel('Delay (microseconds)')
        ax3.set_title('Polarization 1')
        ax4.plot(sample,np.abs(np.fft.fftshift(np.fft.fft(spec[:,2]))),'.')
        ax4.set_xlabel('Delay (microseconds)')
        ax4.set_title('Polarization 2')
        plt.suptitle('Fourier Transforms of Spectra')
        plt.show()
        
        # Treat cross products as complex number and fourier transform
        complexSpec=spec[:,1]+1j*spec[:,2]
        plt.figure()
        plt.plot(sample,np.abs(np.fft.fftshift(np.fft.fft(complexSpec))),'.')
        plt.xlabel('Delay (microseconds)')
        plt.ylabel('abs(Fourier Amplitude)')
        plt.title('Fourier transform of complex cross-spectra')
        plt.show()
        amplitude=np.abs(np.fft.fftshift(np.fft.fft(complexSpec)))
        amplitude=(amplitude-np.mean(amplitude))/np.std(amplitude)
        delay=int(round(sample[np.argmax(amplitude)]*1000))
        if np.amax(amplitude)<5 or delay==0:
            print "No apparent offset!"
        else:
            print "Offset (R-L): "
            print "\t"+str(-delay)+ " ns ~=",
            print str(-int(round(delay/60.)))+' bytes ~=',
            print str(int(round(-299792458*delay*1e-9)))+' m / c.'
    else:
        plt.figure()
        sample=np.linspace(0,spec.shape[0]/(freqBand[1]-freqBand[0])/2,spec.shape[0]/2)
        plt.plot(sample,np.abs(np.fft.fft(spec)[:spec.shape[0]/2]),'.')
        plt.ylabel('abs(Fourier Amplitude)')
        plt.title('Fourier Transform of Spectrum')
        plt.xlabel('Delay (microseconds)')
        plt.show()
