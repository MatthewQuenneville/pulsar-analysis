#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pylab as plt
import pulsarAnalysis.GPs.pulseFinder as pf

# Time to display before pulse peak in seconds
leadWidth=0.0001

# Time to display after pulse peak in seconds
trailWidth=0.0003

# Resolution to use for searching in seconds. Must be larger than or
# equal to phase bin size.
searchRes=1.0/10000

def dynSpec(w,indices=None,normChan=False):
    # Finds the dynamic spectrum for foldspec and icounts arrays 'f'
    # and 'ic', over phase indices 'indices'. If 'normChan', then the
    # flux is normalized by the median in each frequency bin.

    # Get indices to use, defaulting to all bins
    if indices==None:
        indices=range(w.shape[1])

    if w.shape[-1]==4:
        Tsys=w[:,:,(0,3)].sum(-1).mean()
    else:
        Tsys=w.mean()
    n=w[:,indices,...]/Tsys

    # Normalize flux by noise in each frequency bin
    if normChan:        
        n_median=np.median(n,axis=1)
        if w.shape[-1]==4:
            n_median[...,(1,2)]=1
        n/=n_median[:,np.newaxis,...]
        
    return n

def getRFIFreeBins(nChan,telescope):
    # Returns a list of bins that should be roughly free of RFI. This
    # is only approximate, and should only be used for things such as
    # selecting vmax and vmin for plotting.

    freqBand=pf.getFrequencyBand(telescope)
    if telescope=="Jodrell Bank":
        RFI=[(605.,606.5), (614.,615.)]
    if telescope=="GMRT":
        RFI=[(freqBand[0],freqBand[0]+1.0/2048)]
    chanWidth=(freqBand[1]-freqBand[0])/nChan
    bandMin=freqBand[0]
    cleanChan=range(nChan)
    for iBand in RFI:
        for chan in range(nChan):
            if not chan in cleanChan:
                continue
            binRng=(chan*chanWidth+bandMin,(chan+1)*chanWidth+bandMin)
            if binRng[0]<=iBand[0]<=binRng[1] or iBand[0]<=binRng[1]<=iBand[1]:
                cleanChan.remove(chan)
    return cleanChan

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
    leadBins=int(leadWidth/binWidth)
    trailBins=int(trailWidth/binWidth)       
    pulseRange=range(largestPulse-leadBins,largestPulse+trailBins)
    pulseRange_BG=range(largestPulse-2*leadBins-trailBins,largestPulse-leadBins)
    
    # Add entries to dynamic spectra and frequency band dictionaries
    dynamicSpec=dynSpec(w,indices=pulseRange,normChan=False)
    dynamicSpec_BG=dynSpec(w,indices=pulseRange_BG,normChan=False)

    # Get minimum and maximum intensity to plot, ignoring RFI channels
    cleanChans=getRFIFreeBins(w.shape[0],telescope)
    vmin=np.amin(np.amin(dynamicSpec[cleanChans,...],axis=1),axis=0)
    vmax=np.amax(np.amax(dynamicSpec[cleanChans,...],axis=1),axis=0)
    
    # Get min and max of time range:
    tmin=-leadBins*binWidth
    tmax=trailBins*binWidth

    # Plot each polarization if data is present
    if w.shape[-1]==4:
        f,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex='col',sharey='row')
        im=ax1.imshow(dynamicSpec[:,:,0],aspect='auto',origin='lower',
                   interpolation='nearest',cmap=plt.get_cmap('Greys'),
                   extent=[tmin*1e6,tmax*1e6,freqBand[0],freqBand[1]],
                   vmin=vmin[0], vmax=vmax[0])
        plt.colorbar(im,ax=ax1)
        ax1.set_ylabel('Frequency (MHz)')
        ax1.set_title('Polarization 0')

        im=ax2.imshow(dynamicSpec[:,:,3],aspect='auto',origin='lower',
                   interpolation='nearest',cmap=plt.get_cmap('Greys'),
                   extent=[tmin*1e6,tmax*1e6,freqBand[0],freqBand[1]],
                   vmin=vmin[3], vmax=vmax[3])
        plt.colorbar(im,ax=ax2)
        ax2.set_title('Polarization 3')

        im=ax3.imshow(dynamicSpec[:,:,1],aspect='auto',origin='lower',
                   interpolation='nearest',cmap=plt.get_cmap('Greys'),
                   extent=[tmin*1e6,tmax*1e6,freqBand[0],freqBand[1]],
                   vmin=vmin[1], vmax=vmax[1])
        plt.colorbar(im,ax=ax3)
        ax3.set_xlabel('Time (microseconds)')
        ax3.set_ylabel('Frequency (MHz)')
        ax3.set_title('Polarization 1')

        im=ax4.imshow(dynamicSpec[:,:,2],aspect='auto',origin='lower',
                   interpolation='nearest',cmap=plt.get_cmap('Greys'),
                   extent=[tmin*1e6,tmax*1e6,freqBand[0],freqBand[1]],
                   vmin=vmin[2], vmax=vmax[2])
        plt.colorbar(im,ax=ax4)
        ax4.set_xlabel('Time (microseconds)')
        ax4.set_title('Polarization 2')

        plt.suptitle('Dynamic Spectra',size=16)
        plt.show()

    # Plot intensity if no polarization data is present
    else:
        plt.imshow(dynamicSpec[:,:],aspect='auto',origin='lower',
                   interpolation='nearest',cmap=plt.get_cmap('Greys'),
                   extent=[tmin*1e6,tmax*1e6,freqBand[0],freqBand[1]],
                   vmin=vmin, vmax=vmax)
        plt.title('Dynamic Spectrum')
        plt.xlabel('Time (microseconds)')
        plt.ylabel('Frequency (MHz)') 
        plt.show()
