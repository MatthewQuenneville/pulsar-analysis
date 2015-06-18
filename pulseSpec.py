#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pylab as plt
import pulseFinder as pf

# Time to display before pulse peak in seconds
leadWidth=0.0001

# Time to display after pulse peak in seconds
trailWidth=0.0003

# Resolution to use for searching in seconds. Must be larger than or
# equal to phase bin size.
searchRes=1.0/10000

# Normalize intensity in frequency channels
normChan=True

def dynSpec(w,indices=None,normChan=False):
    # Finds the dynamic spectrum for foldspec and icounts arrays 'f'
    # and 'ic', over phase indices 'indices'. If 'normChan', then the
    # flux is normalized by the median in each frequency bin.

    # Get indices to use, defaulting to all bins
    if indices==None:
        indices=range(ic.shape[2])

    n=w[:,indices,...]

    # Normalize flux by noise in each frequency bin
    if normChan:        
        n_median=np.median(n,axis=1)
        if w.shape[-1]==4:
            n_median[...,(1,2)]=1
        n/=n_median[:,np.newaxis,...]
        
    return n

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: %s foldspec" % sys.argv[0]
        # Run the code as eg: ./pulseSpec.py foldspec.npy.
        sys.exit(1)

    # Get run information
    deltat=pf.getDeltaT(sys.argv[1])
    telescope=pf.getTelescope(sys.argv[1])
    startTime=pf.getStartTime(sys.argv[1])
    freqBand=pf.getFrequencyBand(telescope)

    # Folded spectrum axes: time, frequency, phase, pol=4 (XX, XY, YX, YY).
    if 'foldspec' in sys.argv[1]:
        f = np.load(sys.argv[1])
        ic = np.load(sys.argv[1].replace('foldspec', 'icount'))

        # Collapse time axis
        f=f[0,...]
        ic=ic[0,...]

        # Find populated bins
        fullList=np.flatnonzero(ic.sum(0).sum(0))
        w=f/ic[...,np.newaxis]

        binWidth=deltat/f.shape[1]

    elif 'waterfall' in sys.argv[1]:
        w=np.load(sys.argv[1])
        w=np.swapaxes(w,0,1)
        fullList=range(w.shape[1])
        binWidth=pf.getWaterfallBinWidth(telescope,w.shape[0])

    else:
        print "Error, unrecognized file type."
        sys.exit()

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
        sys.exit

    # Find range of pulse to plot
    leadBins=int(leadWidth/binWidth)
    trailBins=int(trailWidth/binWidth)       
    pulseRange=range(largestPulse-leadBins,largestPulse+trailBins)
        
    # Add entries to dynamic spectra and frequency band dictionaries
    dynamicSpec=dynSpec(w,indices=pulseRange,normChan=normChan)

    # Plot each polarization if data is present
    if w.shape[-1]==4:
        for i in range(4):
            plt.imshow(dynamicSpec[:,:,i],aspect='auto',origin='lower',
                       interpolation='nearest',cmap=plt.get_cmap('Greys'),
                       extent=[-leadWidth*1e6,trailWidth*1e6,
                                freqBand[0],freqBand[1]])
            plt.title('Dynamic Spectrum (Polarization '+str(i)+')')
            plt.xlabel('Time (ns)')
            plt.ylabel('Frequency (MHz)')                      
            plt.show()

    # Plot intensity if no polarization data is present
    else:
        plt.imshow(dynamicSpec[:,:],aspect='auto',origin='lower',
                   interpolation='nearest',cmap=plt.get_cmap('Greys'),
                   extent=[-leadBins,trailBins-1,freqBand[0],freqBand[1]])
        plt.title('Dynamic Spectrum')
        plt.xlabel('Time')
        plt.ylabel('Frequency (MHz)') 
        plt.show()
