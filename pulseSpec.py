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

def dynSpec(f,ic,indices=None,normChan=False):
    # Finds the dynamic spectrum for foldspec and icounts arrays 'f'
    # and 'ic', over phase indices 'indices'. If 'normChan', then the
    # flux is normalized by the median in each frequency bin.

    # Get indices to use, defaulting to all bins
    if indices==None:
        indices=range(ic.shape[2])

    # Normalize foldspec for 2 polarisations
    if f.shape[-1]==4:
        n=f[:,:,indices,:].sum(0)/(ic[:,:,indices].sum(0)[:,:,np.newaxis])

    # Normalize foldspec for 1 polarisation
    else:
        n=f[:,:,indices].sum(0)/(ic[:,:,indices].sum(0)[:,:])

    # Normalize flux by noise in each frequency bin
    if normChan:        
        n_noise=np.std(n,axis=1)
        n/=n_noise[:,np.newaxis,...]
        
    return n

def getFrequencyBand(telescope):
    # Returns frequency band of a given telescope

    if telescope=="Jodrell Bank":
        return (605.,615.)
    elif telescope=="GMRT":
        return (601.66666,618.33334)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: %s foldspec" % sys.argv[0]
        # Run the code as eg: ./plotspec.py foldspec.npy.
        sys.exit(1)

    foldspec=sys.argv[1]
    icounts=foldspec.replace('foldspec', 'icount')
        
    # Folded spectrum axes: time, frequency, phase, pol=4 (XX, XY, YX, YY).
    f = np.load(foldspec)
    ic = np.load(icounts)
        
    # Get run information
    deltat=pf.getDeltaT(foldspec)
    telescope=pf.getTelescope(foldspec)
    startTime=pf.getStartTime(foldspec)
    binWidth=float(deltat)/f.shape[2]
    freqBand=getFrequencyBand(telescope)

    # Rebin to find giant pulses, then resolve pulses with finer binning
    nSearchBins=min(int(round(deltat/binWidth)),
                    int(round(deltat/searchRes)))
    f_rebin,ic_rebin=pf.rebin(f,ic,nSearchBins)
    timeSeries_rebin=pf.getTimeSeries(f_rebin,ic_rebin)
    timeSeries=pf.getTimeSeries(f,ic)
    pulseList=pf.getPulses(timeSeries_rebin,searchRes)
    pulseList=[(pf.resolvePulse(
                timeSeries,int(pos*searchRes/binWidth),
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
    dynamicSpec=dynSpec(f,ic,indices=pulseRange,
                                     normChan=normChan)

    # Plot each polarization if data is present
    if f.shape[-1]==4:
        for i in range(4):
            plt.imshow(dynamicSpec[:,:,i],aspect='auto',origin='lower',
                       interpolation='nearest',cmap=plt.get_cmap('Greys'),
                       extent=[-leadBins,trailBins-1,freqBand[0],freqBand[1]])
            plt.title('Dynamic Spectrum (Polarization '+str(i)+')')
            plt.xlabel('Time')
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
