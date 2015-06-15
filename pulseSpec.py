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
searchRes=1.0/100000

# Use same intensity color scale
sameColorScale=False

# Scale Jodrell Bank to same intensity as GMRT
scaleJB=True

def dynspec(f,ic,indices=None,normChan=False):
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

    # Normalize flux by median in each frequency bin
    if normChan:
        n_median=np.median(n,axis=1)
        if f.shape[-1]==4:
            n_median[:,(1,2)]=1
        n/=n_median[:,np.newaxis,...]

    return n

def getFrequencyBand(telescope):
    # Returns frequency band of a given telescope

    if telescope=="Jodrell Bank":
        return (605.,615.)
    elif telescope=="GMRT":
        return (601.66666,618.33334)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: %s foldspec_jb foldspec_gmrt" % sys.argv[0]
        # Run the code as eg: ./plotspec.py foldspec_jb.npy foldspec_gmrt.py.
        sys.exit(1)
    
    # Declare frequency band and dynamic spectra dictionaries
    dynamicSpec={}
    freqBand={}

    # Loop through JB and GMRT files
    for ifoldspec in sys.argv[1:]:
        
        # Folded spectrum axes: time, frequency, phase, pol=4 (XX, XY, YX, YY).
        f = np.load(ifoldspec)
        ic = np.load(ifoldspec.replace('foldspec', 'icount'))
        
        # Get run information
        deltat=pf.getDeltaT(ifoldspec)
        telescope=pf.getTelescope(ifoldspec)
        binWidth=float(deltat)/f.shape[2]

        # Check for polarization data
        if not f.shape[-1]==4:
            print "Error, polarization data is missing for "+telescope+"."
            sys.exit()

        # Rebin to find giant pulses
        f_rebin,ic_rebin=pf.rebin(f,ic,int(round(deltat/searchRes)))
        timeSeries=pf.getTimeSeries(f_rebin,ic_rebin)
        pulseList=pf.getPulses(timeSeries,searchRes)
        try:
            largestPulse=int(pulseList[0][0]*searchRes/binWidth)
        except IndexError:
            print "Error, no giant pulse found in "+telescope+"."
            sys.exit

        # Find range of pulse to plot
        leadBins=int(leadWidth/binWidth)
        trailBins=int(trailWidth/binWidth)       
        pulseRange=range(largestPulse-leadBins,largestPulse+trailBins)
        
        # Add entries to dynamic spectra and frequency band dictionaries
        dynamicSpec[telescope]=dynspec(f,ic,indices=pulseRange,normChan=True)
        freqBand[telescope]=getFrequencyBand(telescope)

    # Determine aspect ratio for plotting
    aspect=(leadBins+trailBins-1)/(
        freqBand["Jodrell Bank"][1]-freqBand["Jodrell Bank"][0])

    # Loop through polarisations to plot
    for i in range(4):

        fig = plt.figure()
        if scaleJB:
            scaleFactor=np.amax(dynamicSpec["GMRT"][:,:,(0,3)].sum(-1))/np.amax(
                dynamicSpec["Jodrell Bank"][:,:,(0,3)].sum(-1))
            dynamicSpec["Jodrell Bank"]*=scaleFactor
            

        # Loop through telescopes
        for j,jtel in enumerate(["Jodrell Bank","GMRT"]):
            
            # Declare max and min z axis values for plotting
            if sameColorScale:
                vmin=np.amin(dynamicSpec["GMRT"][:,:,i])
                vmax=np.amax(dynamicSpec["GMRT"][:,:,i])
            else:
                vmin=np.amin(dynamicSpec[jtel][:,:,i])
                vmax=np.amax(dynamicSpec[jtel][:,:,i])

            # Plot dynamic spectra on subplot, with color bars and labels
            fig.add_subplot(121+j)
            plt.imshow(dynamicSpec[jtel][:,:,i],interpolation='nearest',
                       origin='lower',cmap=plt.get_cmap('Greys'),
                       extent=[-leadBins,trailBins-1,
                                freqBand[jtel][0],freqBand[jtel][1]],
                       aspect=aspect,vmin=vmin,vmax=vmax)
            plt.title(jtel+' ( Pol '+str(i)+' )')
            plt.xlabel('Time')
            plt.ylabel('Frequency')
            if not sameColorScale:
                plt.colorbar()
                
        if sameColorScale:
            plt.colorbar()

        # Prevent overlapping and show figure
        fig.tight_layout()
        plt.show()
