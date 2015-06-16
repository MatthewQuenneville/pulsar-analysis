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
sameColorScale=True

# Scale data sets to have same intensity
scaleData=True

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
        print "Usage: %s foldspec1 foldspec2" % sys.argv[0]
        # Run the code as eg: ./plotspec.py foldspec1.npy foldspec2.py.
        sys.exit(1)
    
    # Declare frequency band and dynamic spectra dictionaries
    dynamicSpec={}
    freqBand={}
    obsList=[]
    pulseTimes={}

    # Loop through JB and GMRT files
    for ifoldspec in sys.argv[1:]:
        
        # Folded spectrum axes: time, frequency, phase, pol=4 (XX, XY, YX, YY).
        f = np.load(ifoldspec)
        ic = np.load(ifoldspec.replace('foldspec', 'icount'))
        
        # Get run information
        deltat=pf.getDeltaT(ifoldspec)
        telescope=pf.getTelescope(ifoldspec)
        startTime=pf.getStartTime(ifoldspec)
        obsList.append((startTime,telescope))
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
        dynamicSpec[obsList[-1]]=dynspec(f,ic,indices=pulseRange,
                                         normChan=True)
        freqBand[obsList[-1]]=getFrequencyBand(telescope)
        pulseTimes[obsList[-1]]=pf.getTime(pulseList[0][1],binWidth,
                                           startTime).iso[:-3]

    # Determine aspect ratio for plotting
    freqRange=[b-a for (a,b) in freqBand.values()]
    maxFreqRange=max(freqRange)
    aspect=2*(leadBins+trailBins-1)/maxFreqRange

    # Loop through polarisations to plot
    for i in range(4):

        fig = plt.figure()
        
        # Apply scaling factor to second data set
        if scaleData:
            scaleFactor=np.amax(
                dynamicSpec[obsList[0]][:,:,(0,3)].sum(-1))/np.amax(
                dynamicSpec[obsList[1]][:,:,(0,3)].sum(-1))
            dynamicSpec[obsList[1]]*=scaleFactor
            

        # Loop through telescopes
        for j,jobs in enumerate(obsList):
            
            # Declare max and min z axis values for plotting
            if sameColorScale:
                vmin=min(
                    np.amin(dynamicSpec[obsList[0]][:,:,i]),
                    np.amin(dynamicSpec[obsList[1]][:,:,i]))
                vmax=max(
                    np.amax(dynamicSpec[obsList[0]][:,:,i]),
                    np.amax(dynamicSpec[obsList[1]][:,:,i]))
            else:
                vmin=np.amin(dynamicSpec[jobs][:,:,i])
                vmax=np.amax(dynamicSpec[jobs][:,:,i])

            # Plot dynamic spectra on subplot, with color bars and labels
            fig.add_subplot(121+j)
            plt.imshow(dynamicSpec[jobs][:,:,i],origin='lower',
                       interpolation='nearest',cmap=plt.get_cmap('Greys'),
                       extent=[-leadBins,trailBins-1,
                                freqBand[jobs][0],freqBand[jobs][1]],
                       aspect=aspect,vmin=vmin,vmax=vmax)
            plt.title(pulseTimes[jobs]+'\n'+jobs[1]+' ( Pol '+str(i)+' )')
            plt.xlabel('Time')
            plt.ylabel('Frequency')
            if not sameColorScale:
                plt.colorbar()
                
        if sameColorScale:
            plt.colorbar()

        # Prevent overlapping and show figure
        fig.tight_layout()
        plt.show()
