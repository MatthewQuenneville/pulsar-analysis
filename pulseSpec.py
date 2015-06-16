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

# Scale data sets to have same intensity
scaleData=False

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
    if len(sys.argv) < 3:
        print "Usage: %s foldspec1 foldspec2" % sys.argv[0]
        # Run the code as eg: ./plotspec.py foldspec1.npy foldspec2.py.
        sys.exit(1)
    
    # Declare observation list, and dynamic spectra, frequency band,
    # and pulse time dictionaries
    obsList=[]
    dynamicSpec={}
    freqBand={}
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
        freqBand[obsList[-1]]=getFrequencyBand(telescope)

        # Check for polarization data
        if not f.shape[-1]==4:
            print "Error, polarization data is missing for "+telescope+"."
            sys.exit()

        # Rebin to find giant pulses
        nSearchBins=min(int(round(deltat/binWidth)),
                        int(round(deltat/searchRes)))
        f_rebin,ic_rebin=pf.rebin(f,ic,nSearchBins)
        timeSeries=pf.getTimeSeries(f_rebin,ic_rebin)
        pulseList=pf.getPulses(timeSeries,searchRes)
        try:
            largestPulse=int(pulseList[0][0]*searchRes/binWidth)
        except IndexError:
            print "Error, no giant pulse found in "+telescope+" for start time:"
            print startTime.iso
            sys.exit

        # Find range of pulse to plot
        leadBins=int(leadWidth/binWidth)
        trailBins=int(trailWidth/binWidth)       
        pulseRange=range(largestPulse-leadBins,largestPulse+trailBins)
        
        # Add entries to dynamic spectra and frequency band dictionaries
        dynamicSpec[obsList[-1]]=dynSpec(f,ic,indices=pulseRange,
                                         normChan=normChan)
        pulseTimes[obsList[-1]]=(pf.getTime(pulseList[0][0],searchRes,
                                           startTime).iso[:-3]).split()[-1]

    # Determine aspect ratio for plotting
    freqRange=[b-a for (a,b) in freqBand.values()]
    maxFreqRange=max(freqRange)
    aspect=2*(leadBins+trailBins-1)/maxFreqRange

    # Apply scaling factor to second data set
    if scaleData:
        scaleFactor=np.amax(dynamicSpec[obsList[0]][:,:,(0,3)].sum(-1))/np.amax(
            dynamicSpec[obsList[1]][:,:,(0,3)].sum(-1))
        dynamicSpec[obsList[1]]*=scaleFactor

    # Declare max and min z axis values for plotting
    vmin=[[],[]]
    vmax=[[],[]]

    for i in range(4):
        if sameColorScale:
            vmin[0].append(min(np.amin(dynamicSpec[obsList[0]][:,:,i]),
                np.amin(dynamicSpec[obsList[1]][:,:,i])))
            vmax[0].append(max(np.amax(dynamicSpec[obsList[0]][:,:,i]),
                np.amax(dynamicSpec[obsList[1]][:,:,i])))
            vmin[1].append(min(np.amin(dynamicSpec[obsList[0]][:,:,i]),
                np.amin(dynamicSpec[obsList[1]][:,:,i])))
            vmax[1].append(max(np.amax(dynamicSpec[obsList[0]][:,:,i]),
                np.amax(dynamicSpec[obsList[1]][:,:,i])))
        else:
            vmin[0].append(np.amin(dynamicSpec[obsList[0]][:,:,i]))
            vmax[0].append(np.amax(dynamicSpec[obsList[0]][:,:,i]))
            vmin[1].append(np.amin(dynamicSpec[obsList[1]][:,:,i]))
            vmax[1].append(np.amax(dynamicSpec[obsList[1]][:,:,i]))

    # Find difference in frequency range and channel widths
    upperDiff=freqBand[obsList[1]][0]-freqBand[obsList[0]][0]
    lowerDiff=freqBand[obsList[1]][1]-freqBand[obsList[0]][1]
    chanWidth=[]
    chanWidth.append(
        (freqBand[obsList[0]][1]-freqBand[obsList[0]][0])/
        dynamicSpec[obsList[0]].shape[0])
    chanWidth.append(
        (freqBand[obsList[0]][1]-freqBand[obsList[0]][0])/
        dynamicSpec[obsList[0]].shape[0])

    if not chanWidth[0]==chanWidth[1]:
        print "Warning, channel widths are not equal."
    
    # Pad smaller image with constant values to maintain aspect ratio
    if upperDiff>0:
        nBinsUpper=int(round(upperDiff/chanWidth[1]))
        dynamicSpec[obsList[1]]=np.pad(dynamicSpec[obsList[1]],
                                       ((0,nBinsUpper),(0,0),(0,0)),
                                       mode='constant',
                                       constant_values=
                                       ((0,vmin[1]),(0,0),(0,0)))
    elif upperDiff<0:
        nBinsUpper=int(round(-upperDiff/chanWidth[0]))
        dynamicSpec[obsList[0]]=np.pad(dynamicSpec[obsList[0]],
                                       ((0,nBinsUpper),(0,0),(0,0)),
                                       mode='constant',
                                       constant_values=
                                       ((0,vmin[0]),(0,0),(0,0)))
    if lowerDiff<0:
        nBinsLower=int(round(-lowerDiff/chanWidth[1]))
        dynamicSpec[obsList[1]]=np.pad(dynamicSpec[obsList[1]],
                                       ((nBinsLower,0),(0,0),(0,0)),
                                       mode='constant',
                                       constant_values=
                                       ((vmin[1],0),(0,0),(0,0)))
    elif lowerDiff>0:
        nBinsLower=int(round(lowerDiff/chanWidth[0]))
        dynamicSpec[obsList[0]]=np.pad(dynamicSpec[obsList[0]],
                                       ((nBinsLower,0),(0,0),(0,0)),
                                       mode='constant',
                                       constant_values=
                                       ((vmin[0],0),(0,0),(0,0)))

    # Loop through polarisations to plot
    for i in range(4):

        fig,axes = plt.subplots(nrows=1,ncols=2)
            
        ymin=min([ifreq[0] for ifreq in freqBand.values()])
        ymax=max([ifreq[1] for ifreq in freqBand.values()])

        # Loop through telescopes
        for j,jobs in enumerate(obsList):

            # Plot image and set titles
            im=axes.flat[j].imshow(dynamicSpec[jobs][:,:,i],origin='lower',
                       interpolation='nearest',cmap=plt.get_cmap('Greys'),
                       extent=[-leadBins,trailBins-1,ymin,ymax],
                       aspect=aspect,vmin=vmin[j][i],vmax=vmax[j][i])
            axes.flat[j].set_title(pulseTimes[jobs]+'\n'+jobs[1]+
                                   ' ( Pol '+str(i)+' )')
            axes.flat[j].set_xlabel('Time')
            axes.flat[j].set_ylabel('Frequency')
            # Plot color bar on each subplot if color scales are different
            if not sameColorScale:
                plt.colorbar(im,ax=axes.flat[j])

        # Add single color bar if color scales are equal, otherwise
        # prevent overlap of colorbars with labels
        if sameColorScale:
            cax = fig.add_axes([0.9,0.1,0.03,0.8])
            fig.colorbar(im,cax=cax)
        else:
            plt.tight_layout()
        
        plt.show()
        """
            # Plot dynamic spectra on subplot, with color bars and labels
            fig.add_subplot(121+j)
            im=plt.imshow(dynamicSpec[jobs][:,:,i],origin='lower',
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
            #fig.add_subplot(155)
            plt.colorbar()

        # Prevent overlapping and show figure
        #fig.tight_layout()
        plt.show()
        """
