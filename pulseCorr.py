#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import pulseFinder as pf
import pulseSpec as ps
import itertools
import os

# Time to display before pulse peak in seconds
leadWidth=0.0001

# Time to display after pulse peak in seconds
trailWidth=0.0003

# Resolution to use for searching in seconds. Must be larger than or
# equal to phase bin size.
searchRes=1.0/10000

# Only take the brightest pulse in each file
onePulsePerFile=True

# Find correlation coefficient between two spectra
def corr(spec1,spec2):
    spec1_norm=(spec1-np.mean(spec1))/np.std(spec1)
    spec2_norm=(spec2-np.mean(spec2))/np.std(spec2)
    normFactor1=np.correlate(spec1_norm,spec1_norm)
    normFactor2=np.correlate(spec2_norm,spec2_norm)
    normFactor=np.sqrt(normFactor1*normFactor2)
    return np.correlate(spec1_norm,spec2_norm)/normFactor

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: %s waterfallDir/" % sys.argv[0]
        # Run the code as eg: ./pulseCorr.py waterfallDir/ .
        sys.exit(1)
    
    # Declare observation list, and dynamic spectra, frequency band, pulse time, and RFI-free channel dictionaries
    pulseDict={}
    freqBand={}
    print "Finding Pulses..."
    # Loop through JB and GMRT files
    for ifilename in os.listdir(sys.argv[1]):
        
        # Folded spectrum axes: time, frequency, phase, pol=4 (XX, XY, YX, YY).

        # Get run information
        deltat=pf.getDeltaT(ifilename)
        telescope=pf.getTelescope(ifilename)
        startTime=pf.getStartTime(ifilename)

        if 'foldspec' in ifilename:
            f = np.load(sys.argv[1]+ifilename)
            ic = np.load(sys.argv[1]+ifilename.replace('foldspec', 'icount'))

            # Collapse time axis
            f=f[0,...]
            ic=ic[0,...]

            # Find populated bins
            fullList=np.flatnonzero(ic.sum(0).sum(0))
            w=f/ic[...,np.newaxis]

            binWidth=deltat/f.shape[1]

        elif 'waterfall' in ifilename:
            w=np.load(sys.argv[1]+ifilename)
            w=np.swapaxes(w,0,1)
            fullList=range(w.shape[1])
            binWidth=pf.getWaterfallBinWidth(telescope,w.shape[0])

        else:
            print "Error, unrecognized file type."
            sys.exit()

        # Check for polarization data
        if not w.shape[-1]==4:
            print "Error, polarization data is missing for "+startTime+"."
            sys.exit()

        # Rebin to find giant pulses
        nSearchBins=min(w.shape[1],int(round(deltat/searchRes)))
    
        w_rebin=pf.rebin(w,nSearchBins)
        timeSeries_rebin=pf.getTimeSeries(w_rebin)
        timeSeries=pf.getTimeSeries(w)
        pulseList=pf.getPulses(timeSeries_rebin,binWidth=searchRes,threshold=20)
        if nSearchBins<w.shape[1]:
            pulseList=[(pf.resolvePulse(
                        timeSeries,int(pos*w.shape[1]/nSearchBins),
                        binWidth=binWidth,searchRadius=1.0/10000),height) 
                       for (pos,height) in pulseList]
        if len(pulseList)==0:
            print "Warning, no pulse found for start time:"
            print startTime.iso
        elif onePulsePerFile:
            pulseList=[pulseList[0]]

        # Find range of pulse to examine
        leadBins=int(leadWidth/binWidth)
        trailBins=int(trailWidth/binWidth)

        # Assign pulse values in dictionary
        for i in pulseList:
            pulseTime=pf.getTime(i[0],binWidth,startTime)
            pulseRange=range(i[0]-leadBins,i[0]+trailBins)
            offRange=range(i[0]-leadBins-(leadBins+trailBins),i[0]-leadBins)

            pulseDict[pulseTime]={}

            pulseDict[pulseTime]['dynSpec']=ps.dynSpec(w,indices=pulseRange,
                                                       normChan=False)
            pulseDict[pulseTime]['spec']=pulseDict[pulseTime]['dynSpec'].sum(1)
            pulseDict[pulseTime]['prof']=pulseDict[pulseTime]['dynSpec'].sum(0)

            pulseDict[pulseTime]['dynSpec_bg']=ps.dynSpec(w,indices=offRange,
                                                       normChan=False)
            pulseDict[pulseTime]['spec_bg']=pulseDict[pulseTime][
                'dynSpec_bg'].sum(1)
            pulseDict[pulseTime]['prof_bg']=pulseDict[pulseTime][
                'dynSpec_bg'].sum(0)
       
    print "Complete!\n"
    print "Calculating correlations..."
    
    # Find all pairs of pulses
    pulsePairs=itertools.combinations(pulseDict.keys(),2)

    # Create lists of time difference and correlation coefficient for
    # each pair of pulses
    dt=[]
    corrCoef=[]
    for (i,j) in pulsePairs:
        dt.append(abs((i-j).sec))
        spec1=(pulseDict[i]['spec'][:,(0,3)]-pulseDict[i][
                'spec_bg'][:,(0,3)]).sum(-1)
        spec2=(pulseDict[j]['spec'][:,(0,3)]-pulseDict[j][
                'spec_bg'][:,(0,3)]).sum(-1)
        corrCoef.append(corr(spec1,spec2))
    print "Complete!"
    
    # Plot results
    plt.figure()
    plt.scatter(dt,corrCoef)
    plt.xlim(0.01,100)
    plt.xscale('log')
    plt.ylabel('Correlation Coefficient')
    plt.xlabel('Time lag (s)')
    plt.show()
    
