#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pylab as plt
import pulsarAnalysis.GPs.pulseFinder as pf
import pulsarAnalysis.GPs.pulseSpec as ps

# Time to display before pulse peak in seconds
leadWidth=0.0003

# Time to display after pulse peak in seconds
trailWidth=0.0009

# Resolution to use for searching in seconds. Must be larger than or
# equal to phase bin size.
searchRes=1.0/10000

ignoreRFI=False

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
    dynamicSpec=ps.dynSpec(w,indices=pulseRange,normChan=False)
    dynamicSpec_BG=ps.dynSpec(w,indices=pulseRange_BG,normChan=False)

    if ignoreRFI:
        noRFI=ps.getRFIFreeBins(w.shape[0],telescope)
        profile=dynamicSpec[noRFI,...].sum(0)-dynamicSpec_BG[noRFI,...].sum(0)
    else:
        profile=dynamicSpec.sum(0)-dynamicSpec_BG.sum(0)

    # Plot profiles
    timeList=[(i-largestPulse)*binWidth*1e6 for i in pulseRange]
    if profile.shape[-1]==4:
        f,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex='col',sharey='row')
        ax1.plot(timeList,profile[:,0])
        ax1.set_ylabel('Intensity (Tsys)')
        ax1.set_title('Polarization 0')

        ax2.plot(timeList,profile[:,3])
        ax2.set_title('Polarization 3')

        ax3.plot(timeList,profile[:,1])
        ax3.set_xlabel('Time (microseconds)')
        ax3.set_ylabel('Intensity (Tsys)')
        ax3.set_title('Polarization 1') 

        ax4.plot(timeList,profile[:,2])
        ax4.set_xlabel('Time (microseconds)')
        ax4.set_title('Polarization 2')

        plt.suptitle('Profiles',size=16)
        plt.show()
        
    else:
        plt.plot(timeList,profile,'.')
        plt.title('Profile')
        plt.ylabel('Intensity (Tsys)')
        plt.xlabel('Frequency Channel')
        plt.show()
    
    # Calculate and plot correlation between polarizations
    normProfile0=(profile[:,0]-np.mean(profile[:,0]))/np.std(profile[:,0])
    normProfile3=(profile[:,3]-np.mean(profile[:,3]))/np.std(profile[:,3])

    #corr=np.correlate(profile[:,0],profile[:,3],mode='same')
    corr=np.correlate(normProfile0,normProfile3,mode='same')
    corr_x=np.arange(-(len(corr)-1)/2,(len(corr)-1)/2+1)
    corr=np.array([icorr/(profile.shape[0]-abs(corr_x[i])) 
          for i,icorr in enumerate(corr)])
    corr_time=binWidth*corr_x*1e6
    sortedCorr=sorted(corr,reverse=True)
    if sortedCorr[0]-sortedCorr[1]<0.05:
        fitRange=np.arange(np.argmax(corr)-int(0.00002/binWidth),
                       np.argmax(corr)+int(0.00002/binWidth))
        fitParams=np.polyfit(corr_time[fitRange],corr[fitRange],2)

        fit=[fitParams[2]+fitParams[1]*i+fitParams[0]*i*i 
             for i in corr_time[fitRange]]
        delay=fitParams[1]/2/fitParams[0]*1e3
    else:
        delay=-corr_time[np.argmax(corr)]*1e3
        
    print "Offset (R-L): "
    print "\t"+str(-delay)+ " ns ~=",
    print str(-int(round(delay/60.)))+' bytes ~=',
    print str(int(round(-299792458*delay*1e-9)))+' m / c.'

    plt.figure()
    plt.xlim(min(corr_time),max(corr_time))
    plt.plot(corr_time,corr)
    if sortedCorr[0]-sortedCorr[1]<0.05:
        plt.plot(corr_time[fitRange],fit)
    plt.ylim(-0.2,1.0)
    plt.ylabel('Correlation')
    plt.xlabel('Time Offset (microseconds)')
    plt.show()
