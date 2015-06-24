#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pylab as plt
import pulseFinder as pf
from math import factorial
import warnings
doFits=True
try:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from lmfit.models import ExponentialGaussianModel
except ImportError:
    print "lmfit module not found. Exponentially modified gaussian fits will not be performed."
    doFits=False

# Time to display before pulse peak in seconds
leadWidth=0.0001

# Time to display after pulse peak in seconds
trailWidth=0.0003

# Resolution to use for searching in seconds. Must be larger than or
# equal to phase bin size.
searchRes=1.0/10000

# Normalize intensity in frequency channels
normChan=False

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
        sys.exit()

    # Find range of pulse to plot
    leadBins=int(leadWidth/binWidth)
    trailBins=int(trailWidth/binWidth)       
    pulseRange=range(largestPulse-leadBins,largestPulse+trailBins)
    pulseRange_BG=range(largestPulse-2*leadBins-trailBins,largestPulse-leadBins)
        
    # Add entries to dynamic spectra and frequency band dictionaries
    dynamicSpec=dynSpec(w,indices=pulseRange,normChan=normChan)
    dynamicSpec_BG=dynSpec(w,indices=pulseRange_BG,normChan=normChan)

    # Get minimum and maximum intensity to plot, ignoring RFI channels
    cleanChans=getRFIFreeBins(w.shape[0],telescope)
    vmin=np.amin(np.amin(dynamicSpec[cleanChans,...],axis=1),axis=0)
    vmax=np.amax(np.amax(dynamicSpec[cleanChans,...],axis=1),axis=0)

    # Plot each polarization if data is present
    if w.shape[-1]==4:
        f,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex='col',sharey='row')
        im=ax1.imshow(dynamicSpec[:,:,0],aspect='auto',origin='lower',
                   interpolation='nearest',cmap=plt.get_cmap('Greys'),
                   extent=[-leadWidth*1e6,trailWidth*1e6,
                            freqBand[0],freqBand[1]],
                   vmin=vmin[0], vmax=vmax[0])
        plt.colorbar(im,ax=ax1)
        ax1.set_ylabel('Frequency (MHz)')
        ax1.set_title('Polarization 0')

        im=ax2.imshow(dynamicSpec[:,:,3],aspect='auto',origin='lower',
                   interpolation='nearest',cmap=plt.get_cmap('Greys'),
                   extent=[-leadWidth*1e6,trailWidth*1e6,
                            freqBand[0],freqBand[1]],
                   vmin=vmin[3], vmax=vmax[3])
        plt.colorbar(im,ax=ax2)
        ax2.set_title('Polarization 3')

        im=ax3.imshow(dynamicSpec[:,:,1],aspect='auto',origin='lower',
                   interpolation='nearest',cmap=plt.get_cmap('Greys'),
                   extent=[-leadWidth*1e6,trailWidth*1e6,
                            freqBand[0],freqBand[1]],
                   vmin=vmin[1], vmax=vmax[1])
        plt.colorbar(im,ax=ax3)
        ax3.set_xlabel('Time (ns)')
        ax3.set_ylabel('Frequency (MHz)')
        ax3.set_title('Polarization 1')

        im=ax4.imshow(dynamicSpec[:,:,2],aspect='auto',origin='lower',
                   interpolation='nearest',cmap=plt.get_cmap('Greys'),
                   extent=[-leadWidth*1e6,trailWidth*1e6,
                            freqBand[0],freqBand[1]],
                   vmin=vmin[2], vmax=vmax[2])
        plt.colorbar(im,ax=ax4)
        ax4.set_xlabel('Time (ns)')
        ax4.set_title('Polarization 2')

        plt.suptitle('Dynamic Spectra',size=16)
        plt.show()

    # Plot intensity if no polarization data is present
    else:
        plt.imshow(dynamicSpec[:,:],aspect='auto',origin='lower',
                   interpolation='nearest',cmap=plt.get_cmap('Greys'),
                   extent=[-leadBins,trailBins-1,freqBand[0],freqBand[1]],
                   vmin=vmin, vmax=vmax)
        plt.title('Dynamic Spectrum')
        plt.xlabel('Time')
        plt.ylabel('Frequency (MHz)') 
        plt.show()
        
    # Plot spectra
    spec=dynamicSpec.sum(1)-dynamicSpec_BG.sum(1)
    if spec.shape[-1]==4:
        f,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex='col',sharey='row')
        ax1.plot(spec[:,0],'.')
        ax1.set_ylabel('Intensity (Tsys)')
        ax1.set_title('Polarization 0')

        ax2.plot(spec[:,3],'.')
        ax2.set_title('Polarization 3')

        im=ax3.plot(spec[:,1],'.')
        ax3.set_xlabel('Frequency Channel')
        ax3.set_ylabel('Intensity (Tsys)')
        ax3.set_title('Polarization 1') 

        im=ax4.plot(spec[:,2],'.')
        ax4.set_xlabel('Frequency Channel')
        ax4.set_title('Polarization 2')

        plt.suptitle('Spectra',size=16)
        plt.show()
        
    else:
        plt.plot(spec,'.')
        plt.title('Spectrum')
        plt.ylabel('Intensity (Tsys)')
        plt.xlabel('Frequency Channel')
        plt.show()

    if spec.shape[-1]==4:
        # Normalize intensity
        spec1=2*spec[:,0]/np.mean(spec[:,0])
        spec2=2*spec[:,3]/np.mean(spec[:,3])
                               
        # Plot histograms and plot if lmfit is found
        xmin=min(min(spec2),min(spec1))
        xmax=max(max(spec2),max(spec1))
        
        bins=np.linspace(np.floor(xmin),np.ceil(xmax),50)
        binCenters=np.array([(bins[i+1]+bins[i])/2. for i in range(len(bins)-1)])
        f,(ax1,ax2) = plt.subplots(1,2,sharex='col',sharey='row')
        specHist1=ax1.hist(spec1,normed=True,bins=bins)
        ax1.set_yscale('log')
        if doFits:
            expGaussMod=ExponentialGaussianModel()
            pars=expGaussMod.guess(specHist1[0],x=binCenters)
            pars['amplitude'].set(value=1.0,vary=False)
            pars['gamma'].set(value=0.5,vary=False)
            pars['center'].set(value=0.0,vary=False)
            out=expGaussMod.fit(specHist1[0],pars,x=binCenters)
            print "===================="
            print "Polarization 0"
            print "Constrained Parameters:"
            print "Amplitude: 1.0"
            print "Lambda: 0.5"
            print "Center: 0.0\n"
            print "Free Parameters:"
            print "Sigma: ", out.params.valuesdict()['sigma']
            print "====================\n"
            ax1.plot(binCenters, out.best_fit)
        ax1.set_xlim(min(bins),max(bins))
        ax1.set_xlabel("Intensity")
        specHist2=ax2.hist(spec2,normed=True,bins=bins)
        if doFits:
            expGaussMod=ExponentialGaussianModel()
            pars=expGaussMod.guess(specHist2[0],x=binCenters)
            pars['amplitude'].set(value=1.0,vary=False)
            pars['gamma'].set(value=0.5,vary=False)
            pars['center'].set(value=0.0,vary=False)
            out=expGaussMod.fit(specHist2[0],pars,x=binCenters)
            print "===================="
            print "Polarization 3"
            print "Constrained Parameters:"
            print "Amplitude: 1.0"
            print "Lambda: 0.5"
            print "Center: 0.0\n"
            print "Free Parameters:"
            print "Sigma: ", out.params.valuesdict()['sigma']
            print "====================\n"
            ax2.plot(binCenters, out.best_fit)
        ax2.set_yscale('log')
        ax2.set_xlim(min(bins),max(bins))
        ax2.set_xlabel("Intensity")
        plt.show()
    else:
        # Normalize intensity
        spec=spec/np.mean(spec)

        # Plot histograms and fit if lmfit is found
        xmin=min(spec)
        xmax=max(spec)
        
        bins=np.linspace(np.floor(xmin),np.ceil(xmax),50)
        binCenters=np.array([(bins[i+1]+bins[i])/2. for i in range(len(bins)-1)])
        specHist=plt.hist(spec,normed=True,bins=bins)
        if doFits:
            expGaussMod=ExponentialGaussianModel()
            pars=expGaussMod.guess(specHist[0],x=binCenters)
            pars['amplitude'].set(value=1.0,vary=False)
            pars['gamma'].set(value=0.5,vary=False)
            pars['center'].set(value=0.0,vary=False)
            out=expGaussMod.fit(specHist[0],pars,x=binCenters)
            print "===================="
            print "Polarization 3"
            print "Constrained Parameters:"
            print "Amplitude: 1.0"
            print "Lambda: 0.5"
            print "Center: 0.0\n"
            print "Free Parameters:"
            print "Sigma: ", out.params.valuesdict()['sigma']
            print "====================\n"
            ax2.plot(binCenters, out.best_fit)
        plt.yscale('log')
        plt.xlim(min(bins),max(bins))
        plt.xlabel("Intensity")
        plt.show()
