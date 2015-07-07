#!/usr/bin/env python

import sys
import os
import numpy as np
import matplotlib.pylab as plt
import string
from astropy.time import Time,TimeDelta
import warnings

# Crab frequency 
# Should implement an ephemeris/polynomial based frequency, but this
# is just used for plotting, so approximate is fine
crabFreq=29.946923 

# Pulse width in seconds 
pulseWidth=0.0001 

# Set number of pulses to find
nPulses=3

# Set noise threshold in number of standard deviations
threshold=5

# Set number of bins over which to normalize by noise. Use 0 for one
# noise bin per second
nNoiseBins=0

# Polarizations to sum over
pol_select = (0, 3)

def loadFiles(pathList,folded=False):
    if len(pathList)==0:
        print "Usage: %s foldspec" % sys.argv[0]
        # Run the code as: ./script.py data_foldspec.npy.
        sys.exit(1) 
    runInfo={}
    for i,iPath in enumerate(pathList):
        if os.path.isdir(iPath):
            iFileList=os.listdir(iPath)
        else:
            iFileList=[iPath]
        for j,jFile in enumerate(iFileList):
            if os.path.isdir(iPath):
                iFile=iPath+'/'+jFile
            else:
                iFile=jFile
            if i==0 and j==0:
                deltat=getDeltaT(iFile)
                telescope=getTelescope(iFile)
                startTime=getStartTime(iFile)
            if folded:
                f=np.load(iFile)
                ic=np.load(iFile.replace('foldspec', 'icount'))
                if i==0:
                    binWidth=deltat/f.shape[2]
                if f.shape[-1]==4:
                    w=f/ic[...,np.newaxis]
                else:
                    w=f/ic
            else:
                if 'foldspec' in iFile:
                    f=np.load(iFile)
                    ic=np.load(iFile.replace('foldspec', 'icount'))
                    f=f.sum(0)
                    ic=ic.sum(0)
                    if i==0:
                        binWidth=deltat/f.shape[1]
                    if f.shape[-1]==4:
                        w=f/ic[...,np.newaxis]
                    else:
                        w=f/ic
                elif 'waterfall' in iFile:
                    w=np.load(iFile)
                    w=np.swapaxes(w,0,1)
                    if i==0:
                        binWidth=getWaterfallBinWidth(telescope,w.shape[0])
                else:
                    print "Error, the following file name is not recognized:"
                    print iFile
            if i==0 and j==0:
                n=w
            elif n.shape==w.shape:
                n=n+w
            else:
                print "Error, shape mismatch in file:"
                print iFile

    if folded:
        if n.shape[-1]==4:
            isNotNan=~np.isnan(n.sum(-1).sum(1).sum(0))
            fullList=np.flatnonzero(isNotNan)
            n=n[:,:,isNotNan,:]
        else:
            isNotNan=~np.isnan(n.sum(1).sum(0))
            fullList=np.flatnonzero(isNotNan)
            n=n[:,:,isNotNan]
    else:
        if n.shape[-1]==4:
            isNotNan=~np.isnan(n.sum(-1).sum(0))
            fullList=np.flatnonzero(isNotNan)
            n=n[:,isNotNan,:]
        else:
            isNotNan=~np.isnan(n.sum(0))
            fullList=np.flatnonzero(isNotNan)
            n=n[:,isNotNan]
    runInfo['binWidth']=binWidth
    runInfo['telescope']=telescope
    runInfo['deltat']=deltat
    runInfo['startTime']=startTime
    runInfo['fullList']=fullList
    return n, runInfo

def getTelescope(fileName):
    # Gets telescope name based on file name

    if 'jb' in string.lower(fileName):
        return 'Jodrell Bank'
    elif 'gmrt' in string.lower(fileName):
        return 'GMRT'
    else:
        print "Telescope not recognized."
        return 'Unknown'

def getFrequencyBand(telescope):
    # Returns frequency band of a given telescope

    if telescope=="Jodrell Bank":
        return (605.,615.)
    elif telescope=="GMRT":
        return (602.0,618.66667)
        #return (601.66666,618.33334)

def getWaterfallBinWidth(telescope,nChan):
    # Returns waterfall bin width for a given telescope observation
    # with 'nChan' channels

    freqBand=getFrequencyBand(telescope)
    chanWidth=1e6*(freqBand[1]-freqBand[0])/nChan
    return 1/chanWidth

def rebin(w,nBins=10000):
    # Rebin 'w' to 'nBins' phase bins

    # Check that requested rebin is to coarser resoution
    nBinsOld=w.shape[1]
    if nBins>nBinsOld:
        print "Error, can't rebin to larger number of bins."
        print "Keeping current dimensions."
        return w

    nBinsCombine=float(nBinsOld)/nBins
    binEdges=[i*nBinsCombine for i in range(nBins+1)]
    w=np.array([w[:,np.floor(binEdges[i]):np.floor(binEdges[i+1]),...].sum(1) for i in range(nBins)])
    w=np.swapaxes(w,0,1)
    
    return w

def rms(sequence):
    # Gets root-mean square of sequence

    return np.sqrt(np.mean(np.power(sequence,2)))

def getStartTime(fileName):
    # Gets start time of trial based on file name
    if 'foldspec' in fileName:
        dateString=fileName.split('foldspec')[1][1:24]
    else:
        dateString=fileName.split('waterfall')[1][1:24]
    return Time(dateString,format='isot',scale='utc',precision=6)

def getDeltaT(fileName):
    # Gets duration of trial based on file name

    timeString=fileName.split('+')[-1].split('sec')[0]
    return float(timeString)

def nanMedian(numArray):
    # Equivalent to np.median(numArray,axis=1) for 2D array numArray,
    # but ignores any 'nan' present in the array

    medianList=[]
    for i in range(numArray.shape[0]):
        cleanedRow=numArray[i,~np.isnan(numArray[i,:])]
        medianList.append(np.median(cleanedRow))
    return np.array(medianList)

def getTime(index,binWidth,startTime):
    # Get time corresponding to 'index' based on 'startTime' and 'binWidth'
    
    return (startTime+TimeDelta(index*binWidth,format='sec'))

def getPeriod(pulseLocation,binWidth,endIndex):
    # Get approximate period worth of bins centered at 'pulseLocation'
    
    crabPeriod=1/crabFreq
    nBins=int(crabPeriod/binWidth)
    Interval=range(pulseLocation-nBins/2, pulseLocation+nBins/2)
    return [i for i in Interval if 0<=i<=endIndex]

def getTimeSeries(w,nNoiseBins=1):
    # Gets profile time series over which to search for giant pulses
    if w.shape[-1] == 4:
        n = w[..., pol_select].sum(-1)
    else:
        n = w
        
    # Normalize by median flux in each frequency bin
    n_median = nanMedian(n)
    nn = n / n_median[:,np.newaxis] - 1.

    # Sum over frequency and remove Nan entries
    timeSeries = nn.sum(0)
    timeSeries = timeSeries[~np.isnan(timeSeries)]

    # Find noise bins, and normlize by noise in each
    noiseBins=[int(i) for i in np.linspace(0,len(timeSeries),nNoiseBins+1)]
    noise=[rms(timeSeries[noiseBins[i]:noiseBins[i+1]]) 
           for i in range(nNoiseBins)]

    for i in range(nNoiseBins):
        for j in range(noiseBins[i],noiseBins[i+1]):
            timeSeries[j]/=noise[i]

    return timeSeries

def resolvePulse(timeSeries,pulseIndex,binWidth=None,searchRadius=1.0/10000):
    # Get pulse center to higher resolution

    nBins=len(timeSeries)
    if binWidth==None:
        binWidth=1.0/nBins
    
    # Ensure there are enough bins to further resolve pulse
    binRadius=int(round(searchRadius/binWidth))
    if binRadius==0:
        print "Warning, not enough bins available to further resolve pulse."
        return pulseIndex
    
    # Get range of bins to search
    binRange=range(pulseIndex-binRadius,pulseIndex+binRadius+1)  
    binRange=[i for i in binRange if i<len(timeSeries)]

    return np.argmax(timeSeries[binRange])+pulseIndex-binRadius

def getPulses(timeSeries,threshold=5,binWidth=None):
    # Gets a list of all pulses higher than the noise threshold

    nBins=len(timeSeries)

    if binWidth==None:
        binWidth=1.0/nBins
    # Define a 'mask' such that previously found pulses are ignored 
    mask=np.ones(len(timeSeries))

    pulseDict={}

    while True:
        # Find largest pulse in masked time series
        currentMax=np.amax(timeSeries*mask)

        # If larger than lowerbound, count as new pulse. Otherwise,
        # break and notify user that not enough giant pulses exist
        if currentMax>threshold:
            
            # Find location of current max
            currentMaxLoc=np.argmax(timeSeries*mask) 
            pulseDict[currentMaxLoc]=currentMax
            # Mask pulse from future searching
            mask[getPeriod(currentMaxLoc,binWidth,nBins-1)]=0.
            
        else:
            break
    pulseList=sorted(pulseDict.keys(),key=lambda x: -pulseDict[x])

    return [(i,pulseDict[i]) for i in pulseList]

if __name__ == "__main__":
    # Load files
    w,runInfo=loadFiles(sys.argv[1:])

    # Get basic run attributes
    telescope=runInfo['telescope']
    startTime=runInfo['startTime']
    deltat=runInfo['deltat']
    binWidth=runInfo['binWidth']
    fullList=runInfo['fullList']

    # Set nNoiseBins to one per second if invalid value is given
    if nNoiseBins<1:
        nNoiseBins=int(np.ceil(deltat))

    # Get timeseries to search for pulses
    timeSeries=getTimeSeries(w,nNoiseBins)

    # Calculate additional information about run. Update start time
    # ignoring empty bins at beginning.
    nBins=len(timeSeries)
    startTime=getTime(fullList[0],binWidth,startTime)

    # Define lower bound of noise to use for pulse finding
    endTime=getTime(nBins-1,binWidth,startTime)

    # Print all run information
    print "\nRun information:"
    if len(fullList)<w.shape[1]:
        emptyTime=(w.shape[1]-len(fullList))*binWidth
        print "\tIgnoring "+str(emptyTime)+" s of empty data."
    print "\tTelescope: ", telescope
    print "\tnBins: ", nBins
    print "\tResolution: ", binWidth, "s"
    print "\tStart time: ", startTime.iso
    print "\tEnd time: ", endTime.iso
    print "\tLower bound: ", threshold, 'sigma \n'
    print "Looking for "+str(nPulses)+" brightest giant pulses."
    print "Pulses: \n"

    # Find pulses with 10000 or fewer bins. Further resolve pulses if
    # finer binning is present.
    if w.shape[1]>10000:
        w_rebin=rebin(w,nBins=10000)
        timeSeries_rebin=getTimeSeries(w_rebin,nNoiseBins)
        pulseList=getPulses(timeSeries_rebin,threshold=threshold)
        pulseList=[(resolvePulse(
                    timeSeries,pos*w.shape[1]/10000,binWidth=binWidth,
                    searchRadius=1.0/10000),height) 
                   for (pos,height) in pulseList]

    else:
        pulseList=getPulses(timeSeries,threshold=threshold)

    # Find pulses
    #pulseList=getPulses(timeSeries,threshold=threshold)

    # Define the axis maximum to use for plotting
    ymax=1.3*np.amax(timeSeries) 

    # Loop through pulses to disply output
    for j in xrange(nPulses):

        # Check if there are still found pulses
        try:
            currentMaxLoc,currentMax=pulseList[j]
        except IndexError:
            print "Not enough giant pulses found!\n"
            break

        # Get pulse time 
        pulseTime=getTime(currentMaxLoc,binWidth,startTime)

        # Output pulse information
        print str(j+1)+'.\tIndex = '+str(currentMaxLoc)      
        print '\tTime = '+str(pulseTime.iso)
        print '\tTime (mjd) = '+str(pulseTime.mjd)
        print '\tPeak pulse height = '+str(round(currentMax,1))+' sigma\n'
        
        ###  Begin plotting ###
        plt.figure()
        indexList=getPeriod(currentMaxLoc,binWidth,nBins-1)
        timeList=indexList
        
        # Plot pulse profile
        plt.plot(timeList,timeSeries[indexList],'k-',label='Pulse Profile')

        # Plot legend, write labels, set limits, and show plot
        plt.ylim(0,ymax)
        plt.xlim(min(timeList),max(timeList))
        plt.xlabel('Time')
        plt.ylabel('Intensity')
        plt.show()
        ### End plotting ###
