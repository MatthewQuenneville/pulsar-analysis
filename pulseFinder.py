#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pylab as plt
import string
from astropy.time import Time,TimeDelta
import warnings

# Crab frequency 
# Should implement an ephemeris/polynomial based frequency, but this
# is just used for plotting, so approximate is fine
crabFreq=29.946923 

# Set number of pulses to find
nPulses=3

# Set noise threshold in number of standard deviations
threshold=5

# Set number of bins over which to normalize by noise. Use 0 for one
# noise bin per second
nNoiseBins=1

# Polarizations to sum over
pol_select = (0, 3)

def getTelescope(fileName):
    # Gets telescope name based on file name

    if 'jb' in string.lower(fileName):
        return 'Jodrell Bank'
    elif 'gmrt' in string.lower(fileName):
        return 'GMRT'
    else:
        print "Telescope not recognized."
        return 'Unknown'

def rms(sequence):
    # Gets root-mean square of sequence

    return np.sqrt(np.mean(np.power(sequence,2)))

def getStartTime(fileName):
    # Gets start time of trial based on file name

    dateString=fileName.split('foldspec')[1][1:24]
    return Time(dateString,format='isot',scale='utc')

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

def getTimeSeries(f,ic,nNoiseBins):
    # Gets profile time series over which to search for giant pulses

    if f.shape[-1] == 4:
        n = f[..., pol_select].sum(-1)[0,:,:]
    else:
        n = f[0,:,:]
        
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        n /= ic[0,:,:]

    n_median = nanMedian(n)
    nn = n / n_median[:,np.newaxis] - 1.

    timeSeries = nn.sum(0)
    timeSeries = timeSeries[~np.isnan(timeSeries)]
    noiseBins=[int(i) for i in np.linspace(0,len(timeSeries),nNoiseBins+1)]
    noise=[rms(timeSeries[noiseBins[i]:noiseBins[i+1]]) 
           for i in range(nNoiseBins)]

    for i in range(nNoiseBins):
        for j in range(noiseBins[i],noiseBins[i+1]):
            timeSeries[j]/=noise[i]

    return timeSeries

def getPulses(timeSeries,threshold=5):
    # Gets a list of all pulses higher than the noise threshold

    nBins=len(timeSeries)

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
    if len(sys.argv) < 2:
        print "Usage: %s foldspec icounts" % sys.argv[0]
        # Run the code as: ./pulseFinder.py data_foldspec.npy data_icounts.py.
        sys.exit(1)  

    # Folded spectrum axes: time, frequency, phase, pol=4 (XX, XY, YX, YY).
    f = np.load(sys.argv[1])
    ic = np.load(sys.argv[2] if len(sys.argv) == 3 else
                 sys.argv[1].replace('foldspec', 'icount'))
    
    # Get basic run attributes
    telescope=getTelescope(sys.argv[1])
    startTime=getStartTime(sys.argv[1])
    deltat=getDeltaT(sys.argv[1])

    # Set nNoiseBins to one per second if invalid value is given
    if nNoiseBins<1:
        nNoiseBins=int(deltat)

    # Find populated bins
    fullList=np.flatnonzero(ic.sum(0).sum(0))

    # Get timeseries to search for pulses
    timeSeries=getTimeSeries(f,ic,nNoiseBins)

    # Calculate additional information about run. Update start time
    # ignoring empty bins at beginning.
    binWidth=float(deltat)/f.shape[2]
    nBins=len(timeSeries)
    startTime=getTime(fullList[0],binWidth,startTime)

    # Define lower bound of noise to use for pulse finding
    endTime=getTime(nBins-1,binWidth,startTime)

    # Print all run information
    print "\nRun information:"
    if len(fullList)<len(ic[0,0,:]):
        emptyTime=(len(ic[0,0,:])-len(fullList))*binWidth
        print "\tIgnoring "+str(emptyTime)+" s of empty data."
    print "\tTelescope: ", telescope
    print "\tnBins: ", nBins
    print "\tResolution: ", binWidth, "s"
    print "\tStart time: ", startTime.iso
    print "\tEnd time: ", endTime.iso
    print "\tLower bound: ", threshold, 'sigma \n'
    print "Looking for "+str(nPulses)+" brightest giant pulses."
    print "Pulses: \n"

    # Find pulses
    pulseList=getPulses(timeSeries,threshold=threshold)

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

        # Plot line for peak intensity
        plt.plot(timeList,[currentMax for i in indexList],'r--',
                 label='Peak Intensity')
        
        # Plot pulse profile
        plt.plot(timeList,timeSeries[indexList],'k-',label='Pulse Profile')

        # Plot line for intensity threshold
        plt.plot(timeList,[threshold for i in indexList],'b--',
                 label='Intensity Threshold')

        # Plot legend, write labels, set limits, and show plot
        plt.legend()
        plt.ylim(0,ymax)
        plt.xlim(min(timeList),max(timeList))
        plt.xlabel('Time')
        plt.ylabel('Intensity')
        plt.show()
        ### End plotting ###
