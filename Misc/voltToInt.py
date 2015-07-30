import sys
import numpy as np
from pulsarAnalysis.Misc import GMRTNaming,GMRTDelay
import os
import pulsarAnalysis.GPs.pulseFinder as pf

def rotateVoltage(w,theta):
    return w*np.exp(1J*theta)

def getIntensity(w,keepdims=False):
    return np.abs(w)**2

def getDishes(filename):
    dishlist=[]
    for i,_ in GMRTNaming.telList:
        if i in filename:
            dishlist.append(i)
            
    return dishlist

def getVoltage(path):
    runInfo={}
    print "Opening file:", path
    if 'voltage' in path:
        w=np.load(path)
        w=np.swapaxes(w,0,1)
        deltat=pf.getDeltaT(path)
        telescope=pf.getTelescope(path)
        startTime=pf.getStartTime(path)
        binWidth=pf.getWaterfallBinWidth(telescope,w.shape[0])
    else:
        print "Error, the following file name is not recognized:"
        print path
        return None

    #isNotNan=~np.isnan(w.sum(0))
    #fullList=np.flatnonzero(isNotNan)
    #n=w[:,isNotNan]
    
    runInfo['dishList']=getDishes(path)
    runInfo['binWidth']=binWidth
    runInfo['telescope']=telescope
    runInfo['deltat']=deltat
    runInfo['startTime']=startTime
    #runInfo['fullList']=fullList
    return w, runInfo

def getSummedVoltages(pathList,rotation=True):
    if len(pathList)==0:
        print "Usage: %s foldspec" % sys.argv[0]
        # Run the code as: ./script.py data_foldspec.npy.
        sys.exit(1) 
    runInfo={}
    fulldishlist=[]
    goodFileFound=False
    print "Opening files..."
    for iPath in pathList:
        if os.path.isdir(iPath):
            iFileList=os.listdir(iPath)
        else:
            iFileList=[iPath]
        for jFile in iFileList:
            if os.path.isdir(iPath):
                iFile=os.path.join(iPath,jFile)
            else:
                iFile=jFile
            print "Reading",iFile 
            if 'voltage' in iFile:
                w=np.load(iFile)
                if not goodFileFound:
                    deltat=pf.getDeltaT(iFile)
                    telescope=pf.getTelescope(iFile)
                    startTime=pf.getStartTime(iFile)
                    binWidth=pf.getWaterfallBinWidth(telescope,w.shape[0])
                    outfileName=jFile.replace('voltage','waterfall')
            else:
                print "Error, the following file name is not recognized:"
                print iFile
                continue
            dishlist=getDishes(jFile)
            if len(dishlist)==1:
                fulldishlist=fulldishlist+dishlist
            else:
                print "Error, multiple dishes found in file:"
                print iFile
                continue

            phase=GMRTDelay.getPhase(dishlist[-1]) if rotation else 0.
            w=rotateVoltage(w,phase)
            if not goodFileFound:
                n=w
                outfileName=outfileName.replace(dishlist[0],'Phased')
                goodFileFound=True
            elif n.shape==w.shape:
                n=n+w
            else:
                print "Error, shape mismatch in file:"
                print iFile

    #isNotNan=~np.isnan(n.sum(0))
    #fullList=np.flatnonzero(isNotNan)
    #n=n[:,isNotNan]
    runInfo['binWidth']=binWidth
    runInfo['telescope']=telescope
    runInfo['deltat']=deltat
    runInfo['startTime']=startTime
    #runInfo['fullList']=fullList
    runInfo['dishList']=dishlist
    runInfo['outfileName']=outfileName
    return n, runInfo

if __name__ == "__main__":
    w,runInfo=getSummedVoltages(sys.argv[1:],rotation=True)
    I=getIntensity(w,keepdims=True)
    print "Saving to:"
    print runInfo['outfileName']
    np.save(runInfo['outfileName'],I)
