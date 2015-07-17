import sys
import numpy as np
from pulsarAnalysis.Misc import GMRTNaming,GMRTDelay
import os
import pulsarAnalysis.GPs.pulseFinder as pf

def rotateVoltage(w,theta):
    w_new=np.zeros(w.shape)
    w_new[...,0]=w[...,0]*np.cos(theta)-w[...,1]*np.sin(theta)
    w_new[...,1]=w[...,0]*np.sin(theta)+w[...,1]*np.cos(theta)
    return w_new

def getIntensity(w,keepdims=False):
    w_new=np.zeros(w.shape[:-1])
    w_new=w[...,0]**2+w[...,1]**2
    if keepdims:
        w_new=w_new[...,np.newaxis]
    return w_new

def getDish(filename):
    dishName=''
    for i,_ in GMRTNaming.telList:
        if i in filename and dishName=='':
            dishName=i
        elif i in filename:
            dishName=None
    if dishName==None:
        print "Error, multiple dish names found for file:"
        print filename
        return ''
    elif dishName=='':
        print "Error, no dish name found for file:"
        print filename
        return ''
    else:
        return dishName    


def getSummedVoltages(pathList,rotation=True):
    if len(pathList)==0:
        print "Usage: %s foldspec" % sys.argv[0]
        # Run the code as: ./script.py data_foldspec.npy.
        sys.exit(1) 
    runInfo={}
    dishList=[]
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
            if 'foldspec' in iFile:
                f=np.load(iFile)
                ic=np.load(iFile.replace('foldspec', 'icount'))
                f=f.sum(0)
                ic=ic.sum(0)
                if f.shape[-1]==2:
                    w=f/ic[...,np.newaxis]
                    w=np.swapaxes(w,0,1)
                    if not goodFileFound:
                        deltat=pf.getDeltaT(iFile)
                        telescope=pf.getTelescope(iFile)
                        startTime=pf.getStartTime(iFile)
                        binWidth=deltat/f.shape[1]
                        outfileName=jFile.replace('foldspec','waterfall')
                        outFileName=outFileName.replace('VoltageDump','')
                else:
                    print "File does not contain exactly 2 time series:"
                    print iFile
                    print "Skipping...\n"
                    continue

            elif 'waterfall' in iFile:
                w=np.load(iFile)
                if not w.shape[-1]==2:
                    print "File does not contain exactly 2 time series:"
                    print iFile
                    print "Skipping...\n"
                    continue
                if not goodFileFound:
                    deltat=pf.getDeltaT(iFile)
                    telescope=pf.getTelescope(iFile)
                    startTime=pf.getStartTime(iFile)
                    binWidth=pf.getWaterfallBinWidth(telescope,w.shape[0])
                    outfileName=jFile.replace('VoltageDump','')
            else:
                print "Error, the following file name is not recognized:"
                print iFile
                continue
            dishName=getDish(jFile)
            dishList.append(dishName)
            phase=GMRTDelay.getPhase(dishName) if rotation else 0.
            w=rotateVoltage(w,phase)
            if not goodFileFound:
                n=w
                outfileName=outfileName.replace(dishName,'VoltSum')
                goodFileFound=True
            elif n.shape==w.shape:
                n=n+w
            else:
                print "Error, shape mismatch in file:"
                print iFile

    isNotNan=~np.isnan(n.sum(-1).sum(0))
    fullList=np.flatnonzero(isNotNan)
    n=n[:,isNotNan,:]
    runInfo['binWidth']=binWidth
    runInfo['telescope']=telescope
    runInfo['deltat']=deltat
    runInfo['startTime']=startTime
    runInfo['fullList']=fullList
    runInfo['dishList']=dishList
    runInfo['outfileName']=outfileName
    return n, runInfo

if __name__ == "__main__":
    w,runInfo=getSummedVoltages(sys.argv[1:])
    I=getIntensity(w,keepdims=True)
    print "Saving to:"
    print runInfo['outfileName']
    np.save(runInfo['outfileName'],I)
