import sys
from pulsarAnalysis.Misc.GMRTNaming import getNodeVolt
delayDict={
    'C00R':None, 'C00L':94,
    'C01R':None, 'C01L':39,
    'C02R':None, 'C02L':0,
    'C03R':None, 'C03L':-14,
    'C04R':0, 'C04L':0,
    'C05R':36, 'C05L':36,
    'C06R':46, 'C06L':46,
    'C08R':78, 'C08L':78,
    'C09R':23, 'C09L':23,
    'C10R':58, 'C10L':58,
    'C11R':25, 'C11L':25,
    'C12R':None, 'C12L':109,
    'C13R':94, 'C13L':94,
    'C14R':43, 'C14L':43,
    'W01R':26, 'W01L':26,
    'W02R':121, 'W02L':121,
    'W03R':None, 'W03L':None,
    'W04R':None, 'W04L':480,
    'W05R':743, 'W05L':743,
    'W06R':1028, 'W06L':1028,
    'E02R':374, 'E02L':375,
    'E03R':664, 'E03L':665,
    'E04R':1153, 'E04L':1153,
    'E05R':1643, 'E05L':None,
    'E06R':1979, 'E06L':1979,
    'S01R':419, 'S01L':420,
    'S02R':605, 'S02L':605,
    'S03R':935, 'S03L':None,
    'S04R':1285, 'S04L':1285,
    'S06R':1824, 'S06L':1824
    }

phaseDict={
    'C00R':None, 'C00L':None,
    'C01R':None, 'C01L':None,
    'C02R':None, 'C02L':None,
    'C03R':None, 'C03L':None,
    'C04R':1.35324132657, 'C04L':None,
    'C05R':-1.9325800295, 'C05L':None,
    'C06R':-0.426757793185, 'C06L':None,
    'C08R':0.867799992303, 'C08L':None,
    'C09R':-0.261623438433, 'C09L':None,
    'C10R':0.002098952482, 'C10L':None,
    'C11R':1.7350886217, 'C11L':None,
    'C12R':None, 'C12L':None,
    'C13R':1.09505939177, 'C13L':None,
    'C14R':2.82525095758, 'C14L':None,
    'W01R':0.254046050469, 'W01L':None,
    'W02R':-2.56150186701, 'W02L':None,
    'W03R':None, 'W03L':None,
    'W04R':None, 'W04L':None,
    'W05R':-1.09183057758, 'W05L':None,
    'W06R':-0.0, 'W06L':None,
    'E02R':1.27555565589, 'E02L':None,
    'E03R':0.0816068036961, 'E03L':None,
    'E04R':-0.0735378475527, 'E04L':None,
    'E05R':2.74345654823, 'E05L':None,
    'E06R':-1.40434619073, 'E06L':None,
    'S01R':0.260075860994, 'S01L':None,
    'S02R':-0.898483592391, 'S02L':None,
    'S03R':-1.8038799855, 'S03L':None,
    'S04R':2.30091730097, 'S04L':None,
    'S06R':-2.0510389114, 'S06L':None
    }

delayDict.update({getNodeVolt(i):delayDict[i] for i in delayDict.keys()})
phaseDict.update({getNodeVolt(i):phaseDict[i] for i in phaseDict.keys()})
def getDelay(dish):
    return delayDict[dish]
def getPhase(dish):
    return phaseDict[dish]

if __name__ == '__main__':
    if len(sys.argv)==2:
        try:
            delay=getDelay(sys.argv[1])
            print delay, "bytes"
            sys.exit()
        except KeyError:
            print "Telescope name not recognized.\n"
    elif len(sys.argv)==3:
        try:
            node=int(sys.argv[1])
            voltStream=int(sys.argv[2])
            try:
                delay=getDelay((node,voltStream))
                print delay, 'bytes'
                sys.exit()
            except KeyError:
                print "Invalid node and voltage stream combination."
                print "Telescope does not exist.\n"
        except ValueError:
            print "Invalid node or voltage stream indices.\n"
    else:
        print "Invalid number of inputs provided.\n"
    print "Usage: "
    print "python %s telescopeName" % sys.argv[0]
    print "or"
    print "python %s node voltageStream\n" % sys.argv[0]
    print "eg.\n"
    print "$ python GMRTNaming.py C01L"
    print "\t37 bytes"
    print "$ python GMRTNaming.py 8 2"
    print "\t37 bytes"
