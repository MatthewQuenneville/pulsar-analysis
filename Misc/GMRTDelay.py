import sys
from pulsarAnalysis.Misc.GMRTNaming import getNodeVolt
delayDict={
    'C00R':0, 'C00L':92,
    'C01R':0, 'C01L':37,
    'C02R':0, 'C02L':-2,
    'C03R':0, 'C03L':-16,
    'C04R':0, 'C04L':-2,
    'C05R':36, 'C05L':34,
    'C06R':46, 'C06L':44,
    'C08R':78, 'C08L':76,
    'C09R':23, 'C09L':21,
    'C10R':58, 'C10L':56,
    'C11R':25, 'C11L':23,
    'C12R':0, 'C12L':107,
    'C13R':94, 'C13L':92,
    'C14R':43, 'C14L':41,
    'W01R':26, 'W01L':24,
    'W02R':121, 'W02L':119,
    'W03R':0, 'W03L':0,
    'W04R':0, 'W04L':478,
    'W05R':743, 'W05L':741,
    'W06R':1028, 'W06L':1026,
    'E02R':374, 'E02L':373,
    'E03R':664, 'E03L':663,
    'E04R':1153, 'E04L':1151,
    'E05R':1643, 'E05L':0,
    'E06R':1979, 'E06L':1977,
    'S01R':419, 'S01L':418,
    'S02R':605, 'S02L':603,
    'S03R':935, 'S03L':0,
    'S04R':1285, 'S04L':1283,
    'S06R':1824, 'S06L':1822
    }
delayDict.update({getNodeVolt(i):delayDict[i] for i in delayDict.keys()})
def getDelay(dish):
    return delayDict[dish]

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
