import sys

telList=[
    ('C00R',(0,1)), ('C00L',(8,1)),
    ('C01R',(0,2)), ('C01L',(8,2)),
    ('C02R',(0,3)), ('C02L',(8,3)),
    ('C03R',(0,4)), ('C03L',(8,4)),
    ('C04R',(1,1)), ('C04L',(9,1)),
    ('C05R',(1,2)), ('C05L',(9,2)),
    ('C06R',(1,3)), ('C06L',(9,3)),
    ('C08R',(1,4)), ('C08L',(9,4)),
    ('C09R',(2,1)), ('C09L',(10,1)),
    ('C10R',(2,2)), ('C10L',(10,2)),
    ('C11R',(2,3)), ('C11L',(10,3)),
    ('C12R',(2,4)), ('C12L',(10,4)),
    ('C13R',(3,1)), ('C13L',(11,1)),
    ('C14R',(3,2)), ('C14L',(11,2)),
    ('W01R',(3,3)), ('W01L',(11,3)),
    ('W02R',(3,4)), ('W02L',(11,4)),
    ('W03R',(4,1)), ('W03L',(12,1)),
    ('W04R',(4,2)), ('W04L',(12,2)),
    ('W05R',(4,3)), ('W05L',(12,3)),
    ('W06R',(4,4)), ('W06L',(12,4)),
    ('E02R',(5,1)), ('E02L',(13,1)),
    ('E03R',(5,2)), ('E03L',(13,2)),
    ('E04R',(5,3)), ('E04L',(13,3)),
    ('E05R',(5,4)), ('E05L',(13,4)),
    ('E06R',(6,1)), ('E06L',(14,1)),
    ('S01R',(6,2)), ('S01L',(14,2)),
    ('S02R',(6,3)), ('S02L',(14,3)),
    ('S03R',(6,4)), ('S03L',(14,4)),
    ('S04R',(7,1)), ('S04L',(15,1)),
    ('S06R',(7,2)), ('S06L',(15,2))
    ]

statDict={i[0]:'g' for i in telList}

# Broken telescopes
statDict['C00R']='m'
statDict['C01R']='m'
statDict['C02R']='m'
statDict['C03R']='m'
statDict['C12R']='m'
statDict['E05L']='m'

statDict['W03R']='b'
statDict['W04R']='b'
statDict['W03L']='b'
statDict['S03L']='b'

statDict.update({i[1]:statDict[i[0]] for i in telList})

def getStatus(dish):
    return statDict[dish]

if __name__ == '__main__':
    if len(sys.argv)==2:
        try:
            status=getStatus(sys.argv[1])
            if status=='b':
                print "Broken"
            elif status=='m':
                print "Missing"
            else:
                print "Good"
            sys.exit()
        except KeyError:
            print "Telescope name not recognized.\n"
    elif len(sys.argv)==3:
        try:
            node=int(sys.argv[1])
            voltStream=int(sys.argv[2])
            try:
                status=getStatus((node,voltStream))
                if status=='b':
                    print "Broken"
                elif status=='m':
                    print "Missing"
                else:
                    print "Good"
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
    print "$ python GMRTNaming.py C01R"
    print "\tMissing"
    print "$ python GMRTNaming.py 0 2"
    print "\tMissing"
