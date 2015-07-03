import sys

def getFileName(telName):
    if not len(telName)==4:
        print "Invalid telescope name. Use format 'C01L' or 'E02R'"
        return ''
    if not telname.lower()[0] in ['c','e','w','s']:
        print "Invalid telescope name. Use format 'C01L' or 'E02R'"
        return ''
    else:
        region=telname.lower()[0]
    try:
        index=int(telname[1:3])
    except ValueError:
        print "Invalid telescope name. Use format 'C01L' or 'E02R'"
        return ''
    if not telname.lower()[-1] in ['l','r']:
        print "Invalid polarization. Use L or R."
        return ''
    else:
        pol=telname.lower()[-1]
        
    if pol=='l':
        sampler=32
    else:
        sampler=0

    if region=='c':
        if not(0<=index<=14) or index==7:
            print "Invalid telescope index. For C, 0<=index<=14 and index~=7."
            return ''
        if index>7:
            index=index-1
    if region=='w':
        if not 1<=index<=6:
            print "Invalid telescope index. For W, 1<=index<=6."
            return ''
        sampler+=13
    if region=='e':
        if not 2<=index<=6:
            print "Invalid telescope index. For E, 2<=index<=6."
            return ''
        sampler+=18
    if region=='s':
        if not(1<=index<=6) or index==5:
            print "Invalid telescope index. For C, 1<=index<=6 and index ~=5."
            return ''
        if index==6:
            index=5
        sampler+=24
    sampler+=index
    node=sampler/4
    rawvolt=sampler%4+1
    citanode=node/2+1
    if node%2==0:
        disk='a'
    else:
        disk='b'
    filename='raw_voltage'+str(rawvolt)+'.B0531+21_1.node'+str(node)+'.scan0'
    print "Nominal path: cita"+str(citanode)+':/mnt/'+disk+'/gsbuser/'
    print "File name: "+filename
telname=sys.argv[1]
getFileName(telname)
