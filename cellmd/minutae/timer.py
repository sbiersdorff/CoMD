#!/usr/bin/python
# replicate for afv

import sys
import usubs

outfile = "tmp.txt"
timingfile = "obj/spu/spuForce.s.timing"
verbose = 1

input = usubs.input()

timingfile = input.usetc("  Timing file?", timingfile)
outfile = timingfile+'.annotated'
outfile = input.usetc("  output file?", outfile)
#verbose = input.useti("  verbosity level (does nothing right now)?",verbose);

try:
    timingfp = open(timingfile,"r")
except:
    print ''
    print 'unable to open timing file:',timingfile
    print ''
    sys.exit(1)


try:
    outfp = open(outfile,"w")
except:
    print ''
    print 'unable to open outputfile:',timingfile
    print ''
    sys.exit(1)

outlines = [];

timinglines=timingfp.readlines();
timingfp.close()

maxlen = 10;
for l in timinglines:
    if ( len(l) > maxlen): maxlen = len(l)

maxlen = maxlen + 5
outfmt = '\n----------------\n    %06d.%02d %s'

k = timinglines[0].split();
if ( k[0] !=  ".file"):
    k = timinglines[1].split();
    idx = 1
    if ( k[0] !=  ".file"):
        print 'this file was not compiled with -g added to your compile line'
        sys.exit(2)
else:
    idx = 2;
print
print
print 'k is ',k
print
print
#
# open the source file
#srcfile = k[2].split('"')[1]
srcfile = k[1].split('"')[1]
try:
    fsrc = open(srcfile,"r")
except:
    print ''
    print 'unable to open source file:',srcfile
    print ''
    sys.exit(3)
    
srcLines = fsrc.readlines();
fsrc.close()
srcLines.insert(0,'dummy line for lining up with output\n')

nowTiming=[]
icount = 0
nowLine = ''
firstLine=1
for l in timinglines:

    nowTiming += [l]
    tmpSplit = l.split()
    if ( tmpSplit[0] == '.loc'):
        nowTiming.pop()
    elif(tmpSplit[0].startswith('1')):
        continue
    elif(tmpSplit[0].startswith('0')):
        continue
    elif(tmpSplit[0] == '.file'):
        continue
    elif(tmpSplit[0] == '.string'):
        continue
    elif(tmpSplit[0] == '.ident'):
        continue
    elif(tmpSplit[0] == '.type'):
        continue
    elif(tmpSplit[0] == '.size'):
        continue
    elif (tmpSplit[0] == '####'):
            continue
    elif(tmpSplit[0].startswith('.loc')):
        nowTiming.pop()
    elif(tmpSplit[0].startswith('.LS_')):
        nowTiming.pop()
        if(firstLine==1):
            outfp.writelines(nowTiming)
            nowTiming = []
        tmp = tmpSplit[0].split('_')
        ifile = int(tmp[2].split('f')[1])
        iline = int(tmp[3].split('l')[1].split(':')[0])
        if (ifile == 1):
            nowLine = outfmt%(iline,icount,srcLines[iline])
        else:
            nowLine = ''
    elif(tmpSplit[0].startswith('.LE_')):
        nowTiming.pop()
        nowTiming[0] = '%s%s'%(nowLine,nowTiming[0])
        outfp.writelines(nowTiming)
        #for x in nowTiming:print x,
        nowTiming = []
        nowLine = ''
        icount = 0
        iline = 0
    elif(len(tmpSplit)>2):
        nowTiming.pop()
        #print tmpSplit
        if((tmpSplit[4] == 'count')):
            icount = int(tmpSplit[6])
            
if(len(nowTiming)>0): outfp.writelines(nowTiming)
outfp.close();




