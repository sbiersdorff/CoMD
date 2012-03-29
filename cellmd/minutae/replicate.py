#!/usr/bin/python
# replicate for afv

import sys
import usubs

mx = 2
my = 2
mz = 1
outfile = "out.txt"
infile = "trimers1.start"

input = usubs.input()

infile = input.usetc("Input file?", infile)
outfile = input.usetc("Input file?", outfile)
mx = input.useti("ex?",mx);
my = input.useti("ey?",my);
mz = input.useti("ez?",mz);


fp = open(infile,"r")
lines=fp.readlines();
fp.close()

fp = open(outfile,"w")

firstwords = lines.pop(0).split()
natom = int(firstwords[0])
ntitle = int(firstwords[1])

comment = []
for i in range(ntitle):
    comment.append(lines.pop(0))
if ( ntitle == 0):
    comment.append('\n');
    lines.pop(0);
taxes = [ float(x) for x in lines.pop(0).split()]

fp.write("%5d%5d\n"%(natom*mx*my*mz,ntitle))

for i in range(ntitle):
    fp.write(comment[i])
if(ntitle == 0): fp.write('\n');

fp.write("%20.8f%20.8f%20.8f\n"%(float(mx)*taxes[0],float(my)*taxes[1],float(mz)*taxes[2]))
for l in lines:
    k = l.split()
    s = "%5d%20.8f%20.8f%20.8f\n"%(int(k[3]),float(k[4]),float(k[5]),float(k[6]))
    for iz in range(mz):
        dz = taxes[2]*float(iz)
        zn = float(k[2])+dz
        for iy in range(my):
            dy = float(iy)*taxes[1]
            yn = float(k[1])+dy
            for ix in range(mx):
                dx = float(ix)*taxes[0]
                xn = float(k[0])+dx
                fp.write("%20.8f%20.8f%20.8f%s"%(xn,yn,zn,s))

fp.close()






