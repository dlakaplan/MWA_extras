import numpy
from astropy.table import Table

fi=open('/Users/dlk/Desktop/pystaratlas/data/constellationship.fab')
Constellations={}
for l in fi.readlines():
    d=l.split()
    name=d[0]
    n=int(d[1])
    data=numpy.array(map(int,d[2:]))
    Constellations[name]=[n,data]
fi.close()      

HIP=Table.read('/Users/dlk/Desktop/pystaratlas/data/hip2.dat',
               format='ascii')
ConstellationStars=[]
I=[]
for c in Constellations.keys():
    for i in xrange(0,len(Constellations[c][1])):
        star=HIP[HIP['HIP']==Constellations[c][1][i]]
        if not Constellations[c][1][i] in ConstellationStars:
            ConstellationStars.append(Constellations[c][1][i])
            I.append(numpy.where(HIP['HIP']==Constellations[c][1][i])[0][0])

HIP_touse=HIP[numpy.array(sorted(I))]
HIP_touse.write('HIP_constellations.dat',format='ascii.commented_header')

