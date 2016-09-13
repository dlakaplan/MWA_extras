import numpy
import casac
import sys,os
from matplotlib import pyplot as plt

tb=casac.casac.table()

for filename in sys.argv[1:]:
    tb.open(filename)
    G=tb.getcol('CPARAM')
    f=plt.figure()
    ax=f.add_subplot(1,1,1)
    antennae=numpy.arange(G.shape[2])
    h=ax.plot(antennae,numpy.abs(G[0]).mean(axis=0),
              antennae,numpy.abs(G[1]).mean(axis=0),'--')
    ax.legend(h,('X','Y'))
    plt.xlabel('Antenna Number')
    plt.ylabel('Gain Mean Absolute Value')
    plt.savefig(os.path.splitext(filename)[0] + '.png')
    print "Wrote %s" % (os.path.splitext(filename)[0] + '.png')
