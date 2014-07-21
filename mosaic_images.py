import os,numpy,re
from optparse import OptionParser
try:
    import astropy.io.fits as pyfits
except ImportError:
    print "Astropy not found...importing pyfits"
    import pyfits

Istring='bm1.0_I'
Vstring='bp1.0_V'

######################################################################

def main():
    usage="Usage: %prog [options] <images>\n"

    parser = OptionParser(usage=usage)
    parser.add_option('-i','--image',dest='image',default=None,
                      help='Template image base')
    parser.add_option('-o','--output',dest='out',default=None,
                      help='Output image name')
    parser.add_option('-t','--threshold',dest='threshold',default=0.05,
                      type='float',
                      help='Beam threshold [default=%default]')
    parser.add_option('--nobeam',dest='beam',default=False,
                      action='store_true',
                      help='Do not use primary beams')
    parser.add_option('-v','--verbose',action="store_true",
                      dest="verbose",default=False,
                      help="Increase verbosity of output")
    (options, args) = parser.parse_args()
    images=args
    dobeam=not options.beam
    if options.image is None:
        image0=images[0]
    else:
        image0=options.image

    if options.out is None:
        print 'Must supply output name'
        sys.exit(1)

    new=images
    if dobeam:
        #newbeam=[i.replace('_regrid.fits','_beamI_regrid.fits') for i in images]
        newbeam=[re.sub('_regrid(\S+).fits','_beamI_regrid\\1.fits',i) for i in images]

    f=pyfits.open(image0)
    D=numpy.zeros(f[0].data.shape)
    W=numpy.zeros(f[0].data.shape)
    N=numpy.zeros(f[0].data.shape)
    for i in xrange(len(new)):
        if options.verbose:
            print 'Reading %d/%d %s' % (i+1,len(new),new[i])
        im=new[i]    
        if dobeam:
            if options.verbose:
                print 'Reading beam %s' % (newbeam[i])

            try:
                fb=pyfits.open(newbeam[i])
            except:
                #newbeam[i]=newbeam[i].replace('-v_','-i_')
                newbeam[i]=newbeam[i].replace(Vstring,Istring)
                print 'Actually reading beam %s' % (newbeam[i])
                fb=pyfits.open(newbeam[i])
            beam=(fb[0].data)
        fi=pyfits.open(im)
        if dobeam:
            beam[0,0][numpy.isnan(fi[0].data[0,0])]=0
            beam[0,0][beam[0,0]<=options.threshold*beam[0,0].max()]=0
        fi[0].data[0,0][numpy.isnan(fi[0].data[0,0])]=0
        if dobeam:
            D[0,0]+=fi[0].data[0,0]*beam[0,0]**2
            W[0,0]+=beam[0,0]**2
            N[0,0]+=beam[0,0]
        else:
            D[0,0]+=fi[0].data[0,0]
            W[0,0]+=(fi[0].data[0,0]>0)
            N[0,0]+=W[0,0]

    D[0,0]/=W[0,0]
    D[0,0][W[0,0]==0]=numpy.nan
    f[0].data=D
    f[0].header['TEMPLATE']=image0
    f[0].header['BMTHRESH']=options.threshold
    f[0].header['NCOMBINE']=len(new)    
    f.append(pyfits.ImageHDU(N, name='WEIGHT'))
    if os.path.exists(options.out):
        os.remove(options.out)
    f.writeto(options.out)
######################################################################

if __name__=="__main__":
    main()
