from astropy.io import fits as pyfits
import sys,os
from optparse import OptionParser

suffix='.2d.fits'
origsuffix='_beamI.fits'
newsuffix='.weight.fits'

usage="Usage: %prog [options] <images>\n"
parser = OptionParser(usage=usage)
parser.add_option('-t','--threshold',dest='threshold',default=0.05,
                  type='float',
                  help='Beam threshold [default=%default]')
parser.add_option('--reduce',dest='reduce',default=False,
                  action='store_true',
                  help='Reduce dimensionality of input images?')
(options, args) = parser.parse_args()



for file in args:
    root,ext=os.path.splitext(file)
    beamfile=root + origsuffix

    if not os.path.exists(file):
        print 'Cannot find image file %s' % file
        sys.exit(0)

    if options.reduce:
        f=pyfits.open(file)
        f[0].data=f[0].data[0,0]
        if f[0].header['CRVAL1']<0:
            f[0].header['CRVAL1']+=360
        f.verify('fix')
        newfile=root + suffix
        if os.path.exists(newfile):
            os.remove(newfile)
        f.writeto(newfile)
        f.close()
        print '%s -> %s' % (file, newfile)
        newroot,newext=os.path.splitext(newfile)
    else:
        newroot, newext=root, ext

    weightfile=newroot + newsuffix
    if not os.path.exists(beamfile):
        print 'Beam file %s does not exist' % beamfile
        sys.exit(1)        
    f=pyfits.open(beamfile)
    f[0].data[f[0].data<options.threshold]=0
    f[0].data=f[0].data[0,0]**2
    if f[0].header['CRVAL1']<0:
        f[0].header['CRVAL1']+=360
    f[0].header['BMTHRESH']=(options.threshold,'Threshold for primary beam correction')
    f.verify('fix')
    if os.path.exists(weightfile):
        os.remove(weightfile)
    f.writeto(weightfile)
    f.close()
    print '%s -> %s' % (beamfile,weightfile)


 
