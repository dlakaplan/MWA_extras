import os,sys,re
from optparse import OptionParser
import aplpy,numpy
import pylab,os
from astropy import wcs
from astropy.io import fits as pyfits


Istring='bm1.0_I'
Vstring='bp1.0_V'


######################################################################

def main():
    usage="Usage: %prog [options] <images>\n"

    parser = OptionParser(usage=usage)
    parser.add_option('-o','--out',default='mosaic',
                      dest='outroot',
                      help='Root name for output')
    parser.add_option('-t','--threshold',dest='threshold',default=0.05,
                      type='float',
                      help='Beam threshold [default=%default]')
    parser.add_option('-v','--verbose',action="store_true",
                      dest="verbose",default=False,
                      help="Increase verbosity of output")
    (options, args) = parser.parse_args()
    images=args

    newbeam=[re.sub('_regrid(\S+).fits','_beamI_regrid\\1.fits',i) for i in images]


    ff0=pyfits.open('MWATS_201404_I_MOL.fits')
    mask=ff0[1].data<=2

    f=pyfits.open(images[0])
    D=numpy.zeros(f[0].data.shape)
    W=numpy.zeros(f[0].data.shape)
    N=numpy.zeros(f[0].data.shape)

    pylab.clf()
    f=pylab.gcf()
    f.set_figwidth(12)
    f.set_figheight(6)

    for i in xrange(len(images)):
        #imagelist=' '.join(images[:i+1])
        #command='python /localhome/kaplan/python/mosaic_images.py -i %s -t %f -o %s_%03d.fits %s' % (
        #images[0],
        #    options.threshold,
        #    options.outroot,
        #    i,
        #    imagelist)
        #print command
        fi=pyfits.open(images[i])

        try:
            fb=pyfits.open(newbeam[i])
        except:
            #newbeam[i]=newbeam[i].replace('-v_','-i_')
            newbeam[i]=newbeam[i].replace(Vstring,Istring)
            fb=pyfits.open(newbeam[i])
        beam=(fb[0].data)
        beam[0,0][numpy.isnan(fi[0].data[0,0])]=0
        beam[0,0][beam[0,0]<=options.threshold*beam[0,0].max()]=0
        fi[0].data[0,0][numpy.isnan(fi[0].data[0,0])]=0
        D[0,0]+=fi[0].data[0,0]*beam[0,0]**2
        W[0,0]+=beam[0,0]**2
        N[0,0]+=beam[0,0]

        DD=D/W
        DD[0,0][W[0,0]==0]=numpy.nan
        DD[mask]=numpy.nan
        fi[0].data=DD
        if os.path.exists('temp.fits'):
            os.remove('temp.fits')
        fi.writeto('temp.fits')
        w=wcs.WCS(fi[0].header,naxis=(1,2))

        f.clf()

        gc=aplpy.FITSFigure('temp.fits',figure=f)
        gc.show_colorscale(vmin=-0.1,vmax=0.1,cmap=pylab.cm.copper)
        gc.recenter(14*15,-22.5, width=160, height=95)
        gc.axis_labels.hide_x()
        gc.axis_labels.hide_y()
        gc.tick_labels.hide()
        gc.add_label(0.1,0.9,'%03d' % i, relative=True,
                     size=16)
        ra=numpy.linspace(8,20,50)*15
        for dec in xrange(-80,30,10):
            x,y=w.wcs_world2pix(ra,ra*0+dec,0)
            pylab.plot(x,y,'k')
            if dec < 20 and dec>-80:
                gc.add_label(15*16,dec,'$\\delta=%d^\\ocirc$' % dec,
                             size=12,verticalalignment='bottom')
        dec=numpy.linspace(-80,20,50)
        for ra in xrange(8,22,2):
            ra_show=ra
            if ra_show < 0:
                ra_show+=24
            x,y=w.wcs_world2pix(0*dec+ra_show*15,dec,0)
            pylab.plot(x,y,'k')
            gc.add_label(ra*15,-45,'$\\alpha=%d^h$' % ra_show,
                         size=12,verticalalignment='bottom',
                         horizontalalignment='left')
        ax=pylab.gca()
        ax.axis('off')
        pylab.box('off')
        pylab.savefig('MWATS_201404_I_MOL_%03d.png' % i)
        print 'MWATS_201404_I_MOL_%03d.png' % i

    if os.path.exists('temp.fits'):
        os.remove('temp.fits')


######################################################################

if __name__=="__main__":
    main()
