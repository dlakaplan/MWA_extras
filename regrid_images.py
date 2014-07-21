"""
Usage:

/usr/local/casapy-stable-42.0.26465-001-64b/casapy --nogui -c ~/python/regrid_images.py -x 4 -y 2 -i 1063373456-image-i.fits 1063373456-image-i.fits 1063373816-image-i.fits 1063374176-image-i.fits 1063374536-image-i.fits 1063374896-image-i.fits
"""

import os,numpy,pyfits,sys,inspect,shutil,time,math
from optparse import OptionParser
try:
    from mwapy import get_observation_info
except:
    pass

Istring='bm1.0_I'
Vstring='bp1.0_V'

######################################################################

def regrid(infits, csys, shape, projection):
    root,ext=os.path.splitext(infits)
    outfile='%s_regrid%s.image' % (root,projection)
    outfits='%s_regrid%s.fits' % (root,projection)

    if os.path.exists('%s.image' % root):
        #os.system('rm -r %s.image' % root)
        shutil.rmtree('%s.image' % root)
    ia.fromfits(root + '.image',infits)
    if os.path.exists(outfile):
        #os.system('rm -r %s' % outfile)    
        shutil.rmtree('%s' % outfile)    
    ia.regrid(outfile=outfile,csys=csys.torecord(),
              shape=shape,
              axes=[0,1])
    ia.close()
    ia.open(outfile)
    if os.path.exists(outfits):
        os.remove(outfits)
    ia.tofits(outfits)
    ia.close()
    #os.system('rm -r %s' % outfile)
    shutil.rmtree('%s' % outfile, ignore_errors=True)
    if os.path.exists(outfile):
        shutil.rmtree('%s' % outfile, ignore_errors=True)        
    #os.system('rm -r %s.image' % root)
    shutil.rmtree('%s.image' % root)
    return outfits


######################################################################

def main():
    usage="Usage: %prog [options] <images>\n"
    usage+='\tRegrids all images onto a common coordinate frame\n'
    usage+='\tFrame is that of template image (if specified), otherwise first image\n'
    usage+='\tImage size is expanded by a given factor if desired\n'
    parser = OptionParser(usage=usage)
    parser.add_option('-i','--image',dest='image',default=None,
                      help='Template image base')
    parser.add_option('-x','--expandx',dest='factorx',default=2,
                      type='int',
                      help='X Expansion factor [default=%default]')
    parser.add_option('-y','--expandy',dest='factory',default=2,
                      type='int',
                      help='Y Expansion factor [default=%default]')
    parser.add_option('--projection',dest='projection',type='choice',
                      default='ZEA',
                      choices=['ZEA','SIN','AIT','HPX','MOL'],
                      help='Image projection [default=%default]')
    parser.add_option('--deczero',dest='deczero',default=False,
                      action='store_true',
                      help='Set reference Dec to 0 for AIT and MOL?')
    parser.add_option('--nobeam',dest='beam',default=False,
                      action='store_true',
                      help='Do not compute primary beams')
    parser.add_option('-v','--verbose',action="store_true",
                      dest="verbose",default=False,
                      help="Increase verbosity of output")

    db=None
    execname=inspect.stack()[0][1]
    imin=1
    while sys.argv[imin] != execname:
        imin+=1
    (options, args) = parser.parse_args(args=sys.argv[imin+1:])
    images=args
    dobeam=not options.beam
    if options.image is None:
        image0=images[0]
    else:
        image0=options.image
    template='_template.fits'
    factorx=options.factorx
    factory=options.factory

    f=pyfits.open(image0)
    shape0=f[0].data.shape
    f[0].data=numpy.zeros((shape0[0],shape0[1],
                           factory*shape0[2],factorx*shape0[3]))
    f[0].header['CRPIX1']=shape0[3]/2*factorx+1
    f[0].header['CRPIX2']=shape0[2]/2*factory+1

    if os.path.exists(template):
        os.remove(template)
    f.writeto(template)
    print 'Created template (%dx%d pixels) from %s' % (factorx*shape0[3],
                                                       factory*shape0[2],
                                                       image0)

    templateroot,ext=os.path.splitext(template)
    if os.path.exists('%s.image' % templateroot):
        os.system('rm -r %s.image' % templateroot)
    ia.fromfits(templateroot + '.image',template)
    csys=ia.coordsys()
    if not options.projection==csys.projection()['type']:
        csys.setprojection(options.projection)
        if options.projection in ['MOL','AIT'] and options.deczero:
            # set the Dec reference to 0
            ddec=csys.increment()['numeric'][1]*180/math.pi
            ref=csys.referencevalue(format='n')
            dec0=ref['numeric'][1]*180/math.pi
            ref['numeric'][1]=0            
            csys.setreferencevalue(ref)
            refpix=csys.referencepixel()
            refpix['numeric'][1]+=int(math.fabs(dec0/ddec))
            csys.setreferencepixel(value=refpix['numeric'])

        if os.path.exists('%s_%s.image' % (templateroot,options.projection)):
            os.system('rm -r %s_%s.image' % (templateroot,options.projection))
        im2=ia.regrid(templateroot + '_' + options.projection + '.image',
                      csys=csys.torecord(),
                      overwrite=True)    
        ia.close()
        ia.open(templateroot + '_' + options.projection + '.image')
        if os.path.exists(templateroot + '_' + options.projection + '.fits'):
            os.remove(templateroot + '_' + options.projection + '.fits')

            print 'Created %s template %s' % (options.projection,
                                              templateroot + '_' + options.projection + '.fits')
            ia.tofits(templateroot + '_' + options.projection + '.fits')
    shape=ia.shape()
    ia.close()

    new=[]
    newbeam=[]
    for im in images:
        if options.verbose:
            print 'Working on %s...' % im
        fi=pyfits.open(im)
        try:
            obsid=fi[0].header['GPSTIME']
            delays=fi[0].header['DELAYS']
            delays=map(int,delays.split(','))
        except:
            if db is None:
                from mwapy.obssched.base import schedule
                db = schedule.getdb()
            d,f=os.path.split(im)
            obsid=int(f.split('-')[0])            
            info=get_observation_info.MWA_Observation(obsid,db=db)
            delays=info.delays

        if dobeam:
            f,e=os.path.splitext(im)
            if Vstring in f and not (os.path.exists(f+'_beamXX'+e) or os.path.exists(f+'_beamI'+e)):
                f=f.replace(Vstring,Istring)
            if not os.path.exists(f+'_beamI'+e):
                if not os.path.exists(f+'_beamXX'+e):
                    # need to create the beam
                    from mwapy.pb import make_beam
                    beamfiles=make_beam.make_beam(im, delays=delays)
                else:
                    print 'XX beam %s already exists' % (f+'_beamXX'+e)
                    beamfiles=[f+'_beamXX'+e,
                               f+'_beamYY'+e]
                if not os.path.exists(f+'_beamI'+e):
                    fx=pyfits.open(beamfiles[0])
                    fy=pyfits.open(beamfiles[1])
                    fx[0].data=0.5*(fx[0].data+fy[0].data)
                    beamfileI=beamfiles[0].replace('beamXX','beamI')
                    if os.path.exists(beamfileI):
                        os.remove(beamfileI)
                    fx.writeto(beamfileI)
                    print 'Made beam for %s: %s' % (im,beamfileI)
            else:
                print 'I beam %s already exists' % (f+'_beamI'+e)
                beamfileI=f+'_beamI'+e
        try:
            outfits=regrid(im, csys, shape, options.projection)
            print 'Regridded %s -> %s' % (im,outfits)
            new.append(outfits)
            if dobeam:
                newbeam.append(regrid(beamfileI, csys, shape, 
                                      options.projection))
                print 'Regridded %s -> %s' % (beamfileI,newbeam[-1])
        except IndexError:
            print 'Failed to regrid %s' % im
    #os.remove(template)
    #os.system('rm -r %s.image' % templateroot)

######################################################################

if __name__=="__main__":
    main()
