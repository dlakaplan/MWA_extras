#!/usr/bin/env python

import logging, sys, os, glob, string, re, urllib, math, time
from optparse import OptionParser
import numpy
import astropy.io.fits as pyfits
from astropy import units as u
from astropy.coordinates import Angle
import time, datetime
import subprocess


import mwapy
from mwapy import get_observation_info, make_metafiles
from mwapy.pb import primary_beam,make_beam

                    
# configure the logging
logging.basicConfig(format='# %(levelname)s:%(name)s: %(message)s')
logger=logging.getLogger('mosaic_driftscan')
logger.setLevel(logging.WARNING)

try:
    from mwapy.obssched.base import schedule
except:
    logger.error("Unable to open connection to database")
    sys.exit(1)
# open up database connection
try:
    db = schedule.getdb()
except:
    logger.error("Unable to open connection to database")
    sys.exit(1)
                    
suffix='.2d.fits'
newsuffix='.weight.fits'
swarpconfig='mosaic.swarp'

######################################################################
class MWAImage():

    ##################################################
    def __init__(self, obsid=None):
        self.obsid=obsid

        self.ximage=None
        self.yimage=None
        self.iimage=None
        self.qimage=None
        self.uimage=None
        self.vimage=None                        
        self.metafits=None
        self.xbeam=None
        self.ybeam=None
        self.ibeam=None
        self.weight=None

######################################################################
def copy_metafits(image, metafits, keys=[]):
    """
    result=copy_metafits(image, metafits, keys=[])
    """
    fm=pyfits.open(metafits)
    try:
        fi=pyfits.open(image,mode='update')
    except IOError, err:
        print 'Unable to open FITS file %s\n\t%s' % (image, err)
        return False
    try:
        for k in fm[0].header.keys():
            if len(keys)==0 or k in keys:
                fi[0].header[k]=fm[0].header[k]
        fi.flush()
        return True
    except:
        logger.error('Error updating FITS header of %s' % image)
        return False
    

######################################################################
def copy_driftscanbeam(template, image,
                       XX='_beamXX',
                       YY='_beamYY'):
    """
    beamXX, beamYY=copy_driftscanbeam(template, image,
    XX='_beamXX',YY='_beamYY')
    """
    base,ext=os.path.splitext(template)
    beamXX=base + XX + ext
    beamYY=base + YY + ext
    if not os.path.exists(beamXX):
        logger.error('XX beam %s does not exist for %s' % (beamXX, template))
        return None
    if not os.path.exists(beamYY):
        logger.error('YY beam %s does not exist for %s' % (beamYY, template))
        return None
    fX=pyfits.open(beamXX)
    fY=pyfits.open(beamYY)
    if not os.path.exists(image):
        logger.error('Image %s does not exist' % image)
        return None
    f=pyfits.open(image)
    if not f[0].header['NAXIS1'] == fX[0].header['NAXIS1'] and f[0].header['NAXIS2'] == fX[0].header['NAXIS2']:
        logger.error('Image sizes do not match: template is %dx%d, while input image is %dx%d' % (
            fX[0].header['NAXIS1'],fX[0].header['NAXIS2'],
            f[0].header['NAXIS1'],f[0].header['NAXIS2']))
        return None
    delays0=numpy.array(map(int,fX[0].header['DELAYS'].split(',')))
    delays=numpy.array(map(int,f[0].header['DELAYS'].split(',')))
    if not (delays==delays0).all():
        logger.error('Delays do not match:\n\ttemplate is %s\n\tinput image is %s' % (delays0,delays))
        return None
    for k in f[0].header.keys():
        if k=='HISTORY':
            continue
        try:
            fX[0].header[k]=f[0].header[k]
        except:
            pass
        try:
            fY[0].header[k]=f[0].header[k]
        except:
            pass
    base,ext=os.path.splitext(image)
    beamXXout=base + XX + ext
    beamYYout=base + YY + ext
    fX[0].header['ORIGNAME']=beamXX
    if os.path.exists(beamXXout):
        os.remove(beamXXout)
    fX.writeto(beamXXout)
    fY[0].header['ORIGNAME']=beamYY
    if os.path.exists(beamYYout):
        os.remove(beamYYout)
    fY.writeto(beamYYout)
    logger.info('Wrote %s, %s' % (beamXXout, beamYYout))
    fX.close()
    fY.close()
    return beamXXout, beamYYout


######################################################################
def circular_mean(x):
    """
    http://en.wikipedia.org/wiki/Mean_of_circular_quantities
    x should be in radians or astropy.Angle etc.
    value returned will be in same units
    """
    c=numpy.cos(x)
    s=numpy.sin(x)
    return numpy.arctan2(s.mean(), c.mean())

######################################################################        

def main():


    localswarpconfig=os.path.join(os.path.split(__file__)[0], swarpconfig)
    
    usage="Usage: %prog [options] <files>\n"
    usage+="\tMosaic together MWA driftscan data\n"
    usage+="\tImages should be either raw XX,YY or primary-beam corrected I (or Q or U or V)\n"
    usage+="\tXX,YY data will be converted into Stokes I\n"
    usage+="\tStokes I beams will be computed for optimal weighting\n"
    usage+="\tswarp is used for mosaicing\n"
    parser = OptionParser(usage=usage,version=mwapy.__version__ + ' ' + mwapy.__date__)
    parser.add_option('--metafits',action="store_true",dest="metafits",default=False,
                      help="Force creation of metafits files?")
    parser.add_option('--copybeam',action="store_true",dest="copy_beam",default=False,
                      help="Copy the beam for the driftscan instead of creating separately?")
    parser.add_option('--analytic',action="store_true",dest="analytic_model",default=False,
                      help="Use the old analytic dipole model, instead of the default Sutinjo 2014 model.")
    parser.add_option('--noprecess',action='store_false',
                      dest='precess',default=True,
                      help='Do not precess coordinates to current epoch (faster but less accurate) [default=False]')
    parser.add_option('-t','--threshold',dest='threshold',default=0.05,
                      type='float',
                      help='Beam threshold [default=%default]')
    parser.add_option('--reduce',dest='reduce',default=False,
                      action='store_true',
                      help='Reduce dimensionality of input images?')
    parser.add_option('--oversample',dest='oversample', default=2,
                      type='int',
                      help='Swarp oversampling [default=%default]')
    parser.add_option('--ra',dest='ra', default=None,
                      help='Center RA of mosaic in hours or "middle" [default=auto]')
    parser.add_option('--dec',dest='dec', default=None,
                      help='Center Dec of mosaic in degrees [default=auto]')
    parser.add_option('--proj',dest='projection',default='MOL',
                      type='choice',
                      choices=['MOL','ZEA','SIN','AIT'],
                      help='Image projection [default=%default]')    
    parser.add_option('--arguments',dest='swarpcommands',default='',
                      help='Additional arguments for Swarp')
    parser.add_option('--out',dest='root', default='mosaic',
                      help='SWARP output root [default=%default]')
    parser.add_option('--swarp',dest='swarp',default=None,
                      help='Path to Swarp executable [default=search path]')
    parser.add_option('--config',dest='config',default=None,
                      help='Path to separate Swarp config file [default=%s]' % localswarpconfig)

    parser.add_option('--clean',action="store_true",dest="clean",default=False,
                      help="Clean temporary files?")
    parser.add_option('-v','--verbose',action="store_true",dest="verbose",default=False,
                      help="Increase verbosity of output")
    (options, args) = parser.parse_args()

        
    if (options.verbose):
        logger.setLevel(logging.INFO)

    if options.config is not None:
        localswarpconfig=options.config
    if not os.path.exists(localswarpconfig):
        logger.error('Cannot find swarp config file %s' % localswarpconfig)
        sys.exit(1)
    logger.info('Will use config file %s...' % localswarpconfig)

    inputfiles=args
    Imagelist={}

    if options.swarp is None:
        # python 2.6 compatibility
        try:
            swarp=subprocess.check_output(['which','swarp']).strip()
        except AttributeError:
            swarp=subprocess.Popen(['which', 'swarp'],
                                   stdout=subprocess.PIPE).communicate()[0].strip()                              
    else:
        swarp=options.swarp

    # go through the images
    # figure out which is for which polarization
    # and see if there is also a metafits file
    times={}
    times['start']=time.time()
    for file in inputfiles:
        polaxis=None
        results=file.split('_')
        obsid=int(results[0])

        if not Imagelist.has_key(obsid):
            Imagelist[obsid]=MWAImage(obsid)
        mwaimage=Imagelist[obsid]            
        currentstokes=None
        try:
            f=pyfits.open(file)
            if f[0].header['CTYPE3']=='STOKES':
                polaxis=3
            if f[0].header['CTYPE4']=='STOKES':
                polaxis=4
            if not polaxis is None:
                if f[0].header['CRVAL%d' % polaxis]==-5:
                    mwaimage.ximage=file
                    currentstokes='XX'
                elif f[0].header['CRVAL%d' % polaxis]==-6:
                    mwaimage.yimage=file
                    currentstokes='YY'
                elif f[0].header['CRVAL%d' % polaxis]>=1 and f[0].header['CRVAL%d' % polaxis]<=4:
                    mwaimage.iimage=file
                    currentstokes='I'
                else:
                    logger.warning('Do not know how to process Stokes=%d for %s; skipping...' % (
                        f[0].header['CRVAL%d' % polaxis], file))
                    continue
            else:
                if 'XX' in file:
                    mwaimage.ximage=file
                    currentstokes='XX'
                elif 'YY' in file:
                    mwaimage.yimage=file                
                    currentstokes='YY'
                
            # see if beams already exist
            # would be identified through header keywords
            if 'BEAM' in f[0].header.keys():
                if currentstokes=='I' and os.path.exists(f[0].header['BEAM']):
                    mwaimage.ibeam=f[0].header['BEAM']
                elif currentstokes=='XX' and os.path.exists(f[0].header['BEAM']):
                    mwaimage.xbeam=f[0].header['BEAM']
                elif currentstokes=='YY' and os.path.exists(f[0].header['BEAM']):
                    mwaimage.ybeam=f[0].header['BEAM']
                    
        except:
            if 'XX' in file:
                mwaimage.ximage=file
                currentstokes='XX'
            elif 'YY' in file:
                mwaimage.yimage=file
                currentstokes='YY'
                
        if os.path.exists('%d.metafits' % obsid):
            mwaimage.metafits='%d.metafits' % obsid

    logger.info('Found %d input images for %d observations\n' % (
        len(inputfiles), len(Imagelist.keys())))

    times['beams']=time.time()
    ##############################
    # now loop through and do what needs to be done
    ##############################
    beams={}
    tempfiles=[]
    for obsid in sorted(Imagelist.keys()):
        logger.info('Processing %d...' % obsid)
        if Imagelist[obsid].ximage is None and Imagelist[obsid].yimage is not None:
            logger.error('Obsid %d has YY image but no XX image' % obsid)
            sys.exit(1)
        if Imagelist[obsid].yimage is None  and Imagelist[obsid].ximage is not None:
            logger.error('Obsid %d has XX image but no YY image' % obsid)
            sys.exit(1)
        s='%d: ' % obsid
        if Imagelist[obsid].ximage is not None:
            s+='(XX,YY)=(%s, %s)' % (Imagelist[obsid].ximage,
                                     Imagelist[obsid].yimage)
        if Imagelist[obsid].iimage is not None:
            s+=' I=%s' % Imagelist[obsid].iimage
        if Imagelist[obsid].metafits is not None:
            s+=' meta=%s' % Imagelist[obsid].metafits
        logger.info(s)

        if Imagelist[obsid].xbeam is not None or Imagelist[obsid].ibeam is not None:
            s='%d: ' % obsid
            if Imagelist[obsid].xbeam is not None:
                s+='(XXbeam,YYbeam)=(%s, %s)' % (Imagelist[obsid].xbeam,
                                                 Imagelist[obsid].ybeam)
            if Imagelist[obsid].ibeam is not None:
                s+=' Ibeam=%s' % Imagelist[obsid].ibeam
            logger.info(s)

        if Imagelist[obsid].iimage is None:
            # no I image exists: deal with XX, YY
            fx=pyfits.open(Imagelist[obsid].ximage)
            fy=pyfits.open(Imagelist[obsid].yimage)
            try:
                result=fx[0].header['DELAYS'] is not None and fx[0].header['DATE-OBS'] is not None
                xheader=True
            except KeyError:
                xheader=False
            try:
                result=fy[0].header['DELAYS'] is not None and fy[0].header['DATE-OBS'] is not None
                yheader=True
            except KeyError:
                yheader=False
            if not (xheader or yheader) or options.metafits:
                # does not have the necessary keys
                # need to make metafits
                if Imagelist[obsid].metafits is None:                
                    logger.warning('FITS header keywords not present for obsid %d; creating metafits...' % obsid)
                    hh=make_metafiles.Corr2UVFITSHeader(obsid, db=db)
                    T=make_metafiles.instrument_configuration(gpstime=obsid, duration=hh.obs.duration,
                                                              db=db)
                    result=T.make_instr_config(quick=True)
                    if result is None:
                        logger.error('Error making instr_config file for %d' % obsid)
                        sys.exit(1)
                    Imagelist[obsid].metafits='%d.metafits' % obsid
                    hh.make_header()
                    T.corr2uvfitsheader=hh
                    h=T.make_metafits(quick=True)                                    
                    if os.path.exists(Imagelist[obsid].metafits):
                        os.remove(Imagelist[obsid].metafits)
                    try:
                        h.writeto(Imagelist[obsid].metafits)
                        logger.info('Metafits for %d written to %s' % (obsid, Imagelist[obsid].metafits))
                    except Exception, e:
                        logger.error('Unable to write metafits file for %d:\n%s' % (obsid,e))
                        sys.exit(1)
            if not xheader:
                fx.close()
                if not Imagelist[obsid].metafits is None:
                    logger.warning('Copying metafits from %s to %s' % (
                        Imagelist[obsid].metafits,
                        Imagelist[obsid].ximage))
                    result=copy_metafits(Imagelist[obsid].ximage,
                                         Imagelist[obsid].metafits)
                else:
                    # assume we can get it from the YY image?
                    logger.warning('Copying metafits from %s to %s' % (
                        Imagelist[obsid].yimage,
                        Imagelist[obsid].ximage))
                    result=copy_metafits(Imagelist[obsid].ximage,
                                         Imagelist[obsid].yimage,
                                         keys=['DELAYS', 'DATE-OBS'])
                
                if not result:
                    logger.error('Error updating FITS header of %s' % Imagelist[obsid].ximage)
                    sys.exit(1)
            if not yheader:
                fy.close()
                if not Imagelist[obsid].metafits is None:
                    logger.warning('Copying metafits from %s to %s' % (
                        Imagelist[obsid].metafits,
                        Imagelist[obsid].yimage))               
                    result=copy_metafits(Imagelist[obsid].yimage,
                                         Imagelist[obsid].metafits)
                else:
                    # assume we can get it from the YY image?
                    logger.warning('Copying metafits from %s to %s' % (
                        Imagelist[obsid].ximage,
                        Imagelist[obsid].yimage))               
                    result=copy_metafits(Imagelist[obsid].yimage,
                                         Imagelist[obsid].ximage,
                                         keys=['DELAYS', 'DATE-OBS'])
                
                if not result:
                    logger.error('Error updating FITS header of %s' % Imagelist[obsid].yimage)
                    sys.exit(1)

            ##############################
            # now create the beams
            ##############################
            fx=pyfits.open(Imagelist[obsid].ximage,'update')
            fy=pyfits.open(Imagelist[obsid].yimage,'update')
            delays=[int(x) for x in fx[0].header['DELAYS'].split(',')]
            if Imagelist[obsid].xbeam is not None and Imagelist[obsid].ybeam is not None:
                logger.info('XX,YY beams already exists for %d; using %s,%s...' % (
                    obsid, Imagelist[obsid].xbeam, Imagelist[obsid].ybeam))

            elif not options.copy_beam or (not beams.has_key(fx[0].header['DELAYS'])):            
                logger.info('Creating primary beams for %d' % obsid)
                out=make_beam.make_beam(Imagelist[obsid].ximage, ext=0, delays=delays,
                                        analytic_model=options.analytic_model,
                                        jones=False,
                                        precess=options.precess)
                if out is None:
                    logger.error('Problem creating primary beam for %s' % obsid)
                    sys.exit(1)
                Imagelist[obsid].xbeam, Imagelist[obsid].ybeam=out
                beams[fx[0].header['DELAYS']]=obsid
                logger.info('Wrote %s,%s' % (out[0], out[1]))
                fx[0].header['BEAM']=Imagelist[obsid].xbeam
                fy[0].header['BEAM']=Imagelist[obsid].ybeam
                fx.flush()
                fy.flush()
                tempfiles+=out
            else:
                logger.info('Copying beam from %s to %s' % (Imagelist[beams[fx[0].header['DELAYS']]].ximage,
                                                            Imagelist[obsid].ximage))
                result=copy_driftscanbeam(Imagelist[beams[fx[0].header['DELAYS']]].ximage,
                                          Imagelist[obsid].ximage)
                if result is None:
                    sys.exit(1)
                Imagelist[obsid].xbeam, Imagelist[obsid].ybeam=result
                fx[0].header['BEAM']=Imagelist[obsid].xbeam
                fy[0].header['BEAM']=Imagelist[obsid].ybeam
                fx.flush()
                fy.flush()
                tempfiles+=result
            fx.close()
            fy.close()
            fx=pyfits.open(Imagelist[obsid].ximage)
            fy=pyfits.open(Imagelist[obsid].yimage)

            # and primary-beam correct
            fxB=pyfits.open(Imagelist[obsid].xbeam)
            fyB=pyfits.open(Imagelist[obsid].ybeam)
            fx[0].data=(fx[0].data/fxB[0].data+fy[0].data/fyB[0].data)/2.0
            fx[0].header['CRVAL%d' % polaxis]=1
            Imagelist[obsid].iimage=Imagelist[obsid].ximage.replace('_XX','_I')
            if os.path.exists(Imagelist[obsid].iimage):
                os.remove(Imagelist[obsid].iimage)
            fx.writeto(Imagelist[obsid].iimage)
            logger.info('Wrote I image %s' % Imagelist[obsid].iimage)
            tempfiles.append(Imagelist[obsid].iimage)
            fxB[0].data=(fxB[0].data+fyB[0].data)/2.0
            fxB[0].header['CRVAL%d' % polaxis]=1            
            Imagelist[obsid].ibeam=Imagelist[obsid].xbeam.replace('_beamXX','_beamI')
            Imagelist[obsid].ibeam=Imagelist[obsid].ibeam.replace('_XX','_I')
            if os.path.exists(Imagelist[obsid].ibeam):
                os.remove(Imagelist[obsid].ibeam)
            fxB.writeto(Imagelist[obsid].ibeam)
            logger.info('Wrote I beam %s' % Imagelist[obsid].ibeam)
            fi=pyfits.open(Imagelist[obsid].iimage, 'update')
            fi[0].header['BEAM']=Imagelist[obsid].ibeam
            fi.flush()
            fi.close()
            tempfiles.append(Imagelist[obsid].ibeam)
        else:
            # Stokes I (corrected) exists
            # now just make I beam
            fi=pyfits.open(Imagelist[obsid].iimage)
            try:
                result=fi[0].header['DELAYS'] is not None and fi[0].header['DATE-OBS'] is not None
                iheader=True
            except KeyError:
                iheader=False
            if not (iheader) or options.metafits:
                # does not have the necessary keys
                # need to make metafits
                if Imagelist[obsid].metafits is None:                
                    logger.warning('FITS header keywords not present for obsid %d; creating metafits...' % obsid)
                    hh=make_metafiles.Corr2UVFITSHeader(obsid, db=db)
                    T=make_metafiles.instrument_configuration(gpstime=obsid, duration=hh.obs.duration,
                                                              db=db)
                    result=T.make_instr_config(quick=True)
                    if result is None:
                        logger.error('Error making instr_config file for %d' % obsid)
                        sys.exit(1)
                    Imagelist[obsid].metafits='%d.metafits' % obsid
                    hh.make_header()
                    T.corr2uvfitsheader=hh
                    h=T.make_metafits(quick=True)                                    
                    if os.path.exists(Imagelist[obsid].metafits):
                        os.remove(Imagelist[obsid].metafits)
                    try:
                        h.writeto(Imagelist[obsid].metafits)
                        logger.info('Metafits for %d written to %s' % (obsid, Imagelist[obsid].metafits))
                    except Exception, e:
                        logger.error('Unable to write metafits file for %d:\n%s' % (obsid,e))
                        sys.exit(1)
            if not iheader:
                fi.close()
                if not Imagelist[obsid].metafits is None:
                    logger.warning('Copying metafits from %s to %s' % (
                        Imagelist[obsid].metafits,
                        Imagelist[obsid].iimage))
                    result=copy_metafits(Imagelist[obsid].iimage,
                                         Imagelist[obsid].metafits)
            if not result:
                logger.error('Error updating FITS header of %s' % Imagelist[obsid].iimage)
                sys.exit(1)
            ##############################
            # now create the beams
            ##############################            
            
            fi=pyfits.open(Imagelist[obsid].iimage, 'update')
            delays=[int(x) for x in fi[0].header['DELAYS'].split(',')]

            if Imagelist[obsid].ibeam is not None:
                logger.info('I beam already exists for %d; using %s...' % (
                    obsid, Imagelist[obsid].ibeam))
            else:
                if not options.copy_beam or (not beams.has_key(fi[0].header['DELAYS'])):
                    logger.info('Creating primary beams for %d' % obsid)
                    out=make_beam.make_beam(Imagelist[obsid].iimage, ext=0, delays=delays,
                                            analytic_model=options.analytic_model,
                                            jones=False,
                                            precess=options.precess)
                    if out is None:
                        logger.error('Problem creating primary beam for %s' % obsid)
                        sys.exit(1)
                    Imagelist[obsid].xbeam, Imagelist[obsid].ybeam=out
                    beams[fi[0].header['DELAYS']]=obsid
                    logger.info('Wrote %s,%s' % (out[0], out[1]))
                    tempfiles+=out
                else:
                    logger.info('Copying beam from %s to %s' % (Imagelist[beams[fi[0].header['DELAYS']]].iimage,
                                                                Imagelist[obsid].iimage))
                    result=copy_driftscanbeam(Imagelist[beams[fi[0].header['DELAYS']]].iimage,
                                              Imagelist[obsid].iimage)
                    if result is None:
                        sys.exit(1)
                    Imagelist[obsid].xbeam, Imagelist[obsid].ybeam=result
                    tempfiles+=result
                # create the I beam
                fxB=pyfits.open(Imagelist[obsid].xbeam)
                fyB=pyfits.open(Imagelist[obsid].ybeam)
                fxB[0].data=(fxB[0].data+fyB[0].data)/2.0
                fxB[0].header['CRVAL%d' % polaxis]=1
                Imagelist[obsid].ibeam=Imagelist[obsid].xbeam.replace('_beamXX','_beamI')
                Imagelist[obsid].ibeam=Imagelist[obsid].ibeam.replace('_XX','_I')
                if os.path.exists(Imagelist[obsid].ibeam):
                    os.remove(Imagelist[obsid].ibeam)
                fxB.writeto(Imagelist[obsid].ibeam)
                logger.info('Wrote I beam %s' % Imagelist[obsid].ibeam)
                # update the header so that the beam can be re-used in the future
                fi[0].header['BEAM']=Imagelist[obsid].ibeam
                fi.flush()
                fi.close()
                tempfiles.append(Imagelist[obsid].ibeam)
        logger.info('\n')

    times['weight']=time.time()
    logger.info('Creating weight images...')
    ##############################
    # now we should have correct I and I beams for everything
    # make 2D images and weights
    ##############################
    for obsid in sorted(Imagelist.keys()):
        root,ext=os.path.splitext(Imagelist[obsid].iimage)    
        if options.reduce:
            f=pyfits.open(Imagelist[obsid].iimage)
            f[0].data=f[0].data[0,0]
            if f[0].header['CRVAL1']<0:
                f[0].header['CRVAL1']+=360
            f.verify('fix')
            newfile=root + suffix
            if os.path.exists(newfile):
                os.remove(newfile)
            f.writeto(newfile)
            f.close()
            logger.info('Wrote 2D version of %s to %s' % (Imagelist[obsid].iimage, newfile))
            newroot,newext=os.path.splitext(newfile)
            Imagelist[obsid].iimage=newfile
            tempfiles.append(newfile)
        else:
            newroot, newext=root, ext

        weightfile=newroot + newsuffix
        f=pyfits.open(Imagelist[obsid].ibeam)
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
        logger.info('Wrote weight file to %s' % weightfile)
        Imagelist[obsid].weight=weightfile
        tempfiles.append(weightfile)

    times['swarp']=time.time()
    logger.info('\n')
    logger.info('Swarping...')
    # and now swarp together
    # get image centers
    RAcenters=[]
    Deccenters=[]
    images=[]
    for obsid in sorted(Imagelist.keys()):
        f=pyfits.open(Imagelist[obsid].iimage)
        RAcenters.append(f[0].header['CRVAL1'])
        Deccenters.append(f[0].header['CRVAL2'])
        images.append(Imagelist[obsid].iimage)
    RAcenters=numpy.array(RAcenters)
    RAcenter=numpy.degrees(circular_mean(numpy.radians(RAcenters)))
    if options.ra is None:
        # make integer hours
        RAcenterhours=Angle(int(round(RAcenter/15.0)), unit=u.hour)
    else:
        if options.ra=='middle':
            RAcenter=RAcenters[len(RAcenters)/2]
            # make integer hours
            RAcenterhours=Angle(int(round(RAcenter/15.0)), unit=u.hour)
        else:
            try:
                RAcenterhours=Angle(float(options.ra), unit=u.hour)
            except:
                logger.warning('Unable to interpret ra center %s' % options.ra)
                RAcenterhours=Angle(int(round(RAcenter/15.0)), unit=u.hour)

    # for a single Dec strip these would all be the same
    # but allow for multiple Dec's
    Deccenter=numpy.unique(Deccenters).mean()
    # command-line option Dec overwrites automatic Dec
    if options.dec is None:
        Deccenter=Angle(Deccenter, unit=u.deg)
    else:
        try:
            Deccenter=Angle(float(options.dec), unit=u.deg)
        except:
            logger.warning('Unable to interpret dec center %s' % options.dec) 
            Deccenter=Angle(Deccenter, unit=u.deg)
            
    # And the image width and height
    width=(f[0].header['NAXIS1'])
    height=(f[0].header['NAXIS2'])
    
    center='%s,%s' % (RAcenterhours.to_string(unit=u.hour, sep=':'),
                                     Deccenter.to_string(unit=u.deg, sep=':'))
    logger.info('Determined center RA=%sh, Dec=%sd' % (RAcenterhours.to_string(unit=u.hour, sep=':'),
                                                      Deccenter.to_string(unit=u.deg, sep=':')))

    command='%s -c %s %s -OVERSAMPLING %d -CENTER %s -IMAGE_SIZE 0' % (swarp,
                                                                       localswarpconfig,
                                                                       ','.join(images),
                                                                       options.oversample,
                                                                       center)
    command+=' -PROJECTION_TYPE %s' % options.projection
    if len(options.swarpcommands)>0:
        command+=' ' + options.swarp
    outname='%s.swarp.fits' % options.root
    outweightname='%s.weight.fits' % options.root
    
    command+=' -IMAGEOUT_NAME %s -WEIGHTOUT_NAME %s' % (outname,
                                                        outweightname)
    logger.info(command)
    p=subprocess.Popen(command.split(),
                       stderr=subprocess.PIPE)
    imagenum=0
    if os.path.exists('swarp.log'):
        os.remove('swarp.log')
    log=open('swarp.log','w')
    for line in p.stderr:
        if not 'line:' in line:
            log.write(line)
        if '-------------- File' in line:
            currentimage=line.split()[-1].replace(':','')
            imagenum+=1
            logger.info('Working on image %04d/%04d...' % (imagenum, len(images)))
        if 'Co-adding frames' in line:
            logger.info('Coadding...')
    log.close()
    
    times['end']=time.time()

    logger.info('Execution time: %d s (setup), %d s (beams), %d s (weights), %d s (swarp) = %d s (total)' % (
        times['beams']-times['start'],
        times['weight']-times['beams'],
        times['swarp']-times['weight'],
        times['end']-times['swarp'],
        times['end']-times['start']))
    
    if options.clean:
        for f in tempfiles:
            os.remove(f)
            
    sys.exit(0)
    
################################################################################

if __name__=="__main__":
    main()
                
