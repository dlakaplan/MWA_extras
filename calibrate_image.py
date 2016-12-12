"""
Todo:

-quality metrics:
  - # of sources detected
-sub-bands & sub-exposures together

To debug:
-table output
-table output of null strings
-table output of filenames
-output of CASA calibrator file in a different directory

"""

import logging,logging.handlers,datetime,math,sys,socket,os,shutil,io,re
try:
    import matplotlib
    matplotlib.use('PDF')
    from matplotlib import pyplot as plt
    _matplotlib=True
except:
    _matplotlib=False
from optparse import OptionParser,OptionGroup
import time,tempfile
import subprocess,fcntl
from astropy.table import Table,Column
import collections,glob,numpy
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy
from astropy import constants as c, units as u
from astropy.wcs import WCS

import extra_utils

logger=extra_utils.makelogger('calibrate_image')
    

import mwapy
from mwapy import metadata
from mwapy.pb import make_beam

# try to use Peter Williams' casac interface
try:
    import drivecasa
    _OLDCASA=True
except ImportError:
    #logger.error('Unable to import drivecasa')
    _OLDCASA=False

try:
    import casac
    _CASA=True
    import pwkit.environments.casa.tasks as tasklib
    import pwkit.environments.casa.util as casautil
except ImportError:
    logger.warning('Unable to import casacore')
    _CASA=False

try:
    from AegeanTools.regions import Region
    from AegeanTools import source_finder, BANE, MIMAS
    # make the logging output from aegean more reasonable
    logging.getLogger("Aegean").setLevel(logging.INFO)
    _aegean=True
except ImportError:
    logger.error('Unable to import aegean')
    _aegean=False



##################################################
# default search paths for ANOKO/mwa-reduce
# the ANOKO path should contain calibrate and applysolutions
##################################################
anokopath=['~kaplan/mwa/anoko/mwa-reduce/build/', 
           '/usr/physics/mwa/pkg/anoko/mwa-reduce/build/']



catalogdir=os.path.join(os.path.split(os.path.abspath(__file__))[0],'catalogs')
#calmodelfile=os.path.join(catalogdir,'model_a-team.txt')
anokocatalog=os.path.join(catalogdir,'model-catalogue_new.txt')
calmodelfile=anokocatalog
anoko=None

# this defines aliases for source names in the anoko
# calibrator model file
# e.g., if the metafits says "PictorA", also look for "PicA" in the model file
# can be many -> one mapping
calmodelaliases={'PicA': ['PictorA'],
                 'PKS2331-41': ['J2334-4'],
                 'HydA': ['HydraA']}

# invert the alias dictionary for faster lookup
calmodelaliases_inverse={}
for k in calmodelaliases.keys():
    for v in calmodelaliases[k]:
        calmodelaliases_inverse[v]=k

# these are sources that don't calibrate well with Anoko
badanokosources=['3C444','PKS2356-61']

if not os.path.exists(calmodelfile):
    logger.warning('Unable to find calibrator model file %s' % calmodelfile)
if not os.path.exists(anokocatalog):
    logger.warning('Unable to find anoko autoprocess catalog file %s' % anokocatalog)

brightsources={'3C353': SkyCoord('17h20m28.1s','-00d58m47s'),
               '3C409': SkyCoord('20h14m27.6s','23d34m53s'),
               '3C444': SkyCoord('-01h45m34.4s','-17d01m30s'),
               'CasA': SkyCoord('23h23m27.9s','58d48m42s'),
               'CenA': SkyCoord('13h25m27.6s','-43d01m09s'),
               'Crab': SkyCoord('05h34m32.0s','22d00m52s'),
               'CygA': SkyCoord('19h59m28.4s','40d44m02s'),
               'ESO362-G021': SkyCoord('05h22m58.0s','-36d27m31s'),
               'ForA': SkyCoord('3h22m41.7s','-37d12m30s'),
               'HerA': SkyCoord('-07h08m52.1s','04d59m16s'),
               'HydA': SkyCoord('09h18m05.1s','-12d05m42s'),
               'NGC253': SkyCoord('00h47m34.7s','-25d17m30s'),
               'PicA': SkyCoord('05h19m30.7s','-45d45m35s'),
               'PKS2153-69': SkyCoord('21h57m06.0s','-69d41m24s'),
               'PKS2331-41': SkyCoord('23h34m26.2s','-41d25m24s'),
               'PKS2356-61': SkyCoord('23h59m04.3s','-60d54m59s'),
               'PKSJ0130-2610': SkyCoord('01h30m27.8s','-26d09m56s'),
               'VirA': SkyCoord('12h30m49.4s','12d23m28s')}

##################################################
def non_block_read(output):
    """
    allows for non-blocking read during a subprocess call
    https://gist.github.com/sebclaeys/1232088
    """
    fd = output.fileno()
    fl = fcntl.fcntl(fd, fcntl.F_GETFL)
    fcntl.fcntl(fd, fcntl.F_SETFL, fl | os.O_NONBLOCK)
    try:
        return output.read()
    except:
        return ""

##################################################
def checkforlock(ms):
    """
    checks for a CASA lock file
    if it finds it, it will delete
    """
    for root, dirs, files in os.walk(ms):
        for name in files:
            if name == 'table.lock': 
                logger.debug('Removing CASA table lock from %s' % os.path.join(root,name))
                os.remove(os.path.join(root,name))


##################################################
def getstatus(p, output=None, error=None):
    """
    pollresults=getstatus(p, output=None, error=None)
    p is a list of subprocess instances
    this checks if they are finished
    and returns a list of False if still running or True if stopped
    """
    poll=[]
    for j in xrange(len(p)):
        try:
            poll.append(p[j].poll() is not None)
        except:
            poll.append(None)
        if output is not None:
            try:
                output[j].flush()
            except:
                pass
        if error is not None:
            try:
                error[j].flush()
            except:
                pass
    return poll
##################################################
def stat_measure(image, fraction=0.5, nonan=True):
    """
    median, rms = stat_measure(image, fraction=0.5, nonan=True)
    uses the central fraction of the image to compute
    the median and rms (using inner-quartile range)
    """

    if fraction>1:
        fraction=1
    if isinstance(image,str):
        try:
            f=fits.open(image)
        except Exception,e:
            logger.error('Unable to open FITS image %s:\n\t%s' % (image,e))
            return None,None
    else:
        f=image

    if 'RA' in f[0].header['CTYPE1']:
        data=f[0].data
    else:
        data=f[0].data.T
    # now it should be Stokes, Freq, Dec, RA
    ny,nx=data.shape[2],data.shape[3]
    ytouse=ny*fraction
    ystart=(ny/2)-ytouse/2
    ystop=(ny/2)+ytouse/2
    xtouse=nx*fraction
    xstart=(nx/2)-xtouse/2
    xstop=(nx/2)+xtouse/2
    d=(data[:,:,int(ystart):int(ystop),int(xstart):int(xstop)]).flatten()
    if nonan:
        d=d[~numpy.isnan(d)]
    q1,m,q2=numpy.percentile(d, [25,50,75])
    return m,(q2-q1)/1.35

##################################################
def match_aegean_sources(sourcelist1, sourcelist2):
    matchradius=10*u.arcsec
    I1=[]
    I2=[]

    coords1=SkyCoord([s.ra for s in sourcelist1],
                     [s.dec for s in sourcelist1],unit=(u.deg,u.deg))
    coords2=SkyCoord([s.ra for s in sourcelist2],
                     [s.dec for s in sourcelist2],unit=(u.deg,u.deg))
    S1=[]
    S2=[]
    for i in xrange(len(coords1)):
        d=coords1[i].separation(coords2)
        if d.min() < matchradius:
            if (not i in I1) and (not numpy.argmin(d) in I2):
                I1.append(i)
                I2.append(numpy.argmin(d))
    for i in xrange(len(I1)):
        S1.append(sourcelist1[I1[i]])
        S2.append(sourcelist2[I2[i]])
    return S1,S2

##################################################
def circle2mimas(ra, dec, radius, filename):
    """
    circle2mimas(ra, dec, radius)
    ra,dec,radius in decimal regrees
    """
    r=Region()
    r.add_circles(numpy.radians(ra),
                  numpy.radians(dec),
                  numpy.radians(radius))
    MIMAS.save_region(r, filename)
    return filename


##################################################
class ANOKOfinder():
    """
    """
    path=anokopath

    def __init__(self, *args):
        """
        """

        if len(args)>0:
            self.paths= list(args) + ANOKOfinder.path
        else:
            self.paths=ANOKOfinder.path

    def find(self):
        actualanoko=None
        if os.environ.has_key('ANOKO'):
            actualanoko=os.environ['ANOKO']
            logger.debug('Identified ANOKO %s from $ANOKO' % actualanoko)
        else:
            for path in self.paths:
                if path is None:
                    continue
                path=os.path.expanduser(path)
                if os.path.exists(path) and os.path.exists(os.path.join(path,'calibrate')) and os.path.exists(os.path.join(path,'applysolutions')):
                    actualanoko=path
                    logger.debug('Identified ANOKO %s' % actualanoko)
                    break
                
        if actualanoko is None:
            logger.error('No ANOKO identified')
            return None

        if not os.path.exists(actualanoko):
            logger.error('ANOKO %s does not exist' % actualanoko)
            return None
        if not (os.path.exists(os.path.join(actualanoko,'calibrate')) and os.path.exists(os.path.join(actualanoko,'applysolutions'))):
            logger.error('ANOKO executables %s and %s do not exist' % ('calibrate','applysolutions'))
            return None

        return actualanoko


##################################################
def makemetafits(obsid, directory=None):        
    """
    makemetafits(obsid)
    returns metafits name on success
    or None on failure
    """
    m=metadata.instrument_configuration(int(obsid))
    h=m.make_metafits()
    if directory is None:
        metafits='%s.metafits' % obsid
    else:
        metafits=os.path.join(directory,'%s.metafits' % obsid)
    if os.path.exists(metafits):
        os.remove(metafits)
    try:
        h.writeto(metafits)
        logger.info('Metafits written to %s' % (metafits))
        return metafits
    except Exception, e:
        logger.error('Unable to write metafits file %s:\n%s' % (metafits,e))
        return None


##################################################
def get_msinfo(msfile):
    """
    chanwidth, inttime, otherkeys=get_msinfo(msfile)
    chanwidth is channel width in kHz
    inttime is integration time in s
    otherkeys is a dictionary containing other header keywords
    """
    if not _CASA:
        logger.error("CASA operation not possible")
        return None

    
    ms=casac.casac.ms()
    ms.open(msfile)
    chanwidth=ms.getspectralwindowinfo()["0"]['ChanWidth']/1e3
    nchans=ms.getspectralwindowinfo()["0"]['NumChan']
    inttime=ms.getscansummary()['1']['0']['IntegrationTime']
    reffreq=ms.getspectralwindowinfo()["0"]['RefFreq']
    ms.close()
    t=casac.casac.table()
    t.open(msfile)
    otherkeys=t.getkeywords()
    t.close()

    logger.debug('Determined channel width of %d kHz for %s' % (chanwidth,msfile))
    logger.debug('Determined integration time of %.1f s for %s' % (inttime,msfile))
    logger.debug('Determined %d channels for %s' % (nchans,msfile))
    logger.debug('Determined reference frequency of %.1f MHz for %s' % (reffreq/1e6,msfile))

    return chanwidth, inttime, nchans, reffreq, otherkeys

    

##################################################
def check_calibrated(msfile):
    """
    check_calibrated(msfile)
    checks for presence of CORRECTED_DATA column in a ms file
    """
    if not _CASA:
        logger.error("CASA operation not possible")
        return None

    t=casac.casac.table()
    t.open(msfile)
    return "CORRECTED_DATA" in t.colnames()

##################################################
# from mwapy/casac_ft_beam
# based on mwapy/ft_beam
# Generate automatic calibration model and form a bandpass solution

# Natasha Hurley-Walker 10/07/2013
# Updated 08/08/2013 to scale the YY and XX beams separately
# Updated 01/10/2013 Use the field name as the calibrator name if the calibrator wasn't filled in properly during scheduling
# Updated 21/11/2013 Added sub-calibrators to complex fields (but didn't find much improvement)
# Updated 02/12/2013 Added a spectral beam option; turned subcalibrators off by default
# Updated 10/03/2014 Try to get the calibrator information from the metafits file
# Updated 18/08/2014 Improved astropy/pyfits compatibility and added option to switch off beam correction
# updated 18/08/2016 to use Peter William's casa python bindings instead of being called from within casa environment (DLK)



##################################################
def importfits(fitsimage, imagename):
   """
   A simplified version of task_importfits
   """
   ia = casautil.tools.image ()
   ia.fromfits(imagename, fitsimage)

##################################################

def ft_beam(vis=None,refant='Tile012',clobber=True,correct_beam=True,
            spectral_beam=False,
            subcalibrator=False,uvrange='>0.03klambda',
            outdir='./'):
    """
    def ft_beam(vis=None,refant='Tile012',clobber=True,correct_beam=True
    spectral_beam=False,
    subcalibrator=False,uvrange='>0.03klambda',outdir='./'):
    
    # Reference antenna
    refant='Tile012'
    # Overwrite files
    clobber=True
    # Option to correct for the attenuation effects of the primary beam
    correct_beam=True
    # Option to include the spectral index of the primary beam
    spectral_beam=False
    # Option to add more sources to the field
    """
    modeldir=os.environ['MWA_CODE_BASE']+'/MWA_Tools/Models/'
    if not os.path.exists(modeldir):
        logger.error('Model directory %s does not exist' % modeldir)
        return None

    qa = casautil.tools.quanta()

    # output calibration solution
    caltable=os.path.join(outdir,re.sub('ms','cal',vis))
    if vis is None or len(vis)==0 or not os.path.exists(vis):
        logger.error('Input visibility must be defined')
        return None
   
    # Get the frequency information of the measurement set
    ms=casac.casac.ms()
    ms.open(vis)
    rec = ms.getdata(['axis_info'])
    df,f0 = (rec['axis_info']['freq_axis']['resolution'][len(rec['axis_info']['freq_axis']['resolution'])/2],rec['axis_info']['freq_axis']['chan_freq'][len(rec['axis_info']['freq_axis']['resolution'])/2])
    F =rec['axis_info']['freq_axis']['chan_freq'].squeeze()/1e6
    df=df[0]*len(rec['axis_info']['freq_axis']['resolution'])
    f0=f0[0]
    rec_time=ms.getdata(['time'])
    sectime=qa.quantity(rec_time['time'][0],unitname='s')
    midfreq=f0
    bandwidth=df
    if isinstance(qa.time(sectime,form='fits'),list):
        dateobs=qa.time(sectime,form='fits')[0]
    else:
        dateobs=qa.time(sectime,form='fits')

    if spectral_beam:
        # Start and end of the channels so we can make the spectral beam image
        startfreq=f0-df/2
        endfreq=f0+df/2
        freq_array=[midfreq,startfreq,endfreq]
    else:
        freq_array=[midfreq]

    # Get observation number directly from the measurement set
    tb=casac.casac.table()
    tb.open(vis+'/OBSERVATION')
    obsnum=int(tb.getcol('MWA_GPS_TIME'))
    tb.close

    # Try getting the calibrator information from the metafits file

    metafits=str(obsnum)+'.metafits'
    calibrator=""
    if os.path.exists(metafits):
        hdu_in=fits.open(metafits)
        try:
            calibrator=hdu_in[0].header['CALIBSRC']
            str_delays=hdu_in[0].header['DELAYS']
            delays=[int(x) for x in str_delays.split(',')]
        except:
            logger.warning('Unable to retrieve calibrator from metafits file.')

    if not calibrator:
        logger.warning('Could not use the metafits file: trying the observation database.')
        info=metadata.MWA_Observation(obsnum)
        logger.info('Retrieved observation info for %d...\n%s\n' % (obsnum,info))

        # Calibrator information
        if info.calibration:
            calibrator=info.calibrators
        else:
            # Observation wasn't scheduled properly so calibrator field is missing: try parsing the fieldname
            # assuming it's something like 3C444_81
            calibrator=info.filename.rsplit('_',1)[0]
        # Delays
        delays=info.delays
        str_delays=','.join(map(str,delays))

    logger.info('Calibrator is %s...' % calibrator)
    logger.info('Delays are: %s' % str_delays)

    # subcalibrators not yet improving the calibration, probably due to poor beam model
    if subcalibrator and calibrator=='PKS0408-65':
        subcalibrator='PKS0410-75'
    elif subcalibrator and calibrator=='HerA':
        subcalibrator='3C353'
    else:
        subcalibrator=False

    # Start models are 150MHz Jy/pixel fits files in a known directory
    model=modeldir+calibrator+'.fits'
    # With a corresponding spectral index map
    spec_index=modeldir+calibrator+'_spec_index.fits'

    if not os.path.exists(model):
        logger.error('Could not find calibrator model %s' % model)
        return None
   
    # load in the model FITS file as a template for later
    ftemplate=fits.open(model)

    # do this for the start, middle, and end frequencies
    for freq in freq_array:
        freqstring=str(freq/1.0e6) + 'MHz'
        # We'll generate images in the local directory at the right frequency for this ms
        outname=os.path.join(outdir,calibrator+'_'+freqstring)
        outnt2=os.path.join(outdir,calibrator+'_'+freqstring+'_nt2')
        # import model, edit header so make_beam generates the right beam in the right place
        if os.path.exists(outname + '.fits') and clobber:
            os.remove(outname + '.fits')
        shutil.copy(model,outname + '.fits')
        fp=fits.open(outname + '.fits','update')
        try:
            fp[0].header['CRVAL3']=freq
            fp[0].header['CDELT3']=bandwidth
            fp[0].header['DATE-OBS']=dateobs
        except KeyError:
            fp[0].header.update('CRVAL3',freq)
            fp[0].header.update('CDELT3',bandwidth)
            fp[0].header.update('DATE-OBS',dateobs)
        fp.flush()
  
        logger.info('Creating primary beam models...')
        beamarray=make_beam.make_beam(outname + '.fits',
                                      model='analytic',
                                      delays=delays)
        # delete the temporary model
        os.remove(outname + '.fits')

        beamimage={}

        for stokes in ['XX','YY']:
            beamimage[stokes]=os.path.join(outdir, calibrator + '_' + freqstring + '_beam' + stokes + '.fits')

    # scale by the primary beam
    # Correct way of doing this is to generate separate models for XX and YY
    # Unfortunately, ft doesn't really understand cubes
    # So instead we just use the XX model, and then scale the YY solution later

    freq=midfreq
    freqstring=str(freq/1.0e6)+'MHz'
    outname=os.path.join(outdir,calibrator+'_'+freqstring)
    outnt2=os.path.join(outdir,calibrator+'_'+freqstring+'_nt2')
    if isinstance(outname, unicode):      
        outname=outname.encode('ascii')
    if isinstance(outnt2, unicode):      
        outnt2=outnt2.encode('ascii')

    # divide to make a ratio beam, so we know how to scale the YY solution later
    fbeamX=fits.open(beamimage['XX'])
    fbeamY=fits.open(beamimage['YY'])
    if correct_beam:
        ratiovalue=(fbeamX[0].data/fbeamY[0].data).mean()
        logger.info('Found <XX/YY>=%.2f' % ratiovalue)
    else:
        ratiovalue=1.0
        logger.info('Not correcting for the beam: using <XX/YY>=%.2f' % ratiovalue)

    # Models are at 150MHz
    # Generate scaled image at correct frequency
    if os.path.exists(outname + '.fits') and clobber:
        os.remove(outname + '.fits')
    # Hardcoded to use the XX beam in the model
    fbeam=fbeamX
    fmodel=fits.open(model)
    fspec_index=fits.open(spec_index)

    if correct_beam:
        ftemplate[0].data=fbeam[0].data * fmodel[0].data/((150000000/f0)**(fspec_index[0].data))
    else:
        ftemplate[0].data=fmodel[0].data/((150000000/f0)**(fspec_index[0].data))
       
    try:
        ftemplate[0].header['CRVAL3']=freq
        ftemplate[0].header['CDELT3']=bandwidth
        ftemplate[0].header['DATE-OBS']=dateobs
        ftemplate[0].header['CRVAL4']=1
    except KeyError:
        ftemplate[0].header.update('CRVAL3',freq)
        ftemplate[0].header.update('CDELT3',bandwidth)
        ftemplate[0].header.update('DATE-OBS',dateobs)
        ftemplate[0].header.update('CRVAL4',1)
        
    ftemplate.writeto(outname + '.fits')
    logger.info('Wrote scaled model to %s' % (outname + '.fits'))
    foutname=fits.open(outname + '.fits')
   
    # Generate 2nd Taylor term
    if os.path.exists(outnt2 + '.fits') and clobber:
        os.remove(outnt2 + '.fits')

    if correct_beam:
        if spectral_beam:
            # Generate spectral image of the beam
            fcalstart=fits.open(calibrator+'_'+str(startfreq/1.0e6)+'MHz_beamXX.fits')
            fcalend=fits.open(calibrator+'_'+str(endfreq/1.0e6)+'MHz_beamXX.fits')   
            ftemplate[0].data=(numpy.log(fcalstart[0].data/fcalend[0].data)/
                               numpy.log((f0-df/2)/(f0+df/2)))
            beam_spec='%s_%sMHz--%sMHz_beamXX.fits' % (calibrator,
                                                       str(startfreq/1.0e6),
                                                       str(endfreq/1.0e6))
            if os.path.exists(beam_spec):
                os.remove(beam_spec)
            ftemplate.writeto(beam_spec)
            fbeam_spec=fits.open(beam_spec)
       
            ftemplate[0].data=foutname[0].data * fbeam[0].data * (fspec_index[0].data+fbeam_spec[0].data)
        else:
            ftemplate[0].data=foutname[0].data * fbeam[0].data * fspec_index[0].data
    else:
        ftemplate[0].data=foutname[0].data

    try:
        ftemplate[0].header['DATE-OBS']=dateobs
    except KeyError:
        ftemplate[0].header.update('DATE-OBS',dateobs)
      
    ftemplate.writeto(outnt2 + '.fits')
    logger.info('Wrote scaled Taylor term to %s' % (outnt2 + '.fits'))

    # import as CASA images
    if os.path.exists(outname + '.im') and clobber:
        shutil.rmtree(outname + '.im')
    if os.path.exists(outnt2 + '.im') and clobber:
        shutil.rmtree(outnt2 + '.im')
    # use the local versions of this 
    # that I defined above
    importfits(outname + '.fits',outname + '.im')
    importfits(outnt2 + '.fits',outnt2 + '.im')
    if not os.path.exists(outname + '.im'):
        logger.error('Cannot find %s' % (outname + '.im'))
        return None
    if not os.path.exists(outnt2 + '.im'):
        logger.error('Cannot find %s' % (outnt2 + '.im'))
        return None
   
    logger.info('Fourier transforming model...')
    cfg=tasklib.FtConfig()
    cfg.vis=vis
    cfg.model=[outname + '.im',outnt2+'.im']
    cfg.usescratch=True
    tasklib.ft(cfg)
    
    logger.info('Calibrating...')
    cfg=tasklib.GaincalConfig()
    cfg.vis=vis
    cfg.caltable=caltable
    cfg.refant=refant
    cfg.uvrange=uvrange
    cfg.gaintype = 'B'
    cfg.combine = ['scan']
    cfg.solint = 'inf'
    cfg.solnorm = False
    tasklib.gaincal(cfg)

    logger.info('Scaling YY solutions by beam ratio...')
    # Scale YY solution by the ratio
    tb.open(caltable)
    G = tb.getcol('CPARAM')
    tb.close()
    
    new_gains = numpy.empty(shape=G.shape, dtype=numpy.complex128)
    
    # XX gains stay the same
    new_gains[0,:,:]=G[0,:,:]
    # YY gains are scaled
    new_gains[1,:,:]=ratiovalue*G[1,:,:]

    tb.open(caltable,nomodify=False)
    tb.putcol('CPARAM',new_gains)
    tb.putkeyword('MODEL',model)
    tb.putkeyword('SPECINDX',spec_index)
    tb.putkeyword('BMRATIO',ratiovalue)
    try:
        tb.putkeyword('MWAVER',mwapy.__version__)
    except:
        pass
    tb.close()
    logger.info('Created %s!' % caltable)
    return caltable

##################################################
def calibrate_casa(obsid, directory=None, minuv=60):
    """
    calibrate_casa(obsid, directory=None, minuv=60)
    minuv in meters
    returns name of calibration (gain) file on success
    returns None on failure
    """
    if not _CASA:
        logger.error("CASA operation not possible")
        return None

    basedir=os.path.abspath(os.curdir)
    if directory is None:
        directory=basedir
    if minuv is not None:
        result=ft_beam(vis='%s.ms' % obsid,
                       uvrange='>%dm' % minuv,
                       outdir='%s' % directory)
    else:
        result=ft_beam(vis='%s.ms' % obsid,
                       outdir='%s' % directory)
    if result is None:
        logger.error('Unable to create calibration table')
        return None
    else:
        outfile=result
    # that file should be the same as the expected output
    if not os.path.split(outfile)[-1] == '%s.cal' % obsid:
        logger.error('CASA calibration produced %s, but expected %s.cal' % (outfile,
                                                                            obsid))
        return None
    return outfile



##################################################
def selfcalibrate_casa(obsid, suffix=None, directory=None, minuv=60):
    """
    selfcalibrate_casa(obsid, directory=None, minuv=60)
    minuv in meters
    returns name of calibration (gain) file on success
    returns None on failure
    """
    if not _CASA:
        logger.error("CASA operation not possible")
        return None
    if suffix is None:
        calfile='%s.cal' % obsid
    else:
        calfile='%s_%s.cal' % (obsid,suffix)
    if directory is not None:
        calfile=os.path.join(directory, calfile)

    basedir=os.path.abspath(os.curdir)
        
    if directory is None:
        directory=basedir
    if minuv is not None:
        uvrange='>%dm' % minuv
    else:
        uvrange=''
    cfg=tasklib.GaincalConfig()
    cfg.vis='%s.ms' % obsid
    cfg.caltable=calfile
    cfg.refant='Tile012'
    cfg.uvrange=uvrange
    cfg.gaintype = 'B'
    cfg.combine = ['scan']
    cfg.solint = 'inf'
    cfg.solnorm = False
    try:
        tasklib.gaincal(cfg)
    except:
        return None
    if not os.path.exists(calfile):
        logger.error('Expected CASA output %s does not exist' % calfile)
        return None
    return calfile


##################################################
def applycal_casa(obsid, calfile):
    """
    applycal_casa(obsid, calfile)
    returns True on success, None on failure
    """
    if not _CASA:
        logger.error("CASA operation not possible")
        return None
    basedir=os.path.abspath(os.curdir)

    cfg=tasklib.ApplycalConfig()
    cfg.vis='%s.ms' % obsid
    cfg.gaintable=[calfile]
    try:
        tasklib.applycal(cfg)
    except:
        return None
    return True


##################################################
def extract_calmodel(filename, sourcename):
    """
    extract_calmodel(filename, sourcename)
    extracts just the data for given <sourcename> from the
    master model file <filename>
    returns the corresponding lines or None on failure
    """
    logger.debug('Extracting calibration model for %s from file %s' % (sourcename,
                                                                       filename))
    try:
        f=open(filename)
    except Exception, e:
        logger.error('Unable to read calibrator model file %s:\n\t%s' % (filename,e))
        return None
    lines=f.readlines()
    i=0
    istart=None
    iend=None
    while i < len(lines):
        if lines[i].startswith('source'):
            d=lines[i+1].split()
            if d[0]=='name' and (sourcename in d[1] or (sourcename in calmodelaliases_inverse.keys() and calmodelaliases_inverse[sourcename] in d[1])):
                istart=i
                i+=1
                while i < len(lines):
                    if lines[i][0]=='}':
                        break
                    i+=1
                iend=i
                break
        i+=1
    if istart is None or iend is None:
        logger.error('Unable to find source "%s" in calibrator model file %s' % (sourcename,
                                                                                 filename))
        return None
    return [lines[0]]+lines[istart:iend+1]

##################################################
def write_calmodelfile(calibrator_name, directory=None):
    """
    write_calmodelfile(calibrator_name, directory=None)
    extracts the relevant info from the master calibrator model file
    and writes it to a single-source file
    returns the name of that file on success, or None on failure
    """
    
    # get a model file
    calibrator_modeldata=extract_calmodel(calmodelfile,
                                          calibrator_name)
    if calibrator_modeldata is None:
        logger.error('No calibrator data found')
        return None
    outputcalmodelfile='%s_%s.model' % (calibrator_name,
                                        datetime.datetime.now().strftime('%Y%m%d'))
    if directory is not None:
        outputcalmodelfile=os.path.join(directory,outputcalmodelfile)
    try:
        f=open(outputcalmodelfile,'w')
    except Exception, e:
        logger.error('Unable to open file %s for writing:\n\t%s' % (outputcalmodelfile,e))
        return None
    for line in calibrator_modeldata:
        f.write(line)
    f.close()
    logger.debug('Wrote calibration model file %s' % outputcalmodelfile)
    return outputcalmodelfile
            
##################################################
def calibrate_anoko(obsid,outputcalmodelfile=None, 
                    minuv=60,
                    suffix=None, 
                    corrected=False,
                    directory=None,
                    ncpus=32):
    """
    calibrate_anoko(obsid,outputcalmodelfile=None, minuv=60,
    suffix=None, 
    corrected=False,
    directory=None,
    ncpus=32)
    returns name of calibration (gain) file on success
    """

    if suffix is None:
        calfile='%s.cal' % obsid
    else:
        calfile='%s_%s.cal' % (obsid,suffix)
    if directory is not None:
        calfile=os.path.join(directory, calfile)
    calibratecommand=[os.path.join(anoko,'calibrate')]
    if minuv is not None and minuv > 0:
        calibratecommand+=['-minuv',
                           str(minuv)]
                                    
    if corrected:
        calibratecommand+=['-datacolumn',
                           'CORRECTED_DATA']
    calibratecommand+=['-j',
                       str(ncpus),
                       '-a',
                       '0.001',
                       '0.0001']
    if outputcalmodelfile is not None:
        calibratecommand+=['-m',
                           outputcalmodelfile]
    calibratecommand+=['%s.ms' % obsid,
                       calfile]
    logger.info('Will run:\n\t%s' % ' '.join(calibratecommand))
    p=subprocess.Popen(' '.join(calibratecommand),
                       stderr=subprocess.PIPE,
                       stdout=subprocess.PIPE,
                       shell=True,
                       close_fds=False)
    while True:
        p.stdout.flush()
        p.stderr.flush()
        for l in p.stdout.readlines():
            logger.debug(l.rstrip())
        for l in p.stderr.readlines():
            logger.error(l.rstrip())
        returncode=p.poll()
        if returncode is not None:
            break
        time.sleep(1)
    return calfile
        

##################################################
def applycal_anoko(obsid, calfile, corrected=False):
    """
    applycal_anoko(obsid, calfile, corrected=False):
    returns True on success
    """
    applycalcommand=[os.path.join(anoko,'applysolutions')]
    if corrected:
        applycalcommand+=['-datacolumn',
                          'CORRECTED_DATA']
    applycalcommand+=['-copy',
                      '%s.ms' % obsid,
                      calfile]
    logger.info('Will run:\n\t%s' % ' '.join(applycalcommand))
    p=subprocess.Popen(' '.join(applycalcommand),
                       stderr=subprocess.PIPE,
                       stdout=subprocess.PIPE,
                       shell=True,
                       close_fds=False)
    while True:
        p.stdout.flush()
        p.stderr.flush()
        for l in p.stdout.readlines():
            logger.debug(l.rstrip())
        for l in p.stderr.readlines():
            logger.error(l.rstrip())
        returncode=p.poll()
        if returncode is not None:
            break
        time.sleep(1)
    return True

######################################################################
def identify_calibrators(observation_data):
    calibrators=numpy.where(observation_data['iscalibrator'])[0]
    notcalibrators=numpy.where(~observation_data['iscalibrator'])[0]
    if len(calibrators)==0:
        logger.warning('No calibrators identified')
        return calibrators, notcalibrators, {}
    
    cal_observations={}
    for i in notcalibrators:
        # find which match in freq and are not calibrators
        good=(observation_data['channel']==observation_data[i]['channel']) & observation_data['iscalibrator']
        if good.sum() > 0:
            logger.info('For observation %d (channel=%d) identified %d matching calibrator observations' % (observation_data[i]['obsid'],
                                                                                                            observation_data[i]['channel'],
                                                                                                            good.sum()))
            if good.sum()>1:
                # find the closest in time
                dt=numpy.abs(observation_data[i]['obsid']-observation_data[good]['obsid'])
                closest=observation_data[good][dt==dt.min()]
                logger.info('Will use %d (separation=%d s) for calibration of %s' % (closest['obsid'],
                                                                                     numpy.abs(closest['obsid']-observation_data[i]['obsid']),
                                                                                     observation_data[i]['obsid']))
                cal_observations[observation_data[i]['obsid']]=closest['obsid']
            else:
                logger.info('Will use %d for calibration of %s' % (observation_data[good]['obsid'],
                                                                   observation_data[i]['obsid']))
                cal_observations[observation_data[i]['obsid']]=observation_data[good]['obsid']
                
        else:
            logger.info('For observation %d (channel=%d) identified no matching calibrator observations' % (observation_data[i]['obsid'],
                                                                                                            observation_data[i]['channel']))
            cal_observations[observation_data[i]['obsid']]=None
                             

    return calibrators, notcalibrators, cal_observations


######################################################################
class Observation(metadata.MWA_Observation):

    ##############################
    def __init__(self, obsid, 
                 caltype='anoko',
                 outputdir='./', ncpus=32, 
                 memfraction=50,
                 clobber=False, delete=False):
        """
        Observation(obsid, 
        caltype='anoko',
        outputdir='./', ncpus=32, 
        memfraction=50,
        clobber=False, delete=False)
        """
        assert caltype in ['anoko','casa']

        self.obsid=obsid
        self.basedir=os.path.abspath(outputdir)
        self.clobber=clobber
        self.delete=delete
        self.ncpus=ncpus
        self.memfraction=memfraction

        self.metafits=None
        self.rawfiles=[]
        self.beamfiles=[]
        self.corrfiles=[]

        # this is the obsid of the calibrator to be used                
        self.calibrator=None
        # this is the model file used for anoko calibration
        self.calmodelfile=None
        # this is the type of calibration used (anoko, casa)
        self.caltype=caltype
        # this is the actual calibration (gain) file
        self.calibratorfile=None
        self.calminuv=None
        self.inttime=0
        self.chanwidth=0
        self.nchans=0
        self.reffreq=0
        self.otherkeys={}
        self.autoprocesssources='None'
        self.centerchanged=False
        self.selfcal=False
        self.clean_iterations_selfcal=4000
        self.nimages=1
        self.clean_threshold=None
        self.sources=[]
        self.aegean_region=None

        self.mincpus=1
        self.minmem=3

        if os.path.exists(os.path.join(self.basedir,'%s.metafits' % self.obsid)):
            self.metafits=os.path.join(self.basedir,'%s.metafits' % self.obsid)
        else:
            self.metafits=makemetafits(self.obsid, directory=self.basedir)
            
        self.observation=metadata.MWA_Observation(self.obsid)

        if not os.path.exists('%s.ms' % self.obsid):
            logger.error('MS file %s.ms does not exist' % self.obsid)
        result=get_msinfo('%s.ms' % self.obsid)
        if result is None:
            logger.error('Cannot get info for MS file %s.ms' % self.obsid)
        else:
            self.chanwidth, self.inttime, self.nchans, self.reffreq, self.otherkeys=result
        self.filestodelete=[]

    ##############################
    def __del__(self):
        """
        This desctructor deletes temporary files
        if self.delete is True
        """
        if self.delete:
            for file in self.filestodelete:
                if os.path.exists(file):
                    logger.debug('Deleting %s' % file)
                    try:
                        if os.path.isdir(file):
                            shutil.rmtree(file)
                        else:
                            os.remove(file)
                    except Exception,e:
                        logger.error('Unable to delete file %s:\n\t%s' % (file,e))


    ##############################
    def __getattr__(self, name):
        """
        overloads getattr to lookup info in the self.observation
        class if needed
        """
        if self.__dict__.has_key(name):
            return self.__dict__[name]
        else:
            return self.observation.__dict__[name]

    ##############################
    def checkforlock(self):
        checkforlock('%s.ms' % self.obsid)

    ##############################
    def make_cal(self, minuv=60, selfcal=False):
        """
        make_cal(self, minuv=60, selfcal=False)
        make a calibration solution
        returns the name of the calibration solution on success
        or None on failure
        also sets self.calibratorfile to that file
        """
        assert self.caltype in ['anoko','casa']
        if not selfcal:
            if not self.calibration:
                return None
            if self.caltype=='anoko' and self.calibratorsource in badanokosources:
                logger.warning('Calibrator source %s does not work well with Anoko; switching to CASA calibration...' % self.calibratorsource)
                self.caltype='casa'

            if self.caltype=='anoko':
                calibrator_name=self.calibratorsource
                self.calmodelfile=write_calmodelfile(calibrator_name, directory=self.basedir)
                if self.calmodelfile is None:
                    return None
                self.filestodelete.append(self.calmodelfile)
            self.calibratorfile=os.path.join(self.basedir,'%s.cal' % self.obsid)
            if os.path.exists(self.calibratorfile):
                if self.clobber:
                    logger.warning('Calibration file %s exists but clobber=True; deleting...' % self.calibratorfile)
                    if os.path.isdir(self.calibratorfile):
                        shutil.rmtree(self.calibratorfile)
                    else:
                        os.remove(self.calibratorfile)
                else:
                    logger.warning('Calibration file %s exists and clobber=False; continuing...' % self.calibratorfile)
                    return self.calibratorfile
        if self.caltype=='anoko':
            if not selfcal:
                result=calibrate_anoko(self.obsid,
                                       self.calmodelfile,
                                       minuv=minuv,    
                                       directory=self.basedir,
                                       ncpus=self.ncpus)
            else:
                result=calibrate_anoko(self.obsid,
                                       suffix='selfcal',
                                       corrected=True,
                                       minuv=minuv,                                   
                                       directory=self.basedir,
                                       ncpus=self.ncpus)

            self.calminuv=minuv
        elif self.caltype=='casa':        
            if not selfcal:
                result=calibrate_casa(self.obsid, directory=self.basedir, minuv=minuv)
            else:
                result=selfcalibrate_casa(self.obsid,
                                          suffix='selfcal',
                                          directory=self.basedir,
                                          minuv=minuv)

            self.calminuv=minuv
        if result is not None:
            self.calibratorfile=result
            logger.info('Wrote %s' % self.calibratorfile)
            return self.calibratorfile
        else:
            return None

    ##############################
    def calibrate(self, 
                  selfcal=False,
                  recalibrate=True):
        """
        calibrate(self, selfcal=False,
        recalibrate=True)
        apply a calibration solution
        return True on success, False on error
        """
        assert self.caltype in ['anoko','casa']
        #if self.caltype=='casa' and selfcal:
        #    logger.error('Cannot do selfcal with CASA calibration')
        #    return False
        if not selfcal:
            if check_calibrated('%s.ms' % self.obsid) and not recalibrate:
                logger.info('%s.ms is already calibrated and recalibrate is False' % self.obsid)
                return True
        if self.calibratorfile is None:
            logger.error('No calibration file present for %s.ms' % self.obsid)
            return False
                
        if not os.path.exists(self.calibratorfile):
            logger.error('Cannot find calibration file %s for observation %s' % (self.calibratorfile,
                                                                                 self.obsid))
            return False
        logger.info('Will calibrate %s.ms with %s' % (self.obsid,
                                                      self.calibratorfile))

        if self.caltype=='anoko':
            if os.path.isdir(self.calibratorfile):
                logger.warning('Calibrator file %s appears to be CASA output but Anoko calibration is selected; using CASA applycal...' % self.calibratorfile)
                self.caltype='casa'
        elif self.caltype=='casa':
            if not os.path.isdir(self.calibratorfile):
                logger.warning('Calibrator file %s appears to be Anoko output but CASA calibration is selected; using Anoko applycal...' % self.calibratorfile)
                self.caltype='anoko'


        if self.caltype=='anoko':
            if not selfcal:
                result=applycal_anoko(self.obsid,self.calibratorfile)
            else:
                result=applycal_anoko(self.obsid,self.calibratorfile,
                                      corrected=True)
        elif self.caltype=='casa':
            result=applycal_casa(self.obsid,self.calibratorfile)
        if result is None:
            return False

        if not check_calibrated('%s.ms' % self.obsid):
            logger.error('%s.ms does not appear to be calibrated; no CORRECTED_DATA file is present' % self.obsid)
            return False
        return True


    ##############################
    def autoprocess(self,imagesize=2048,
                    pixelscale=0.015,
                ):
        """
        autoprocess(self,imagesize=2048,
                    pixelscale=0.015,
                ):
        run autoprocess, peeling only
        if the source to be peeled is within the primary FOV 
        (as determined by pixelscale and imagesize) then will not be peeled
        returns True on success
        also sets self.autoprocesssources to the source peeled
        """
        autoprocesscommand=['autoprocess',
                            '-noselfcal',
                            anokocatalog,
                            '%s.ms' % self.obsid]
        logger.info('Will run:\n\t%s' % ' '.join(autoprocesscommand))
        p=subprocess.Popen(' '.join(autoprocesscommand),
                           stderr=subprocess.PIPE,
                           stdout=subprocess.PIPE,
                           shell=True,
        close_fds=False)
        output=[]
        while True:
            p.stdout.flush()
            p.stderr.flush()
            for l in p.stdout.readlines():
                output.append(l)
                logger.debug(l.rstrip())
            for l in p.stderr.readlines():
                if 'dUTC(Double)' in l:
                    continue
                logger.error(l.rstrip())
            returncode=p.poll()
            if returncode is not None:
                break
            time.sleep(1)
        tasks=[]
        distances={}
        peelsource=None
        for l in output:
            if 'distance' in l:
                source=l.split()[0]
                distance=float(l.split()[-2].replace('distance=',''))
                distances[source]=distance
            if l.startswith('- '):
                # if not '- No' in l:
                if 'Peel out source' in l:
                    tasks.append(l[2:])
                    peelsource=l.split()[-1]
    
        if len(tasks)==0:
            logger.info('No autoprocess required')
            self.autoprocesssources='None'
            return None
        for i in xrange(len(tasks)):
            logger.info('autoprocess will do: %s' % tasks[i])
            if distances[peelsource] < (imagesize/2)*pixelscale:
                logger.info('Actually, %s will be in the FOV so no need to peel' % peelsource)
                self.autoprocesssources='None'
                return None

        autoprocesscommand=['autoprocess',
                            '-noselfcal',
                            '-go',
                            anokocatalog,
                            '%s.ms' % self.obsid]
        logger.info('Will run:\n\t%s' % ' '.join(autoprocesscommand))
        p=subprocess.Popen(' '.join(autoprocesscommand),
                           stderr=subprocess.PIPE,
                           stdout=subprocess.PIPE,
                           shell=True,
        close_fds=False)
        while True:
            p.stdout.flush()
            p.stderr.flush()
            for l in p.stdout.readlines():
                logger.debug(l.rstrip())
            for l in p.stderr.readlines():
                if 'dUTC(Double)' in l:
                    continue
                logger.error(l.rstrip())
            returncode=p.poll()
            if returncode is not None:
                break
            time.sleep(1)
        
        self.autoprocesssources=peelsource
        return True


    ##############################
    def chgcentre(self):
        """
        chgcentre()
        returns True on success
        also sets self.centerchanged to True if run
        """
        chgcentrecommand=['chgcentre',
                          '-minw',
                          '-shiftback',
                          '%s.ms' % self.obsid]
        logger.info('Will run:\n\t%s' % ' '.join(chgcentrecommand))
        p=subprocess.Popen(' '.join(chgcentrecommand),
                           stderr=subprocess.PIPE,
                           stdout=subprocess.PIPE,
                           shell=True,
        close_fds=False)
        while True:
            p.stdout.flush()
            p.stderr.flush()
            for l in p.stdout.readlines():
                logger.debug(l.rstrip())
            for l in p.stderr.readlines():
                logger.error(l.rstrip())
            returncode=p.poll()
            if returncode is not None:
                break
            time.sleep(1)
        
        self.centerchanged=True
        return True


    ##############################
    def image(self, 
              predict=False,
              suffix=None,
              clean_weight='uniform',
              imagesize=2048,
              pixelscale=0.015,
              clean_iterations=100,
              clean_gain=0.1,
              clean_mgain=1.0,
              clean_minuv=0,
              clean_maxuv=0,
              clean_threshold=0,
              fullpolarization=True,
              smallinversion=True,
              cleanborder=5,
              fitbeam=True,
              wsclean_arguments='',
              subbands=1,
              updateheader=True):
        """
        image(prefict=False,
              suffix=None,
              clean_weight='uniform',
              imagesize=2048,
              pixelscale=0.015,
              clean_iterations=100,
              clean_gain=0.1,
              clean_mgain=1.0,
              clean_minuv=0,
              clean_maxuv=0,
              clean_threshold=0,
              fullpolarization=True,
              smallinversion=True,
              cleanborder=5,
              fitbeam=True,
              wsclean_arguments='',
              subbands=1,
              updateheader=True)
        images with wsclean
        creates output like <obsid>_<suffix>-<pol>-image.fits
        returns the list of files on success, None on failure
        
        other files (model, etc) are put into list of files to delete
        sets:
        self.clean_weight
        self.imagesize
        self.pixelscale
        self.clean_iterations
        self.clean_gain
        self.clean_mgain
        self.clean_minuv
        self.clean_maxuv
        self.clean_threshold
        self.fullpolarization
        self.smallinversion
        self.cleanborder
        self.fitbeam
        self.wsclean_arguments

        will update the header of the output if desired to include info from 
        this task and from metafits

        """
        self.clean_weight=clean_weight
        self.imagesize=imagesize
        self.pixelscale=pixelscale
        self.clean_iterations=clean_iterations
        self.clean_gain=clean_gain
        self.clean_mgain=clean_mgain
        self.clean_minuv=clean_minuv
        self.clean_maxuv=clean_maxuv
        self.fullpolarization=fullpolarization
        self.smallinversion=smallinversion
        self.cleanborder=cleanborder
        self.fitbeam=fitbeam
        self.wsclean_arguments=wsclean_arguments
        if clean_threshold is not None and clean_threshold > 0:
            self.clean_threshold=clean_threshold

        # self.checkforlock()

        if suffix is not None:
            name='%s_%s' % (self.obsid,suffix)
        else:
            name=str(self.obsid)
        name=os.path.join(self.basedir, name)
        wscleancommand=['wsclean',
                        '-j',
                        str(self.ncpus),
                        '-mem',
                        str(self.memfraction),
                        '-weight',
                        self.clean_weight,
                        '-size',
                        '%d %d'% (self.imagesize,self.imagesize),
                        '-scale',
                        '%.4fdeg' % self.pixelscale]
        if clean_threshold is not None and clean_threshold > 0:
            wscleancommand+=['-threshold',
                             str(clean_threshold)]
        if subbands > 1:
            wscleancommand+=['-channelsout',
                             str(subbands)]
        if self.clean_iterations>0:
            wscleancommand+=['-niter',
                             str(self.clean_iterations),
                             '-gain',
                             str(self.clean_gain),
                             '-mgain',
                             str(self.clean_mgain)]
        if smallinversion:
            wscleancommand+=['-smallinversion']
        else:
            wscleancommand+=['-nosmallinversion']
        wscleancommand+=['-cleanborder %d' % self.cleanborder]
        if fitbeam:
            wscleancommand+=['-fitbeam']
        else:
            wscleancommand+=['-nofitbeam']
        if predict:
            wscleancommand.append('-predict')        
        wscleancommand+=['-name',
                         name]
        expected_output=[]
        if self.clean_minuv > 0:
            wscleancommand.append('-minuv-l')
            wscleancommand.append(str(self.clean_minuv))
        if self.clean_maxuv > 0:
            wscleancommand.append('-maxuv-l')
            wscleancommand.append(str(self.clean_maxuv))
        if self.fullpolarization:
            wscleancommand.append('-pol')
            wscleancommand.append('xx,xy,yx,yy')
            if not predict:
                wscleancommand.append('-joinpolarizations')
        else:
            wscleancommand.append('-pol')
            wscleancommand.append('xx,yy')
        if len(self.wsclean_arguments)>0:
            wscleancommand+=self.wsclean_arguments.split()
        wscleancommand.append('%s.ms' % self.obsid)
        
        if subbands==1:
            if self.fullpolarization:
                for pol in ['XX','YY','XY','XYi']:
                    expected_output.append('%s-%s-image.fits' % (name,pol))
            else:
                for pol in ['XX','YY']:
                    expected_output.append('%s-%s-image.fits' % (name,pol))                
        else:
            for subband in xrange(subbands):
                if self.fullpolarization:
                    for pol in ['XX','YY','XY','XYi']:
                        expected_output.append('%s-%04d-%s-image.fits' % (name,subband,pol))
                else:
                    for pol in ['XX','YY']:
                        expected_output.append('%s-%04d-%s-image.fits' % (name,subband,pol))                
            if self.fullpolarization:
                for pol in ['XX','YY','XY','XYi']:
                    expected_output.append('%s-MFS-%s-image.fits' % (name,pol))
            else:
                for pol in ['XX','YY']:
                    expected_output.append('%s-MFS-%s-image.fits' % (name,pol))                
            

        if predict:
            expected_output=[]

        self.wscleancommand=wscleancommand
        logger.info('Will run:\n\t%s' % ' '.join(wscleancommand))
        p=subprocess.Popen(' '.join(wscleancommand),
                           stderr=subprocess.PIPE,
                           stdout=subprocess.PIPE,
                           shell=True,
                           close_fds=True)
        while True:
            s=non_block_read(p.stdout)
            if len(s)>0:
                logger.debug(s)
            s=non_block_read(p.stderr)
            if len(s)>0:
                logger.error(s)
            returncode=p.poll()
            if returncode is not None:
                break
            time.sleep(1)


        for j in xrange(len(expected_output)):
            file=expected_output[j]
            if not os.path.exists(file):
                logger.error('Expected wsclean output %s does not exist' % file)
                return None
            logger.info('Created wsclean output %s' % file)
            # add other output to cleanup list
            for ext in ['dirty','model','residual']:
                self.filestodelete.append(file.replace('-image','-%s' % ext))

            self.rawfiles.append(file)

            if not updateheader:
                continue
            try:
                f=fits.open(file,'update')
            except Exception,e:
                logger.error('Unable to open image output %s for updating:\n\t%s' % (file,e))
                return None

            f[0].header['CALTYPE']=(self.caltype,'Calibration type')
            #f[0].header['CALIBRAT']=(observation_data[i]['calibrator'],
            #                         'Calibrator observation')
            f[0].header['CALFILE']=(self.calibratorfile,
                                    'Calibrator file')
            f[0].header['CALMINUV']=(self.calminuv,
                                     '[m] Min UV for calibration')
            f[0].header['AUTOPEEL']=self.autoprocesssources
            f[0].header['CHGCENTR']=self.centerchanged
            f[0].header['SELFCAL']=self.selfcal
            med,rms=stat_measure(f)
            if numpy.isnan(rms):
                logger.error('RMS of %s is NaN' % file)
                return None

            f[0].header['IMAGERMS']=(rms,'[Jy/beam] Image RMS')
            try:
                f[0].header.add_history(' '.join(self.wscleancommand))
            except:
                f[0].header.add_history(self.wscleancommand)
            try:
                fm=fits.open(self.metafits)
            except Exception,e:
                logger.open('Unable to open metafits file %s:\n\t%s' % (self.metafits,
                                                                        e))
                return None
            for k in fm[0].header.keys():
                try:
                    f[0].header[k]=(fm[0].header[k],
                                    fm[0].header.comments[k])
                except:
                    pass
            prev_inttime=f[0].header['INTTIME']
            prev_finechan=f[0].header['FINECHAN']
            prev_nchan=f[0].header['NCHANS']
            f[0].header['INTTIME']=(self.inttime,'[s] Integration time')
            f[0].header['FINECHAN']=(self.chanwidth,'[kHz] Fine channel width')
            f[0].header['NCHANS']=(self.nchans,
                                   'Number of fine channels in spectrum')
            if prev_inttime != f[0].header['INTTIME']:
                f[0].header['NSCANS']=int(f[0].header['NSCANS']*prev_inttime/f[0].header['INTTIME'])
                logger.info('New integration time (%.1f s) differs from previous integration time (%.1f s): changing number of scans to %d' % (f[0].header['INTTIME'],
                                                                                                                                               prev_inttime,
                                                                                                                                               f[0].header['NSCANS']))
                
            f[0].header['BANDWDTH']=(self.nchans*self.chanwidth/1000,
                                     '[MHz] Total bandwidth')
            f[0].header['FREQCENT']=(self.reffreq/1e6,
                                     '[MHz] Center frequency of observation')
            if subbands > 1:
                try:
                    f[0].header['FREQCENT']=f[0].header['CRVAL3']/1e6
                    f[0].header['BANDWDTH']=(f[0].header['CDELT3']/1e6,'[MHz] Total bandwidth')
                    f[0].header['NCHANS']=(int(round(f[0].header['CDELT3']/self.chanwidth)),'Number of fine channels in spectrum')
                    # startchannel=f[0].header['WSCCHANS']
                    # stopchannel=f[0].header['WSCCHANE']                    
                    # f[0].header['NCHANS']=(int(stopchannel-startchannel), 'Number of fine channels in spectrum')
                    # f[0].header['BANDWDTH']=((stopchannel-startchannel)*self.chanwidth/1e3, 
                    #                         '[MHz] Total bandwidth')
                    
                    # f[0].header['FREQCENT']=f[0].header['FREQCENT']+(self.chanwidth/1e3)*((stopchannel-startchannel)/2-self.nchans/2)
                except:
                    pass
            f.verify('fix')
            f.flush()

        if not predict:
            self.filestodelete.append('%s-psf.fits' % name)
            return self.rawfiles
        else:
            return True


    ##############################
    def multiimage_time(self, 
                        suffix=None,
                        imagetime=2,
                        clean_weight='uniform',
                        imagesize=2048,
                        pixelscale=0.015,
                        clean_iterations=100,
                        clean_gain=0.1,
                        clean_mgain=1.0,
                        clean_minuv=0,
                        clean_maxuv=0,
                        clean_threshold=0,
                        fullpolarization=True,
                        wsclean_arguments='',
                        updateheader=True,
                        ntorun=8):
        """
        multiimage_time(self, 
                        suffix=None,
                        imagetime=2,
                        clean_weight='uniform',
                        imagesize=2048,
                        pixelscale=0.015,
                        clean_iterations=100,
                        clean_gain=0.1,
                        clean_mgain=1.0,
                        clean_minuv=0,
                        clean_maxuv=0,
                        clean_threshold=0,
                        fullpolarization=True,
                        wsclean_arguments='',
                        updateheader=True,
                        ntorun=8)
        images with wsclean
        creates output like <obsid>_<suffix>_t<i>-<j>-<pol>-image.fits
        does this for many sub-images in parallel
        returns the list of files on success, None on failure
        
        other files (model, etc) are put into list of files to delete
        sets:
        self.clean_weight
        self.imagesize
        self.pixelscale
        self.clean_iterations
        self.clean_gain
        self.clean_mgain
        self.clean_minuv
        self.clean_maxuv
        self.fullpolarization
        self.wsclean_arguments
        self.clean_threshold

        will update the header of the output if desired to include info from 
        this task and from metafits

        """
        self.clean_weight=clean_weight
        self.imagesize=imagesize
        self.pixelscale=pixelscale
        if clean_threshold is not None and clean_threshold > 0:
            self.clean_threshold=clean_threshold
        self.clean_iterations=clean_iterations
        self.clean_gain=clean_gain
        self.clean_mgain=clean_mgain
        self.clean_minuv=clean_minuv
        self.clean_maxuv=clean_maxuv
        self.fullpolarization=fullpolarization
        self.wsclean_arguments=wsclean_arguments

        if suffix is not None:
            name0='%s_%s' % (self.obsid,suffix)
        else:
            name0=str(self.obsid)
        name0=os.path.join(self.basedir, name0)
        
        if imagetime < self.inttime:
            logger.warning('Requested image time of %.1f s is < integration time of %.1f s; setting to integration time' % (imagetime,
                                                                                                                            self.inttime))
            imagetime=self.inttime

        # determine number of images to make
        nimages = int(self.duration / imagetime)
        self.nimages=nimages
        integrationsperimage=int(imagetime/self.inttime)
        logger.info('Will make %d images of %.1f s each' % (nimages, imagetime))
        commands=[]
        expected_outputs=[]
        tempdirs=[]
        for i in xrange(nimages):
            expected_output=[]
            startindex=i*integrationsperimage
            stopindex=(i+1)*integrationsperimage
            if stopindex > int(self.duration / self.inttime):
                stopindex=int(self.duration / self.inttime)
            name='%s_t%03d-%03d' % (name0, startindex, stopindex)

            ncpus=max(self.mincpus,self.ncpus/ntorun)
            mem=max(self.minmem,self.memfraction/ntorun)
            wscleancommand=['wsclean',
                            '-j',
                            str(ncpus),
                            '-mem',
                            str(mem),
                            '-weight',
                            self.clean_weight,
                            '-size',
                            '%d %d'% (self.imagesize,self.imagesize),
                            '-scale',
                            '%.4fdeg' % self.pixelscale,
                            '-interval',
                            '%d %d' % (startindex, stopindex)]
            if clean_threshold is not None and clean_threshold > 0:
                wscleancommand+=['-threshold',
                                 str(clean_threshold)]
            if self.clean_iterations>0:
                wscleancommand+=['-niter',
                                 str(self.clean_iterations),
                                 '-gain',
                                 str(self.clean_gain),
                                 '-mgain',
                                 str(self.clean_mgain)]
                if self.clean_mgain < 1 and ntorun>1:
                    logger.warning('Unable to run multiple wscleans in parallel because major cycles are requested.  Setting nprocess to 1...')
                    #ntorun=1
            # if we run multiple in parallel, 
            # need to use tempdirs for reordering
            if ntorun > 1:
                #wscleancommand+=['-no-reorder']
                tempdir=tempfile.mkdtemp(dir=self.basedir)
                tempdirs.append(tempdir)
                wscleancommand+=['-tempdir',
                                 tempdir,
                                 '-no-update-model-required']

            wscleancommand+=['-name',
                             name]
            if self.clean_minuv > 0:
                wscleancommand.append('-minuv-l')
                wscleancommand.append(str(self.clean_minuv))
            if self.clean_maxuv > 0:
                wscleancommand.append('-maxuv-l')
                wscleancommand.append(str(self.clean_maxuv))
            if self.fullpolarization:
                wscleancommand.append('-pol')
                wscleancommand.append('xx,xy,yx,yy')
                for pol in ['XX','YY','XY','XYi']:
                    expected_output.append('%s-%s-image.fits' % (name,pol)) 
            else:
                wscleancommand.append('-pol')
                wscleancommand.append('xx,yy')
                for pol in ['XX','YY']:
                    expected_output.append('%s-%s-image.fits' % (name,pol)) 
            if len(self.wsclean_arguments)>0:
                wscleancommand+=self.wsclean_arguments.split()
            wscleancommand.append('%s.ms' % self.obsid)
            commands.append(' '.join(wscleancommand))
            expected_outputs.append(expected_output)
            self.wscleancommand=wscleancommand

        p=[None]*ntorun
        o=[None]*ntorun
        e=[None]*ntorun
        i=0
        idone=0
        while (i<len(commands) or idone<len(commands)):
            if (None in p):
                for j in xrange(len(p)):
                    if (p[j] is None and i < len(commands)):
                        logger.info("Will run (%d/%d) in %d:\n\t%s" % (i,len(commands),j,commands[i]))
                        op,ep=subprocess.PIPE,subprocess.PIPE
                        p[j]=subprocess.Popen(commands[i], shell=True,
                                              stderr=ep,
                                              stdout=op,
                                              close_fds=True)
                        i+=1
            poll=getstatus(p,output=o,error=e)

            while (not (any(poll))):
                time.sleep(5)
                poll=getstatus(p,output=o,error=e)

            for j in xrange(len(p)):
                if (poll[j]):
                    idone+=1
                    logger.info("\t%d finished" % j)
                    p[j]=None

        for i in xrange(len(expected_outputs)):
            expected_output=expected_outputs[i]
            for j in xrange(len(expected_output)):
                file=expected_output[j]
                if not os.path.exists(file):
                    logger.error('Expected wsclean output %s does not exist' % file)
                    return None
                logger.info('Created wsclean output %s' % file)
                # add other output to cleanup list
                for ext in ['dirty','model','residual']:
                    self.filestodelete.append(file.replace('-image','-%s' % ext))

                self.rawfiles.append(file)

                if not updateheader:
                    continue
                try:
                    f=fits.open(file,'update')
                except Exception,e:
                    logger.error('Unable to open image output %s for updating:\n\t%s' % (file,e))
                    return None

                f[0].header['CALTYPE']=(self.caltype,'Calibration type')
                #f[0].header['CALIBRAT']=(observation_data[i]['calibrator'],
                #                         'Calibrator observation')
                f[0].header['CALFILE']=(self.calibratorfile,
                                        'Calibrator file')
                f[0].header['CALMINUV']=(self.calminuv,
                                         '[m] Min UV for calibration')
                f[0].header['AUTOPEEL']=self.autoprocesssources
                f[0].header['CHGCENTR']=self.centerchanged
                f[0].header['SELFCAL']=self.selfcal
                med,rms=stat_measure(f)
                if numpy.isnan(rms):
                    logger.error('RMS of %s is NaN' % file)
                    return None

                f[0].header['IMAGERMS']=(rms,'[Jy/beam] Image RMS')

                f[0].header.add_history(commands[i])
                try:
                    fm=fits.open(self.metafits)
                except Exception,e:
                    logger.open('Unable to open metafits file %s:\n\t%s' % (self.metafits,
                                                                        e))
                    return None
                for k in fm[0].header.keys():
                    if k in ['DATE-OBS']:
                        continue
                    try:                        
                        f[0].header[k]=(fm[0].header[k],
                                        fm[0].header.comments[k])
                    except:
                        pass

                prev_inttime=f[0].header['INTTIME']
                prev_finechan=f[0].header['FINECHAN']
                prev_nchan=f[0].header['NCHANS']   

                f[0].header['INTTIME']=(self.inttime,'[s] Integration time')
                f[0].header['FINECHAN']=(self.chanwidth,'[kHz] Fine channel width')
                f[0].header['NCHANS']=(self.nchans,
                                       'Number of fine channels in spectrum')

                f[0].header['NSCANS']=int(integrationsperimage)
                f[0].header['EXPOSURE']=(integrationsperimage*self.inttime,'[s] duration of observation')
                f[0].header['BANDWDTH']=(self.nchans*self.chanwidth/1000,
                                         '[MHz] Total bandwidth')
                f[0].header['FREQCENT']=(self.reffreq/1e6,
                                         '[MHz] Center frequency of observation')
            
                f.verify('fix')
                f.flush()

        if len(tempdirs)>0:
            for tempdir in tempdirs:
                shutil.rmtree(tempdir)

        return self.rawfiles







    ##############################
    def makepb(self, suffix=None, ifile=0, beammodel='2014i'):
        """
        makepb(self, suffix=None, ifile=0, beammodel='2014i')
        make primary beam using anoko/beam
        
        creates files like beam-<obsid>_<suffix>-<pol>.fits
        
        returns list of new beam files on success, None on failure
        """
        assert beammodel in ['2013','2014','2014i']

        name=self.rawfiles[ifile].replace('-XX-image.fits','').replace('-XX-image.fits','').replace('-XY-image.fits','').replace('-XYi-image.fits','').replace('-YY-image.fits','')
        if '.fits' in name:
            logger.error('Cannot construct name base from output %s' % self.rawfiles[ifile])
            return None
        name=os.path.split(name)[1]
        if suffix is not None:
            name+='_%s' % suffix
        beamname=os.path.join(self.basedir,'beam-%s' % name)
        name=os.path.join(self.basedir, name)

        self.beammodel=beammodel
        beamcommand=['beam',
                     '-' + beammodel,
                     '-name',
                     beamname,
                     '-proto',
                     self.rawfiles[ifile],
                     '-ms',
                     '%s.ms' % self.obsid]
        
        expected_beamoutput=[]
        for pol in ['xxr','xxi','yyr','yyi','xyr','xyi','yxr','yxi']:
            expected_beamoutput.append('%s-%s.fits' % (beamname,
                                                       pol))

        remake=False
        # check to see if the expected file already exists
        # if so, make sure it has the right beammodel 
        # and size/position/frequency
        for file in expected_beamoutput:
            if not os.path.exists(file):
                remake=True
            else:
                fb=fits.open(file)
                f=fits.open(self.rawfiles[ifile])
                if not 'PBMODEL' in fb[0].header.keys() and fb[0].header['PBMODEL']==self.beammodel:
                    remake=True
                keys=['NAXIS1','NAXIS2',
                      'CRVAL1','CRVAL2',
                      'CDELT1','CDELT2',
                      'CRVAL3']
                for k in keys:
                    if fb[0].header[k] != f[0].header[k]:
                        remake=True
        if not remake:
            if not self.clobber:
                logger.info('All expected beam outputs already exist and match requirements')
                return expected_beamoutput
            else:
                logger.info('All expected beam outputs already exist and match requirements but clobber=True; remaking...')
                remake=True

        logger.info('Will run:\n\t%s' % ' '.join(beamcommand))
        p=subprocess.Popen(' '.join(beamcommand),
                           stderr=subprocess.PIPE,
                           stdout=subprocess.PIPE,
                           shell=True,
                           close_fds=False)
        while True:
            p.stdout.flush()
            p.stderr.flush()
            for l in p.stdout.readlines():
                logger.debug(l.rstrip())
            for l in p.stderr.readlines():
                # I get a lot of errors about UTC-UT being
                # possible wrong, but I don't think it matters
                # here
                if 'dUTC(Double)' in l:
                    continue
                logger.error(l.rstrip())
            returncode=p.poll()
            if returncode is not None:
                break
            time.sleep(1)

        for j in xrange(len(expected_beamoutput)):
            if not os.path.exists(expected_beamoutput[j]):
                logger.error('Expected beam output %s does not exist' % expected_beamoutput[j])
                return None
            
            self.filestodelete.append(expected_beamoutput[j])
            file=expected_beamoutput[j]
            try:
                f=fits.open(file,'update')
            except Exception,e:
                logger.open('Unable to open PB output %s for updating:\n\t%s' % (file,e))
                return None            
            f[0].header['PBMODEL']=(self.beammodel,'Primary beam model')
            f.verify('fix')
            f.flush()
            self.beamfiles.append(expected_beamoutput[j])

        return expected_beamoutput

    ##############################
    def pbcorrect(self, suffix=None, model=False, uncorrect=False, ifile=0):
        """
        pbcorrect(self, suffix=None, model=False, uncorrect=False, ifile=0)
        primary beam correct using anoko/pbcorrect

        creates output like <obsid>_<suffix>-<pol>.fits
        if doing fullpol, will make all Stokes
        otherwise Q,U,V get put into list of files to delete
        
        returns list of corrected images on success, None on failure
        """


        name=self.rawfiles[ifile].replace('-XX-image.fits','').replace('-XY-image.fits','').replace('-XYi-image.fits','').replace('-YY-image.fits','')
        if '.fits' in name:
            logger.error('Cannot construct name base from output %s' % self.rawfiles[ifile])
            return None
        name=os.path.split(name)[1]
        if suffix is not None:
            name+='_%s' % suffix
        beamname=os.path.join(self.basedir,'beam-%s' % name)
        name=os.path.join(self.basedir, name)
        outname=name
        if not model:
            imagetype='image.fits'
        else:
            imagetype='model.fits'
            outname+='_model'

        pbcorrectcommand=['pbcorrect']
        if uncorrect:
            pbcorrectcommand.append('-uncorrect')
            name+='_model'
        pbcorrectcommand+=[name,
                           imagetype,
                           beamname,
                           outname]
        expected_pbcorroutput=[]
        if uncorrect:
            for pol in ['XX','YY','XY','XYi']:
                expected_pbcorroutput.append('%s-%s-%s' % (name,
                                                           pol,
                                                           imagetype))
        elif self.fullpolarization:
            for pol in ['I','Q','U','V']:
                expected_pbcorroutput.append('%s-%s.fits' % (outname,
                                                             pol))
        else:
            for pol in ['I']:
                expected_pbcorroutput.append('%s-%s.fits' % (outname,
                                                             pol))

        logger.info('Will run:\n\t%s' % ' '.join(pbcorrectcommand))
        p=subprocess.Popen(' '.join(pbcorrectcommand),
                           stderr=subprocess.PIPE,
                           stdout=subprocess.PIPE,
                           shell=True,
                           close_fds=False)
        while True:
            p.stdout.flush()
            p.stderr.flush()
            for l in p.stdout.readlines():
                logger.debug(l.rstrip())
            for l in p.stderr.readlines():
                logger.error(l.rstrip())
            returncode=p.poll()
            if returncode is not None:
                break
            time.sleep(1)

        for j in xrange(len(expected_pbcorroutput)):
            if not os.path.exists(expected_pbcorroutput[j]):
                logger.error('Expected pbcorrect output %s does not exist' % expected_pbcorroutput[j])
                return None
            logger.info('Created pbcorrect output %s' % expected_pbcorroutput[j])
            self.corrfiles.append(expected_pbcorroutput[j])
            file=expected_pbcorroutput[j]
            try:
                f=fits.open(file,'update')
            except Exception,e:
                logger.error('Unable to open image output %s for updating:\n\t%s' % (file,e))
                return None
            f[0].header['PBMODEL']=(self.beammodel,'Primary beam model')


            # update header from raw file
            try:
                fraw=fits.open(self.rawfiles[ifile])
            except Exception,e:
                logger.error('Unable to open file %s:\n\t%s' % (self.rawfiles[ifile],e))

            med,rms=stat_measure(f)
            if numpy.isnan(rms):
                logger.error('RMS of %s is NaN' % file)
                return None
            f[0].header['IMAGERMS']=(rms,'[Jy/beam] Image RMS')
            for k in fraw[0].header.keys():
                if not k in f[0].header.keys():
                    f[0].header[k]=(fraw[0].header[k],
                                    fraw[0].header.comments[k])
            f.verify('fix')
            f.flush()

        if not self.fullpolarization:
            for pol in ['Q','U','V']:
                if os.path.exists('%s-%s.fits' % (outname,pol)):
                    self.filestodelete.append('%s-%s.fits' % (outname,pol))


            
        return self.corrfiles

    ##############################
    def multipbcorrect_time(self, suffix=None):
        """
        pbcorrect(self, suffix=None)
        primary beam correct using anoko/pbcorrect

        creates output like <obsid>_<suffix>-<pol>.fits
        if doing fullpol, will make all Stokes
        otherwise Q,U,V get put into list of files to delete
        
        returns list of corrected images on success, None on failure
        """
        
        beamname=None
        for i in xrange(self.nimages):
            if self.fullpolarization:
                rawfile=self.rawfiles[i*4]
            else:
                rawfile=self.rawfiles[i*2]

            name=rawfile.replace('-XX-image.fits','').replace('-XY-image.fits','').replace('-XYi-image.fits','').replace('-YY-image.fits','')
            if '.fits' in name:
                logger.error('Cannot construct name base from output %s' % rawfile)
                return None
            name=os.path.split(name)[1]
            if suffix is not None:
                name+='_%s' % suffix
            if beamname is None:
                beamname=os.path.join(self.basedir,'beam-%s' % name)
            name=os.path.join(self.basedir, name)
            imagetype='image.fits'

            pbcorrectcommand=['pbcorrect']
            pbcorrectcommand+=[name,
                               imagetype,
                               beamname,
                               name]
            expected_pbcorroutput=[]
            if self.fullpolarization:
                for pol in ['I','Q','U','V']:
                    expected_pbcorroutput.append('%s-%s.fits' % (name,
                                                                 pol))
            else:
                for pol in ['I']:
                    expected_pbcorroutput.append('%s-%s.fits' % (name,
                                                                 pol))

            logger.info('Will run:\n\t%s' % ' '.join(pbcorrectcommand))
            p=subprocess.Popen(' '.join(pbcorrectcommand),
                               stderr=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               shell=True,
                               close_fds=False)
            while True:
                p.stdout.flush()
                p.stderr.flush()
                for l in p.stdout.readlines():
                    logger.debug(l.rstrip())
                for l in p.stderr.readlines():
                    logger.error(l.rstrip())
                returncode=p.poll()
                if returncode is not None:
                    break
                time.sleep(1)

            for j in xrange(len(expected_pbcorroutput)):
                if not os.path.exists(expected_pbcorroutput[j]):
                    logger.error('Expected pbcorrect output %s does not exist' % expected_pbcorroutput[j])
                    return None
                logger.info('Created pbcorrect output %s' % expected_pbcorroutput[j])
                self.corrfiles.append(expected_pbcorroutput[j])
                file=expected_pbcorroutput[j]
                try:
                    f=fits.open(file,'update')
                except Exception,e:
                    logger.open('Unable to open image output %s for updating:\n\t%s' % (file,e))
                    return None
                f[0].header['PBMODEL']=(self.beammodel,'Primary beam model')

                # update header from raw file
                try:
                    fraw=fits.open(rawfile)
                except Exception,e:
                    logger.error('Unable to open file %s:\n\t%s' % (rawfile,e))
                    
                for k in fraw[0].header.keys():
                    if not k in f[0].header.keys():
                        f[0].header[k]=(fraw[0].header[k],
                                        fraw[0].header.comments[k])
                f.verify('fix')
                f.flush()

            if not self.fullpolarization:
                for pol in ['Q','U','V']:
                    if os.path.exists('%s-%s.fits' % (name,pol)):
                        self.filestodelete.append('%s-%s.fits' % (name,pol))


            
        return self.corrfiles

            
######################################################################
# from fluxmatch.py
######################################################################
def get_rms_background(imagename, cores=16):
    outbase=os.path.splitext(imagename)[0]
    rmsimage=outbase + '_rms.fits'
    bgimage=outbase + '_bkg.fits'
    if not (os.path.exists(rmsimage) and os.path.exists(bgimage)):
        BANE.filter_image(imagename,
                          outbase,
                          mask=False, 
                          cores=cores)
    else:
        logger.info('Using existing background and rms images %s, %s' % (bgimage,rmsimage))
    return rmsimage, bgimage
    
######################################################################
def find_sources_in_image(imagename, max_summits=5, nsigma=10, usebane=True,
                          region=None, cores=16):
    """
    sources,rmsimage,bgimage=find_sources_in_image(imagename, max_summits=5, nsigma=10, usebane=True, region=None, cores=16)
    runs aegean.find_sources_in_image
    but first runs BANE to get the BG/rms estimates
    if region is supplied (.mim format) only sources inside that will be identified
    """
    if usebane:
        rmsimage, bgimage=get_rms_background(imagename, cores=cores)
    else:
        rmsimage,bgimage=None,None

    sources=source_finder.SourceFinder().find_sources_in_image(imagename,
                                                               max_summits=max_summits,
                                                               innerclip=nsigma,
                                                               rmsin=rmsimage,
                                                               bkgin=bgimage,
                                                               mask=region,
                                                               cores=cores)
    return sources,rmsimage,bgimage

######################################################################
def aegean2table(sources):
    """
    sourcesTable=aegean2table(sources)
    return astropy table corresponding to aegean source list
    """
    ra=Column([s.ra for s in sources],name='RA')
    dec=Column([s.dec for s in sources],name='Dec')
    peakflux=Column([s.peak_flux for s in sources],name='PeakFlux')
    peakfluxerr=Column([s.err_peak_flux for s in sources],name='PeakFluxErr')
    intflux=Column([s.int_flux for s in sources],name='IntFlux')
    intfluxerr=Column([s.err_int_flux for s in sources],name='IntFluxErr')
    rms=Column([s.local_rms for s in sources],name='RMS')

    sourcesTable=Table([ra,dec,peakflux,peakfluxerr,intflux,intfluxerr,rms])
    return sourcesTable
######################################################################
def find_beam(image):
    f=fits.open(image)
    if 'BEAM' in f[0].header.keys():
        if os.path.exists(f[0].header['BEAM']):
            logger.info('Found existing primary beam in header: %s' % f[0].header['BEAM'])
            return f[0].header['BEAM']
    if not 'DELAYS' in f[0].header.keys():
        return None
    delays=[int(x) for x in f[0].header['DELAYS'].split(',')]
    out=make_beam.make_beam(image, ext=0,
                            delays=delays,
                            model='analytic',
                            jones=False,
                            precess=False)
    if f[0].header['CRVAL4']==-5:
        # XX
        return out[0]
    elif f[0].header['CRVAL4']==-6:
        # YY
        return out[1]
    elif f[0].header['CRVAL4']==1:
        # I
        fXX=fits.open(out[0])
        fYY=fits.open(out[1])
        fXX[0].data=0.5*(fXX[0].data+fYY[0].data)
        fXX[0].header['CRVAL4']=1
        out=out[0].replace('_beamXX','_beamI')
        if os.path.exists(out):
            os.remove(out)
        fXX.writeto(out)
        return out            

######################################################################
def fluxmatch(image,
              catalog='GLEAMIDR3.fits',
              nsigma=10,
              rmsfactor=3,
              matchradius=120,
              rejectsigma=3,
              maxdistance=20,
              minbeam=0.5,
              psfextent=1.1,
              refineposition=False,
              maxscale=10,
              update=False,
              plot=True,
              cores=1):
    """
    catalog='GLEAMIDR3.fits',
    signal-to-noise for source finding
    nsigma=10,
    ratio of local rms to minimum image RMS that is OK
    rmsfactor=3,
    distance between catalog source and source in image (arcsec)
    matchradius=120,
    rejection sigma for flux ratio outliers
    rejectsigma=3,
    distance from the image center that is OK (deg)
    maxdistance=20,
    minimum beam power compared to max
    minbeam=0.5,
    area of source/area of psf threshold
    psfextent=1.1,
    """

    if not isinstance(matchradius,astropy.units.quantity.Quantity):
        matchradius=matchradius*u.arcsec
    if not isinstance(maxdistance,astropy.units.quantity.Quantity):
        maxdistance=maxdistance*u.deg        

    if not _matplotlib:
        logger.warning('matplotlib is not enabled; unable to plot')

    if not os.path.exists(image):
        logger.error('Cannot find input image %s' % image)
        return None
    if not os.path.exists(catalog):
        logger.error('Cannot find GLEAM catalog %s' % catalog)
        return None
    beam=find_beam(image)
    if beam is None:
        logger.warning('Did not generate primary beam: will ignore')
        minbeam=None
    if not os.path.exists(beam):
        logger.warning('Cannot find primary beam %s: will ignore' % beam)
        minbeam=None                
        beam=None
    outbase=os.path.splitext(image)[0]           
    sources, rmsimage, bgimage=find_sources_in_image(image,
                                                     nsigma=nsigma,
                                                     cores=cores)
    logger.info('Found %d sources above %d sigma in %s' % (len(sources),
                                                           nsigma,
                                                           image))
    logger.info('Wrote %s and %s' % (rmsimage, bgimage))
    # convert to astropy table
    sourcesTable=aegean2table(sources)
    
    fimage=fits.open(image)
    frequency=fimage[0].header['CRVAL3']
    w=WCS(fimage[0].header,naxis=2)
    frmsimage=fits.open(rmsimage)
    minrms=frmsimage[0].data.min()
    logger.info('Minimum RMS in image is %.1f mJy' % (minrms*1e3))

    if beam is not None:
        fbeam=fits.open(beam)    

    x,y=w.wcs_world2pix(sourcesTable['RA'],sourcesTable['Dec'],0)
    sourcesTable.add_column(Column(x,name='X'))
    sourcesTable.add_column(Column(y,name='Y'))
    pointingcenter=SkyCoord(fimage[0].header['RA'],fimage[0].header['DEC'],
                            unit=('deg','deg'))
    coords=SkyCoord(sourcesTable['RA'],sourcesTable['Dec'],unit=(u.deg,u.deg))
    sourcesTable.add_column(Column(coords.separation(pointingcenter).to(u.deg),
                                   name='SOURCEDIST'))
    if beam is not None:
        sourcesTable.add_column(Column(fbeam[0].data[0,0,numpy.int16(y),
                                                     numpy.int16(x)],
                                       name='BEAM'))
    else:
        sourcesTable.add_column(Column(0*x,
                                       name='BEAM'))
        
    fcatalog=fits.open(catalog)
    catalogTable=Table(fcatalog[1].data)
    bandfrequencies=numpy.array([int(s.split('_')[-1]) for s in numpy.array(catalogTable.colnames)[numpy.nonzero(numpy.array([('int_flux' in c) and not ('deep' in c) for c in catalogTable.colnames]))[0]]])

    # find the indices of the bands just above and below the observation
    # linearly weight the fluxes just above and below to match
    # the observation frequency
    indexplus=(bandfrequencies>=frequency/1e6).nonzero()[0].min()
    indexminus=(bandfrequencies<frequency/1e6).nonzero()[0].max()
    logger.info('Observation frequency of %.1f MHz: interpolating between %d MHz and %d MHz' % (frequency/1e6,bandfrequencies[indexminus],bandfrequencies[indexplus]))
    
    weightplus=(frequency/1e6-bandfrequencies[indexminus])/(bandfrequencies[indexplus]-bandfrequencies[indexminus])
    weightminus=1-weightplus
    gleamflux=catalogTable['int_flux_%03d' % bandfrequencies[indexminus]]*weightminus+catalogTable['int_flux_%03d' % bandfrequencies[indexplus]]*weightplus
    gleamfluxerr=numpy.sqrt((catalogTable['err_fit_flux_%03d' % bandfrequencies[indexminus]]*weightminus)**2+(catalogTable['err_fit_flux_%03d' % bandfrequencies[indexplus]]*weightplus)**2)

    catalogcoords=SkyCoord(catalogTable['RAJ2000'],
                           catalogTable['DECJ2000'],unit=(u.deg,u.deg))
    # match the catalog to the data
    idx,sep2d,sep3d=coords.match_to_catalog_sky(catalogcoords)
    # add the matched columns to the soure table
    sourcesTable.add_column(Column(catalogTable['Name'][idx],
                                   name='Name'))
    sourcesTable.add_column(Column(catalogTable['RAJ2000'][idx],
                                   name='GLEAMRA'))
    sourcesTable.add_column(Column(catalogTable['DECJ2000'][idx],
                                   name='GLEAMDEC'))
    sourcesTable.add_column(Column(sep2d.to(u.arcsec),
                                   name='GLEAMSep'))
    sourcesTable.add_column(Column(gleamflux[idx],
                                   name='GLEAMFlux'))
    sourcesTable.add_column(Column(gleamfluxerr[idx],
                                   name='GLEAMFluxErr'))
    sourcesTable.add_column(Column(catalogTable['psf_a_%03d' % bandfrequencies[indexplus]][idx] * catalogTable['psf_b_%03d' % bandfrequencies[indexplus]][idx],
                                   name='PSFAREA'))
    sourcesTable.add_column(Column(catalogTable['a_%03d' % bandfrequencies[indexplus]][idx] * catalogTable['b_%03d' % bandfrequencies[indexplus]][idx],
                                   name='SOURCEAREA'))

    dRA=(sourcesTable['RA']-sourcesTable['GLEAMRA'])
    dDEC=(sourcesTable['Dec']-sourcesTable['GLEAMDEC'])
    iterations=1
    if refineposition:
        iterations=2
        

    for iter in xrange(iterations):
        # determine the good matches
        # first criterion is separation
        good=(sourcesTable['GLEAMSep']<matchradius)
        logger.info('%04d/%04d sources are within %.1f arcsec' % (good.sum(),
                                                                  len(good),
                                                                  matchradius.to(u.arcsec).value))
        # only point sources
        if psfextent is not None and psfextent>0:
            good=good & (sourcesTable['SOURCEAREA']<=psfextent*sourcesTable['PSFAREA'])
            logger.info('%04d/%04d sources also have source a*b < %.1f * psf a*b' % (good.sum(),
                                                                                     len(good),
                                                                                     psfextent))
        # cut on the local rms compared to the minimum in the image
        if rmsfactor is not None and rmsfactor>0:
            good=good & (sourcesTable['RMS']<=rmsfactor*minrms)
            logger.info('%04d/%04d sources also have RMS < %.1f mJy' % (good.sum(),
                                                                        len(good),
                                                                        rmsfactor*minrms*1e3))
        # distance from pointing center
        if maxdistance is not None and maxdistance>0:
            good=good & (sourcesTable['SOURCEDIST'] < maxdistance)
            logger.info('%04d/%04d sources also are within %.1f deg of pointing center' % (good.sum(),
                                                                                           len(good),
                                                                                           maxdistance.to(u.deg).value))
        # primary beam power
        if minbeam is not None and minbeam>0:
            good=good & (sourcesTable['BEAM']>minbeam*fbeam[0].data.max())
            logger.info('%04d/%04d sources also are at primary beam power > %.2f' % (good.sum(),len(good),minbeam))

        # require that all sources are > 5 sigma detections
        # and that flux uncertainties are > 0
        good=good & (sourcesTable['IntFluxErr']<0.2*sourcesTable['IntFlux']) & (sourcesTable['IntFluxErr']>0) & (sourcesTable['GLEAMFluxErr']>0) & (sourcesTable['GLEAMFluxErr']<0.2*sourcesTable['GLEAMFlux'])
        logger.info('%04d/%04d sources match all cuts' % (good.sum(),
                                                          len(good)))

        fitres=numpy.polyfit(sourcesTable['GLEAMFlux'][good],
                             sourcesTable['IntFlux'][good],
                             deg=1,
                             w=1/sourcesTable['IntFluxErr'][good]**2)
        ratio=sourcesTable['IntFlux']/sourcesTable['GLEAMFlux']
        ratioerr=numpy.sqrt((sourcesTable['IntFluxErr']/sourcesTable['GLEAMFlux'])**2+(sourcesTable['IntFlux']*sourcesTable['GLEAMFluxErr']/sourcesTable['GLEAMFlux']**2)**2)
        if rejectsigma is not None:
            # do a bit of sigma clipping just in case
            good=(good) & (numpy.abs(ratio-numpy.median(ratio[good]))<=ratioerr*rejectsigma)
        fittedratio=(ratio[good]/ratioerr[good]**2).sum()/(1/ratioerr[good]**2).sum()
        fittedratioerr=numpy.sqrt(1/(1/ratioerr[good]**2).sum())
        chisq=(((ratio[good]-fittedratio)/ratioerr[good])**2).sum()
        ndof=good.sum()-1
        logger.info('Found ratio of %s / %s = %.3f +/- %.3f' % (image,
                                                                catalog,
                                                                fittedratio,
                                                                fittedratioerr))
        if refineposition and iter==0:
            sourcesTable['RA']-=dRA[good].mean()
            sourcesTable['Dec']-=dDEC[good].mean()
            logger.info('Applied shift of (%.1f sec, %.1f arcsec)' % (dRA[good].mean()*3600,
                                                                      dDEC[good].mean()*3600))
            coords=SkyCoord(sourcesTable['RA'],sourcesTable['Dec'],unit=(u.deg,u.deg))
            idx,sep2d,sep3d=coords.match_to_catalog_sky(catalogcoords)
            sourcesTable['GLEAMSep']=sep2d.to(u.arcsec)

    sourcesTable.add_column(Column(good,name='GOOD'))
    sourcesTable.meta['ratio']=fittedratio
    sourcesTable.meta['ratio_err']=fittedratioerr
    sourcesTable.meta['chisq']=chisq
    sourcesTable.meta['ndof']=ndof
    sourcesTable.meta['slope']=fitres[0]
    sourcesTable.meta['intercept']=fitres[1]
    sourcesTable.meta['rashift']=0
    sourcesTable.meta['decshift']=0
    if refineposition:
        sourcesTable.meta['rashift']=dRA[good].mean()*3600
        sourcesTable.meta['decshift']=dDEC[good].mean()*3600
    if os.path.exists(outbase + '_fluxmatch.hdf5'):
        os.remove(outbase + '_fluxmatch.hdf5')
    sourcesTable.write(outbase + '_fluxmatch.hdf5',path='data')
    logger.info('Wrote %s_fluxmatch.hdf5' % outbase)


    if update:
        if fittedratio > maxscale or fittedratio < 1.0/maxscale:
            logger.warning('Ratio exceeds reasonable limits; skipping...')
        else:
            fimage=fits.open(image,'update')
            fimage[0].data/=fittedratio
            fimage[0].header['FLUXSCAL']=(fittedratio,'Flux scaling relative to catalog')
            fimage[0].header['FLUX_ERR']=(fittedratioerr,'Flux scaling uncertainty relative to catalog')
            fimage[0].header['FLUXCAT']=(catalog,'Flux scaling catalog')
            fimage[0].header['NFLUXSRC']=(good.sum(),'Number of sources used for flux scaling')
            fimage[0].header['FLUXCHI2']=(chisq,'Flux scaling chi-squared')
            fimage[0].header['FLUXSLOP']=(fitres[0],'Flux scaling slope')
            if refineposition:
                fimage[0].header['RASHIFT']=(dRA[good].mean()*3600,'[s] RA Shift for catalog match')
                fimage[0].header['DECSHIFT']=(dDEC[good].mean()*3600,'[arcsec] DEC Shift for catalog match')
                fimage[0].header['CRVAL1']-=dRA[good].mean()
                fimage[0].header['CRVAL2']-=dDEC[good].mean()
            if 'IMAGERMS' in fimage[0].header.keys():
                fimage[0].header['IMAGERMS']/=fittedratio
            fimage.flush()
            logger.info('Scaled %s by %.3f' % (image,fittedratio))
        

    if plot:

        plt.clf()
        xx=numpy.logspace(-2,10)
        plt.loglog(xx,xx*fittedratio,'r')
        plt.loglog(xx,numpy.polyval(fitres,xx),
                 'r--')
        plt.errorbar(sourcesTable[good]['GLEAMFlux'],
                     sourcesTable[good]['IntFlux'],
                     xerr=sourcesTable[good]['GLEAMFluxErr'],
                     yerr=sourcesTable[good]['IntFluxErr'],
                     fmt='b.')
        #plt.gca().set_xscale('log')
        #plt.gca().set_yscale('log')
        plt.axis([0.1,2,0.1,2])
        plt.xlabel('Flux Density in %s (Jy)' % catalog,fontsize=16)
        plt.ylabel('Flux Density in %s (Jy)' % image,fontsize=16)
        plt.gca().tick_params(labelsize=16)
        plt.savefig('%s_fluxflux.pdf' % outbase)
        logger.info('Wrote %s_fluxflux.pdf' % outbase)

        plt.clf()
        plt.hist(ratio[good],30)
        plt.xlabel('Flux Density in %s / Flux Density in %s' % (image,catalog),
                   fontsize=16)
        plt.ylabel('Number of Sources',fontsize=16)
        plt.plot(fittedratio*numpy.array([1,1]),
                 plt.gca().get_ylim(),'r-')
        plt.gca().tick_params(labelsize=16)
        plt.savefig('%s_hist.pdf' % outbase)
        logger.info('Wrote %s_hist.pdf' % outbase)
        
        plt.clf()
        plt.plot(x,y,'k.')
        h=plt.scatter(x[good],y[good],s=60,
                      c=ratio[good],
                      norm=matplotlib.colors.LogNorm(vmin=0.5,vmax=2),
                      cmap=plt.cm.BrBG)
        plt.xlabel('X',fontsize=16)
        plt.ylabel('Y',fontsize=16)
        cbar = plt.gcf().colorbar(h,ticks=[0.5,1,2])
        plt.gca().tick_params(labelsize=16)
        plt.savefig('%s_scatter.pdf' % outbase)
        logger.info('Wrote %s_scatter.pdf' % outbase)    

        plt.clf()
        plt.plot((sourcesTable['RA'][good]-sourcesTable['GLEAMRA'][good])*3600,
                 (sourcesTable['Dec'][good]-sourcesTable['GLEAMDEC'][good])*3600,
                 'ro')
        plt.plot(plt.gca().get_xlim(),[0,0],'k--')
        plt.plot([0,0],plt.gca().get_ylim(),'k--')
        plt.xlabel('$\\alpha$(%s)-$\\alpha$(%s)' % (image,catalog),fontsize=16)
        plt.ylabel('$\\delta$(%s)-$\\delta$(%s)' % (image,catalog),fontsize=16)
        plt.gca().tick_params(labelsize=16)
        plt.savefig('%s_position.pdf' % outbase)
        logger.info('Wrote %s_position.pdf' % outbase)    

        plt.clf()
        xx=numpy.linspace(0,300,50)
        plt.hist(sourcesTable['GLEAMSep'].to(u.arcsec).value[~good],
                 xx,color='b',alpha=0.5)
        plt.hist(sourcesTable['GLEAMSep'].to(u.arcsec).value[good],
                 xx,color='r',alpha=0.5)
        plt.plot(matchradius.to(u.arcsec).value*numpy.array([1,1]),
                 plt.gca().get_ylim(),
                 'k--')
        plt.xlabel('Separation %s vs. %s (arcsec)' % (image,catalog),
                   fontsize=16)
        plt.ylabel('Number of sources',fontsize=16)
        plt.gca().tick_params(labelsize=16)
        plt.savefig('%s_separation.pdf' % outbase)
        logger.info('Wrote %s_separation.pdf' % outbase)    

        
    return fittedratio, fittedratioerr, chisq, ndof, fitres[0], fitres[1], sourcesTable.meta['rashift'], sourcesTable.meta['decshift']

                

##################################################                              
def main():
    
    times=collections.OrderedDict()
    times['start']=time.time()

    usage="Usage: %prog [options] <msfiles>\n"
    parser = OptionParser(usage=usage)

    calibration_parser=OptionGroup(parser, 'Calibration options')
    imaging_parser=OptionGroup(parser, 'Imaging options')
    control_parser=OptionGroup(parser, 'Control options')
    
    calibration_parser.add_option('--caltype',dest='caltype',default='casa',
                                  type='choice',
                                  choices=['anoko','casa'],
                                  help='Type of calibration [default=%default]')
    calibration_parser.add_option('--recalibrate',dest='recalibrate',default=False,
                                  action='store_true',
                                  help='Recalibrate existing output?')
    calibration_parser.add_option('--calminuv',dest='calminuv',default=60,
                                  type='float',
                                  help='Minimum UV distance in m for calibration [default=%default]')
    control_parser.add_option('--autoprocess',dest='autoprocess',default=False,
                              action='store_true',
                              help='Run autoprocess?')
    control_parser.add_option('--chgcentre',dest='chgcentre',default=False,
                              action='store_true',
                              help='Run chgcentre?')
    control_parser.add_option('--quick',dest='quick',default=False,
                              action='store_true',
                              help='Run quick initial clean?')    
    control_parser.add_option('--selfcal',dest='selfcal',default=False,
                              action='store_true',
                              help='Run selfcal?')    
    imaging_parser.add_option('--selfcalthreshold',dest='selfcal_clean_threshold',default=0.01,
                              type='float',
                              help='Clean threshold during initial (pre-selfcal) clean (Jy) [default=%default]')
    #control_parser.add_option('--nofluxscale',dest='fluxscale',default=True,
    #                          action='store_false',
    #                          help='Do not scale the fluxes after selfcal using source finding?')    
    #control_parser.add_option('--fluxscaleregion',dest='fluxscaleregion',default=None,type='str',
    #                          help='Region for flux scaling. Default is HPBW.  Can also be radius [deg] around pointing center or MIMAS region file')
    control_parser.add_option('--fluxscale',dest='fluxscale',default=False,
                              action='store_true',
                              help='Do not scale the fluxes with catalog matching')
    control_parser.add_option('--catalog',dest='catalog', default='/work/MWA/GLEAMIDR3.fits',
                              help='GLEAM catalog for flux scaling [default=%default]')
    
    control_parser.add_option('--autosize',default=1,type='float',
                              help='Maximum possible size increase factor to include a bright source [default=%default]')
    imaging_parser.add_option('--size','--imagesize',dest='imagesize',default=2048,type='int',
                              help='Image size in pixels [default=%default]')
    imaging_parser.add_option('--scale',dest='pixelscale',default=0.015,type='float',
                              help='Pixel scale in deg [default=%default]')
    imaging_parser.add_option('--niter',dest='clean_iterations',default=100,
                              type='int',
                              help='Clean iterations [default=%default]')
    imaging_parser.add_option('--threshold',dest='clean_threshold',default=0,
                              type='float',
                              help='Clean threshold (Jy) [default=%default]')
    imaging_parser.add_option('--thresholdsigma',dest='clean_threshold_sigma',default=3,
                              type='float',
                              help='Clean threshold (sigma) [default=%default]')
    imaging_parser.add_option('--mgain',dest='clean_mgain',default=1,
                              type='float',
                              help='Clean major cycle gain [default=%default]')
    imaging_parser.add_option('--gain',dest='clean_gain',default=0.1,
                              type='float',
                              help='Clean gain [default=%default]')
    imaging_parser.add_option('--weight',dest='clean_weight',default='briggs -1.0',
                              help='Clean weighting [default=%default]')
    imaging_parser.add_option('--minuv',dest='clean_minuv',default=0,
                              type='float',
                              help='Minimum UV distance in wavelengths [default=%default]')
    imaging_parser.add_option('--maxuv',dest='clean_maxuv',default=0,
                              type='float',
                              help='Maximum UV distance in wavelengths [default=%default]')
    imaging_parser.add_option('--fullpol','--fullpolarization',dest='fullpolarization',default=False,
                              action='store_true',
                              help='Process full polarization (including cross terms)?')
    imaging_parser.add_option('--nosmallinversion',dest='smallinversion',default=True,
                              action='store_false',
                              help='Do not perform small inversion (sacrifice aliasing for speed)?')
    imaging_parser.add_option('--cleanborder',dest='cleanborder',default=5,
                              type='int',
                              help='Clean border in percent of total image size [default=%default]')
    imaging_parser.add_option('--nofitbeam',dest='fitbeam',default=True,
                              action='store_false',
                              help='Determine beam shape from longest projected baseline (not from fitting the PSF)')
    imaging_parser.add_option('--wscleanargs',dest='wsclean_arguments',default='',
                              help='Additional wsclean arguments')
    imaging_parser.add_option('--subexptime',dest='subexptime',default=None,
                              type='float',
                              help='Exposure time for sub-exposures (s) [default=None]')
    imaging_parser.add_option('--subbands',dest='subbands',default=1,
                              type='int',
                              help='Number of subbands to image [default=%default]')
    parser.add_option('--BANE',dest='BANE',default=False,action='store_true',
                      help='Run BANE for background and RMS estimation?')
    parser.add_option('--beam',dest='beammodel',default='2014i',type='choice',
                      choices=['2014i','2014','2013'],
                      help='Primary beam model [default=%default]')
    parser.add_option('--summaryname',dest='summaryname',default='summary.dat',
                      help='Name for output summary table [default=%default]')
    parser.add_option('--summaryformat',dest='summaryformat',default='ascii.commented_header',
                      help='Format for output summary table (from astropy.table) [default=%default]')   
    parser.add_option('--email',dest='notifyemail', default=None,
                      help='Notification email [default=no email]')
    control_parser.add_option('--nocal',dest='nocal',default=False,
                              action='store_true',
                              help='Do not image the calibrator observations?')
    parser.add_option('-o','--out',dest='out',default='./',
                      help='Output destination directory [default=%default]')
    parser.add_option('--clobber',dest='clobber',default=False,
                      action='store_true',
                      help='Clobber existing output?')
    parser.add_option('--noclobber',dest='clobber',default=False,
                      action='store_false',
                      help='Do not clobber existing output?')
    parser.add_option('--delete',dest='delete',default=False,
                      action='store_true',
                      help='Delete unneeded output?')
    parser.add_option('--nodelete',dest='delete',default=False,
                      action='store_false',
                      help='Do not delete unneeded output?')
    parser.add_option('--cpus',dest='ncpus',default=32,
                      type='int',
                      help='Number of CPUs to be used [default=%default]')
    parser.add_option('--memory',dest='memfraction',default=50,
                      type='int',
                      help='Percentage of total RAM to be used [default=%default]')   
    parser.add_option('--nprocess',dest='nprocess',default=8,
                      type='int',
                      help='Number of possible subprocesses to run [default=%default]')   
    parser.add_option('--casapy',dest='casapy',default=None,
                      help='Path to CASAPY directory (can override with $CASAPY) [default=%default]')
    parser.add_option('--anoko',dest='anoko',default=None,
                      help='Path to anoko/mwa-reduce executable directory (calibrate, applysolutions) (can override with $ANOKO) [default=%default]')
    parser.add_option("-v", "--verbose", dest="loudness", default=0, action="count",
                      help="Each -v option produces more informational/debugging output")
    parser.add_option("-q", "--quiet", dest="quietness", default=0, action="count",
                      help="Each -q option produces less error/warning/informational output")

    parser.add_option_group(control_parser)
    parser.add_option_group(calibration_parser)
    parser.add_option_group(imaging_parser)
    
    command=' '.join(sys.argv)
    (options, args) = parser.parse_args()
    loglevels = {0: [logging.DEBUG, 'DEBUG'],
                 1: [logging.INFO, 'INFO'],
                 2: [logging.WARNING, 'WARNING'],
                 3: [logging.ERROR, 'ERROR'],
                 4: [logging.CRITICAL, 'CRITICAL']}
    logdefault = 2    # WARNING
    level = max(min(logdefault - options.loudness + options.quietness,4),0)
    logging.getLogger('').handlers[1].setLevel(loglevels[level][0])
    #logger.info('Log level set: messages that are %s or higher will be shown.' % loglevels[level][1])

    eh=extra_utils.ExitHandler(os.path.split(sys.argv[0])[-1], email=options.notifyemail)
    logger.info('**************************************************')
    logger.info('%s starting at %s on host %s with user %s' % (sys.argv[0],
                                                               datetime.datetime.now(),
                                                               socket.gethostname(),
                                                               os.environ['USER']))
    logger.debug('In %s' % os.path.abspath(os.curdir))
    logger.debug('Command was:\n\t%s' % command)
    logger.info('**************************************************')


    if options.subbands > 1 and options.subexptime is not None:
        logger.error('Cannot do subbands>1 and subexptime<total together')
        eh.exit(1)

    ##################################################
    # figure out where anoko executables are
    ##################################################
    global anoko
    if options.anoko is not None:
        anoko=ANOKOfinder(options.anoko).find()
    else:
        anoko=ANOKOfinder().find()
    if anoko is None:
        if options.anoko is not None:
            searchpath=ANOKOfinder(options.anoko).paths
        else:
            searchpath=ANOKOfinder().paths
        logger.warning('Unable to find ANOKO in search path:\n\t%s' % ('\n\t'.join(searchpath)))


    files=args
    if len(files)==0:
        logger.error('Must supply >=1 ms files to process')
        eh.exit(1)

    if options.caltype=='casa' and not os.environ.has_key('MWA_CODE_BASE'):
        logger.error('Environment variable $MWA_CODE_BASE is not set; please set and re-run')
        eh.exit(1)

    logger.debug('Using mwapy version %s' % mwapy.__version__)
    result=subprocess.Popen(['wsclean','-version'],
                            stdout=subprocess.PIPE).communicate()[0].strip()
    wscleanversion=result.split('\n')[0].split('version')[1].strip()
    logger.debug('Using wsclean version %s' % wscleanversion)
    #logger.debug('Using casapy version %s' % casapyversion)

    # figure out which if any is a calibrator
    # and which sources it would apply to
    observations=collections.OrderedDict()
    observation_data=numpy.zeros((len(files),),
                                 dtype=[('obsid','i4'),
                                        ('metafits','a40'),
                                        ('inttime','f4'),
                                        ('chanwidth','f4'),
                                        ('ms_cotterversion','a20'),
                                        ('ms_metadataversion','a20'),
                                        ('ms_mwapyversion','a20'),
                                        ('iscalibrator',numpy.bool),
                                        ('channel','i4'),
                                        ('calibrator','i4'),                                        
                                        ('caltype','a10'),
                                        ('calibratorsource','a20'),
                                        ('calibratorfile','a30'),
                                        ('calminuv','f4'),
                                        ('chgcentre',numpy.bool),
                                        ('autoprocess',numpy.bool),
                                        ('selfcal',numpy.bool),
                                        ('clean_threshold','i4'),
                                        ('clean_iterations','i4'),
                                        ('imagesize','i4'),
                                        ('pixelscale','f4'),
                                        ('clean_weight','a20'),
                                        ('clean_minuv','f4'),
                                        ('clean_maxuv','f4'),
                                        ('clean_gain','f4'),
                                        ('clean_mgain','f4'),
                                        ('fullpolarization',numpy.bool),
                                        ('subbands','i4'),
                                        ('subband_bandwidth','f4'),
                                        ('subexposures','i4'),
                                        ('exposure_time','f4'),
                                        ('fluxscale','f4'),
                                        ('fluxscale_nsources','f4'),
                                        ('fluxscale_chisq','f4'),
                                        ('rashift_sec','f4'),
                                        ('decshift_arcsec','f4'),
                                        ('fluxscale_catalog','a200'),
                                        ('wsclean_arguments','a60'),
                                        ('wsclean_command','a200'),
                                        #('rawimages','a30',(6,)),               
                                        ('beammodel','a10'),
                                        #('corrimages','a30',(4,)), 
                                    ])

    ##################################################
    # gather initial data
    ##################################################
    logger.info('**************************************************')
    logger.info('Gathering initial information...')
    if not (os.path.exists(options.out) and os.path.isdir(options.out)):
        logger.warning('Requested output directory %s does not exist; creating...' % options.out)
        try:
            os.mkdir(options.out)
        except Exception,e:
            logger.error('Problem making output directory %s:\n\t%s' % (options.out,e))
            eh.exit(1)            
    increasesize=False  
    for i in xrange(len(files)):
        if not os.path.exists(files[i]):
            logger.error('MS %s does not exist' % files[i])
            eh.exit(1)
        file=files[i]
        obsid=int(file.split('.')[0])
        observations[obsid]=Observation(obsid, 
                                        outputdir=options.out,
                                        caltype=options.caltype,
                                        clobber=options.clobber,
                                        ncpus=options.ncpus,
                                        memfraction=options.memfraction,
                                        delete=options.delete)
        # set a per-observation threshold
        # so that it can be changed later
        observations[obsid].clean_threshold=options.clean_threshold
        if observations[obsid].inttime==0:
            eh.exit(1)        
        observation_data[i]['obsid']=observations[obsid].observation_number
        observation_data[i]['inttime']=observations[obsid].inttime
        observation_data[i]['chanwidth']=observations[obsid].chanwidth
        observation_data[i]['ms_cotterversion']=observations[obsid].otherkeys['MWA_COTTER_VERSION']
        observation_data[i]['ms_metadataversion']=observations[obsid].otherkeys['MWA_METADATA_VERSION']
        observation_data[i]['ms_mwapyversion']=observations[obsid].otherkeys['MWA_MWAPY_VERSION']

        observation_data[i]['iscalibrator']=observations[obsid].calibration
        observation_data[i]['channel']=observations[obsid].center_channel
        if observation_data[i]['iscalibrator']:
            observations[obsid].calibrator=observation_data[i]['obsid']
            observation_data[i]['calibrator']=observation_data[i]['obsid']
            observations[obsid].calibratorsource=''.join(observations[obsid].calibrators)
            observation_data[i]['calibratorsource']=''.join(observations[obsid].calibrators)
            
        else:
            observation_data[i]['calibratorsource']=''

        observation_data[i]['metafits']=observations[obsid].metafits
        pointingcenter=SkyCoord(observations[obsid].RA*u.degree,
                                observations[obsid].Dec*u.degree)

        if not observation_data[i]['iscalibrator']:
            for source in brightsources.keys():
                distance=pointingcenter.separation(brightsources[source])
                if distance > 0.9*(options.imagesize/2)*options.pixelscale*u.degree and distance < options.autosize*(options.imagesize/2)*options.pixelscale*u.degree:
                    logger.warning('Source %s is %.1f deg from field center of %d but outside imaged area; recommend increasing imaged area' % (source,distance.value,observations[obsid].observation_number))
                    increasesize=True

        if False:
            # turn off the self-cal flux scaling
            # do it wrt GLEAM now
            if options.fluxscaleregion is None:
                # default is HPBW
                regionradius=15*(150/(1.28*observations[obsid].center_channel))
            else:
                try:
                    regionradius=float(options.fluxscaleregion)
                except ValueError:
                    if os.path.exists(options.fluxscaleregion):
                        logger.info('Will use MIMAS region %s' % options.fluxscaleregion)
                        observations[obsid].aegean_region=options.fluxscaleregion
                        regionradius=None
                    else:
                        logger.error('Unable to interpret region %s' % options.fluxscaleregion)
                        eh.exit(1)
            if regionradius is not None and observations[obsid].aegean_region is None:
                observations[obsid].aegean_region=circle2mimas(observations[obsid].RA,observations[obsid].Dec,
                                                               regionradius,
                                                               '%d.mim' % obsid)
                logger.info('Will use MIMAS region %s: circle(%.3f,%.3f,%.2f)' % (observations[obsid].aegean_region,
                                                                                  observations[obsid].RA,observations[obsid].Dec,
                                                                                  regionradius))
                                   

    if increasesize and options.autosize>1:
        logger.warning('Increasing imaged area to %dx%d (%.1f deg)' % (options.imagesize*options.autosize,
                                                                       options.imagesize*options.autosize,
                                                                       options.imagesize*options.autosize*options.pixelscale))
        options.imagesize*=options.autosize

    calibrators, notcalibrators, cal_observations=identify_calibrators(observation_data)
    for i in notcalibrators:
        try:
            observation_data[i]['calibrator']=cal_observations[observation_data[i]['obsid']]
            observations[observation_data[i]['obsid']].calibrator=cal_observations[observation_data[i]['obsid']]
        except:
            observation_data[i]['calibrator']=-1
            observations[observation_data[i]['obsid']].calibrator=-1
    times['init']=time.time()

    ##################################################
    # generate calibration solutions
    ##################################################
    logger.info('**************************************************')
    logger.info('Generate calibration solutions...')
    for i in calibrators:
        result=observations[observation_data[i]['obsid']].make_cal(minuv=options.calminuv)
        if result is None:
            eh.exit(1)
        observation_data[i]['calibratorfile']=result
        observation_data[i]['caltype']=observations[observation_data[i]['obsid']].caltype
    times['makecal']=time.time()        

    ##################################################
    # apply calibration
    ##################################################
    logger.info('**************************************************')
    logger.info('Apply calibration solutions...')
    for i in xrange(len(observation_data)):
        if not observation_data[i]['iscalibrator']:
            try:
                cal_touse=numpy.where(observation_data[i]['calibrator']==observation_data['obsid'])[0][0]
            except:
                cal_touse=-1
        else:
            cal_touse=i
        if cal_touse>=0:
            observations[observation_data[i]['obsid']].calibratorfile=observations[observation_data[cal_touse]['obsid']].calibratorfile
        
        result=observations[observation_data[i]['obsid']].calibrate(recalibrate=options.recalibrate)
        if result is False or result is None:
            eh.exit(1)
        observation_data[i]['calibratorfile']=observations[observation_data[i]['obsid']].calibratorfile

    times['applycal']=time.time()                       

    ##################################################
    # autoprocess
    ##################################################
    if options.autoprocess:
        logger.info('**************************************************')
        logger.info('autoprocess...')
        for i in xrange(len(observation_data)):
            if observations[observation_data[i]['obsid']].calibration and options.nocal:
                logger.debug('Skipping autoprocess of calibrator observation %s' % observation_data[i]['obsid'])
                continue
            results=observations[observation_data[i]['obsid']].autoprocess(imagesize=options.imagesize,
                                                                       pixelscale=options.pixelscale)
            observation_data[i]['autoprocess']=True
    times['autoprocess']=time.time()                       

    ##################################################
    # chgcentre
    ##################################################
    if options.chgcentre:
        logger.info('**************************************************')
        logger.info('chgcentre...')
        for i in xrange(len(observation_data)):
            if observations[observation_data[i]['obsid']].calibration and options.nocal:
                logger.debug('Skipping chgcentre of calibrator observation %s' % observation_data[i]['obsid'])
                continue
            results=observations[observation_data[i]['obsid']].chgcentre()
            observation_data[i]['chgcentre']=True
    times['chgcentre']=time.time()                       


    ##################################################
    # quick clean
    ##################################################
    if options.quick and not options.selfcal:
        logger.info('**************************************************')
        logger.info('Quick clean...')    
        for i in xrange(len(observation_data)):
            if observations[observation_data[i]['obsid']].calibration and options.nocal:
                logger.debug('Skipping quick clean of calibrator observation %s' % observation_data[i]['obsid'])
                continue
            # do an initial clean
            results=observations[observation_data[i]['obsid']].image(
                suffix='quick',
                clean_weight=options.clean_weight,
                imagesize=options.imagesize,
                pixelscale=options.pixelscale,
                clean_iterations=observations[observation_data[i]['obsid']].clean_iterations_selfcal,
                clean_gain=options.clean_gain,
                clean_mgain=options.clean_mgain,
                clean_minuv=options.clean_minuv,
                clean_maxuv=options.clean_maxuv,
                clean_threshold=0.01,
                fullpolarization=True,
                smallinversion=options.smallinversion,
                cleanborder=options.cleanborder,
                fitbeam=options.fitbeam,
                wsclean_arguments='-stopnegative',
                updateheader=True)
            if results is None:
                eh.exit(1)
            # we don't want these in the end
            observations[observation_data[i]['obsid']].filestodelete+=results

            # figure out the rms of the residual image
            residimage=results[0].replace('-image','-residual')
            med,rms=stat_measure(residimage)
            logger.info('Measured rms of %d mJy in %s' % (rms*1e3, residimage))
            logger.info('Setting clean threshold to %d*rms=%d mJy' % (options.clean_threshold_sigma,
                                                                      options.clean_threshold_sigma*rms*1e3))
            observations[observation_data[i]['obsid']].clean_threshold=options.clean_threshold_sigma*rms
            observations[observation_data[i]['obsid']].rawfiles=[]
            

    ##################################################
    # selfcal
    ##################################################
    if options.selfcal:
        logger.info('**************************************************')
        logger.info('Selfcal...')    
        for i in xrange(len(observation_data)):
            if observations[observation_data[i]['obsid']].calibration and options.nocal:
                logger.debug('Skipping selfcal of calibrator observation %s' % observation_data[i]['obsid'])
                continue
            # do an initial clean
            results=observations[observation_data[i]['obsid']].image(
                suffix='selfcal',
                clean_weight=options.clean_weight,
                imagesize=options.imagesize,
                pixelscale=options.pixelscale,
                clean_iterations=observations[observation_data[i]['obsid']].clean_iterations_selfcal,
                clean_gain=options.clean_gain,
                clean_mgain=options.clean_mgain,
                clean_minuv=options.clean_minuv,
                clean_maxuv=options.clean_maxuv,
                clean_threshold=options.selfcal_clean_threshold,
                fullpolarization=True,
                smallinversion=options.smallinversion,
                fitbeam=options.fitbeam,
                cleanborder=options.cleanborder,
                wsclean_arguments='-stopnegative',
                updateheader=False)
            if results is None:
                eh.exit(1)
            # we don't want these in the end
            observations[observation_data[i]['obsid']].filestodelete+=results

            # figure out the rms of the residual image
            residimage=results[0].replace('-image','-residual')
            med,rms=stat_measure(residimage)
            logger.info('Measured rms of %d mJy in %s' % (rms*1e3, residimage))
            logger.info('Setting clean threshold to %d*rms=%.1f mJy' % (options.clean_threshold_sigma,
                                                                        options.clean_threshold_sigma*rms*1e3))
            observations[observation_data[i]['obsid']].clean_threshold=options.clean_threshold_sigma*rms

                
            # determine the primary beam
            results=observations[observation_data[i]['obsid']].makepb(beammodel=options.beammodel)
            # we don't want these in the end
            observations[observation_data[i]['obsid']].filestodelete+=results
                                                                      
            if results is None:
                eh.exit(1)

            if True:
                # correct the preliminary image
                results=observations[observation_data[i]['obsid']].pbcorrect()
                if results is None:
                    eh.exit(1)
                Iimage=None
                for image in results:
                    if '-I.fits' in image:
                        Iimage=image
                if False:
                    # turn this off for now
                    if Iimage is None:
                        logger.warning('Could not find I image for source finding')
                    elif _aegean and options.fluxscale:
                        try:
                            observations[observation_data[i]['obsid']].sources,rmsimage,bgimage=find_sources_in_image(Iimage, 
                                                                                                                      max_summits=5,
                                                                                                                      nsigma=10,
                                                                                                                      cores=observations[observation_data[i]['obsid']].ncpus,
                                                                                                                      region=observations[observation_data[i]['obsid']].aegean_region)
                            # if rmsimage is not None:
                            #    observations[observation_data[i]['obsid']].filestodelete.append(rmsimage)
                            # if bgimage is not None:
                            #    observations[observation_data[i]['obsid']].filestodelete.append(bgimage)

                        except Exception,e:
                            logger.error('Unable to run aegean on %s:\n\t%s' % (Iimage,e))
                            eh.exit(1)

                    if len(observations[observation_data[i]['obsid']].sources)==0:
                        logger.warning('Aegean found %d sources in %s' % (len(observations[observation_data[i]['obsid']].sources), Iimage))
                    else:
                        logger.info('Aegean found %d sources in %s' % (len(observations[observation_data[i]['obsid']].sources), Iimage))


            # correct the model image
            results=observations[observation_data[i]['obsid']].pbcorrect(model=True)
            # we don't want these in the end
            observations[observation_data[i]['obsid']].filestodelete+=results

            if results is None:
                eh.exit(1)
            # delete the other Stokes images
            for result in results:
                if '-Q' in result or '-U' in result or '-V' in result:
                    try:
                        os.remove(result)
                    except:
                        logger.warning('Could not delete %s' % result)
                        pass
                else:
                    f=fits.open(result)
                    if f[0].data.max()==0:
                        logger.error('Model image %s has max of 0; you may consider decreasing --selfcalthreshold' % result)
                        eh.exit(1)
            # uncorrect the beam
            results=observations[observation_data[i]['obsid']].pbcorrect(model=True,
                                                                         uncorrect=True)
            if results is None:
                eh.exit(1)
            # we don't want these in the end
            observations[observation_data[i]['obsid']].filestodelete+=results

            results=observations[observation_data[i]['obsid']].image(
                predict=True,
                suffix='selfcal_model',
                clean_weight=options.clean_weight,
                imagesize=options.imagesize,
                pixelscale=options.pixelscale,
                clean_iterations=0,
                clean_gain=options.clean_gain,
                clean_mgain=options.clean_mgain,
                clean_minuv=options.clean_minuv,
                clean_maxuv=options.clean_maxuv,
                clean_threshold=0,
                fullpolarization=True,
                smallinversion=options.smallinversion,
                fitbeam=options.fitbeam,
                cleanborder=options.cleanborder,
                updateheader=False)

            results=observations[observation_data[i]['obsid']].make_cal(minuv=options.calminuv,
                                                                        selfcal=True)
            if results is None:
                eh.exit(1)
            
            result=observations[observation_data[i]['obsid']].calibrate(selfcal=True)
            if result is False or result is None:
                eh.exit(1)
            observation_data[i]['calibratorfile']=observations[observation_data[i]['obsid']].calibratorfile
            observations[observation_data[i]['obsid']].selfcal=True

            observations[observation_data[i]['obsid']].rawfiles=[]
            observations[observation_data[i]['obsid']].beamfiles=[]
            observations[observation_data[i]['obsid']].corrfiles=[]
            observation_data[i]['selfcal']=True

    times['selfcal']=time.time()                       


    ##################################################
    # imaging
    ##################################################
    logger.info('**************************************************')
    logger.info('Image...')    
    for i in xrange(len(observation_data)):
        if observations[observation_data[i]['obsid']].calibration and options.nocal:
            logger.debug('Skipping imaging of calibrator observation %s' % observation_data[i]['obsid'])
            continue

        if options.subexptime is not None and options.subexptime > 0:
            # make multiple images for different subexposures
            results=observations[observation_data[i]['obsid']].multiimage_time(
                imagetime=options.subexptime,
                clean_weight=options.clean_weight,
                imagesize=options.imagesize,
                pixelscale=options.pixelscale,
                clean_iterations=options.clean_iterations,
                clean_gain=options.clean_gain,
                clean_mgain=options.clean_mgain,
                clean_minuv=options.clean_minuv,
                clean_maxuv=options.clean_maxuv,
                clean_threshold=observations[observation_data[i]['obsid']].clean_threshold,
                fullpolarization=options.fullpolarization,
                wsclean_arguments=options.wsclean_arguments,
                ntorun=options.nprocess)
            if results is None:
                eh.exit(1)

            if options.fullpolarization:
                observation_data[i]['subexposures']=len(results)/4
            else:
                observation_data[i]['subexposures']=len(results)/2

        else:
            results=observations[observation_data[i]['obsid']].image(
                clean_weight=options.clean_weight,
                imagesize=options.imagesize,
                pixelscale=options.pixelscale,
                clean_iterations=options.clean_iterations,
                clean_gain=options.clean_gain,
                clean_mgain=options.clean_mgain,
                clean_minuv=options.clean_minuv,
                clean_maxuv=options.clean_maxuv,
                clean_threshold=observations[observation_data[i]['obsid']].clean_threshold,
                fullpolarization=options.fullpolarization,
                smallinversion=options.smallinversion,
                fitbeam=options.fitbeam,
                cleanborder=options.cleanborder,
                subbands=options.subbands,                
                wsclean_arguments=options.wsclean_arguments)
            if results is None:
                eh.exit(1)
            
            observation_data[i]['subbands']=options.subbands

        for j in xrange(len(results)):
            file=results[j]
            try:
                f=fits.open(file,'update')
            except Exception,e:
                logger.open('Unable to open image output %s for updating:\n\t%s' % (file,e))
                eh.exit(1)

            if j==0 and observation_data[i]['subbands']>1:
                observation_data[i]['subband_bandwidth']=f[0].header['BANDWDTH']
            if j==0:
                observation_data[i]['exposure_time']=f[0].header['EXPOSURE']

            #f[0].header['CASAVER']=(casapyversion,'CASAPY version')
            f[0].header['MWAPYVER']=(mwapy.__version__,'MWAPY version')
            f[0].header['WSCLNVER']=(wscleanversion,'WSCLEAN version')
            f[0].header['COTTRVER']=(observation_data[i]['ms_cotterversion'],'Cotter version')
            f[0].header.add_history(command)
            f.flush()
            #observation_data[i]['rawimages'][j]=results[j]

        observation_data[i]['clean_iterations']=options.clean_iterations
        observation_data[i]['clean_threshold']=observations[observation_data[i]['obsid']].clean_threshold
        observation_data[i]['clean_weight']=options.clean_weight
        observation_data[i]['clean_gain']=options.clean_gain
        observation_data[i]['clean_mgain']=options.clean_mgain
        observation_data[i]['imagesize']=options.imagesize
        observation_data[i]['pixelscale']=options.pixelscale
        observation_data[i]['clean_minuv']=options.clean_minuv
        observation_data[i]['clean_maxuv']=options.clean_maxuv
        observation_data[i]['fullpolarization']=options.fullpolarization
        observation_data[i]['wsclean_arguments']=options.wsclean_arguments
        try:
            observation_data[i]['wsclean_command']=' '.join(observations[observation_data[i]['obsid']].wscleancommand)
        except:
            pass

        if options.subbands==1:
            results=observations[observation_data[i]['obsid']].makepb(beammodel=options.beammodel)
            if results is None:
                eh.exit(1)
        else:
            for subband in xrange(options.subbands+1):
                filesperband=2+2*options.fullpolarization
                results=observations[observation_data[i]['obsid']].makepb(beammodel=options.beammodel,
                                                                          ifile=subband*filesperband)
                if results is None:
                    eh.exit(1)
        observation_data[i]['beammodel']=options.beammodel

        if observations[observation_data[i]['obsid']].nimages==1:
            if options.subbands==1:
                results=observations[observation_data[i]['obsid']].pbcorrect()
                if results is None:
                    eh.exit(1)
            else:
                for subband in xrange(options.subbands+1):
                    filesperband=2+2*options.fullpolarization
                    results=observations[observation_data[i]['obsid']].pbcorrect(ifile=subband*filesperband)
                    if results is None:
                        eh.exit(1)
 
        else:
            results=observations[observation_data[i]['obsid']].multipbcorrect_time()
            if results is None:
                eh.exit(1)


        if options.fluxscale:
            if not os.path.exists(options.catalog):
                logger.warning('Could not find GLEAM catalog %s' % options.catalog)
            else:
                observation_data[i]['fluxscale_catalog']=options.catalog
                # need to find the I image
                # if subbands, need to make sure it is MFS
                Iimage=None
                for image in results:
                    if options.subbands==1:
                        if '-I.fits' in image:
                            Iimage=image
                    else:
                        if '-MFS-I.fits' in image:
                            Iimage=image
                    out=fluxmatch(Iimage, catalog=options.catalog,
                                  nsigma=10,
                                  rmsfactor=4,
                                  matchradius=60,
                                  rejectsigma=3,
                                  maxdistance=20,
                                  minbeam=0.4,
                                  psfextent=1.1,
                                  refineposition=True,
                                  plot=True,
                                  update=True,
                                  cores=1)
                    if out is None:
                        logger.warning('Flux match with %s returned None for %s; skipping...' % (options.catalog,
                                                                                                 Iimage))
                    else:
                        fittedratio, fittedratioerr, chisq, ndof, slope, intercept, rashift, decshift=out
                        observation_data[i]['fluxscale']=fittedratio
                        observation_data[i]['fluxscale_nsources']=ndof+1
                        observation_data[i]['fluxscale_chisq']=chisq
                        observation_data[i]['rashift_sec']=rashift
                        observation_data[i]['decshift_arcsec']=decshift
                    

        if options.BANE:
            for image in results:
                rmsimage,bgimage=get_rms_background(image)
                logger.info('Created rms image %s, background image %s' % (rmsimage,bgimage))

        # do another round of source finding if we did selfcal
        if options.selfcal:
            fluxratio=None
            # need to find the I image
            # if subbands, need to make sure it is MFS
            Iimage=None
            for image in results:
                if options.subbands==1:
                    if '-I.fits' in image:
                        Iimage=image
                else:
                    if '-MFS-I.fits' in image:
                        Iimage=image
            if False:
                if Iimage is None:
                    logger.warning('Could not find I image for source finding')
                elif _aegean and options.fluxscale:
                    try:
                        newsources,newrmsimage,newbgimage=find_sources_in_image(Iimage, 
                                                                                max_summits=5,
                                                                                nsigma=10,
                                                                                cores=observations[observation_data[i]['obsid']].ncpus,
                                                                                region=observations[observation_data[i]['obsid']].aegean_region)
                        if newrmsimage is not None:
                            observations[observation_data[i]['obsid']].filestodelete.append(newrmsimage)
                        if newbgimage is not None:
                            observations[observation_data[i]['obsid']].filestodelete.append(newbgimage)

                    except Exception,e:
                        logger.warning('Unable to run aegean on %s:\n\t%s' % (Iimage,e))
                        # eh.exit(1)
                        newsources=[]
                        newrmsimage,newbgimage=None,None
                    if len(newsources)==0:
                        logger.warning('Aegean found %d sources in %s' % (len(newsources), Iimage))
                    else:
                        logger.info('Aegean found %d sources in %s' % (len(newsources), Iimage))
                        origsources, newsources=match_aegean_sources(observations[observation_data[i]['obsid']].sources,
                                                                     newsources)
                        if len(origsources)==0:
                            logger.warning('Could not match sources between images')
                        else:
                            origfluxes=numpy.array([s.peak_flux for s in origsources])
                            newfluxes=numpy.array([s.peak_flux for s in newsources])
                            fluxratio=numpy.median(newfluxes/origfluxes)
                            logger.info('Determined median flux ratio of %.2f from %d matched sources' % (fluxratio,
                                                                                                          len(origfluxes)))
                            if fluxratio < 0.5 or fluxratio > 2:
                                logger.warning('Flux ratio exceeds acceptable limits; skipping scaling...')
                            else:
                                for image in results:
                                    try:
                                        f=fits.open(image,'update')
                                        f[0].data/=fluxratio
                                        f[0].header['FLUXSCAL']=(fluxratio,'Ratio of flux compared to pre-selfcal image')
                                        f.flush()
                                        logger.info('Divided %s by %.2f to correct flux scale' % (image,fluxratio))
                                    except:
                                        logger.warning('Unable to correct flux scale for %s' % image)

        # for j in xrange(len(results)):
        #    observation_data[i]['corrimages'][j]=results[j]

   
    times['image']=time.time()
    try:
        observation_data_table=Table(observation_data)
    except Exception,e:
        logger.error('Unable to create Table for summary data:\n\t%s' % e)
        eh.exit(1)

    try:
        observation_data_table.write(options.summaryname,
                                     delimiter='|',
                                     format=options.summaryformat)
        logger.info('Summary table written to %s' % options.summaryname)
    except Exception, e:
        print observation_data[0]
        logger.error('Unable to write summary table %s with format %s:\n%s' % (options.summaryname,
                                                                               options.summaryformat,
                                                                               e))
            
    times['end']=time.time()


    # get rid of extraneous log files
    observations[observation_data[i]['obsid']].filestodelete+=glob.glob('casapy-*.log')
    observations[observation_data[i]['obsid']].filestodelete+=glob.glob('ipython-*.log')

    logger.info('Execution time: %d s (setup), %d s (make cal), %d s (apply cal), %d s (autoprocess), %d s (chgcentre), %d s (selfcal), %d s (image) %d s (tidy) = %d s (total)' % (
        times['init']-times['start'],
        times['makecal']-times['init'],
        times['applycal']-times['makecal'],
        times['autoprocess']-times['applycal'],        
        times['chgcentre']-times['autoprocess'],        
        times['selfcal']-times['chgcentre'],        
        times['image']-times['selfcal'],        
        times['end']-times['image'],
        times['end']-times['start']))




    eh.exit(0)

                                                                
######################################################################

if __name__=="__main__":

    main()
