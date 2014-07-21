from mwapy import ephem_utils
import sys
import numpy,math,numpy.ma
import os,time,datetime
import ephem
import tpipe,mwa_pipe
import pyfits,pywcs
import logging
from optparse import OptionParser
#import pylab

# /usr/local/lib/python2.6/dist-packages/:/localhome/kaplan/python/
# PYTHONPATH=/usr/local/lib/python2.6/dist-packages/:/localhome/kaplan/python/:/localhome/kaplan/lib/python2.6/site-packages/; export PYTHONPATH

logging.basicConfig(format='# %(levelname)s:%(name)s: %(message)s')
logger=logging.getLogger('make_phasedbeam')
logger.setLevel(logging.WARNING)

_maxnint=32

################################################################################
def getpath(satellite, MJD, dt=0, ds=0.0):
    """
    ra,dec,alt,speed,angle=getpath(satellite, MJD, dt=0, ds=0)
    all are returned as degrees
    dt (sec) adjusts the position forward/backward along the path
    ds (deg) adjusts the position perpendicular to the track
    """
    observer=ephem.Observer()
    # make sure no refraction is included
    observer.pressure=0
    mwa=ephem_utils.Obs[ephem_utils.obscode['MWA']]
    observer.long=mwa.long/ephem_utils.DEG_IN_RADIAN
    observer.lat=mwa.lat/ephem_utils.DEG_IN_RADIAN
    observer.elevation=mwa.elev
    ra=[]
    dec=[]
    alt=[]
    speed=[]
    angle=[]
    t=ephem_utils.MWATime()
    for mjd in MJD: 
        t.MJD=mjd
        tuse=t.datetime+datetime.timedelta(seconds=dt)
        #print tuse
        observer.date=tuse
        #print float(observer.date),observer.date,observer.date.tuple()
        # make it deal with fractional seconds
        observer.date=float(observer.date)+tuse.microsecond/1.0e6*ephem.second
        date0=float(observer.date)
        satellite.compute(observer)
        ra.append(satellite.ra)
        dec.append(satellite.dec)
        alt.append(satellite.alt)
        #print '%.7f %s %f %f' % (observer.date,observer.date,numpy.degrees(satellite.ra),numpy.degrees(satellite.dec))
        #print ''
        observer.date=date0-0.5*ephem.second
        satellite.compute(observer)
        ra1=satellite.ra
        dec1=satellite.dec
        observer.date=date0+0.5*ephem.second
        satellite.compute(observer)
        ra2=satellite.ra
        dec2=satellite.dec
        ddec=dec2-dec1
        dra=ra2-ra1
        dracosdec=dra*numpy.cos(dec[-1])
        speed.append(numpy.sqrt(dracosdec**2+ddec**2))
        angle.append(numpy.arctan2(ddec,dracosdec))
        
    ra=numpy.degrees(numpy.array(ra))
    dec=numpy.degrees(numpy.array(dec))
    alt=numpy.degrees(numpy.array(alt))
    speed=numpy.degrees(numpy.array(speed))
    angle=numpy.degrees(numpy.array(angle))
    dx=ds*numpy.sin(numpy.radians(angle))
    dy=-ds*numpy.cos(numpy.radians(angle))
    ra+=dx/numpy.cos(numpy.radians(dec))
    dec+=dy

    return ra,dec,alt,speed,angle


################################################################################
def main():

    try:
        x=tpipe.miriad
    except AttributeError:
        logging.error('tpipe.miriad is not available')
        sys.exit(0)

    usage="Usage: %prog [options]\n"
    parser = OptionParser(usage=usage)
    parser.add_option('-f','--filename',dest="filename",default=None,
                      help="Miriad filename")
    parser.add_option('-o','--output',dest="output",default=None,
                      help="Output filename")
    parser.add_option('--channels',dest="channels",default=768,type='int',
                      help="Number of channels")
    parser.add_option('-i','--integrations',dest="integrations",default=0,type='int',
                      help="Number of integrations to process")
    parser.add_option('-s','--skip',dest="skip",default=0,type='int',
                      help="Number of integrations to skip")
    parser.add_option('--tle','--satellite',dest="tlefile",default=None,
                      help="Satellite TLE file",metavar="TLEFILE")
    parser.add_option('--dt',dest="dt",default=0,type='float',
                      help="Number of integrations to shift in absolute WCS")
    parser.add_option('--ds',dest="ds",default=0,type='float',
                      help="Distance in degrees to shift in absolute WCS")
    parser.add_option('--dra',dest="dra",default=0,type='float',
                      help="RA distance in degrees to shift")
    parser.add_option('--ddec',dest="ddec",default=0,type='float',
                      help="Dec distance in degrees to shift")
    parser.add_option('--filter',action="store_true",dest="filter",default=False,
                      help="Filter the phased-array beam according to the expected shape")
    parser.add_option('--nofilter',action="store_false",dest="filter",default=False,
                      help="Do not filter the phased-array beam according to the expected shape")
    parser.add_option('--width',default=3.0,type='float',dest='width',
                      help="Width for filter (arcmin)")
    parser.add_option('--ra',default=None,type='float',dest='ra',
                      help="Right Ascension (deg) for phasing")
    parser.add_option('--dec',default=None,type='float',dest='dec',
                      help="Declination (deg) for phasing")
    parser.add_option('--maxbuf',default=_maxnint,type='int',dest='maxbuf',
                      help="Maximum number of integrations to process at once")
    parser.add_option('--flip',action="store_true",dest="flip",default=False)
    parser.add_option('--flag',dest="flag",default='',type='str',
                      help='List of antennas to flag (1 origin)')
    parser.add_option('--imag',action="store_true",dest="imag",default=False)


    (options, args) = parser.parse_args()

    nchan=options.channels
    nskip=options.skip
    nints=options.integrations
    outputfile=options.output

    if options.filename is None:
        logging.error('Must supply Miriad file')
        sys.exit(0)

    if not os.path.exists(options.filename):
        logging.error('Miriad file %s does not exist' % options.filename)
        sys.exit(1)
        
    satellite=None
    if options.tlefile is not None:
        if not os.path.exists(options.tlefile):
            logging.error('TLE file %s does not exist' % options.tlefile)
            sys.exit(1)
        f=open(options.tlefile)
        tlelines=f.readlines()
        satellite_label=tlelines[0].replace('_','\_').replace('\n','')
        satellite=ephem.readtle(tlelines[0],
                                tlelines[1],
                                tlelines[2])

    toflag=[]
    if len(options.flag)>0:
        toflag=[int(x) for x in options.flag.split(',')]
        print 'Will flag antennas %s' % toflag
    
    dataph=numpy.zeros((nints,nchan,4))
    if options.imag:
        dataph_imag=numpy.zeros((nints,nchan,4))
    MJD=numpy.zeros((nints))
    RA=numpy.zeros_like(MJD)
    Dec=numpy.zeros_like(MJD)
    Altitude=numpy.zeros_like(MJD)
    Angle=numpy.zeros_like(MJD)
    Speed=numpy.zeros_like(MJD)
    DL=numpy.zeros_like(MJD)
    DM=numpy.zeros_like(MJD)
    if nints <= options.maxbuf:
        try:
            pipe=mwa_pipe.pipe_mwa(options.filename,profile='mwa',
                                   chans=numpy.arange(nchan),
                                   nints=nints,dmarr=[0.],nskip=nskip,
                                   selectpol=['XX','YY','XY','YX'])        
        except:
            logging.error('Error reading data')
            sys.exit(1)
        print 'Data read'
        if len(toflag)>0:
            pipe.flag_antennas(numpy.array(toflag)-1)

        #dataph=pipe.dataph.data
        dataph=(pipe.data.mean(axis=1)).real
        if options.imag:
            dataph_imag=(pipe.data.mean(axis=1)).imag        
        MJD=pipe.time-2400000.5

        if options.ra is not None:
            dl,dm=pipe.get_shift(options.ra,options.dec)
            pipe.phaseshift(dl,dm)
            dataph=(pipe.data.mean(axis=1)).real
            #dataph=pipe.dataph
            if options.imag:
                dataph_imag=(pipe.data.mean(axis=1)).imag        


    else:
        nint_remaining=nints
        nstart=0
        #if isinstance(dataph_out,numpy.ma.core.MaskedArray):
        #    dataph_imag=numpy.ma.zeros((nints,nchan))            

        while nint_remaining>0:
            ntoread=min(options.maxbuf,nint_remaining)
            print 'Reading %d integrations starting with %d' % (ntoread,
                                                                nskip)
            pipe=mwa_pipe.pipe_mwa(options.filename,profile='mwa',
                                   chans=numpy.arange(nchan),
                                   nints=ntoread,
                                   dmarr=[0.],nskip=nskip,
                                   selectpol=['XX','YY','XY','YX']) 
            if len(toflag)>0:
                pipe.flag_antennas(numpy.array(toflag)-1)

            #dataph_out=pipe.dataph
            dataph_out=(pipe.data.mean(axis=1)).real
            if options.imag:
                dataph_out_imag=(pipe.data.mean(axis=1)).imag                            

            if options.ra is not None:
                dl,dm=pipe.get_shift(options.ra,options.dec)
                pipe.phaseshift(dl,dm)
                #dataph_out=pipe.dataph
                dataph_out=(pipe.data.mean(axis=1)).real
                if options.imag:
                    dataph_out_imag=(pipe.data.mean(axis=1)).imag        

            dataph[nstart:nstart+min(options.maxbuf,nint_remaining)]=dataph_out
            if options.imag:
                dataph_imag[nstart:nstart+min(options.maxbuf,nint_remaining)]=dataph_out_imag
                
            MJD[nstart:nstart+min(options.maxbuf,nint_remaining)]=pipe.time-2400000.5
            try:
                RA[nstart:nstart+min(options.maxbuf,nint_remaining)]=ra
                Dec[nstart:nstart+min(options.maxbuf,nint_remaining)]=dec
                Altitude[nstart:nstart+min(options.maxbuf,nint_remaining)]=alt
                Speed[nstart:nstart+min(options.maxbuf,nint_remaining)]=speed
                Angle[nstart:nstart+min(options.maxbuf,nint_remaining)]=angle
                DL[nstart:nstart+min(options.maxbuf,nint_remaining)]=dl
                DM[nstart:nstart+min(options.maxbuf,nint_remaining)]=dm
            except:
                pass

            nint_remaining-=ntoread
            nskip+=ntoread
            nstart+=ntoread
            
    tint=numpy.round((MJD[1]-MJD[0])*86400,1)

    wcs=pywcs.WCS(naxis=3)
    # see Greisen & Calabretta 2002, 395, 1061
    # Table 7
    # I=1
    # XX=-5
    # YY=-6
    # XY=-7
    # YX=-8

    wcs.wcs.ctype=['FREQ-LSR','UTC','STOKES']
    wcs.wcs.cunit=['MHz','s','']
    wcs.wcs.crval=[numpy.round(pipe.freq[0]*1e3,2),0.0,-5]
    wcs.wcs.crpix=[1.0,0.5,1]
    wcs.wcs.cdelt=[numpy.round((pipe.freq[1]-pipe.freq[0])*1e3,2),tint,1]
    f=pyfits.PrimaryHDU(header=wcs.to_header())
    f.header['TIMESYS']='UTC'
    t=ephem_utils.MWATime()
    t.MJD=MJD[0]
    f.header['DATEREF']=t.datetime.strftime('%Y-%m-%dT%H:%M:%S')
    if isinstance(dataph,numpy.ma.core.MaskedArray):
        dataph.data[dataph.mask]=0
        f.data=dataph.data
    else:
        f.data=dataph
    if os.path.exists(outputfile):
        os.remove(outputfile)
    f.header['FILENAME']=(options.filename,'Input miriad file')
    f.header['NSKIP']=(options.skip,'Number of integrations that were skipped')
    f.header['NINT']=(nints,'Number of integrations that were processed')
    f.header['INTTIME']=(tint,'[s] Integration time')
    f.header['CHANNEL']=(int(((pipe.freq[1]-pipe.freq[0])*1e6)),
                         '[kHz] Channel width')
    if len(toflag)>0:
        f.header['FLAGGED']=(','.join([str(x) for x in toflag]),'Flagged antennas')
    else:
        f.header['FLAGGED']=('NONE','Flagged antennas')
    if satellite is not None:
        f.header['DT']=(options.dt,'[s] Time offset for absolute WCS')
        f.header['DS']=(options.ds,'[deg] Perpendicular offset for absolute WCS')
        f.header['DRA']=(options.dra,'[deg] RA offset')
        f.header['DDEC']=(options.ddec,'[deg] Dec offset')
        f.header['TLEFILE']=(options.tlefile,'Satellite TLE filename')
        if options.filter:
            f.header['FILTER']=(pyfits.TRUE,'Filtered phased-array data?')
        else:
            f.header['FILTER']=(pyfits.FALSE,'Filtered phased-array data?')        
        c1=pyfits.Column(name='MJD',format='D',unit='day',array=MJD)
        c2=pyfits.Column(name='RA',format='E',unit='deg',array=RA)
        c3=pyfits.Column(name='Dec',format='E',unit='deg',array=Dec)
        c4=pyfits.Column(name='Altitude',format='E',unit='deg',array=Altitude)
        c5=pyfits.Column(name='Speed',format='E',unit='deg/int',array=Speed)
        c6=pyfits.Column(name='Angle',format='E',unit='deg',array=Angle)
        c7=pyfits.Column(name='DL',format='E',unit='deg',array=DL)
        c8=pyfits.Column(name='DM',format='E',unit='deg',array=DM)
        coldefs=pyfits.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8])
        tbhdu=pyfits.new_table(coldefs)
        ftemp=f
        f=pyfits.HDUList([ftemp,tbhdu])
    if options.ra is not None:
        f.header['RAPHASE']=(options.ra,'[deg] RA of phase center')
        f.header['DECPHASE']=(options.dec,'[deg] Dec of phase center')

    if options.imag:
        z=numpy.zeros(dataph_imag.shape)
        if isinstance(dataph_imag,numpy.ma.core.MaskedArray):
            dataph_imag.data[dataph_imag.mask]=0
            z=dataph_imag.data
        else:
            z=dataph_imag
        if isinstance(f,pyfits.hdu.hdulist.HDUList):            
            f.append(pyfits.ImageHDU(data=z))
        else:
            ftemp=f
            f=pyfits.HDUList([ftemp,pyfits.ImageHDU(data=dataph_imag.data)])
            #f.data=dataph_imag.data
        
    
    f.writeto(outputfile)
    

    sys.exit(0)
    
################################################################################
    
if __name__=="__main__":
    main()
