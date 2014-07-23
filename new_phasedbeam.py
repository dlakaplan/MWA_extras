import numpy,tpipe,os,sys
import astropy.io.fits as pyfits
import astropy.wcs as pywcs
from mwapy import ephem_utils

minuv=None
sourcename=None
ra,dec=None,None

#filename='1063400456_corr.ms'
filename=sys.argv[5]
#filename='1070377184.ms'
#ra,dec=69.3158945833333, -47.2523944444444
#ra,dec=63.2009, -56.0136
#sourcename='J0412b'
#sourcename='PSRJ0437'
if len(sys.argv)>6:
    minuv=int(sys.argv[6])


#filename='1070380064_autopeel_PKS0436.ms'
#ra,dec=None,None

nchannels=None
maxbuf=16
nints=None
if minuv is not None:
    outfilename=filename.replace('.ms',
                                 '_cut%d_phased.fits' % minuv)
else:
    outfilename=filename.replace('.ms',
                                 '_phased.fits')
if sourcename is not None:
    outfilename=outfilename.replace('_phased.fits','_%s_phased.fits' % sourcename)

print 'Reading %s...' % filename
print 'Will write %s...' % outfilename

if nchannels is None or nints is None:
    ms.open(filename)
    spwinfo = ms.getspectralwindowinfo()
    scansummary = ms.getscansummary()
    if nchannels is None:
        nchannels=spwinfo['0']['NumChan']
        print 'Determined %d channels' % nchannels
    if nints is None:
        BeginTime=scansummary['1']['0']['BeginTime']
        EndTime=scansummary['1']['0']['EndTime']
        IntegrationTime=scansummary['1']['0']['IntegrationTime']
        nints=int((EndTime-BeginTime)*86400/IntegrationTime)
        print 'Determined %d integrations' % nints

XX=numpy.ma.zeros((nints,nchannels),dtype=numpy.complex128)
YY=numpy.ma.zeros((nints,nchannels),dtype=numpy.complex128)
XY=numpy.ma.zeros((nints,nchannels),dtype=numpy.complex128)
YX=numpy.ma.zeros((nints,nchannels),dtype=numpy.complex128)
MJD=numpy.zeros((nints))


nint_remaining=nints
nstart=0
nskip=0
data=None
dataC=None
uvw=None
badbaselines=None
while nint_remaining>0:
    ntoread=min(maxbuf, nint_remaining)
    print 'Reading %d integrations starting with %d' % (ntoread, nskip)
    pipe=tpipe.pipe_msint(filename,nints=ntoread,nskip=nskip,
                          selectpol=['XX','XY','YX','YY'],
                          dmarr=[0.],chans=numpy.arange(nchannels),
                          datacol='data')            
    if uvw is None and minuv is not None:
        uvw=numpy.sqrt(pipe.u[0]**2+pipe.v[0]**2+pipe.w[0]**2)                 
        # flag short baselines
        badbaselines=uvw<=minuv
        print 'Flagging %d%% of baselines that are <%d lambda' % (badbaselines.sum()/(1.0*len(badbaselines))*100,
                                                                  minuv)
        
    if data is None:
        data=numpy.ma.zeros(pipe.data.shape)
        dataC=numpy.ma.zeros(pipe.data.shape)

    if ra is not None:
        # rotate the phase center
        dl,dm=pipe.get_shift(ra, dec)
        pipe.phaseshift(dl,dm)      

    if badbaselines is not None:
        pipe.data.mask[:,badbaselines,:,:]=True

    data=pipe.data
    dataC=pipe.data.conj()

    # add the conjugate of the data
    # V(u,v,w) = V*(-u,-v,-w)
    phased_data=data.mean(axis=1) + dataC.mean(axis=1)
    #phased_data=data.mean(axis=1)
    XX[nstart:nstart+ntoread]=phased_data[:,:,0]
    XY[nstart:nstart+ntoread]=phased_data[:,:,1]
    YX[nstart:nstart+ntoread]=phased_data[:,:,2]
    YY[nstart:nstart+ntoread]=phased_data[:,:,3]
    MJD[nstart:nstart+ntoread]=pipe.time-2400000.5
    nint_remaining-=ntoread
    nskip+=ntoread
    nstart+=ntoread
            
tint=numpy.round((MJD[1]-MJD[0])*86400,1)

wcs=pywcs.WCS(naxis=2)
wcs.wcs.ctype=['FREQ-LSR','UTC']
wcs.wcs.cunit=['MHz','s']
wcs.wcs.crval=[numpy.round(pipe.freq[0]*1e3,2),0.0]
wcs.wcs.crpix=[1.0,0.5]
wcs.wcs.cdelt=[numpy.round((pipe.freq[1]-pipe.freq[0])*1e3,2),tint]
f=pyfits.PrimaryHDU(header=wcs.to_header())


f.header['TIMESYS']='UTC'
t=ephem_utils.MWATime()
t.MJD=MJD[0]
f.header['DATEREF']=t.datetime.strftime('%Y-%m-%dT%H:%M:%S')
f.header['FILENAME']=(filename,'Input file')
f.header['NSKIP']=(0,'Number of integrations that were skipped')
f.header['NINT']=(nints,'Number of integrations that were processed')
f.header['INTTIME']=(tint,'[s] Integration time')
f.header['CHANNEL']=(int(((pipe.freq[1]-pipe.freq[0])*1e6)),
                     '[kHz] Channel width')
if ra is not None:
    f.header['RAPHASE']=(ra,'[deg] RA of new phase center')
    f.header['DECPHASE']=(dec,'[deg] DEC of new phase center')


XX.data[XX.mask]=0
YY.data[YY.mask]=0
XY.data[XY.mask]=0
YX.data[YX.mask]=0

l=[f]
l.append(pyfits.ImageHDU(data=XX.data.real,
                         header=f.header))
l[-1].header['EXTNAME']='XX-REAL'
l.append(pyfits.ImageHDU(data=XX.data.imag,
                         header=f.header))
l[-1].header['EXTNAME']='XX-IMAG'
l.append(pyfits.ImageHDU(data=YY.data.real,
                         header=f.header))
l[-1].header['EXTNAME']='YY-REAL'
l.append(pyfits.ImageHDU(data=YY.data.imag,
                         header=f.header))
l[-1].header['EXTNAME']='YY-IMAG'
l.append(pyfits.ImageHDU(data=XY.data.real,
                         header=f.header))
l[-1].header['EXTNAME']='XY-REAL'
l.append(pyfits.ImageHDU(data=XY.data.imag,
                         header=f.header))
l[-1].header['EXTNAME']='XY-IMAG'
l.append(pyfits.ImageHDU(data=YX.data.real,
                         header=f.header))
l[-1].header['EXTNAME']='YX-REAL'
l.append(pyfits.ImageHDU(data=YX.data.imag,
                         header=f.header))
l[-1].header['EXTNAME']='YX-IMAG'

if os.path.exists(outfilename):
    os.remove(outfilename)
hl=pyfits.HDUList(l)
hl.writeto(outfilename)
