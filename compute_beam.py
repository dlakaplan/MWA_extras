import math,numpy,sys,numpy.linalg,os
import ephem
from mwapy import ephem_utils,get_observation_info
from mwapy.pb import primary_beam
from mwapy.obssched.base import schedule
import astropy.io.fits as pyfits

# open up database connection
try:
	db = schedule.getdb()
except:
	print "Unable to open connection to database"
	sys.exit(1)

RA,Dec=69.3158945833333, -47.2523944444444
gpstime=1063400456
filename='1063400456_corr_phased.fits'

o=get_observation_info.MWA_Observation(gpstime, db=db)

mwatime=ephem_utils.MWATime(gpstime=gpstime+o.duration/2)
frequency=o.center_channel*1.28e6
delays=o.delays

RAnow,Decnow=ephem_utils.precess(RA,Dec,2000,mwatime.epoch)
HA=float(mwatime.LST)-RAnow
mwa=ephem_utils.Obs[ephem_utils.obscode['MWA']]
Az,Alt=ephem_utils.eq2horz(HA,Decnow,mwa.lat)
theta=(90-Alt)*math.pi/180
phi=Az*math.pi/180


def MWA_Jones_Analytic(theta, phi, ha, dec,
                       freq=100.0e6, delays=None,
                       zenithnorm=True):
    c=2.998e8
    # wavelength in meters
    lam=c/freq

    dip_sep=primary_beam._DIPOLE_SEPARATION
    delay_int=primary_beam._DELAY_INT
    dipheight=primary_beam._DIPOLE_HEIGHT

    if (delays is None):
        delays=0

    if (isinstance(delays,float) or isinstance(delays,int)):
        delays=delays*numpy.ones((16))
    if (isinstance(delays,numpy.ndarray) and len(delays)==1):
        delays=delays[0]*numpy.ones((16))        

    # direction cosines (relative to zenith) for direction az,za
    projection_east=numpy.sin(theta)*numpy.sin(phi)
    projection_north=numpy.sin(theta)*numpy.cos(phi)
    projection_z=numpy.cos(theta)

    # dipole position within the tile
    dipole_north=dip_sep*numpy.array([1.5,1.5,1.5,1.5,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5,-1.5,-1.5,-1.5,-1.5])
    dipole_east=dip_sep*numpy.array([-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5])
    dipole_z=dip_sep*numpy.zeros(dipole_north.shape)
    
    # loop over dipoles
    array_factor=0.0

    for i in xrange(4):
        for j in xrange(4):
            k=4*j+i
            # relative dipole phase for a source at (theta,phi)
            phase=numpy.exp((1j)*2*math.pi/lam*(dipole_east[k]*projection_east + dipole_north[k]*projection_north +
                                                dipole_z[k]*projection_z-delays[k]*c*delay_int))
            array_factor+=phase/16.0

    ground_plane=2*numpy.sin(2*math.pi*dipheight/lam*numpy.cos(theta))
    # make sure we filter out the bottom hemisphere
    ground_plane*=(theta<=math.pi/2)
    # normalize to zenith
    if (zenithnorm):
        ground_plane/=2*numpy.sin(2*math.pi*dipheight/lam)


    # this is a guess
    haAntennaZenith=0
    decAntennaZenith=dec-numpy.radians(mwa.lat)

    sinDecAntennaZenith=numpy.sin(decAntennaZenith)
    cosDecAntennaZenith=numpy.cos(decAntennaZenith)
    sinDec,cosDec=numpy.sin(dec), numpy.cos(dec)
    sinHa=numpy.sin(ha-haAntennaZenith)
    cosHa=numpy.cos(ha-haAntennaZenith)
    
    rot=numpy.array([cosHa,
                     sinDec*sinHa,
                     -sinDecAntennaZenith*sinHa,
                     cosDecAntennaZenith*cosDec+sinDecAntennaZenith*sinDec*cosHa])
    
    gain=rot*array_factor*ground_plane
    return gain


#gainXX,gainYY=primary_beam.MWA_Tile_analytic(theta, phi,
#                                             freq=frequency*1e6,
#                                             delays=delays)

gain=MWA_Jones_Analytic(theta, phi,
			numpy.radians(HA), numpy.radians(Decnow),
			freq=frequency,
			delays=delays)
print gain
gain=gain.reshape((2,2))
B=numpy.linalg.inv(gain)
Bprime=B.transpose().conj()

data=pyfits.open(filename)

j=numpy.sqrt(numpy.complex(-1))
XX=data['XX-REAL'].data+j*data['XX-IMAG'].data
YY=data['YY-REAL'].data+j*data['YY-IMAG'].data
XY=data['XY-REAL'].data+j*data['XY-IMAG'].data
YX=data['YX-REAL'].data+j*data['YX-IMAG'].data
I=numpy.zeros(XX.shape)
Q=numpy.zeros(XX.shape)
U=numpy.zeros(XX.shape)
V=numpy.zeros(XX.shape)

for t in xrange(XX.shape[0]):
	for f in xrange(XX.shape[1]):
		#D=numpy.reshape(numpy.array([[XX[t,f].real,XY[t,f],
		#XY[t,f].conj(),YY[t,f].real]]),(2,2))
		D=numpy.reshape(numpy.array([[XX[t,f].real,XY[t,f],
					      YX[t,f],YY[t,f].real]]),(2,2))
		Out=numpy.dot(numpy.dot(B,D),Bprime)
		if t==5 and f==0:
			out=Out
		I[t,f]=0.5*(Out[1,1].real+Out[0,0].real)
		Q[t,f]=0.5*(Out[1,1].real-Out[0,0].real)
		U[t,f]=0.5*(Out[0,1].real+Out[1,0].real)
		V[t,f]=0.5*(-Out[1,0].imag+Out[0,1].imag)


outfilename=filename.replace('.fits','_corrected.fits')
f=pyfits.PrimaryHDU(header=data[0].header)
l=[f]
l.append(pyfits.ImageHDU(data=I,
			 header=data[0].header))
l[-1].header['EXTNAME']='I'
l.append(pyfits.ImageHDU(data=Q,
			 header=data[0].header))
l[-1].header['EXTNAME']='Q'
l.append(pyfits.ImageHDU(data=U,
			 header=data[0].header))
l[-1].header['EXTNAME']='U'
l.append(pyfits.ImageHDU(data=V,
			 header=data[0].header))
l[-1].header['EXTNAME']='V'
if os.path.exists(outfilename):
	os.remove(outfilename)
hl=pyfits.HDUList(l)
hl.writeto(outfilename)

