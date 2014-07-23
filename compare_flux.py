import pyfits,numpy,numpy.ma

D=[]
D2=[]
for i in xrange(16):      
    f=pyfits.open('1070380064_autopeel_PKS0436-%04d-XX-image.fits' % i)
    D.append(f[0].data[0,0,65,65])

f=pyfits.open('1070380064_autopeel_PKS0436_phased.fits')
d=f['XX-REAL'].data
d2=numpy.ma.array(d,mask=d==0)
for i in xrange(16):      
    D2.append(d2[:,(i*12):((i+1)*12)].mean())

x=numpy.arange(16)
