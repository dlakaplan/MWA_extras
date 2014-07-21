import tpipe
import numpy as n

class pipe_mwa(tpipe.MiriadReader, tpipe.ProcessByDispersion2):
    """ Create pipeline object for reading MWA data in (as Miriad data) and doing dispersion-based analysis.
    This version uses optimized code for dedispersion, which also has some syntax changes.
    nints is the number of integrations to read.
    nskip is the number of integrations to skip before reading.
    nocal,nopass are options for applying calibration while reading Miriad data.
    Can also set some parameters as key=value pairs.
    """

    def __init__(self, file, profile='default', nints=1024, nskip=0, nocal=False,
                 nopass=False, selectpol=[], **kargs):
        self.set_profile(profile=profile)
        self.set_params(**kargs)
        self.read(file=file, nints=nints, nskip=nskip, nocal=nocal,
                  nopass=nopass, 
                  selectpol=selectpol)
        self.prep(deleteraw=True)


    def flag_antenna(self, antenna):
        """
        flag_antenna(self,antenna)
        """
        baseline_to_flag=n.any(self.blarr==antenna,axis=1)
        old_flagged=self.data.mask.sum()
        self.data.mask[:,baseline_to_flag]=True
        new_flagged=self.data.mask.sum()
        print "Flagged antenna %d: %d visibilities flagged" % (antenna, new_flagged-old_flagged)
        return new_flagged-old_flagged

    def flag_antennas(self, antennas):
        """
        flag_antennas(self,antennas)
        """
        flagged=0
        for antenna in antennas:
            flagged+=self.flag_antenna(antenna)
        return flagged
            

    def get_uv_kernel(self, width, length, angle):
        """
        Takes parameters of a elliptical Gaussian that represents model of spatial distribution of source.
        Returns kernel in uv plane, that is, the evaluation of Fourier inverse of spatial model at uv sampling points.
        Added by DLK 2013-04-09
        """

        sigmax=length/2.0
        sigmay=width/2.0
        # this is the correlation matrix in real space
        # [[a,b],[b,c]]        
        a=n.cos(angle)**2/2/sigmax**2+n.sin(angle)**2/2/sigmay**2
        c=n.sin(angle)**2/2/sigmax**2+n.cos(angle)**2/2/sigmay**2
        b=-n.sin(2*angle)/4/sigmax**2+n.sin(2*angle)/4/sigmay**2
        # invert to get UV
        det=a*c-b**2
        ainv=c/det
        binv=-b/det
        cinv=a/det
        #kernel=n.zeros((self.nints,self.nbl,self.nchan))
        kernel=n.zeros_like(self.data)
        for i in xrange(self.nints):
            # form channel dependent uvw
            u_ch = n.outer(self.u[i], self.freq/self.freq_orig[0])
            v_ch = n.outer(self.v[i], self.freq/self.freq_orig[0])
            if not isinstance(sigmax,n.ndarray):            
                kernel[i]=(n.exp(-(ainv*u_ch**2+2*binv*u_ch*v_ch+cinv*v_ch**2))).reshape(kernel[i].shape)
            else:
                kernel[i]=(n.exp(-(ainv[i]*u_ch**2+2*binv[i]*u_ch*v_ch+cinv[i]*v_ch**2))).reshape(kernel[i].shape)                        
            
        return kernel

    def make_phasedbeam(self, kernel=None):
        """ Integrates data at dmtrack for each pair of elements in dmarr, time.
        Not threaded.  Uses dmthread directly.
        Stores mean of detected signal after dmtrack, effectively forming beam at phase center.
        Ignores zeros in any bl, freq, time.
        """

        self.phasedbeam = n.zeros((len(self.dmarr),len(self.reltime)), dtype='float64')

        for i in xrange(len(self.dmarr)):
            if kernel is None:
                self.phasedbeam[i] = self.dedisperse(dmbin=i).mean(axis=3).mean(axis=2).mean(axis=1).real               # dedisperse and mean
            else:
                # average over pol, then multiply by kernel, then average over other axes
                self.phasedbeam[i] = (self.dedisperse(dmbin=i).mean(axis=3)*kernel).mean(axis=2).mean(axis=1).real               # dedisperse and mean                
            print 'dedispersed for ', self.dmarr[i]

    def general_dedisperse(self, data, dmbin):
        """ Creates dedispersed visibilities integrated over frequency.
        Uses ur track for each dm, then shifts by tint. Faster than using n.where to find good integrations for each trial, but assumes int-aligned pulse.

        will do this for generic array
        """

        dddata = data.copy()
        twidth = self.twidths[dmbin]
        delay = self.delay[dmbin]

        # dedisperse by rolling time axis for each channel
        for i in xrange(len(self.chans)):
            if len(data.shape)==4:
                dddata[:,:,i,:] = n.roll(data[:,:,i,:], -delay[i], axis=0)
            elif len(data.shape)==3:
                dddata[:,:,i] = n.roll(data[:,:,i], -delay[i], axis=0)

        return dddata

    def dmtrack(self, dm = 0., t0 = 0., show=0, middle=True):
        """ Takes dispersion measure in pc/cm3 and time offset from first integration in seconds.
        t0 defined at first (unflagged) channel. Need to correct by flight time from there to freq=0 for true time.
        Returns an array of (timebin, channel) to select from the data array.
        """


        reltime = self.reltime
        chans = self.chans
        tint = self.inttime

        # given freq, dm, dfreq, calculate pulse time and duration
        pulset_firstchan = 4.2e-3 * dm * self.freq[len(self.chans)-1]**(-2)   # used to start dmtrack at highest-freq unflagged channel
        pulset_midchan = 4.2e-3 * dm * self.freq[len(self.chans)/2]**(-2)   # used to start dmtrack at highest-freq unflagged channel. fails to find bright j0628 pulse
        if middle:
            pulset = 4.2e-3 * dm * self.freq**(-2) + t0 - pulset_midchan  # time in seconds referenced to some frequency (first, mid, last)
        else:
            pulset = 4.2e-3 * dm * self.freq**(-2) + t0 - pulset_firstchan  # time in seconds referenced to some frequency (first, mid, last)
        pulsedt = n.sqrt( (8.3e-6 * dm * (1000*self.sdf) * self.freq**(-3))**2 + self.pulsewidth**2)   # dtime in seconds

        chantimes=n.outer((pulset+pulsedt).reshape((self.nchan,1)),n.ones((1,len(reltime))))
        inttimes=n.outer(n.ones((self.nchan,1)),self.reltime.reshape((1,len(reltime))))
        I=(chantimes>=(inttimes+tint/2.0))
        
        chanbin=self.chans
        timebin=I.sum(axis=1)

        track = (list(timebin), list(chanbin))

        if show:
            p.plot(track[0], track[1], 'w*')

        return track


    def filter_and_dedisperse(self, kernel=None):
        """
        phased_data, phased_dd_data=filter_and_dedisperse(self, kernel=None)
        returns phased-array data (with or without kernel)
        and a copy that has been dedispersed
        dedispersion is carried out after summing over baselines
        """
        
        # compute the tied-array beam
        # with or without a kernel

        if kernel is not None:
            dataph=((self.data*kernel).mean(axis=3).mean(axis=1)).real
        else:
            dataph=((self.data.mean(axis=3)).mean(axis=1)).real
        dddata=n.ma.zeros((len(self.time),len(self.chans),len(self.dmarr)))
        for dm in xrange(len(self.dmarr)):
            for i in xrange(len(self.chans)):
                dddata[:,i,dm] = n.roll(dataph[:,i],-self.delay[dm][i],axis=0)
        return dataph, dddata.mean(axis=1)
    
    def calc_sefd(self, eta=1.0, i1=0, i2=None):
        """
        sefd,sigma=pipe.calc_sefd(eta=1.0, i1=0, i2=None)
        calculates the SEFD and uncertainty for each antenna, channel, and polarization
        averaging over the time range i1:i2
        assumes overall efficiency eta

        based on Lu Feng's calculate_sefd.py
        updated to do more vectorized math

        output has dimensions [antennas, channels, polarizations]
        """
        if i2 is None:
            i2=-1

        # get the channel width and integration time
        t_int=(self.time[1]-self.time[0])*86400
        chan_width=(self.freq[1]-self.freq[0])*1e9
        
        # initialize output arrays
        sefd=n.ma.zeros((self.nants,self.nchan,self.npol))
        sigma=n.ma.zeros((self.nants,self.nchan,self.npol))
        # set svd cutoff based on machine precision
        svd_cond=self.nbl*n.finfo(n.float).eps  # svd cutoff (1.10134124043e-13)
        sefd.mask=(sefd.data>0)
        sigma.mask=(sefd.data>0)
        
        corr_matrix=n.zeros((self.nbl,self.nants))
        for b in xrange(self.nbl):
            ant1=self.blarr[b][0]
            ant2=self.blarr[b][1]
            # make sure we do not include autos
            if ant1 == ant2:
                continue
            i=n.where(self.ants==ant1)[0][0]
            j=n.where(self.ants==ant2)[0][0]
            corr_matrix[b,i]=1
            corr_matrix[b,j]=1                    

        # use single value decomposition (svd)
        # the code below reproduces linalg.pinv
        u,w,v=n.linalg.svd(corr_matrix,full_matrices=False)
        vT=v.transpose()
        uT=u.transpose()
        w[w < (n.max(w)*svd_cond)]=0  # svd cutoff
        w_inv=1/w
        w_inv[n.isinf(w_inv)]=0  # svd: replace 1/0 by 0

        for pol in xrange(self.npol):
            for ch in xrange(self.nchan):
                calvis_imag=self.data[i1:i2,:,ch,pol].imag

                stddev=calvis_imag.std(axis=0,ddof=1)
                stddev_vector=2*n.log(eta*stddev*n.sqrt(2*chan_width*t_int))
        
                # calculate sefd/error for all antennas
                udotb=n.dot(uT,stddev_vector)
                log_sefd=n.inner(vT,w_inv*udotb)
                log_sigma_sq=((vT*w_inv)**2).sum(axis=1)
                # system is linear in log-space
                # propagate error back into normal space
                sefd[:,ch,pol]=n.exp(log_sefd)
                sigma[:,ch,pol]=sefd[:,ch,pol]*n.sqrt(log_sigma_sq)
        
        sefd.mask[(sefd.data<10) | (sigma.data < 1e-3)]=True
        sigma.mask[(sefd.data<10) | (sigma.data < 1e-3)]=True
        return sefd, sigma
    
