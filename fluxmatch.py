import matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot as plt
import numpy,os,sys,logging
from optparse import OptionParser
import astropy
from astropy import constants as c, units as u
from astropy.table import Table,Column
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import aegean,BANE,MIMAS
from AegeanTools.regions import Region

import mwapy
from mwapy.pb import make_beam
from mwapy import metadata

# configure the logging
logging.basicConfig(format='# %(levelname)s:%(name)s: %(message)s')
logger=logging.getLogger('fluxmatch')
logger.setLevel(logging.WARNING)

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
def find_sources_in_image(imagename, max_summits=5, innerclip=10, usebane=True,
                          region=None, cores=16):
    """
    sources,rmsimage,bgimage=find_sources_in_image(imagename, max_summits=5, innerclip=10, usebane=True, region=None, cores=16)
    runs aegean.find_sources_in_image
    but first runs BANE to get the BG/rms estimates
    if region is supplied (.mim format) only sources inside that will be identified
    """
    if usebane:
        rmsimage, bgimage=get_rms_background(imagename, cores=cores)
    else:
        rmsimage,bgimage=None,None

    sources=aegean.find_sources_in_image(imagename,
                                         max_summits=max_summits,
                                         innerclip=innerclip,
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
    logger.info('Creating primary beam for %s' % image)
    out=make_beam.make_beam(image, ext=0,
                            delays=delays,
                            analytic_model=True,
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
              fluxcolumn=None,
              fluxerrcolumn=None,
              nsigma=10,
              rmsfactor=3,
              matchradius=120,
              rejectsigma=3,
              maxdistance=20,
              minbeam=0.5,
              psfextent=1.1,
              limit=10,
              refineposition=False,
              update=False,
              plot=True,
              region=True,
              cores=1):
    """
    catalog='GLEAMIDR3.fits',
    fluxcolumn=None,
    fluxerrcolumn=None,
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
    max ratio of new to old fluxes (or reciprocal)
    limit=10,
    """

    if not isinstance(matchradius,astropy.units.quantity.Quantity):
        matchradius=matchradius*u.arcsec
    if not isinstance(maxdistance,astropy.units.quantity.Quantity):
        maxdistance=maxdistance*u.deg        

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
                                                     innerclip=nsigma,
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
    if 'RA' in fimage[0].header.keys():
        pointingcenter=SkyCoord(fimage[0].header['RA'],fimage[0].header['DEC'],
                                unit=('deg','deg'))
    else:
        # get the pointing center from the metadata
        logger.warning('Pointing metadata not present in header; retrieving...')
        try:
            obs=metadata.MWA_Observation(fimage[0].header['GPSTIME'])
            logger.info('Found pointing center %f,%f' % (obs.RA,obs.Dec))
            pointingcenter=SkyCoord(obs.RA,obs.Dec,
                                    unit=('deg','deg'))
        except:
            logger.warning('Using CRVAL1/CRVAL2 for pointing center')
            pointingcenter=SkyCoord(fimage[0].header['CRVAL1'],fimage[0].header['CRVAL2'],
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

    if '.fits' in catalog:
        # this seems to be faster than going straight to the Table.read()
        try:
            fcatalog=fits.open(catalog)
        except:
            logger.error('Unable to open FITS catalog %s' % catalog)
            return None
        catalogTable=Table(fcatalog[1].data)
    else:
        try:
            catalogTable=Table.read(catalog)
        except:
            logger.error('Unable to read catalog %s' % catalog)
            return None        
    bandfrequencies=numpy.array([int(s.split('_')[-1]) for s in numpy.array(catalogTable.colnames)[numpy.nonzero(numpy.array([('int_flux' in c) and not ('deep' in c) and not ('wide' in c) for c in catalogTable.colnames]))[0]]])
    
    if len(bandfrequencies)>0:
        # find the indices of the bands just above and below the observation
        # linearly weight the fluxes just above and below to match
        # the observation frequency
        indexplus=(bandfrequencies>=frequency/1e6).nonzero()[0].min()
        indexminus=(bandfrequencies<frequency/1e6).nonzero()[0].max()
        logger.info('Observation frequency of %.1f MHz: interpolating between %d MHz and %d MHz' % (frequency/1e6,bandfrequencies[indexminus],bandfrequencies[indexplus]))
    
        weightplus=(frequency/1e6-bandfrequencies[indexminus])/(bandfrequencies[indexplus]-bandfrequencies[indexminus])
        weightminus=1-weightplus
        gleamflux=catalogTable['int_flux_%03d' % bandfrequencies[indexminus]]*weightminus+catalogTable['int_flux_%03d' % bandfrequencies[indexplus]]*weightplus
        try:
            gleamfluxerr=numpy.sqrt((catalogTable['err_fit_flux_%03d' % bandfrequencies[indexminus]]*weightminus)**2+(catalogTable['err_fit_flux_%03d' % bandfrequencies[indexplus]]*weightplus)**2)
        except KeyError:
            gleamfluxerr=numpy.sqrt((catalogTable['err_int_flux_%03d' % bandfrequencies[indexminus]]*weightminus)**2+(catalogTable['err_int_flux_%03d' % bandfrequencies[indexplus]]*weightplus)**2)
    else:
        logger.warning('Could not identify GLEAM band fluxes')
        if fluxcolumn is None:
            logger.error('Could not identify flux columns to use')
            return None            
        if fluxcolumn in catalogTable.colnames and fluxerrcolumn in catalogTable.colnames:
            logger.warning('Using %s and %s columns' % (fluxcolumn,fluxerrcolumn))
            gleamflux=catalogTable[fluxcolumn]
            gleamfluxerr=catalogTable[fluxerrcolumn]
        else:
            logger.error('Could not identify flux columns to use')
            return None

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
    try:
        sourcesTable.add_column(Column(catalogTable['psf_a_%03d' % bandfrequencies[indexplus]][idx] * catalogTable['psf_b_%03d' % bandfrequencies[indexplus]][idx],
                                       name='PSFAREA'))
        sourcesTable.add_column(Column(catalogTable['a_%03d' % bandfrequencies[indexplus]][idx] * catalogTable['b_%03d' % bandfrequencies[indexplus]][idx],
                                       name='SOURCEAREA'))
    except:
        pass
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
        good=good & (sourcesTable['GLEAMFlux']>=sourcesTable['IntFlux'][good].min())
        logger.info('%04d/%04d sources match all cuts' % (good.sum(),
                                                          len(good)))
        if good.sum()<5:
            logger.error('Insufficient sources for flux scaling')
            return None

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
    if refineposition:
        sourcesTable.meta['rashift']=dRA[good].mean()*3600
        sourcesTable.meta['decshift']=dDEC[good].mean()*3600
    if os.path.exists(outbase + '_fluxmatch.hdf5'):
        os.remove(outbase + '_fluxmatch.hdf5')
    sourcesTable.write(outbase + '_fluxmatch.hdf5',path='data')
    logger.info('Wrote %s_fluxmatch.hdf5' % outbase)

    if region:
        outreg=outbase + '_fluxmatch.reg'
        if os.path.exists(outreg):
            os.remove(outreg)
        foutreg=open(outreg,'w')
        for i in xrange(len(sourcesTable)):
            if sourcesTable[i]['GOOD']:
                foutreg.write('icrs;circle(%f,%f,60") # text={%03d} color={green}\n' % (sourcesTable[i]['RA'],
                                                                                        sourcesTable[i]['Dec'],
                                                                                        i))
            else:
                foutreg.write('icrs;box(%f,%f,60",60",0) # text={%03d} color={red}\n' % (sourcesTable[i]['RA'],
                                                                                         sourcesTable[i]['Dec'],
                                                                                         i))
        logger.info('Wrote %s' % outreg)
        foutreg.close()

    if update:
        if fittedratio > limit or fittedratio < 1.0/limit:
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

        imagename=image.replace('_','\_')

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
        plt.axis([0.1,100,0.1,100])
        plt.xlabel('Flux Density in %s (Jy)' % catalog,fontsize=16)
        plt.ylabel('Flux Density in %s (Jy)' % imagename,fontsize=16)
        plt.gca().tick_params(labelsize=16)
        plt.savefig('%s_fluxflux.pdf' % outbase)
        logger.info('Wrote %s_fluxflux.pdf' % outbase)

        plt.clf()
        plt.hist(ratio[good],30)
        plt.xlabel('Flux Density in %s / Flux Density in %s' % (imagename,catalog),
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
        plt.xlabel('$\\alpha$(%s)-$\\alpha$(%s)' % (imagename,catalog),fontsize=16)
        plt.ylabel('$\\delta$(%s)-$\\delta$(%s)' % (imagename,catalog),fontsize=16)
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
        plt.xlabel('Separation %s vs. %s (arcsec)' % (imagename,catalog),
                   fontsize=16)
        plt.ylabel('Number of sources',fontsize=16)
        plt.gca().tick_params(labelsize=16)
        plt.savefig('%s_separation.pdf' % outbase)
        logger.info('Wrote %s_separation.pdf' % outbase)    

        
    return fittedratio, fittedratioerr, chisq, ndof, fitres[0], fitres[1]

######################################################################        

def main():
    
    usage="Usage: %prog [options] <files>\n"
    usage+="\tCompare MWA image against GLEAM catalog and determine flux scaling\n"
    usage+="\tSelection of comparison sources depends on:\n"
    usage+="\t\tSeparation\n\t\tPoint source extent\n\t\tLocal RMS\n"
    usage+="\t\tPrimary beam power\n\t\tDistance from pointing center\n"
    usage+="\tCan optionally refine positions, update image, make diagnostic plots\n"

    parser = OptionParser(usage=usage,version=mwapy.__version__ + ' ' + mwapy.__date__)
    parser.add_option('-c','--catalog',dest='catalog',default='GLEAMIDR3.fits',
                      help='Location of GLEAM catalog file [default=%default]')
    parser.add_option('--fluxcol',dest='fluxcolumn',default='FLUX',
                      help='Column for flux density if not GLEAM standard [default=%default]')
    parser.add_option('--fluxerrcol',dest='fluxerrcolumn',default='FLUXERR',
                      help='Column for flux density errors if not GLEAM standard [default=%default]')
    parser.add_option('--nsigma',dest='nsigma',default=10,type='float',
                      help='Threshold in sigma for source finding [default=%default]')
    parser.add_option('--match',dest='matchradius',default=60,type='float',
                      help='Matching radius in arcsec [default=%default]')
    parser.add_option('--rmsfactor',dest='rmsfactor',default=3,type='float',
                      help='Max ratio of local RMS to min RMS [default=%default]')
    parser.add_option('--maxdistance',dest='maxdistance',default=20,type='float',
                      help='Max distance from pointing center in degrees [default=%default]')
    parser.add_option('--minbeam',dest='minbeam',default=0.5,type='float',
                      help='Minimum primary beam power [default=%default]')
    parser.add_option('--psfextent',dest='psfextent',default=1.1,type='float',
                      help='Max source / PSF extent [default=%default]')
    parser.add_option('--rejectsigma',dest='rejectsigma',default=3,type='float',
                      help='Sigma clipping threshold [default=%default]')
    parser.add_option('--limit',dest='ratiolimit',default=10,type='float',
                      help='Max ratio of new to old fluxes (or old to new) [default=%default]')
    parser.add_option('--update',action="store_true",dest='update',default=False,
                      help="Update original FITS image?")
    parser.add_option('--refineposition',dest='refineposition',default=False,action='store_true',
                      help="Refine positions?")
    parser.add_option('--plot',action="store_true",dest="plot",default=False,
                      help="Save diagnostic plots?")
    parser.add_option('--region',action="store_true",dest="region",default=False,
                      help="Save ds9 region?")
    parser.add_option('-m','--cores',default=1,type='int',dest='cores',
                      help='Number of cores for Aegean [default=%default]')
    parser.add_option('-v','--verbose',action="store_true",dest="verbose",default=False,
                      help="Increase verbosity of output")
    (options, args) = parser.parse_args()

        
    if (options.verbose):
        logger.setLevel(logging.INFO)

    for file in args:
        out=fluxmatch(file,
                      catalog=options.catalog,
                      fluxcolumn=options.fluxcolumn,
                      fluxerrcolumn=options.fluxerrcolumn,
                      nsigma=options.nsigma,
                      matchradius=options.matchradius,
                      rmsfactor=options.rmsfactor,
                      maxdistance=options.maxdistance,
                      minbeam=options.minbeam,
                      psfextent=options.psfextent,
                      limit=options.ratiolimit,
                      rejectsigma=options.rejectsigma,
                      update=options.update,
                      refineposition=options.refineposition,
                      plot=options.plot,
                      region=options.region)


        
    sys.exit(0)
    
################################################################################

if __name__=="__main__":
    main()
                
