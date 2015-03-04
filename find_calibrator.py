import logging,datetime,math
import optparse

# configure the logging                                                                                                                                                        
logging.basicConfig(format='# %(levelname)s:%(name)s: %(message)s')
logger=logging.getLogger('metadata')
logger.setLevel(logging.WARNING)

import mwapy
from mwapy import metadata
import sys,os
from astropy.time import Time,TimeDelta
from astropy import constants as c, units as u
from astropy.coordinates import SkyCoord
import numpy

##################################################
def issamenight(obsid1, obsid2):
    """
    issamenight(obsid1, obsid2)
    checks whether they are the same night
    assumes AWST=UTC+8
    """
    t1=Time(obsid1, format='gps', scale='utc')
    t2=Time(obsid2, format='gps', scale='utc')
    if int(t1.jd)==int(t2.jd) and t1.datetime.hour>=10 and t2.datetime.hour>=10 and t1.datetime.hour<=22 and t2.datetime.hour<=22:
        return True
    return False

##################################################


def find_calibrator(obsid, 
                    maxtimediff=TimeDelta(1*u.d),
                    maxseparation=180*u.deg,
                    sourcename=None,
                    matchproject=True,
                    matchnight=True,
                    priority='time',
                    all=False):
    """
    obsid, timediff, distance=find_calibrator(obsid, 
                    maxtimediff=TimeDelta(1*u.d),
                    maxseparation=180*u.deg,
                    sourcename=None,
                    matchproject=True,
                    matchnight=True,
                    priority='time',
                    all=False):

    timediff is in days
    distance is in degrees
    """
    assert priority in ['time','distance']

    starttime=Time(obsid, format='gps',scale='utc')-maxtimediff
    stoptime=Time(obsid, format='gps',scale='utc')+maxtimediff

    try:
        baseobs=metadata.MWA_Observation(obsid)
    except Exception,e:
        logger.error('Cannot fetch info for observation %d:\n\t%s' % (obsid, e))
        return None

    if baseobs is None:
        logger.error('Cannot fetch info for observation %d:\n\t%s' % (obsid, e))
        return None

    try:
        basepointing=SkyCoord(baseobs.ra_phase_center, baseobs.dec_phase_center, 
                              unit='deg', frame='icrs')
    except:
        basepointing=SkyCoord(baseobs.RA, baseobs.Dec,
                          unit='deg', frame='icrs')

    basetime=Time(obsid, format='gps',scale='utc')
    
    if sourcename is not None and len(sourcename)>0:
        sourcename='%' + sourcename + '%'
    else:
        sourcename='%'
    if matchproject:
        results=metadata.fetch_observations(mintime=int(starttime.gps)-1,
                                            maxtime=int(stoptime.gps)-1,
                                            projectid=baseobs.projectid,
                                            calibration=1,
                                            obsname=sourcename,
                                            cenchan=baseobs.center_channel)
    else:
        results=metadata.fetch_observations(mintime=int(starttime.gps)-1,
                                            maxtime=int(stoptime.gps)-1,
                                            calibration=1,
                                            obsname=sourcename,
                                            cenchan=baseobs.center_channel)
    if len(results)==0:
        return None

    goodlist=numpy.zeros((len(results),),
                         dtype=[('obsid','i4'),
                                ('distance','f4'),
                                ('timediff','f4'),
                                ('good',numpy.bool)])

    for i in xrange(len(results)):
        result=results[i]
        obs=metadata.MWA_Observation_Summary(result)
        pointing=SkyCoord(obs.ra,obs.dec,unit='deg', frame='icrs')
        separation=basepointing.separation(pointing)
        good=True
        good=good and separation < maxseparation
        good=good and Time(obs.obsid,format='gps',scale='utc')-basetime<maxtimediff
        good=good and (not matchnight or issamenight(baseobs.observation_number, obs.obsid))
        goodlist[i]['obsid']=obs.obsid
        goodlist[i]['good']=good
        goodlist[i]['distance']=separation.deg
        goodlist[i]['timediff']=numpy.abs((Time(obs.obsid,format='gps',scale='utc')-basetime).jd)
        
    if priority=='time':
        goodlist=numpy.sort(goodlist, order='timediff')
    elif priority=='distance':
        goodlist=numpy.sort(goodlist, order='distance')
    if not all:
        return (goodlist[0]['obsid'],goodlist[0]['timediff'],goodlist[0]['distance'])
    else:
        return (list(goodlist['obsid']),list(goodlist['timediff']),list(goodlist['distance']))

##################################################
def main():

    usage="Usage: %prog [options] <obsid>\n"
    o = optparse.OptionParser(usage=usage,version=mwapy.__version__ + ' ' + mwapy.__date__)
    o.add_option('--separation',dest='separation',default=180,
                 type='float',
                 help='Maximum separation (deg) [default=%default]')
    o.add_option('--timediff',dest='timediff',default='1d',
                 help='Maximum time difference (d, h, m, or s) [default=%default]')
    o.add_option('--matchproject',dest='matchproject',default=False,
                 action='store_true',
                 help='Match project to original ObsID?')
    o.add_option('--matchnight',dest='matchnight',default=False,
                 action='store_true',
                 help='Match night to original ObsID?')
    o.add_option('--source',dest='source',default=None,
                 help='Calibrator source name [default=None]')
    o.add_option('--priority',dest='priority',default='time',
                 type='choice',
                 choices=['time','distance'],
                 help='Return the closest in time or distance? [default=%default]')
    o.add_option('--all',dest='all',default=False,
                 action='store_true',
                 help='Return all possible matches?')
    o.add_option('-v','--verbose',dest='verbose',default=False,
                 action='store_true',
                 help='Give verbose output?')

    options, args = o.parse_args()

    if len(args)==0:
        logger.error('Must specify >1 obsids')
        sys.exit(1)

    if 'd' in options.timediff:
        maxtimediff=TimeDelta(float(options.timediff[:-1])*u.d)
    elif 'h' in options.timediff:
        maxtimediff=TimeDelta(float(options.timediff[:-1])*u.h)
    elif 'm' in options.timediff:
        maxtimediff=TimeDelta(float(options.timediff[:-1])*u.m)
    elif 's' in options.timediff:
        maxtimediff=TimeDelta(float(options.timediff[:-1])*u.s)
    maxseparation=options.separation*u.deg
    matchproject=options.matchproject
    matchnight=options.matchnight
    priority=options.priority

    for obsid in args:
        result=find_calibrator(int(obsid),
                               maxtimediff=maxtimediff,
                               maxseparation=maxseparation,
                               matchproject=matchproject,
                               matchnight=matchnight,
                               sourcename=options.source,
                               priority=priority,
                               all=options.all)
        if result is None:
            print "None"
        elif not options.all:
            print '# Obsid\t\tCalObsid\tTimediff(day)\tDistance(deg)'
            print '%s\t%s\t%.3f\t\t%.1f' % (obsid,result[0],result[1],result[2])
            if options.verbose:
                o=metadata.MWA_Observation(int(result[0]))
                print o
        else:
            print '# Obsid\t\tCalObsid\tTimediff(day)\tDistance(deg)'
            for i in xrange(len(result[0])):
                print '%s\t%s\t%.3f\t\t%.1f' % (obsid,result[0][i],
                                                result[1][i],result[2][i])
                if options.verbose:
                    o=metadata.MWA_Observation(int(result[0][i]))
                    print o
            
    sys.exit(0)
################################################################################                                                                                               

if __name__=="__main__":
    main()


