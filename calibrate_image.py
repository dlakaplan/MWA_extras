import logging,datetime,math,sys,socket,os,json,shutil,io
from optparse import OptionParser,OptionGroup
import threading
import urllib2, urllib
import base64
import time
import subprocess
from astropy.table import Table,Column
import collections,glob,numpy

logging.basicConfig(format='# %(levelname)s:%(name)s: %(message)s')
logger=logging.getLogger('')
logger.setLevel(logging.INFO)


import mwapy
from mwapy import metadata

try:
    import drivecasa
    _CASA=True
except ImportError:
    _CASA=False


    

##################################################
def calibrate_casa(obsid):
    if not _CASA:
        logger.error("CASA operation not possible")
        return None
    basedir=os.path.abspath(os.curdir)
    try:
        casa = drivecasa.Casapy(casa_dir=casapy,
                                working_dir=basedir,
                                timeout=1200,
                                gui=False)
    except Exception, e:
        logger.error('Unable to instantiate casa:\n%s' % e)
        return None

    command=['from mwapy import ft_beam',
             """ft_beam.ft_beam(vis='%s.ms')""" % obsid]
    logger.info('Will run in casa:\n\t%s' % '\n\t'.join(command))
    result=casa.run_script(command)
    if len(result[1])>0:
        logger.error('CASA/ft_beam returned some errors:\n\t%s' % '\n\t'.join(result[1]))
        return None
    outfile=None
    for line in result[0]:
        if 'Created' in line:
            outfile=line.split()[1].replace('!','')
    if outfile is None:
        logger.error('No output created')
        return None
    # that file should be the same as the expected output
    if not outfile == '%s.cal' % obsid:
        logger.error('CASA calibration produced %s, but expected %s.cal' % (outfile,
                                                                            obsid))
        return None
    if not os.path.exists('%s.cal' % obsid):
        logger.error('CASA calibration command produced no ouptut')
        return None
    else:
        return '%s.cal' % obsid

##################################################
def applycal_casa(obsid, calfile):
    if not _CASA:
        logger.error("CASA operation not possible")
        return None
    basedir=os.path.abspath(os.curdir)
    try:
        casa = drivecasa.Casapy(casa_dir=casapy,
                                working_dir=basedir,
                                timeout=1200,
                                gui=False)
    except Exception, e:
        logger.error('Unable to instantiate casa:\n%s' % e)
        return None

    command=['applycal(vis="%s.ms", gaintable="%s")' % (obsid,
                                                        calfile)]
    logger.info('Will run in casa:\n\t%s' % '\n\t'.join(command))
    result=casa.run_script(command)
    # this doesn't really produce output
    if len(result[1])>0:
        logger.error('CASA/applycal returned some errors:\n\t%s' % '\n\t'.join(result[1]))
        return None
    return True

##################################################
def extract_calmodel(filename, sourcename):
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
            if d[0]=='name' and sourcename in d[1]:
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
def write_calmodelfile(calibrator_name):
    
    # get a model file
    calibrator_modeldata=extract_calmodel(calmodelfile,
                                          calibrator_name)
    if calibrator_modeldata is None:
        logger.error('No calibrator data found')
        return None
    outputcalmodelfile='%s_%s.model' % (calibrator_name,
                                        datetime.datetime.now().strftime('%Y%m%d'))
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
def calibrate_anoko(obsid,outputcalmodelfile):
    calfile='%s.cal' % obsid
    calibratecommand=[anokocalibrate,
                      '-j',
                      str(ncpus),
                      '-m',
                      outputcalmodelfile,
                      '%s.ms' % obsid,
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
def applycal_anoko(obsid, calfile):
    applycalcommand=[anokoapplycal,
                     '-copy',
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
        
    

##################################################                                                                                                                             
casapy='/usr/local/casapy'
calmodelfile='/home/kaplan/MWA_extras/model_a-team.txt'
anokocalibrate='~kaplan/mwa/anoko/mwa-reduce/build/calibrate'
anokoapplycal='~kaplan/mwa/anoko/mwa-reduce/build/applysolutions'
ncpus=16

if not os.environ.has_key('MWA_CODE_BASE'):
    logger.error('Environment variable $MWA_CODE_BASE is not set; please set and re-run')
    sys.exit(1)


clobber=False
calibration='anoko'
files=sorted(glob.glob('*.ms'))

# figure out which if any is a calibrator
# and which sources it would apply to
observations=[]
observation_data=numpy.zeros((len(files),),
                             dtype=[('obsid','i4'),
                                    ('iscalibrator',numpy.bool),
                                    ('channel','i4'),
                                    ('calibrator','i4'),
                                    ('calibratorsource','a20'),
                                    ('calibratorfile','a30')])

for i in xrange(len(files)):
    file=files[i]
    observations.append(metadata.MWA_Observation(int(file.split('.')[0])))
    observation_data[i]['obsid']=observations[-1].observation_number
    observation_data[i]['iscalibrator']=observations[-1].calibration
    observation_data[i]['channel']=observations[-1].center_channel
    if observation_data[i]['iscalibrator']:
        observation_data[i]['calibrator']=observation_data[i]['obsid']
        observation_data[i]['calibratorsource']=''.join(observations[-1].calibrators)
    

calibrators=numpy.where(observation_data['iscalibrator'])[0]
notcalibrators=numpy.where(~observation_data['iscalibrator'])[0]
if len(calibrators)==0:
    logger.error('No calibrators identified')
    sys.exit(1)

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
            closest=observation_data[dt==dt.min()]
            logger.info('Will use %d (separation=%d s) for calibration' % (closest['obsid'],
                                                           numpy.abs(closest['obsid']-observation_data[i]['obsid'])))
            observation_data[i]['calibrator']=closest['obsid']
        else:
            logger.info('Will use %d for calibration' % observation_data[good]['obsid'])
            observation_data[i]['calibrator']=observation_data[good]['obsid']
            
    else:
        logger.info('For observation %d (channel=%d) identified no matching calibrator observations' % (observation_data[i]['obsid'],
                                                                                                        observation_data[i]['channel']))
                                                                                                            
for i in calibrators:
    if calibration=='anoko':
        calibrator_name=observation_data[i]['calibratorsource']
        outputcalmodelfile=write_calmodelfile(calibrator_name)
        if outputcalmodelfile is None:
            sys.exit(1)
        
    calfile='%d.cal' % observation_data[i]['obsid']
    observation_data[i]['calibratorfile']=calfile
    if os.path.exists(calfile):
        if clobber:
            logger.warning('Calibration file %s exists but clobber=True; deleting...' % calfile)
            shutil.rmtree(calfile)
        else:
            logger.warning('Calibration file %s exists and clobber=False; continuing...' % calfile)
            continue
    if calibration=='anoko':
        result=calibrate_anoko(observation_data[i]['obsid'],outputcalmodelfile)
    elif calibration=='casa':        
        result=calibrate_casa(observation_data[i]['obsid'])
    if result is not None:
        calfile=result
        logger.info('Wrote %s' % calfile)
    else:
        sys.exit(1)
                       
for i in xrange(len(observation_data)):
    if not observation_data[i]['iscalibrator']:
        cal_touse=numpy.where(observation_data[i]['calibrator']==observation_data['obsid'])[0][0]
    else:
        cal_touse=i
    if not os.path.exists(observation_data[cal_touse]['calibratorfile']):
        logger.error('Cannot find calibration file %s for observation %s' % (observation_data[cal_touse]['calibratorfile'],
                                                                             observation_data[i]['obsid']))
        sys.exit(1)
    calibratorfile=observation_data[cal_touse]['calibratorfile']
    logger.info('Will calibrate %s.ms with %s' % (observation_data[i]['obsid'],
                                                  observation_data[cal_touse]['calibratorfile']))

    if calibration=='anoko':
        result=applycal_anoko(observation_data[i]['obsid'],calibratorfile)
    elif calibration=='casa':
        result=applycal_casa(observation_data[i]['obsid'],calibratorfile)
    if result is None:
        sys.exit(1)

