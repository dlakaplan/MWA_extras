"""
Todo:

-quality metrics:
  - image rms
  - # of sources detected
-selfcal

To debug:

-table output of null strings
-table output of filenames

"""

import logging,datetime,math,sys,socket,os,json,shutil,io
from optparse import OptionParser,OptionGroup
import threading
import urllib2, urllib
import base64
import time
import collections
import subprocess
from astropy.table import Table,Column
import collections,glob,numpy
from astropy.io import fits

##############################
# Custom formatter
# http://stackoverflow.com/questions/1343227/can-pythons-logging-format-be-modified-depending-on-the-message-log-level
class MyFormatter(logging.Formatter):
    err_fmt  = '\n\n# %(asctime)s %(levelname)s:%(name)s: %(message)s\n\n'
    warning_fmt  = '\n# %(asctime)s %(levelname)s:%(name)s: %(message)s\n'
    other_fmt  = '# %(asctime)s %(levelname)s:%(name)s: %(message)s'

    def __init__(self, fmt=other_fmt):
        logging.Formatter.__init__(self, fmt)

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.WARNING:
            self._fmt = MyFormatter.warning_fmt

        elif record.levelno == logging.ERROR:
            self._fmt = MyFormatter.err_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._fmt = format_orig

        return result
    

# configure the logging
logFormatter=logging.Formatter('# %(asctime)s %(levelname)s:%(name)s: %(message)s')
logging.basicConfig(format='# %(asctime)s %(levelname)s:%(name)s: %(message)s')
logger = logging.getLogger('calibrate_image')

fmt = MyFormatter()
hdlr = logging.StreamHandler(sys.stdout)

hdlr.setFormatter(fmt)
logger.addHandler(hdlr)
logger.setLevel(logging.WARNING)


import mwapy
from mwapy import metadata

try:
    import drivecasa
    _CASA=True
except ImportError:
    _CASA=False

casapy='/usr/local/casapy'
calmodelfile='/home/kaplan/MWA_extras/model_a-team.txt'
anokocalibrate='~kaplan/mwa/anoko/mwa-reduce/build/calibrate'
anokoapplycal='~kaplan/mwa/anoko/mwa-reduce/build/applysolutions'
anokocatalog='~kaplan/mwa/anoko/mwa-reduce/catalogue/model-catalogue_new.txt'

##################################################
def makemetafits(obsid):        
    m=metadata.instrument_configuration(int(obsid))
    h=m.make_metafits()
    metafits='%s.metafits' % obsid
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

    if not _CASA:
        logger.error('requires drivecasa')
        return None

    try:
        casa = drivecasa.Casapy(casa_dir=casapy,
                                working_dir=os.path.abspath(os.curdir),
                                timeout=1200,
                                gui=False)
    except Exception, e:
        logger.error('Unable to instantiate casa:\n%s' % e)
        return None

    command=['import casac',
             'ms=casac.casac.ms()',
             'ms.open("%s")' % msfile,
             'print "chanwidth=%d" % (ms.getspectralwindowinfo()["0"]["ChanWidth"]/1e3)',
             'print "inttime=%f" % (ms.getscansummary()["1"]["0"]["IntegrationTime"])',
             't=casac.casac.table()',
             't.open("%s")' % msfile,
             'keys=t.getkeywords()',
             'print "\\n".join(["%s=%s" % (k,keys[k]) for k in keys.keys()])']
             
    logger.debug('Will run in casa:\n\t%s' % '\n\t'.join(command))
    result=casa.run_script(command)
    if len(result[1])>0:
        logger.error('CASA returned some errors:\n\t%s' % '\n\t'.join(result[1]))
        return None
    chanwidth=None
    inttime=None
    otherkeys={}
    for l in result[0]:
        if 'chanwidth' in l:
            chanwidth=float(l.split('=')[1])
        elif 'inttime' in l:
            inttime=float(l.split('=')[1])
        elif '=' in l:
            k,v=l.split('=')
            otherkeys[k]=v

    return chanwidth, inttime, otherkeys
    
##################################################
def check_calibrated(msfile):
    if not _CASA:
        logger.error('requires drivecasa')
        return None

    try:
        casa = drivecasa.Casapy(casa_dir=casapy,
                                working_dir=os.path.abspath(os.curdir),
                                timeout=1200,
                                gui=False)
    except Exception, e:
        logger.error('Unable to instantiate casa:\n%s' % e)
        return None

    command=['import casac',
             't=casac.casac.table()',
             't.open("%s")' % msfile,
             'print "CORRECTED_DATA" in t.colnames()']
             
    logger.debug('Will run in casa:\n\t%s' % '\n\t'.join(command))
    result=casa.run_script(command)
    if len(result[1])>0:
        logger.error('CASA returned some errors:\n\t%s' % '\n\t'.join(result[1]))
        return None
    for l in result[0]:
        if len(l)>0 and l=='True':
            return True
        if len(l)>0 and l=='False':
            return False
    return None
##################################################
def getcasaversion():
    if not _CASA:
        logger.error('requires drivecasa')
        return None

    try:
        casa = drivecasa.Casapy(casa_dir=casapy,
                                working_dir=os.path.abspath(os.curdir),
                                timeout=1200,
                                gui=False)
    except Exception, e:
        logger.error('Unable to instantiate casa:\n%s' % e)
        return None

    command=['pass']
             
    logger.debug('Will run in casa:\n\t%s' % '\n\t'.join(command))
    result=casa.run_script(command)
    if len(result[1])>0:
        logger.error('CASA returned some errors:\n\t%s' % '\n\t'.join(result[1]))
        return None
    return result[0]


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
    logger.debug('Will run in casa:\n\t%s' % '\n\t'.join(command))
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
    logger.debug('Will run in casa:\n\t%s' % '\n\t'.join(command))
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

######################################################################
def identify_calibrators(observation_data):
    calibrators=numpy.where(observation_data['iscalibrator'])[0]
    notcalibrators=numpy.where(~observation_data['iscalibrator'])[0]
    if len(calibrators)==0:
        logger.error('No calibrators identified')
        return None
    
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
                closest=observation_data[dt==dt.min()]
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
def autoprocess(self):
    autoprocesscommand=['autoprocess',
                        anokocatalog,
                        '%s.ms' % obsid]
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
            logger.error(l.rstrip())
        returncode=p.poll()
        if returncode is not None:
            break
        time.sleep(1)
    tasks=[]
    for l in output:
        if l.startswith('- '):
            if not '- No' in l:
                tasks.append(l[2:])
    
    if len(tasks)==0:
        logger.info('No autoprocess required')
        self.autoprocess='None'
        return None
    for i in xrange(len(tasks)):
        logger.info('autoprocess will do: %s' % tasks[i])

    autoprocesscommand=['autoprocess',
                        '-go',
                        anokocatalog,
                        '%s.ms' % obsid]
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
            logger.error(l.rstrip())
        returncode=p.poll()
        if returncode is not None:
            break
        time.sleep(1)

    self.autoprocess=';'.join(tasks)
    return True

        
######################################################################
class Observation(metadata.MWA_Observation):

    ##############################
    def __init__(self, obsid, outputdir='./', clobber=False, delete=False):
        self.obsid=obsid
        self.basedir=os.path.abspath(outputdir)
        self.clobber=clobber
        self.delete=delete

        self.metafits=None
        self.rawfiles=[]
        self.beamfiles=[]
        self.corrfiles=[]

        self.calibrator=None
        self.calmodelfile=None
        self.caltype=None
        self.calibratorfile=None
        self.inttime=0
        self.chanwidth=0
        self.otherkeys={}
        self.autoprocess=False

        if os.path.exists('%s.metafits' % self.obsid):
            self.metafits='%s.metafits' % self.obsid
        else:
            self.metafits=makemetafits(self.obsid)
            
        self.observation=metadata.MWA_Observation(self.obsid)

        if not os.path.exists('%s.ms' % self.obsid):
            logger.error('MS file %s.ms does not exist' % self.obsid)
        self.chanwidth, self.inttime, self.otherkeys=get_msinfo('%s.ms' % self.obsid)
        self.filestodelete=[]

    ##############################
    def __del__(self):
        if self.delete:
            for file in self.filestodelete:
                if os.path.exists(file):
                    try:
                        if os.path.isdir(file):
                            shutil.rmtree(file)
                        else:
                            os.remove(file)
                    except Exception,e:
                        logger.error('Unable to delete file %s:\n\t%s' % (file,e))


    ##############################
    def __getattr__(self, name):
        if self.__dict__.has_key(name):
            return self.__dict__[name]
        else:
            return self.observation.__dict__[name]

    ##############################
    def make_cal(self, caltype):
        if not self.calibration:
            return None
        if caltype=='anoko':
            calibrator_name=self.calibratorsource
            self.calmodelfile=write_calmodelfile(calibrator_name)
            if self.calmodelfile is None:
                return None
            self.filestodelete.append(self.calmodelfile)
        self.calibratorfile='%s.cal' % self.obsid
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
        if caltype=='anoko':
            result=calibrate_anoko(self.obsid,
                                   self.calmodelfile)
        elif caltype=='casa':        
            result=calibrate_casa(self.obsid)
        if result is not None:
            self.calibratorfile=result
            self.caltype=caltype
            logger.info('Wrote %s' % calibratorfile)
            return self.calibratorfile
        else:
            return None

    ##############################
    def calibrate(self, caltype='anoko', recalibrate=True):
        if not os.path.exists(self.calibratorfile):
            logger.error('Cannot find calibration file %s for observation %s' % (self.calibratorfile,
                                                                                 self.obsid))
            return False
        logger.info('Will calibrate %s.ms with %s' % (self.obsid,
                                                      self.calibratorfile))

        if check_calibrated('%s.ms' % self.obsid) and not recalibrate:
            logger.info('%s.ms is already calibrated and recalibrate is False' % self.obsid)
            return True
        if self.caltype=='anoko':
            result=applycal_anoko(self.obsid,self.calibratorfile)
        elif self.caltype=='casa':
            result=applycal_casa(self.obsid,self.calibratorfile)
        if result is None:
            return False

        if not check_calibrated('%s.ms' % self.obsid):
            logger.error('%s.ms does not appear to be calibrated; no CORRECTED_DATA file is present' % self.obsid)
            return False
        return True

    ##############################
    def image(self, 
              clean_weight='uniform',
              imagesize=2048,
              pixelscale=0.015,
              clean_iterations=100,
              clean_gain=0.1,
              clean_mgain=1.0,
              clean_minuv=0,
              clean_maxuv=0,
              fullpolarization=True,
              wsclean_arguments='',
              ncpus=32,
              memfraction=50):

        self.clean_weight=clean_weight
        self.imagesize=imagesize
        self.pixelscale=pixelscale
        self.clean_iterations=clean_iterations
        self.clean_gain=clean_gain
        self.clean_mgain=clean_mgain
        self.clean_minuv=clean_minuv
        self.clean_maxuv=clean_maxuv
        self.fullpolarization=fullpolarization
        self.wsclean_arguments=wsclean_arguments

        wscleancommand=['wsclean',
                        '-j',
                        str(ncpus),
                        '-mem',
                        str(memfraction),
                        '-name',
                        str(self.obsid),
                        '-weight',
                        self.clean_weight,
                        '-size',
                        '%d %d'% (self.imagesize,self.imagesize),
                        '-scale',
                        '%.4fdeg' % self.pixelscale,
                        '-niter',
                        str(self.clean_iterations),
                        '-gain',
                        str(self.clean_gain),
                        '-mgain',
                        str(self.clean_mgain)]
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
            wscleancommand.append('-joinpolarizations')
            for pol in ['XX','YY','XY','XYi']:
                expected_output.append('%s-%s-image.fits' % (self.obsid,pol))
        else:
            wscleancommand.append('-pol')
            wscleancommand.append('xx,yy')
            for pol in ['XX','YY']:
                expected_output.append('%s-%s-image.fits' % (self.obsid,pol))
        if len(self.wsclean_arguments)>0:
            wscleancommand+=self.wsclean_arguments.split()
        wscleancommand.append('%s.ms' % self.obsid)
        self.wscleancommand=wscleancommand
        logger.info('Will run:\n\t%s' % ' '.join(wscleancommand))
        p=subprocess.Popen(' '.join(wscleancommand),
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
        for j in xrange(len(expected_output)):
            file=expected_output[j]
            if not os.path.exists(file):
                logger.error('Expected wsclean output %s does not exist' % file)
                return None
            logger.info('Created wsclean output %s' % file)
            self.rawfiles.append(file)

            try:
                f=fits.open(file,'update')
            except Exception,e:
                logger.open('Unable to open image output %s for updating:\n\t%s' % (file,e))
                return None

            f[0].header['CALTYPE']=(self.calibration,'Calibration type')
            #f[0].header['CALIBRAT']=(observation_data[i]['calibrator'],
            #                         'Calibrator observation')
            f[0].header['CALFILE']=(self.calibratorfile,
                                    'Calibrator file')
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
            f[0].header['INTTIME']=(self.inttime,'[s] Integration time')
            f[0].header['FINECHAN']=(self.chanwidth,'[kHz] Fine channel width')
            if prev_inttime != f[0].header['INTTIME']:
                f[0].header['NSCANS']=int(f[0].header['NSCANS']*prev_inttime/f[0].header['INTTIME'])
                logger.info('New integration time (%.1f s) differs from previous integration time (%.1f s): changing number of scans to %d' % (f[0].header['INTTIME'],
                                                                                                                                               prev_inttime,
                                                                                                                                               f[0].header['NSCANS']))
            if prev_finechan != f[0].header['FINECHAN']:
                f[0].header['NCHANS']=int(f[0].header['NCHANS']*prev_finechan/f[0].header['FINECHAN'])
                logger.info('New final channel (%d kHz) differs from previous fine channel (%d kHz): changing number of channels to %d' % (f[0].header['FINECHAN'],
                                                                                                                                           prev_finechan,
                                                                                                                                           f[0].header['NCHANS']))
                
            
            f.verify('fix')
            f.flush()

            # add other output to cleanup list
            for ext in ['dirty','model','residual']:
                self.filestodelete.append(file.replace('-image','-%s' % ext))

        return self.rawfiles

    ##############################
    def makepb(self, beammodel='2014i'):
        assert beammodel in ['2013','2014','2014i']

        self.beammodel=beammodel
        beamcommand=['beam',
                     '-' + beammodel,
                     '-name',
                     'beam-%s' % self.obsid,
                     '-proto',
                     self.rawfiles[0],
                     '-ms',
                     '%s.ms' % self.obsid]
        
        expected_beamoutput=[]
        for pol in ['xxr','xxi','yyr','yyi','xyr','xyi','yxr','yxi']:
            expected_beamoutput.append('beam-%s-%s.fits' % (self.obsid,
                                                            pol))

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
        return expected_beamoutput

    ##############################
    def pbcorrect(self):
        pbcorrectcommand=['pbcorrect',
                          str(self.obsid),
                          'image.fits',
                          'beam-%s' % str(self.obsid),
                          str(self.obsid)]
        expected_pbcorroutput=[]
        if self.fullpolarization:
            for pol in ['I','Q','U','V']:
                expected_pbcorroutput.append('%s-%s.fits' % (self.obsid,
                                                             pol))
        else:
            for pol in ['I']:
                expected_pbcorroutput.append('%s-%s.fits' % (self.obsid,
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
                fraw=fits.open(self.rawfiles[0])
            except Exception,e:
                logger.error('Unable to open file %s:\n\t%s' % (self.rawfiles[0],e))

            for k in fraw[0].header.keys():
                if not k in f[0].header.keys():
                    f[0].header[k]=(fraw[0].header[k],
                                    fraw[0].header.comments[k])
            f.verify('fix')
            f.flush()

        if not self.fullpolarization:
            for pol in ['Q','U','V']:
                if os.path.exists('%s-%s.fits' % (self.obsid,pol)):
                    self.filestodelete.append('%s-%s.fits' % (self.obsid,pol))


            
        return self.corrfiles
            

##################################################                              
def main():
    
    times=collections.OrderedDict()
    times['start']=time.time()

    usage="Usage: %prog [options] <msfiles>\n"
    parser = OptionParser(usage=usage)

    correct_primarybeam=True
    parser.add_option('--caltype',dest='caltype',default='anoko',
                      type='choice',
                      choices=['anoko','casa'],
                      help='Type of calibration [default=%default]')
    parser.add_option('--recalibrate',dest='recalibrate',default=False,
                      action='store_true',
                      help='Recalibrate existing output?')
    parser.add_option('--autoprocess',dest='autoprocess',default=False,
                      action='store_true',
                      help='Run autoprocess?')
    parser.add_option('--size',dest='imagesize',default=2048,type='int',
                      help='Image size in pixels [default=%default]')
    parser.add_option('--scale',dest='pixelscale',default=0.015,type='float',
                      help='Pixel scale in deg [default=%default]')
    parser.add_option('--niter',dest='clean_iterations',default=100,
                      type='int',
                      help='Clean iterations [default=%default]')
    parser.add_option('--mgain',dest='clean_mgain',default=1,
                      type='float',
                      help='Clean major cycle gain [default=%default]')
    parser.add_option('--gain',dest='clean_gain',default=0.1,
                      type='float',
                      help='Clean gain [default=%default]')
    parser.add_option('--weight',dest='clean_weight',default='uniform',
                      help='Clean weighting [default=%default]')
    parser.add_option('--minuv',dest='clean_minuv',default=0,
                      type='float',
                      help='Minimum UV distance in wavelengths [default=%default]')
    parser.add_option('--maxuv',dest='clean_maxuv',default=0,
                      type='float',
                      help='Maximum UV distance in wavelengths [default=%default]')
    parser.add_option('--fullpol',dest='fullpolarization',default=False,
                      action='store_true',
                      help='Process full polarization (including cross terms)?')
    parser.add_option('--wscleanargs',dest='wsclean_arguments',default='',
                      help='Additional wsclean arguments')
    parser.add_option('--beam',dest='beammodel',default='2014i',type='choice',
                      choices=['2014i','2014','2013'],
                      help='Primary beam model [default=%default]')
    parser.add_option('--summaryname',dest='summaryname',default='summary.dat',
                      help='Name for output summary table [default=%default]')
    parser.add_option('--summaryformat',dest='summaryformat',default='ascii.commented_header',
                      help='Format for output summary table (from astropy.table) [default=%default]')   
    parser.add_option('--nocal',dest='nocal',default=False,
                      action='store_true',
                      help='Do not image the calibrator observations?')
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
    parser.add_option("-v", "--verbose", dest="loudness", default=0, action="count",
                      help="Each -v option produces more informational/debugging output")
    parser.add_option("-q", "--quiet", dest="quietness", default=0, action="count",
                      help="Each -q option produces less error/warning/informational output")

    (options, args) = parser.parse_args()
    loglevels = {0: [logging.DEBUG, 'DEBUG'],
                 1: [logging.INFO, 'INFO'],
                 2: [logging.WARNING, 'WARNING'],
                 3: [logging.ERROR, 'ERROR'],
                 4: [logging.CRITICAL, 'CRITICAL']}
    logdefault = 2    # WARNING
    level = max(min(logdefault - options.loudness + options.quietness,4),0)
    logger.setLevel(loglevels[level][0])
    logger.info('Log level set: messages that are %s or higher will be shown.' % loglevels[level][1])

    files=args
    if len(files)==0:
        logger.error('Must supply >=1 ms files to process')
        sys.exit(1)

    if not os.environ.has_key('MWA_CODE_BASE'):
        logger.error('Environment variable $MWA_CODE_BASE is not set; please set and re-run')
        sys.exit(1)

    logger.debug('Using mwapy version %s' % mwapy.__version__)
    result=subprocess.Popen(['wsclean','-version'],
                            stdout=subprocess.PIPE).communicate()[0].strip()
    wscleanversion=result.split('\n')[0].split('version')[1].strip()
    logger.debug('Using wsclean version %s' % wscleanversion)
    casapyversion='-'.join(os.path.realpath(casapy).split('-')[1:])
    logger.debug('Using casapy version %s' % casapyversion)

    # figure out which if any is a calibrator
    # and which sources it would apply to
    observations=collections.OrderedDict()
    observation_data=numpy.zeros((len(files),),
                                 dtype=[('obsid','i4'),
                                        ('metafits','a20'),
                                        ('inttime','f4'),
                                        ('chanwidth','f4'),
                                        ('ms_cotterversion','a20'),
                                        ('ms_metadataversion','a20'),
                                        ('ms_mwapyversion','a20'),
                                        ('iscalibrator',numpy.bool),
                                        ('channel','i4'),
                                        ('calibrator','i4'),
                                        ('calibratorsource','a20'),
                                        ('calibratorfile','a30'),
                                        ('clean_iterations','i4'),
                                        ('imagesize','i4'),
                                        ('pixelscale','f4'),
                                        ('clean_weight','a20'),
                                        ('clean_minuv','f4'),
                                        ('clean_maxuv','f4'),
                                        ('clean_gain','f4'),
                                        ('clean_mgain','f4'),
                                        ('fullpolarization',numpy.bool),
                                        ('wsclean_arguments','a60'),
                                        ('wsclean_command','a200'),
                                        ('rawimages','a30',(6,)),               
                                        ('beammodel','a10'),
                                        ('corrimages','a30',(4,)), 
                                    ])

    ##################################################
    # gather initial data
    ##################################################
    for i in xrange(len(files)):
        file=files[i]
        obsid=int(file.split('.')[0])
        observations[obsid]=Observation(obsid, clobber=options.clobber,
                                        delete=options.delete)
        if observations[obsid].inttime==0:
            sys.exit(1)        
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
                             
    calibrators, notcalibrators, cal_observations=identify_calibrators(observation_data)
    for i in notcalibrators:
        observation_data[i]['calibrator']=cal_observations[observation_data[i]['obsid']]
        observations[observation_data[i]['obsid']].calibrator=cal_observations[observation_data[i]['obsid']]
    times['init']=time.time()

    ##################################################
    # generate calibration solutions
    ##################################################
    for i in calibrators:
        result=observations[observation_data[i]['obsid']].make_cal(options.caltype)
        if result is None:
            sys.exit(1)
        observation_data[i]['calibratorfile']=result

    times['makecal']=time.time()        

    ##################################################
    # apply calibration
    ##################################################
    for i in xrange(len(observation_data)):
        if not observation_data[i]['iscalibrator']:
            cal_touse=numpy.where(observation_data[i]['calibrator']==observation_data['obsid'])[0][0]
        else:
            cal_touse=i
        observations[observation_data[i]['obsid']].calibratorfile=observation_data[cal_touse]['calibratorfile']
        
        result=observations[observation_data[i]['obsid']].calibrate(options.caltype,
                                                                    options.recalibrate)
        if result is False or result is None:
            sys.exit(1)
        observation_data[i]['calibratorfile']=observations[observation_data[i]['obsid']].calibratorfile

    times['applycal']=time.time()                       

    if options.autoprocess:
        results=observations[observation_data[i]['obsid']].autoprocess()

    times['autoprocess']=time.time()                       
    ##################################################
    # imaging
    ##################################################
    for i in xrange(len(observation_data)):
        if observations[observation_data[i]['obsid']].calibration and options.nocal:
            continue
        results=observations[observation_data[i]['obsid']].image(
            clean_weight=options.clean_weight,
            imagesize=options.imagesize,
            pixelscale=options.pixelscale,
            clean_iterations=options.clean_iterations,
            clean_gain=options.clean_gain,
            clean_mgain=options.clean_mgain,
            clean_minuv=options.clean_minuv,
            clean_maxuv=options.clean_maxuv,
            fullpolarization=options.fullpolarization,
            wsclean_arguments=options.wsclean_arguments,
            ncpus=options.ncpus,
            memfraction=options.memfraction)
        if results is None:
            sys.exit(1)
        for j in xrange(len(results)):
            file=results[j]
            try:
                f=fits.open(file,'update')
            except Exception,e:
                logger.open('Unable to open image output %s for updating:\n\t%s' % (file,e))
                sys.exit(1)

            f[0].header['CASAVER']=(casapyversion,'CASAPY version')
            f[0].header['MWAPYVER']=(mwapy.__version__,'MWAPY version')
            f[0].header['WSCLNVER']=(wscleanversion,'WSCLEAN version')
            f[0].header['COTTRVER']=(observation_data[i]['ms_cotterversion'],'Cotter version')
            f.flush()
            observation_data[i]['rawimages'][j]=results[j]

        observation_data[i]['clean_iterations']=options.clean_iterations
        observation_data[i]['clean_weight']=options.clean_weight
        observation_data[i]['clean_gain']=options.clean_gain
        observation_data[i]['clean_mgain']=options.clean_mgain
        observation_data[i]['imagesize']=options.imagesize
        observation_data[i]['pixelscale']=options.pixelscale
        observation_data[i]['clean_minuv']=options.clean_minuv
        observation_data[i]['clean_maxuv']=options.clean_maxuv
        observation_data[i]['fullpolarization']=options.fullpolarization
        observation_data[i]['wsclean_arguments']=options.wsclean_arguments
        observation_data[i]['wsclean_command']=' '.join(observations[observation_data[i]['obsid']].wscleancommand)

        results=observations[observation_data[i]['obsid']].makepb(options.beammodel)
        if results is None:
            sys.exit(1)
        observation_data[i]['beammodel']=options.beammodel

        results=observations[observation_data[i]['obsid']].pbcorrect()
        if results is None:
            sys.exit(1)
        for j in xrange(len(results)):
            observation_data[i]['corrimages'][j]=results[j]

   
    times['image']=time.time()
    try:
        observation_data_table=Table(observation_data)
    except Exception,e:
        logger.error('Unable to create Table for summary data:\n\t%s' % e)
        sys.exit(1)

    try:
        observation_data_table.write(options.summaryname,
                                     format=options.summaryformat)
        logger.info('Summary table written to %s' % options.summaryname)
    except Exception, e:
        logger.error('Unable to write summary table %s with format %s:\n%s' % (options.summaryname,
                                                                               options.summaryformat,
                                                                               e))
        sys.exit(1)

    times['end']=time.time()


    logger.info('Execution time: %d s (setup), %d s (make cal), %d s (apply cal), %d s (autoprocess), %d s (image) %d s (tidy) = %d s (total)' % (
        times['init']-times['start'],
        times['makecal']-times['init'],
        times['applycal']-times['makecal'],
        times['autoprocess']-times['applycal'],        
        times['image']-times['autoprocess'],        
        times['end']-times['image'],
        times['end']-times['start']))


    sys.exit(0)

                                                                
######################################################################

if __name__=="__main__":

    main()
