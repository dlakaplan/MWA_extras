import logging,logging.handlers,datetime,math,sys,socket,os,json,shutil,io,re
from optparse import OptionParser,OptionGroup
import threading
import urllib2, urllib
import base64
import time
import subprocess
import ssl
import uuid
import stat
from astropy.table import Table,Column
from astropy.time import Time
from astropy.coordinates import SkyCoord
import collections

if sys.version_info[0] >= 3:
    from http.client import HTTPSConnection
else:
    from httplib import HTTPSConnection

import extra_utils
logger=extra_utils.makelogger('preprocess')

import mwapy
from mwapy import metadata
import find_calibrator

################################################################################
# from obsdownload.py
# 2016-08-17
################################################################################

KEY_FILE = os.path.join(os.getenv('HOME'),'obskey.key')
CERT_FILE = os.path.join(os.getenv('HOME'),'obscrt.crt')

######################################################################
class HTTPSClientAuthConnection(HTTPSConnection):
    """Class to make a HTTPS connection, with support for full client-based SSL Authentication"""

    def __init__(self, host, key_file, cert_file, timeout=None):
        HTTPSConnection.__init__(self, host, key_file=key_file, cert_file=cert_file)
        self.key_file = key_file
        self.cert_file = cert_file
        self.ca_file = None
        self.timeout = timeout

    def connect(self):
        """ Connect to a host on a given (SSL) port.
            If ca_file is pointing somewhere, use it to check Server Certificate.

            Redefined/copied and extended from httplib.py:1105 (Python 2.6.x).
            This is needed to pass cert_reqs=ssl.CERT_REQUIRED as parameter to ssl.wrap_socket(),
            which forces SSL to check server certificate against our client certificate.
        """
        sock = socket.create_connection((self.host, self.port), self.timeout)
        if self._tunnel_host:
            self.sock = sock
            self._tunnel()
        # If there's no CA File, don't force Server Certificate Check
        if self.ca_file:
            self.sock = ssl.wrap_socket(sock, self.key_file, self.cert_file,
                                        ca_certs=self.ca_file, cert_reqs=ssl.CERT_REQUIRED)
        else:
            self.sock = ssl.wrap_socket(sock, self.key_file, self.cert_file, cert_reqs=ssl.CERT_NONE)



class FileStatus():
    def __init__(self, numfiles):
        self.status = {}
        self.lock = threading.RLock()
        self.errors = []
        self.filesComplete = 0
        self.totalfiles = numfiles;

    def file_error(self, err):
        with self.lock:
            self.errors.append(err)
            print(err)

    def file_starting(self, filename):
        with self.lock:
            print('%s [INFO] Downloading %s' % (time.strftime('%c'), filename))

    def file_complete(self, filename):
        with self.lock:
            self.filesComplete = self.filesComplete + 1
            print('%s [INFO] %s complete [%d of %d]' % (time.strftime('%c'), filename,
                                                        self.filesComplete, self.totalfiles))

# def delete_key_file():
#     try:
#         os.remove(KEY_FILE)
#     except:
#         pass
#     try:
#         os.remove(CERT_FILE)
#     except:
#         pass

# def create_key_file():
#     tmpkey = str(uuid.uuid1())
#     tmpcert = str(uuid.uuid1())

#     with open(tmpkey, 'wb') as out:
#         out.write(KEY)
#         out.flush()

#     with open(tmpcert, 'wb') as out:
#         out.write(CERT)
#         out.flush()

#     os.rename(tmpkey, KEY_FILE)
#     os.rename(tmpcert, CERT_FILE)
#     os.chmod(KEY_FILE, stat.S_IRUSR)
#     os.chmod(CERT_FILE, stat.S_IRUSR)

######################################################################
def query_observation(obs=None, host=None, flagonly=False, metaonly=False, uvfits=False, chstart=None, chcount=False):

    if not os.path.exists(KEY_FILE):
        logger.error('Unable to find HTTPS key file %s' % KEY_FILE)
        return None
    if not os.path.exists(CERT_FILE):
        logger.error('Unable to find HTTPS certificate file %s' % CERT_FILE)
        return None


    response = None
    try:
        conn = HTTPSClientAuthConnection(host, key_file=KEY_FILE, cert_file=CERT_FILE)
        conn.request('GET', '/QUERY?query=files_like&like=%s%%&format=json' % obs)
        response = conn.getresponse()
        if response.status != 200:
            raise Exception('HTTP error: %s Reason: %s' \
                            % (response.status, response.reason))

        length = int(response.getheader('content-length', 0))
        buffsize = 4096
        readin = 0
        resultbuffer = []
        while readin < length:
            left = length - readin
            buff = response.read(buffsize if left >= buffsize else left)
            if not buff:
                raise Exception('socket read error')
            resultbuffer.append(buff)
            readin += len(buff)
        filemap = {}
        json_files = json.loads(b''.join(resultbuffer).decode('utf8'))
        files = json_files.get('files_like', [])

        for f in files:
            add = False

            filename = f['col3']
            filesize = f['col6']

            if flagonly or uvfits or metaonly:
                if 'flags.zip' in filename and flagonly:
                    add = True
                if '.uvfits' in filename and uvfits:
                    add = True
                if '_ppds' in filename and metaonly:
                    add = True
            else:
                if chstart is not None:
                    # it must be a valid MWA visibility fits file which contains
                    # the coarse channel. If valid extract the coarse channel
                    # and check with the range requested.
                    try:
                        num = int(filename.split('_')[2].split('gpubox')[1])
                        if (num >= chstart) and (num < chstart+chcount):
                            add = True
                    except:
                        pass
                else:
                    add = True

            if add is True:
                filemap[filename] = filesize

        return filemap

    finally:
        if response:
            response.close()


######################################################################
def check_complete(filename, size, dir):
    path = dir + filename
    if os.path.isfile(path) is True:
        filesize = os.stat(path).st_size
        if filesize == int(size):
            return True

    return False


######################################################################
def download_worker(url, host, size, filename, sem, out, stat, bufsize):

    resp = None

    try:
        stat.file_starting(filename)

        conn = HTTPSClientAuthConnection(host, key_file = KEY_FILE,
                                        cert_file = CERT_FILE, timeout = 5000)
        conn.request('GET', '/RETRIEVE?file_id=%s' % filename)
        response = conn.getresponse()
        if response.status != 200:
            raise Exception('HTTP error: %s Reason: %s' \
                            % (response.status, response.reason))

        length = int(response.getheader('content-length', 0))

        with open(out + filename, 'wb') as f:
            readin = 0
            while readin < length:
                try:
                    response.fp._sock.settimeout(120)
                except:
                    pass
                left = length - readin
                buff = response.read(bufsize if left >= bufsize else left)
                if not buff:
                    raise Exception('socket read error')
                f.write(buff)
                readin += len(buff)

        stat.file_complete(filename)

    except Exception as exp:
        stat.file_error('%s [ERROR] %s: %s' % (time.strftime('%c'), filename,
                                                    str(exp)))
    finally:
        if resp:
            resp.close()

        sem.release()


######################################################################
def run_obsdownload(obs,ngashost='fe1.pawsey.org.au:8777',
                    numdownload=4):
    try:
        # get system TCP buffer size
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        bufsize = s.getsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF)
        s.close()

        if numdownload <= 0 or numdownload > 6:
            logger.error('Number of simultaneous downloads must be > 0 and <= 6')
            return None

        logger.info('Finding observation %s' % (obs,))
            
        fileresult = query_observation(obs, ngashost)
        if len(fileresult) <= 0:
            logger.error('No file(s) found for observation %s' % (obs,))
            return None

        logger.info('Found %s files' % (str(len(fileresult)),))

        outputdir=os.path.join(os.curdir,
                               str(obs),'')
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)

        urls = []
        stat = FileStatus(len(fileresult))

        for key, value in fileresult.items():
            url = 'https://%s/RETRIEVE?file_id=%s' % (ngashost, key)
            if check_complete(key, int(value), outputdir) is False:
                urls.append((url, value, key))
            else:
                stat.file_complete(key)

        threads = []
        sem = threading.BoundedSemaphore(value = numdownload)
        for u in urls:
            while True:
                if sem.acquire(blocking = False):
                    t = threading.Thread(target = download_worker,
                                         args = (u[0], ngashost, u[1], u[2],
                                                 sem, outputdir, stat, int(bufsize)))
                    t.setDaemon(True)
                    t.start()
                    threads.append(t)
                    break
                else:
                    time.sleep(1)

        # wait for all the threads to finish
        while len(threads) > 0:
            for t in list(threads):
                t.join(1)
                if not t.isAlive():
                    threads.remove(t)

        logger.info('File Transfer Complete.')

        # check if we have errors
        if stat.errors:
            logger.warning('File Transfer Error Summary:')
            for errs in stat.errors:
                logger.warning(errs)
                raise Exception()
        else:
            logger.info('File Transfer Success.')
        return fileresult

    except KeyboardInterrupt as k:
        logger.error('Interrupted; terminating download')
        return None
    
    except Exception as e:
        logger.error(e)
        return None



################################################################################
# from obsdownload.py
# 2015-02-05
################################################################################
# username = 'ngas'
# password = base64.decodestring('bmdhcw==')
    
# ######################################################################
# class PrintStatus():
    
#     def __init__(self, numfiles):
#         self.status = {}
#         self.lock = threading.RLock()
#         self.currentbytes = 0
#         self.totalbytes = 0
#         self.runtime = 0
#         self.errors = []
#         self.files = 0
#         self.filesComplete = 0
#         self.totalfiles = numfiles;

    
#     def fileError(self, err):
#         with self.lock:
#             self.errors.append(err)
#             logger.error(err)
        
#     def fileStarting(self, filename):
#         with self.lock:
#             logger.info('Downloading %s' % (filename))
            
#     def fileComplete(self, filename):
#         with self.lock:
#             self.filesComplete = self.filesComplete + 1
#             logger.info('%s complete [%d of %d]' % (filename,
#                                                     self.filesComplete,
#                                                     self.totalfiles))
            
# ######################################################################
# def queryObs(obs, host, flagonly):
#     #http://fe1.pawsey.ivec.org:7777/QUERY?query=files_like&like=1061316296%&format=json
    
#     base64string = base64.encodestring('%s:%s' % (username, password)).replace('\n', '')
    
#     url = 'http://%s/QUERY?' % host
#     data = "query=files_like&like=" + str(obs) + "%&format=json"
    
#     response = None
#     try:
#         # Send HTTP POST request
#         request = urllib2.Request(url + data)
#         request.add_header("Authorization", "Basic %s" % base64string)   
        
#         response = urllib2.urlopen(request)
        
#         resultbuffer = ''
#         while True:
#             result = response.read(2048)
#             if result:
#                 resultbuffer = resultbuffer + result
#             else:
#                 break

#         filemap = {}
#         files = json.loads(resultbuffer)['files_like']
                
#         for f in files:
#             # retrieve flag files only
#             if flagonly:
#                 if 'flags.zip' in f['col3']:
#                     filemap[f['col3']] = f['col6']
#             # else add all the files: having a map removes diplicates
#             else:
#                 filemap[f['col3']] = f['col6']
        
#         return filemap
    
#     finally:
#         if response:
#             response.close()

 
# ######################################################################
# def checkFile(filename, size, dir):
#     path = dir + filename
        
#     # check the file exists
#     if os.path.isfile(path) is True:
#         #check the filesize matches
#         filesize = os.stat(path).st_size

#         if filesize == int(size):
#             return True
        
#     return False

# ######################################################################
# def worker(url, size, filename, s, out, stat, bufsize):
#     u = None
#     f = None
    
#     try:
#         stat.fileStarting(filename)

#         # open file URL
#         request = urllib2.Request(url)
#         base64string = base64.encodestring('%s:%s' % (username, password)).replace('\n', '')
#         request.add_header("Authorization", "Basic %s" % base64string)   
        
#         u = urllib2.urlopen(request, timeout=1400)        
#         u.fp.bufsize = bufsize
        
#         # get file size
#         meta = u.info()
#         file_size = int(meta.getheaders("Content-Length")[0])
        
#         # open file for writing
#         f = open(out + filename, 'wb')
        
#         current = 0
#         file_size_dl = 0
#         block_sz = bufsize
        
#         while file_size_dl < file_size:
#             buffer = u.read(block_sz)
#             if buffer:
#                 f.write(buffer)

#                 current = len(buffer)
#                 file_size_dl += current
                
#             else:
#                 if file_size_dl != file_size:
#                     raise Exception("size mismatch %s %s" % str(file_size), str(file_size_dl))
                
#                 break

#         stat.fileComplete(filename)
        
#     except urllib2.HTTPError as e:
#         stat.fileError("%s [ERROR] %s %s" % (time.strftime("%c"), filename, str(e.read()) ))
    
#     except urllib2.URLError as urlerror:
#         if hasattr(urlerror, 'reason'):
#             stat.fileError("%s [ERROR] %s %s" % (time.strftime("%c"), filename, str(urlerror.reason) ))
#         else:
#             stat.fileError("%s [ERROR] %s %s" % (time.strftime("%c"), filename, str(urlerror) ))
    
#     except Exception as exp:
#         stat.fileError("%s [ERROR] %s %s" % (time.strftime("%c"), filename, str(exp) ))
        
#     finally:
#         if u != None:
#             u.close()
            
#         if f != None:
#             f.flush()
#             f.close()
            
#         s.release()



# ######################################################################
# def old.run_obsdownload(obs,ngashost='fe1.pawsey.ivec.org:7777',
#                     numdownload=4):
#     if numdownload <= 0 or numdownload>6:
#         logger.error('Number of downloads must be between 0 and 6')
#         return None

#     try:
#         # get system TCP buffer size
#         s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
#         bufsize = s.getsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF)
#         s.close()

#         logger.info('Finding observation %s' % obs)
        
#         fileresult = queryObs(obs, ngashost, False)
#         #fileresult = queryObs(obs, ngashost, True)
        
#         if len(fileresult) <= 0:
#             logger.error('No files found for observation %s' % obs)
#             return None
        
#         logger.info('Found %s files' % len(fileresult))
#         outputdir=os.path.join(os.curdir,
#                                str(obs),
#                                '')
    
#         if not os.path.exists(outputdir):
#             os.makedirs(outputdir)
        
#         stat = PrintStatus(len(fileresult))
#         urls = []
        
#         for key, value in fileresult.iteritems():
#             url = "http://" + ngashost + "/RETRIEVE?file_id=" + key
#             if checkFile(key, int(value), outputdir) is False:
#                 urls.append((url, value, key))
#                 # stat.update(url, r[2], 0, 0, 0, 0)    
#             else:
#                 stat.fileComplete(key)

#         s = threading.BoundedSemaphore(value=numdownload)        
#         for u in urls:
#             while(1):
#                 if s.acquire(blocking=False):
#                     t = threading.Thread(target=worker, args=(u[0],u[1],u[2],s,outputdir,stat,int(bufsize)))
#                     t.setDaemon(True)
#                     t.start()
#                     break
#                 else:
#                     time.sleep(1)

#         while (1):
#             main_thread = threading.currentThread()
#             # if the main thread and print thread are left then we must be done; else wait join
#             if len(threading.enumerate()) == 1:
#                 break;
            
#             for t in threading.enumerate():
#                 # if t is main_thread or t is status_thread:
#                 if t is main_thread:
#                     continue
                
#                 t.join(1)
#                 if t.isAlive():
#                     continue

#         logger.info('File Transfer Complete')
        
#         # check if we have errors
#         if len(stat.errors) > 0:
#             logger.warning('File Transfer Error Summary:')
#             for i in stat.errors:
#                 logger.warning(i)
                
#             #raise Exception()
#         else:
#             logger.info('File Transfer Success')

#         return fileresult
#     except KeyboardInterrupt as k:
#         logger.error('Interrupted; terminating download')
#         return None
    
#     except Exception as e:
#         logger.error(e)
#         return None

##################################################
def tokenize(input):
    """
    web queries use different wildcards than normal people
    turn ? to _
    turn % to *
    """
    
    if input is not None:
        input=input.replace('_','\_')
        return input.replace('?','_').replace('*','%')
    else:
        return input

##################################################
def parse_time(input):
    """
    parse input time (GPStime or datetime string)
    """
    try:
        return int(input)
    except ValueError:
        # try to parse it as a date
        try:
            t=Time(input, scale='utc')
        except ValueError,e:
            try:
                t=Time(datetime.datetime.strptime(input, '%Y%m%d%H%M%S'), scale='utc')
            except ValueError,e:                                
                logger.error('Unable to parse input time %s: %s' % (input,
                                                                    e))
                return None
        return int(t.gps)

######################################################################
def get_size(start_path = '.'):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)
    return total_size

######################################################################
def get_free(path=None):
    if path is None:
        path=os.path.abspath(os.cwd())

    if ':' in path:
        # assume a remote location
        # need ssh
        host,path=path.split(':')
        command='ssh %s "df -k %s"' % (host,path)
        try:
            result=subprocess.check_output(command, shell=True)
            if len(result.split('\n'))==3:
                return int(result.split('\n')[1].split()[3])*1024
            elif len(result.split('\n'))==4:
                print result.split('\n')[2]
                # an extra line-break
                return int(result.split('\n')[2].split()[2])*1024
        except Exception, e:
            logger.error('Unable to check free space on %s:%s\n%s' % (host,path,e))
            return 0
    else:
        # it is local
        try:
            stat=os.statvfs(path)
        except OSError:
            return 0
        return stat.f_bavail*stat.f_bsize


######################################################################
class Observation():

    ##############################
    def __init__(self, obsid, basedir='./', clobber=False, 
                 cleanfits=False,
                 cleanms=False,
                 cleanflag=False,
                 checkfree=True,
                 checksize=True,
                 copycommand='rsync -aruvP'):
        self.obsid=obsid
        self.basedir=basedir
        self.downloadedfiles=None
        self.metafits=None
        self.cottercommand=None
        self.msfile=None
        self.clobber=clobber
        self.copycommand=copycommand
        self.cleanfits=cleanfits
        self.cleanms=cleanms
        self.cleanflag=cleanflag
        self.destination=None
        self.rawsize=None
        self.fitssize=None
        self.mssize=None
        self.checkfree=checkfree
        self.checksize=checksize

        self.outputdir=os.path.join(self.basedir,
                                    str(self.obsid),
                                    '')
        os.chdir(basedir)
        logger.info('##################################################')
        logger.info('%d' % self.obsid)
        self.observation=metadata.MWA_Observation(self.obsid)
        logger.debug('\n' + str(self.observation))
        

        obsinfo=metadata.fetch_obsinfo(self.obsid)
        self.rawsize=sum([obsinfo['files'][f]['size'] for f in obsinfo['files'].keys()])
        self.fitssize=sum([obsinfo['files'][f]['size']*(obsinfo['files'][f]['filetype']==8) for f in obsinfo['files'].keys()])
        logger.info('Raw data size is %.1f MB' % (self.rawsize/1024./1024.))



    ##############################
    def __del__(self):
        if self.cleanfits:
            logger.warning('Deleting downloaded FITS files %s' % ','.join(self.downloadedfiles))
            for f in self.downloadedfiles:
                os.remove(os.path.join(self.outputdir,f))
        if self.cleanms:
            if self.destination is None:
                logger.warning('Requesting delete of MS file %s, but it has not been copied elsewhere' % self.msfile)
            else:
                logger.warning('Deleting MS file %s' % self.msfile)
                shutil.rmtree(os.path.join(self.outputdir,self.msfile))
        if self.cleanflag:
            logger.warning('Deleting downloaded flag files %s/*.mwaf' % self.outputdir)
            for f in glob.glob('%s/*.mwaf' % self.outputdir):
                os.remove(f)


    ##############################
    def download(self, numdownload=4):
        if self.checkfree and (self.rawsize is not None and self.rawsize > 0 and self.rawsize > get_free(self.basedir)):
            logger.error('Data download will require %.1f MB, but only %.1f MB free space on %s' % (self.rawsize/1024./1024,
                                                                                                    get_free(self.basedir)/1024./1024,
                                                                                                    self.basedir))
            return None
        self.downloadedfiles=run_obsdownload(self.obsid, numdownload=numdownload)
        #downloadedfiles2=run_obsdownload(self.obsid, numdownload=numdownload)
        #for k in downloadedfiles2.keys():
        #    self.downloadedfiles[k]=downloadedfiles2[k]
        
        os.system(s)

        return self.downloadedfiles
    ##############################
    def makemetafits(self, min_bad_dipoles=5):        
        logger.info('Making metafits with min_bad_dipoles=%d' % min_bad_dipoles)
        m=metadata.instrument_configuration(self.obsid, min_bad_dipoles=min_bad_dipoles)
        h=m.make_metafits()
        self.metafits=os.path.join(self.outputdir,'%s.metafits' % self.obsid)
        if os.path.exists(self.metafits):
                os.remove(self.metafits)
        try:
            h.writeto(self.metafits)
            logger.info('Metafits written to %s' % (self.metafits))
            return True
        except Exception, e:
            logger.error('Unable to write metafits file %s:\n%s' % (self.metafits,e))
            return None
    ##############################
    def cotter(self, cottercpus=4,
               timeres=4,
               freqres=40,
               phasecenter=None,
               allowmissing=True):


        if timeres < self.observation.inttime:
            logger.warning('Desired time resolution %.1f s is < native time resolution %.1f s' % (
                timeres,self.observation.inttime))
            timeres=self.observation.inttime
            logger.warning('Changing to %.1f s' % timeres)
        if freqres < self.observation.fine_channel:
            logger.error('Desired frequency resolution %d kHz is < native frequency resolution %d kHz' % (
                freqres,self.observation.fine_channel))
            return None
        compressionfactor_time=int(timeres/self.observation.inttime)
        compressionfactor_freq=int(freqres/self.observation.fine_channel)
        logger.info('Expect compression of %d in time and %d in frequency' % (compressionfactor_time,
                                                                              compressionfactor_freq))
        self.mssize=self.fitssize/compressionfactor_time/compressionfactor_freq
        if self.checkfree and (self.mssize is not None and self.mssize > 0 and self.mssize > get_free(self.basedir)):
            logger.error('Cotter will require %.1f MB, but only %.1f MB free space on %s' % (self.mssize/1024./1024,
                                                                                             get_free(self.basedir)/1024./1024,
                                                                                             self.basedir))
            return None

        os.chdir(self.outputdir)
        self.useflagfiles=False
        # check to see if a .zip file for flags was included
        for f in self.downloadedfiles:
            if '.zip' in f:
                logger.info('Unzipping flags from %s' % f)
                subprocess.Popen(['unzip','-j',f],stderr=subprocess.PIPE,
                                 stdout=subprocess.PIPE)
                self.useflagfiles=True
        # cotter gets angry if only one timeslice
        newfiles=[]
        for f in self.downloadedfiles:
            if re.match('\d+_\d+_gpubox\d\d\.fits',f) is not None:
                logger.warning('Format for %s is wrong; linking %s_00.fits to %s...' % (f,os.path.splitext(f)[0],f))
                try:
                    os.symlink(f,'%s_00.fits' % os.path.splitext(f)[0])
                except OSError:
                    logger.warning('%s_00.fits exists' % os.path.splitext(f)[0])
                newfiles.append('%s_00.fits' % os.path.splitext(f)[0])
        for f in newfiles:
            self.downloadedfiles[f]=0
                
        self.cottercommand=['cotter','-j',str(cottercpus),
                            '-m',self.metafits,
                            '-timeres',str(timeres),
                            '-freqres',str(freqres),
                            '-o','%s.ms' % self.obsid]
        if self.useflagfiles:
            self.cottercommand+=['-flagfiles','%s_%%%%.mwaf' % self.obsid]
        #self.cottercommand+=['-norfi']
        if phasecenter is not None:
            self.cottercommand+=['-centre',phasecenter.to_string('hmsdms').split()[0],
                                 phasecenter.to_string('hmsdms').split()[1]]
                                 
        if allowmissing:
            self.cottercommand+=['-allowmissing']
        self.cottercommand+=['*gpubox*_*.fits']        
        self.msfile='%s.ms' % self.obsid
        if os.path.exists(self.msfile):
            if not self.clobber:
                logger.warning('MSfile %s exists; will not clobber' % self.msfile)
                return True
            else:
                logger.warning('MSfile %s exists; will clobber' % self.msfile)
                shutil.rmtree(self.msfile)
                
        logger.info('Will run:\n\t%s' % ' '.join(self.cottercommand))
        p=subprocess.Popen(' '.join(self.cottercommand),
                           stderr=subprocess.PIPE,
                           stdout=subprocess.PIPE,
                           shell=True,
                           close_fds=False)
        waserror=False
        while True:
            p.stdout.flush()
            p.stderr.flush()
            for l in p.stdout.readlines():
                logger.debug(l.rstrip())
            for l in p.stderr.readlines():
                logger.error(l.rstrip())
                waserror=True
            returncode=p.poll()
            if returncode is not None:
                break
            time.sleep(1)
        if waserror:
            logger.error('Cotter had an error')
            return None

        if not os.path.exists(self.msfile):
            logger.error('MSfile %s was not created' % self.msfile)
            return None

        newsize=get_size(self.msfile)
        if math.fabs(newsize-self.mssize)/float(self.mssize) > 0.75 and self.checksize:
            logger.error('Expected a size of %d bytes for MSfile %s, but actual file has %d bytes' % (self.mssize,
                                                                                                      self.msfile,
                                                                                                      newsize))
            return None
        self.mssize=get_size(self.msfile)
        logger.info('MSfile %s was created with size %d bytes' % (self.msfile,self.mssize))
        return self.mssize
    
    ##############################
    def copy(self, destination):
        if destination is not None:
            if self.checkfree and (self.mssize is not None and self.mssize > 0 and self.mssize > get_free(destination)):
                logger.error('Copy will require %.1f MB, but only %.1f MB free space on %s' % (self.mssize/1024./1024,
                                                                                               get_free(destination)/1024./1024,
                                                                                               destination))
                return None

            copycommand=[self.copycommand,
                         self.msfile,
                         destination + '/']
            logger.info('Will run:\n\t%s' % ' '.join(copycommand))
            p=subprocess.Popen(' '.join(copycommand),
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
            if returncode > 0:
                logger.error('Error copying %s to %s' % (self.msfile, destination))
                return None

            logger.info('%s copied to %s' % (self.msfile,destination))

            copycommand=[self.copycommand,
                         self.metafits,
                         destination + '/']
            logger.info('Will run:\n\t%s' % ' '.join(copycommand))
            p=subprocess.Popen(' '.join(copycommand),
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
            if returncode > 0:
                logger.error('Error copying %s to %s' % (self.metafits, destination))
                return None

            logger.info('%s copied to %s' % (self.metafits,destination))

            self.destination=destination
            return True


######################################################################
def main():

    observation_num=None

    usage="Usage: %prog [options] <obsids>\n"
    parser = OptionParser(usage=usage)
    parser.add_option('--start',dest='starttime',default=None,type='str',
                      help='Starting time for processing (GPSseconds or ISO UT date)')
    parser.add_option('--stop',dest='stoptime',default=None,type='str',
                      help='Final time for processing (GPSseconds or ISO UT date)')
    parser.add_option('--proj',dest='project',default=None,
                      help='Project ID')
    parser.add_option('--destination',dest='destination',default=None,
                      help='Final destination to (remote) copy the MS file')
    parser.add_option('--summaryname',dest='summaryname',default=None,
                      help='Name for output summary table')
    parser.add_option('--downloads',dest='downloads',default=4,type='int',
                      help='Number of simultaneous NGAS downloads [default=%default]')
    parser.add_option('--nomissing',dest='allowmissing',default=True,
                      action='store_false',
                      help='Do not allow missing GPU box files')
    parser.add_option('--timeres',dest='timeres',default=4,type='float',
                      help='Output time resolution (s) [default=%default]')
    parser.add_option('--freqres',dest='freqres',default=40,type='int',
                      help='Output frequency resolution (kHz) [default=%default]')
    parser.add_option('--caltimeres',dest='caltimeres',default=0,type='float',
                      help='Output time resolution for calibrators (s) [default=<timeres>]')
    parser.add_option('--calfreqres',dest='calfreqres',default=0,type='int',
                      help='Output frequency resolution (kHz) [default=<freqres>]')
    parser.add_option('--ra',dest='ra',default=None,
                      help='Output RA phase center (deg) [default=metafits]')
    parser.add_option('--dec',dest='dec',default=None,
                      help='Output Dec phase center (deg) [default=metafits]')
    parser.add_option('--baddipoles',dest='minbaddipoles',default=2,type='int',
                      help='Number of bad dipoles above which tile is flagged [default=%default]')
    parser.add_option('--cal',dest='cal',default=False,
                      action='store_true',
                      help='Always include a calibrator observation?')
    parser.add_option('--cpus',dest='cottercpus',default=4,type='int',
                      help='Number of CPUs for cotter [default %default]')
    parser.add_option('--clobber',dest='clobber',default=False,action='store_true',
                      help='Clobber existing MS file? [default=False]')
    parser.add_option('--cleanall',dest='cleanall',default=False,action='store_true',
                      help='Delete downloaded FITS files, flag files, and original MS file when done? [default=False]')
    parser.add_option('--cleanfits',dest='cleanfits',default=False,action='store_true',
                      help='Delete downloaded FITS files when done? [default=False]')
    parser.add_option('--cleanms',dest='cleanms',default=False,action='store_true',
                      help='Delete original MS file when done? [default=False]')
    parser.add_option('--copycommand',dest='copycommand',default='rsync -aruvP',
                      help='Command for (remote) copying [default=%default]')
    parser.add_option('--summaryformat',dest='summaryformat',default='ascii.commented_header',
                      help='Format for output summary table (from astropy.table) [default=%default]')
    parser.add_option('--email',dest='notifyemail', default=None,
                      help='Notification email [default=no email]')
    parser.add_option('--ignorespace',dest='ignorespace',default=False,
                      action='store_true',
                      help='Ignore free-space for file download, processing, and copying?')
    parser.add_option('--ignoresize',dest='ignoresize',default=False,
                      action='store_true',
                      help='Ignore size of the .ms file?')
    parser.add_option("-v", "--verbose", dest="loudness", default=0, action="count",
                      help="Each -v option produces more informational/debugging output")
    parser.add_option("-q", "--quiet", dest="quietness", default=0, action="count",
                      help="Each -q option produces less error/warning/informational output")


    (options, args) = parser.parse_args()
    if options.cleanall:
        options.cleanfits=True
        options.cleanms=True
        options.cleanflag=True

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
    logger.info('%s starting at %s UT on host %s with user %s' % (sys.argv[0],
                                                                  datetime.datetime.now(),
                                                                  socket.gethostname(),
                                                                  os.environ['USER']))
    logger.info('**************************************************')

    if options.ra is not None and options.dec is not None:
        phasecenter=SkyCoord(options.ra, options.dec,frame='icrs',unit='deg')
    else:
        phasecenter=None

    if len(args)==0:
        if options.starttime is None:
            logger.error('Must supply starttime')
            eh.exit(1)
        if options.stoptime is None:
            logger.error('Must supply stoptime')
            eh.exit(1)
    
    logger.debug('Using mwapy version %s' % mwapy.__version__)
    result=subprocess.Popen(['cotter','-version'],
                            stdout=subprocess.PIPE).communicate()[0].strip()
    cotterversion=result.split('\n')[0].split('version')[1].strip()
    AOFlaggerversion=result.split('\n')[1].split('AOFlagger')[1].strip()
    logger.debug('Using cotter version %s' % cotterversion)
    logger.debug('Using AOFlagger version %s' % AOFlaggerversion)
    results=[]
    havecalibrator={}

    if not (options.starttime is None and options.stoptime is None):
        GPSstart=parse_time(options.starttime)
        GPSstop=parse_time(options.stoptime)
        if options.project is not None:
            logger.info('Will preprocess observations from %s to %s with project=%s' % (GPSstart,
                                                                                        GPSstop,
                                                                                        options.project))
        else:
            logger.info('Will preprocess observations from %s to %s' % (GPSstart,
                                                                        GPSstop))


        if options.project is None:
            results=metadata.fetch_observations(mintime=GPSstart-1,
                                                maxtime=GPSstop+1)
        else:
            results=metadata.fetch_observations(mintime=GPSstart-1,
                                                maxtime=GPSstop+1,
                                                projectid=tokenize(options.project))
    if len(args)>0:
        # add in some additional obsids
        for arg in args:
            results+=metadata.fetch_observations(mintime=int(arg)-1,
                                                 maxtime=int(arg)+1)

    if results is None or len(results)==0:
        logger.error('No observations found')
        eh.exit(1)
    logger.info(metadata.MWA_Observation_Summary.string_header())
    # check if there is a calibrator present
    observations=[]
    for item in results:            
        o=metadata.MWA_Observation_Summary(item)
        logger.info(str(o))
        observations.append(metadata.MWA_Observation(item[0]))
        if havecalibrator.has_key(observations[-1].center_channel):
            havecalibrator[observations[-1].center_channel] = havecalibrator[observations[-1].center_channel] or observations[-1].calibration
        else:
            havecalibrator[observations[-1].center_channel] = observations[-1].calibration

    if not all(havecalibrator.values()):
        for k in havecalibrator.keys():
            if not havecalibrator[k]:
                logger.warning('No calibrator observation found for center channel %d' % k)

        if options.cal:
            cals=[]
            for i in xrange(len(results)):
                channel=observations[i].center_channel
                if not havecalibrator[channel]:
                    cal=find_calibrator.find_calibrator(observations[i].observation_number,
                                                        matchproject=False,
                                                        priority='time',
                                                        all=False)
                    if cal is None:
                        logger.error('Unable to find an appropriate calibrator scan')
                        eh.exit(1)
                    logger.info('Found calibrator observation %d for center channel %d' % (cal[0],channel))
                    cals+=metadata.fetch_observations(mintime=cal[0]-1,
                                                      maxtime=cal[0]+1)
                    logger.info(str(metadata.MWA_Observation_Summary(cals[-1])))
                    havecalibrator[channel]=True
            results+=cals
    data=None

    logger.info('Processing observations...\n\n')
    basedir=os.path.abspath(os.curdir)
    for obs in results:
        timeres,freqres=None,None
        processstart=datetime.datetime.now()
        observation=Observation(obs[0], basedir=basedir,
                                copycommand=options.copycommand,
                                clobber=options.clobber,
                                cleanms=options.cleanms,
                                cleanfits=options.cleanfits,
                                checkfree=not options.ignorespace,
                                checksize=not options.ignoresize)
        downloadedfiles=observation.download(numdownload=options.downloads)
        if downloadedfiles is None:
            break

        if observation.makemetafits(min_bad_dipoles=options.minbaddipoles) is None:
            break
            
        if observation.observation.calibration is not None and observation.observation.calibration:
            if options.caltimeres==0:
                timeres=options.timeres
            else:
                timeres=options.caltimeres
            if options.calfreqres==0:
                freqres=options.freqres
            else:
                freqres=options.calfreqres
        else:
            timeres,freqres=options.timeres,options.freqres
        msfilesize=observation.cotter(cottercpus=options.cottercpus,
                                      timeres=timeres,
                                      freqres=freqres,
                                      phasecenter=phasecenter,
                                      allowmissing=options.allowmissing)
        if msfilesize is None:
            break

        if options.destination is not None and observation.copy(options.destination) is None:
            break

        row=collections.OrderedDict({'date': processstart.strftime('%Y-%m-%dT%H:%M:%S'),
                                     'obsid': obs[0],
                                     'user': os.environ['USER'],
                                     'host': socket.gethostname(),
                                     'cotter': cotterversion.replace(' ','_'),
                                     'mwapy': mwapy.__version__,
                                     'rawfiles': len(observation.downloadedfiles),
                                     'msfilesize': msfilesize})
        if options.destination is not None:
            row['msfile']=os.path.join(options.destination,observation.msfile)
        else:
            row['msfile']=os.path.join(observation.outputdir,observation.msfile)
        if data is None:
            data=Table([row])
        else:
            data.add_row(row)
    os.chdir(basedir)
    if data is not None and options.summaryname is not None:
        data=data[('date','obsid','user','host',
                   'cotter','mwapy','rawfiles','msfilesize','msfile')]
        try:
            data.write(options.summaryname,format=options.summaryformat)
            logger.info('Summary table written to %s' % options.summaryname)
        except Exception, e:
            logger.error('Unable to write summary table %s with format %s:\n%s' % (options.summaryname,
                                                                                   options.summaryformat,
                                                                                   e))

    eh.exit(0)

        
                                                                   
######################################################################

if __name__=="__main__":

    main()
