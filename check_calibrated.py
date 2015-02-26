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
logger.setLevel(logging.WARNING)

try:
    import drivecasa
    _CASA=True
except ImportError:
    _CASA=False

casapy='/usr/local/casapy'

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
             
    logger.info('Will run in casa:\n\t%s' % '\n\t'.join(command))
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

######################################################################                                                                                                        
def main():

    usage="Usage: %prog [options] <msfile>\n"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", dest="verbose",default=False,action='store_true',
                      help='More verbose output?')

    (options, args) = parser.parse_args()

    if options.verbose:
        logger.setLevel(logging.INFO)

    for msfile in args:
        result=check_calibrated(msfile)
        if result is not None:
            if result:
                print '%s has CORRECTED_DATA column' % msfile
            else:
                print '%s has no CORRECTED_DATA column' % msfile
        else:
            print 'Unable to determine if %s has CORRECTED_DATA column' % msfile
    

######################################################################                                                                                                        
if __name__=="__main__":

    main()

