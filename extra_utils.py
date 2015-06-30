import smtplib
from email.mime.text import MIMEText
import os,datetime,socket,sys
import logging,logging.handlers

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

        elif record.levelno >= logging.ERROR:
            self._fmt = MyFormatter.err_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._fmt = format_orig

        return result

######################################################################
def makelogger(name):

    fmt = MyFormatter()
    # set up logging to file - see previous section for more details
    logging.basicConfig(level=logging.DEBUG,
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filename='/dev/null')


    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.WARNING)
    console.setFormatter(fmt)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)
    
    filehandler = logging.handlers.RotatingFileHandler(name + '.log',
                                                       backupCount=10)
    filehandler.setLevel(logging.DEBUG)
    filehandler.setFormatter(fmt)
    filehandler.doRollover()
    logging.getLogger('').addHandler(filehandler)

    return logging.getLogger(name)



######################################################################
class ExitHandler():
    """
    a handler for exiting from a script
    will do various cleanup and notification items on exit
    including sending an email if desired
    """
    def __init__(self,
                 name,
                 email=None):
        self.name=name
        self.email=email

    def sendemail(self):
        if self.email is None:
            return
        fromaddr='%s@%s' % (self.name,socket.gethostname())
        toaddr=self.email
        # check if we are in a screen
        if os.environ.has_key('STY'):
            subject='%s (PID=%s, STY=%s) finished at %s with exit code %d' % (self.name,
                                                                              os.getpid(),
                                                                              os.environ['STY'],
                                                                              datetime.datetime.now(),
                                                                              self.code)
        else:
            subject='%s (PID=%s) finished at %s with exit code %d' % (self.name,
                                                                      os.getpid(),
                                                                      datetime.datetime.now(),
                                                                      self.code)
            
        logfile=logging.getLogger('').handlers[2].baseFilename
        f=open(logfile)
        lines=f.readlines()
        f.close()
        # filter out some extra aegean output if needed
        goodlines=[]
        for line in lines:
            if not 'DEBUG:root' in line:
                goodlines.append(line)
        msg=MIMEText(''.join(goodlines))
        msg['Subject']=subject
        msg['To']=toaddr
        msg['From']=fromaddr
        server = smtplib.SMTP('localhost')
        try:
            server.sendmail(fromaddr, [toaddr], msg.as_string())
            logging.getLogger('').info('Sent email to %s' % toaddr)
        except Exception,e:
            logging.getLogger('').error('Cannot send notification email to %s:\n\t%s' % (toaddr,e))
            msg=MIMEText('[message text too large]')
            msg['Subject']=subject
            msg['To']=toaddr
            msg['From']=fromaddr
            try:
                server.sendmail(fromaddr, [toaddr], msg.as_string())
                logging.getLogger('').info('Sent truncated email to %s' % toaddr)
            except Exception,e:
                logging.getLogger('').error('Cannot send truncated notification email to %s:\n\t%s' % (toaddr,e))            
        server.quit()
        

    def exit(self, code):        
        self.code=code
        self.sendemail()
        sys.exit(code)
