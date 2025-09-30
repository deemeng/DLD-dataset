import logging

class MyHandler(logging.StreamHandler):

    def __init__(self):
        logging.StreamHandler.__init__(self)
        fmt = '%(asctime)s %(filename)-18s %(levelname)-8s: %(message)s'
        fmt_date = '%Y-%m-%dT%T%Z'
        formatter = logging.Formatter(fmt, fmt_date)
        self.setFormatter(formatter)

def logger(logname):

    '''
    config log infomation.
    params:
        logname - file name
    '''
    logging.basicConfig(filename=logname,
                        filemode='a',
                        level=logging.DEBUG)
    log = logging.getLogger('root')
    log.addHandler(MyHandler())

    # log.debug('debug message') 
    # log.info('info message')
    # log.warning('warning message')
    # log.error('error message')
    return log
