import MySQLdb
import pandas as pd
import pandas.io.sql as pdsql
import scipy.io as sio
import numpy as np

# connect to database
db = MySQLdb.connect(host='localhost',
    user='root', passwd='', db='bartc')

# create cursor
cur = db.cursor()

# make a test query
qstr = """
    SELECT * FROM spikes WHERE time < 1 
    """

# read data
dat = pdsql.read_frame(qstr, db)

# write data
pdsql.write_frame(dat, con=db, name='supp',
    if_exists='replace', flavor='mysql')

########### tests ##################
patient = 18
dset = 1
chan = 1
unit = 2

ftup = (patient, dset, chan, unit)

ddir = '/home/jmp33/data/bartc/plexdata/'
pdir = 'patient' + str(ftup[0]).zfill(3)
fname = ('times_' + str(ftup[0]) + '.' + str(ftup[1]) + '.plx' + 
    str(ftup[2]) + '.mat')
fullname = ddir + pdir + '/' + fname

dat = sio.loadmat(fullname)['cluster_class']
unit = dat[:,0].astype('int')
times = dat[:,1]/1000. #times are originally in ms
sortord = np.argsort(times) #spikes aren't always time-sorted
times = times[sortord]
unit = unit[sortord]
ddict = {'patient': ftup[0], 'dataset': ftup[1], 
'channel': ftup[2], 'unit': unit, 'time': times}
df = pd.DataFrame(ddict)

db = MySQLdb.connect(host='localhost',
    user='root', passwd='', db='bartc')
pdsql.write_frame(df, con=db, name='spikes',
    if_exists='append', flavor='mysql')