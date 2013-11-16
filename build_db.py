import MySQLdb
import pandas as pd
import pandas.io.sql as pdsql
import scipy.io as sio
import numpy as np

def SetupDB():
    # connect to database server
    db = MySQLdb.connect(host='localhost',
        user='root', passwd='')

    # create cursor
    cur = db.cursor()

    # create database and switch to it
    setupstr = """
        CREATE DATABASE bartc;
        USE bartc;
        CREATE TABLE spikes (patient TINYINT, dataset TINYINT, 
            channel TINYINT, unit TINYINT, time DECIMAL(7,3));
        """
    cur.execute(setupstr)

def ImportSpikes(ftup, datadir):
    pdir = 'patient' + str(ftup[0]).zfill(3)
    fname = ('times_' + str(ftup[0]) + '.' + str(ftup[1]) + '.plx' + 
        str(ftup[2]) + '.mat')
    fullname = datadir + pdir + '/' + fname

    dat = sio.loadmat(fullname)['cluster_class']
    
    unit = dat[:,0].astype('int')
    times = np.around(dat[:,1]/1000., decimals=3) #times are originally in ms
    sortord = np.argsort(times) #spikes aren't always time-sorted
    times = times[sortord]
    unit = unit[sortord]
    
    # now restrict to valid units:
    valid = unit == ftup[3]
    times = times[valid]
    unit = unit[valid]

    ddict = {'patient': ftup[0], 'dataset': ftup[1], 
    'channel': ftup[2], 'unit': unit, 'time': times}
    df = pd.DataFrame(ddict)

    db = MySQLdb.connect(host='localhost',
        user='root', passwd='', db='bartc')
    pdsql.write_frame(df, con=db, name='spikes',
        if_exists='append', flavor='mysql')
    db.close()


# read data
# dat = pdsql.read_frame(qstr, db)

# # write data
# pdsql.write_frame(dat, con=db, name='supp',
#     if_exists='replace', flavor='mysql')

########### tests ##################
# patient = 18
# dset = 1
# chan = 1
# unit = 2

# ftup = (patient, dset, chan, unit)


if __name__ == '__main__':

    # build database and tables
    SetupDB()

    # locations of relevant files
    ddir = '/home/jmp33/data/bartc/plexdata/'
    spkfile = '/home/jmp33/code/hephys/valid_units.csv'

    # get list of tuples with valid channels
    ulist = []

    with open(spkfile) as infile:
        for line in infile:
            ulist.append(tuple(map(int, line.split(','))))

    # iterate through files, loading data
    for ftup in ulist:
        ImportSpikes(ftup, ddir)

