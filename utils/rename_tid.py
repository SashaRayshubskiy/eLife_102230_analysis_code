#! /usr/bin/python

import glob
import sys
import string
import os

pathtodata = '/data/drive_fast/sasha/151224_bfly_2/2p/'

files = glob.glob(pathtodata + '*.tif')

for filename in files:
    print 'Before: ' + filename

    fs = filename.split('_')
    # print fs[7]
    fs[7] = str(int(fs[7]) - 1)
    print 'After: ' + string.join( fs, '_' )
    new_filename = string.join( fs, '_' )
    os.rename(filename, new_filename)
