#!/usr/bin/env python
# coding: utf-8

# In[48]:


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table
from copy import deepcopy
import os
import pdb
import glob
import warnings

# In[31]:
homeDir = os.environ['HOME']
relPath = 'Documents/jwst/flight_data/01442_reproc/jw01442001001_04103_00001-seg006_nrcalong_uncal.fits'
def_uncal_path = os.path.join(homeDir,relPath)

max_nints = 10

def mini_split_one_seg(uncal_path=def_uncal_path,diagnostics=False,
                       reDo=False):
    print("Making mini segments for {}".format(uncal_path))
    
    HDUList = fits.open(uncal_path)
    #HDUList.info()
    
    origHeader = HDUList[0].header
    
    origDir = os.path.split(uncal_path)[0]
    outDir = os.path.join(origDir,'miniseg')
    if os.path.exists(outDir) == False:
        os.mkdir(outDir)
    
    orig_name = os.path.basename(uncal_path)
    
    nint_orig = origHeader['INTEND'] - origHeader['INTSTART'] + 1
    
    num_mini = int(np.ceil(float(nint_orig)/float(max_nints)))
    
    mini_starts_rel = np.arange(0,nint_orig,max_nints)
    mini_ends_rel = mini_starts_rel + max_nints
    mini_ends_rel[-1] = nint_orig
    
    if diagnostics == True:
        print("intstart={}".format(origHeader['INTSTART']))
        print("intend={}".format(origHeader['INTEND']))
        print("mini_starts_rel=")
        print(mini_starts_rel)
        print("mini_ends_rel=")
        print(mini_ends_rel)
        
    table = Table(HDUList['INT_TIMES'].data)
        
    for oneMini in np.arange(num_mini):
        tmp = orig_name
        outName = "_".join(tmp.split('_')[0:3]) + "_miniseg{:03d}_".format(oneMini) + "_".join(tmp.split('_')[3:])
        outPath = os.path.join(outDir,outName)
        
        if (os.path.exists(outPath) == True) & (reDo == False):
            print("Already found, skipping {}".format(outName))
        else:
            shortHDU_prim = deepcopy(fits.PrimaryHDU(None,origHeader))
            relstart = mini_starts_rel[oneMini]
            relend =  mini_ends_rel[oneMini]
            shortHDU_prim.header['INTSTART'] = origHeader['INTSTART'] + relstart
            shortHDU_prim.header['INTEND'] = origHeader['INTSTART'] + relend - 1
            #shortHDU_prim.header['INTSTART'] + nint_use - 1
            
            #HDUList['SCI'].data.shape
            shortHDU_sci = fits.ImageHDU(HDUList['SCI'].data[relstart:relend],HDUList['SCI'].header)
            
            shortHDU_group = HDUList['GROUP']
            shortHDU_asdf = HDUList['ASDF']
            
            int_timesHDU = fits.BinTableHDU(table[relstart:relend],header=HDUList['INT_TIMES'].header)
            
            shortHDUList = fits.HDUList([shortHDU_prim,shortHDU_sci,HDUList['GROUP'],
                                         int_timesHDU,HDUList['ASDF']])
            
            if 'ZEROFRAME' in HDUList:
                shortHDU_zero = fits.ImageHDU(HDUList['ZEROFRAME'].data[relstart:relend],HDUList['ZEROFRAME'].header)
                shortHDUList.append(shortHDU_zero)
            
            if 'REFOUT' in HDUList:
                shortHDU_refout = fits.ImageHDU(HDUList['REFOUT'].data[relstart:relend],HDUList['REFOUT'].header)
                shortHDUList.append(shortHDU_refout)
            
            shortHDUList.writeto(outPath,overwrite=True)
            
            shortHDUList.close()

    HDUList.close()

# In[ ]:

def loop_minisegments(fileSearch):
    fileList = np.sort(glob.glob(fileSearch))
    for oneFile in fileList:
        if os.path.exists(oneFile) == False:
            ## for symbolic links, check if they exist
            warnings.warn("File {} does not exist".format(oneFile))
        else:
            mini_split_one_seg(oneFile)


