import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table
from copy import deepcopy
import os
#from jwst import datamodels
from copy import deepcopy
import glob
import pdb
from astropy.time import Time
import astropy.units as u
import tqdm


def make_pairs(fileSearch = 'proc_nrca3_p012_flight_singleBias/jw*_uncal_jumpstep.fits'):

    fileList = np.sort(glob.glob(fileSearch))

    origDir = os.path.split(fileList[0])[0]
    outDir = os.path.join(origDir,'cds_split')
    if os.path.exists(outDir) == False:
        os.makedirs(outDir)

    for oneFile in fileList:
        #oneFile = fileList[-1]
    
    
        HDUList = fits.open(oneFile)
        origHead = HDUList[0].header
        dat4D = HDUList['SCI'].data
        nint, ngroup, ny, nx = dat4D.shape
    
        texp_start = Time(origHead["EXPSTART"],format='mjd')
    
        print("Splitting file {} of {}".format(oneFile,len(fileList)))
        for oneInt in tqdm.tqdm(np.arange(nint)):
            intNum = origHead['INTSTART'] + oneInt
            intHDU = HDUList['INT_TIMES']
    
            intStart_appr = texp_start + (origHead['TFRAME'] + origHead['EFFINTTM']) * intNum * u.second
            intStart_better_UT = Time(HDUList['INT_TIMES'].data['int_start_MJD_UTC'][oneInt],format='mjd')
            intStart_better_BJD = Time(HDUList['INT_TIMES'].data['int_start_BJD_TDB'][oneInt],format='mjd')
    
            for oneGroup in np.arange(ngroup-1):
                diff = dat4D[oneInt,oneGroup+1,:,:] - dat4D[oneInt,oneGroup,:,:]
                cds = diff / origHead['TGROUP']
                primHDU = fits.PrimaryHDU(None,deepcopy(HDUList[0].header))
                primHDU.header['CDSB'] = (oneGroup+2,'minuend group (B) in B-A, 1-based')
                primHDU.header['CDSA'] = (oneGroup+1,'subtrahend group (A) in B-A, 1-based')
                primHDU.header['ON_INT'] = (intNum,'Current integration (1-based counting)')
    
                ## approximate time becuase the INT_TIMES were off in early data
                ## also doesn't include the fast frame resets but those are really small corrections
                groupOffset = (oneGroup + 0.5) * origHead['TGROUP'] * u.second
                thisT = intStart_appr + groupOffset
                primHDU.header['DATE-OBS'] = (thisT.fits,'approximate UT time,backup to INT_TIMES')
            
                ## time from INT_TIMES extension
                thisT2 = intStart_better_UT + groupOffset
                primHDU.header['DATEGRUT'] = (thisT2.fits,'UT time from INT_TIMES extension + group times')
                thisT3 = intStart_better_BJD + groupOffset
                primHDU.header['DATEGBJD'] = (thisT3.fits,'BJD time from INT_TIMES extension + group times')
            
                hduSCI = fits.ImageHDU(cds)
                hduSCI.header['BUNIT'] = ('DN/s','count rate')
                hduSCI.name = 'SCI'
                err_diff = np.sqrt(HDUList['ERR'].data[oneInt,oneGroup+1,:,:]**2 +
                                   HDUList['ERR'].data[oneInt,oneGroup,:,:]**2)
                errHDU = fits.ImageHDU(err_diff/origHead['TGROUP'])
                errHDU.header['BUNIT'] = ('DN/s','count rate')
                errHDU.name = 'ERR'
            
                dqCDS = (HDUList['GROUPDQ'].data[oneInt,oneGroup+1,:,:] |
                         HDUList['GROUPDQ'].data[oneInt,oneGroup,:,:])
                dqCDScomb = HDUList['PIXELDQ'].data[:,:] | dqCDS
                dqHDU = fits.ImageHDU(dqCDS)
                dqHDU.name = 'DQ'
    
                outHDUList = fits.HDUList([primHDU,hduSCI,errHDU,dqHDU])
    
                outName_pref = os.path.splitext(os.path.basename(oneFile))[0]
                outName = "{}_cds_int_{:04d}_grp_{:03d}.fits".format(outName_pref,intNum,oneGroup+1)
                outPath = os.path.join(outDir,outName)
                outHDUList.writeto(outPath,overwrite=True)
            
        HDUList.close(oneFile)
    



