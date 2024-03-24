from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import csv
import numpy as np
import asdf
import astropy.units as u
import glob
import time
import yaml
import pdb
from scipy import ndimage
from copy import deepcopy
import pkg_resources
import os
import glob
import crds
import tqdm
from astropy.table import Table

def quick_ff_divide(searchPath,customFlat=None,
                    outputDir='ff_div'):
    """
    Divide the files by a flat field

    outputDir: str
        output directory name
    
    """
    
    fileList = np.sort(glob.glob(searchPath))
    
    for oneFile in tqdm.tqdm(fileList):
        head = fits.getheader(oneFile)
        
        ## make output path
        outDir = os.path.join(os.path.split(oneFile)[0],outputDir)
        if os.path.exists(outDir) == False:
            os.mkdir(outDir)

        if head['INSTRUME'] == 'NIRSPEC':
            instrumName = 'nirspec'
        else:
            instrumName = 'nircam'
        
        crds_path = os.path.join(os.environ['CRDS_PATH'],'references','jwst',instrumName)
        if 'PUPIL' in head:
            pupil = head['PUPIL']
        else:
            pupil = "none"
        
        if (head['DETECTOR'] == 'NRCALONG') & ("GRISM" in pupil):
            if head['FILTER'] == 'F444W':
                flatName = 'jwst_nircam_flat_0313.fits'
            elif head['FILTER'] == 'F322W2':
                flatName = 'jwst_nircam_flat_0266.fits'
            else:
                raise NotImplementedError("Have to add this filter {}".format(head['FILTER']))
        else:
            recs = crds.getrecommendations(head)
            if instrumName == 'nirspec':
                flatKey = 'dflat'
                flatName = recs['dflat']
            else:
                flatName = recs['flat']
                flatKey = 'flat'

        if customFlat is None:
            flatPath = os.path.join(crds_path,flatName)
            if os.path.exists(flatPath) == False:
                crds.getreferences(head,reftypes=[flatKey])
        else:
            flatPath = customFlat
        
        if instrumName == 'nirspec':
            flatCube = fits.getdata(flatPath)
            nPlanes = flatCube.shape[0]
            flatData = flatCube[nPlanes // 2]
            if head['SUBARRAY'] == 'SUB2048':
                if head['DETECTOR'] == 'NRS1':
                    customMask = True
                    customPx_x = [929 , 1861]
                    customPx_y = [ 18 , 9   ]
                else:
                    customMask = True
                    customPx_x = [1974, 1975, 1976, 1976, 1976, 1975, 1974, 1974]
                    customPx_y = [ 25,    25,  25 ,   24,   23,   23,   23,   24]
            else:
                customMask = False
        else:
            flatData = fits.getdata(flatPath)
            customMask = False
        
        if customMask == True:
            customBadPxTable = Table()
            customBadPxTable['x'] = customPx_x
            customBadPxTable['y'] = customPx_y
            pxHDU = fits.BinTableHDU(customBadPxTable)
            pxHDU.name = 'PXMASK'
        
        outName = os.path.basename(oneFile).replace('.fits','_ff.fits')
        outPath = os.path.join(outDir,outName)
        
        HDUList = fits.open(oneFile)
        xStart = HDUList[0].header['SUBSTRT1'] - 1
        xEnd = xStart + HDUList[0].header['SUBSIZE1']
        yStart = HDUList[0].header['SUBSTRT2'] - 1
        yEnd = yStart + HDUList[0].header['SUBSIZE2']

        

        subFlat = flatData[yStart:yEnd,xStart:xEnd]
        HDUList['SCI'].data = HDUList['SCI'].data / subFlat
        badpt = (HDUList['DQ'].data & 2**0) > 0
        HDUList['SCI'].data[badpt] = np.nan
        if customMask == True:
            for badpxInd in np.arange(len(customPx_x)):
                x1, y1 = customPx_x[badpxInd], customPx_y[badpxInd]
                HDUList['SCI'].data[y1,x1] = np.nan
            HDUList.append(pxHDU)
        
        HDUList['SCI'].header['FFNAME'] = (os.path.basename(flatPath),
                                           'manual flat field file')
        HDUList['SCI'].header['FFDIV'] = (True,
                                          'Is the science frame divided by a flat?')
        HDUList['SCI'].header['CUSTMSK'] = (customMask,
                                          'Is the science frame divided by a flat?')

        HDUList.writeto(outPath,overwrite=True)
        HDUList.close()
    
