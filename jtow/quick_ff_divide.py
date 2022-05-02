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

def quick_ff_divide(searchPath):
    """
    Divide the files by a flat field
    """
    
    fileList = np.sort(glob.glob(searchPath))
    
    for oneFile in tqdm.tqdm(fileList):
        head = fits.getheader(oneFile)
        
        ## make output path
        outDir = os.path.join(os.path.split(oneFile)[0],'ff_div')
        if os.path.exists(outDir) == False:
            os.mkdir(outDir)
        
        nc_path = os.path.join(os.environ['CRDS_PATH'],'references','jwst','nircam')
        if (head['DETECTOR'] == 'NRCALONG') & ("GRISM" in head['PUPIL']):
            if head['FILTER'] == 'F444W':
                flatName = 'jwst_nircam_flat_0313.fits'
            elif head['FILTER'] == 'F322W2':
                flatName = 'jwst_nircam_flat_0266.fits'
            else:
                raise NotImplementedError("Have to add this filter {}".format(head['FILTER']))
        else:
            recs = crds.getrecommendations(head)
            flatName = recs['flat']
        
        flatPath = os.path.join(nc_path,flatName)
        flatData = fits.getdata(flatPath)
        
        outName = os.path.basename(oneFile).replace('.fits','_ff.fits')
        outPath = os.path.join(outDir,outName)
        
        HDUList = fits.open(oneFile)
        xStart = HDUList[0].header['SUBSTRT1'] - 1
        xEnd = xStart + HDUList[0].header['SUBSIZE1']
        yStart = HDUList[0].header['SUBSTRT2'] - 1
        yEnd = yStart + HDUList[0].header['SUBSIZE2']
        
        subFlat = flatData[yStart:yEnd,xStart:xEnd]
        HDUList['SCI'].data = HDUList['SCI'].data / subFlat
        badpt = (HDUList['DQ'].data & 2**0)
        HDUList['SCI'].data[badpt] = np.nan
        HDUList.writeto(outPath,overwrite=True)
        HDUList.close()
    