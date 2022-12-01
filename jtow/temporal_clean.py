#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import interp1d
from copy import deepcopy
import os
import pdb
import tqdm

defSearchPath = os.path.join(os.environ['HOME'],
                             ('Documents/jwst/flight_data/proc/01185/'+
                             'gj3470_obs016_nrca3_proc_001_example/*.fits'))

savePath = os.path.join(os.environ['HOME'],
                       ('Documents/jwst/flight_data/proc/01185/'+
                       'gj3470_obs016_nrca3_proc_001_cleaned/'))

fileList = np.sort(glob.glob(defSearchPath))

def clean_data(searchPath,savePath,thresholdSize=20,
               chunkSize=150):
    """
    Find a series of images using a the search path.
    Find points beyond the rejection threshold and replace with 
    linear interpolation along the temporal axis
                        
    Inputs
    ------
    fileList: list of str
        List of file names to go through. Must all be the same size
    savePath: str
        A directory where to save files
    thresholdSize: float
        The rejection threshold in the temporal direction (sigma)
    chunkSize: int
        How many images to examine at a time? Uses np.array_split so may not be exact
    """
    fileList = np.sort(glob.glob(searchPath))
    
    numChunks = len(fileList) // chunkSize + 1
    fileChunks = np.array_split(fileList,numChunks)
    
    for chunk in tqdm.tqdm(fileChunks):
        clean_from_fileList(chunk,savePath=savePath,thresholdSize=thresholdSize)
    

def clean_from_fileList(fileList,savePath,
                        thresholdSize = 20):
    """
    Clean a file list of images using a temporal rejection and 
    linear interpolation replacement scheme along the temporal axis
                        
    Inputs
    ------
    fileList: list of str
        List of file names to go through. Must all be the same size
    savePath: str
        A directory where to save files
    thresholdSize: float
        The rejection threshold in the temporal direction (sigma)
    """
    
    oneImg = fits.getdata(fileList[0])
    oneErr = fits.getdata(fileList[0],extname='ERR')
    nFile = len(fileList)
    
    ## Read in a cube of images
    imgCube = np.zeros([nFile,oneImg.shape[0],oneImg.shape[1]])    
    for ind,oneFile in enumerate(fileList):
        imgCube[ind] = fits.getdata(oneFile)
    
    ## find outliers
    medImg = np.nanmedian(imgCube,axis=0)
    
    diffCube = np.abs(imgCube - medImg) / oneErr
    
    outliers = diffCube > thresholdSize
    
    numOutliers = np.sum(outliers)
    numOutliers
    
    outlierMap = np.sum(outliers,axis=0)
    numOutlierUniquePixels = np.sum(outlierMap > 0)
    
    outY, outX = np.where(outlierMap > 0)
    
    ## make a cleaned cube with interpolation
    cleanedCube = deepcopy(imgCube)
    
    
    intNumArr = np.arange(nFile)
    for ptLook in np.arange(numOutlierUniquePixels):
    #    pltData = plt.plot(intNumArr,diffCube[:,outY[ptLook],outX[ptLook]])
    #    thisCol = pltData[0].get_color()
        outPt = outliers[:,outY[ptLook],outX[ptLook]]
        goodPt = (outPt == False)
    #    plt.plot(intNumArr[outPt],
    #             diffCube[outPt,outY[ptLook],outX[ptLook]],'x',color=thisCol)
        yGood = imgCube[goodPt,outY[ptLook],outX[ptLook]]
        if np.sum(goodPt) > 10: ## need a least 10 points for interpolation
            f_interp = interp1d(intNumArr[goodPt],yGood,
                                bounds_error=False,fill_value=np.nanmedian(yGood))
            cleanedCube[outPt,outY[ptLook],outX[ptLook]] = f_interp(intNumArr[outPt])
    #    plt.plot(intNumArr[outPt],cleanedCube[outPt,outY[ptLook],outX[ptLook]],'o',
    #             color=thisCol)
    
    if os.path.exists(savePath) == False:
        os.makedirs(savePath)
    
    for ind,oneFile in enumerate(fileList):
        baseName = os.path.splitext(os.path.basename(oneFile))[0]
        outName = "{}_cln.fits".format(baseName)
        outPath = os.path.join(savePath,outName)
        HDUList = fits.open(oneFile)
        HDUList[0].header['TCLEANED'] = (True, 'Is the data temporally cleaned')
        HDUList[0].header['CLTHRESH'] = (thresholdSize, 'Temporal Cleaning threshold (sigma)')
        HDUList['SCI'].data = cleanedCube[ind]
        HDUList.writeto(outPath,overwrite=True)
        HDUList.close()
    
    # In[16]:

