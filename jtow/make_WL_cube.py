# coding: utf-8
import glob
import numpy as np
from astropy.io import fits, ascii
import tqdm
import os
import matplotlib.pyplot as plt
import pdb
from copy import deepcopy


relFileSearch = 'Documents/jwst/flight_data/proc/01274/nrca3_fenrir_proc_002_groebak/split_output/*.fits'
defFileSearch = os.path.join(os.environ['HOME'],relFileSearch)

def make_WL_cube(fileSearch=defFileSearch,
                 saveFits=False,xLim=[950,1150]):
    
    fileList = np.sort(glob.glob(fileSearch))
    nImg = len(fileList)
    
    outDir = os.path.join(os.path.split(fileSearch)[0],'WL_cube')
    firstPath = fileList[0]
    firstHead = fits.getheader(firstPath)
    
    if os.path.exists(outDir) == False:
        os.mkdir(outDir)
    
    excpMessage = "No default coord set for {} {} {}".format(firstHead['FILTER'],firstHead['SUBARRAY'],
                                                             firstHead['DETECTOR'])
    if firstHead['FILTER'] == 'WLP4':
        if firstHead['SUBARRAY'] == 'SUBGRISM64':
            if firstHead['DETECTOR'] == 'NRCA3':
                x1, x2, y1, y2 = 1025, 1105, 0, 64
            elif firstHead['DETECTOR'] == 'NRCA1':
                x1, x2, y1, y2 = 1756, 1836, 0, 64
            else:
                raise Exception(excpMessage)
        else:
            raise Exception(excpMessage)
    else:
        raise Exception(excpMessage)
    
    wlCube = np.zeros([nImg,y2-y1,x2-x1])
    
    for ind in tqdm.tqdm(np.arange(nImg)):
        oneFile = fileList[ind]
        
        HDUList = fits.open(oneFile)
        thisImg = HDUList['SCI'].data
        
        wlCube[ind] = thisImg[y1:y2,x1:x2]
        
    
    outName = os.path.basename(firstPath).replace('.fits','_WL_cube.fits')
    
    HDUList_out = fits.PrimaryHDU(wlCube,firstHead)
    HDUList_out.writeto(os.path.join(outDir,outName),overwrite=True)
