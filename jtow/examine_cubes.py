import os
from astropy.io import fits, ascii
import numpy as np
import glob


#fileSearch = '../obs001_hd189_nrca1_proc_001/split_output/ff_cleaned/*.fits'
#fileList = np.sort(glob.glob(fileSearch))

def avg_image(fileList,outName='data/avg_img.fits'):
    oneImg = fits.getdata(fileList[0])
    firstHead = fits.getheader(fileList[0])
    oneErr = fits.getdata(fileList[0],extname='ERR')
    nFile = len(fileList)
    
    ## Read in a cube of images
    imgCube = np.zeros([nFile,oneImg.shape[0],oneImg.shape[1]])
    for ind,oneFile in enumerate(fileList):
        imgCube[ind] = fits.getdata(oneFile)
        
    ## find average
    avgImg = np.nanmean(imgCube,axis=0)

    ## save result
    outHDU = fits.PrimaryHDU(avgImg,firstHead)
    outHDU.writeto(outName,overwrite=True)

def do_avg(fileSearch,startInd=1007,endInd=1017,saveDir='avg_img_data'):
    """
    find the average image for a give start and end index
    """
    fileList = np.sort(glob.glob(fileSearch))
    
    imgSubset = fileList[startInd:endInd]
    baseName = os.path.splitext(os.path.basename(fileList[0]))[0]
    outName = 'avg_{}_in_{:04d}_{:04d}.fits'.format(baseName,startInd,endInd)
    
    outPath = os.path.join(saveDir,outName)
    if os.path.exists(saveDir) == False:
        os.makedirs(saveDir)
    
    avg_image(imgSubset,outName=outPath)


def do_specific_avg():
    startIndList = [5,1007,1500]
    endIndList=[1000,1494,2743]

    for i in np.arange(3):
        do_avg(startIndList[i],endIndList[i])
    
            
