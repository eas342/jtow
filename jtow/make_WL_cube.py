# coding: utf-8
import glob
import numpy as np
from astropy.io import fits, ascii
import tqdm
import os
import matplotlib.pyplot as plt
import pdb
from copy import deepcopy
from sklearn.decomposition import PCA

relFileSearch = 'Documents/jwst/flight_data/proc/01274/nrca3_fenrir_proc_002_groebak/split_output/*.fits'
defFileSearch = os.path.join(os.environ['HOME'],relFileSearch)

class wlcubeMake(object):
    def __init__(self,fileSearch=defFileSearch):
        """
        Initialize WL cube object
        """
        self.fileList = np.sort(glob.glob(fileSearch))
        self.nImg = len(self.fileList)
        
        self.outDir = os.path.join(os.path.split(fileSearch)[0],'WL_cube')
        self.firstPath = self.fileList[0]
        self.firstHead = fits.getheader(self.firstPath)
        
        if os.path.exists(self.outDir) == False:
            os.mkdir(outDir)

        outName = os.path.basename(self.firstPath).replace('.fits','_WL_cube.fits')
        self.outPath = os.path.join(self.outDir,outName)
        pcaName = os.path.basename(self.firstPath).replace('.fits','PCA_results.fits')
        self.pcaPath = os.path.join(self.outDir,pcaName)
        
        

    def make_WL_cube(self):
        firstHead = self.firstHead
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
        
        wlCube = np.zeros([self.nImg,y2-y1,x2-x1])
        
        for ind in tqdm.tqdm(np.arange(self.nImg)):
            oneFile = self.fileList[ind]
            
            HDUList = fits.open(oneFile)
            thisImg = HDUList['SCI'].data
            
            wlCube[ind] = thisImg[y1:y2,x1:x2]
            HDUList.close()
            
        


        HDUList_out = fits.PrimaryHDU(wlCube,self.firstHead)
        HDUList_out.writeto(self.outPath,overwrite=True)

    def run_pca(self,n_components=8,hideImg=[]):
        cube = fits.getdata(self.outPath)

        ## replace NaN with 0
        nanpt = np.isfinite(cube) == False
        cube[nanpt] = 0.0
        ## hide some points (from cosmic rays)
        ## replace them with the median image
        for onePt in hideImg:
            cube[onePt,:,:] = np.median(cube,axis=0)

        nz, ny, nx = cube.shape
        dat2D = np.reshape(cube,[nz,ny * nx])
        pca = PCA(n_components=n_components)
        pca.fit(dat2D)
        principalComponents = pca.fit_transform(dat2D)
        pcaImgCube = np.reshape(pca.components_,[n_components,ny,nx])
        
        pcaHDU = fits.PrimaryHDU(pcaImgCube,self.firstHead)
        tserHDU = fits.ImageHDU(principalComponents)
        outHDUList = fits.HDUList([pcaHDU,tserHDU])
        print("Writing PCA results to {}".format(self.pcaPath))
        outHDUList.writeto(self.pcaPath,overwrite=True)

    def run_all(self):
        self.make_WL_cube()
    

def make_WL_cube(fileSearch=defFileSearch):
    wlCubeObj = wlCubeMake()
    wlCubeObj.run_all()
