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
            os.mkdir(self.outDir)

        outName = os.path.basename(self.firstPath).replace('.fits','_WL_cube.fits')
        self.outPath = os.path.join(self.outDir,outName)
        pcaName = os.path.basename(self.firstPath).replace('.fits','PCA_results.fits')
        self.pcaPath = os.path.join(self.outDir,pcaName)
        self.pcaPlotDir = os.path.join(self.outDir,'plots')
        if os.path.exists(self.pcaPlotDir) == False:
            os.mkdir(self.pcaPlotDir)
        
        if self.nImg < 13:
            self.default_n_components = 6
        else:
            self.default_n_components = 13

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
            if firstHead['SUBARRAY'] == 'SUBGRISM256':
                if firstHead['DETECTOR'] == 'NRCA3':
                    x1, x2, y1, y2 = 994, 1128, 87, 247
                elif firstHead['DETECTOR'] == 'NRCA1':
                    x1, x2, y1, y2 = 1727, 1861, 82, 242
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

    def run_pca(self,n_components=None,hideImg=[]):
        cube = fits.getdata(self.outPath)

        if n_components is None:
            n_components = self.default_n_components

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

    def plot_pca_tser(self,renorm=False):
        HDUList = fits.open(self.pcaPath)
        tser2D = HDUList[1].data
        fig, ax = plt.subplots(figsize=(10,18))
        if renorm == True:
            offsetVal = 7
            ylabel = 'Normalized PC Value + offset'
            labelOffset =  0.5 * offsetVal
        else:
            offsetVal = 3000
            labelOffset =  0.3 * offsetVal
            ylabel = "PC Value + offset"
        for ind in np.arange(tser2D.shape[1]):
            offset = ind * offsetVal
            yShow = tser2D[:,ind]
            if renorm == True:
                yShow = yShow / np.std(yShow)
            
            lineDat = plt.plot(yShow - offset)
            ax.text(0,-offset + labelOffset,'Comp {}'.format(ind),
                    color=lineDat[0].get_color())
        ax.set_xlabel("Time (Int Number)")
        ax.set_ylabel(ylabel)
        ax.legend()
        figPath = os.path.join(self.pcaPlotDir,'eigenpca_norm_{}.pdf'.format(renorm))
        fig.savefig(figPath,bbox_inches='tight')

        HDUList.close()


    def plot_eigenimages(self):
        HDUList = fits.open(self.pcaPath)
        EigenImages = HDUList[0].data
        
        fig, axArr2D = plt.subplots(4,4,gridspec_kw={'wspace':0.2,'hspace':0.4},
                            figsize=(12,12))
        axArr = axArr2D.ravel()
        for ind,oneImg in enumerate(EigenImages):

            ax = axArr[ind]
            ax.imshow(oneImg,origin='lower')
            ax.set_title('Comp {}'.format(ind))
        
        figPath = os.path.join(self.pcaPlotDir,'eigenimages.pdf')
        print('Saving PCA figure to {}'.format(figPath))
        fig.savefig(figPath,bbox_inches='tight')
        HDUList.close()

    def plot_pca(self):
        self.plot_eigenimages()
        for renormalize in [False,True]:
            self.plot_pca_tser(renorm=renormalize)


    def run_all(self):
        self.make_WL_cube()
        self.run_pca()
        self.plot_pca()
    

def make_WL_cube(fileSearch=defFileSearch):
    wlCubeObj = wlcubeMake(fileSearch=fileSearch)
    wlCubeObj.run_all()
