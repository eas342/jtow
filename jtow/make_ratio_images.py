# coding: utf-8
import glob
import numpy as np
from astropy.io import fits, ascii
import tqdm
import os
import matplotlib.pyplot as plt
import pdb
from copy import deepcopy

defRefImg = '/fenrirdata1/es_tso/flight_data/01366/nrca3_proc_001/split_output/jw01366002001_04103_00001-seg001_nrca3_uncal_1_rampfitstep_I00117.fits'
#defRefImg = 'split_output/jw01442001001_04103_00001-seg001_nrca3_uncal_1_rampfitstep_I00005.fits'

#defFileSearch = 'split_output/jw01442001001_04103_00001-seg00*.fits'
defFileSearch = '/fenrirdata1/es_tso/flight_data/01366/nrca3_proc_001/split_output/*.fits'
defOutDir = '/fenrirdata1/es_tso/flight_data/01366/nrca3_proc_001/ratio_images'

def make_ratio_images(refImage=defRefImg,fileSearch=defFileSearch,outDir=defOutDir,
                      saveFits=False,xLim=[950,1150]):
    refImg = fits.getdata(refImage)
    fileList = np.sort(glob.glob(fileSearch))
    #os.mkdir('rel_images')
    
    if os.path.exists(outDir) == False:
        os.mkdir(outDir)
    if saveFits == True:
        fitsOutDir = os.path.join(outDir,'ratio_data')
        if os.path.exists(fitsOutDir) == False:
            os.mkdir(fitsOutDir)
    
    for ind in tqdm.tqdm(np.arange(len(fileList))):
        oneFile = fileList[ind]
        fig,ax = plt.subplots()
        HDUList = fits.open(oneFile)
        thisImg = HDUList['SCI'].data
        
        ratioImg = thisImg / refImg
        pdata = ax.imshow(ratioImg,vmin=0.98,vmax=1.02,origin='lower')
        
        ax.set_xlim(xLim[0],xLim[1])
        ax.set_title("Image {}".format(ind))
        
        outName = os.path.basename(oneFile).replace('.fits','.png')
        fig.colorbar(pdata)
        fig.savefig(os.path.join(outDir,outName))
        plt.close(fig)
        

        
        if saveFits == True:
            outNameFits = 'ratio_img_{}'.format(os.path.basename(oneFile))
            HDUList_out = deepcopy(HDUList)
            HDUList_out['SCI'].data = ratioImg
            HDUList_out.writeto(os.path.join(fitsOutDir,outNameFits),overwrite=True)
        
        HDUList.close()
    
