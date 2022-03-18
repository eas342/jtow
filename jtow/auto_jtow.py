from . import jtow
import glob
import numpy as np
import warnings
import os
import pkg_resources
from astropy.io import fits, ascii
from copy import deepcopy
import pdb

path_to_defaults = "params/default_params.yaml"
defaultParamPath = pkg_resources.resource_filename('jtow',path_to_defaults)
defaultParam = jtow.read_yaml(defaultParamPath)

class auto_jtow(object):
    def __init__(self,searchString):
        """
        An object to run an automatic pipeline wrapper run
        
        """
        self.gather_files(searchString)
        

    
    def gather_files(self,searchString):
        """
        Gather the files to run jtow on
        """
        self.fileList = np.sort(glob.glob(searchString))
        
        if len(self.fileList) == 0:
            warnings.warn("No files found at {}".format(searchString))
    
    def find_rate_file(self,oneFile):
        """
        Find the default rate file but check if it exists
        """
        rate_file_guess = oneFile.replace('_uncal.fits','_rate.fits')
        if os.path.exists(rate_file_guess):
            rate_file = rate_file_guess
        else:
            rate_file = None
            warnings.warn("No rate file found in the same place as _uncal.fits")
            warnings.warn("Proceed at your own risk")
        
        self.rate_file = rate_file
        return rate_file
    
    def get_roeba_threshold(self,rate_file):
        """
        Find the ROEBA threshold
        """
        if rate_file is None:
            rateThreshold = 5.
        else:
            err_est = fits.getdata(rate_file,extname='ERR')
            valid_pt = np.isfinite(err_est)
            nonzero = err_est[valid_pt] > 0
            
            rateThreshold = 5. * np.min(err_est[valid_pt][nonzero])
            
        return rateThreshold
        
    
    def set_up_parameters(self,oneFile):
        """
        Set up parameters
        """
        directParam = deepcopy(defaultParam)
        directParam['rawFileSearch'] = oneFile
        outPath_initial = os.path.split(oneFile)[0]
        if outPath_initial == '':
            directParam['outputDir'] = '.'
        else:
            directParam['outputDir'] = outPath_initial
        directParam['add_noutputs_keyword'] = False
        directParam['ROEBAmaskfromRate'] = self.find_rate_file(oneFile)
        directParam['ROEBAmaskfromRateThreshold'] = self.get_roeba_threshold(self.rate_file)
        directParam['custBias'] = 'selfBias' ## for now until a good set is prepped
        directParam['jumpRejectionThreshold'] = 3.0
        directParam['ROEBAmaskGrowthSize'] = 15
        directParam['maxCores'] = 'none'
        
        return directParam
    
    def run_jtow(self):
        """ Run the jtow in a for loop """
        
        for ind,oneFile in enumerate(self.fileList):
            directParam = self.set_up_parameters(oneFile)
            self.run_jtow_one(directParam)
    
    def run_jtow_one(self,directParam):
        jw = jtow.jw(directParam=directParam)
        jw.run_jw()

def run_auto_jtow(searchString):
    aj = auto_jtow(searchString)
    aj.run_jtow()
    
    