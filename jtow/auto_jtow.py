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
    def __init__(self,searchString,iteration=1,
                 skipExisting=False):
        """
        An object to run an automatic pipeline wrapper run
        
        """
        self.gather_files(searchString)
        self.iteration = iteration
        self.skipExisting = skipExisting

    
    def gather_files(self,searchString):
        """
        Gather the files to run jtow on
        """
        self.fileList = np.sort(glob.glob(searchString))
        
        if len(self.fileList) == 0:
            warnings.warn("No files found at {}".format(searchString))
    
    def predicted_output_name(self,uncalName):
        """
        Predicted the output name from the jtow stage 1
        """
        return uncalName.replace('_uncal.fits','_0_rampfitstep.fits')

    def find_rate_file(self,oneFile):
        """
        Find the default rate file but check if it exists
        """
        if self.iteration == 1:
            rate_file_guess = oneFile.replace('_uncal.fits','_rate.fits')
        else:
            rate_file_guess = self.predicted_output_name(oneFile)
        
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
            rate_img = fits.getdata(rate_file,extname='SCI')
            valid_pt = np.isfinite(err_est)
            nonzero = err_est[valid_pt] > 0
            
            self.backg_est = np.nanmedian(rate_img)
            self.err_est = np.percentile(err_est[valid_pt][nonzero],1)
            if self.iteration == 1:
                rateThreshold = 5. * self.err_est + self.backg_est * 1.3
            else:
                rateThreshold = 12. * self.err_est
        
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
        directParam['ROEBApreserveBackg'] = True
        #directParam['custBias'] = 'selfBias' ## for now until a good set is prepped
        directParam['jumpRejectionThreshold'] = 6.0
        directParam['ROEBAmaskGrowthSize'] = 7
        directParam['maxCores'] = 'quarter'
        
        
        return directParam
    
    def run_jtow(self):
        """ Run the jtow in a for loop """
        
        for ind,oneFile in enumerate(self.fileList):
            if self.skipExisting == True:
                predicted_output = self.predicted_output_name(oneFile)
                if os.path.exists(predicted_output):
                    skipFile = True
                    print("Alread found {}".format(predicted_output))
                else:
                    skipFile = False
            else:
                skipFile = False
            if skipFile == True:
                pass
            else:
                directParam = self.set_up_parameters(oneFile)
                self.run_jtow_one(directParam)
    
    def run_jtow_one(self,directParam):
        jw = jtow.jw(directParam=directParam)
        jw.run_jw()

def run_auto_jtow(searchString,iterations=1,skipExisting=False):
    for iteration in np.arange(iterations) + 1:
        aj = auto_jtow(searchString,iteration=iteration,
                       skipExisting=skipExisting)
        aj.run_jtow()
    
    
