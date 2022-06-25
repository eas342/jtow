#Download the super short GRISM256 data and its supporting SuperBias Files.

#Import Library

#For Zenodo downloading
import urllib.request
import shutil

#General
import os
import yaml
import glob

#unittest
import unittest
from jtow import jtow

#Get the raw uncalibrated data for NIRCam short-wavelength detector

rawdata_filename_ext = 'jw88888001001_01101_00002_nrca5_uncal.fits' #raw data Name

cwd = os.getcwd() #Get the current working directory 
data_dir = cwd+'/jtow_test_results/jtow_nrca5_custBias_selfBias/Raw_Data/' #Output directroy for download data

#Make a directory to store the test data
if (os.path.exists(data_dir) == False):
    os.makedirs(data_dir)

rawdata_file = data_dir+rawdata_filename_ext #File for the raw data 
if (os.path.exists(rawdata_file) == False):
    raw_data = urllib.request.urlopen('https://zenodo.org/record/6688451/files/{}'.format(rawdata_filename_ext)) #Zenodo link
    
    with open('{}'.format(rawdata_filename_ext),'wb') as rawdata_out:
         rawdata_out.write(raw_data.read())
            
    shutil.move(str(rawdata_filename_ext), str(data_dir + '/' + rawdata_filename_ext)) #Save the raw data 
    
#Make a bias file
#uncalfile = data_dir+rawdata_filename_ext
#HDUList = fits.open(file)

#image2D = HDUList[1].data[0][0] #Acting as Bias File
jtow_Params_file = cwd+'/jtow_test_results/jtow_nrca5_custBias_selfBias/jtow_nrca5_custBias_selfBias_params.yaml' #File for the param file 
if (os.path.exists(jtow_Params_file) == False):
    jtow_Params = {'rawFileSearch': rawdata_file,
                   'outputDir': cwd+'/jtow_test_results/jtow_nrca5_custBias_selfBias/',
                   'photParam': None,
                   'noutputs': None,
                   'add_noutputs_keyword': False, 
                   'ROEBACorrection': True,
                   'autoROEBAmasks': True,
                   'maxCores': 'quarter',
                   'ROEBAmaskfromRate': None,
                   'ROEBAmaskfromRateThreshold': 0.5,
                   'custBias': 'selfBias',
                   'saveBiasStep': False,
                   'saveROEBAdiagnostics': False,
                   'jumpRejectionThreshold': 15.0,
                   'ROEBAmaskGrowthSize': None,
                   'saveJumpStep': False,
                   'biasCycle':  None,
                   'biasCycleSearch': None}
                
    #Write the test parameter file
    with open(cwd+'/jtow_test_results/jtow_nrca5_custBias_selfBias/jtow_nrca5_custBias_selfBias_params.yaml', 'w') as file:
        paramfile = yaml.dump(jtow_Params, file)
    
#Test jtow
class Test_JTOW(unittest.TestCase):
   
    def setUp(self):
        """Set up test fixtures, if any."""
        
        self.examParamPath = cwd+'/jtow_test_results/jtow_nrca5_custBias_selfBias/jtow_nrca5_custBias_selfBias_params.yaml' #Path to YAML file
       
        self.jtow_obj = jtow.jw(paramFile = self.examParamPath) #Assign YAML file to JTOW object
    
    def test_run_jw(self):
        self.jtow_obj.run_jw()