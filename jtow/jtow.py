import os
#os.environ['CRDS_PATH'] = '/fenrirdata1/kg_data/crds_cache/' #These pathways should be defined in your ~./bash profile. If not, you can set them within the notebook.
#os.environ['CRDS_SERVER_URL']= 'https://jwst-crds.stsci.edu'
#os.environ['CRDS_CONTEXT']='jwst_0756.pmap' #Occasionally, the JWST CRDS pmap will be updated. Updates may break existing code. Use this command to revert to an older working verison until the issue is fixed.


import jwst
print(jwst.__version__) #Print what version of the pipeline you are using.

from jwst.pipeline.calwebb_detector1 import Detector1Pipeline #Stage 1
from jwst.pipeline.calwebb_image2 import Image2Pipeline #Stage 2
from jwst.pipeline.calwebb_tso3 import Tso3Pipeline #Stage 3
from jwst.associations.asn_from_list import asn_from_list #Association file imports
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase

#General
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

# Individual steps that make up calwebb_detector1
from jwst.dq_init import DQInitStep
from jwst.saturation import SaturationStep
from jwst.superbias import SuperBiasStep
from jwst.ipc import IPCStep                                                                                    
from jwst.refpix import RefPixStep                                                                
from jwst.linearity import LinearityStep
from jwst.persistence import PersistenceStep
from jwst.dark_current import DarkCurrentStep
from jwst.jump import JumpStep
from jwst.ramp_fitting import RampFitStep
from jwst import datamodels


# In[359]:


## ES custom pipeline
from tshirt.pipeline import phot_pipeline
from tshirt.pipeline.instrument_specific import rowamp_sub
import tqdm

path_to_defaults = "params/default_params.yaml"
defaultParamPath = pkg_resources.resource_filename('jtow',path_to_defaults)


def read_yaml(filePath):
    with open(filePath) as yamlFile:
        yamlStructure = yaml.safe_load(yamlFile)
    return yamlStructure

class jw(object):
    def __init__(self,paramFile=defaultParamPath,directParam=None):
        """
        An wrapper object to run the jwst pipeline (with some modifications)
        
        Parameters
        -----------
        paramFile: str
            Location of the YAML file containing the parameters.
        
        directParam: dict
            A dictionary of parameters. Overrides the paramFile
        """
        
        self.get_parameters(paramFile,directParam=directParam)
        
        defaultParams = read_yaml(defaultParamPath)
        
        for oneKey in defaultParams.keys():
            if oneKey not in self.param:
                self.param[oneKey] = defaultParams[oneKey]
        
        self.set_up_dirs()
        self.get_files()
        
        self.make_descrip()
        
        self.max_cores = self.param['maxCores']
        
        self.make_roeba_masks()
    

    
    def get_parameters(self,paramFile,directParam=None):
        if directParam is None:
            self.paramFile = paramFile
            self.param = read_yaml(paramFile)
        else:
            self.paramFile = 'direct dictionary'
            self.param = directParam
        
        if self.param['add_noutputs_keyword'] == True:
            warnings.warn("This code will modify the uncal file NOUTPUTS. This is DANGEROUS. Only use for older mirage simulations that lacked NOUTPUTS keyword")
        
    
    def set_up_dirs(self):
        """
        Set up directories
        """
        self.output_dir = self.param["outputDir"]
        
        self.diagnostic_dir = os.path.join(self.param["outputDir"],'diagnostics')
        
        ## make sure the directories exist
        for oneDir in [self.output_dir,self.diagnostic_dir]:
            if os.path.exists(oneDir) == False:
                os.makedirs(oneDir)
        
    
    def get_files(self):
        
        all_uncal_files = [] 
        #output_dir = '/fenrirdata1/es_tso/sim_data/mirage_035_hatp14_short_for_pipe_tests/stsci_proc/'
        #output_dir = '/fenrirdata1/es_tso/sim_data/mirage_032_hatp14_car33_no_backg/stsci_proc/'
        #output_dir = '/fenrirdata1/es_tso/sim_data/mirage_032_hatp14_car33_no_backg/stsci_proc_003_es_refcor/'
        #rawFileSearch = "/fenrirdata1/es_tso/sim_data/mirage_029_hd189733b_transit/raw/*nrca1_uncal.fits"
        #output_dir = '/fenrirdata1/es_tso/sim_data/mirage_029_hd189733b_transit/proc_roeba_nrca1/'
        #rawFileSearch = "/fenrirdata1/es_tso/sim_data/mirage_029_hd189733b_transit/raw/*nrca1_uncal.fits"
        #rawFileSearch = "/fenrirdata1/es_tso/sim_data/mirage_037_hatp14_lower_well_frac/raw/*nrca3_uncal.fits"
        #output_dir = "/fenrirdata1/es_tso/sim_data/mirage_037_hatp14_lower_well_frac/proc_roeba_nrca3"

        
        rawList = np.sort(glob.glob(self.param['rawFileSearch']))
        
        for fitsName in rawList: #Grabbing only these files from the directory
            if self.param['add_noutputs_keyword'] = True
                HDUList = fits.open(fitsName, 'update')
                #This was not input at the time of the simulation. Therefore, we manually must input this information.
                HDUList[0].header['NOUTPUTS'] = (self.param['noutputs'], 'Number of output amplifiers') 
                HDUList.close()
            all_uncal_files.append(fitsName)
    
        self.all_uncal_files = sorted(all_uncal_files) #sort files alphabetically.
    
    def make_descrip(self):
        """
        Make a description for diagnostics and saving info
        """
        self.firstUncal = os.path.basename(self.all_uncal_files[0])
        self.descrip = self.firstUncal.replace('_uncal.fits','')
    
    def make_roeba_masks(self):
        """
        Make masks for Row-by-row, odd-even by amplifier correction (ROEBA)
        
        """
        if self.param['autoROEBAmasks'] == True:
            firstHead = fits.getheader(self.all_uncal_files[0])
            firstHead_sci = fits.getheader(self.all_uncal_files[0],extname='SCI')
            Nx = firstHead_sci['NAXIS1']
            Ny = firstHead_sci['NAXIS2']
            if self.param['noutputs'] is None:
                if 'NOUTPUTS' in firstHead:
                    self.param['noutputs'] = firstHead['NOUTPUTS']
                else:
                    raise Exception("NOUTPUTS not found in first header. Try setting it manually with noutputs")
            
            if self.param['ROEBAmaskfromRate'] != None:
                HDUList = fits.open(self.param['ROEBAmaskfromRate'])
                rateDat = HDUList['SCI'].data
                
                self.photParam = None
                ROEBAmask = (rateDat < self.param['ROEBAmaskfromRateThreshold'])
                
                self.bad_dq_mask = HDUList['DQ'].data > 0
                
                HDUList.close()
                
            elif firstHead['PUPIL'] == 'GRISMR':
                grismsFilterList = ['F322W2','F444W']
                if firstHead['FILTER'] in grismsFilterList:
                    self.photParam = None
                    mask1 = np.ones([Ny,Nx],dtype=bool)
                    
                    #mask1[0:4,:] = False
                    mask1[:,0:4] = False
                    mask1[:,-4:] = False
                    if firstHead['FILTER'] == 'F444W':
                        mask1[12:57,794:2048] = False
                    elif firstHead['FILTER'] == 'F322W2':
                        mask1[12:57,28:1794] = False
                    else:
                        raise NotImplementedError
                    
                    ROEBAmask = mask1
                
                    self.bad_dq_mask = np.zeros_like(ROEBAmask,dtype=bool)
                
                else:
                    raise NotImplementedError
            ## NOTE TO SELF: I SHOULD CHECK ALL HEADERS, NOT JUST ONE!! Will fix this later
            
            elif firstHead['EXP_TYPE'] == 'NRC_TSIMAGE':
                if firstHead['PUPIL'] == 'WLP8':
                    backRadii = [100,101]
                elif (firstHead['PUPIL'] == 'CLEAR') & (firstHead['FILTER'] == 'WLP4'):
                    backRadii = [49,50]
                else:
                    backRadii = [12,13]
                
                xLoc = firstHead_sci['XREF_SCI']
                yLoc = firstHead_sci['YREF_SCI']
                
                #photParam = {'refStarPos': [[67.0 - 1.0,30.0 - 1.0]],'backStart':49,'backEnd': 50,
                self.photParam = {'refStarPos': [[xLoc-1,yLoc-1]],'backStart':backRadii[0],'backEnd': backRadii[1],
                                  'FITSextension': 1,
                                  'isCube': True,'cubePlane':0,'procFiles':'*.fits'}
                refpixMask = np.ones([Ny,Nx],dtype=bool)
                refpixMask[:,0:4] = False
                refpixMask[:,-4:] = False
                ROEBAmask = refpixMask
                
                self.bad_dq_mask = np.zeros_like(ROEBAmask,dtype=bool)
            else:
                raise Exception("Unrecognized header metadata to create an automatic ROEBA mask")
        else:
            self.photParam = None
            ROEBAmask = None
            self.bad_dq_mask = None
        
        if (self.param['ROEBAmaskGrowthSize'] is None) | (ROEBAmask is None):
            self.ROEBAmask = ROEBAmask
        else:
            self.ROEBAmask = self.grow_mask(ROEBAmask)
            
        if self.param['saveROEBAdiagnostics'] == True:
            self.save_roeba_masks()
        
    
    def save_diagnostic_img(self,diagnostic_img,suffix):
        """
        Save a diagnostic file
        """
        primHDU = fits.PrimaryHDU(np.array(diagnostic_img,dtype=int))
        outPath = os.path.join(self.diagnostic_dir,'{}_{}.fits'.format(self.descrip,suffix))
        print("Saving {} to {}".format(suffix,outPath))
        primHDU.writeto(outPath,overwrite=True)
        
        del primHDU
    
    def grow_mask(self,img):
        """
        Grow the mask to extend into the wings
        
        Parameters
        ----------
        img: numpy 2D array
            Mask image to be grown
        """
        
        ## construct a round tophat kernel, rounded to the nearest pixel
        growth_r = self.param['ROEBAmaskGrowthSize']
        ksize = int(growth_r * 2 + 4)
        y, x= np.mgrid[0:ksize,0:ksize]
        midptx = ksize/2. - 0.5
        midpty = midptx
        r = np.sqrt(((x-midptx)**2 + (y-midpty)**2))
        k = r < growth_r
        
        if self.param['saveROEBAdiagnostics'] == True:
            self.save_diagnostic_img(k,'roeba_mask_growth_kernel')
            self.save_diagnostic_img(img,'roeba_mask_before_growth')
        
        ## Keep in mind, the source=False and Backg=True (as of 2022-03-03)
        
        ## the source pixels are 0 = False, so we want to grow those
        ## but don't grow the bad pixels
        arr_to_convolve = (img == 0) & (self.bad_dq_mask == False)
        grown = (ndimage.convolve(np.array(arr_to_convolve,dtype=int), k, mode='constant', cval=0.0))
        
        # Now that we've grown the source pixels, we want to find the background pixels again
        # Maybe the mask should have been with the source=False, background=True from the start
        #, but this is the way it works currently (2022-03-03)
        initialROEBAmask = (grown == 0)
        # Have to add the original bad DQ mask in
        finalROEBAmask = initialROEBAmask & (img > 0)
        
        return finalROEBAmask 
        
    
    def save_roeba_masks(self):
        """
        Save the background mask used by ROEBA
        """
        if self.ROEBAmask is None:
            print("No ROEBA mask found, nothing to save")
        else:
            self.save_diagnostic_img(self.ROEBAmask,'roeba_mask')
            
    
    def run_jw(self):
        """
        Run the JWST pipeline for all uncal files
        """
        
        startTime = time.time() #Time how long this step takes
        
        for uncal_file in self.all_uncal_files:
                # Using the run() method. Instantiate and set parameters
            dq_init_step = DQInitStep()
            dq_init = dq_init_step.run(uncal_file)
            
            
            # ## Saturation Flagging
            # Using the run() method
            saturation_step = SaturationStep()
            # Call using the the output from the previously-run dq_init step
            saturation = saturation_step.run(dq_init)
            del dq_init ## try to save memory
    
            # Using the run() method
            superbias_step = SuperBiasStep()
            
            if self.param['custBias'] is None:
                pass
            elif self.param['custBias'] == 'selfBias':
                superbias_step.skip = True
                saturation.data = saturation.data - saturation.data[0][0]
            else:
                superbias_step.override_superbias = self.param['custBias']
            
            if self.param['saveBiasStep'] == True:
                superbias_step.output_dir = self.output_dir
                superbias_step.save_results = True
                if self.param['custBias'] == 'selfBias':
                    ## Have to save it manually if this step is skipped because of self bias subtraction
                    origName = saturation.meta.filename
                    if '_uncal.fits' in origName:
                        outName = origName.replace('_uncal.fits','_superbiasstep.fits')
                    else:
                        outName = 'cust_superbiasstep.fits'
                    
                    outPath = os.path.join(self.output_dir,outName)
                    saturation.to_fits(outPath,overwrite=True)
                
    
            # Call using the the output from the previously-run saturation step
            superbias = superbias_step.run(saturation)
            
            
            del saturation ## try to save memory
            
            
            if self.param['ROEBACorrection'] == True:
                # try using a copy of the bias results as the refpix output
                # refpix = refpix_step.run(superbias)
                # refpix_res = deepcopy(refpix)
                # the old way was to run the refpix and then replace it
                refpix_res = deepcopy(superbias)
    
                ngroups = superbias.meta.exposure.ngroups
                nints = superbias.data.shape[0] ## use the array size because segmented data could have fewer ints
                ## (instead of 
    
                # First, make sure that the aperture looks good. Here I have cheated and used a final rampfit result.
    
                # In[389]:
    
                if self.photParam is None:
                    phot = None
                else:
                    phot = phot_pipeline.phot(directParam=self.photParam)
    
    
                # In[390]:
    
    
                #phot.showStamps(showPlot=True,boxsize=200,vmin=0,vmax=1)
    
    
                # Everything inside the larger blue circle will be masked when doing reference pixel corrections
    
                # In[391]:
                
                
                for oneInt in tqdm.tqdm(np.arange(nints)):
                    for oneGroup in np.arange(ngroups):
            
                        rowSub, modelImg = rowamp_sub.do_backsub(superbias.data[oneInt,oneGroup,:,:],
                                                                 phot,amplifiers=self.param['noutputs'],
                                                                 backgMask=self.ROEBAmask,
                                                                 saveDiagnostics=self.param['saveROEBAdiagnostics'])
                        refpix_res.data[oneInt,oneGroup,:,:] = rowSub
                
                # # Linearity Step
    
                # In[328]:
                if self.param['saveROEBAdiagnostics'] == True:
                    origName = deepcopy(refpix_res.meta.filename)
                    if '.fits' in origName:
                        outName = origName.replace('.fits','_refpixstep.fits')
                    else:
                        outName = 'ROEBAstep.fits'
                
                    outPath = os.path.join(self.output_dir,outName)
                    refpix_res.to_fits(outPath,overwrite=True)
                
                
            else:
                # Instantiate and set parameters
                refpix_step = RefPixStep()
                #refpix_step.output_dir = output_dir
                #refpix_step.save_results = True
                refpix_res = refpix_step.run(superbias)
            
            del superbias ## try to save memory
    
            # Using the run() method
            linearity_step = LinearityStep()
            # refpix step
            linearity = linearity_step.run(refpix_res)
            del refpix_res ## try to save memory
    
            # # Persistence Step
    
            # Using the run() method
            persist_step = PersistenceStep()

            # skip for now since ref files are zeros
            persist_step.skip = True
    
            persist = persist_step.run(linearity)
            del linearity ## try to save memory
    
            # # Dark current step
    
            # Using the run() method
            dark_step = DarkCurrentStep()
    
            # There was a CRDS error so I'm skipping
            dark_step.skip = True
    
            # Call using the persistence instance from the previously-run
            # persistence step
            dark = dark_step.run(persist)

            del persist ## try to save memory
    
            # # Jump Step
    
            # In[335]:
    
    
            # Using the run() method
            jump_step = JumpStep()
            #jump_step.output_dir = output_dir
            #jump_step.save_results = True
            jump_step.rejection_threshold = self.param['jumpRejectionThreshold']

            jump_step.maximum_cores = self.max_cores
            # Call using the dark instance from the previously-run
            # dark current subtraction step
            jump = jump_step.run(dark)
            del dark ## try to save memory
    
            # # Ramp Fitting
    
            # In[344]:
    
    
            # Using the run() method
            ramp_fit_step = RampFitStep()
            ramp_fit_step.maximum_cores = self.max_cores
    
            ramp_fit_step.output_dir = self.output_dir
            ramp_fit_step.save_results = True
    
            # Let's save the optional outputs, in order
            # to help with visualization later
            #ramp_fit_step.save_opt = True
    
            # Call using the dark instance from the previously-run
            # jump step
            ramp_fit = ramp_fit_step.run(jump)
            del jump ## try to save memory
            del ramp_fit ## try to save memory
    
    
        executionTime = (time.time() - startTime)
        print('Stage 1 Execution Time in Seconds: ' + str(executionTime)) #Time how long this step takes
