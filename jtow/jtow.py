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
from configparser import ConfigParser

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

import warnings

# In[359]:
import gc


## ES custom pipeline
from tshirt.pipeline import phot_pipeline
from tshirt.pipeline.instrument_specific import rowamp_sub
import tqdm
from splintegrate import splintegrate
from . import quick_ff_divide
from . import temporal_clean

path_to_defaults = "params/default_params.yaml"
defaultParamPath = pkg_resources.resource_filename('jtow',path_to_defaults)


def read_yaml(filePath):
    with open(filePath) as yamlFile:
        yamlStructure = yaml.safe_load(yamlFile)
    return yamlStructure

def log_output(TargetName):
    """
    Output the JWST pipeline log infromation to a seperate file.
    
    Parameters
    ---------
    TargetName: str
        Name of the target
    """
    config_object = ConfigParser()
    
    #Required sections for the configuration file
    config_object["*"] = {"handler": "file:{}_pipeline.log".format(TargetName), "level": "INFO"}
    
    #Write the above sections to stpipe-log.cfg file
    pwd = os.getcwd() #current working directory
    
    #Write the file to the working directory
    with open(pwd+'/stpipe-log.cfg'.format(TargetName), 'w') as conf:
        config_object.write(conf)
        
    print("A configuration file stpipe-log.cfg and a log output file {}_pipeline.log for the JWST Pipeline will be created in the working directory".format(TargetName, TargetName))

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
        
        ## check that there are no unexpected parameters
        for oneKey in self.param:
            if oneKey not in defaultParams.keys():
                warnings.warn("{} not an expected parameter".format(oneKey))
        
        self.set_up_dirs()
        self.get_files()
        
        self.make_descrip()
        
        self.max_cores = self.param['maxCores']
        
        self.make_roeba_masks()
        
        if self.param['custBias'] == 'cycleBias':
            self.check_biasCycle()
        
        if self.param['simpleSlopes'] == None:
            self.do_simple_ramp_fit = False
            self.do_full_ramp_fit = True
        elif self.param['simpleSlopes'] == 'Both':
            self.do_simple_ramp_fit = True
            self.do_full_ramp_fit = True
        elif self.param['simpleSlopes'] == 'Only':
            self.do_simple_ramp_fit = True
            self.do_full_ramp_fit = False
        else:
            raise Exception("Unrecognized simpleSlopes option {}. Options are None, Both and Only".format(simpelSlopes))
        
        
    
    def get_parameters(self,paramFile,directParam=None):
        if directParam is None:
            self.paramFile = paramFile
            self.param = read_yaml(paramFile)
        else:
            self.paramFile = 'direct dictionary'
            self.param = directParam
        
        if self.param['add_noutputs_keyword'] == True:
            warnings.warn("This code will modify the uncal file NOUTPUTS. This is DANGEROUS. Only use for older mirage simulations that lacked NOUTPUTS keyword")
        
    
    def check_biasCycle(self):
        self.biasCycleSearch = self.param['biasCycleSearch']
        searchResult = np.sort(glob.glob(self.biasCycleSearch))
        if len(searchResult) != 2:
            warnings.warn('Found something other than 2 files for the biasCycleSearch')
        
    
    def set_up_dirs(self):
        """
        Set up directories
        """
        self.output_dir = self.param["outputDir"]
        
        self.diagnostic_dir = os.path.join(self.param["outputDir"],'diagnostics')
        self.splitDir = os.path.join(self.param['outputDir'],'split_output')
        
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
            if self.param['add_noutputs_keyword'] == True:
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
            if self.param['forceHeaderChange'] is None:
                pass
            else:
                for oneKey in self.param['forceHeaderChange'].keys():
                    firstHead[oneKey] = self.param['forceHeaderChange'][oneKey]
                
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
            elif self.param['photParam'] != None:
                self.photParam = self.param['photParam']
                ROEBAmask = None
            elif firstHead['PUPIL'] == 'GRISMR':
                grismsFilterList = ['F322W2','F444W']
                if firstHead['FILTER'] in grismsFilterList:
                    self.photParam = None
                    mask1 = np.ones([Ny,Nx],dtype=bool)
                    
                    #mask1[0:4,:] = False
                    useSideRef = self.param['useGrismRefpx']
                    
                    if (useSideRef == True) & (firstHead['FILTER'] == 'F322W2'):
                        pass ## keep mask1[:,0:4] = True
                    else:
                        mask1[:,0:4] = False
                    
                    if (useSideRef == True) & (firstHead['FILTER'] == 'F444W'):
                        pass ## keep mask1[:,-4:] = True
                    else:
                        mask1[:,-4:] = False
                    
                    if firstHead['FILTER'] == 'F444W':
                        mask1[4:,637:2045] = False
                    elif firstHead['FILTER'] == 'F322W2':
                        mask1[4:,4:1846] = False
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
            grown_mask = self.grow_mask(ROEBAmask)
            good_rows = np.sum(np.sum(grown_mask,axis=1) >= 4)
            
            if good_rows != grown_mask.shape[0]:
                warnMessage = 'grown ROEBA mask has too few rows to fit for {}. Skipping the growth'.format(self.descrip)
                print(warnMessage)
                warnings.warn(warnMessage)
                self.ROEBAmask = ROEBAmask
            else:
                self.ROEBAmask = grown_mask
        
        if self.ROEBAmask is None:
            pass
        else:
            self.good_rows = np.sum(np.sum(self.ROEBAmask,axis=1) >= 4)
            
            if self.good_rows != self.ROEBAmask.shape[0]:
                warnMessage = 'final ROEBA mask has too few rows to fit for {}. Setting it to None and turning off ROEBA'.format(self.descrip)
                print(warnMessage)
                warnings.warn(warnMessage)
                self.ROEBAmask = None
                self.param['ROEBACorrection'] = False
        
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
        ## but don't grow the bad pixels or isolated pixels
        border = np.array([[1,1,1],[1,0,1],[1,1,1]])
        has_neighbors = ndimage.convolve(np.array(img == 0),border,mode='constant',cval=0.0)
        
        arr_to_convolve = (img == 0) & (self.bad_dq_mask == False) & has_neighbors
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
            
    def lineInterceptBias(self,stepResult):
        """
        Fit a line to all pixels and use the intercept
        """
        data = stepResult.data
        
        result = deepcopy(data)
        
        for ind,oneInt in enumerate(data):
            nz, ny, nx = oneInt.shape
            x = np.arange(nz)
    
            flatDat = np.reshape(oneInt,[nz,nx * ny])
            pfit = np.polyfit(x,flatDat,1)
            intercept2D = np.reshape(pfit[1],[ny,nx])
            slope2D = np.reshape(pfit[0],[ny,nx])
            result[ind] = oneInt - intercept2D
        
        return result 
    
        
    def cycleBiasSub(self,stepResult):
        """
        Cycle through the bias pattern defined by biasCycle
        """
        int_start = stepResult.meta.exposure.integration_start
        cycleLen = len(self.param['biasCycle'])
        cycler_counter = np.mod(0 + int_start - 1,cycleLen)
        ngroups = stepResult.meta.exposure.ngroups
        data = stepResult.data
        
        result = deepcopy(data)
        
        for ind,oneInt in enumerate(data):
            biasType = self.param['biasCycle'][cycler_counter]
            biasPath = self.biasCycleSearch.replace('?',biasType)
            
            
            datBias = fits.getdata(biasPath)
            tiledBias = np.tile(datBias,[ngroups,1,1])
            
            result[ind] = oneInt - tiledBias
            
            cycler_counter = np.mod(cycler_counter + 1,cycleLen)
            
        return result
    
    def simple_ramp_fit(self,ramp4D,uncal_name):
        nint, ngroup, ny, nx = ramp4D.shape
        ramp4D.meta.exposure
        
        tgroup = ramp4D.meta.exposure.group_time
        x = np.arange(ngroup) * tgroup
        
        rampfit3D_slope = np.zeros([nint,ny,nx])
        rampfit3D_intercept = np.zeros([nint,ny,nx])
        for intNum in tqdm.tqdm(np.arange(nint)):
            oneInt = ramp4D.data[intNum]
            
            flatDat = np.reshape(oneInt,[ngroup,nx * ny])
            pfit = np.polyfit(x,flatDat,1)
            intercept2D = np.reshape(pfit[1],[ny,nx])
            slope2D = np.reshape(pfit[0],[ny,nx])
            rampfit3D_slope[intNum] = slope2D
            rampfit3D_intercept[intNum] = intercept2D
        
        ## Save the result
        origHead = fits.getheader(uncal_name)
        
        ## save the DQ
        if hasattr(ramp4D,'groupdq'):
            saveDQ = True
            DQ_res = np.bitwise_or.reduce(ramp4D.groupdq,axis=1)
            dqHDU = fits.ImageHDU(DQ_res)
            dqHDU.name = 'DQ'
        else:
            saveDQ = False
        
        for oneOutput in ['slopes','intercepts']:
            if oneOutput == 'slopes':
                outSuffix = '_simple_slopes.fits'
                result3D = rampfit3D_slope
            else:
                outSuffix = '_simple_intercepts.fits'
                result3D = rampfit3D_intercept
            
            outName_result = ramp4D.meta.filename.replace('.fits',outSuffix).replace('_uncal','')
            outPath_result = os.path.join(self.output_dir,outName_result)
            primHDU = fits.PrimaryHDU(None,origHead)
            resultHDU = fits.ImageHDU(result3D)
            resultHDU.name = "SCI"
            HDUList_result = fits.HDUList([primHDU,resultHDU])
            ## Also save the DQ
            if saveDQ == True:
                HDUList_result.append(dqHDU)
            
            print("Saving result to {}".format(outPath_result))
            HDUList_result.writeto(outPath_result,overwrite=True)
            
            del HDUList_result
        
    
    def run_roeba(self,superbias):
        """
        Do the ROEBA (row-by-row, odd/even by amplifier algorithm)
        
        Parameters
        -----------
        superbias: jwst cube model (I think)
            The result of the superbias subtraction to run w/ ROEBA
        
        
        """
        # try using a copy of the bias results as the refpix output
        # refpix = refpix_step.run(superbias)
        # refpix_res = deepcopy(refpix)
        # the old way was to run the refpix and then replace it
        refpix_res = deepcopy(superbias)
        
        
        ## (instead of 

        # First, make sure that the aperture looks good. Here I have cheated and used a final rampfit result.

        # In[389]:

        if self.photParam is None:
            phot = None
        else:
            phot = phot_pipeline.phot(directParam=self.photParam)


        # In[390]:

        nints,ngroups,ny,nx = superbias.data.shape
        #phot.showStamps(showPlot=True,boxsize=200,vmin=0,vmax=1)


        # Everything inside the larger blue circle will be masked when doing reference pixel corrections

        # In[391]:
        
        
        ## have to transpose NIRISS and NIRSpec, but not NIRCam
        if superbias.meta.instrument.name == "NIRCAM":
            transposeForROEBA = False
            backgMask = self.ROEBAmask
        else:
            transposeForROEBA = True
            if self.ROEBAmask is None:
                backgMask = self.ROEBAmask
            else:
                backgMask = self.ROEBAmask.T
        
        for oneInt in tqdm.tqdm(np.arange(nints)):
            if self.param['ROEBAK'] == True:
                iterations = 2
                fastRead = [False,True]
                intermediate_result = np.zeros([ngroups,ny,nx])
                slowReadModelCube = np.zeros([ngroups,ny,nx])
            else:
                iterations=1
                fastRead = [True]
            
            roeba_one_int = np.zeros([ngroups,ny,nx])
            
            for oneIteration in np.arange(iterations):
                doFastRead = fastRead[oneIteration]
                
                if (self.param['ROEBAK'] == True) & (oneIteration == 2-1):
                    kTCadjustment = np.median(intermediate_result,axis=0)
                    cubeToCorrect = intermediate_result - kTCadjustment
                    
                    if self.param['saveROEBAdiagnostics'] == True:
                        origName = deepcopy(refpix_res.meta.filename)
                        if '.fits' in origName:
                            outName = origName.replace('.fits','_kTC_int_{:04d}.fits'.format(oneInt))
                        else:
                            outName = 'ROEBA_kTC_int_{:04d}.fits'.format(oneInt)
        
                        outPath = os.path.join(self.output_dir,'diagnostics',outName)
                        HDUList = fits.PrimaryHDU(kTCadjustment)
                        HDUList.writeto(outPath,overwrite=True)
                    
                else:
                    cubeToCorrect = superbias.data[oneInt]
                
                for oneGroup in np.arange(ngroups):
                    
                    if transposeForROEBA == True:
                        imgToCorrect = cubeToCorrect[oneGroup,:,:].T
                    else:
                        imgToCorrect = cubeToCorrect[oneGroup,:,:]
                    
                    if self.param['ROEBACorrection'] == 'GROEBA':
                        GROEBA = True
                    else:
                        GROEBA = False
                    
                    rowSub, slowImg, fastImg = rowamp_sub.do_backsub(imgToCorrect,
                                                             phot,amplifiers=self.param['noutputs'],
                                                             backgMask=backgMask,
                                                             saveDiagnostics=self.param['saveROEBAdiagnostics'],
                                                             returnFastSlow=True,
                                                             colByCol=self.param['colByCol'],
                                                             smoothSlowDir=self.param['smoothSlowDir'],
                                                             GROEBA=GROEBA)
                    
                    if (self.param['ROEBAK'] == True) & (oneIteration == 2-1):
                        groupResult = intermediate_result[oneGroup] - fastImg
                    else:
                        groupResult = rowSub
                    
                    ## Save on the last iteration
                    if oneIteration == iterations - 1:
                        if transposeForROEBA == True:
                            roeba_one_int[oneGroup,:,:] = groupResult.T
                        else:
                            roeba_one_int[oneGroup,:,:] = groupResult
                    else:
                        intermediate_result[oneGroup,:,:] = rowSub
                        slowReadModelCube[oneGroup,:,:] = slowImg
            
            refpix_res.data[oneInt,:,:,:] = roeba_one_int
        
        
        # In[328]:
        if self.param['saveROEBAdiagnostics'] == True:
            origName = deepcopy(refpix_res.meta.filename)
            if '.fits' in origName:
                outName = origName.replace('.fits','_refpixstep.fits')
            else:
                outName = 'ROEBAstep.fits'
        
            outPath = os.path.join(self.output_dir,outName)
            refpix_res.to_fits(outPath,overwrite=True)
        return refpix_res
    
    def delete_object(self,obj,step=None):
        """
        Try to delete an object's data to free up memory
        """
        # if hasattr(obj,'data'):
        #     del obj.data
        # if hasattr(obj,'groupdq'):
        #     del obj.groupdq
        # if hasattr(obj,'err'):
        #     del obj.err
        # if hasattr(obj,'dq'):
        #     del obj.dq
        # if hasattr(obj,'refout'):
        #     del obj.refout
        #
        if step is None:
            del obj
        else:
            if step.skip == True:
                pass
            else:
                del obj
        
        gc.collect()
        pass
    
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
            self.delete_object(dq_init) ## try to save memory
    
            # Using the run() method
            superbias_step = SuperBiasStep()
            
            if self.param['custBias'] is None:
                pass
            elif self.param['custBias'] == 'selfBias':
                superbias_step.skip = True
                saturation.data = saturation.data - saturation.data[0][0]
            elif self.param['custBias'] == 'cycleBias':
                superbias_step.skip = True
                saturation.data = self.cycleBiasSub(saturation)
            elif self.param['custBias'] == 'lineIntercept':
                superbias_step.skip = True
                saturation.data = self.lineInterceptBias(saturation)
            else:
                superbias_step.override_superbias = self.param['custBias']
            
            if self.param['saveBiasStep'] == True:
                superbias_step.output_dir = self.output_dir
                superbias_step.save_results = True
                if self.param['custBias'] in ['selfBias','cycleBias','lineIntercept']:
                    ## Have to save it manually if this step is skipped because of self bias subtraction
                    origName = deepcopy(saturation.meta.filename)
                    if '_uncal.fits' in origName:
                        outName = origName.replace('_uncal.fits','_superbiasstep.fits')
                    else:
                        outName = 'cust_superbiasstep.fits'
                    
                    outPath = os.path.join(self.output_dir,outName)
                    saturation.to_fits(outPath,overwrite=True)
                    
                    ## try to return filename back to original
                    saturation.meta.filename = origName
    
            # Call using the the output from the previously-run saturation step
            superbias = superbias_step.run(saturation)
            
            ngroups = superbias.meta.exposure.ngroups
            nints = superbias.data.shape[0] ## use the array size because segmented data could have fewer ints
            
            if self.param['custGroupDQfile'] is not None:
                custGroupDQ = fits.getdata(self.param['custGroupDQfile'])
                tiled_custGroup = np.tile(custGroupDQ,[nints,1,1,1])
                superbias.groupdq = (superbias.groupdq | tiled_custGroup)
                
            self.delete_object(saturation,step=superbias_step) ## try to save memory
            
            
            if (self.param['ROEBACorrection'] == True) | (self.param['ROEBACorrection'] == "GROEBA"):
                refpix_res = self.run_roeba(superbias)
                
                
                
            else:
                # Instantiate and set parameters
                refpix_step = RefPixStep()
                refpix_step.output_dir = self.output_dir
                if self.param['saveROEBAdiagnostics'] == True:
                    refpix_step.save_results = True
                
                refpix_step.side_smoothing_length=self.param['side_smoothing_length']
                refpix_res = refpix_step.run(superbias)
            
                self.delete_object(superbias,step=refpix_step) ## try to save memory
            
            # # Linearity Step   
            # Using the run() method
            linearity_step = LinearityStep()
            
            if self.param['doLincor'] == True:
                linearity_step.skip = False
            elif self.param['doLincor'] == False:
                linearity_step.skip = True
            else:
                raise Exception("Unrecognized doLincor value {}".format(self.param['doLinearity']))
                
            
            linearity = linearity_step.run(refpix_res)
            
            self.delete_object(refpix_res,step=linearity_step)
            
            # # Persistence Step
    
            # Using the run() method
            #persist_step = PersistenceStep()
            #
            ## skip for now since ref files are zeros
            #persist_step.skip = True
            #
            #persist = persist_step.run(linearity)
            #self.delete_object(linearity) ## try to save memory
    
            # # Dark current step
    
            # Using the run() method
            #dark_step = DarkCurrentStep()
            #
            # There was a CRDS error so I'm skipping
            #dark_step.skip = True
    
            # Call using the persistence instance from the previously-run
            # persistence step
            #dark = dark_step.run(persist)
            
            #self.delete_object(persist)

            ## to save memory, just move on without running a step with skip=True
            dark_result = linearity
            # # Jump Step
            # In[335]:
    
    
            # Using the run() method
            jump_step = JumpStep()
            #jump_step.output_dir = output_dir
            #jump_step.save_results = True
            jump_step.rejection_threshold = self.param['jumpRejectionThreshold']
            jump_step.skip = self.param['skipJumpDet']
            
            jump_step.maximum_cores = self.max_cores
            
            jump_step.flag_4_neighbors = self.param['jumpStepFlag4Neighbors']
            
            if self.param['saveJumpStep'] == True:
                jump_step.output_dir = self.output_dir
                jump_step.save_results = True
            else:
                pass
            
            # Call using the dark instance from the previously-run
            # dark current subtraction step
            jump = jump_step.run(dark_result)
            
            self.delete_object(dark_result,step=jump_step)
            
            # # Ramp Fitting
    
            # In[344]:
    
            ## Do a simple ramp fit if parameters are set
            if self.do_simple_ramp_fit == True:
                self.simple_ramp_fit(jump,uncal_file)
            
            if self.do_full_ramp_fit == True:
                # Using the run() method
                ramp_fit_step = RampFitStep()
                ramp_fit_step.weighting = self.param['rampFitWeighting']
                
                ramp_fit_step.maximum_cores = self.max_cores
    
                ramp_fit_step.output_dir = self.output_dir
                ramp_fit_step.save_results = True
                
                if hasattr(ramp_fit_step,'suppress_one_group'):
                    ramp_fit_step.suppress_one_group = self.param['suppressOneGroup']
                
                # Let's save the optional outputs, in order
                # to help with visualization later
                #ramp_fit_step.save_opt = True
            

            
                # Call using the dark instance from the previously-run
                # jump step
                ramp_fit0, ramp_fit1 = ramp_fit_step.run(jump)
                
                self.delete_object(ramp_fit0)
                self.delete_object(ramp_fit1)
            
            self.delete_object(jump) ## try to save memory
            
    
    
        executionTime = (time.time() - startTime)
        print('Stage 1 Execution Time in Seconds: ' + str(executionTime)) #Time how long this step takes

    def splintegrate(self):
        """
        Split up the rateints files into individual ones
        """
        filesToSplit = os.path.join(self.param['outputDir'],'*1_rampfitstep.fits')
        splintegrate.run_on_multiple(inFiles=filesToSplit,
                                    outDir=self.splitDir,overWrite=True,
                                    detectorName=None,
                                    flipToDet=False,
                                    mirageSeedFile=False)
        
    def do_flat(self):
        """
        Flat field the split up (ie. splintegrated) rateints results
        """
        filesToFlat = os.path.join(self.splitDir,'*.fits')
        quick_ff_divide.quick_ff_divide(filesToFlat)
    
    def temporal_clean(self):
        """
        Do a sigma clipping clean
        """
        filesToClean = os.path.join(self.splitDir,'ff_div','*.fits')
        cleanedPath = os.path.join(self.splitDir,'ff_cleaned')
        temporal_clean.clean_data(searchPath=filesToClean,
                                  savePath=cleanedPath)
    
    def run_all(self):
        self.run_jw()
        self.splintegrate()
        self.do_flat()
    
