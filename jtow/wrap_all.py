from . import jtow
from . import make_minisegments
import os
import glob
import yaml
import numpy as np
from astropy.table import Table
from astropy.io import fits, ascii
import pdb
from tshirt.pipeline import phot_pipeline, spec_pipeline
import pkg_resources

path_to_defaults_tshirt_phot = "params/default_tshirt_phot_params.yaml"
defaultParamPath_tshirt_phot = pkg_resources.resource_filename('jtow',path_to_defaults_tshirt_phot)
defaultParamPath_jtow_nrcalong = pkg_resources.resource_filename('jtow',
                                                                 'params/default_jtow_nrcalong.yaml')

tshirt_baseDir = phot_pipeline.get_baseDir()

class wrap(object):
    """
    Wrapper to run everything
    """

    def __init__(self,progID,obsNum):
        """
        wrapper to run everything
        """
        self.progID = progID
        self.obsNum = obsNum
        self.mast_path = os.environ['JWSTDOWNLOAD_OUTDIR']
        self.prog_dir = os.path.join(self.mast_path,"{:05d}".format(self.progID))
        self.obs_dir = os.path.join(self.prog_dir,"obsnum{:02d}".format(self.obsNum))

    def run_all(self):
        self.organize_files()
        self.make_miniseg()
        #self.make_jtow_param()
        #self.run_jtow()
        self.make_tshirt_phot_param()
        #self.make_tshirt_spec_param()
        #self.run_tshirt()

    def organize_files(self):
        ta_search = 'jw{:05d}{:03d}???_02102*'.format(self.progID,self.obsNum)
        ta_search_path = os.path.join(self.obs_dir,ta_search)
        ta_dir = os.path.join(self.obs_dir,'ta_files')
        
        move_files(ta_search_path,ta_dir)
        specFileTable = make_fileTable(os.path.join(self.obs_dir,'jw*'))
        unique_descriptors = np.unique(specFileTable['suffix'])
        for oneSuffix in unique_descriptors:
            descriptor_path = os.path.join(self.obs_dir,oneSuffix.replace('.','_'))
            fileSearch = os.path.join(self.obs_dir,'*{}'.format(oneSuffix))
            move_files(fileSearch,descriptor_path)
                                      
    def make_miniseg(self):
        spec_uncal_dir = os.path.join(self.obs_dir,'nrcalong_uncal_fits')
        uncal_search = os.path.join(spec_uncal_dir,'*uncal.fits')
        make_minisegments.loop_minisegments(uncal_search)
        first_uncal = np.sort(glob.glob(uncal_search))[0]
        firstHead = fits.getheader(first_uncal)
        self.LWFilter = firstHead['FILTER']
        self.LWPupil = firstHead['PUPIL']
        if firstHead['FILTER'] == 'F444W':
            self.SWdetSearch = 'nrca1_uncal_fits'
        else:
            self.SWdetSearch = 'nrca3_uncal_fits'

        self.SWprocDir = self.SWdetSearch.replace('uncal_fits','proc')
            
        self.SWdetSearchPath = os.path.join(self.obs_dir,self.SWdetSearch,'*uncal.fits')
        make_minisegments.loop_minisegments(self.SWdetSearchPath)
        
    def make_tshirt_phot_param(self): 
        photParams = jtow.read_yaml(defaultParamPath_tshirt_phot)
        
        first_sw_uncal = np.sort(glob.glob(self.SWdetSearchPath))[0]
        firstHead = fits.getheader(first_sw_uncal)
        self.SWFilter = firstHead['FILTER']
        self.SWPupil = firstHead['PUPIL']
        
        photParams['procFiles'] = os.path.join(self.obs_dir,self.SWprocDir,
                                               'split_output',
                                               'ff_cleaned','*.fits')
        photParams['srcName'] = firstHead['TARGPROP']
        photParams['srcNameShort'] = "auto_params_001"
        srcFileName = photParams['srcName'].strip().replace(' ','_')
        photParams['nightName'] = "prog{}_{}_{}".format(firstHead['VISIT_ID'],srcFileName,self.LWFilter)
        if self.LWFilter == 'F444W':
            if (firstHead['SUBARRAY'] == 'SUBGRISM256') | (firstHead['SUBARRAY'] == 'FULL'):
                starPos = [1794.27,161.54]
            elif firstHead['SUBARRAY'] == 'SUBGRISM64':
                starPos = [1796.0,35.5]
            else:
                raise NotImplementedError
        else:
            if (firstHead['SUBARRAY'] == 'SUBGRISM256') | (firstHead['SUBARRAY'] == 'FULL'):
                starPos = [1060.7, 165.9]
            elif firstHead['SUBARRAY'] == 'SUBGRISM64':
                starPos = [1064.4, 30.9]
            else:
                raise NotImplementedError
        photParams['refStarPos'] = [starPos]
        if self.SWPupil == 'WLP8':
            apertures = [79,79,100]
        elif self.SWPupil == 'WLP4':
            apertures = [31.5,32,60]
        else:
            raise NotImplementedError
        photParams['apRadius'] = apertures[0]
        photParams['backStart'] = apertures[1]
        photParams['backEnd'] = apertures[2]
        
        tshirt_photDirPath = os.path.join(tshirt_baseDir,
                                          'parameters',
                                          'phot_params',
                                          'jwst_flight_data',
                                          'prog01185'.format(firstHead['PROGRAM']))
        tshirt_photName = "phot_param_{}_autoparam_001.yaml".format(photParams['nightName'])
        tshirt_photPath = os.path.join(tshirt_photDirPath,tshirt_photName)
        print("Writing photom auto parameter file to {}".format(tshirt_photPath))
        with open(tshirt_photPath,'w') as outFile:
            yaml.dump(photParams,outFile,default_flow_style=False)


    def make_jtow_nrcalong(self): 
        jtowParams = jtow.read_yaml(defaultParamPath_jtow_nrcalong)

        rawFileSearch = os.path.join(self.obs_dir,'nrcalong_uncal_fits',
                                     'miniseg','*uncal.fits')
        jtowParams['rawFileSearch'] = rawFileSearch
        
        procFilePath = os.path.join(self.obs_dir,'nrcalong_proc')
        jtowParams['outputDir'] = procFilePath
        first_lw_uncal = np.sort(glob.glob(rawFileSearch))[0]
        
        firstHead = fits.getheader(first_lw_uncal)

        srcFileName = firstHead['TARGPROP'].strip().replace(' ','_')
        
        jtow_paramName = "flight_{}_nrcalong_{}_autoparam_001.yaml".format(firstHead['VISIT_ID'],
                                                                           srcFileName)
        
        print("Writing photom auto parameter file to {}".format(jtow_paramName))
        with open(jtow_paramName,'w') as outFile:
            yaml.dump(jtowParams,outFile,default_flow_style=False)
        

def make_fileTable(searchPath):
    t = Table()
    fileList = np.sort(glob.glob(searchPath))
    detectorDescriptors = []
    file_suffixes = []
    for oneFile in fileList:
        detector_and_descrip = oneFile.split('_')[-2:]
        file_suffix = '_'.join(detector_and_descrip)
        detectorDescriptor = file_suffix.replace('.','_')
        detectorDescriptors.append(detectorDescriptor)
        file_suffixes.append(file_suffix)
    t['name'] = fileList
    t['det+descrip'] = detectorDescriptors
    t['suffix'] = file_suffixes
    return t
        
def ensure_directory(path):
    if os.path.exists(path) == False:
        os.makedirs(path)

def move_files(searchPath,destinationDir):
    fileList = np.sort(glob.glob(searchPath))
    ensure_directory(destinationDir)
    for oneFile in fileList:
        baseName = os.path.basename(oneFile)
        outName = os.path.join(destinationDir,baseName)
        os.rename(oneFile,outName)


    
