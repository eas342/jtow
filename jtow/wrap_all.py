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
defaultParamPath_jtow_nrc_SW = pkg_resources.resource_filename('jtow',
                                                               'params/default_jtow_nrc_short.yaml')
defaultParamPath_tshirt_spec = pkg_resources.resource_filename('jtow',
                                                               'params/default_tshirt_spec_params.yaml')

defaultParamPath_jtow_nrs_grating = pkg_resources.resource_filename('jtow',
                                                               'params/default_nrs_grating.yaml')

defaultParamPath_tshirt_nrs_grating = pkg_resources.resource_filename('jtow',
                                                               'params/default_tshirt_nrs_grating.yaml')

defaultParamPath_jtow_miri = pkg_resources.resource_filename('jtow',
                                                             'params/default_miri_lrs.yaml')
defaultParamPath_tshirt_lrs = pkg_resources.resource_filename('jtow',
                                                             'params/default_tshirt_miri_lrs.yaml')

tshirt_baseDir = phot_pipeline.get_baseDir()

class wrap(object):
    """
    Wrapper to run everything
    """

    def __init__(self,progID,obsNum,
                 recenteredNIRCamGrism=False,
                 crdsContext=None):
        """
        wrapper to run everything
        recenteredNIRCam grism allows for offset position
        """
        self.progID = progID
        self.obsNum = obsNum
        self.mast_path = os.environ['JWSTDOWNLOAD_OUTDIR']
        self.prog_dir = os.path.join(self.mast_path,"{:05d}".format(self.progID))
        self.obs_dir = os.path.join(self.prog_dir,"obsnum{:02d}".format(self.obsNum))
        self.recenteredNIRCamGrism = recenteredNIRCamGrism
        self.crdsContext = crdsContext

    def run_all(self):
        self.lookup_configuration()
        self.organize_files()
        self.make_miniseg()
        if self.instrument == 'NIRCAM':
            self.make_jtow_nrcalong()
            self.make_jtow_nrc_SW()
            self.run_jtow_nrcalong()
            self.make_tshirt_spec_param()
            self.run_jtow_nrc_SW()
            self.make_tshirt_phot_param()
        elif self.instrument == 'NIRSPEC':
            if self.grating == 'PRISM':
                self.make_jtow_prism()
            else:
                for detector in ['nrs1','nrs2']:
                    nrs_paramFile = self.make_jtow_nrs_grating(detector=detector)
                    jw = jtow.jw(nrs_paramFile)
                    jw.run_all()
                    tshirt_param =self.make_tshirt_spec_param(detector=detector)
                    spec = spec_pipeline.spec(tshirt_param)
                    spec.showStarChoices(showPlot=False)
                    spec.do_extraction(useMultiprocessing=True)
        elif self.instrument == 'MIRI':
            miri_paramFile = self.make_jtow_miri_lrs()
            jw = jtow.jw(miri_paramFile)
            jw.run_all()
            tshirt_param = self.make_tshirt_spec_param(detector='mirimage')
                
        #self.run_tshirt()

    def lookup_configuration(self):
        all_files = np.sort(glob.glob(os.path.join(self.obs_dir,'*.fits')))
        oneHead = fits.getheader(all_files[-1])
        self.instrument = oneHead['INSTRUME']
        if self.instrument == 'NIRSPEC':
            self.grating = oneHead['GRATING']
        
        
    def organize_files(self):
        if self.instrument == 'NIRCAM':
            ta_search = 'jw{:05d}{:03d}???_02102*'.format(self.progID,self.obsNum)
        else:
            ta_search = 'jw{:05d}{:03d}???_02101*'.format(self.progID,self.obsNum)
        
        ta_search_path = os.path.join(self.obs_dir,ta_search)
        ta_dir = os.path.join(self.obs_dir,'ta_files')
        move_or_link_files(ta_search_path,ta_dir,operation='link')

        if self.instrument == 'MIRI':
            dirimg_search = 'jw{:05d}{:03d}???_04102*'.format(self.progID,self.obsNum)
            dirimg_search_path = os.path.join(self.obs_dir,dirimg_search)
            dirimg_dir = os.path.join(self.obs_dir,'dirimg_files')
            move_or_link_files(dirimg_search_path,dirimg_dir,operation='link')
        else:
            dirimg_search_path = ''
        
        specFileTable = make_fileTable(os.path.join(self.obs_dir,'jw*'))
        unique_descriptors = np.unique(specFileTable['suffix'])
        for oneSuffix in unique_descriptors:
            descriptor_path = os.path.join(self.obs_dir,oneSuffix.replace('.','_'))
            fileSearch = os.path.join(self.obs_dir,'*{}'.format(oneSuffix))
            move_or_link_files(fileSearch,descriptor_path,
                               operation='link',
                               excludeSearch=ta_search_path,
                               excludeSearch2=dirimg_search_path)
            
    def make_miniseg(self,reDo=False):
        if self.instrument == 'NIRCAM':
            spec_uncal_dir = os.path.join(self.obs_dir,'nrcalong_uncal_fits')
            uncal_search = os.path.join(spec_uncal_dir,'*uncal.fits')
            self.LWdetSearchPath = uncal_search
            make_minisegments.loop_minisegments(uncal_search,reDo=reDo)
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
            make_minisegments.loop_minisegments(self.SWdetSearchPath,reDo=reDo)
        elif self.instrument == 'NIRSPEC':
            spec_uncal_dir = os.path.join(self.obs_dir,'nrs1_uncal_fits')
            uncal_search = os.path.join(spec_uncal_dir,'*uncal.fits')
            make_minisegments.loop_minisegments(uncal_search,reDo=reDo)
            if self.grating != 'PRISM':
                spec_uncal_dir2 = os.path.join(self.obs_dir,'nrs2_uncal_fits')
                uncal_search2 = os.path.join(spec_uncal_dir2,'*uncal.fits')
                make_minisegments.loop_minisegments(uncal_search2,reDo=reDo)
        elif self.instrument == 'MIRI':
            spec_uncal_dir = os.path.join(self.obs_dir,'mirimage_uncal_fits')
            uncal_search = os.path.join(spec_uncal_dir,'*uncal.fits')
            make_minisegments.loop_minisegments(uncal_search,reDo=reDo)
        else:
            raise NotImplementedError
    
    def get_SW_starPos(self,firstHead):
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
        return starPos
    
        
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
        starPos = self.get_SW_starPos(firstHead)
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
                                          'prog{}'.format(firstHead['PROGRAM']))
        if os.path.exists(tshirt_photDirPath) == False:
            os.makedirs(tshirt_photDirPath)
        tshirt_photName = "phot_param_{}_autoparam_001.yaml".format(photParams['nightName'])
        tshirt_photPath = os.path.join(tshirt_photDirPath,tshirt_photName)
        print("Writing photom auto parameter file to {}".format(tshirt_photPath))
        with open(tshirt_photPath,'w') as outFile:
            yaml.dump(photParams,outFile,default_flow_style=False)

    def make_tshirt_spec_param(self,detector='nrcalong'): 
        if (self.instrument == 'NIRCAM'):
            specParams = jtow.read_yaml(defaultParamPath_tshirt_spec)
            instrument_abbrev = 'nrc'
        elif (self.instrument == 'NIRSPEC'):
            if self.grating == 'PRISM':
                raise NotImplementedError
            else:
                specParams = jtow.read_yaml(defaultParamPath_tshirt_nrs_grating)
            instrument_abbrev = 'nrs'
        elif (self.instrument == 'MIRI'):
            specParams = jtow.read_yaml(defaultParamPath_tshirt_lrs)
            instrument_abbrev = 'mir'
        else:
            raise NotImplementedError

        spec_uncal_dir = os.path.join(self.obs_dir,'{}_uncal_fits'.format(detector))
        uncal_search = os.path.join(spec_uncal_dir,'*uncal.fits')
        first_lw_uncal = np.sort(glob.glob(uncal_search))[0]
        firstHead = fits.getheader(first_lw_uncal)
        
        specParams['procFiles'] = os.path.join(self.obs_dir,'{}_proc'.format(detector),
                                               'split_output',
                                               'ff_cleaned','*.fits')
        specParams['srcName'] = firstHead['TARGPROP']
        specParams['srcNameShort'] = "auto_params_001"
        srcFileName = specParams['srcName'].strip().replace(' ','_')
        
        if (self.instrument == 'NIRCAM'):
            if self.LWFilter == 'F444W':
                starPos = 31
                bkgRegionsY = [[5,21],[41,64]]
                dispPixels = [750,2040]
            elif self.LWFilter == 'F322W2':
                starPos = 34
                bkgRegionsY = [[5,24],[44,65]]
                dispPixels = [4,1747]
            else:
                raise NotImplementedError
            
            specParams['starPositions'] = [starPos]
            specParams['bkgRegionsY'] = bkgRegionsY
            filterDescrip = self.LWFilter
        elif (self.instrument == 'NIRSPEC'):
            if detector == 'nrs1':
                detectorGain = 1.420 ## NRS1 median
            else:
                detectorGain =  1.614 ## NRS2 median
            
            specParams['detectorGain'] = detectorGain
            if (firstHead['GRATING'] == 'G395H') & (firstHead['FILTER'] == 'F290LP'):
                if detector == 'nrs1':
                    dispPixels = [550,2044]
                else:
                    dispPixels = [4,2044]
            elif (firstHead['GRATING'] == 'G395M') & (firstHead['FILTER'] == 'F290LP'):
                if detector == 'nrs1':
                    dispPixels = [680,2044]
                else:
                    dispPixels = [None,None]
            else:
                raise NotImplementedError("Grating and filter not implemented.")

            filterDescrip = '{}_{}'.format(firstHead['GRATING'],detector)
        elif (self.instrument == 'MIRI'):
            filterDescrip = 'lrs'
            dispPixels = specParams['dispPixels'] ## just copy default, no filter change
        else:
            raise NotImplementedError
        
        specParams['nightName'] = "prog{}_{}_{}".format(firstHead['VISIT_ID'],srcFileName,filterDescrip)
        specParams['dispPixels'] = dispPixels
        
        tshirt_specDirPath = os.path.join(tshirt_baseDir,
                                          'parameters',
                                          'spec_params',
                                          'jwst',
                                          'prog_{}'.format(firstHead['PROGRAM']))
        if os.path.exists(tshirt_specDirPath) == False:
            os.makedirs(tshirt_specDirPath)
        
        tshirt_specName = "spec_{}_{}_autoparam_001.yaml".format(instrument_abbrev,
                                                                 specParams['nightName'])
        tshirt_specPath = os.path.join(tshirt_specDirPath,tshirt_specName)
        print("Writing spec auto parameter file to {}".format(tshirt_specPath))
        with open(tshirt_specPath,'w') as outFile:
            yaml.dump(specParams,outFile,default_flow_style=False,sort_keys=False)
        return tshirt_specPath

    def make_jtow_nrcalong(self):
        defaultParamPath = defaultParamPath_jtow_nrcalong
        jtow_paramName = self.make_jtow_spec(defaultParamPath,detName='nrcalong',
                                             recenteredNIRCamGrism=self.recenteredNIRCamGrism)
        self.jtow_nrcalong_paramfile = jtow_paramName

    def make_jtow_nrs_grating(self,detector='nrs1'):
        defaultParamPath = defaultParamPath_jtow_nrs_grating
        return self.make_jtow_spec(defaultParamPath,detName=detector)
    
    def make_jtow_miri_lrs(self):
        defaultParamPath = defaultParamPath_jtow_miri
        return self.make_jtow_spec(defaultParamPath,detName='mirimage')
            
    def make_jtow_spec(self,defaultParamPath,detName,
                       recenteredNIRCamGrism=False): 
        jtowParams = jtow.read_yaml(defaultParamPath)
        
        rawFileSearch = os.path.join(self.obs_dir,'{}_uncal_fits'.format(detName),
                                     'miniseg','*uncal.fits')
        jtowParams['rawFileSearch'] = rawFileSearch
        
        procFilePath = os.path.join(self.obs_dir,'{}_proc'.format(detName))
        jtowParams['outputDir'] = procFilePath
        first_spec_uncal = np.sort(glob.glob(rawFileSearch))[0]
        
        firstHead = fits.getheader(first_spec_uncal)

        srcFileName = firstHead['TARGPROP'].strip().replace(' ','_')
        
        if self.instrument == 'NIRSPEC':
            """
            Use rate files for mask
            """
            origFileName = firstHead['FILENAME']
            rate_file_use = os.path.join(self.obs_dir,origFileName.replace('uncal.fits','rate.fits'))
            jtowParams['ROEBAmaskfromRate'] = rate_file_use

        if "autoParamVersion" in jtowParams:
            autoParamVersion = jtowParams["autoParamVersion"]
        else:
            autoParamVersion = 1

        if recenteredNIRCamGrism == True:
            jtowParams['recenteredNIRCamGrism'] = recenteredNIRCamGrism
        
        jtowParams= self.set_crdsContext(jtowParams)

        jtow_paramName = "flight_{}_{}_{}_autoparam_{:03d}.yaml".format(firstHead['VISIT_ID'],
                                                                     detName,
                                                                     srcFileName,
                                                                     autoParamVersion)
        print("Writing photom auto parameter file to {}".format(jtow_paramName))
        with open(jtow_paramName,'w') as outFile:
            yaml.dump(jtowParams,outFile,default_flow_style=False)
        return jtow_paramName

    def run_jtow_nrcalong(self):
        jw = jtow.jw(self.jtow_nrcalong_paramfile)
        jw.run_all()
    
    def set_crdsContext(self,jtowParams):
        if self.crdsContext is not None:
            jtowParams['crdsContext'] = self.crdsContext
        return jtowParams

    def make_jtow_nrc_SW(self): 
        jtowParams = jtow.read_yaml(defaultParamPath_jtow_nrc_SW)
        
        rawFileSearch = os.path.join(self.obs_dir,self.SWdetSearch,
                                     'miniseg','*uncal.fits')
        jtowParams['rawFileSearch'] = rawFileSearch
        
        procFilePath = os.path.join(self.obs_dir,self.SWprocDir)
        jtowParams['outputDir'] = procFilePath
        first_lw_uncal = np.sort(glob.glob(rawFileSearch))[0]
        
        jtowParams= self.set_crdsContext(jtowParams)

        firstHead = fits.getheader(first_lw_uncal)

        srcFileName = firstHead['TARGPROP'].strip().replace(' ','_')
        detName = self.SWdetSearch.split('_')[0]
        jtow_paramName = "flight_{}_{}_{}_autoparam_001.yaml".format(firstHead['VISIT_ID'],
                                                                     detName,
                                                                     srcFileName)
        starPos = self.get_SW_starPos(firstHead)
        jtowParams['photParam']['refStarPos'] = [starPos]
        print("Writing photom auto parameter file to {}".format(jtow_paramName))
        self.jtow_SW_paramfile = jtow_paramName
        with open(jtow_paramName,'w') as outFile:
            yaml.dump(jtowParams,outFile,default_flow_style=False)
    
    def run_jtow_nrc_SW(self):
        jw = jtow.jw(self.jtow_SW_paramfile)
        jw.run_all()

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


def move_or_link_files(searchPath,destinationDir,excludeSearch='',
                       operation='link',
                       excludeSearch2=''):
    fileList = np.sort(glob.glob(searchPath))
    ensure_directory(destinationDir)

    excludeList1 = glob.glob(excludeSearch)
    excludeList2 = glob.glob(excludeSearch2)
    excludeList = excludeList1 + excludeList2
    for oneFile in fileList:
        baseName = os.path.basename(oneFile)
        outName = os.path.join(destinationDir,baseName)
        if (oneFile in excludeList) | (baseName in excludeList):
            excluded = True
        else:
            excluded = False
        
        if (os.path.exists(outName) == False) & (excluded == False):
            if operation == 'link':
                os.symlink(oneFile,outName)
            elif operation == 'move':
                os.rename(oneFile,outName)
            else:
                raise NotImplementedError
            
                
