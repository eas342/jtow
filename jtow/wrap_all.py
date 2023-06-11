from . import jtow
from . import make_minisegments
import os
import glob
import numpy as np
from astropy.table import Table
import pdb

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

    def run_all():
        self.organize_files()
        self.make_miniseg()
        self.make_jtow_param()
        self.run_jtow()
        self.make_tshirt_param()
        self.run_tshirt()

    def organize_files(self):
        self.prog_dir = os.path.join(self.mast_path,"{:05d}".format(self.progID))
        self.obs_dir = os.path.join(self.prog_dir,"obsnum{:02d}".format(self.obsNum))
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


    
