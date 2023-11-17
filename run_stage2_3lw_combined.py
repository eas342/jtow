import asdf
import copy
import os
import shutil

# Numpy library:
import numpy as np

# For downloading data
import requests
import glob

from jwst.associations.asn_from_list import asn_from_list
from jwst.pipeline.calwebb_image3 import Image3Pipeline

from jwst.associations.lib.rules_level3_base import DMS_Level3_Base

from jwst.pipeline.calwebb_image2 import Image2Pipeline

prefix = 'jw01160022001_02103'
fileList = np.sort(glob.glob('{}*nrc?long_0_rampfitstep.fits'.format(prefix)))
descrip = '{}_lw_f410m_asn'.format(prefix)
calSearch = '{}*nrc?long*0_cal.fits'.format(prefix)

output_dir = '.'

for oneFile in fileList:
    stage2 = Image2Pipeline()
    stage2.flat_field.skip = False
    stage2.save_results = True
    stage2.output_dir = output_dir
    stage2.resample.skip = True ## this only a diagnostic image
    stage2.run(oneFile)
    del stage2


cal_files = np.sort(glob.glob(calSearch))

#cal_files = np.sort(glob.glob('jw01059119001_02104_00001_nrc?long_cal.fits'))
#descrip = 'jw01059119_001_02104lw_asn'  

asn_file = descrip + '.json'

stage3_asn = asn_from_list(cal_files,product_name=descrip)


with open(asn_file,'w') as fh:
    name, serialized = stage3_asn.dump(format='json')
    fh.write(serialized)

stage3 = Image3Pipeline()
stage3.tweakreg.skip = False
stage3.skymatch.skip = True ## sky is subtracted by ROEBA
stage3.save_results = True
stage3.run(asn_file)
