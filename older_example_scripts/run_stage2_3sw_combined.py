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

fileList = np.sort(glob.glob('jw01059119*_nrc*_rate.fits'))
descrip = 'jw01059119_maskipr_pupil_sw_comb_asn'  
calSearch = 'jw01059119*_nrc??comb_cal.fits'

output_dir = '.'

for oneFile in fileList:
    stage2 = Image2Pipeline()
    stage2.flat_field.skip = True ## skip because we want to see the flat
    stage2.save_results = True
    stage2.output_dir = output_dir
    stage2.resample.skip = True
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
stage3.tweakreg.skip = True
stage3.skymatch.skip = True
stage3.save_results = True
stage3.run(asn_file)
