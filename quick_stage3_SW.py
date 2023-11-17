from jwst.associations.asn_from_list import asn_from_list #Association file imports
import os
from jwst.pipeline.calwebb_image2 import Image2Pipeline #Stage 2
from jwst.pipeline.calwebb_tso3 import Tso3Pipeline #Stage 3
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
import glob
import numpy as np
from astropy.io import fits, ascii

all_calints_files = np.sort(glob.glob('*nrca3*calints.fits'))

outDir = '.'

if os.path.exists(outDir) == False:
    os.makedirs(outDir)


#Generate an association file required for stage 3
level3_asn = (os.path.join(outDir, 'nrca3_level3_asn.json')) #Name the stage 3 association file and give it a path.
asn_stage3 = asn_from_list(all_calints_files, product_name ='nrca3_level3_asn') #The rateints files; Name the output..
with open(level3_asn, 'w') as fh: #Write an association file.
       fh.write(asn_stage3.dump()[1])


#The file to use is the stage 3 association file defined above.

# Instantiate the class. Do not provide a configuration file.
pipeline_stage3 = Tso3Pipeline()

pipeline_stage3.outlier_detection.skip = True

# Specify that you want results saved to a file
pipeline_stage3.save_results = True
pipeline_stage3.output_dir = outDir

# Execute the pipeline using the run method
result_stage3 = pipeline_stage3.run(level3_asn)
