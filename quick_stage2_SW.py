from jwst.associations.asn_from_list import asn_from_list #Association file imports
import os
from jwst.pipeline.calwebb_image2 import Image2Pipeline #Stage 2
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
import glob
import numpy as np
from astropy.io import fits, ascii

all_rateints_files = np.sort(glob.glob('proc_roeba/*nrca3_1_rampfitstep.fits'))

outDir = 'proc_roeba_stage2/'

asn_dir = outDir #Name the association file's directory.
level2_asn = (os.path.join(asn_dir, 'nrca3_level2_asn.json')) #Name the stage 2 association file and give it a path.
asn_stage2 = asn_from_list(all_rateints_files,rule=DMSLevel2bBase) #The rateints files; DMSLevel2bBase indicates that a Level2 association is to be created.
with open(level2_asn, 'w') as fh: #Write an association file.
    fh.write(asn_stage2.dump()[1])



# Instantiate the class. Do not provide a configuration file.
pipeline_stage2 = Image2Pipeline()

# Specify that you want results saved to a file
pipeline_stage2.save_results = True
pipeline_stage2.output_dir = outDir

# Execute the pipeline using the run method
result_stage2 = pipeline_stage2.run(level2_asn)


## correct the header to make sure the photometric aperture is centered

all_calints = np.sort(glob.glob(os.path.join(outDir,'*nrca3*calints.fits')))
for fitsName in all_calints: #Grabbing only nrca3 files from the directory
    print("Updating XREF in header for {}".format(fitsName))
    HDUList = fits.open(fitsName, 'update')
    HDUList[1].header['XREF_SCI'] = (2048. - 992.0, 'Aperture X reference point in SCI frame') #Fix x-position centering
    HDUList[1].header['YREF_SCI'] = (166.0, 'Aperture Y reference point in SCI frame') #Fix the y-position centering
    HDUList.close()
