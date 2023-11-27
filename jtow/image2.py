## Some tweaks to stage 2
from jwst.pipeline import Image2Pipeline
import os
import glob
import numpy as np

def run_image2(searchPaths):
    """
    quick fun to loop over many files
    """
    searchRes = np.sort(glob.glob(searchPaths))
    for oneFile in searchRes:
        proc_one_image(oneFile)

def predicted_output(oneFile):
    splitR = oneFile.split('_')
    splitR[-1] = 'cal.fits'
    return '_'.join(splitR)

        
def proc_one_image(oneFile):
    """
    Quick wrapper for image2 w/ some step tweaks
    """
    outputPredict = predicted_output(oneFile)
    if os.path.exists(outputPredict) == True:
        print("already found {}".format(outputPredict))
    else:
        res = Image2Pipeline.call(oneFile,
                                  save_results=True,
                                  steps={'resample': {'skip':True}})
    
