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
#from scipy import ndimage
from copy import deepcopy
#import pkg_resources
import os
import glob
import crds
import tqdm
from astropy.table import Table


class custGroupDQ(object):
    def __init__(self,rateFile=None,
                 fracSat=0.9,custSatVal=None):
        """
        Make an object to make a custom Group DQ

        The original intention is to create a column saturated DQ
        for TSO observations so that
            1) the number of included pixels is constant
            and doesn't switch depending on the flux level or noise
            above/below the stauration threshold
            2) the background and sky pixels use the same number of groups

        Parameters
        -----------
        rateFile: str
            Path to rate file
        fracSat: float
            Fraction of the saturation level
        custSatVal: int or float
            Custom saturation value. Otherwise, the saturation reference
            file is used
        """
        self.rateFile = rateFile
        self.head = fits.getheader(self.rateFile)
        self.sciHead = fits.getheader(self.rateFile,extname='SCI')
        self.get_ref_files()
        self.fracSat = fracSat
        self.custSatVal = custSatVal

        self.outName = self.rateFile.replace('_rate.fits','_rate_groupdq.fits')
        self.outPlotName = self.rateFile.replace('_rate.fits','_groupdq_image.pdf')

    def get_ref_files(self):
        """
        Get the default saturation file
        """
        reftypes = ['saturation','superbias']
        self.refs = crds.getreferences(self.head,reftypes=reftypes)
    
    def calc_satpoints(self,useBias=True):
        """
        Calculate the points at which a pixel saturates
        """
        self.rate = fits.getdata(self.rateFile)

        groupTime = self.head['TGROUP']
        ngroup = self.head['NGROUPS']
        groupArr = np.arange(ngroup)
        ny = self.sciHead['NAXIS2']
        nx = self.sciHead['NAXIS1']
        groupNum3D = np.transpose(np.tile(groupArr,[nx,ny,1]),[2,1,0])

        if useBias == True:
            bias = fits.getdata(self.refs['superbias'])
            biasHead = fits.getheader(self.refs['superbias'])

            xStart = self.head['SUBSTRT1'] - biasHead['SUBSTRT1']
            xEnd = xStart + self.head['SUBSIZE1']
            yStart = self.head['SUBSTRT2'] - biasHead['SUBSTRT2']
            yEnd = yStart + self.head['SUBSIZE2']
            useBias = bias[yStart:yEnd,xStart:xEnd]
            
        else:
            useBias = 0.

        if self.custSatVal == True:
            useSat = self.custSatVal
        else:
            sat = fits.getdata(self.refs['saturation'])
            satHead = fits.getheader(self.refs['saturation'])

            xStart2 = self.head['SUBSTRT1'] - satHead['SUBSTRT1']
            xEnd2 = xStart2 + self.head['SUBSIZE1']
            yStart2 = self.head['SUBSTRT2'] - satHead['SUBSTRT2']
            yEnd2 = yStart2 + self.head['SUBSIZE2']
            useSat = sat[yStart2:yEnd2,xStart2:xEnd2]

        self.calcCube = groupNum3D * self.rate * groupTime + useBias
        self.isSat = self.calcCube > self.fracSat * useSat
        satInColumn = np.sum(self.isSat > 0,axis=1) > 0
        tmpTile = np.tile(satInColumn,[ny,1,1])
        self.declareSat = np.transpose(tmpTile,[1,0,2])
        
        custGroupDQ = 2**1 * self.declareSat
        groupDQHDU = fits.PrimaryHDU(custGroupDQ)
        copyTerms = ['SUBSTRT1','SUBSTRT2','SUBSIZE1','SUBSIZE2']
        for oneTerm in copyTerms:
            groupDQHDU.header[oneTerm] = self.head[oneTerm]
        groupDQHDU.writeto(self.outName,overwrite=True)

        fig, ax = plt.subplots()
        ax.imshow(np.sum(self.declareSat,axis=0))
        fig.savefig(self.outPlotName)

        

