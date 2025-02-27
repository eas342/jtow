#!/usr/bin/env python
# coding: utf-8

from astroquery.mast import Observations
import http
import os

def get_by_obsnum(obsList,obsNumString):
    obsnumList = []
    for oneRow in obsList:
        obsID = oneRow['obs_id']
        obsnum = obsID[7:10]
        obsnumList.append(obsnum)
    obsList['obsnum'] = obsnumList
    pts_use = obsList['obsnum'] == obsNumString
    return pts_use

def do_download(propID=1185,obsNum=103,downloadAll=True,
                products=['RATE','UNCAL'],
                max_attempts=5):
    """
    Repeat the download process a few times in case of connection errors
    """
    attempt = 1
    while attempt <= max_attempts:
        try:
            res = do_download1(propID=propID,
                               obsNum=obsNum,
                               downloadAll=downloadAll,
                               products=products)
            attempt = max_attempts + 1
        except:
            print("Download failed. Trying again")
            attempt = attempt + 1
            res = []
    return res



def do_download1(propID=1185,obsNum=103,downloadAll=True,
                products=['RATE','UNCAL']):
    """
    Automatically download products from propID and obsNum (1 attempt)

    Parameters
    ----------
    propID: int
        Proposal ID
    obsNum: int
        Observation number
    downloadAll: bool
        Download all files? Otherwise, uses a subset for testing
    """


    my_session = Observations.login(token=os.environ['MAST_API_TOKEN'])
    obsNumString02 = "{:02d}".format(obsNum) ## for legacy directory organization
    obsNumString03 = "{:03d}".format(obsNum)
    propIDString = "{:05d}".format(propID)
    outPath = os.path.join(os.environ['JWSTDOWNLOAD_OUTDIR'],
                        propIDString,
                        "obsnum"+obsNumString02)
    if os.path.exists(outPath) == False:
        os.makedirs(outPath)

    ## For now, this is kind of slow so I will use Observations instead

    # missions = MastMissions(mission='jwst')
    # #my_session = missions.login(token=os.environ['MAST_API_TOKEN'])

    # # Use query_criteria method to use selected form search conditions for making missions_mast search API call
    # results = missions.query_criteria(select_cols=[
    # ],
    #     program='1281',
    #     observtn='1')

    # print(results)


    observation = Observations.query_criteria(proposal_id = propIDString)
    data_products = Observations.get_product_list(observation)

    # Filter them to get ramps and rateints; only for the fourth segment of data:
    downloadResList = []
    if downloadAll == True:
        download_products = products
    else:
        download_products = [products[0]]
    
    for oneProduct in download_products:
        files = Observations.filter_products(data_products, productType = 'SCIENCE', 
                                             productSubGroupDescription = oneProduct)
        pts_use = get_by_obsnum(files,obsNumString03)
        downloadRes = Observations.download_products(files[pts_use],flat=True,
                                                     download_dir=outPath)
        downloadResList.append(downloadRes)

    return downloadResList
