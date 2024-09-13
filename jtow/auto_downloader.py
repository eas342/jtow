#!/usr/bin/env python
# coding: utf-8

from astroquery.mast import Observations
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

def do_download(propID=1185,obsNum=103,downloadAll=True):
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
    rates = Observations.filter_products(data_products, productType = 'SCIENCE', productSubGroupDescription = 'RATE')
    uncal = Observations.filter_products(data_products, productType = 'SCIENCE', productSubGroupDescription = 'UNCAL')

    pts_use_rate = get_by_obsnum(rates,obsNumString03)
    pts_use_uncal = get_by_obsnum(uncal,obsNumString03)

    downloadRes = Observations.download_products(rates[pts_use_rate],flat=True,download_dir=outPath)
    
    if downloadAll == True:
        uncal_to_download = uncal[pts_use_uncal]
    else:
        ## for testing, only download first uncal (probly a TA image)
        uncal_to_download = uncal[pts_use_uncal][0]
    downloadRes2 = Observations.download_products(uncal_to_download,flat=True,download_dir=outPath)

    return [downloadRes,downloadRes2]
