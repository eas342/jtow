import crds
from astropy.io import fits, ascii
import os


example_file = os.path.join(os.environ['JWSTDOWNLOAD_OUTDIR'],
                            '01185/uncal_obs_017_gj3470',
                            'jw01185017001_02102_00001-seg001_nrcalong_uncal.fits')


def compare_recs(file1=example_file,
                 context1=1093,context2=1137):
    
    pmap1 = "jwst_{}.pmap".format(context1)
    pmap2 = "jwst_{}.pmap".format(context2)
    

    head = fits.getheader(file1)
    rec1 = crds.getrecommendations(head,context=pmap1)
    rec2 = crds.getrecommendations(head,context=pmap2)

    counter = 0
    for oneKey in rec1.keys():
        if rec1[oneKey] == rec2[oneKey]:
            pass
        else:
            print("Difference for {}. {} vs {}".format(oneKey,
                                                       rec1[oneKey],
                                                       rec2[oneKey]))
            counter = counter + 1
    if counter == 0:
        print("No differences found")

if __name__ == "__main__":
    compare_recs(example_file,1093,1137)
