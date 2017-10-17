#!/usr/bin/env/python

##  NODDI microstructural model (Zhang et al, 2012).  Uses AMICO toolbox (Daducci et al, 2015; https://github.com/daducci/AMICO)

import os,sys

def NODDI_fxn(subj):
	import amico
	amico.core.setup()
	ae = amico.Evaluation("/N/dc2/projects/lifebid/Concussion/concussion_real/NODDI",subj)
	ae.load_data(dwi_filename = "data_acpc.nii.gz", scheme_filename = "bvals.scheme", mask_filename = "nodif_brain_mask.nii.gz", b0_thr = 1)
	ae.set_model("NODDI")
	ae.generate_kernels()
	ae.load_kernels()
	ae.fit()
	ae.save_results()
	print("NODDI model generation complete")

if __name__ == '__main__':
    import sys
    function = getattr(sys.modules[__name__], sys.argv[1])
    subj = sys.argv[2]
    NODDI_fxn(subj)
