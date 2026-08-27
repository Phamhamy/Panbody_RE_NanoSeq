#!/usr/bin/env python

print("Entred python...")

from SigProfilerAssignment import Analyzer as Analyze
print("Starting decompose fit...")


Analyze.decompose_fit(samples="sigpro_sbs96_sample_t100.txt", output="Sigpro_deconvolution_hdp96_sub",
                       input_type="matrix", signatures="hdp_sigs.txt",
                      signature_database = "./cosmic_SBS_signatures.txt", 
                      genome_build="GRCh37", cosmic_version=3.4,
                    verbose=False, make_plots= False,
                       new_signature_thresh_hold=0.9,
                       exome=False)
