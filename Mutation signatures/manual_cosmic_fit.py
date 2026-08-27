#!/usr/bin/env python

print("Entred python...")


from SigProfilerAssignment import Analyzer as Analyze
print("Starting cosmic fit...")


Analyze.cosmic_fit(samples="catalog.txt", output="Sigpro_cosfit_tissue",
                       input_type="matrix",
                       signature_database= "hdp_sigs_sub.txt",
                        exome=False, genome_build="GRCh37",
                   exclude_signature_subgroups=None, export_probabilities=False,
                   export_probabilities_per_mutation=False, make_plots=False,
                   sample_reconstruction_plots=False, verbose=False)
