#!/usr/bin/env python

import os
import multiprocessing
import shutil
import argparse
import pandas as pd
import numpy as np
from joblib import parallel_config


TARGET_CORES = 16
os.environ["LOKY_MAX_CPU_COUNT"] = str(TARGET_CORES)
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["JOBLIB_TEMP_FOLDER"] = "/tmp"

# Monkey-patch system CPU counters so SigProfilerAssignment respects our limits
os.cpu_count = lambda: TARGET_CORES
multiprocessing.cpu_count = lambda: TARGET_CORES

from SigProfilerAssignment import Analyzer as Analyze


RANDOM_SEED = 123
rng = np.random.default_rng(RANDOM_SEED)


# 1. Configuration

parser = argparse.ArgumentParser(
    description="Bootstrap signature deconvolution with 95% CIs."
)
parser.add_argument("table1_path", help="Path to donor SBS96 mutation counts table")
parser.add_argument("table2_path", help="Path to list of signatures per sample table")
parser.add_argument("table3_path", help="Path to COSMIC SBS signatures profiles table")
parser.add_argument("bootstrap_number", type=int, help="Number of bootstrap resamples to run")
args = parser.parse_args()
 
TABLE1_PATH = args.table1_path
TABLE2_PATH = args.table2_path
TABLE3_PATH = args.table3_path
BOOTSTRAPS = args.bootstrap_number
 
COUNTS_OUT = 'signature_presence_counts.tsv'
ACTIVITY_OUT = 'final_bootstrapped_activities.tsv'
CI_OUT = 'signature_confidence_intervals.tsv'
 
MUTATION_COL = 'Type' 
THRESHOLD_VAL = 0.95 * BOOTSTRAPS 
 
print("Loading data...")
mutations_df = pd.read_csv(TABLE1_PATH, sep='\t')
sample_sigs_df = pd.read_csv(TABLE2_PATH, sep='\t')
signatures_df = pd.read_csv(TABLE3_PATH, sep='\t')
 
presence_results = []
activity_results = []
ci_results = []
 

# 2. Iterate Through Samples

for index, row in sample_sigs_df.iterrows():
    sample_name = row['Sample']
    target_sigs = [sig.strip() for sig in str(row['Signatures']).split(',')]
    
    if sample_name not in mutations_df.columns:
        continue
    
    print(f"\nProcessing Sample: {sample_name}")
    
    counts = mutations_df[sample_name].values
    total_n = counts.sum()
    if total_n == 0: continue

    probabilities = counts / total_n
    bootstrapped_matrix = rng.multinomial(total_n, probabilities, size=BOOTSTRAPS).T
    
    bs_df = pd.DataFrame(bootstrapped_matrix, columns=[f"boot_{i}" for i in range(BOOTSTRAPS)])
    bs_df.insert(0, MUTATION_COL, mutations_df[MUTATION_COL])
    
    valid_sigs = [sig for sig in target_sigs if sig in signatures_df.columns]
    sample_sig_df = signatures_df[[MUTATION_COL] + valid_sigs]
    
    temp_dir = f"temp_{sample_name}"
    os.makedirs(temp_dir, exist_ok=True)
    bs_file = os.path.join(temp_dir, "bs_matrix.txt")
    sig_file = os.path.join(temp_dir, "bs_sigs.txt")
    
    bs_df.to_csv(bs_file, sep='\t', index=False)
    sample_sig_df.to_csv(sig_file, sep='\t', index=False)
    
    # Run Deconvolution
    try:
        with parallel_config(backend='loky', n_jobs=TARGET_CORES):
            Analyze.cosmic_fit(
                samples=bs_file, 
                output=temp_dir,
                input_type="matrix",
                signature_database=sig_file,
                add_background_signatures=False,
                nnls_add_penalty=0.015,
                nnls_remove_penalty=0.01,
                exome=False, genome_build="GRCh37",
                make_plots=False, verbose=False
            )
        
        # Process Results
        res_path = os.path.join(temp_dir, 'Assignment_Solution', 'Activities', 'Assignment_Solution_Activities.txt')
        
        if os.path.exists(res_path):
            res_df = pd.read_csv(res_path, sep='\t')
            sample_presence = {'Sample': sample_name}
            sample_activity = {'Sample': sample_name}
            
            for sig in valid_sigs:
                sig_data = res_df[sig]
                count = (sig_data > 0.01).sum()
                sample_presence[sig] = count
                
                median_val = sig_data.median()
                sample_activity[sig] = median_val if count >= THRESHOLD_VAL else 0
                
                ci_results.append({
                    'sequenceid': sample_name,
                    'signatures': sig,
                    'lower_95_CI': np.percentile(sig_data, 2.5),
                    'upper_95_CI': np.percentile(sig_data, 97.5),
                    'mean_attribution': sig_data.mean()
                })
                    
            presence_results.append(sample_presence)
            activity_results.append(sample_activity)
            
    except Exception as e:
        print(f"Error processing sample {sample_name}: {e}")
    finally:
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)


# 3. Final Table Generation

def finalize_table(results_list, out_path):
    if not results_list: return
    df = pd.DataFrame(results_list).fillna(0)
    cols = ['Sample'] + sorted([c for c in df.columns if c != 'Sample'])
    df[cols].to_csv(out_path, sep='\t', index=False)
    print(f"Saved: {out_path}")

finalize_table(presence_results, COUNTS_OUT)
finalize_table(activity_results, ACTIVITY_OUT)
pd.DataFrame(ci_results).to_csv(CI_OUT, sep='\t', index=False)
print("\nAll tasks completed.")
