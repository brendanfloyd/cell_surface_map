import argparse
import csv
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import scipy.stats as st
import math

# Set up argument parser to handle command-line inputs
parser = argparse.ArgumentParser(description='Make DIA-NN output file appropriate for Perseus')
parser.add_argument('file_path', metavar='p', type=str,
                    help='DIA-NN protein group .pg file path')
# parser.add_argument('proteintsv_file_path', metavar='u', type=str,
#                     help='Fragpipe protein.tsv file that has unique peptide info')
parser.add_argument('out_path', metavar='o', type=str,
                    help='Output file path')
args = parser.parse_args()

# Read DIA-NN output file and replace column names
exp_data = pd.read_csv(args.file_path, sep='\t')

# Log transform intensity values to normalize distribution
exp_data['target_1 log Intensity'] = np.log2(exp_data['target_1 Intensity'])
exp_data['target_2 log Intensity'] = np.log2(exp_data['target_2 Intensity'])
exp_data['target_3 log Intensity'] = np.log2(exp_data['target_3 Intensity'])
exp_data['control_1 log Intensity'] = np.log2(exp_data['control_1 Intensity'])
exp_data['control_2 log Intensity'] = np.log2(exp_data['control_2 Intensity'])
exp_data['control_3 log Intensity'] = np.log2(exp_data['control_3 Intensity'])

# Extract intensity values into lists for processing
exp_uin_1 = list(exp_data['target_1 Intensity'])
exp_uin_2 = list(exp_data['target_2 Intensity'])
exp_uin_3 = list(exp_data['target_3 Intensity'])
ctl_uin_1 = list(exp_data['control_1 Intensity'])
ctl_uin_2 = list(exp_data['control_2 Intensity'])
ctl_uin_3 = list(exp_data['control_3 Intensity'])

# Compute mean intensity values for target and control groups
exp_mean_list = []
ctl_mean_list = []
i = 0
for item in exp_uin_1:
    mean = np.nanmean([item, exp_uin_2[i], exp_uin_3[i]])  # Compute mean while ignoring NaNs
    exp_mean_list.append(mean)
    i += 1

i = 0
for item in ctl_uin_1:
    mean = np.nanmean([item, ctl_uin_2[i], ctl_uin_3[i]])  # Compute mean while ignoring NaNs
    ctl_mean_list.append(mean)
    i += 1

# Store computed mean intensity values in the dataframe
exp_data['target_mean_intensity'] = exp_mean_list
exp_data['control_mean_intensity'] = ctl_mean_list

# Compute fold change between target and control groups
exp_data['fold_change'] = ((exp_data['target_mean_intensity'] + 1) - (exp_data['control_mean_intensity'] + 1)) / \
                          (exp_data['control_mean_intensity'] + 1)

# Compute log2 fold change for better visualization in downstream analysis
exp_data['log_fold_change'] = np.log2(exp_data['fold_change'] + 1)

# Extract protein IDs from dataset
all_prots = list(exp_data['Protein ID'])

# Print number of proteins processed
print(len(exp_data))

# Save the processed data to the specified output path
out = args.out_path
exp_data.to_csv(out)