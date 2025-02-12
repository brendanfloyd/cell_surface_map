import argparse
import csv
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import scipy.stats as st
import math
import matplotlib.pyplot as plt
import seaborn as sns

# Set up argument parser to take file paths as input parameters
parser = argparse.ArgumentParser(
    description='Make Perseus output file appropriate for volcano plot and network addition')
parser.add_argument('-c', '--clean_file_path', type=str,
                    help='Output from diann_out_to_perseus_in file path',
                    )
parser.add_argument('-p', '--perseus_file_path', type=str,
                    help='Perseus out file .tsv',
                    )
parser.add_argument('-o', '--out_path', type=str,
                    help='Output file path',
                    )
args = parser.parse_args()

# Read input data files
exp_data = pd.read_csv(args.clean_file_path)  # Read cleaned DIA-NN output file
perseus = pd.read_csv(args.perseus_file_path, sep='\t')  # Read Perseus output file

# Extract q-value from Perseus data
perseus['qvalue'] = perseus["Student's T-test q-value target log Intensity_control log Intensity"]

# Remove first two rows, assuming they contain metadata rather than data
perseus = perseus.iloc[2:]

# Initialize lists to store fold change values
fc_list = list(perseus['fold_change'])
prot_list = list(perseus['Protein ID'])
final_fc_list = []
i = 0

# Iterate through fold change list and impute missing values
for item in fc_list:
    if str(item) in 'nan':  # Check for NaN values
        tmp_df = perseus.loc[perseus['Protein ID'].isin([prot_list[i]])]  # Filter Perseus data for the protein

        # Convert log intensities back to normal intensities using base-2 exponentiation
        tmp_df['target_1 Intensity'] = np.exp2(float(list(tmp_df['target_1 log Intensity'])[0]))
        tmp_df['target_2 Intensity'] = np.exp2(float(list(tmp_df['target_2 log Intensity'])[0]))
        tmp_df['target_3 Intensity'] = np.exp2(float(list(tmp_df['target_3 log Intensity'])[0]))
        tmp_df['control_1 Intensity'] = np.exp2(float(list(tmp_df['control_1 log Intensity'])[0]))
        tmp_df['control_2 Intensity'] = np.exp2(float(list(tmp_df['control_2 log Intensity'])[0]))
        tmp_df['control_3 Intensity'] = np.exp2(float(list(tmp_df['control_3 log Intensity'])[0]))

        # Calculate mean intensities for target and control groups
        exp_mean_list = [np.nanmean([tmp_df['target_1 Intensity'].iloc[0],
                                     tmp_df['target_2 Intensity'].iloc[0],
                                     tmp_df['target_3 Intensity'].iloc[0]])]
        ctl_mean_list = [np.nanmean([tmp_df['control_1 Intensity'].iloc[0],
                                     tmp_df['control_2 Intensity'].iloc[0],
                                     tmp_df['control_3 Intensity'].iloc[0]])]

        tmp_df['target_mean_intensity'] = exp_mean_list
        tmp_df['control_mean_intensity'] = ctl_mean_list

        # Compute fold change with pseudocount to avoid division by zero
        tmp_df['fold_change'] = ((tmp_df['target_mean_intensity'] + 1) - (
                tmp_df['control_mean_intensity'] + 1)) / \
                                (tmp_df['control_mean_intensity'] + 1)

        final_fc_list.append(list(tmp_df['fold_change'])[0])
    else:
        final_fc_list.append(item)
    i += 1

# Update Perseus data with new fold change values
perseus['fold_change'] = final_fc_list
perseus['fold_change'] = perseus['fold_change'].astype(float)
perseus['log_fold_change'] = np.log2(perseus['fold_change'] + 1)

# Rename intensity columns for clarity
perseus['perseus target_1 log Intensity'] = perseus['target_1 log Intensity']
perseus['perseus target_2 log Intensity'] = perseus['target_2 log Intensity']
perseus['perseus target_3 log Intensity'] = perseus['target_3 log Intensity']
perseus['perseus control_1 log Intensity'] = perseus['control_1 log Intensity']
perseus['perseus control_2 log Intensity'] = perseus['control_2 log Intensity']
perseus['perseus control_3 log Intensity'] = perseus['control_3 log Intensity']

# Retain only relevant columns
perseus = perseus.filter(items=['Protein ID', 'fold_change', 'log_fold_change', 'qvalue',
                                'perseus target_1 log Intensity', 'perseus target_2 log Intensity',
                                'perseus target_3 log Intensity', 'perseus control_1 log Intensity',
                                'perseus control_2 log Intensity', 'perseus control_3 log Intensity'])

# Merge Perseus data with experiment data
exp_data = exp_data.drop(columns=['fold_change', 'log_fold_change'], errors='ignore')
new_exp_data = exp_data.merge(perseus, on='Protein ID', how='left')

# Compute -log10 of q-values for visualization purposes
new_exp_data['qvalue'] = new_exp_data['qvalue'].astype(float)
new_exp_data['log_pvalue'] = -np.log10(new_exp_data['qvalue'])

# Filter out specific unwanted genes from the dataset
gene_list = list(new_exp_data['Gene'])
new_gene_list = [item for item in gene_list if
                 not any(x in str(item) for x in ['IGKV', 'IGHV', 'IGHG', 'IGHA', 'IGLV', 'IGL', 'IGKC', 'KRT', 'ALB'])]
new_exp_data = new_exp_data.loc[new_exp_data['Gene'].isin(new_gene_list)]

# Save the processed data to the specified output file
out = args.out_path
print(out)
new_exp_data.to_csv(args.out_path)