import argparse
import pandas as pd
import numpy as np

# Suppress chained assignment warnings
pd.options.mode.chained_assignment = None

bait_results_file_1 = 'path_to_bait_experiment_1'
bait_results_file_2 = 'path_to_bait_experiment_2'

def protein_fc_matrix_assembler(target_id, data_file):
    """
    Reads experimental data, filters for valid fold changes, and retains proteins present in a given network file.

    Parameters:
    - target_id (str): Identifier for the target protein.
    - data_file (str): Path to the experimental data file.
    - network_file (str): Path to the network file containing valid protein IDs.

    Returns:
    - pd.DataFrame: Filtered dataframe containing protein IDs and fold changes.
    """
    exp_data = pd.read_csv(data_file)
    exp_data = exp_data.dropna(subset=['fold_change'])

    # Add bait column
    exp_data['Bait'] = target_id

    # Filter out unwanted genes
    unwanted_genes = ['IGKV', 'IGHV', 'IGHG', 'IGHA', 'IGLV', 'IGL', 'KRT', 'ALB']
    exp_data = exp_data[~exp_data['Gene'].astype(str).str.contains('|'.join(unwanted_genes), na=False)]

    # Rename fold change column
    new_name = f"{target_id}_fold_change"
    exp_data[new_name] = exp_data['fold_change']

    # Retain only relevant cell surface proteome proteins if you have a known surfaceome. Uncomment if needed.
    '''
    cell_surface_prots = [list_of_cell_surface_uniprot_ids]
    exp_data = exp_data.loc[exp_data['Protein ID'].isin(cell_surface_prots)]   
    '''
    return exp_data[['Protein ID', new_name]]

def fc_matrix_assembler(output_file):
    """
    Assembles a fold change matrix for multiple bait proteins, filters based on presence in experiments,
    and adds gene annotations.
    Note that you can add as many baits as you have sampled. In the case of Jurkats, 65 bait experiments were used
    Parameters:
    - results_file_1 (str): Path to perseus clean output results file for 1 bait experiment.
    - results_file_2 (str): Path to perseus clean output results file for 1 bait experiment.
    - output_file (str): Path to save the assembled fold change matrix.
    """
    # Generate fold change matrices for target proteins. Add more as needed
    bait_result_1 = protein_fc_matrix_assembler('Bait gene name 1', bait_results_file_1, network_file)
    bait_result_2 = protein_fc_matrix_assembler('Bait gene name 1', tfrc_results_file, network_file)

    # Merge results
    complete_df = bait_result_1.merge(bait_result_2, on='Protein ID', how='outer').fillna(0)

    # Filter proteins appearing in at least 10 experiments. Change this to desired experimental observation threshold.
    protein_counts = (complete_df.set_index('Protein ID') != 0).sum(axis=1)
    complete_df = complete_df[protein_counts >= 10]

    # Save output
    complete_df.to_csv(output_file, index=False)

fc_matrix_assembler(output_file)