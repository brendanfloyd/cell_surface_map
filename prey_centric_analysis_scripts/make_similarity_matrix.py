import pandas as pd
import scipy.stats
import numpy as np
import itertools
import matplotlib.pyplot as plt
from scipy.spatial import distance
from sklearn import metrics

pd.options.mode.chained_assignment = None  # Suppress chained assignment warning

# This script calculates pairwise correlations between protein intensity profiles to identify interactors.
# The correlations are computed using Spearman and Pearson coefficients, along with various distance metrics.
def make_similarity_matrix():
    # Load fold change matrix file containing protein intensity data
    fc_matrix_file = 'path_to_fold_change_matrix'
    fc_matrix = pd.read_csv(fc_matrix_file)
    og_fc_matrix = fc_matrix.filter(items=['Protein ID', 'Gene'])  # Keep original protein ID and gene mappings

    # Remove 'Gene' column if it exists to facilitate numerical analysis
    try:
        fc_matrix = fc_matrix.drop('Gene', axis=1)
    except:
        pass

    fc_matrix = fc_matrix.fillna(0)  # Fill missing values with zeros
    prot_list = list(fc_matrix['Protein ID'])  # List of protein IDs
    combo_list = itertools.combinations(prot_list, 2)  # Generate all unique protein pairs
    corr_list = []  # Store correlation results

    # Iterate over each protein pair and compute correlation metrics
    for item in combo_list:
        # Extract intensity vectors for the two proteins
        tmp_df_1 = fc_matrix.loc[fc_matrix['Protein ID'].isin([item[0]])].set_index(list(fc_matrix)[0])
        tmp_df_2 = fc_matrix.loc[fc_matrix['Protein ID'].isin([item[1]])].set_index(list(fc_matrix)[0])

        vector_1 = list(tmp_df_1.loc[item[0], :])
        vector_2 = list(tmp_df_2.loc[item[1], :])

        # Remove zero values to ensure meaningful correlation calculations
        nonzero_vector_1 = [val for val in vector_1 if val != 0]
        nonzero_vector_2 = [val for val in vector_2 if val != 0]

        # Only compute correlations if both proteins have sufficient nonzero values
        if len(nonzero_vector_1) >= 5 and len(nonzero_vector_2) >= 5:
            spearman = scipy.stats.spearmanr(vector_1, vector_2)
            pearson = scipy.stats.pearsonr(vector_1, vector_2)
            euclid = distance.euclidean(vector_1, vector_2)
            bc = distance.braycurtis(vector_1, vector_2)
            manhattan = metrics.pairwise.manhattan_distances([vector_1], [vector_2])
            chebyshev = distance.chebyshev(vector_1, vector_2)

            # Store computed correlation and distance metrics
            corr_list.append(
                [item[0], item[1], pearson.statistic, spearman.statistic, euclid, bc, manhattan[0][0], chebyshev])

    # Create a DataFrame to store correlation results
    corr_df = pd.DataFrame(columns=['Protein ID 1', 'Protein ID 2', 'pearson', 'spearman', 'euclidean',
                                    'bray_curtis', 'manhattan', 'chebyshev'], data=corr_list)

    # Extract and merge gene names for each protein pair
    tmp_corr_df_1 = corr_df.filter(items=['Protein ID 1']).rename(columns={'Protein ID 1': 'Protein ID'})
    tmp_corr_df_2 = corr_df.filter(items=['Protein ID 2']).rename(columns={'Protein ID 2': 'Protein ID'})

    tmp_corr_df_1 = tmp_corr_df_1.merge(og_fc_matrix, on='Protein ID', how='left')
    tmp_corr_df_2 = tmp_corr_df_2.merge(og_fc_matrix, on='Protein ID', how='left')

    corr_df['Gene 1'] = list(tmp_corr_df_1['Gene'])
    corr_df['Gene 2'] = list(tmp_corr_df_2['Gene'])

    # Generate a list of unique gene pairs
    pair_list = [
        f"{sorted([corr_df['Gene 1'][i], corr_df['Gene 2'][i]])[0]}:{sorted([corr_df['Gene 1'][i], corr_df['Gene 2'][i]])[1]}"
        for i in range(len(corr_df))]
    corr_df['Paired'] = pair_list

    # Plot histogram of Pearson correlation values
    plt.hist(corr_df['pearson'], bins=20, color='dimgray', edgecolor='black')

    # Customize axis labels and font sizes
    xticks_positions = [-0.5, 0, 0.5]  # Specify tick positions
    xticks_labels = [f"{x:.2f}" for x in xticks_positions]  # Format tick labels
    plt.xticks(xticks_positions, xticks_labels, fontsize=24)
    plt.xlabel('Pearson Coefficient', fontsize=36)
    plt.ylabel('Count', fontsize=36)

    # Save and display the plot
    png_out = 'path_to_png_outfile'
    plt.savefig(png_out)
    plt.show()

make_similarity_matrix()
