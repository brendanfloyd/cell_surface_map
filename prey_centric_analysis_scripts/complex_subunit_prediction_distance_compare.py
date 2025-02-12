import pandas as pd
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from adjustText import adjust_text


# Function to calculate pairwise distances between protein subunits
# Uses a CIF file to extract structural information
# Returns a dataframe with subunit distances

def get_complex_subunit_distance():
    cif_file = 'path_to_protein_structure_cif_file.cif'
    distances, subunit_centers = calculate_pairwise_distances_from_cif(cif_file)

    # Dictionary mapping PDB chain identifiers to gene names
    name_dict = {'a': 'CD247', 'b': 'CD247', 'd': 'CD3D', 'e': 'CD3E', 'f': 'CD3E', 'g': 'CD3G', 'm': 'TRAC',
                 'n': 'TRBC1'}
    dist_list = []

    for (chain1, chain2), distance in distances.items():
        print(f"Distance between subunit {name_dict[chain1]} and subunit {name_dict[chain2]}: {distance:.2f} Å")
        pair = sorted([name_dict[chain1], name_dict[chain2]])
        paired = str(pair[0]) + ':' + str(pair[1])
        dist_list.append([paired, pair[0], pair[1], distance, distance * 0.1])
        print([paired, distance])

    dist_df = pd.DataFrame(columns=['Paired', 'Gene 1', 'Gene 2', 'angstrom_distance', 'nm_distance'], data=dist_list)
    return dist_df


# Function to compare subunit distances with predicted interaction scores
# Merges experimental distances with computational predictions
# Generates a scatter plot comparing Pearson correlation with subunit distances

def complex_subunit_prediction_distance_compare():
    dist_df = get_complex_subunit_distance()

    # Averaging distances for complexes with multiple copies of a subunit
    new_dist_list = []
    for item in set(list(dist_df['Paired'])):
        tmp_df = dist_df.loc[dist_df['Paired'].isin([item])]
        gene1 = list(tmp_df['Gene 1'])[0]
        gene2 = list(tmp_df['Gene 2'])[0]

        if len(tmp_df) == 1:
            new_dist_list.append(tmp_df.values.tolist()[0])
        else:
            mean_ang = np.mean(list(tmp_df['angstrom_distance']))
            mean_nm = np.mean(list(tmp_df['nm_distance']))
            new_dist_list.append([item, gene1, gene2, mean_ang, mean_nm])

    dist_df = pd.DataFrame(columns=['Paired', 'Gene 1', 'Gene 2', 'angstrom_distance', 'nm_distance'],
                           data=new_dist_list)

    # Load predicted interaction scores
    prediction_file = 'path_to_similarity_matrix'
    pred_df = pd.read_csv(prediction_file)
    pred_df = pred_df.filter(
        items=['Paired', 'pearson', 'spearman', 'euclidean', 'bray_curtis', 'manhattan', 'chebyshev'])

    # Merge experimental distances with predictions
    merge_df = dist_df.merge(pred_df, on='Paired', how='left')
    print(merge_df)

    # Save merged data to CSV
    out = 'output_path_for_csv_with_distances_and_predictions.csv'
    merge_df.to_csv(out)
    merge_df = merge_df.dropna()

    # Compute correlation statistics
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(list(merge_df['nm_distance']),
                                                                         list(merge_df['pearson']))
    print(r_value, p_value)
    spearman = scipy.stats.spearmanr(list(merge_df['nm_distance']), list(merge_df['pearson']))
    print(spearman.statistic)
    rho = spearman.statistic
    spear_pvalue = spearman.pvalue

    # Generate scatter plot comparing distance with Pearson correlation
    plt.figure(figsize=(8, 6))
    plt.scatter(merge_df['nm_distance'], merge_df['pearson'], color='black', marker='o', s=150)

    line_values = [slope * x + intercept for x in list(merge_df['nm_distance'])]
    plt.plot(merge_df['nm_distance'], line_values, 'r-', label=f'Line of best fit: y={slope:.2f}x+{intercept:.2f}',
             color='black')

    # Annotate each point with the 'Paired' label
    texts = []
    for x, y, label in zip(merge_df['nm_distance'], merge_df['pearson'], merge_df['Paired']):
        texts.append(plt.text(x, y, label, fontsize=28, ha='center', va='bottom', color='black'))

    adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=2))

    # Label axes
    plt.xlabel('Subunit pairwise distance (nm)', fontsize=28)
    plt.ylabel('Pearson r', fontsize=28)

    # Annotate plot with Spearman correlation results
    plt.annotate(f'Spearman rho: {rho:.2f}\np-value: {spear_pvalue:.4f}', xy=(1.1, 0.95), xycoords='axes fraction',
                 fontsize=26, ha='right')

    # Set font sizes for axis ticks
    plt.xticks(fontsize=28)
    plt.yticks(fontsize=28)
    plt.tight_layout()

    # Save plot
    png_out = 'out_file_for_comparison_plot.png'
    plt.savefig(png_out)
    plt.show()

# Run comparison function
complex_subunit_prediction_distance_compare()