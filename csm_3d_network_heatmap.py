import argparse
import csv
import pandas as pd
pd.options.mode.chained_assignment = None  # Suppress pandas chained assignment warnings
import numpy as np
import scipy.stats as st
import math
import itertools
import statsmodels.stats.multitest as sm
import sklearn as sk
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve
from sklearn import metrics
import networkx as nx
import matplotlib.pyplot as plt
import netwulf as nw
import plotly.graph_objects as go
import plotly.io as pio
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib import cm
from scipy.interpolate import griddata
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import json
import ffmpeg

# Function to normalize data using the min-max scaling method
def min_max_normalize(data):
    min_val = min(data)
    max_val = max(data)
    normalized_data = [(x - min_val) / (max_val - min_val) for x in data]
    return normalized_data


# Function to process data for creating a network heatmap.
# The input is a protein of interest (POI) gene name to visualize its distribution on the cell surface.
def network_3d_heatmap_df_formatting(POI):
    print('Processing data files')

    # Load the network file (replace with actual file path)
    network_file = 'path_to_network_file.csv'
    network = pd.read_csv(network_file)

    # Load the fold change matrix file (replace with actual file path)
    fc_matrix_file = 'path_to_fc_matrix_file.csv'
    fc_matrix = pd.read_csv(fc_matrix_file)

    # Fill missing values in the fold change matrix with zero
    fc_matrix = fc_matrix.fillna(0)

    # Select rows where the Gene is in the list of genes corresponding to the POI
    poi_fc_matrix = fc_matrix.loc[fc_matrix['Gene'].isin([POI])]

    # Extract column names (except for 'Gene' column)
    fc_names = list(poi_fc_matrix.columns.values)
    fc_names = fc_names[1:len(fc_names) + 1]

    pure_fc_genes_dict = {}
    pure_fc_genes = []

    # Process fold change data for each gene in the POI fold change matrix
    for item in fc_names:
        tmp_fc = list(poi_fc_matrix[item])
        tmp_fc = tmp_fc[0]
        new_item = item[0:item.find('_f')]  # Extract gene name from column name
        pure_fc_genes_dict[new_item] = tmp_fc
        pure_fc_genes.append(new_item)

    # Filter only bait proteins from the network
    list_o_baits = set(list(network['Bait']))
    print(list_o_baits)

    # Combine network data with bait protein information
    combined_result = network.loc[network['Gene'].isin(list_o_baits)]

    # Adding reciprocal nodes (edges)
    cr_gene_list = list(combined_result['Gene'])
    cr_bait_list = list(combined_result['Bait'])
    add_gene_list = []
    add_bait_list = []
    i = 0
    for gene in cr_gene_list:
        if gene not in cr_bait_list:
            tmpdf = combined_result.loc[combined_result['Gene'].isin([gene])]
            tmpdf = tmpdf.loc[tmpdf['Bait'].isin([cr_bait_list[i]])]
            tmpdf['Gene'] = cr_bait_list[i]
            tmpdf['Bait'] = gene
            tmplist = tmpdf.values.tolist()[0]
            combined_result.loc[len(combined_result)] = tmplist  # Add the new row

        i += 1

    # Extract final bait proteins and their corresponding fold change values
    final_bait_list = list(combined_result['Bait'])
    final_fc_list = []

    for item in final_bait_list:
        final_fc_list.append(pure_fc_genes_dict[item])

    # Normalize the fold change values
    normal_fc_list = min_max_normalize(final_fc_list)
    print(set(final_fc_list))
    print(normal_fc_list)

    # Add new columns for fold change values to the network results
    combined_result['poi_log_fold_change'] = final_fc_list
    combined_result['poi_normalized_fold_change'] = normal_fc_list

    # Output file (replace with actual output file path)
    combined_result_out = 'path_to_output_file.csv'

    return (list_o_baits, combined_result, fc_matrix)


# Function to create a 3D network heatmap for the given POI
def network_3d_heatmap(POI):
    list_o_baits, combined_result, fc_matrix = network_3d_heatmap_df_formatting(POI)

    # Set the random seed for reproducibility
    np.random.seed(42)

    # Process the fold change matrix by dropping unnecessary columns
    new_fc_matrix = fc_matrix.drop('Gene', axis=1)
    new_fc_matrix = new_fc_matrix.drop('Protein ID', axis=1)

    # Find the min and max fold change values
    max_fc_value = new_fc_matrix.to_numpy().max()
    min_fc_value = new_fc_matrix.to_numpy().min()
    print("Plotting...")

    # Create a graph using networkx
    G = nx.Graph()
    edge_list = [(row['Bait'], row['Gene'], {'fc': row['poi_normalized_fold_change']}) for _, row in
                 combined_result.iterrows()]
    G.add_edges_from(edge_list)

    # Add nodes to the graph with fold change values
    for node in combined_result['Bait']:
        G.add_node(node,
                   fc=combined_result.loc[combined_result['Bait'] == node, 'poi_normalized_fold_change'].values[0])

    # Generate positions for the nodes using a layout algorithm (Kamarda-Kawai)
    pos = nx.kamada_kawai_layout(G, dim=3)
    radius = 1.0

    # Normalize node positions
    for node, position in pos.items():
        norm = np.linalg.norm(position)
        pos[node] = position / norm * radius

    # Convert node positions to spherical coordinates
    points = np.array([pos[node] for node in G])
    spherical_coords = np.array(
        [[np.arctan2(np.sqrt(p[0] ** 2 + p[1] ** 2), p[2]), np.arctan2(p[1], p[0])] for p in points])
    values = np.array([G.nodes[node]['fc'] for node in G])
    print(len(values))

    # Identify the nearest nodes to the North and South poles (theta = 0 and pi)
    north_pole_nearest = np.argmin(spherical_coords[:, 0])  # Closest to theta=0
    south_pole_nearest = np.argmax(spherical_coords[:, 0])  # Closest to theta=pi

    # Set values for the poles based on nearest nodes
    north_pole_value = 1 * values[north_pole_nearest]
    south_pole_value = 1 * values[south_pole_nearest]

    # Adding artificial nodes at the poles with zero values
    pole_positions = {
        'north_pole': [0, 0, 1 * radius],
        'south_pole': [0, 0, -1 * radius]
    }

    # Add north and south pole nodes
    i = 0
    for pole in pole_positions:
        if i == 0:
            G.add_node(pole, fc=north_pole_value)
            pos[pole] = pole_positions[pole]
            i += 1
        else:
            G.add_node(pole, fc=south_pole_value)
            pos[pole] = pole_positions[pole]

    # Prepare spherical grid for visualizing the network
    theta = np.linspace(0, np.pi, 100)
    phi = np.linspace(0, 2 * np.pi + np.pi / 50, 100)  # Extend beyond 2*pi
    theta, phi = np.meshgrid(theta, phi)
    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)

    # Interpolate fold change values over the spherical grid
    grid_x, grid_y, grid_z = np.ravel(x), np.ravel(y), np.ravel(z)
    grid_theta = np.arctan2(np.sqrt(grid_x ** 2 + grid_y ** 2), grid_z)
    grid_phi = np.mod(np.arctan2(grid_y, grid_x), 2 * np.pi + np.pi / 50)
    grid_values = griddata(spherical_coords, values, (grid_theta, grid_phi), method='cubic')

    # Create the plot
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    grid_values = grid_values.reshape(x.shape)

    # Set color normalization
    norm = Normalize(vmin=0, vmax=1)

    # Plot the surface of the sphere with heatmap colors
    ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=cm.plasma(norm(grid_values)), shade=False)

    # Plot the nodes on the surface
    for node, position in pos.items():
        ax.scatter(*position, color=cm.plasma(norm(G.nodes[node]['fc'])), s=0)

    # Add colorbar for fold change values
    mappable = cm.ScalarMappable(cmap=cm.plasma, norm=norm)
    mappable.set_array(grid_values)
    fig.colorbar(mappable, shrink=0.5, aspect=5)

    ax.set_box_aspect([1, 1, 1])  # Set equal aspect ratio

    # Display the plot
    plt.show()


# Calling the function with a specified protein of interest (POI)
POI = 'Your_Protein_Of_Interest'  # Replace with actual protein name
network_3d_heatmap(POI)
