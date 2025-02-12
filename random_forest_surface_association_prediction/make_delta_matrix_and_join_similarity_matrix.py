import argparse
import csv
import pandas as pd
import scipy.stats
import sklearn.metrics
pd.options.mode.chained_assignment = None
import numpy as np
import scipy.stats as st
import math
import itertools
import statsmodels.stats.multitest as sm
import sklearn as sk
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve
from sklearn import metrics
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import openai
import umap.umap_ as umap
from sklearn.cluster import KMeans
from scipy.cluster import hierarchy
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import cdist
from scipy.spatial import distance
import networkx as nx
from Bio.PDB import MMCIFParser
from adjustText import adjust_text
from Bio import SeqIO
from scipy.stats import mannwhitneyu
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from scipy.stats import entropy
import scipy.cluster.hierarchy as sch
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind
import json
from Bio import pairwise2
from Bio.Seq import Seq

# Function to normalize data using min-max scaling
def min_max_normalize(data):
    min_val = min(data)
    max_val = max(data)
    normalized_data = [(x - min_val) / (max_val - min_val) for x in data]
    return normalized_data

# Function to create a delta matrix by calculating pairwise differences between protein expression profiles
def delta_matrix_formation():
    matrix_file = 'path_to_matrix_file.csv'
    matrix = pd.read_csv(matrix_file)
    og_matrix = matrix.filter(items=['Protein ID', 'Gene'])
    try:
        matrix = matrix.drop('Gene', axis=1)
    except:
        pass

    matrix = matrix.fillna(0)
    prot_list = list(matrix['Protein ID'])

    combo_list = itertools.combinations(prot_list, 2)
    exp_matrix = matrix.set_index(list(matrix)[0])
    #exp_matrix = exp_matrix.drop('Gene', axis=1)
    exp_list = exp_matrix.columns.values.tolist()
    delta_list = []
    normalized_delta_list = []
    prot1 = []
    prot2 = []
    # Iterating over each protein pair
    for item in combo_list:
        tmp_df_1 = matrix.loc[matrix['Protein ID'].isin([item[0]])]
        tmp_df_1 = tmp_df_1.set_index(list(tmp_df_1)[0])

        tmp_df_2 = matrix.loc[matrix['Protein ID'].isin([item[1]])]
        tmp_df_2 = tmp_df_2.set_index(list(tmp_df_2)[0])

        # Creating lists of expression values for each protein
        vector_1 = list(tmp_df_1.loc[item[0], :])
        nonzero_vector_1 = []
        [nonzero_vector_1.append(item) for item in vector_1 if item != 0]
        normal_vector_1 = min_max_normalize(vector_1)

        vector_2 = list(tmp_df_2.loc[item[1], :])
        nonzero_vector_2 = []
        [nonzero_vector_2.append(item) for item in vector_2 if item != 0]
        normal_vector_2 = min_max_normalize(vector_2)

        # Proceed only if both proteins have enough non-zero values
        if len(nonzero_vector_1) >= 20 and len(nonzero_vector_2) >= 20:
            tmp_delta = []
            tmp_norm_delta = []
            i=0
            for value in vector_1:
                tmp_delta.append(value - vector_2[i])
                i+=1
            i=0
            for value in normal_vector_1:
                tmp_norm_delta.append(value-normal_vector_2[i])
                i+=1
            prot1.append(item[0])
            prot2.append(item[1])
            delta_list.append(tmp_norm_delta)
        else:
            continue

    # Creating a DataFrame from the delta list
    delta_df = pd.DataFrame(columns = exp_list, data=delta_list)

    delta_df['Protein ID 1'] = prot1
    delta_df['Protein ID 2'] = prot2

    # Merging original gene information with the delta matrix
    tmp_delta_df_1 = delta_df.filter(items=['Protein ID 1'])
    tmp_delta_df_2 = delta_df.filter(items=['Protein ID 2'])
    tmp_delta_df_1 = tmp_delta_df_1.rename(columns={'Protein ID 1': 'Protein ID'})
    tmp_delta_df_2 = tmp_delta_df_2.rename(columns={'Protein ID 2': 'Protein ID'})

    tmp_delta_df_1 = tmp_delta_df_1.merge(og_matrix, on='Protein ID', how='left')
    tmp_delta_df_2 = tmp_delta_df_2.merge(og_matrix, on='Protein ID', how='left')

    # Adding Gene information to the delta matrix
    delta_genes_1 = list(tmp_delta_df_1['Gene'])
    delta_genes_2 = list(tmp_delta_df_2['Gene'])

    delta_df['Gene'] = delta_genes_1
    delta_df['Gene 2'] = delta_genes_2

    # Creating a paired column with sorted gene names
    pair_list = []
    i = 0
    for item in delta_genes_1:
        pair = sorted([item, delta_genes_2[i]])
        pair_list.append(str(pair[0]) + ':' + str(pair[1]))
        i += 1
    delta_df['Paired'] = pair_list

    return(delta_df)

#this little bit will combine the delta matrix and similarity matrix files
def join_delta_similarity_matrices():
    similarity_matrix_file = 'path_to_similarity_matrix.csv'
    similarity_matrix = pd.read_csv(similarity_matrix_file)

    # Normalizing and inverting distance measures (Euclidean, Bray-Curtis, Manhattan, Chebyshev)
    inverted_list = ['euclidean','bray_curtis', 'manhattan', 'chebyshev']
    new_euclidean = []
    new_bc = []
    new_manhattan = []
    new_chebyshev = []
    i=0
    # Processing each similarity measure
    for item in inverted_list:
        tmp_list = list(similarity_matrix[item])
        i+=1
        for value in tmp_list:
            new_value = float(1/value) # Inverting the values
            if i==1:
                new_euclidean.append(new_value)
            elif i==2:
                new_bc.append(new_value)
            elif i==3:
                new_manhattan.append(new_value)
            elif i==4:
                new_chebyshev.append(new_value)
            else:
                print('Ahhhh')
    new_euclidean = min_max_normalize(new_euclidean)
    new_bc = min_max_normalize(new_bc)
    new_manhattan = min_max_normalize(new_manhattan)
    new_chebyshev = min_max_normalize(new_chebyshev)

    # Updating the similarity matrix with the new values
    similarity_matrix['euclidean']= new_euclidean
    similarity_matrix['bray_curtis'] = new_bc
    similarity_matrix['manhattan'] = new_manhattan
    similarity_matrix['chebyshev'] = new_chebyshev

    # Getting the delta matrix
    delta_matrix = delta_matrix_formation()

    # Dropping unnecessary columns
    delta_matrix = delta_matrix.drop('Protein ID 1', axis=1)
    delta_matrix = delta_matrix.drop('Protein ID 2', axis=1)
    delta_matrix = delta_matrix.drop('Gene', axis=1)
    delta_matrix = delta_matrix.drop('Gene 2', axis=1)

    # Merging the delta and similarity matrices
    combined_matrix = delta_matrix.merge(similarity_matrix, on='Paired', how='inner')

    # Outputting the combined matrix to a CSV file
    out = 'output_path_for_delta_matrix_joined_with_similarity_matrix.csv'
    combined_matrix.to_csv(out)
join_delta_similarity_matrices()