import csv
import pandas as pd
import numpy as np
import random
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier


def delta_matrix_label(delta_matrix_file, ppi_file, string_ppi_file):
    """
    Loads and processes the delta similarity matrix and reference protein-protein interaction (PPI) datasets.

    Args:
        delta_matrix_file (str): Path to the delta similarity matrix CSV file.
        ppi_file (str): Path to the gold standard PPI CSV file.
        string_ppi_file (str): Path to the STRING PPI CSV file.

    Returns:
        tuple: Processed training data matrix, training labels, withheld data matrix, full dataset matrix.
    """
    delta_matrix = pd.read_csv(delta_matrix_file)
    pairs = list(delta_matrix['Paired'])

    gold_ppi = pd.read_csv(ppi_file)
    complex_pairs = list(gold_ppi['Paired'])
    delta_complex_pairs = delta_matrix.loc[delta_matrix['Paired'].isin(complex_pairs)]
    hc_complex_pairs = delta_complex_pairs.loc[delta_complex_pairs['pearson'] >= 0.3]

    perc = 0.75  # Percentage of high-confidence pairs used for training
    hc_pairs = list(hc_complex_pairs['Paired'])
    shuffled_hc_pairs = random.sample(hc_pairs, len(hc_pairs))
    withheld_hc_pairs = shuffled_hc_pairs[round(len(shuffled_hc_pairs) * perc):]

    string_ppi = pd.read_csv(string_ppi_file)
    full_train_pairs = string_ppi.loc[string_ppi['Paired'].isin(pairs)]
    lc_string_ppi = delta_matrix.drop_duplicates('Paired')
    lc_string_ppi = lc_string_ppi.loc[~lc_string_ppi['Paired'].isin(hc_pairs)]

    lc_pairs = list(lc_string_ppi['Paired'])
    shuffled_lc_pairs = random.sample(lc_pairs, len(lc_pairs))
    final_lc_pairs = shuffled_lc_pairs[:10000]
    withheld_lc_pairs = shuffled_lc_pairs[10000:20000]

    withheld_pairs = withheld_hc_pairs + withheld_lc_pairs
    all_pairs = hc_pairs + final_lc_pairs

    train_delta_matrix = delta_matrix.loc[delta_matrix['Paired'].isin(all_pairs)]
    non_train_delta_matrix = delta_matrix.loc[delta_matrix['Paired'].isin(withheld_pairs)]
    full_delta_matrix = delta_matrix.drop(['Protein ID 1', 'Protein ID 2', 'Gene 1', 'Gene 2'], axis=1).set_index(
        'Paired')

    training_labels = [1 if item in hc_pairs else 0 for item in all_pairs]

    return train_delta_matrix.set_index('Paired'), training_labels, non_train_delta_matrix.set_index(
        'Paired'), full_delta_matrix


def rf_train_test_split(delta_matrix_file, ppi_file, string_ppi_file):
    """
    Splits the dataset into training and testing sets for Random Forest classification.

    Args:
        delta_matrix_file (str): Path to the delta similarity matrix CSV file.
        ppi_file (str): Path to the gold standard PPI CSV file.
        string_ppi_file (str): Path to the STRING PPI CSV file.

    Returns:
        tuple: Scaled training and testing data, labels, and full dataset matrix.
    """
    X, y, non_train_delta_matrix, full_delta_matrix = delta_matrix_label(delta_matrix_file, ppi_file, string_ppi_file)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, train_size=0.7, test_size=0.3, random_state=42)

    return X_train, X_test, y_train, y_test, non_train_delta_matrix, X, full_delta_matrix, y


def run_random_forest_with_cv(delta_matrix_file, ppi_file, string_ppi_file):
    """
    Trains a Random Forest classifier using cross-validation.

    Args:
        delta_matrix_file (str): Path to the delta similarity matrix CSV file.
        ppi_file (str): Path to the gold standard PPI CSV file.
        string_ppi_file (str): Path to the STRING PPI CSV file.
    """
    X_train, X_test, y_train, y_test, non_train_delta_matrix, X, full_delta_matrix, y = rf_train_test_split(
        delta_matrix_file, ppi_file, string_ppi_file)

    rf_classifier = RandomForestClassifier(n_estimators=100, max_depth=10, random_state=42)
    rf_classifier.fit(X_train, y_train)

    test_score = rf_classifier.score(X_test, y_test)
    print(f"Random Forest Test Accuracy: {test_score:.4f}")
