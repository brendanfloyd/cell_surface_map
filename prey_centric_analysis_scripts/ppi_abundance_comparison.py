import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt

# Suppress chained assignment warning
pd.options.mode.chained_assignment = None

def ppi_abundance_scatter(prot_1_name, prot_2_name):
    # Define the path to the fold change matrix CSV file
    matrix_file = 'path_to_fold_change_matrix'

    # Read the CSV file into a pandas DataFrame
    matrix = pd.read_csv(matrix_file)

    # Replace NaN values with 0
    matrix = matrix.fillna(0)

    # Extract rows corresponding to the given protein names
    prot_1 = matrix.loc[matrix['Gene'].isin([prot_1_name])]
    prot_2 = matrix.loc[matrix['Gene'].isin([prot_2_name])]

    # Remove non-numeric columns and convert to list format
    prot_1 = prot_1.drop(['Gene', 'Protein ID'], axis=1).iloc[0].to_list()
    prot_2 = prot_2.drop(['Gene', 'Protein ID'], axis=1).iloc[0].to_list()

    # Perform linear regression analysis
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(prot_1, prot_2)
    print(p_value)  # Output the p-value

    # Create a scatter plot of protein abundance changes
    plt.figure(figsize=(9, 6))  # Set figure size
    plt.scatter(prot_1, prot_2, color='black', marker='o', s=100)  # Scatter plot

    # Generate line of best fit
    line_values = [slope * x + intercept for x in prot_1]
    plt.plot(prot_1, line_values, 'black', linewidth=3, label=f'Line of best fit: y={slope:.2f}x+{intercept:.2f}')

    # Label axes and title
    plt.xlabel(prot_1_name + ' log fold\nchange values')
    plt.ylabel(prot_2_name + ' log fold\nchange values')
    plt.title('Scatter Plot of ' + prot_1_name + ' vs ' + prot_2_name)

    # Adjust font sizes for clarity
    plt.xticks(fontsize=28)
    plt.yticks(fontsize=28)

    # Customize plot appearance
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2.5)
    ax.spines['bottom'].set_linewidth(2.5)

    # Improve layout
    plt.tight_layout()
    plt.subplots_adjust(left=0.25)

    # Save the figure to file
    png_out = 'file_output_path'
    plt.savefig(png_out)

    # Show the plot
    plt.show()

# Example function call with placeholder protein names
ppi_abundance_scatter('protein_1_name', 'protein_2_name')
