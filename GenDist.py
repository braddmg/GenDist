import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def calculate_distances_per_sample(data, gene1, gene2):
    # Strip extra spaces from column names
    data.columns = data.columns.str.strip()

    # Check if the 'gene' column exists
    if 'gene' not in data.columns:
        raise ValueError("The input file does not contain a 'gene' column.")

    # Filter the data for the two specified genes
    filtered_data = data[data['gene'].isin([gene1, gene2])]

    # Sort data by Sample and start position
    filtered_data = filtered_data.sort_values(by=['Sample', 'start'])

    # Initialize a dictionary to store distances by sample
    sample_distances = {}

    # Group by Sample
    for sample, group in filtered_data.groupby('Sample'):
        group = group.reset_index(drop=True)
        distances = []

        # Find all occurrences of gene1 and gene2
        gene1_indices = group[group['gene'] == gene1].index
        gene2_indices = group[group['gene'] == gene2].index

        # Ensure both genes are ordered in the same sample
        for i in gene1_indices:
            for j in gene2_indices:
                # Calculate distance by subtracting the end of one gene from the start of the other
                if group.loc[i, 'start'] < group.loc[j, 'start']:
                    distance = group.loc[j, 'start'] - group.loc[i, 'end']
                else:
                    distance = group.loc[i, 'start'] - group.loc[j, 'end']

                distances.append(distance)

        # Store the minimum distance for this sample
        if distances:
            sample_distances[sample] = min(distances)

    return sample_distances

def calculate_summary_stats(sample_distances):
    # Calculate summary statistics of the minimum distances across all samples
    distances = list(sample_distances.values())

    if distances:
        median_distance = np.median(distances)
        min_distance = min(distances)
        max_distance = max(distances)

        # Compute mode properly
        mode_result = stats.mode(distances, keepdims=True)
        mode_distance = mode_result.mode[0] if mode_result.mode.size > 0 else None

        return min_distance, median_distance, max_distance, mode_distance
    else:
        return None, None, None, None

def plot_histogram(sample_distances, min_distance, median_distance, max_distance, output_image):
    # Plot KDE plot (smooth histogram) for the minimum distances
    distances = list(sample_distances.values())
    if distances:
        plt.figure(figsize=(10, 6))

        # KDE Plot for smooth distribution
        sns.kdeplot(distances, shade=True, color="#525252", alpha=0.5)

        # Add vertical lines for min, median, and max values
        if min_distance is not None:
            plt.axvline(min_distance, color='blue', linestyle='solid', linewidth=2, label=f'Min: {min_distance:.2f} bp')
        if median_distance is not None:
            plt.axvline(median_distance, color='red', linestyle='dashed', linewidth=2, label=f'Median: {median_distance:.2f} bp')
        if max_distance is not None:
            plt.axvline(max_distance, color='green', linestyle='solid', linewidth=2, label=f'Max: {max_distance:.2f} bp')

        # Title and labels
        plt.title("Smooth Distribution of Minimum Genomic Distances", fontsize=16)
        plt.xlabel("Distance (bp)", fontsize=12)
        plt.ylabel("Density", fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.7)

        # Add legend
        plt.legend()

        # Save the plot as an image
        plt.savefig(output_image, dpi=300, bbox_inches='tight')

        # Show the plot
        plt.show()

def save_results_to_csv(sample_distances, output_csv):
    # Convert dictionary to DataFrame
    df = pd.DataFrame(list(sample_distances.items()), columns=["Sample", "Minimum Distance"])
    
    # Save as CSV file
    df.to_csv(output_csv, index=False)

    print(f"\nResults saved to '{output_csv}'")

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description='Calculate genomic distances between two genes.')
    parser.add_argument('file', type=str, help='The input file with genomic data (CSV format).')
    parser.add_argument('gene1', type=str, help='The first gene to compare.')
    parser.add_argument('gene2', type=str, help='The second gene to compare.')
    parser.add_argument('--output', type=str, default="genomic_distances", help='Base name for output files (CSV and image).')

    # Parse the arguments
    args = parser.parse_args()

    # Load the data (for CSV with commas)
    data = pd.read_csv(args.file, sep=",")  # Changed sep to ',' for CSV

    # Print the column names to debug
    print(f"Columns in the input file: {data.columns}")

    # Calculate distances for each sample
    distances_per_sample = calculate_distances_per_sample(data, args.gene1, args.gene2)

    # Print the distances for each sample
    for sample, distance in distances_per_sample.items():
        print(f"Sample: {sample}")
        print(f"Minimum Distance: {distance} bp")

    # Save results to a CSV file
    output_csv = f"{args.output}_results.csv"
    save_results_to_csv(distances_per_sample, output_csv)

    # Calculate the min, median, max, and mode distances
    min_distance, median_distance, max_distance, mode_distance = calculate_summary_stats(distances_per_sample)

    # Print the summary statistics
    if min_distance is not None and median_distance is not None and max_distance is not None:
        print(f"\nThe minimum genomic distance is {min_distance:.2f} bp.")
        print(f"The median of the minimum genomic distances is {median_distance:.2f} bp.")
        print(f"The maximum genomic distance is {max_distance:.2f} bp.")
        print(f"The mode of the minimum genomic distances is {mode_distance} bp.")
    else:
        print("\nNo distances were calculated.")

    # Plot the histogram and save it
    output_image = f"{args.output}_histogram.png"
    plot_histogram(distances_per_sample, min_distance, median_distance, max_distance, output_image)

if __name__ == '__main__':
    main()
