import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

print("Calculating distances")
def calculate_distances_per_sample(data, gene1, gene2, circular_data=None, all_circular=False):
    data.columns = data.columns.str.strip()
    if 'GENE' not in data.columns:
        raise ValueError("The input file does not contain a 'GENE' column.")

    filtered_data = data[data['GENE'].isin([gene1, gene2])]
    filtered_data = filtered_data.sort_values(by=['SEQUENCE', 'START'])

    sequence_distances = {}

    # Load circularity data if provided
    if circular_data is not None:
        circular_dict = dict(zip(circular_data['SEQUENCE'], circular_data['CIRCULAR']))

    for sequence, group in filtered_data.groupby('SEQUENCE'):
        group = group.reset_index(drop=True)
        distances = []

        gene1_indices = group[group['GENE'] == gene1].index
        gene2_indices = group[group['GENE'] == gene2].index

        for i in gene1_indices:
            for j in gene2_indices:
                # Normal distance calculation
                if group.loc[i, 'START'] < group.loc[j, 'START']:
                    distance = group.loc[j, 'START'] - group.loc[i, 'END']
                else:
                    distance = group.loc[i, 'START'] - group.loc[j, 'END']

                # Check if the sequence is circular
                if circular_data is not None:
                    is_circular = circular_dict.get(sequence, 0) == 1 or all_circular
                    if is_circular:
                        total_length = circular_data[circular_data['SEQUENCE'] == sequence]['LENGTH'].values[0]
                        if distance > total_length / 2:
                            # Adjust the distance when it exceeds half of the total length
                            if group.loc[i, 'START'] < group.loc[j, 'START']:
                                # Gene i comes first
                                distance = (total_length + group.loc[i, 'START']) - group.loc[j, 'END']
                            else:
                                # Gene j comes first
                                distance = (total_length + group.loc[j, 'START']) - group.loc[i, 'END']

                distances.append(distance)

        if distances:
            sequence_distances[sequence] = min(distances)

    return sequence_distances

def calculate_summary_stats(sequence_distances):
    distances = list(sequence_distances.values())
    if distances:
        median_distance = np.median(distances)
        min_distance = min(distances)
        max_distance = max(distances)
        mode_result = stats.mode(distances, keepdims=True)
        mode_distance = mode_result.mode[0] if mode_result.mode.size > 0 else None
        return min_distance, median_distance, max_distance, mode_distance
    else:
        return None, None, None, None

def plot_histogram(sequence_distances, median_distance, output_pdf):
    distances = list(sequence_distances.values())
    if distances:
        plt.figure(figsize=(8, 6))
        sns.histplot(distances, bins=30, color="#525252", alpha=1)

        # Only plot the median
        plt.axvline(median_distance, color='black', linestyle='dashed', linewidth=1.2, label=f'Median: {median_distance:.2f} bp')

        plt.title("Genomic Distance Distribution", fontsize=16)
        plt.xlabel("Distance (bp)", fontsize=12)
        plt.ylabel("Frequency", fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()

        plt.text(median_distance, plt.ylim()[1] * 0.9, f"Median: {median_distance:.2f}", color="black", rotation=90, va="bottom")

        plt.savefig(output_pdf, dpi=600, bbox_inches='tight')
        plt.close()

def save_results_to_csv(sequence_distances, output_csv):
    df = pd.DataFrame(list(sequence_distances.items()), columns=["SEQUENCE", "Minimum Distance"])
    df.to_csv(output_csv, index=False)
    print(f"Results saved to '{output_csv}'")

def main():
    parser = argparse.ArgumentParser(description='Calculate genomic distances between two genes.')
    parser.add_argument('-f', '--file', type=str, required=True, help='The input file with genomic data (TSV format).')
    parser.add_argument('-g1', '--gene1', type=str, required=True, help='The first gene to compare.')
    parser.add_argument('-g2', '--gene2', type=str, required=True, help='The second gene to compare.')
    parser.add_argument('--output', type=str, default="genomic_distances", help='Base name for output files (CSV and PDF).')
    parser.add_argument('--circular', type=str, help='File containing sequence names, topology (0 or 1), and length for circular sequences.')
    parser.add_argument('--all-circular', action='store_true', help='Assume all sequences are circular.')
    args = parser.parse_args()

    data = pd.read_csv(args.file, sep="\t")  # Updated to read TSV
    circular_data = None
    if args.circular:
        circular_data = pd.read_csv(args.circular, sep="\t", names=["SEQUENCE", "CIRCULAR", "LENGTH"])

    distances_per_sequence = calculate_distances_per_sample(data, args.gene1, args.gene2, circular_data, args.all_circular)
    output_csv = f"{args.output}_results.csv"
    save_results_to_csv(distances_per_sequence, output_csv)

    min_distance, median_distance, max_distance, mode_distance = calculate_summary_stats(distances_per_sequence)
    if min_distance is not None:
        print(f"Min: {min_distance:.2f} bp, Median: {median_distance:.2f} bp, Max: {max_distance:.2f} bp, Mode: {mode_distance}")

    output_pdf = f"{args.output}_histogram.pdf"
    plot_histogram(distances_per_sequence, median_distance, output_pdf)

if __name__ == '__main__':
    main()
