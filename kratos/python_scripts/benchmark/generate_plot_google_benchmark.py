import json
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def main(filenames):
    """
    Main function to generate a combined plot from multiple Google Benchmark JSON outputs.

    Args:
        filenames (list of str): List of paths to the JSON files containing benchmark results.

    The function performs the following steps:
    1. Loads the JSON files specified in the filenames list.
    2. Combines the 'benchmarks' section of all JSON data into a single pandas DataFrame.
    3. Plots the 'cpu_time' and 'real_time' for each benchmark with bars grouped by input filename.

    The resulting plot displays:
    - X-axis: Benchmark names (with labels based on input filenames)
    - Y-axis: Time in nanoseconds
    - Bars grouped by CPU Time and Real Time for each file

    NOTE: The JSON files must be generated using Google Benchmark with:
    ./benchmark_name --benchmark_format=json --benchmark_out=filename.json
    """
    combined_data = []

    # Process each JSON file
    for filename in filenames:
        with open(filename) as f:
            data = json.load(f)

        # Extract benchmark information and add filename as a source
        benchmarks = data["benchmarks"]
        for benchmark in benchmarks:
            benchmark["source"] = filename
        combined_data.extend(benchmarks)

    # Convert combined data to a pandas DataFrame
    benchmark_df = pd.DataFrame(combined_data)

    # Plot the results
    plt.figure(figsize=(12, 8))
    bar_width = 0.35
    x_labels = benchmark_df["name"].unique()
    x = range(len(x_labels))
    offset = 0

    # Plot for each input file
    for filename in filenames:
        subset = benchmark_df[benchmark_df["source"] == filename]
        plt.bar(
            [pos + offset for pos in x],
            subset["real_time"],
            bar_width,
            label=f"{filename} - Real Time",
            alpha=0.7
        )
        plt.bar(
            [pos + offset for pos in x],
            subset["cpu_time"],
            bar_width,
            label=f"{filename} - CPU Time",
            alpha=0.7
        )
        offset += bar_width

    # Customize the plot
    plt.title("Benchmark Performance Across Files")
    plt.xlabel("Benchmark Name")
    plt.ylabel("Time (ns)")
    plt.xticks([pos + bar_width for pos in x], x_labels, rotation=45, ha="right")
    plt.legend()
    plt.tight_layout()

    # Show the plot
    plt.show()

if __name__ == "__main__":
    """
    This script can be run from the command line to generate a combined plot from multiple Google Benchmark JSON outputs.
    Example usage:
    python generate_plot_google_benchmark.py --filenames filename1.json filename2.json
    """
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Generate a combined plot from multiple Google Benchmark JSON results.')
    parser.add_argument('--filenames', nargs='+', type=str, help='The list of JSON files containing benchmark results')
    args = parser.parse_args()
    main(args.filenames)