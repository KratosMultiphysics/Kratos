import json
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import argparse
import os
import sys
from datetime import datetime

# ─────────────────────────────────────────────
#  Color palette (colorblind-friendly)
# ─────────────────────────────────────────────
PALETTE = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52",
    "#8172B3", "#937860", "#DA8BC3", "#8C8C8C",
]

def load_json(filename: str) -> dict:
    """Load and validate a Google Benchmark JSON file."""
    if not os.path.exists(filename):
        print(f"[ERROR] File not found: {filename}")
        sys.exit(1)

    with open(filename, "r") as f:
        try:
            data = json.load(f)
        except json.JSONDecodeError as e:
            print(f"[ERROR] Failed to parse JSON in '{filename}': {e}")
            sys.exit(1)

    if "benchmarks" not in data:
        print(f"[ERROR] No 'benchmarks' key found in '{filename}'. "
              "Make sure the file was generated with --benchmark_format=json.")
        sys.exit(1)

    return data

def extract_context(data: dict) -> dict:
    """Extract run context metadata (host, date, CPU info, etc.)."""
    ctx = data.get("context", {})
    return {
        "date":          ctx.get("date", "N/A"),
        "host":          ctx.get("host_name", "N/A"),
        "executable":    ctx.get("executable", "N/A"),
        "num_cpus":      ctx.get("num_cpus", "N/A"),
        "mhz_per_cpu":   ctx.get("mhz_per_cpu", "N/A"),
        "cpu_scaling":   ctx.get("cpu_scaling_enabled", "N/A"),
        "library_build": ctx.get("library_build_type", "N/A"),
    }

def build_dataframe(filenames: list[str]) -> tuple[pd.DataFrame, dict]:
    """
    Load all JSON files, combine benchmark data into one DataFrame,
    and collect per-file context metadata.
    """
    combined_rows = []
    contexts      = {}

    for filename in filenames:
        data               = load_json(filename)
        contexts[filename] = extract_context(data)

        for bm in data["benchmarks"]:
            # Skip aggregate rows (mean / median / stddev)
            if bm.get("run_type", "") == "aggregate":
                continue
            bm["source"] = filename
            combined_rows.append(bm)

    if not combined_rows:
        print("[ERROR] No benchmark data found across the provided files.")
        sys.exit(1)

    df = pd.DataFrame(combined_rows)

    # Ensure numeric columns exist
    for col in ("cpu_time", "real_time", "iterations"):
        if col not in df.columns:
            df[col] = np.nan

    return df, contexts

def smart_time_unit(max_ns: float) -> tuple[float, str]:
    """
    Choose a human-friendly time unit based on the maximum value (in ns).

    Returns:
        (divisor, unit_label)
    """
    if max_ns >= 1e9:
        return 1e9, "s"
    elif max_ns >= 1e6:
        return 1e6, "ms"
    elif max_ns >= 1e3:
        return 1e3, "µs"
    else:
        return 1.0, "ns"

def print_summary_table(df: pd.DataFrame, contexts: dict, divisor: float, unit: str):
    """Print a formatted summary table to the console."""
    print("\n" + "═" * 70)
    print("  GOOGLE BENCHMARK — SUMMARY REPORT")
    print(f"  Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("═" * 70)

    for filename, ctx in contexts.items():
        print(f"\n📄 File      : {filename}")
        print(f"   Host      : {ctx['host']}")
        print(f"   Date      : {ctx['date']}")
        print(f"   CPUs      : {ctx['num_cpus']}  @  {ctx['mhz_per_cpu']} MHz")
        print(f"   Scaling   : {ctx['cpu_scaling']}")
        print(f"   Build     : {ctx['library_build']}")

    print("\n" + "─" * 70)
    print(f"  {'Benchmark':<35} {'File':<20} {'CPU Time':>12} {'Real Time':>12} {'Iterations':>12}")
    print("─" * 70)

    for _, row in df.iterrows():
        src = os.path.basename(str(row["source"]))
        print(
            f"  {str(row['name']):<35} "
            f"{src:<20} "
            f"{row['cpu_time'] / divisor:>11.3f}{unit} "
            f"{row['real_time'] / divisor:>11.3f}{unit} "
            f"{int(row['iterations']) if not pd.isna(row['iterations']) else 'N/A':>12}"
        )

    print("═" * 70 + "\n")

def plot_benchmarks(df: pd.DataFrame, filenames: list[str],
                    contexts: dict, divisor: float, unit: str,
                    output: str | None = None):
    """
    Render a grouped bar chart with:
      - CPU Time & Real Time bars per file
      - Value labels on top of each bar
      - A context info box
      - Grid lines for readability
    """
    benchmark_names = df["name"].unique()
    n_benchmarks    = len(benchmark_names)
    n_files         = len(filenames)

    bars_per_group  = 2 * n_files
    bar_width       = min(0.8 / bars_per_group, 0.25)
    group_centers   = np.arange(n_benchmarks)

    fig, ax = plt.subplots(figsize=(max(12, n_benchmarks * 2.5), 7))
    fig.patch.set_facecolor("#F8F9FA")
    ax.set_facecolor("#F8F9FA")

    bar_index = 0

    for f_idx, filename in enumerate(filenames):
        color_cpu  = PALETTE[f_idx * 2 % len(PALETTE)]
        color_real = PALETTE[(f_idx * 2 + 1) % len(PALETTE)]
        subset     = df[df["source"] == filename]

        # Align subset to benchmark_names order
        subset_indexed = (
            subset.set_index("name")
                  .reindex(benchmark_names)
        )

        cpu_vals  = subset_indexed["cpu_time"].values  / divisor
        real_vals = subset_indexed["real_time"].values / divisor

        short_name = os.path.basename(filename)

        # ── CPU Time bars ──────────────────────────────────────────────
        offset_cpu = (bar_index - bars_per_group / 2 + 0.5) * bar_width
        bars_cpu   = ax.bar(
            group_centers + offset_cpu,
            cpu_vals,
            bar_width,
            label=f"{short_name} — CPU Time",
            color=color_cpu,
            alpha=0.85,
            edgecolor="white",
            linewidth=0.6,
            zorder=3,
        )
        bar_index += 1

        # ── Real Time bars ─────────────────────────────────────────────
        offset_real = (bar_index - bars_per_group / 2 + 0.5) * bar_width
        bars_real   = ax.bar(
            group_centers + offset_real,
            real_vals,
            bar_width,
            label=f"{short_name} — Real Time",
            color=color_real,
            alpha=0.85,
            edgecolor="white",
            linewidth=0.6,
            zorder=3,
            hatch="//",
        )
        bar_index += 1

        # ── Value labels on top of each bar ───────────────────────────
        for bars in (bars_cpu, bars_real):
            for bar in bars:
                h = bar.get_height()
                if not np.isnan(h) and h > 0:
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        h * 1.015,
                        f"{h:.2f}",
                        ha="center", va="bottom",
                        fontsize=7, color="#333333",
                        rotation=45,
                    )

    # ── Axes formatting ────────────────────────────────────────────────
    ax.set_title(
        "Google Benchmark — Performance Comparison",
        fontsize=15, fontweight="bold", pad=18, color="#1A1A2E"
    )
    ax.set_xlabel("Benchmark Name", fontsize=11, labelpad=10, color="#333333")
    ax.set_ylabel(f"Time ({unit})",  fontsize=11, labelpad=10, color="#333333")

    ax.set_xticks(group_centers)
    ax.set_xticklabels(benchmark_names, rotation=40, ha="right", fontsize=9)
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter(f"%.2f {unit}"))

    ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)

    # ── Legend ─────────────────────────────────────────────────────────
    ax.legend(
        loc="upper right",
        fontsize=8,
        framealpha=0.9,
        edgecolor="#CCCCCC",
        title="Source File — Metric",
        title_fontsize=8,
    )

    # ── Context info box (bottom-left) ─────────────────────────────────
    info_lines = [f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"]
    for filename, ctx in contexts.items():
        short = os.path.basename(filename)
        info_lines.append(
            f"{short}: {ctx['num_cpus']} CPUs @ {ctx['mhz_per_cpu']} MHz | "
            f"Host: {ctx['host']} | Build: {ctx['library_build']}"
        )

    ax.text(
        0.01, 0.01, "\n".join(info_lines),
        transform=ax.transAxes,
        fontsize=6.5, color="#555555",
        verticalalignment="bottom",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                  edgecolor="#CCCCCC", alpha=0.8),
    )

    plt.tight_layout()

    if output:
        fig.savefig(output, dpi=180, bbox_inches="tight")
        print(f"[INFO] Plot saved to: {output}")
    else:
        plt.show()

def main(filenames: list[str], output: str | None = None):
    """
    Orchestrate loading, processing, reporting, and plotting
    of Google Benchmark JSON results.
    """
    df, contexts = build_dataframe(filenames)

    # Pick the best time unit automatically
    max_ns        = max(df["cpu_time"].max(), df["real_time"].max())
    divisor, unit = smart_time_unit(max_ns)

    print_summary_table(df, contexts, divisor, unit)
    plot_benchmarks(df, filenames, contexts, divisor, unit, output)

if __name__ == "__main__":
    """
    Example usage:
        python generate_plot_google_benchmark.py \\
            --filenames result1.json result2.json \\
            --output comparison.png
    """
    parser = argparse.ArgumentParser(
        description="Generate an enhanced comparison plot from Google Benchmark JSON results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  python generate_plot_google_benchmark.py "
            "--filenames a.json b.json --output out.png\n\n"
            "JSON files must be produced with:\n"
            "  ./benchmark_binary --benchmark_format=json --benchmark_out=result.json"
        ),
    )
    parser.add_argument(
        "--filenames", nargs="+", required=True,
        help="One or more Google Benchmark JSON result files."
    )
    parser.add_argument(
        "--output", type=str, default=None,
        help="Optional path to save the plot (e.g. chart.png). "
             "If omitted, the plot is displayed interactively."
    )
    args = parser.parse_args()
    main(args.filenames, args.output)