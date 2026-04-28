import json
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import argparse
import os
import sys
from datetime import datetime

plt.style.use("seaborn-v0_8-whitegrid")

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

def _truncate(name: str, max_len: int = 30) -> str:
    """Truncate a benchmark name with ellipsis if too long."""
    return name if len(name) <= max_len else name[:max_len - 1] + "…"


def plot_benchmarks(df: pd.DataFrame, filenames: list[str],
                    contexts: dict, divisor: float, unit: str,
                    output: str | None = None):
    """
    Render a grouped bar chart with:
      - CPU Time & Real Time bars per file
      - Value labels on top of each bar (green = fastest, red = slowest)
      - Bar shadows for depth
      - Log-scale Y-axis
      - Benchmarks sorted by average CPU time (descending)
      - Minimum-value reference line per group
      - Iteration count annotation below each group
      - A context info box
    """
    # ── Sort benchmarks by mean cpu_time descending (slowest first) ────
    mean_cpu = df.groupby("name")["cpu_time"].mean()
    benchmark_names = mean_cpu.sort_values(ascending=False).index.tolist()
    n_benchmarks    = len(benchmark_names)
    n_files         = len(filenames)

    bars_per_group  = 2 * n_files
    bar_width       = min(0.8 / bars_per_group, 0.25)
    group_centers   = np.arange(n_benchmarks)

    fig_height = max(7, 5 + n_files * 0.4)
    fig, ax = plt.subplots(figsize=(max(18, n_benchmarks * 3.0), fig_height))
    fig.patch.set_facecolor("#F8F9FA")
    ax.set_facecolor("#F8F9FA")

    # Pre-collect all bar heights per benchmark group to color-code labels
    all_heights: dict[str, list[float]] = {name: [] for name in benchmark_names}

    bar_index = 0
    # Store (bar_container, x_positions) for label coloring in second pass
    bar_records: list[tuple] = []

    for f_idx, filename in enumerate(filenames):
        color_cpu  = PALETTE[f_idx * 2 % len(PALETTE)]
        color_real = PALETTE[(f_idx * 2 + 1) % len(PALETTE)]
        subset     = df[df["source"] == filename]

        subset_indexed = (
            subset.set_index("name")
                  .reindex(benchmark_names)
        )

        cpu_vals  = subset_indexed["cpu_time"].values  / divisor
        real_vals = subset_indexed["real_time"].values / divisor
        iters     = subset_indexed["iterations"].values

        short_name = os.path.basename(filename)

        # ── CPU Time bars ──────────────────────────────────────────────
        offset_cpu = (bar_index - bars_per_group / 2 + 0.5) * bar_width
        x_cpu      = group_centers + offset_cpu
        # Shadow
        ax.bar(x_cpu + 0.008, cpu_vals, bar_width,
               color="#000000", alpha=0.08, zorder=2)
        bars_cpu = ax.bar(
            x_cpu, cpu_vals, bar_width,
            label=f"{short_name} — CPU Time",
            color=color_cpu, alpha=0.88,
            edgecolor="white", linewidth=0.6, zorder=3,
        )
        bar_index += 1

        # ── Real Time bars ─────────────────────────────────────────────
        offset_real = (bar_index - bars_per_group / 2 + 0.5) * bar_width
        x_real      = group_centers + offset_real
        ax.bar(x_real + 0.008, real_vals, bar_width,
               color="#000000", alpha=0.08, zorder=2)
        bars_real = ax.bar(
            x_real, real_vals, bar_width,
            label=f"{short_name} — Real Time",
            color=color_real, alpha=0.88,
            edgecolor="white", linewidth=0.6, zorder=3, hatch="//",
        )
        bar_index += 1

        # Accumulate heights per benchmark for min/max detection
        for i, name in enumerate(benchmark_names):
            for h in (cpu_vals[i], real_vals[i]):
                if not np.isnan(h) and h > 0:
                    all_heights[name].append(h)

        bar_records.append((bars_cpu,  cpu_vals))
        bar_records.append((bars_real, real_vals))

        # ── Iteration count annotation (below group, once per file) ───
        for i, (name, it) in enumerate(zip(benchmark_names, iters)):
            if not np.isnan(it):
                ax.text(
                    group_centers[i], -0.06, f"{int(it):,} iters",
                    transform=ax.get_xaxis_transform(),
                    ha="center", va="top",
                    fontsize=5.5, color="#888888",
                )

    # ── Color-coded value labels (green=fastest, red=slowest) ──────────
    global_min = min(v for vals in all_heights.values() for v in vals) if all_heights else 0
    global_max = max(v for vals in all_heights.values() for v in vals) if all_heights else 1

    for bars_container, vals in bar_records:
        for bar, h in zip(bars_container, vals):
            if np.isnan(h) or h <= 0:
                continue
            if h == global_min:
                lbl_color = "#2ca02c"   # green  — fastest
            elif h == global_max:
                lbl_color = "#d62728"   # red    — slowest
            else:
                lbl_color = "#444444"   # neutral
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                h * 1.015,
                f"{h:.2f}",
                ha="center", va="bottom",
                fontsize=7, color=lbl_color,
                fontweight="bold" if lbl_color != "#444444" else "normal",
                rotation=45, zorder=5,
            )

    # ── Minimum reference line per benchmark group ─────────────────────
    for i, name in enumerate(benchmark_names):
        if all_heights[name]:
            min_h = min(all_heights[name])
            half  = bar_width * bars_per_group / 2
            ax.hlines(min_h,
                      group_centers[i] - half - 0.02,
                      group_centers[i] + half + 0.02,
                      colors="#2ca02c", linewidths=1.2,
                      linestyles="dashed", zorder=4, alpha=0.7)

    # ── Log scale Y-axis ───────────────────────────────────────────────
    pos_vals = [h for vals in all_heights.values() for h in vals if h > 0]
    if pos_vals and max(pos_vals) / max(min(pos_vals), 1e-12) > 10:
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f"{x:.3g} {unit}")
        )
    else:
        ax.yaxis.set_major_formatter(mticker.FormatStrFormatter(f"%.2f {unit}"))

    # ── Axes formatting ────────────────────────────────────────────────
    ax.set_title(
        "Google Benchmark — Performance Comparison",
        fontsize=15, fontweight="bold", pad=22, color="#1A1A2E"
    )
    # Accent underline beneath title
    ax.annotate("", xy=(1, 1.038), xytext=(0, 1.038),
                xycoords="axes fraction", textcoords="axes fraction",
                arrowprops=dict(arrowstyle="-", color="#4C72B0",
                                lw=2, alpha=0.35))

    ax.set_xlabel("Benchmark Name", fontsize=11, labelpad=14, color="#333333")
    ax.set_ylabel(f"Time ({unit})",  fontsize=11, labelpad=10, color="#333333")

    truncated = [_truncate(n) for n in benchmark_names]
    ax.set_xticks(group_centers)
    ax.set_xticklabels(truncated, rotation=40, ha="right", fontsize=9)

    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)

    # ── Legend (outside the plot, to the right) ──────────────────────
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0,
        fontsize=8,
        framealpha=0.9,
        edgecolor="#CCCCCC",
        title="Source File — Metric",
        title_fontsize=8,
    )

    # ── Context info box (right side, below the legend) ────────────────
    info_lines = [f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"]
    for filename, ctx in contexts.items():
        short = os.path.basename(filename)
        info_lines.append(
            f"{short}:  {ctx['num_cpus']} CPUs @ {ctx['mhz_per_cpu']} MHz  |  "
            f"Host: {ctx['host']}  |  Build: {ctx['library_build']}"
        )

    ax.text(
        1.01, 0.30, "\n\n".join(info_lines),
        transform=ax.transAxes,
        fontsize=5.0, color="#555555",
        verticalalignment="top",
        horizontalalignment="left",
        wrap=False,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                  edgecolor="#CCCCCC", alpha=0.8),
    )


    plt.tight_layout(rect=[0, 0, 0.78, 1])

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

# ─────────────────────────────────────────────────────────────────────────────
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