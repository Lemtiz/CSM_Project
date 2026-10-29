"""Plot VV/r-NAPA execution times and speedup from benchmark.csv."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterMathtext, LogLocator


REQUIRED_COLUMNS = {
    "k",
    "vv_seconds",
    "napa_seconds",
    "speedup",
}


def load_results(path: Path) -> dict[str, list[float]]:
    if not path.is_file():
        raise FileNotFoundError(f"Benchmark file not found: {path}")

    with path.open(newline="", encoding="utf-8") as file:
        reader = csv.DictReader(file)
        columns = set(reader.fieldnames or [])
        missing = REQUIRED_COLUMNS.difference(columns)
        if missing:
            missing_columns = ", ".join(sorted(missing))
            raise ValueError(f"Missing benchmark columns: {missing_columns}")
        rows = list(reader)

    if not rows:
        raise ValueError(f"Benchmark file is empty: {path}")

    results = {
        column: [float(row[column]) for row in rows]
        for column in REQUIRED_COLUMNS
    }
    if any(value <= 0.0 for value in results["k"]):
        raise ValueError("All k values must be positive for the logarithmic axis")

    order = sorted(range(len(results["k"])), key=results["k"].__getitem__)
    return {
        column: [values[index] for index in order]
        for column, values in results.items()
    }


def plot_results(results: dict[str, list[float]], output: Path) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update(
        {
            "axes.titlesize": 17,
            "axes.labelsize": 15,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "legend.fontsize": 13,
        }
    )
    figure, axes = plt.subplots(
        1,
        2,
        figsize=(12, 5.2),
        constrained_layout=True,
    )

    axes[0].plot(
        results["k"],
        results["vv_seconds"],
        marker="o",
        linewidth=2,
        markersize=5,
        color="#202020",
        label="Velocity Verlet",
    )
    axes[0].plot(
        results["k"],
        results["napa_seconds"],
        marker="s",
        linewidth=2,
        markersize=5,
        color="#777777",
        label="r-NAPA",
    )
    axes[0].set_xscale("log")
    axes[0].set_xlabel(r"$k\ [\varepsilon/\sigma^2]$")
    axes[0].set_ylabel("Execution time [s]")
    axes[0].legend(frameon=True)

    axes[1].plot(
        results["k"],
        results["speedup"],
        marker="o",
        linewidth=2,
        markersize=5,
        color="#404040",
        label="Measured speedup",
    )
    log_k = [math.log(value) for value in results["k"]]
    log_speedup = [math.log(value) for value in results["speedup"]]
    mean_log_k = sum(log_k) / len(log_k)
    mean_log_speedup = sum(log_speedup) / len(log_speedup)
    denominator = sum(
        (value - mean_log_k) ** 2 for value in log_k
    )
    exponent = sum(
        (x - mean_log_k) * (y - mean_log_speedup)
        for x, y in zip(log_k, log_speedup)
    ) / denominator
    intercept = mean_log_speedup - exponent * mean_log_k
    residuals = [
        y - (intercept + exponent * x)
        for x, y in zip(log_k, log_speedup)
    ]
    degrees_of_freedom = len(log_k) - 2
    exponent_error = math.sqrt(
        sum(residual ** 2 for residual in residuals)
        / degrees_of_freedom
        / denominator
    )
    fit_k_min = min(results["k"])
    fit_k_max = max(results["k"])
    fit_k = [
        fit_k_min * (fit_k_max / fit_k_min) ** (index / 199)
        for index in range(200)
    ]
    fit_speedup = [
        math.exp(intercept) * value ** exponent
        for value in fit_k
    ]
    axes[1].plot(
        fit_k,
        fit_speedup,
        color="#909090",
        linewidth=1.8,
        linestyle="--",
        label=rf"Fit $S \propto k^a$, $a={exponent:.3f}\pm{exponent_error:.3f}$",
    )
    axes[1].axhline(
        1.0,
        color="#909090",
        linewidth=1.2,
        linestyle="dotted",
    )
    axes[1].set_xscale("log")
    axes[1].xaxis.set_major_locator(
        LogLocator(base=10.0, subs=(1.0,))
    )
    axes[1].xaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
    speedup_tick_max = math.ceil(max(results["speedup"]))
    axes[1].set_yticks(
        list(range(1, speedup_tick_max + 1))
    )
    axes[1].set_xlabel(r"$k\ [\varepsilon/\sigma^2]$")
    axes[1].set_ylabel(r"Speedup $t_{\mathrm{VV}}/t_{\mathrm{r\!-\!NAPA}}$")
    axes[1].legend(frameon=True)

    for axis in axes:
        axis.set_axisbelow(True)
        axis.grid(True, which="major", linewidth=0.8, alpha=0.7)
        axis.grid(True, which="minor", linewidth=0.4, alpha=0.35)

    figure.savefig(output, dpi=220, bbox_inches="tight")
    figure.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    print(f"Saved {output}")
    print(f"Saved {output.with_suffix('.pdf')}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot benchmark execution times and r-NAPA speedup."
    )
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=Path("benchmark.csv"),
        help="CSV generated by the C++ benchmark",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("benchmark_plots.png"),
        help="PNG output path",
    )
    args = parser.parse_args()

    results = load_results(args.input)
    plot_results(results, args.output)


if __name__ == "__main__":
    main()
