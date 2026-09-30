#!/usr/bin/env python3
"""Optidose 96-well growth curves for the three CRISPRi knockdown libraries.

Reads the plate-reader exports for E. coli, E. cloacae and K. pneumoniae,
applies the shared plate map, and draws one small-multiples grid:
rows = organism, columns = imipenem dose, solid = induced, dashed = uninduced.

Usage:  python3 Sequencing/plot_optidose.py
Writes: Sequencing/optidose_96well.png
"""

import csv
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent

# --- plate map (shared by all three plates) -------------------------------
# Rows A/H and columns 1/12 are media blanks; B,C,D are uninduced replicates
# and E,F,G are induced replicates. Column 2 is the no-drug control; columns
# 3..11 are a 2-fold series running 0.125 -> 32 ug/mL.
UNINDUCED_ROWS = ("B", "C", "D")
INDUCED_ROWS = ("E", "F", "G")
DOSES = [0.0, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
DOSE_BY_COLUMN = {2: 0.0} | {col: 0.125 * 2 ** (col - 3) for col in range(3, 12)}

ORGANISMS = [
    ("E. coli", "IMI_Ecoli/ECimiRDW-growth-data.tsv", "#2a78d6"),
    ("E. cloacae", "IMI_Ecloacae/growth-data-2.tsv", "#eb6834"),
    ("K. pneumoniae", "IMI_Kpneumoniae/growth-data2A.tsv", "#1baf7a"),
]

INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#e6e5e1"


def read_plate(path):
    """Return (hours, {well: OD600 array}) from either plate-reader layout.

    Two export shapes appear in this repo: timepoints-as-rows with a leading
    'Seconds' column, and the transposed 'Cycle Nr.' form with wells as rows.
    """
    with open(path, newline="") as handle:
        rows = [row for row in csv.reader(handle, delimiter="\t") if row]

    if rows[0][0].strip() == "Seconds":  # timepoints as rows
        wells = rows[0][2:]
        seconds = np.array([float(r[0].rstrip("s")) for r in rows[1:]])
        cols = list(zip(*(r[2:] for r in rows[1:])))
        traces = {w: np.array([float(v) for v in col]) for w, col in zip(wells, cols)}
    else:  # transposed: 'Cycle Nr.' header, 'Time [s]' row, then one row per well
        by_label = {r[0].strip(): r[1:] for r in rows}
        seconds = np.array([float(v) for v in by_label["Time [s]"]])
        traces = {
            label: np.array([float(v) for v in values])
            for label, values in by_label.items()
            if re.fullmatch(r"[A-H]\d{1,2}", label)
        }

    return seconds / 3600.0, traces


def replicate_stats(traces, rows, column):
    """Mean and standard error across the replicate rows of one dose column."""
    stack = np.vstack([traces[f"{row}{column}"] for row in rows if f"{row}{column}" in traces])
    mean = stack.mean(axis=0)
    sem = stack.std(axis=0, ddof=1) / np.sqrt(stack.shape[0]) if stack.shape[0] > 1 else None
    return mean, sem


def main():
    plates = [(name, color, *read_plate(ROOT / rel)) for name, rel, color in ORGANISMS]

    fig, axes = plt.subplots(
        len(plates), len(DOSES), figsize=(13.5, 5.4), sharex=True, sharey=True
    )
    fig.patch.set_facecolor("#fcfcfb")

    dose_by_column = {dose: col for col, dose in DOSE_BY_COLUMN.items()}

    for r, (name, color, hours, traces) in enumerate(plates):
        for c, dose in enumerate(DOSES):
            ax = axes[r, c]
            column = dose_by_column[dose]
            for style, rows in (("--", UNINDUCED_ROWS), ("-", INDUCED_ROWS)):
                mean, sem = replicate_stats(traces, rows, column)
                if sem is not None:
                    ax.fill_between(
                        hours, mean - sem, mean + sem, color=color, alpha=0.25, linewidth=0
                    )
                ax.plot(hours, mean, style, color=color, linewidth=1.6, solid_capstyle="round")

            ax.set_yscale("log")
            ax.set_ylim(0.07, 1.4)
            ax.set_xlim(0, max(hours))
            ax.grid(True, which="major", color=GRID, linewidth=0.6)
            ax.set_axisbelow(True)
            ax.tick_params(which="both", length=0, labelsize=8, colors=INK_MUTED)
            ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
            ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([0.1, 0.3, 1.0]))
            ax.yaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(["0.1", "0.3", "1.0"]))
            for side in ax.spines.values():
                side.set_visible(False)

            if r == 0:
                ax.set_title(
                    f"{dose:g}" if dose else "0",
                    fontsize=9,
                    color=INK,
                    pad=6,
                )
            if c == 0:
                ax.set_ylabel(name, fontsize=9.5, style="italic", color=INK, labelpad=8)
            if r == len(plates) - 1:
                ax.set_xticks([0, 6, 12, 18])

    fig.suptitle(
        "Imipenem dose response of the CRISPRi knockdown libraries",
        x=0.012,
        y=0.98,
        ha="left",
        fontsize=13,
        color=INK,
    )
    fig.text(
        0.012,
        0.915,
        "Imipenem (µg/mL)",
        ha="left",
        fontsize=9.5,
        color=INK_MUTED,
    )
    fig.text(0.5, 0.045, "Time (h)", ha="center", fontsize=9.5, color=INK_MUTED)
    fig.text(
        0.006,
        0.5,
        "OD600 (log scale)",
        va="center",
        rotation=90,
        fontsize=9.5,
        color=INK_MUTED,
    )

    fig.legend(
        handles=[
            Line2D([], [], color=INK_MUTED, linestyle="-", linewidth=1.6, label="induced"),
            Line2D([], [], color=INK_MUTED, linestyle="--", linewidth=1.6, label="uninduced"),
        ],
        loc="upper right",
        bbox_to_anchor=(0.995, 1.0),
        frameon=False,
        ncols=2,
        fontsize=9,
        labelcolor=INK_MUTED,
        handlelength=2.0,
    )

    fig.subplots_adjust(left=0.075, right=0.995, top=0.855, bottom=0.11, wspace=0.14, hspace=0.18)
    out = ROOT / "optidose_96well.png"
    fig.savefig(out, dpi=200, facecolor=fig.get_facecolor())
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
