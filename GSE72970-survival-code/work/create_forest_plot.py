#!/usr/bin/env python3
"""Render the requested Cox OS forest plot when R is unavailable."""

from pathlib import Path
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

FONT = FontProperties(family="Helvetica")
FONT_BOLD = FontProperties(family="Helvetica", weight="bold")
RED = "#A63D40"
TEXT = "#262626"
GRAY = "#666666"

matplotlib.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "font.family": "Helvetica",
    "font.sans-serif": ["Helvetica", "Arial"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.linewidth": 0.65,
})

fig, ax = plt.subplots(figsize=(3.2, 2.3))

# Reference line and estimates.
ax.axvline(1.0, color="#808080", linestyle=(0, (4, 3)), linewidth=0.65, zorder=1)
ax.errorbar(
    1.724, 1, xerr=[[1.724 - 1.138], [2.612 - 1.724]],
    fmt="o", color=RED, ecolor=RED, elinewidth=1.15,
    capsize=0, markersize=4.5, markeredgewidth=0, zorder=3,
)
ax.plot(1.0, 0, "o", color="#555555", markersize=3.8, markeredgewidth=0, zorder=3)

# Log10 HR axis.
ax.set_xscale("log", base=10)
ax.set_xlim(0.7, 3.2)
ax.set_ylim(-0.55, 1.55)
ticks = [0.75, 1.0, 1.5, 2.0, 3.0]
ax.xaxis.set_major_locator(FixedLocator(ticks))
ax.xaxis.set_major_formatter(FixedFormatter(["0.75", "1", "1.5", "2", "3"]))
ax.xaxis.set_minor_formatter(NullFormatter())

ax.set_yticks([0, 1], labels=["Low", "High"], fontproperties=FONT, fontsize=8.5)
ax.set_xlabel("Hazard ratio (95% CI)", fontproperties=FONT, fontsize=8.5, labelpad=5)
ax.set_title("Overall survival", loc="left", fontproperties=FONT_BOLD, fontsize=10, pad=8)

# Exact requested model annotation, without significance stars.
ax.text(
    0.99, 1.02, "High: HR 1.724 (1.138–2.612)   P = 0.010",
    transform=ax.transAxes, ha="right", va="bottom",
    fontproperties=FONT, fontsize=6.3, color=TEXT,
)

ax.tick_params(axis="x", which="major", labelsize=7.5, colors=TEXT, width=0.65, length=3)
ax.tick_params(axis="x", which="minor", bottom=False)
ax.tick_params(axis="y", left=False, length=0, colors=TEXT)
for label in ax.get_xticklabels():
    label.set_fontproperties(FONT)

# Classic manuscript styling: no grid; only bottom axis retained.
ax.grid(False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_color(TEXT)

fig.subplots_adjust(left=0.20, right=0.96, top=0.78, bottom=0.25)

fig.savefig(OUT / "Cox_OS_forest_plot.pdf", format="pdf")
fig.savefig(OUT / "Cox_OS_forest_plot.svg", format="svg")
fig.savefig(OUT / "Cox_OS_forest_plot_600dpi.png", format="png", dpi=600, facecolor="white")
plt.close(fig)

with (OUT / "Cox_OS_results_table.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle)
    writer.writerow(["Group", "Patients", "Events", "HR (95% CI)", "Cox P"])
    writer.writerow(["Low", 62, 39, "1.00 (Reference)", "—"])
    writer.writerow(["High", 62, 53, "1.72 (1.14–2.61)", "0.010"])
