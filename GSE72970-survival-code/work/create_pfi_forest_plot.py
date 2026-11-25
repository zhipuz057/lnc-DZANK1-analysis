#!/usr/bin/env python3
"""Create the GSE72970 progression-endpoint Cox forest plot.

GSE72970 names this endpoint `pfs`; the output uses the requested PFI label
while retaining that source-field definition in the companion results table.
"""

from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

HR = 1.2321548074864526
LOWER = 0.8505744182822137
UPPER = 1.7849178590135408
COX_P = 0.2695699372838645

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.unicode_minus": False,
})

regular = FontProperties(family="Helvetica", size=8.5)
bold = FontProperties(family="Helvetica", size=10, weight="bold")

fig, ax = plt.subplots(figsize=(3.2, 2.3), facecolor="white")
ax.set_facecolor("white")

# Null-effect line.
ax.axvline(1.0, color="#808080", linewidth=0.5, linestyle=(0, (4, 3)), zorder=1)

# High-expression estimate and 95% CI.
ax.errorbar(
    HR,
    1,
    xerr=[[HR - LOWER], [UPPER - HR]],
    fmt="o",
    color="#A63D40",
    ecolor="#A63D40",
    elinewidth=0.9,
    capsize=0,
    markersize=5.2,
    markeredgewidth=0,
    zorder=3,
)

# Low-expression reference estimate; no CI by design.
ax.plot(1.0, 0, "o", color="#555555", markersize=4.4, markeredgewidth=0, zorder=4)

ax.set_xscale("log", base=10)
ax.set_xlim(0.7, 3.2)
ax.set_xticks([0.75, 1.0, 1.5, 2.0, 3.0])
ax.set_xticklabels(["0.75", "1", "1.5", "2", "3"])
ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
ax.set_ylim(-0.65, 1.55)
ax.set_yticks([1, 0])
ax.set_yticklabels(["High", "Low"])

ax.set_xlabel("Hazard ratio (95% CI)", fontproperties=regular, labelpad=7)
ax.set_title("Progression-free interval", loc="left", fontproperties=bold, pad=8)
ax.text(
    0.99,
    1.02,
    "High: HR 1.232 (0.851–1.785)   P = 0.270",
    transform=ax.transAxes,
    ha="right",
    va="bottom",
    fontproperties=FontProperties(family="Helvetica", size=6.3),
    color="#262626",
)

for label in ax.get_xticklabels():
    label.set_fontproperties(FontProperties(family="Helvetica", size=7.5))
    label.set_color("#262626")
for label in ax.get_yticklabels():
    label.set_fontproperties(regular)
    label.set_color("#262626")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_color("#262626")
ax.spines["bottom"].set_linewidth(0.6)
ax.tick_params(axis="y", length=0)
ax.tick_params(axis="x", which="major", length=3, width=0.6, color="#262626")
ax.grid(False)

fig.subplots_adjust(left=0.20, right=0.97, bottom=0.24, top=0.79)

fig.savefig(OUT / "Cox_PFI_forest_plot.pdf", facecolor="white")
fig.savefig(OUT / "Cox_PFI_forest_plot.svg", facecolor="white")
fig.savefig(OUT / "Cox_PFI_forest_plot_600dpi.png", dpi=600, facecolor="white")
plt.close(fig)

(OUT / "Cox_PFI_results_table.csv").write_text(
    "Group,Patients,PFI events,HR (95% CI),Cox P\n"
    "Low,62,55,1.00 (Reference),—\n"
    "High,62,59,1.23 (0.85–1.78),0.270\n",
    encoding="utf-8",
)

