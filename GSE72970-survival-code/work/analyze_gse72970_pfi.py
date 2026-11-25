#!/usr/bin/env python3
"""Median-split progression-endpoint analysis for Lnc-DZANK1 in GSE72970."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "work"))
from analyze_gse72970 import characteristics, cox_binary_efron, km_curve, logrank, parse_matrix  # noqa: E402


OUT = ROOT / "outputs"
PROBE = "1561307_at"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    metadata, samples, expression = parse_matrix()
    chars = characteristics(metadata)
    titles = metadata["!Sample_title"][0]
    time = np.array([float(v) for v in chars["pfs"]], dtype=float)
    event = np.array([int(float(v)) for v in chars["pfs censored"]], dtype=int)
    median = float(np.median(expression))
    high = expression > median
    group = np.where(high, "High", "Low")

    chisq, logrank_p = logrank(time, event, high)
    _, _, hr, lo, hi, _, cox_p, _ = cox_binary_efron(time, event, high)

    with (OUT / "GSE72970_Lnc-DZANK1_patient_expression_groups_PFI.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "GEO_sample", "patient_sample", "probe_set", "Lnc-DZANK1_expression",
            "median_cutoff", "expression_group", "source_PFS_time_months",
            "source_PFS_event"
        ])
        for row in zip(samples, titles, expression, group, time, event):
            writer.writerow([
                row[0], row[1], PROBE, f"{row[2]:.9f}", f"{median:.9f}",
                row[3], f"{row[4]:.9f}", row[5]
            ])

    with (OUT / "GSE72970_Lnc-DZANK1_PFI_results.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "endpoint_label", "GEO_source_fields", "probe_set", "median_cutoff",
            "group", "n", "events", "HR_vs_low", "CI95_low", "CI95_high",
            "cox_p", "logrank_chisq", "logrank_p"
        ])
        common = ["Progression-free interval", "pfs; pfs censored", PROBE, f"{median:.9f}"]
        writer.writerow(common + [
            "Low (reference)", int(np.sum(~high)), int(np.sum(event[~high])),
            1, "", "", "", f"{chisq:.9g}", f"{logrank_p:.9g}"
        ])
        writer.writerow(common + [
            "High", int(np.sum(high)), int(np.sum(event[high])),
            f"{hr:.9g}", f"{lo:.9g}", f"{hi:.9g}", f"{cox_p:.9g}",
            f"{chisq:.9g}", f"{logrank_p:.9g}"
        ])

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })
    fig, ax = plt.subplots(figsize=(4.6, 3.6), facecolor="white")
    colors = {"Low": "#555555", "High": "#A63D40"}
    styles = {"Low": "-", "High": "-"}
    for label, mask in [("Low", ~high), ("High", high)]:
        x, y = km_curve(time[mask], event[mask])
        ax.step(
            x, y, where="post", color=colors[label], linestyle=styles[label],
            linewidth=1.5,
            label=f"{label} (n={np.sum(mask)}, events={np.sum(event[mask])})"
        )
        for censored_time in np.unique(time[mask & (event == 0)]):
            idx = np.searchsorted(x, censored_time, side="right") - 1
            ax.plot(
                censored_time, y[idx], marker="|", color=colors[label],
                markersize=6, markeredgewidth=0.9
            )

    fig.suptitle(
        "Progression-free interval by Lnc-DZANK1 expression",
        x=0.13, y=0.975, ha="left", fontsize=10, weight="bold"
    )
    fig.text(
        0.13, 0.925,
        f"GSE72970; median cutoff {median:.4f}; log-rank P = {logrank_p:.3f}",
        ha="left", va="top", fontsize=8, color="#404040"
    )
    ax.set_xlabel("Progression-free time (months)")
    ax.set_ylabel("Progression-free probability")
    ax.set_xlim(left=0)
    ax.set_ylim(0, 1.02)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(False)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(OUT / "GSE72970_Lnc-DZANK1_KM_PFI.pdf", facecolor="white")
    fig.savefig(OUT / "GSE72970_Lnc-DZANK1_KM_PFI.svg", facecolor="white")
    fig.savefig(OUT / "GSE72970_Lnc-DZANK1_KM_PFI_600dpi.png", dpi=600, facecolor="white")
    plt.close(fig)

    print(f"median={median:.9f}")
    print(f"low n={np.sum(~high)} events={np.sum(event[~high])}")
    print(f"high n={np.sum(high)} events={np.sum(event[high])}")
    print(f"log-rank p={logrank_p:.12g}")
    print(f"HR={hr:.9f} CI=({lo:.9f}, {hi:.9f}) Cox p={cox_p:.12g}")


if __name__ == "__main__":
    main()
