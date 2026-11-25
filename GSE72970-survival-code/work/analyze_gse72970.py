#!/usr/bin/env python3
"""Median-split OS analysis for lnc-DZANK1-1 in GSE72970.

The expression feature is Affymetrix probe set 1561307_at. UCSC's hg19
affyU133Plus2 alignment places its target sequence across the exons of
LNCipedia lnc-DZANK1-1; the cross-hybridizing companion 1561309_x_at is not
used. GEO's `os censored` value is treated as the event indicator (1=death),
consistent with short observed OS records being marked 1 and long-follow-up
records marked 0.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MATRIX = ROOT / "work/geo/GSE72970_series_matrix.txt.gz"
OUT = ROOT / "outputs"
PROBE = "1561307_at"


def unquote(value: str) -> str:
    return value.strip().strip('"')


def parse_matrix():
    metadata: dict[str, list[list[str]]] = {}
    expr_header = None
    expr_values = None
    in_table = False
    with gzip.open(MATRIX, "rt", newline="") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                continue
            if line == "!series_matrix_table_end":
                break
            if not in_table:
                if line.startswith("!Sample_"):
                    row = next(csv.reader([line], delimiter="\t"))
                    metadata.setdefault(row[0], []).append([unquote(v) for v in row[1:]])
                continue
            row = next(csv.reader([line], delimiter="\t"))
            if expr_header is None:
                expr_header = [unquote(v) for v in row[1:]]
            elif unquote(row[0]) == PROBE:
                expr_values = np.array([float(unquote(v)) for v in row[1:]], dtype=float)
    if expr_header is None or expr_values is None:
        raise RuntimeError(f"Probe {PROBE} not found in {MATRIX}")
    return metadata, expr_header, expr_values


def characteristics(metadata):
    result = {}
    for row in metadata["!Sample_characteristics_ch1"]:
        key, _, _ = row[0].partition(":")
        values = []
        for value in row:
            _, sep, payload = value.partition(":")
            values.append(payload.strip() if sep else value.strip())
        result[key.strip().lower()] = values
    return result


def km_curve(times, events):
    order = np.argsort(times)
    times = np.asarray(times)[order]
    events = np.asarray(events)[order]
    x = [0.0]
    y = [1.0]
    surv = 1.0
    for t in np.unique(times[events == 1]):
        n_risk = np.sum(times >= t)
        d = np.sum((times == t) & (events == 1))
        x.extend([float(t), float(t)])
        y.extend([surv, surv * (1.0 - d / n_risk)])
        surv = y[-1]
    x.append(float(times.max()))
    y.append(surv)
    return np.array(x), np.array(y)


def logrank(times, events, high):
    oe = 0.0
    variance = 0.0
    for t in np.unique(times[events == 1]):
        at_risk = times >= t
        at_t = times == t
        n1 = np.sum(at_risk & high)
        n0 = np.sum(at_risk & ~high)
        n = n1 + n0
        d1 = np.sum(at_t & (events == 1) & high)
        d0 = np.sum(at_t & (events == 1) & ~high)
        d = d1 + d0
        if n > 1:
            oe += d1 - d * n1 / n
            variance += n1 * n0 * d * (n - d) / (n * n * (n - 1))
    chisq = oe * oe / variance
    p = math.erfc(math.sqrt(chisq / 2.0))
    return chisq, p


def cox_binary_efron(times, events, high):
    """Univariate Cox PH fit with Efron handling of tied event times."""
    x = high.astype(float)

    def score_info(beta):
        ebx = np.exp(beta * x)
        score = 0.0
        info = 0.0
        loglik = 0.0
        for t in np.unique(times[events == 1]):
            deaths = (times == t) & (events == 1)
            risk = times >= t
            d = int(np.sum(deaths))
            sx0 = np.sum(ebx[risk])
            sx1 = np.sum(ebx[risk] * x[risk])
            sx2 = np.sum(ebx[risk] * x[risk] * x[risk])
            dx0 = np.sum(ebx[deaths])
            dx1 = np.sum(ebx[deaths] * x[deaths])
            dx2 = np.sum(ebx[deaths] * x[deaths] * x[deaths])
            score += np.sum(x[deaths])
            loglik += beta * np.sum(x[deaths])
            for ell in range(d):
                frac = ell / d
                r0 = sx0 - frac * dx0
                r1 = sx1 - frac * dx1
                r2 = sx2 - frac * dx2
                score -= r1 / r0
                info += r2 / r0 - (r1 / r0) ** 2
                loglik -= math.log(r0)
        return score, info, loglik

    beta = 0.0
    for _ in range(100):
        score, info, _ = score_info(beta)
        step = score / info
        beta += step
        if abs(step) < 1e-12:
            break
    _, info, loglik = score_info(beta)
    se = math.sqrt(1.0 / info)
    z = beta / se
    p = math.erfc(abs(z) / math.sqrt(2.0))
    hr = math.exp(beta)
    lo = math.exp(beta - 1.959963984540054 * se)
    hi = math.exp(beta + 1.959963984540054 * se)
    return beta, se, hr, lo, hi, z, p, loglik


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    metadata, samples, expression = parse_matrix()
    chars = characteristics(metadata)
    titles = metadata["!Sample_title"][0]
    os_time = np.array([float(v) for v in chars["os"]])
    event = np.array([int(float(v)) for v in chars["os censored"]])
    median = float(np.median(expression))
    high = expression > median
    group = np.where(high, "High", "Low")

    if not (len(samples) == len(titles) == len(expression) == len(os_time) == len(event)):
        raise RuntimeError("Sample/phenotype/expression lengths do not match")
    if np.isnan(expression).any() or np.isnan(os_time).any():
        raise RuntimeError("Unexpected missing expression or OS time")
    if not set(event).issubset({0, 1}):
        raise RuntimeError("Unexpected OS event coding")

    chisq, logrank_p = logrank(os_time, event, high)
    beta, se, hr, lo, hi, z, cox_p, loglik = cox_binary_efron(os_time, event, high)

    patient_rows = [
        {
            "GEO_sample": gsm,
            "patient_sample": title,
            "probe_set": PROBE,
            "expression": float(expr),
            "median_cutoff": median,
            "group": grp,
            "OS_months": float(months),
            "OS_event": int(evt),
        }
        for gsm, title, expr, grp, months, evt
        in zip(samples, titles, expression, group, os_time, event)
    ]
    group_rows = [
        {"group": "Low (reference)", "n": int(np.sum(~high)), "events": int(np.sum(event[~high])), "HR": 1.0, "CI_low": None, "CI_high": None, "cox_p": None},
        {"group": "High", "n": int(np.sum(high)), "events": int(np.sum(event[high])), "HR": hr, "CI_low": lo, "CI_high": hi, "cox_p": cox_p},
    ]
    km_rows = []
    for label, mask in [("Low", ~high), ("High", high)]:
        kx, ky = km_curve(os_time[mask], event[mask])
        for tx, sy in zip(kx, ky):
            km_rows.append({
                "OS_months": float(tx), "survival_probability": float(sy), "group": label,
                "group_n": int(np.sum(mask)), "group_events": int(np.sum(event[mask])),
            })

    csv_path = OUT / "GSE72970_Lnc-DZANK1_patient_expression_groups.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "GEO_sample", "patient_sample", "probe_set", "Lnc-DZANK1_expression",
            "median_cutoff", "expression_group", "OS_time_months", "OS_event"
        ])
        for row in zip(samples, titles, expression, group, os_time, event):
            writer.writerow([row[0], row[1], PROBE, f"{row[2]:.9f}", f"{median:.9f}", row[3], f"{row[4]:.9f}", row[5]])

    summary_path = OUT / "GSE72970_Lnc-DZANK1_survival_results.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["endpoint", "probe_set", "median_cutoff", "group", "n", "events", "HR_vs_low", "CI95_low", "CI95_high", "cox_p", "logrank_chisq", "logrank_p"])
        writer.writerow(["Overall survival", PROBE, f"{median:.9f}", "Low (reference)", int(np.sum(~high)), int(np.sum(event[~high])), 1, "", "", "", f"{chisq:.9g}", f"{logrank_p:.9g}"])
        writer.writerow(["Overall survival", PROBE, f"{median:.9f}", "High", int(np.sum(high)), int(np.sum(event[high])), f"{hr:.9g}", f"{lo:.9g}", f"{hi:.9g}", f"{cox_p:.9g}", f"{chisq:.9g}", f"{logrank_p:.9g}"])

    generated_at = datetime.now(timezone.utc).isoformat()
    geo_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72970"
    source = {
        "id": "gse72970",
        "label": "NCBI GEO GSE72970 series matrix",
        "path": "GSE72970_series_matrix.txt.gz",
        "href": geo_url,
    }
    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": "GSE72970 Lnc-DZANK1 median-split survival analysis",
            "description": "Kaplan–Meier, log-rank, and univariate Cox analysis using overall survival.",
            "generatedAt": generated_at,
            "sources": [source],
            "cards": [
                {"id": "hr", "dataset": "headline", "sourceId": "gse72970", "description": "Univariate Cox hazard ratio for High versus Low.", "metrics": [{"label": "HR, High vs Low", "field": "HR", "format": "number"}]},
                {"id": "cox_p", "dataset": "headline", "sourceId": "gse72970", "description": "Two-sided Wald P value from the univariate Cox model.", "metrics": [{"label": "Cox P value", "field": "cox_p", "format": "number"}]},
                {"id": "logrank_p", "dataset": "headline", "sourceId": "gse72970", "description": "Two-sided log-rank test comparing the Kaplan–Meier curves.", "metrics": [{"label": "Log-rank P", "field": "logrank_p", "format": "number"}]},
                {"id": "median", "dataset": "headline", "sourceId": "gse72970", "description": "Cohort median of normalized 1561307_at expression.", "metrics": [{"label": "Median cutoff", "field": "median", "format": "number"}]},
            ],
            "charts": [{
                "id": "km", "title": "Overall survival by Lnc-DZANK1 expression group",
                "subtitle": "GSE72970; median cutoff 4.3199; n=124; repeated x-values encode Kaplan–Meier steps",
                "intent": "trend", "question": "Do overall-survival curves differ between the median-defined expression groups?",
                "rationale": "A two-series survival curve shows the time-to-event separation and censoring context across follow-up.",
                "type": "line", "dataset": "km_curve", "sourceId": "gse72970", "layout": "full",
                "encodings": {"x": {"field": "OS_months", "type": "quantitative", "label": "Overall survival (months)"}, "y": {"field": "survival_probability", "type": "quantitative", "label": "Survival probability"}, "color": {"field": "group", "type": "nominal", "label": "Expression group"}, "lineStyle": {"field": "group", "type": "nominal"}},
                "valueFormat": "number", "unit": "probability",
            }],
            "tables": [
                {"id": "groups", "title": "Group sizes, events, and Cox estimate", "dataset": "groups", "sourceId": "gse72970", "layout": "full", "defaultSort": {"field": "group", "direction": "asc"}, "columns": [{"field": "group", "label": "Group", "type": "text"}, {"field": "n", "label": "Patients", "format": "number"}, {"field": "events", "label": "Events", "format": "number"}, {"field": "HR", "label": "HR vs Low", "format": "number"}, {"field": "CI_low", "label": "95% CI low", "format": "number"}, {"field": "CI_high", "label": "95% CI high", "format": "number"}, {"field": "cox_p", "label": "Cox P", "format": "number"}]},
                {"id": "patients", "title": "Patient-level expression, group, and OS data", "subtitle": "All 124 analyzed GSE72970 samples", "dataset": "patients", "sourceId": "gse72970", "layout": "full", "density": "dense", "defaultSort": {"field": "expression", "direction": "desc"}, "columns": [{"field": "GEO_sample", "label": "GEO sample", "type": "text"}, {"field": "patient_sample", "label": "Patient sample", "type": "text"}, {"field": "expression", "label": "Lnc-DZANK1 expression", "format": "number"}, {"field": "group", "label": "Group", "type": "text"}, {"field": "OS_months", "label": "OS (months)", "format": "number"}, {"field": "OS_event", "label": "OS event", "format": "number"}]},
            ],
            "blocks": [
                {"id": "title", "type": "markdown", "body": "# GSE72970 Lnc-DZANK1 median-split survival analysis"},
                {"id": "summary", "type": "markdown", "sourceId": "gse72970", "body": f"## Technical summary\n\nHigh Lnc-DZANK1 expression was associated with poorer overall survival. The high group had HR **{hr:.2f}** (95% CI **{lo:.2f}–{hi:.2f}**; Cox P **{cox_p:.4f}**) versus the low group. The Kaplan–Meier curves differed by log-rank P **{logrank_p:.4f}**."},
                {"id": "metrics", "type": "metric-strip", "cardIds": ["hr", "cox_p", "logrank_p", "median"]},
                {"id": "finding", "type": "markdown", "sourceId": "gse72970", "body": f"## High expression tracks with shorter overall survival\n\nThe groups were balanced at 62 patients each. There were **53 deaths** in the high-expression group and **39 deaths** in the low-expression group. This is an unadjusted association, not evidence of causality or independence from clinical covariates."},
                {"id": "km_block", "type": "chart", "chartId": "km"},
                {"id": "scope", "type": "markdown", "sourceId": "gse72970", "body": "## Cohort, feature, and endpoint definitions\n\nThe analysis includes all 124 GSE72970 samples with complete expression and OS data. Lnc-DZANK1 is represented by the non-cross-hybridizing Affymetrix probe set **1561307_at**, whose hg19 target alignment overlaps LNCipedia **lnc-DZANK1-1**. OS is measured in months; the GEO field `os censored` is used as the event indicator (1=death, 0=censored)."},
                {"id": "groups_block", "type": "table", "tableId": "groups"},
                {"id": "methods", "type": "markdown", "sourceId": "gse72970", "body": "## Median-only grouping and statistical methods\n\nThe cutoff is the cohort median (4.319903367). High is defined as expression strictly above the median and Low as expression at or below it; there were no values exactly equal to the median. Kaplan–Meier estimates use right censoring, the curve comparison uses the standard two-sample log-rank test, and the univariate Cox proportional hazards model uses Low as reference with Efron handling of tied event times."},
                {"id": "limits", "type": "markdown", "body": "## Limitations and robustness\n\nThis is a single-probe, univariate retrospective analysis. The Cox proportional-hazards assumption was not formally tested here, and clinical covariates were not adjusted. The cross-hybridizing companion probe 1561309_x_at was excluded. The label “Lnc-DZANK1” is interpreted as the LNCipedia lnc-DZANK1-1 locus rather than the protein-coding DZANK1 gene."},
                {"id": "patient_heading", "type": "markdown", "body": "## Patient-level values"},
                {"id": "patients_block", "type": "table", "tableId": "patients"},
                {"id": "next", "type": "markdown", "body": "## Recommended next steps\n\nConfirm the probe-to-lncRNA mapping against the intended Lnc-DZANK1 naming convention and test the proportional-hazards assumption before using the HR in a formal prognostic model."},
                {"id": "questions", "type": "markdown", "body": "## Further questions\n\nWould the association remain after adjustment for age, performance status, tumor location, and treatment regimen? Is it reproducible with an independent transcript-level assay?"},
            ],
        },
        "snapshot": {"version": 1, "generatedAt": generated_at, "status": "ready", "datasets": {"headline": [{"HR": hr, "cox_p": cox_p, "logrank_p": logrank_p, "median": median, "CI_low": lo, "CI_high": hi}], "groups": group_rows, "km_curve": km_rows, "patients": patient_rows}, "accessIssues": []},
        "sources": [source],
    }
    with (OUT / "GSE72970_Lnc-DZANK1_report_artifact.json").open("w") as handle:
        json.dump(artifact, handle, indent=2)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    fig, ax = plt.subplots(figsize=(8.2, 6.1), dpi=180)
    colors = {"Low": "#245B9B", "High": "#D58A20"}
    styles = {"Low": "-", "High": "--"}
    for label, mask in [("Low", ~high), ("High", high)]:
        kx, ky = km_curve(os_time[mask], event[mask])
        ax.step(kx, ky, where="post", color=colors[label], linestyle=styles[label], linewidth=2.2,
                label=f"{label} (n={np.sum(mask)}, events={np.sum(event[mask])})")
        censored = np.unique(os_time[mask & (event == 0)])
        for ct in censored:
            idx = np.searchsorted(kx, ct, side="right") - 1
            ax.plot(ct, ky[idx], marker="|", color=colors[label], markersize=8, markeredgewidth=1.2)
    fig.suptitle("GSE72970: overall survival by Lnc-DZANK1 expression",
                 x=0.13, y=0.975, ha="left", fontsize=14, weight="bold")
    fig.text(0.13, 0.935,
             f"Median split at {median:.4f} (probe {PROBE}); log-rank P = {logrank_p:.3g}",
             ha="left", va="top", color="#46505A")
    ax.set_xlabel("Overall survival (months)")
    ax.set_ylabel("Survival probability")
    ax.set_xlim(left=0)
    ax.set_ylim(0, 1.02)
    ax.grid(axis="y", color="#D9DEE3", linewidth=0.7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="lower left")
    ax.text(1, -0.16, "Ticks indicate censoring. OS_event: 1=death, 0=censored.", transform=ax.transAxes,
            ha="right", va="top", fontsize=8.5, color="#5D6670")
    fig.tight_layout(rect=(0, 0.03, 1, 0.91))
    fig.savefig(OUT / "GSE72970_Lnc-DZANK1_KM_OS.png", bbox_inches="tight")
    fig.savefig(OUT / "GSE72970_Lnc-DZANK1_KM_OS.svg", bbox_inches="tight")
    plt.close(fig)

    print(f"n={len(samples)} median={median:.9f}")
    print(f"low n={np.sum(~high)} events={np.sum(event[~high])}")
    print(f"high n={np.sum(high)} events={np.sum(event[high])}")
    print(f"logrank chisq={chisq:.9f} p={logrank_p:.12g}")
    print(f"cox beta={beta:.9f} se={se:.9f} HR={hr:.9f} CI=({lo:.9f},{hi:.9f}) z={z:.9f} p={cox_p:.12g} loglik={loglik:.9f}")


if __name__ == "__main__":
    main()
