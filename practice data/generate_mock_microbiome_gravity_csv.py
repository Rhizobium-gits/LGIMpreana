#!/usr/bin/env python3
"""generate_mock_microbiome_gravity_csv.py

A reproducible *synthetic* 16S rRNA (genus-level) dataset generator for practice.

It creates:
  1) a counts table (integer read counts)
  2) a relative abundance table (floats summing to 1 per sample)

Design goals
------------
- 3 donors (D1, D2, D3) with distinct, realistic-ish baseline gut microbiome profiles
- 5 gravity conditions: 0g, 1_6g, 1g, 1g_s (shake), 5g
- 3 culture times: 8h, 16h, 24h
- 3 biological/technical replicates per condition-time-donor
- Total: 3 baseline + (5*3*3*3)=135 cultured samples = 138 samples

Notes
-----
- This dataset is **synthetic**. It is built to look like a plausible human gut 16S genus table.
- Network inference on compositional data should use log-ratio transforms (e.g., CLR). This script only
  generates the raw table; analysis choices are up to you.
- A fixed RNG seed makes outputs reproducible.

Usage
-----
python generate_mock_microbiome_gravity_csv.py \
  --seed 20251223 \
  --out-counts gut_microbiome_16S_mock_gravity_culture.csv \
  --out-relative gut_microbiome_16S_mock_gravity_culture_relative_abundance.csv

Dependencies
------------
- numpy
- pandas
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


# ----------------------------
# Configuration / Taxa list
# ----------------------------
TAXA: List[str] = [
  "Bacteroides",
  "Prevotella",
  "Faecalibacterium",
  "Roseburia",
  "Blautia",
  "Agathobacter",
  "Subdoligranulum",
  "Ruminococcus",
  "Ruminococcus_gnavus_group",
  "Ruminococcaceae_UCG-002",
  "Ruminococcaceae_UCG-005",
  "Ruminococcaceae_NK4A214_group",
  "Lachnospira",
  "Lachnospiraceae_NK4A136_group",
  "Lachnospiraceae_UCG-001",
  "Lachnospiraceae_UCG-004",
  "Lachnoclostridium",
  "Eubacterium_hallii_group",
  "Eubacterium_eligens_group",
  "Eubacterium_coprostanoligenes_group",
  "Coprococcus",
  "Anaerostipes",
  "Dorea",
  "Fusicatenibacter",
  "Butyricicoccus",
  "Christensenellaceae_R-7_group",
  "Oscillibacter",
  "Dialister",
  "Veillonella",
  "Megamonas",
  "Phascolarctobacterium",
  "Parabacteroides",
  "Alistipes",
  "Barnesiella",
  "Odoribacter",
  "Butyricimonas",
  "Bilophila",
  "Desulfovibrio",
  "Sutterella",
  "Parasutterella",
  "Akkermansia",
  "Bifidobacterium",
  "Collinsella",
  "Eggerthella",
  "Enterococcus",
  "Streptococcus",
  "Lactobacillus",
  "Enterobacter",
  "Escherichia_Shigella",
  "Klebsiella",
  "Citrobacter",
  "Clostridium_sensu_stricto_1",
  "Clostridioides",
  "Peptostreptococcus",
  "Finegoldia",
  "Paraprevotella",
  "Romboutsia",
  "Turicibacter",
  "Holdemanella",
  "Intestinibacter",
  "Catenibacterium",
  "Anaerotruncus",
  "Parvimonas",
  "Fusobacterium",
  "Haemophilus",
  "Campylobacter",
  "Helicobacter",
  "Methanobrevibacter",
  "Anaerococcus",
  "Porphyromonas",
  "Alloprevotella",
  "Granulicatella",
  "Actinomyces",
  "Bacillus",
  "Staphylococcus",
  "Corynebacterium",
  "Pseudomonas",
  "Acinetobacter",
  "Enterocloster",
  "Flavonifractor",
  "Butyrivibrio",
  "Clostridiales_vadinBB60_group",
  "Family_XIII_AD3011_group",
  "Monoglobus",
  "Bacteroidetes_unclassified",
  "Firmicutes_unclassified",
  "Lachnospiraceae_unclassified",
  "Ruminococcaceae_unclassified",
  "Bacteroidales_unclassified",
  "Proteobacteria_unclassified",
  "Actinobacteria_unclassified",
  "Verrucomicrobia_unclassified",
  "Erysipelotrichaceae_UCG-003",
  "Erysipelotrichaceae_UCG-006",
  "Erysipelatoclostridium",
  "Senegalimassilia",
  "Slackia",
  "Gordonibacter",
  "Anaerobutyricum",
  "Coprobacter",
  "Murimonas",
  "Anaerofustis",
  "Pseudoflavonifractor",
  "Oscillospira",
  "Anaeroplasma",
  "UCG-010",
  "UCG-014",
  "UCG-003",
  "UCG-009",
  "UCG-013",
  "UCG-007",
  "Rikenellaceae_RC9_gut_group",
  "Candidatus_Saccharimonas",
  "Mitsuokella",
  "Succinivibrio",
  "Succiniclasticum",
  "Victivallis",
  "Anaerovorax",
  "Howardella",
  "Butyricimonas_group"
]

META_COLS: List[str] = [
    "SampleID",
    "SampleType",
    "Donor",
    "Gravity",
    "Time",
    "Replicate",
    "TotalReads",
]

GRAVITIES: List[str] = ["0g", "1_6g", "1g", "1g_s", "5g"]
TIMES: List[str] = ["8h", "16h", "24h"]
DONORS: List[int] = [1, 2, 3]
REPLICATES: List[int] = [1, 2, 3]


# ----------------------------
# Helper utilities
# ----------------------------
def _normalize(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    s = float(x.sum())
    if s <= 0:
        raise ValueError("Vector sum must be > 0.")
    return x / s


def _safe_dirichlet(rng: np.random.Generator, alpha: np.ndarray) -> np.ndarray:
    alpha = np.asarray(alpha, dtype=float)
    if np.any(alpha <= 0):
        raise ValueError("Dirichlet alpha must be strictly positive.")
    return rng.dirichlet(alpha)


def _draw_total_reads(rng: np.random.Generator) -> int:
    """Draw a realistic-ish library size around ~24k reads."""
    # lognormal gives a slight right tail; then clamp
    x = int(rng.lognormal(mean=np.log(24000), sigma=0.12))
    return int(np.clip(x, 15000, 45000))


def _make_sample_id(gravity: str, time: str, donor: int, rep: int) -> str:
    return "G" + gravity + "_T" + time + "_D" + str(donor) + "_R" + str(rep)


# ----------------------------
# Domain-ish priors (gut-like)
# ----------------------------

# Taxa groups used to impose "plausible" culture shifts.
FACULTATIVE = {
    "Escherichia_Shigella",
    "Enterobacter",
    "Klebsiella",
    "Enterococcus",
    "Streptococcus",
    "Lactobacillus",
    "Veillonella",
    "Haemophilus",
}

BUTYRATE_PRODUCERS = {
    "Faecalibacterium",
    "Roseburia",
    "Agathobacter",
    "Subdoligranulum",
    "Anaerostipes",
    "Coprococcus",
    "Eubacterium_hallii_group",
    "Eubacterium_eligens_group",
    "Eubacterium_coprostanoligenes_group",
}

BILE_TOLERANT = {
    "Bacteroides",
    "Alistipes",
    "Parabacteroides",
    "Bilophila",
}

MUCIN_ASSOCIATED = {"Akkermansia"}


def _base_target_composition() -> Dict[str, float]:
    """A reasonable-ish 'average' gut genus composition (not a real reference).

    Values are rough targets and are later adjusted per donor and then sampled
    via a Dirichlet distribution.
    """
    # Start with a tiny baseline for everything
    w: Dict[str, float] = {t: 1e-4 for t in TAXA}

    # Dominant / common gut genera & groups
    w.update(
        {
            "Bacteroides": 0.20,
            "Prevotella": 0.08,
            "Faecalibacterium": 0.08,
            "Blautia": 0.06,
            "Roseburia": 0.05,
            "Agathobacter": 0.04,
            "Subdoligranulum": 0.03,
            "Ruminococcus": 0.03,
            "Ruminococcus_gnavus_group": 0.02,
            "Alistipes": 0.03,
            "Parabacteroides": 0.02,
            "Akkermansia": 0.02,
            "Bifidobacterium": 0.02,
            "Lachnospira": 0.015,
            "Lachnoclostridium": 0.01,
            "Anaerostipes": 0.012,
            "Coprococcus": 0.012,
            "Dialister": 0.008,
            "Collinsella": 0.01,
            "Dorea": 0.01,
            "Oscillibacter": 0.008,
            "Bilophila": 0.003,
            "Escherichia_Shigella": 0.002,
            "Streptococcus": 0.002,
            "Enterococcus": 0.0015,
            "Veillonella": 0.002,
            "Lactobacillus": 0.0015,
        }
    )

    # Spread some mass across many low-abundance taxa to look realistic
    for t in TAXA:
        if w[t] <= 1e-4:
            w[t] = 0.0002

    # Normalize to sum=1
    s = float(sum(w.values()))
    return {k: v / s for k, v in w.items()}


def _apply_donor_archetype(target: Dict[str, float], donor: int) -> Dict[str, float]:
    """Make D1/Bacteroides-type, D2/Prevotella-type, D3/Ruminococcaceae-type baselines."""
    w = dict(target)

    def mult(taxon: str, factor: float) -> None:
        if taxon in w:
            w[taxon] *= factor

    if donor == 1:
        # Bacteroides-dominant, moderate Faecalibacterium, some Bifidobacterium
        mult("Bacteroides", 2.2)
        mult("Prevotella", 0.15)
        mult("Faecalibacterium", 1.3)
        mult("Bifidobacterium", 1.7)
        mult("Akkermansia", 0.9)
        mult("Alistipes", 1.2)
        mult("Parabacteroides", 1.2)
    elif donor == 2:
        # Prevotella-dominant, often co-occurs with Dialister/Veillonella
        mult("Prevotella", 2.8)
        mult("Bacteroides", 0.7)
        mult("Dialister", 1.7)
        mult("Veillonella", 1.7)
        mult("Faecalibacterium", 1.1)
        mult("Akkermansia", 0.6)
        mult("Bifidobacterium", 0.6)
    elif donor == 3:
        # More butyrate producers & mucin-associated (Akkermansia)
        mult("Faecalibacterium", 2.0)
        mult("Akkermansia", 2.3)
        mult("Subdoligranulum", 1.6)
        mult("Roseburia", 1.4)
        mult("Ruminococcus", 1.4)
        mult("Prevotella", 0.35)
        mult("Bacteroides", 0.9)
    else:
        raise ValueError("Donor must be 1, 2, or 3.")

    s = float(sum(w.values()))
    return {k: v / s for k, v in w.items()}


def _dirichlet_from_target(
    rng: np.random.Generator, target: Dict[str, float], concentration: float
) -> np.ndarray:
    """Sample a composition around the target with given concentration."""
    p = np.array([target[t] for t in TAXA], dtype=float)
    p = _normalize(p)
    alpha = np.maximum(p * float(concentration), 1e-6)
    return _safe_dirichlet(rng, alpha)


def _time_effect_logfc(time: str) -> Dict[str, float]:
    """Log-fold-change vector for culture time (synthetic)."""
    if time == "8h":
        k = 0.6
    elif time == "16h":
        k = 1.0
    elif time == "24h":
        k = 1.4
    else:
        raise ValueError("Time must be 8h/16h/24h.")

    eff: Dict[str, float] = {}

    # Facultatives bloom under batch culture / mild oxygen etc.
    for t in FACULTATIVE:
        eff[t] = +0.35 * k

    # Butyrate producers can be depleted in some culture settings.
    for t in BUTYRATE_PRODUCERS:
        eff[t] = eff.get(t, 0.0) + (-0.22 * k)

    # Bile tolerant taxa are relatively stable; mild trend only.
    for t in BILE_TOLERANT:
        eff[t] = eff.get(t, 0.0) + (-0.05 * k)

    # Mucin-associated: mild drift
    for t in MUCIN_ASSOCIATED:
        eff[t] = eff.get(t, 0.0) + (-0.03 * k)

    return eff


def _gravity_effect_logfc(gravity: str) -> Dict[str, float]:
    """Small gravity-specific shifts (synthetic)."""
    eff: Dict[str, float] = {}

    if gravity == "1g":
        return eff

    if gravity == "1_6g":
        fac = 0.04
        ana = -0.02
    elif gravity == "0g":
        fac = 0.10
        ana = -0.04
    elif gravity == "1g_s":
        # shaking: more mixing / oxygen ingress
        fac = 0.18
        ana = -0.10
    elif gravity == "5g":
        fac = 0.08
        ana = -0.05
    else:
        raise ValueError("Unknown gravity.")

    for t in FACULTATIVE:
        eff[t] = eff.get(t, 0.0) + fac
    for t in BUTYRATE_PRODUCERS:
        eff[t] = eff.get(t, 0.0) + ana

    # Small tilt under extremes
    if gravity in ("0g", "5g"):
        eff["Bacteroides"] = eff.get("Bacteroides", 0.0) + 0.04
        eff["Prevotella"] = eff.get("Prevotella", 0.0) - 0.04

    return eff


def _apply_logfc(p0: np.ndarray, logfc: Dict[str, float]) -> np.ndarray:
    """Apply log-fold-change to a composition vector and renormalize."""
    p = p0.copy().astype(float)
    idx = {t: i for i, t in enumerate(TAXA)}
    for t, v in logfc.items():
        if t in idx:
            p[idx[t]] *= float(np.exp(v))
    return _normalize(p)


def _sample_culture_composition(
    rng: np.random.Generator, baseline_p: np.ndarray, gravity: str, time: str
) -> np.ndarray:
    """Generate a culture composition from donor baseline by applying effects + Dirichlet noise."""
    logfc: Dict[str, float] = {}
    logfc.update(_time_effect_logfc(time))
    g = _gravity_effect_logfc(gravity)
    for k, v in g.items():
        logfc[k] = logfc.get(k, 0.0) + v

    p_mean = _apply_logfc(baseline_p, logfc)

    # Replicate variability: lower concentration -> more variability.
    if time == "8h":
        conc = 220.0
    elif time == "16h":
        conc = 190.0
    else:
        conc = 160.0

    alpha = np.maximum(p_mean * conc, 1e-6)
    return _safe_dirichlet(rng, alpha)


@dataclass(frozen=True)
class Dataset:
    counts: pd.DataFrame
    rel_abundance: pd.DataFrame


def generate(seed: int = 20251223) -> Dataset:
    rng = np.random.default_rng(seed)

    # 1) Donor baseline compositions
    base_target = _base_target_composition()

    donor_baseline_p: Dict[int, np.ndarray] = {}
    baseline_rows: List[dict] = []

    for donor in DONORS:
        target = _apply_donor_archetype(base_target, donor)
        p0 = _dirichlet_from_target(rng, target, concentration=420.0)
        donor_baseline_p[donor] = p0

        total = _draw_total_reads(rng)
        counts = rng.multinomial(total, p0)

        row = {
            "SampleID": "D" + str(donor),
            "SampleType": "baseline",
            "Donor": donor,
            "Gravity": "baseline",
            "Time": "0h",
            "Replicate": 0,
            "TotalReads": int(total),
        }
        row.update({t: int(c) for t, c in zip(TAXA, counts)})
        baseline_rows.append(row)

    # 2) Culture samples
    culture_rows: List[dict] = []
    for gravity in GRAVITIES:
        for time in TIMES:
            for donor in DONORS:
                for rep in REPLICATES:
                    sid = _make_sample_id(gravity, time, donor, rep)

                    p = _sample_culture_composition(rng, donor_baseline_p[donor], gravity, time)
                    total = _draw_total_reads(rng)
                    counts = rng.multinomial(total, p)

                    row = {
                        "SampleID": sid,
                        "SampleType": "culture",
                        "Donor": donor,
                        "Gravity": gravity,
                        "Time": time,
                        "Replicate": rep,
                        "TotalReads": int(total),
                    }
                    row.update({t: int(c) for t, c in zip(TAXA, counts)})
                    culture_rows.append(row)

    # 3) Assemble counts table
    df = pd.DataFrame(baseline_rows + culture_rows)
    df = df[META_COLS + TAXA]  # column order

    # 4) Relative abundance table
    rel = df.copy()
    rel[TAXA] = rel[TAXA].div(rel["TotalReads"].astype(float), axis=0)

    return Dataset(counts=df, rel_abundance=rel)


def write_outputs(ds: Dataset, out_counts: Path, out_relative: Path) -> None:
    out_counts.parent.mkdir(parents=True, exist_ok=True)
    out_relative.parent.mkdir(parents=True, exist_ok=True)
    ds.counts.to_csv(out_counts, index=False)
    ds.rel_abundance.to_csv(out_relative, index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=20251223, help="RNG seed (default: 20251223)")
    ap.add_argument("--out-counts", type=Path, required=True, help="Output CSV path for counts table")
    ap.add_argument("--out-relative", type=Path, required=True, help="Output CSV path for relative abundance table")
    args = ap.parse_args()

    ds = generate(seed=args.seed)
    write_outputs(ds, args.out_counts, args.out_relative)

    print("Wrote counts   :", args.out_counts)
    print("Wrote relative :", args.out_relative)
    print("Samples:", len(ds.counts), "(expected 138)")
    print("Taxa   :", len(TAXA), "(expected 120)")


if __name__ == "__main__":
    main()
