#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
SUPP_FIG_S5_V07_WITHIN_WORKFLOW_ISOFORM_HEATMAPS.py

Scientific redesign of Supplementary Figure S5.

WHY V07
-------
The previous direct four-workflow shared-row heatmaps were not informative:
AP contained many workflow-specific exact isoforms, while only very few
exact isoforms had interpretable cell-line profiles in multiple workflows.

V07 therefore separates two questions:

S5A
    How many exact isoforms are quantitatively accessible in each workflow?

    Two quantities are shown:
    - "Any quantitative evidence": >=1 observed run
    - "Cell-line profile": >=1 observed replicate in >=2 cell lines

S5B
    What cell-line abundance patterns are recovered WITHIN each workflow?

    Each workflow receives its OWN exact-isoform row set.
    Rows are NOT forced to match across AP/FP/MM/MQ.
    This avoids turning workflow-specific missingness into a giant grey block.

For each workflow:
    exact single_isoform_unique peptides
    -> positive intensity
    -> log2
    -> within-run median centering
    -> median exact-isoform abundance
    -> profile support: >=1 replicate in >=2 cell lines
    -> Perseus-style imputation within supported rows
       downshift 1.8 SD, width 0.3 SD
    -> row-wise z-score
    -> rank by between-cell-line variance
    -> show up to TOP_ROWS_PER_WORKFLOW
    -> row hierarchical clustering: Euclidean + average linkage

Columns are kept fixed:
    C33A_1 C33A_2 C33A_3
    SiHa_1 SiHa_2 SiHa_3
    HeLa_1 HeLa_2 HeLa_3

IMPORTANT
---------
V07 imports the already-working V05 extraction engine from the SAME Code folder:

    SUPP_FIG_S5_ISOFORM_HEATMAP_4WORKFLOW_V05_RELAXED_UNION.py

This prevents a new source-discovery bug and reuses the raw-source extraction
that already produced the real AP / FP / MM / MQ MBR-OFF matrices.

OUTPUTS
-------
FigureS5A_V07_WORKFLOW_ISOFORM_ACCESSIBILITY.png
FigureS5A_V07_WORKFLOW_ISOFORM_ACCESSIBILITY_600dpi.tiff

FigureS5B_V07_WITHIN_WORKFLOW_ISOFORM_HEATMAPS.png
FigureS5B_V07_WITHIN_WORKFLOW_ISOFORM_HEATMAPS_600dpi.tiff

plus QC / selection / matrix audit CSV files.
"""

from __future__ import annotations

import importlib.util
import json
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist


# =============================================================================
# 01. IMPORT THE VALIDATED V05 ENGINE
# =============================================================================

SCRIPT_DIR = Path(__file__).resolve().parent

V05_FILE = (
    SCRIPT_DIR
    / "SUPP_FIG_S5_ISOFORM_HEATMAP_4WORKFLOW_V05_RELAXED_UNION.py"
)

if not V05_FILE.exists():
    raise FileNotFoundError(
        "\nV07 requires the validated V05 script in the SAME Code folder:\n"
        f"{V05_FILE}\n"
    )

spec = importlib.util.spec_from_file_location(
    "s5_v05_engine",
    V05_FILE,
)

if spec is None or spec.loader is None:
    raise RuntimeError(
        f"Could not import V05 engine:\n{V05_FILE}"
    )

v5 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v5)


# =============================================================================
# 02. GLOBAL SETTINGS
# =============================================================================

ROOT = v5.ROOT

OUT = (
    ROOT
    / "MANUSCRIPT_REVISION_20260813"
    / "MAIN_FIGURES_REBUILD_20260815_V01"
    / "SUPPLEMENTARY_FIGURES"
    / "SUPP_FIG_S5_HEATMAP"
    / "V07_WITHIN_WORKFLOW"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

# Redirect all V05 audit output into the V07 folder.
v5.OUT = OUT

PROGRAMS = list(v5.PROGRAMS)
CELL_LINES = list(v5.CELL_LINES)
REPLICATES = list(v5.REPLICATES)
SAMPLES = list(v5.SAMPLES)

RANDOM_SEED = int(v5.RANDOM_SEED)
ZSCORE_LIMIT = float(v5.ZSCORE_LIMIT)

IMPUTE_DOWNSHIFT_SD = float(v5.IMPUTE_DOWNSHIFT_SD)
IMPUTE_WIDTH_SD = float(v5.IMPUTE_WIDTH_SD)

# -------------------------------------------------------------------------
# Biological support criterion
# -------------------------------------------------------------------------

# Very relaxed, but still meaningful for a C33A / SiHa / HeLa comparison:
# one observed replicate in at least two cell lines.
MIN_OBSERVED_REPS_PER_CELL_LINE = 1
MIN_OBSERVED_CELL_LINES = 2

# Show all if <= 24; otherwise select the most cell-line-variable isoforms.
TOP_ROWS_PER_WORKFLOW = 24

# -------------------------------------------------------------------------
# Typography — restrained, around approved supplementary scale
# -------------------------------------------------------------------------

PROGRAM_TITLE_PT = 10.2
CELL_LINE_PT = 7.3
SAMPLE_PT = 6.4
ROW_PT_MAX = 6.0
ROW_PT_MIN = 4.6
LEGEND_PT = 7.0
BAR_LABEL_PT = 7.2
BAR_TICK_PT = 7.5

# -------------------------------------------------------------------------
# Output
# -------------------------------------------------------------------------

S5A_PNG = OUT / "FigureS5A_V07_WORKFLOW_ISOFORM_ACCESSIBILITY.png"
S5A_TIFF = OUT / "FigureS5A_V07_WORKFLOW_ISOFORM_ACCESSIBILITY_600dpi.tiff"

S5B_PNG = OUT / "FigureS5B_V07_WITHIN_WORKFLOW_ISOFORM_HEATMAPS.png"
S5B_TIFF = OUT / "FigureS5B_V07_WITHIN_WORKFLOW_ISOFORM_HEATMAPS_600dpi.tiff"


# =============================================================================
# 03. SMALL HELPERS
# =============================================================================

def save_tiff_lzw(
    fig,
    path: Path,
    dpi: int = 600,
) -> None:
    try:
        fig.savefig(
            path,
            dpi=dpi,
            format="tiff",
            facecolor="white",
            bbox_inches="tight",
            pad_inches=0.04,
            pil_kwargs={
                "compression": "tiff_lzw",
            },
        )
    except Exception as exc:
        warnings.warn(
            f"LZW TIFF save failed ({exc}); saving uncompressed TIFF."
        )

        fig.savefig(
            path,
            dpi=dpi,
            format="tiff",
            facecolor="white",
            bbox_inches="tight",
            pad_inches=0.04,
        )


def dynamic_row_font(
    n_rows: int,
) -> float:
    if n_rows <= 12:
        return ROW_PT_MAX

    if n_rows <= 18:
        return 5.6

    if n_rows <= 24:
        return 5.1

    return ROW_PT_MIN


def observed_count(
    matrix: pd.DataFrame,
    accession: str,
) -> int:

    if accession not in matrix.index:
        return 0

    vals = matrix.loc[
        accession,
        SAMPLES,
    ].to_numpy(
        dtype=float
    )

    return int(
        np.isfinite(
            vals
        ).sum()
    )


def cell_line_observation_count(
    matrix: pd.DataFrame,
    accession: str,
) -> int:

    if accession not in matrix.index:
        return 0

    row = matrix.loc[
        accession
    ]

    n_cell_lines = 0

    for cell in CELL_LINES:

        cols = [
            f"{cell}_{rep}"
            for rep in REPLICATES
        ]

        vals = row[
            cols
        ].to_numpy(
            dtype=float
        )

        if (
            np.isfinite(
                vals
            ).sum()
            >= MIN_OBSERVED_REPS_PER_CELL_LINE
        ):
            n_cell_lines += 1

    return int(
        n_cell_lines
    )


# =============================================================================
# 04. BUILD REAL EXACT-ISOFORM MATRICES USING THE VALIDATED V05 ENGINE
# =============================================================================

def build_all_workflow_matrices():

    (
        frozen,
        exact_map,
        mapping_source,
    ) = v5.load_frozen_mapping()

    gene_lookup = v5.build_gene_lookup(
        exact_map
    )

    iso_matrix_by_program = {}
    iso_long_by_program = {}
    qc_rows = []

    for program in PROGRAMS:

        q = v5.extract_program(
            program
        )

        centered = v5.center_peptide_quant(
            program,
            q,
            frozen,
        )

        (
            iso_long,
            iso_matrix,
        ) = v5.aggregate_exact_isoforms(
            program,
            centered,
            exact_map,
        )

        iso_matrix_by_program[
            program
        ] = iso_matrix

        iso_long_by_program[
            program
        ] = iso_long

        exact_keys = set(
            exact_map.loc[
                exact_map[
                    "Program"
                ].eq(
                    program
                ),
                "PeptideKey",
            ]
        )

        centered[
            centered[
                "PeptideKey"
            ].isin(
                exact_keys
            )
        ].to_csv(
            OUT
            / f"03_{program}_single_unique_peptide_quant_centered.csv",
            index=False,
        )

        iso_long.to_csv(
            OUT
            / f"04_{program}_exact_isoform_abundance_long.csv",
            index=False,
        )

        iso_matrix.to_csv(
            OUT
            / f"05_{program}_exact_isoform_matrix_preimputation.csv"
        )

        any_quant = int(
            iso_matrix.notna().any(
                axis=1
            ).sum()
        )

        profile_supported = 0

        for accession in iso_matrix.index:

            if (
                cell_line_observation_count(
                    iso_matrix,
                    accession,
                )
                >= MIN_OBSERVED_CELL_LINES
            ):
                profile_supported += 1

        qc_rows.append(
            {
                "Program": program,
                "ExactSingleUniquePeptides": int(
                    exact_map.loc[
                        exact_map[
                            "Program"
                        ].eq(
                            program
                        ),
                        "PeptideKey",
                    ].nunique()
                ),
                "ExactIsoformsWithAnyQuant": any_quant,
                "ExactIsoformsWithCellLineProfile": profile_supported,
                "ProfileCriterion": (
                    f">={MIN_OBSERVED_REPS_PER_CELL_LINE} observed replicate "
                    f"in >={MIN_OBSERVED_CELL_LINES} cell lines"
                ),
            }
        )

    qc = pd.DataFrame(
        qc_rows
    )

    qc.to_csv(
        OUT
        / "06_V07_workflow_isoform_accessibility_summary.csv",
        index=False,
    )

    return (
        frozen,
        exact_map,
        mapping_source,
        gene_lookup,
        iso_matrix_by_program,
        iso_long_by_program,
        qc,
    )


# =============================================================================
# 05. S5A — CLEAN ACCESSIBILITY BAR PLOT
# =============================================================================

def draw_s5a_accessibility(
    qc: pd.DataFrame,
) -> None:

    x = np.arange(
        len(
            PROGRAMS
        )
    )

    any_vals = (
        qc.set_index(
            "Program"
        )
        .loc[
            PROGRAMS,
            "ExactIsoformsWithAnyQuant",
        ]
        .to_numpy(
            dtype=float
        )
    )

    profile_vals = (
        qc.set_index(
            "Program"
        )
        .loc[
            PROGRAMS,
            "ExactIsoformsWithCellLineProfile",
        ]
        .to_numpy(
            dtype=float
        )
    )

    fig, ax = plt.subplots(
        figsize=(
            6.0,
            3.7,
        ),
        facecolor="white",
    )

    width = 0.31

    bars1 = ax.bar(
        x
        - width
        / 2,
        any_vals,
        width,
        label="Any quantitative evidence",
        color="#B8C9D8",
        edgecolor="none",
    )

    bars2 = ax.bar(
        x
        + width
        / 2,
        profile_vals,
        width,
        label="Cell-line profile",
        color="#4C78A8",
        edgecolor="none",
    )

    ax.set_xticks(
        x
    )

    ax.set_xticklabels(
        PROGRAMS,
        fontsize=9.2,
        fontweight="bold",
    )

    ax.set_ylabel(
        "Exact isoforms",
        fontsize=8.8,
    )

    ax.tick_params(
        axis="y",
        labelsize=7.5,
        length=3,
    )

    ax.tick_params(
        axis="x",
        length=0,
    )

    ax.spines[
        "top"
    ].set_visible(
        False
    )

    ax.spines[
        "right"
    ].set_visible(
        False
    )

    ax.spines[
        "left"
    ].set_linewidth(
        0.6
    )

    ax.spines[
        "bottom"
    ].set_linewidth(
        0.6
    )

    ax.grid(
        axis="y",
        color="#E5E8EB",
        linewidth=0.55,
        zorder=0,
    )

    ax.set_axisbelow(
        True
    )

    ymax = max(
        1,
        int(
            max(
                any_vals.max(),
                profile_vals.max(),
            )
        ),
    )

    ax.set_ylim(
        0,
        ymax
        * 1.17,
    )

    for bars in [
        bars1,
        bars2,
    ]:

        for bar in bars:

            h = bar.get_height()

            ax.text(
                bar.get_x()
                + bar.get_width()
                / 2,
                h
                + ymax
                * 0.018,
                f"{int(h)}",
                ha="center",
                va="bottom",
                fontsize=BAR_LABEL_PT,
                fontweight="bold",
            )

    ax.legend(
        frameon=False,
        fontsize=7.1,
        loc="upper right",
    )

    fig.tight_layout(
        pad=0.7
    )

    fig.savefig(
        S5A_PNG,
        dpi=600,
        facecolor="white",
        bbox_inches="tight",
        pad_inches=0.04,
    )

    save_tiff_lzw(
        fig,
        S5A_TIFF,
        dpi=600,
    )

    plt.close(
        fig
    )


# =============================================================================
# 06. WORKFLOW-SPECIFIC PROFILE SUPPORT + PERSEUS PROCESSING
# =============================================================================

def make_profile_support_series(
    matrix: pd.DataFrame,
) -> pd.Series:

    support = pd.Series(
        False,
        index=matrix.index,
        dtype=bool,
    )

    for accession in matrix.index:

        support.loc[
            accession
        ] = (
            cell_line_observation_count(
                matrix,
                accession,
            )
            >= MIN_OBSERVED_CELL_LINES
        )

    return support


def cell_line_median_variance(
    z_matrix: pd.DataFrame,
    accession: str,
) -> float:

    if accession not in z_matrix.index:
        return float(
            "-inf"
        )

    medians = []

    for cell in CELL_LINES:

        cols = [
            f"{cell}_{rep}"
            for rep in REPLICATES
        ]

        vals = z_matrix.loc[
            accession,
            cols,
        ].to_numpy(
            dtype=float
        )

        vals = vals[
            np.isfinite(
                vals
            )
        ]

        if len(
            vals
        ) == 0:

            medians.append(
                np.nan
            )

        else:

            medians.append(
                float(
                    np.median(
                        vals
                    )
                )
            )

    m = np.array(
        medians,
        dtype=float,
    )

    m = m[
        np.isfinite(
            m
        )
    ]

    if len(
        m
    ) < 2:
        return float(
            "-inf"
        )

    return float(
        np.var(
            m,
            ddof=0,
        )
    )


def prepare_program_heatmap(
    program: str,
    matrix: pd.DataFrame,
    seed: int,
):

    support = make_profile_support_series(
        matrix
    )

    (
        imputed,
        z,
        flags,
    ) = v5.perseus_impute_and_zscore(
        matrix,
        support,
        seed=seed,
    )

    supported_accessions = [
        accession
        for accession in matrix.index
        if bool(
            support.get(
                accession,
                False,
            )
        )
        and accession
        in z.index
    ]

    ranking_rows = []

    for accession in supported_accessions:

        ranking_rows.append(
            {
                "Accession": accession,
                "BetweenCellLineVariance": cell_line_median_variance(
                    z,
                    accession,
                ),
                "ObservedRuns": observed_count(
                    matrix,
                    accession,
                ),
                "ObservedCellLines": cell_line_observation_count(
                    matrix,
                    accession,
                ),
            }
        )

    ranking = pd.DataFrame(
        ranking_rows
    )

    if ranking.empty:

        raise RuntimeError(
            f"{program}: no exact isoform had an interpretable "
            f">={MIN_OBSERVED_CELL_LINES}-cell-line profile."
        )

    ranking = ranking.sort_values(
        [
            "BetweenCellLineVariance",
            "ObservedRuns",
            "Accession",
        ],
        ascending=[
            False,
            False,
            True,
        ],
    ).reset_index(
        drop=True
    )

    if len(
        ranking
    ) > TOP_ROWS_PER_WORKFLOW:

        ranking = ranking.head(
            TOP_ROWS_PER_WORKFLOW
        ).copy()

    selected = ranking[
        "Accession"
    ].tolist()

    z_selected = z.loc[
        selected,
        SAMPLES,
    ].copy()

    # Any extremely rare non-finite values after the V05 imputation branch
    # are filled with 0 only for distance geometry / display robustness.
    # This should normally not be needed for supported rows.
    z_selected = z_selected.replace(
        [
            np.inf,
            -np.inf,
        ],
        np.nan,
    )

    for accession in z_selected.index:

        row = z_selected.loc[
            accession
        ]

        if row.isna().any():

            finite = row[
                row.notna()
            ]

            fill_value = (
                float(
                    finite.mean()
                )
                if len(
                    finite
                )
                else 0.0
            )

            z_selected.loc[
                accession
            ] = row.fillna(
                fill_value
            )

    if len(
        z_selected
    ) == 1:

        tree = None

        order = z_selected.index.tolist()

    else:

        tree = linkage(
            pdist(
                z_selected.to_numpy(
                    dtype=float
                ),
                metric="euclidean",
            ),
            method="average",
        )

        d = dendrogram(
            tree,
            no_plot=True,
        )

        order = [
            z_selected.index[
                i
            ]
            for i in d[
                "leaves"
            ]
        ]

    ranking[
        "SelectedForHeatmap"
    ] = ranking[
        "Accession"
    ].isin(
        order
    ).astype(
        int
    )

    ranking[
        "RowOrder"
    ] = ranking[
        "Accession"
    ].map(
        {
            accession: i
            + 1
            for i, accession in enumerate(
                order
            )
        }
    )

    ranking.to_csv(
        OUT
        / f"07_{program}_V07_heatmap_isoform_selection.csv",
        index=False,
    )

    imputed.to_csv(
        OUT
        / f"08_{program}_V07_profile_imputed.csv"
    )

    z.to_csv(
        OUT
        / f"09_{program}_V07_profile_row_zscore.csv"
    )

    flags.astype(
        int
    ).to_csv(
        OUT
        / f"10_{program}_V07_profile_imputation_mask.csv"
    )

    z_selected.loc[
        order,
        SAMPLES,
    ].to_csv(
        OUT
        / f"11_{program}_V07_final_heatmap_matrix.csv"
    )

    return {
        "program": program,
        "support": support,
        "imputed": imputed,
        "z": z,
        "flags": flags,
        "ranking": ranking,
        "order": order,
        "tree": tree,
        "matrix": z_selected.loc[
            order,
            SAMPLES,
        ].copy(),
    }


# =============================================================================
# 07. DRAW ONE PROGRAM PANEL INTO A GRID
# =============================================================================

def draw_dendrogram_axis(
    ax,
    tree,
    n_rows: int,
):

    if tree is not None:

        dendrogram(
            tree,
            orientation="left",
            no_labels=True,
            color_threshold=0,
            above_threshold_color="#3A3A3A",
            link_color_func=lambda _: "#3A3A3A",
            ax=ax,
        )

        ax.set_ylim(
            10
            * n_rows,
            0,
        )

    else:

        ax.set_ylim(
            n_rows
            - 0.5,
            -0.5,
        )

    ax.set_xticks(
        []
    )

    ax.set_yticks(
        []
    )

    for spine in ax.spines.values():
        spine.set_visible(
            False
        )


def draw_label_axis(
    ax,
    order,
    gene_lookup,
):

    n_rows = len(
        order
    )

    row_fs = dynamic_row_font(
        n_rows
    )

    ax.set_xlim(
        0,
        1,
    )

    ax.set_ylim(
        n_rows
        - 0.5,
        -0.5,
    )

    ax.axis(
        "off"
    )

    for i, accession in enumerate(
        order
    ):

        ax.text(
            0.99,
            i,
            v5.row_label(
                accession,
                gene_lookup,
            ),
            ha="right",
            va="center",
            fontsize=row_fs,
            clip_on=False,
        )


def draw_heat_axis(
    ax,
    program_result,
    cmap,
):

    matrix = program_result[
        "matrix"
    ]

    arr = matrix.to_numpy(
        dtype=float
    )

    image = ax.imshow(
        arr,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=-ZSCORE_LIMIT,
        vmax=ZSCORE_LIMIT,
        origin="upper",
    )

    ax.axvline(
        2.5,
        color="white",
        linewidth=1.05,
    )

    ax.axvline(
        5.5,
        color="white",
        linewidth=1.05,
    )

    ax.set_xticks(
        range(
            9
        )
    )

    ax.set_xticklabels(
        [
            "1",
            "2",
            "3",
            "1",
            "2",
            "3",
            "1",
            "2",
            "3",
        ],
        fontsize=SAMPLE_PT,
    )

    ax.tick_params(
        axis="x",
        length=0,
        pad=2,
    )

    ax.set_yticks(
        []
    )

    ax.set_title(
        f"{program_result['program']} (n={len(matrix)})",
        fontsize=PROGRAM_TITLE_PT,
        fontweight="bold",
        pad=23,
    )

    for x, cell in zip(
        [
            1,
            4,
            7,
        ],
        CELL_LINES,
    ):

        ax.text(
            x,
            1.012,
            cell,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=CELL_LINE_PT,
            fontweight="bold",
        )

    for side in [
        "top",
        "right",
        "bottom",
        "left",
    ]:

        ax.spines[
            side
        ].set_linewidth(
            0.45
        )

        ax.spines[
            side
        ].set_color(
            "#BDBDBD"
        )

    return image


# =============================================================================
# 08. S5B — 2x2 WORKFLOW-SPECIFIC HEATMAPS
# =============================================================================

def draw_s5b_heatmaps(
    results,
    gene_lookup,
) -> None:

    cmap = LinearSegmentedColormap.from_list(
        "Perseus_GYR_V07",
        [
            "#1A9850",
            "#FEE08B",
            "#D73027",
        ],
        N=256,
    )

    # Panel geometry is set from a fixed physical manuscript size.
    # Each workflow gets its own row set, so cells are never stretched by
    # another workflow's missing rows.
    fig = plt.figure(
        figsize=(
            11.7,
            8.4,
        ),
        facecolor="white",
    )

    outer = fig.add_gridspec(
        nrows=2,
        ncols=3,
        width_ratios=[
            1.0,
            1.0,
            0.045,
        ],
        height_ratios=[
            1.0,
            1.0,
        ],
        left=0.025,
        right=0.965,
        bottom=0.075,
        top=0.925,
        wspace=0.13,
        hspace=0.25,
    )

    positions = {
        "AP": (
            0,
            0,
        ),
        "FP": (
            0,
            1,
        ),
        "MM": (
            1,
            0,
        ),
        "MQ": (
            1,
            1,
        ),
    }

    image = None

    for program in PROGRAMS:

        r, c = positions[
            program
        ]

        inner = outer[
            r,
            c,
        ].subgridspec(
            nrows=1,
            ncols=3,
            width_ratios=[
                0.72,
                1.65,
                3.25,
            ],
            wspace=0.035,
        )

        ax_den = fig.add_subplot(
            inner[
                0,
                0,
            ]
        )

        ax_lab = fig.add_subplot(
            inner[
                0,
                1,
            ]
        )

        ax_heat = fig.add_subplot(
            inner[
                0,
                2,
            ]
        )

        result = results[
            program
        ]

        draw_dendrogram_axis(
            ax_den,
            result[
                "tree"
            ],
            len(
                result[
                    "order"
                ]
            ),
        )

        draw_label_axis(
            ax_lab,
            result[
                "order"
            ],
            gene_lookup,
        )

        image = draw_heat_axis(
            ax_heat,
            result,
            cmap,
        )

    cax = fig.add_subplot(
        outer[
            :,
            2,
        ]
    )

    if image is None:
        raise RuntimeError(
            "S5B heatmap image was not created."
        )

    cb = fig.colorbar(
        image,
        cax=cax,
    )

    cb.set_label(
        "row z-score",
        fontsize=LEGEND_PT,
        labelpad=5,
    )

    cb.set_ticks(
        [
            -2,
            -1,
            0,
            1,
            2,
        ]
    )

    cb.ax.tick_params(
        labelsize=6.4,
        length=2.0,
    )

    cb.outline.set_linewidth(
        0.45
    )

    fig.text(
        0.49,
        0.026,
        "Replicate",
        ha="center",
        va="center",
        fontsize=7.2,
    )

    fig.savefig(
        S5B_PNG,
        dpi=600,
        facecolor="white",
        bbox_inches="tight",
        pad_inches=0.04,
    )

    save_tiff_lzw(
        fig,
        S5B_TIFF,
        dpi=600,
    )

    plt.close(
        fig
    )


# =============================================================================
# 09. MAIN
# =============================================================================

def main():

    print(
        "="
        * 108
    )

    print(
        "SUPPLEMENTARY FIGURE S5 — V07 — WORKFLOW-SPECIFIC ISOFORM HEATMAP REDESIGN"
    )

    print(
        "S5A accessibility counts | S5B separate AP / FP / MM / MQ exact-isoform heatmaps"
    )

    print(
        "="
        * 108
    )

    (
        frozen,
        exact_map,
        mapping_source,
        gene_lookup,
        iso_matrix_by_program,
        iso_long_by_program,
        qc,
    ) = build_all_workflow_matrices()

    print(
        "\nWorkflow exact-isoform accessibility:"
    )

    print(
        qc.to_string(
            index=False
        )
    )

    # ---------------------------------------------------------------------
    # S5A
    # ---------------------------------------------------------------------

    draw_s5a_accessibility(
        qc
    )

    # ---------------------------------------------------------------------
    # S5B
    # ---------------------------------------------------------------------

    results = {}

    for i, program in enumerate(
        PROGRAMS
    ):

        results[
            program
        ] = prepare_program_heatmap(
            program=program,
            matrix=iso_matrix_by_program[
                program
            ],
            seed=(
                RANDOM_SEED
                + 200
                + i
            ),
        )

        print(
            f"\n{program}: "
            f"{int(qc.set_index('Program').loc[program, 'ExactIsoformsWithCellLineProfile'])} "
            f"isoforms pass the >=2-cell-line profile criterion; "
            f"{len(results[program]['order'])} displayed."
        )

    draw_s5b_heatmaps(
        results,
        gene_lookup,
    )

    # ---------------------------------------------------------------------
    # Method record
    # ---------------------------------------------------------------------

    method_record = {
        "mapping_source": str(
            mapping_source
        ),
        "exact_isoform_evidence": (
            "single_isoform_unique only; one explicit suffixed UniProt accession"
        ),
        "quantification": (
            "positive peptide intensity -> log2 -> within-run median centering -> "
            "median across exact single-isoform-unique peptides per isoform/run"
        ),
        "S5A": {
            "panel": (
                "workflow-specific exact-isoform accessibility counts"
            ),
            "any_quantitative_evidence": (
                ">=1 observed quantitative run"
            ),
            "cell_line_profile": (
                f">={MIN_OBSERVED_REPS_PER_CELL_LINE} observed replicate "
                f"in >={MIN_OBSERVED_CELL_LINES} cell lines"
            ),
        },
        "S5B": {
            "row_universe": (
                "workflow-specific; rows are deliberately NOT forced to match across software"
            ),
            "profile_support": (
                f">={MIN_OBSERVED_REPS_PER_CELL_LINE} observed replicate "
                f"in >={MIN_OBSERVED_CELL_LINES} cell lines"
            ),
            "imputation": {
                "style": (
                    "Perseus-style within supported workflow rows"
                ),
                "downshift_sd": IMPUTE_DOWNSHIFT_SD,
                "width_sd": IMPUTE_WIDTH_SD,
                "unsupported_rows_imputed": False,
            },
            "standardization": (
                "row-wise z-score within each workflow"
            ),
            "selection": (
                f"up to {TOP_ROWS_PER_WORKFLOW} isoforms per workflow, "
                "ranked by variance of C33A/SiHa/HeLa median z-scores"
            ),
            "row_clustering": (
                "independent within each workflow; Euclidean distance + average linkage"
            ),
            "column_order": (
                "fixed C33A_1-3, SiHa_1-3, HeLa_1-3"
            ),
            "color_scale": (
                "green=relatively low; yellow=center; red=relatively high"
            ),
            "interpretation": (
                "S5B tests whether each workflow recovers coherent cell-line-dependent "
                "patterns among the exact isoforms it can quantify. It is not a direct "
                "row-by-row cross-software abundance comparison."
            ),
        },
    }

    with open(
        OUT
        / "12_METHOD_RECORD_V07.json",
        "w",
        encoding="utf-8",
    ) as fh:

        json.dump(
            method_record,
            fh,
            indent=2,
            ensure_ascii=False,
        )

    print(
        "\nOUTPUTS"
    )

    print(
        f"  S5A PNG : {S5A_PNG}"
    )

    print(
        f"  S5A TIFF: {S5A_TIFF}"
    )

    print(
        f"  S5B PNG : {S5B_PNG}"
    )

    print(
        f"  S5B TIFF: {S5B_TIFF}"
    )

    print(
        f"  Folder  : {OUT}"
    )

    print(
        "="
        * 108
    )


if __name__ == "__main__":
    main()
