# ==============================================================================
# OPTICO-ACS | pySCENIC regulon inference and differential regulon activity
#
# Panels : Fig. 5i
# Input  : T_reclustering_2025/data/SCENIC/*_RNA_*counts.csv, *_metadata.csv (from Fig5g_oligoclonal_DEG.R)
#          cisTarget hg38 v10 ranking databases and motif annotation (see README)
# Output : T_reclustering_2025/data/SCENIC/SCENIC_differential_regulons.xlsx, SCENIC_inspection_tables.xlsx
# ==============================================================================

import os
OPTICO_DIR = os.environ.get("OPTICO_DIR", os.path.expanduser("~/OPTICO"))  # project folder

# %% [markdown]
# # pySCENIC — OPTICO T cells: Expanded vs Not_expanded
# **Goal:** Identify transcription factor regulons driving clonal expansion
# in CD4 and CD8 T cells; compare between Erosion and Rupture.
#
# **Workflow:**
# 1. GRN inference (GRNBoost2) — per lineage
# 2. Regulon pruning (cisTarget) — per lineage
# 3. AUCell scoring — per lineage
# 4. Differential regulon analysis: Expanded vs Not_expanded
# 5. Secondary comparison: Erosion vs Rupture within Expanded cells
# 6. Inspection table → curated lollipop plot

# %% [markdown]
# ## 0. Imports and paths

# %%
import os
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# ── Directory paths ──────────────────────────────────────────────────────────
SCENIC_DIR  = os.path.expanduser(
    os.path.join(OPTICO_DIR, "T_reclustering_2025/data/SCENIC")
)
PLOT_DIR    = os.path.join(SCENIC_DIR, "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

# ── Shared resource files (reuse from ApoB project) ──────────────────────────
# Update these paths if your resource files live elsewhere
APOB_DIR     = os.path.join(OPTICO_DIR, "external_data/SCENIC")
TF_LIST      = os.path.join(APOB_DIR, "allTFs_hg38.txt")
FEATHER_500  = os.path.join(APOB_DIR, "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
FEATHER_10K  = os.path.join(APOB_DIR, "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
MOTIF_ANNOT  = os.path.join(APOB_DIR, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl.txt")

# ── Per-lineage input files (exported from DEG_enrichGO_upset_expanded_SCENIC.R)
LINEAGES = ["CD8_all", "CD4_all"]

# Colors matching R script
COLOR_ERO = "#4DBBD5"
COLOR_RUP = "#E64B35"
COLOR_EXP = "#2C3E50"
COLOR_NOT = "#BDC3C7"

# %% [markdown]
# ## 1. GRN inference (GRNBoost2) — run once per lineage
# This is the slow step (~20–60 min per lineage depending on cell number).
# Re-use the output TSV on subsequent runs.

# %%
for lineage in LINEAGES:
    norm_csv    = os.path.join(SCENIC_DIR, f"{lineage}_RNA_norm_counts.csv")
    adj_out     = os.path.join(SCENIC_DIR, f"{lineage}_adjacencies.tsv")

    if os.path.exists(adj_out):
        print(f"[{lineage}] GRN output already exists — skipping. Delete to rerun:\n  {adj_out}")
        continue

    cmd = (
        f"pyscenic grn "
        f"--seed 1200 "
        f"--num_workers 12 "
        f"--method grnboost2 "
        f"--output {adj_out} "
        f"{norm_csv} "
        f"{TF_LIST}"
    )
    print(f"[{lineage}] Running GRN:\n  {cmd}\n")
    os.system(cmd)

# %% [markdown]
# ## 2. Regulon pruning (cisTarget) — run once per lineage
# Uses both 500bp and 10kb feather databases.

# %%
for lineage in LINEAGES:
    adj_tsv     = os.path.join(SCENIC_DIR, f"{lineage}_adjacencies.tsv")
    norm_csv    = os.path.join(SCENIC_DIR, f"{lineage}_RNA_norm_counts.csv")
    regulon_out = os.path.join(SCENIC_DIR, f"{lineage}_regulons.csv")

    if os.path.exists(regulon_out):
        print(f"[{lineage}] Regulon file already exists — skipping. Delete to rerun:\n  {regulon_out}")
        continue

    cmd = (
        f"pyscenic ctx "
        f"--nes_threshold 2.5 "
        f"--num_workers 12 "
        f"--mode dask_multiprocessing "
        f"{adj_tsv} "
        f"{FEATHER_500} {FEATHER_10K} "
        f"--expression_mtx_fname {norm_csv} "
        f"--annotations_fname {MOTIF_ANNOT} "
        f"--mask_dropouts "
        f"--output {regulon_out}"
    )
    print(f"[{lineage}] Running ctx:\n  {cmd}\n")
    os.system(cmd)

# %% [markdown]
# ## 3. AUCell scoring — run once per lineage

# %%
for lineage in LINEAGES:
    norm_csv    = os.path.join(SCENIC_DIR, f"{lineage}_RNA_norm_counts.csv")
    regulon_csv = os.path.join(SCENIC_DIR, f"{lineage}_regulons.csv")
    auc_out     = os.path.join(SCENIC_DIR, f"{lineage}_SCENIC_auc.csv")

    if os.path.exists(auc_out):
        print(f"[{lineage}] AUC file already exists — skipping. Delete to rerun:\n  {auc_out}")
        continue

    cmd = (
        f"pyscenic aucell "
        f"{norm_csv} "
        f"{regulon_csv} "
        f"--output {auc_out} "
        f"--num_workers 12 "
        f"--seed 1200 "
        f"--nes_threshold 2.5"
    )
    print(f"[{lineage}] Running AUCell:\n  {cmd}\n")
    os.system(cmd)

# %% [markdown]
# ## 4. Load AUC matrices and metadata

# %%
auc   = {}
meta  = {}

for lineage in LINEAGES:
    auc_path  = os.path.join(SCENIC_DIR, f"{lineage}_SCENIC_auc.csv")
    meta_path = os.path.join(SCENIC_DIR, f"{lineage}_metadata.csv")

    df_auc  = pd.read_csv(auc_path)
    df_auc.index = df_auc["Cell"]
    df_auc  = df_auc.drop(columns=["Cell"])

    df_meta = pd.read_csv(meta_path, index_col=0)

    # Align on shared cell barcodes
    shared = df_auc.index.intersection(df_meta.index)
    auc[lineage]  = df_auc.loc[shared]
    meta[lineage] = df_meta.loc[shared]

    print(f"\n{lineage}")
    print(f"  Cells: {len(shared)} | Regulons: {df_auc.shape[1]}")
    print(meta[lineage]["expansion_group"].value_counts())
    print(meta[lineage]["disease"].value_counts())

# %% [markdown]
# ## 5. Differential regulon analysis — Expanded vs Not_expanded
#
# For each lineage:
#   - Mann-Whitney U test per regulon (Expanded vs Not_expanded)
#   - Effect size: mean delta AUC (Expanded − Not_expanded)
#   - BH correction across all regulons within each lineage
#   - Export ranked table for inspection before plotting

# %%
np.float = np.float64   # compatibility patch for older numpy/pyscenic

def diff_regulons(auc_df, meta_df, group_col, group_a, group_b, label=""):
    """
    Mann-Whitney U test per regulon between two groups.
    Returns a DataFrame sorted by FDR then delta_AUC.
    """
    idx_a = meta_df.index[meta_df[group_col] == group_a]
    idx_b = meta_df.index[meta_df[group_col] == group_b]

    # Intersect with AUC index
    idx_a = idx_a.intersection(auc_df.index)
    idx_b = idx_b.intersection(auc_df.index)

    print(f"  {label}: {group_a} n={len(idx_a)} | {group_b} n={len(idx_b)}")

    results = []
    for reg in auc_df.columns:
        a_vals = auc_df.loc[idx_a, reg].values
        b_vals = auc_df.loc[idx_b, reg].values
        stat, pval = stats.mannwhitneyu(a_vals, b_vals, alternative="two-sided")
        results.append({
            "regulon":       reg,
            "mean_AUC_a":    a_vals.mean(),
            "mean_AUC_b":    b_vals.mean(),
            "delta_AUC":     a_vals.mean() - b_vals.mean(),   # positive = higher in group_a
            "pval":          pval,
            "group_a":       group_a,
            "group_b":       group_b,
            "n_a":           len(idx_a),
            "n_b":           len(idx_b),
        })

    df = pd.DataFrame(results)
    _, df["FDR"], _, _ = multipletests(df["pval"], method="fdr_bh")
    df = df.sort_values(["FDR", "delta_AUC"], key=lambda x: x.abs() if x.name == "delta_AUC" else x,
                        ascending=[True, False])
    return df

diff_results = {}

for lineage in LINEAGES:
    print(f"\n── {lineage} ─────────────────────────────────────────")

    # Primary: Expanded vs Not_expanded (all diseases combined)
    key = f"{lineage}__Expanded_vs_NotExpanded"
    diff_results[key] = diff_regulons(
        auc[lineage], meta[lineage],
        group_col="expansion_group",
        group_a="Expanded",
        group_b="Not_expanded",
        label=key
    )

    # Secondary: within Expanded cells — Erosion vs Rupture
    mask_exp = meta[lineage]["expansion_group"] == "Expanded"
    auc_exp  = auc[lineage].loc[mask_exp]
    meta_exp = meta[lineage].loc[mask_exp]

    key2 = f"{lineage}__Erosion_vs_Rupture_in_Expanded"
    diff_results[key2] = diff_regulons(
        auc_exp, meta_exp,
        group_col="disease",
        group_a="Erosion",
        group_b="Rupture",
        label=key2
    )

# Export all differential results to Excel
with pd.ExcelWriter(os.path.join(SCENIC_DIR, "SCENIC_differential_regulons.xlsx")) as writer:
    for key, df in diff_results.items():
        sheet = key[:31]   # Excel sheet name limit
        df.to_excel(writer, sheet_name=sheet, index=False)

print("\nDifferential regulon tables exported.")

# %% [markdown]
# ## 6. Inspection table — ranked meaningful regulons
#
# Filters for:
#   - FDR < 0.05 AND |delta_AUC| > 0.005 (practical effect size threshold)
#   - Flags known T cell–relevant TFs for priority review
#
# **Review this table and note which regulons to include in the lollipop plot.**

# %%
# Known TFs relevant to T cell biology and clonal expansion
KNOWN_TCELL_TFS = {
    # Exhaustion / dysfunction
    "TOX", "TOX2", "NR4A1", "NR4A2", "NR4A3",
    # Effector / cytotoxicity
    "TBX21", "EOMES", "RUNX3", "ID2",
    # Memory / stemness
    "TCF7", "LEF1", "FOXO1", "KLF2",
    # Helper / Th1
    "GATA3", "RORC", "RORA", "STAT4", "STAT1",
    # Regulatory T cells
    "FOXP3",
    # Co-stimulation / activation
    "BATF", "BATF3", "IRF4", "IRF7", "NFATC1", "NFATC2",
    "REL", "RELA", "NFKB1",
    # Proliferation
    "MYC", "E2F1", "E2F4",
    # Stress / inflammation
    "JUN", "FOS", "FOSL1", "FOSL2", "ATF3",
    # T cell trafficking / adhesion
    "KLF4", "ZEB1", "ZEB2",
}

def make_inspection_table(diff_df, fdr_cut=0.05, delta_cut=0.005):
    """Filter significant regulons and flag known T cell TFs."""
    sig = diff_df[
        (diff_df["FDR"] < fdr_cut) &
        (diff_df["delta_AUC"].abs() > delta_cut)
    ].copy()

    # Strip the (+) suffix to get TF name for lookup
    sig["TF_name"] = sig["regulon"].str.replace(r"\(\+\)$", "", regex=True)
    sig["known_Tcell_TF"] = sig["TF_name"].isin(KNOWN_TCELL_TFS)
    sig["direction"] = sig["delta_AUC"].apply(
        lambda x: f"↑ {sig.iloc[0]['group_a']}" if x > 0 else f"↑ {sig.iloc[0]['group_b']}"
    )
    sig["direction"] = np.where(sig["delta_AUC"] > 0, sig["group_a"], sig["group_b"])
    return sig[["regulon", "TF_name", "known_Tcell_TF", "direction",
                "delta_AUC", "mean_AUC_a", "mean_AUC_b", "pval", "FDR"]].reset_index(drop=True)

print("=" * 70)
print("INSPECTION TABLE — Expanded vs Not_expanded")
print("Review delta_AUC and known_Tcell_TF to select regulons for lollipop")
print("=" * 70)

inspect_tables = {}
for lineage in LINEAGES:
    for contrast in ["Expanded_vs_NotExpanded", "Erosion_vs_Rupture_in_Expanded"]:
        key = f"{lineage}__{contrast}"
        if key not in diff_results:
            continue
        tbl = make_inspection_table(diff_results[key])
        inspect_tables[key] = tbl

        print(f"\n── {key}")
        print(f"   Significant regulons: {len(tbl)} "
              f"(known T cell TFs: {tbl['known_Tcell_TF'].sum()})")
        print(tbl[tbl["known_Tcell_TF"]].to_string(index=False))   # show known TFs first

# Export inspection tables
with pd.ExcelWriter(os.path.join(SCENIC_DIR, "SCENIC_inspection_tables.xlsx")) as writer:
    for key, tbl in inspect_tables.items():
        tbl.to_excel(writer, sheet_name=key[:31], index=False)

print(f"\nInspection tables exported to:\n  {SCENIC_DIR}/SCENIC_inspection_tables.xlsx")
print("\nRegulons shown in the lollipop plot are listed in SELECTED_REGULONS below.")

# %% [markdown]
# ## 7. Lollipop plot — curated regulons
#
# Regulon names selected from the inspection table
# you want to include. One list per contrast. The plot is built automatically.

# %%
# ─────────────────────────────────────────────────────────────────────────────
# Selected from the inspection table (section 6)
# Example format: "TF_NAME(+)"
# ─────────────────────────────────────────────────────────────────────────────
SELECTED = {
    "CD8_all__Expanded_vs_NotExpanded": [
        # e.g. "TBX21(+)", "TOX(+)", "EOMES(+)",
    ],
    "CD4_all__Expanded_vs_NotExpanded": [
        # e.g. "BATF(+)", "IRF4(+)", "FOXP3(+)",
    ],
    "CD8_all__Erosion_vs_Rupture_in_Expanded": [
        # e.g. "TCF7(+)", "NR4A1(+)",
    ],
    "CD4_all__Erosion_vs_Rupture_in_Expanded": [
        # e.g. "GATA3(+)", "RORC(+)",
    ],
}

# ─────────────────────────────────────────────────────────────────────────────
# Color scheme for lollipop dots
# positive delta_AUC = higher in group_a (Expanded or Erosion)
# negative delta_AUC = higher in group_b (Not_expanded or Rupture)
# ─────────────────────────────────────────────────────────────────────────────
CONTRAST_COLORS = {
    "Expanded_vs_NotExpanded":         (COLOR_EXP, COLOR_NOT),
    "Erosion_vs_Rupture_in_Expanded":  (COLOR_ERO, COLOR_RUP),
}

def lollipop_plot(diff_df, selected_regulons, title, group_a, group_b,
                  color_a, color_b, save_path=None):
    """
    Horizontal lollipop plot of delta AUC for selected regulons.
    Dot color indicates which group the regulon is enriched in.
    Dot size scales with -log10(FDR).
    """
    if len(selected_regulons) == 0:
        print(f"No regulons selected for: {title}\nCheck the SELECTED dict above.")
        return

    plot_df = diff_df[diff_df["regulon"].isin(selected_regulons)].copy()
    plot_df = plot_df.sort_values("delta_AUC")
    plot_df["neg_log10_FDR"] = -np.log10(plot_df["FDR"].clip(lower=1e-300))
    plot_df["color"] = np.where(plot_df["delta_AUC"] > 0, color_a, color_b)

    fig, ax = plt.subplots(figsize=(6, 0.4 * len(plot_df) + 1.5), dpi=150)

    # Stems
    for _, row in plot_df.iterrows():
        ax.plot([0, row["delta_AUC"]], [row["regulon"], row["regulon"]],
                color="grey", linewidth=0.8, zorder=1)

    # Dots — size proportional to -log10(FDR)
    scatter = ax.scatter(
        plot_df["delta_AUC"],
        plot_df["regulon"],
        c=plot_df["color"],
        s=plot_df["neg_log10_FDR"] * 15 + 20,
        zorder=2,
        edgecolors="white",
        linewidths=0.4,
    )

    # Zero line
    ax.axvline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)

    # FDR significance stars
    for _, row in plot_df.iterrows():
        stars = "***" if row["FDR"] < 0.001 else ("**" if row["FDR"] < 0.01 else
                ("*"  if row["FDR"] < 0.05  else ""))
        if stars:
            x_pos = row["delta_AUC"] + (0.0003 if row["delta_AUC"] > 0 else -0.0003)
            ha = "left" if row["delta_AUC"] > 0 else "right"
            ax.text(x_pos, row["regulon"], stars, va="center", ha=ha,
                    fontsize=7, color="black")

    # Legend for dot size
    sizes = [20 + 15 * s for s in [1, 2, 3]]
    labels = ["-log10(FDR) = 1", "= 2", "= 3"]
    handles = [plt.scatter([], [], s=s, color="grey", alpha=0.6) for s in sizes]
    leg1 = ax.legend(handles, labels, title="Effect size", loc="lower right",
                     fontsize=7, title_fontsize=7, framealpha=0.5)

    # Legend for color
    patch_a = mpatches.Patch(color=color_a, label=f"↑ {group_a}")
    patch_b = mpatches.Patch(color=color_b, label=f"↑ {group_b}")
    ax.legend(handles=[patch_a, patch_b], loc="upper left",
              fontsize=7, framealpha=0.5)
    ax.add_artist(leg1)

    ax.set_xlabel(f"Δ mean AUC ({group_a} − {group_b})", fontsize=9)
    ax.set_ylabel("")
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Saved: {save_path}")
    plt.show()
    return fig

# ── Plot all four contrasts ──────────────────────────────────────────────────
contrast_meta = {
    "Expanded_vs_NotExpanded":        ("Expanded", "Not_expanded"),
    "Erosion_vs_Rupture_in_Expanded": ("Erosion",  "Rupture"),
}

for lineage in LINEAGES:
    for contrast, (group_a, group_b) in contrast_meta.items():
        key = f"{lineage}__{contrast}"
        if key not in diff_results:
            continue

        colors = CONTRAST_COLORS[contrast]
        lollipop_plot(
            diff_df           = diff_results[key],
            selected_regulons = SELECTED.get(key, []),
            title             = f"{lineage} — {contrast.replace('_', ' ')}",
            group_a           = group_a,
            group_b           = group_b,
            color_a           = colors[0],
            color_b           = colors[1],
            save_path         = os.path.join(
                PLOT_DIR, f"lollipop_{key}.pdf"
            ),
        )

# %% [markdown]
# ## 8. RSS plot — regulon specificity per cell type annotation
# Replicates your ApoB workflow: RSS scores across T cell subtype annotations.
# Useful as a quality check and for the supplementary figure.

# %%
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss

for lineage in LINEAGES:
    if lineage not in auc:
        continue

    rss = regulon_specificity_scores(auc[lineage], meta[lineage]["annotation"])
    rss.to_csv(os.path.join(SCENIC_DIR, f"{lineage}_rss_by_annotation.csv"))

    # Get unique annotations present in this lineage
    annotations = meta[lineage]["annotation"].unique().tolist()
    n = len(annotations)
    ncols = min(4, n)
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 3.5, nrows * 3.5), dpi=100)
    axes = np.array(axes).flatten()

    for i, ann in enumerate(sorted(annotations)):
        plot_rss(rss, ann, ax=axes[i])
        axes[i].set_title(ann, fontsize=8)
        if i % ncols != 0:
            axes[i].set_ylabel("")

    # Hide empty axes
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle(f"RSS by annotation — {lineage}", fontsize=10, y=1.01)
    plt.tight_layout()
    save_rss = os.path.join(PLOT_DIR, f"RSS_by_annotation_{lineage}.pdf")
    fig.savefig(save_rss, dpi=300, bbox_inches="tight")
    plt.show()
    print(f"RSS plot saved: {save_rss}")

# %% [markdown]
# ## 9. RSS plot — regulon specificity by expansion group
# Adapts the RSS approach to score regulons by expansion_group rather than
# cell type. Useful as a high-level QC before the differential analysis.

# %%
for lineage in LINEAGES:
    if lineage not in auc:
        continue

    rss_exp = regulon_specificity_scores(auc[lineage], meta[lineage]["expansion_group"])
    rss_exp.to_csv(os.path.join(SCENIC_DIR, f"{lineage}_rss_by_expansion.csv"))

    fig, axes = plt.subplots(1, 2, figsize=(8, 4), dpi=100)
    for i, grp in enumerate(["Expanded", "Not_expanded"]):
        if grp in rss_exp.index:
            plot_rss(rss_exp, grp, ax=axes[i])
            axes[i].set_title(f"{lineage} — {grp}", fontsize=9)
            if i > 0:
                axes[i].set_ylabel("")

    plt.tight_layout()
    save_rss2 = os.path.join(PLOT_DIR, f"RSS_by_expansion_{lineage}.pdf")
    fig.savefig(save_rss2, dpi=300, bbox_inches="tight")
    plt.show()
    print(f"RSS plot saved: {save_rss2}")
