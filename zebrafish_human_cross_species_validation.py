# %% [markdown]
# # Zebrafish-to-human hematopoietic growth-program validation
#
# This notebook validates the **ortholog-mapped proliferation and apoptosis expression
# programs** from the saved zebrafish analysis. It does **not** rerun GraphVelo, moscot,
# or the previous moscot-like transport code.
#
# Primary human reference:
# [Popescu et al., Nature 2019](https://doi.org/10.1038/s41586-019-1652-y),
# annotated 10x fetal-liver cells from E-MTAB-7407 (14 donors, 7+6 to 17 PCW).
#
# Secondary replication reference:
# [GSE189161](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189161),
# CD34-enriched HSPCs from 26 human samples. Only the nine fetal-liver donors are used
# for the fetal replication test; the full lifespan cohort is retained for sensitivity
# analyses.
#
# ## What can and cannot be concluded
#
# A concordant result supports conservation of an ortholog-mapped **transcriptomic
# program**. It does not prove that a zebrafish transport coupling, fate probability,
# velocity vector, or absolute growth rate is numerically transferable to humans.
# GSE189161 is a strong *independent HSPC replication dataset*, but it is not a broad
# blood/immune atlas because cells were CD34 enriched.

# %%
from __future__ import annotations

import gzip
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
import time
import warnings
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_cross_species_validation")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/xdg_cross_species_validation")

import anndata as ad
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from IPython.display import Markdown, display
from scipy import sparse
from scipy.stats import binomtest, pearsonr, spearmanr
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=FutureWarning)
SEED = 322
RNG = np.random.default_rng(SEED)

mpl.rcParams.update(
    {
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "text.color": "black",
        "axes.edgecolor": "black",
        "axes.labelcolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 10,
    }
)

EXPANDED_MAPPING_NAME = "reviewed_zebrafish_to_human_growth_orthologs_expanded.csv"

ROOT = Path(os.environ.get("ZEBRAFISH_PROJECT_DIR", str(Path.cwd()))).expanduser().resolve()
DATA_ROOT = Path(os.environ.get("ZEBRAFISH_DATA_DIR", str(ROOT))).expanduser()

VALIDATION_DIR = ROOT / "cross_species_validation"
DATA_DIR = VALIDATION_DIR / "data"
CACHE_DIR = VALIDATION_DIR / "cache"
RESULTS_DIR = VALIDATION_DIR / "results"
FIGURES_DIR = VALIDATION_DIR / "figures"
for directory in (CACHE_DIR, RESULTS_DIR, FIGURES_DIR):
    directory.mkdir(parents=True, exist_ok=True)

ZF_CHECKPOINT = (
    ROOT
    / "blood_growth_graphvelo_moscot_v2/checkpoints/"
    / "zebrafish_growth_graphvelo_moscot_v2.h5ad"
)
ZF_FULL = DATA_ROOT / "zf_atlas_hematopoetic_endothelial_v4_release.h5ad"
ORTHOLOG_CSV = ROOT / EXPANDED_MAPPING_NAME
MARKER_SOURCE_MANIFEST = ROOT / "growth_marker_reference/source_manifest.csv"
ZFIN_ORTHOLOG_DUMP = ROOT / "human_orthos.txt"
POPESCU_H5AD = DATA_DIR / "popescu_E_MTAB_7407/fetal_liver_alladata.h5ad"
GSE_DIR = DATA_DIR / "GSE189161"
GSE_MATRIX = GSE_DIR / "GSE189161_matrix.mtx.gz"
GSE_FEATURES = GSE_DIR / "GSE189161_features.txt.gz"
GSE_CELLS = GSE_DIR / "GSE189161_cells.txt.gz"
GSE_METADATA = GSE_DIR / "GSE189161_metadata.csv.gz"
GSE_SERIES = GSE_DIR / "GSE189161_series_matrix.txt.gz"

REQUIRED_FILES = [
    ZF_CHECKPOINT,
    ZF_FULL,
    ORTHOLOG_CSV,
    ZFIN_ORTHOLOG_DUMP,
    MARKER_SOURCE_MANIFEST,
    POPESCU_H5AD,
    GSE_MATRIX,
    GSE_FEATURES,
    GSE_CELLS,
    GSE_METADATA,
    GSE_SERIES,
]
missing = [str(p) for p in REQUIRED_FILES if not p.exists()]
assert not missing, "Missing required files:\n" + "\n".join(missing)

print("Python:", sys.version.split()[0])
print("Executable:", sys.executable)
print("anndata:", ad.__version__)
print("numpy / scipy / pandas:", np.__version__, scipy.__version__, pd.__version__)
print("Root:", ROOT)
print("This notebook imports neither moscot nor GraphVelo.")

# %% [markdown]
# ## 1. Statistical target and formulas
#
# Counts are library-size normalized within each dataset:
#
# $$
# L_{ig}=\log\left(1+10^4\frac{C_{ig}}{\sum_h C_{ih}}\right).
# $$
#
# To avoid comparing incompatible absolute expression scales across platforms and
# species, each gene is standardized *within species*:
#
# $$
# Z_{ig}=\mathrm{clip}\left(\frac{L_{ig}-\mu_g}{\sigma_g},-5,5\right),
# \qquad
# S_i^{(M)}=\frac{1}{|M|}\sum_{g\in M}Z_{ig}.
# $$
#
# Human paralogs that map to multiple zebrafish genes are collapsed by their mean
# normalized expression before scoring. Formal time tests use replicate-by-lineage
# pseudobulks and a lineage-adjusted model:
#
# $$
# S_{rb}=\beta_0+\beta_t\,\widetilde{t}_r+
# \sum_{k=2}^{K}\gamma_k I(b=k)+\varepsilon_{rb},
# $$
#
# with cluster-robust uncertainty at the fish/donor level. Cell-level distributions
# are descriptive only. Multiple tested programs use Benjamini-Hochberg FDR.

# %%
def human_bytes(n: int) -> str:
    value = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024 or unit == "TB":
            return f"{value:.2f} {unit}"
        value /= 1024
    return f"{value:.2f} TB"


def short_sha256(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


audit_rows = []
for path in REQUIRED_FILES:
    audit_rows.append(
        {
            "file": str(path.relative_to(ROOT)),
            "bytes": path.stat().st_size,
            "size": human_bytes(path.stat().st_size),
            "sha256": short_sha256(path),
        }
    )
file_audit = pd.DataFrame(audit_rows)
file_audit.to_csv(RESULTS_DIR / "input_file_audit.csv", index=False)
display(file_audit)

# %% [markdown]
# ## 2. Externally supported ortholog map and phase-resolved marker sets
#
# The complete zebrafish gene universe is **not** treated as a proliferation list.
# The frozen expanded panel requires external cell-cycle/pathway evidence, a ZFIN
# human-zebrafish orthology record, and detection in the zebrafish analysis. S-phase
# and G2/M genes are read from the mapping table rather than duplicated as hard-coded
# lists here.

# %%
ortholog_map = pd.read_csv(ORTHOLOG_CSV)
required_columns = {
    "zebrafish_gene",
    "human_gene",
    "marker_class",
    "marker_subclass",
    "mapping_database",
    "evidence_status",
    "pathway_sources",
    "selection_rule",
}
assert required_columns.issubset(ortholog_map.columns)
assert ortholog_map["evidence_status"].eq("reviewed_external_consensus").all()

ortholog_map["zebrafish_gene"] = ortholog_map["zebrafish_gene"].astype(str)
ortholog_map["human_gene"] = ortholog_map["human_gene"].astype(str).str.upper()

# Independent local cross-check against the downloaded ZFIN human orthology table.
# Columns 1 and 3 (zero-based) are the zebrafish and human gene symbols.
zfin_orthologs = pd.read_csv(ZFIN_ORTHOLOG_DUMP, sep="\t", header=None, dtype=str)
zfin_pairs = set(
    zip(
        zfin_orthologs.iloc[:, 1].str.casefold(),
        zfin_orthologs.iloc[:, 3].str.upper(),
    )
)
ortholog_map["present_in_zfin_dump"] = [
    (zf.casefold(), human.upper()) in zfin_pairs
    for zf, human in zip(ortholog_map.zebrafish_gene, ortholog_map.human_gene)
]
assert ortholog_map["present_in_zfin_dump"].all(), (
    "At least one reviewed mapping is absent from the local ZFIN orthology dump."
)
ortholog_map.to_csv(RESULTS_DIR / "expanded_ortholog_mapping_audit.csv", index=False)

human_genes = ortholog_map["human_gene"].drop_duplicates().tolist()

PROGRAMS = {
    "proliferation": ortholog_map.loc[
        ortholog_map.marker_class.eq("proliferation"), "human_gene"
    ].drop_duplicates().tolist(),
    "S_phase": ortholog_map.loc[
        ortholog_map.marker_subclass.eq("S_phase"), "human_gene"
    ].drop_duplicates().tolist(),
    "G2M_phase": ortholog_map.loc[
        ortholog_map.marker_subclass.eq("G2M_phase"), "human_gene"
    ].drop_duplicates().tolist(),
    "apoptosis": ortholog_map.loc[
        ortholog_map.marker_class.eq("apoptosis"), "human_gene"
    ].drop_duplicates().tolist(),
}

program_sizes = pd.Series({k: len(v) for k, v in PROGRAMS.items()}, name="n_genes")
display(program_sizes.to_frame())
display(ortholog_map)

# %% [markdown]
# ## 3. Memory-safe marker extraction
#
# The annotated Popescu H5AD stores raw counts for 113,063 cells. Loading its entire
# 27,080-gene matrix is unnecessary. The helper below reads only the reviewed marker
# columns from HDF5 CSR matrices in row chunks.

# %%
def _decode_array(values) -> list[str]:
    return [v.decode("utf-8") if isinstance(v, (bytes, np.bytes_)) else str(v) for v in values]


def _axis_names(h5: h5py.File, axis: str) -> list[str]:
    node = h5[axis]
    if isinstance(node, h5py.Dataset):
        assert node.dtype.names and "index" in node.dtype.names
        return _decode_array(node["index"])
    index_key = node.attrs.get("_index", "_index")
    if isinstance(index_key, bytes):
        index_key = index_key.decode()
    return _decode_array(node[index_key][:])


def extract_h5ad_csr_columns(
    path: Path,
    requested_genes: list[str],
    matrix_key: str = "X",
    chunk_rows: int = 2_000,
) -> tuple[np.ndarray, list[str]]:
    """Return dense n_cells x n_requested marker counts without loading all genes."""
    with h5py.File(path, "r") as h5:
        matrix = h5[matrix_key]
        assert isinstance(matrix, h5py.Group), f"{matrix_key} is not a sparse HDF5 group"
        fmt = matrix.attrs.get("encoding-type", matrix.attrs.get("h5sparse_format", ""))
        if isinstance(fmt, bytes):
            fmt = fmt.decode()
        assert "csr" in str(fmt).lower(), f"Expected CSR matrix; found {fmt!r}"

        shape = matrix.attrs.get("shape", matrix.attrs.get("h5sparse_shape"))
        n_rows, n_cols = map(int, shape)
        var_names = _axis_names(h5, "var")
        assert len(var_names) == n_cols

        exact_lookup: dict[str, int] = {}
        casefold_lookup: dict[str, list[int]] = {}
        for index, gene in enumerate(var_names):
            assert gene not in exact_lookup, f"Duplicate exact feature name: {gene}"
            exact_lookup[gene] = index
            casefold_lookup.setdefault(gene.casefold(), []).append(index)

        # Prefer an exact symbol match. Fall back to case-insensitive matching only
        # when it identifies one feature. This resolves the atlas's CKS2/cks2 pair
        # to the ZFIN symbol `cks2` instead of arbitrarily choosing one feature.
        found = []
        resolved_columns = []
        ambiguous = []
        for gene in requested_genes:
            if gene in exact_lookup:
                found.append(gene)
                resolved_columns.append(exact_lookup[gene])
                continue
            candidates = casefold_lookup.get(gene.casefold(), [])
            if len(candidates) == 1:
                found.append(gene)
                resolved_columns.append(candidates[0])
            elif len(candidates) > 1:
                ambiguous.append(gene)
        assert not ambiguous, (
            "Requested genes lack an exact match and have ambiguous case-insensitive "
            f"matches: {ambiguous}"
        )
        columns = np.asarray(resolved_columns, dtype=np.int64)
        out = np.zeros((n_rows, len(found)), dtype=np.float32)
        indptr = matrix["indptr"][:]

        for start_row in range(0, n_rows, chunk_rows):
            stop_row = min(start_row + chunk_rows, n_rows)
            start_nnz = int(indptr[start_row])
            stop_nnz = int(indptr[stop_row])
            local_indptr = indptr[start_row : stop_row + 1].astype(np.int64) - start_nnz
            local = sparse.csr_matrix(
                (
                    matrix["data"][start_nnz:stop_nnz],
                    matrix["indices"][start_nnz:stop_nnz],
                    local_indptr,
                ),
                shape=(stop_row - start_row, n_cols),
            )
            out[start_row:stop_row] = local[:, columns].toarray()
    return out, found


def log_cpm(counts: np.ndarray | sparse.spmatrix, totals=None) -> np.ndarray:
    if sparse.issparse(counts):
        counts = counts.toarray()
    counts = np.asarray(counts, dtype=np.float64)
    if totals is None:
        totals = counts.sum(axis=1)
    totals = np.asarray(totals, dtype=np.float64).reshape(-1)
    totals = np.where(np.isfinite(totals) & (totals > 0), totals, 1.0)
    return np.log1p(10_000.0 * counts / totals[:, None]).astype(np.float32)


def zscore_columns(values: np.ndarray, fit_mask=None, clip: float = 5.0) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    fit = values if fit_mask is None else values[np.asarray(fit_mask)]
    mean = np.nanmean(fit, axis=0)
    std = np.nanstd(fit, axis=0)
    std = np.where(np.isfinite(std) & (std > 1e-8), std, 1.0)
    return np.clip((values - mean) / std, -clip, clip).astype(np.float32)


def module_scores(z_values: np.ndarray, genes: list[str]) -> pd.DataFrame:
    gene_to_col = {g: i for i, g in enumerate(genes)}
    result = {}
    for program, members in PROGRAMS.items():
        columns = [gene_to_col[g] for g in members if g in gene_to_col]
        assert columns, f"No genes available for {program}"
        result[program] = np.nanmean(z_values[:, columns], axis=1)
    return pd.DataFrame(result)


def collapse_zebrafish_paralogs(
    zf_expression: np.ndarray,
    zf_genes: list[str],
    mapping: pd.DataFrame,
) -> tuple[np.ndarray, list[str]]:
    gene_to_col = {g.casefold(): i for i, g in enumerate(zf_genes)}
    collapsed = []
    labels = []
    for human_gene, block in mapping.groupby("human_gene", sort=False):
        columns = [gene_to_col[g.casefold()] for g in block.zebrafish_gene if g.casefold() in gene_to_col]
        if columns:
            collapsed.append(np.mean(zf_expression[:, columns], axis=1))
            labels.append(human_gene)
    return np.column_stack(collapsed).astype(np.float32), labels

# %% [markdown]
# ## 4. Load the saved zebrafish checkpoint without rerunning transport
#
# The saved checkpoint is the authoritative source for cell inclusion, annotations,
# time, fish IDs, and previous growth scores. Its feature space contains only 2,986
# selected genes, so it cannot support the complete expanded proliferation audit. The
# matching cells are therefore aligned to the original full-count H5AD solely to
# recover the reviewed marker counts. No velocity or transport model is refit.

# %%
zf_checkpoint = ad.read_h5ad(ZF_CHECKPOINT, backed="r")
zf_obs = zf_checkpoint.obs.copy()
zf_checkpoint_cells = zf_checkpoint.obs_names.astype(str).tolist()
zf_checkpoint_shape = zf_checkpoint.shape
zf_checkpoint.file.close()

zf_full_counts, zf_found = extract_h5ad_csr_columns(
    ZF_FULL,
    ortholog_map.zebrafish_gene.drop_duplicates().tolist(),
    matrix_key="layers/counts",
)
zf_source = ad.read_h5ad(ZF_FULL, backed="r")
zf_source_obs = zf_source.obs.copy()
zf_full_cells = zf_source.obs_names.astype(str).tolist()
zf_source.file.close()

row_lookup = {cell: i for i, cell in enumerate(zf_full_cells)}
missing_cells = [cell for cell in zf_checkpoint_cells if cell not in row_lookup]
assert not missing_cells, f"Checkpoint cells absent from the full H5AD: {missing_cells[:5]}"
zf_rows = np.asarray([row_lookup[cell] for cell in zf_checkpoint_cells])
zf_full_counts = zf_full_counts[zf_rows]

# Use the library totals from the same full-count object, not the selected-feature
# checkpoint's spliced-size field.
zf_source_totals = zf_source_obs.iloc[zf_rows]["total_counts"].to_numpy()
zf_log = log_cpm(zf_full_counts, zf_source_totals)
zf_human_log, zf_human_genes = collapse_zebrafish_paralogs(zf_log, zf_found, ortholog_map)
assert set(human_genes) == set(zf_human_genes)

ZF_TYPE_MAP = {
    "hematopoietic_stem_and_progenitor_HSPC": "stem_progenitor",
    "erythroblast": "erythroid",
    "macrophage": "macrophage",
    "macrophage_immature": "macrophage",
    "neutrophil": "neutrophil",
    "innate-lymphoid-cell": "innate_lymphoid",
}
zf_meta = zf_obs.copy()
zf_meta["replicate"] = zf_meta["fish_id"].astype(str)
zf_meta["age"] = pd.to_numeric(zf_meta["hpf"], errors="coerce")
zf_meta["broad_type"] = zf_meta["official_fine_cell_type"].map(ZF_TYPE_MAP)
zf_keep = zf_meta["broad_type"].notna().to_numpy()

zf_z = zscore_columns(zf_human_log, fit_mask=zf_keep)
zf_scores = module_scores(zf_z, zf_human_genes)
zf_scores.index = zf_meta.index
for column in zf_scores:
    zf_meta[column] = zf_scores[column]

checkpoint_score_audit = []
for previous, recomputed in [
    ("zf_proliferation_score", "proliferation"),
    ("zf_apoptosis_score", "apoptosis"),
]:
    valid = zf_meta[previous].notna() & zf_meta[recomputed].notna()
    rho, pvalue = spearmanr(zf_meta.loc[valid, previous], zf_meta.loc[valid, recomputed])
    checkpoint_score_audit.append(
        {"saved_score": previous, "recomputed_score": recomputed, "spearman_rho": rho, "pvalue": pvalue}
    )
checkpoint_score_audit = pd.DataFrame(checkpoint_score_audit)
checkpoint_score_audit.to_csv(RESULTS_DIR / "zebrafish_checkpoint_score_audit.csv", index=False)

print("Saved checkpoint shape:", zf_checkpoint_shape)
print("Aligned cells:", len(zf_meta))
print("Recovered zebrafish marker genes:", len(zf_found), "/", ortholog_map.zebrafish_gene.nunique())
display(checkpoint_score_audit)

# %% [markdown]
# ## 5. Load annotated Popescu E-MTAB-7407 fetal-liver cells
#
# The official portal file contains 113,063 annotated 10x fetal-liver cells from 14
# donors. `orig.ident` is used as the donor unit. Gestational age is parsed from
# `fetal.ids`, including day offsets such as `16+2PCW`.

# %%
popescu = ad.read_h5ad(POPESCU_H5AD, backed="r")
pop_obs = popescu.obs.copy()
pop_shape = popescu.shape
popescu.file.close()

pop_counts, pop_found = extract_h5ad_csr_columns(POPESCU_H5AD, human_genes, matrix_key="X")
assert set(pop_found) == set(human_genes)

def parse_pcw(value: str) -> float:
    text = str(value)
    match = re.search(r"_(\d+)(?:\+(\d+))?PCW", text, flags=re.IGNORECASE)
    if not match:
        raise ValueError(f"Could not parse gestational age from {value!r}")
    weeks = float(match.group(1))
    days = float(match.group(2) or 0)
    return weeks + days / 7.0


POP_TYPE_MAP = {
    "HSC_MPP": "stem_progenitor",
    "Early Erythroid": "erythroid",
    "Mid Erythroid": "erythroid",
    "Late Erythroid": "erythroid",
    "MEMP": "erythroid",
    "Kupffer Cell": "macrophage",
    "Mono-Mac": "macrophage",
    "Monocyte": "macrophage",
    "Monocyte precursor": "macrophage",
    "VCAM1+ EI macrophage": "macrophage",
    "Neutrophil-myeloid progenitor": "neutrophil",
    "NK": "innate_lymphoid",
    "ILC precursor": "innate_lymphoid",
}

pop_meta = pop_obs.copy()
pop_meta["replicate"] = pop_meta["orig.ident"].astype(str)
pop_meta["age"] = pop_meta["fetal.ids"].map(parse_pcw)
pop_meta["broad_type"] = pop_meta["cell.labels"].map(POP_TYPE_MAP)
pop_keep = pop_meta["broad_type"].notna().to_numpy()
pop_log = log_cpm(pop_counts, pop_meta["n_counts"].to_numpy())
pop_z = zscore_columns(pop_log, fit_mask=pop_keep)
pop_scores = module_scores(pop_z, pop_found)
pop_scores.index = pop_meta.index
for column in pop_scores:
    pop_meta[column] = pop_scores[column]

print("Popescu shape:", pop_shape)
print("Donors:", pop_meta.replicate.nunique())
print("Age range (PCW):", pop_meta.age.min(), "to", pop_meta.age.max())
print("Mapped hematopoietic cells:", int(pop_keep.sum()))
display(
    pop_meta.loc[pop_keep]
    .groupby(["replicate", "age", "broad_type"], observed=True)
    .size()
    .rename("n_cells")
    .reset_index()
    .head(20)
)

# %% [markdown]
# ## 6. Load GSE189161 marker counts with a reusable cache
#
# GSE189161 contains a very large Matrix Market file (130,073,329 nonzero entries).
# The first execution scans it once with `gzip` and `awk`, retains the current mapped
# human growth-program genes, and saves a compact sparse cache. The cache filename
# includes a marker-panel hash, so an old 38/12 panel can never be reused silently.
# This does not use the network and does not rerun the zebrafish notebook.

# %%
def build_gse_marker_cache(
    matrix_path: Path,
    feature_path: Path,
    target_genes: list[str],
    cache_path: Path,
) -> sparse.csr_matrix:
    features = pd.read_csv(feature_path, header=None, names=["gene"])["gene"].astype(str)
    feature_lookup = {gene.casefold(): i for i, gene in enumerate(features)}
    found = [gene for gene in target_genes if gene.casefold() in feature_lookup]
    assert set(found) == set(target_genes), sorted(set(target_genes) - set(found))

    row_map_path = CACHE_DIR / f"gse189161_marker_row_map_{marker_cache_tag}.tsv"
    triplet_path = CACHE_DIR / f"gse189161_growth_marker_triplets_{marker_cache_tag}.tsv"
    row_map = pd.DataFrame(
        {
            "matrix_row_1based": [feature_lookup[g.casefold()] + 1 for g in target_genes],
            "target_col_1based": np.arange(1, len(target_genes) + 1),
        }
    )
    row_map.to_csv(row_map_path, sep="\t", index=False, header=False)

    awk_program = (
        'NR==FNR {keep[$1]=$2; next} '
        'FNR<=3 {next} '
        '($1 in keep) {print ($2-1) "\\t" (keep[$1]-1) "\\t" $3}'
    )
    print("Scanning GSE189161 Matrix Market file once; this can take several minutes...")
    started = time.time()
    with triplet_path.open("wb") as output:
        unzip = subprocess.Popen(["gzip", "-cd", str(matrix_path)], stdout=subprocess.PIPE)
        assert unzip.stdout is not None
        awk = subprocess.run(
            ["awk", awk_program, str(row_map_path), "-"],
            stdin=unzip.stdout,
            stdout=output,
            stderr=subprocess.PIPE,
            check=False,
        )
        unzip.stdout.close()
        unzip_code = unzip.wait()
    if awk.returncode != 0 or unzip_code != 0:
        raise RuntimeError(
            f"GSE marker extraction failed (gzip={unzip_code}, awk={awk.returncode}): "
            + awk.stderr.decode(errors="replace")
        )

    triplets = pd.read_csv(
        triplet_path,
        sep="\t",
        header=None,
        names=["cell", "gene", "count"],
        dtype={"cell": np.int32, "gene": np.int16, "count": np.float32},
    )
    cells = pd.read_csv(GSE_CELLS, header=None).shape[0]
    matrix = sparse.csr_matrix(
        (triplets["count"], (triplets["cell"], triplets["gene"])),
        shape=(cells, len(target_genes)),
        dtype=np.float32,
    )
    sparse.save_npz(cache_path, matrix, compressed=True)
    pd.Series(target_genes, name="human_gene").to_csv(
        CACHE_DIR / f"gse189161_growth_marker_genes_{marker_cache_tag}.csv", index=False
    )
    triplet_path.unlink()
    print(f"GSE marker cache built in {(time.time() - started) / 60:.1f} minutes")
    return matrix


marker_cache_tag = hashlib.sha256("\n".join(human_genes).encode()).hexdigest()[:12]
gse_cache = CACHE_DIR / f"gse189161_growth_marker_counts_{marker_cache_tag}.npz"
gse_gene_cache = CACHE_DIR / f"gse189161_growth_marker_genes_{marker_cache_tag}.csv"
if gse_cache.exists():
    gse_counts = sparse.load_npz(gse_cache).tocsr()
    cached_genes = pd.read_csv(gse_gene_cache)["human_gene"].tolist()
    assert cached_genes == human_genes
else:
    gse_counts = build_gse_marker_cache(GSE_MATRIX, GSE_FEATURES, human_genes, gse_cache)

gse_meta = pd.read_csv(GSE_METADATA, index_col=0)
gse_cells = pd.read_csv(GSE_CELLS, header=None)[0].astype(str)
assert gse_cells.tolist() == gse_meta["Cell"].astype(str).tolist()
assert gse_counts.shape == (len(gse_meta), len(human_genes))

GSE_SAMPLE_AGE = {
    "Samp21": ("fetal_liver", 10.0), "GR97": ("fetal_liver", 11.0),
    "Samp22": ("fetal_liver", 13.0), "Samp23": ("fetal_liver", 14.0),
    "Samp41": ("fetal_liver", 15.0), "Samp42": ("fetal_liver", 20.0),
    "GR72": ("fetal_liver", 22.0), "GR73": ("fetal_liver", 23.0),
    "GR74": ("fetal_liver", 23.0), "CB1": ("cord_blood", np.nan),
    "CB2": ("cord_blood", np.nan), "GR95": ("bone_marrow", 2.0 * 52),
    "Samp44": ("bone_marrow", 4.0 * 52), "GR96": ("bone_marrow", 10.0 * 52),
    "Samp43": ("bone_marrow", 12.0 * 52), "GR90": ("bone_marrow", 17.0 * 52),
    "BM1": ("bone_marrow", 25.0 * 52), "BM3": ("bone_marrow", 25.0 * 52),
    "BM2": ("bone_marrow", 32.0 * 52), "Samp24": ("bone_marrow", 35.0 * 52),
    "Samp25": ("bone_marrow", 45.0 * 52), "Samp26": ("bone_marrow", 53.0 * 52),
    "GR91": ("bone_marrow", 62.0 * 52), "GR94": ("bone_marrow", 62.0 * 52),
    "GR92": ("bone_marrow", 76.0 * 52), "GR93": ("bone_marrow", 77.0 * 52),
}
GSE_TYPE_MAP = {
    "HSC": "stem_progenitor", "MPP-1": "stem_progenitor",
    "MPP-2": "stem_progenitor", "MPP-3": "stem_progenitor",
    "E-Prog-1": "erythroid", "E-Prog-2": "erythroid", "E-Prog-3": "erythroid",
    "Mk/E-MPP": "erythroid", "Mk-Prog": "erythroid",
    "My-MPP": "myeloid", "G-Prog-1": "myeloid", "G-Prog-2": "myeloid",
    "G-Prog-3": "myeloid", "G/Mono-Prog": "myeloid",
    "Mono/cDC-Prog": "myeloid", "pDC-Prog": "myeloid",
    "Baso/Mast-Prog": "myeloid",
    "Ly-Prog-1": "lymphoid", "Ly-Prog-2": "lymphoid", "Ly-Prog-3": "lymphoid",
    "Ly-Prog-4": "lymphoid", "Ly-Prog-5": "lymphoid",
}

gse_meta["replicate"] = gse_meta["orig.ident"].astype(str)
gse_meta["tissue"] = gse_meta["replicate"].map(lambda x: GSE_SAMPLE_AGE[x][0])
gse_meta["age"] = gse_meta["replicate"].map(lambda x: GSE_SAMPLE_AGE[x][1])
gse_meta["broad_type"] = gse_meta["cluster_name"].map(GSE_TYPE_MAP)
gse_fetal = gse_meta["tissue"].eq("fetal_liver").to_numpy()

gse_log = log_cpm(gse_counts, gse_meta["nCount_RNA"].to_numpy())
gse_z = zscore_columns(gse_log, fit_mask=gse_fetal)
gse_scores = module_scores(gse_z, human_genes)
gse_scores.index = gse_meta.index
for column in gse_scores:
    gse_meta[column] = gse_scores[column]

print("GSE189161 cells:", len(gse_meta))
print("Samples:", gse_meta.replicate.nunique())
print("Fetal-liver donors/cells:", gse_meta.loc[gse_fetal, "replicate"].nunique(), int(gse_fetal.sum()))
display(gse_meta.loc[gse_fetal].groupby(["replicate", "age", "broad_type"]).size().unstack(fill_value=0))

# %% [markdown]
# ## 7. Coverage and harmonization audit
#
# Cell types are harmonized conservatively. Ambiguous zebrafish
# `non_specific_hematopoietic`, brain microglia, and human non-haematopoietic liver
# cells are excluded from the primary lineage comparison. GSE189161 labels denote
# progenitor biases—not mature immune-cell identities—so that cohort is not used as
# a substitute for the broad Popescu atlas.

# %%
coverage_rows = []
for dataset, available in {
    "zebrafish_full_counts": set(g.upper() for g in zf_found),
    "Popescu_E-MTAB-7407": set(g.upper() for g in pop_found),
    "GSE189161": set(human_genes),
}.items():
    for marker_class, block in ortholog_map.groupby("marker_class"):
        target = (
            set(block.zebrafish_gene.str.upper())
            if dataset == "zebrafish_full_counts"
            else set(block.human_gene.str.upper())
        )
        coverage_rows.append(
            {
                "dataset": dataset,
                "marker_class": marker_class,
                "found": len(target.intersection(available)),
                "required": len(target),
                "fraction": len(target.intersection(available)) / len(target),
                "missing": ";".join(sorted(target - available)),
            }
        )
coverage = pd.DataFrame(coverage_rows)
coverage.to_csv(RESULTS_DIR / "marker_coverage.csv", index=False)

type_audit = pd.concat(
    [
        pd.DataFrame({"dataset": "zebrafish", "original_label": list(ZF_TYPE_MAP), "broad_type": list(ZF_TYPE_MAP.values())}),
        pd.DataFrame({"dataset": "Popescu", "original_label": list(POP_TYPE_MAP), "broad_type": list(POP_TYPE_MAP.values())}),
        pd.DataFrame({"dataset": "GSE189161", "original_label": list(GSE_TYPE_MAP), "broad_type": list(GSE_TYPE_MAP.values())}),
    ],
    ignore_index=True,
)
type_audit.to_csv(RESULTS_DIR / "cell_type_harmonization.csv", index=False)
display(coverage)
display(type_audit)

# %% [markdown]
# ## 8. Replicate-by-lineage pseudobulks
#
# Pseudobulk here means the mean log-normalized expression or mean standardized
# program score for a fish/donor and broad cell class. A minimum cell count is applied
# before inference. This prevents a donor with thousands of cells from masquerading as
# thousands of biological replicates.

# %%
def make_pseudobulk(
    meta: pd.DataFrame,
    expression: np.ndarray,
    genes: list[str],
    keep_mask: np.ndarray,
    min_cells: int,
    dataset: str,
) -> pd.DataFrame:
    assert len(meta) == expression.shape[0]
    positions = np.flatnonzero(keep_mask)
    grouping = meta.iloc[positions][["replicate", "age", "broad_type"]].copy()
    grouping["position"] = positions
    score_columns = list(PROGRAMS)
    rows = []
    for (replicate, age, broad_type), block in grouping.groupby(
        ["replicate", "age", "broad_type"], observed=True, dropna=False
    ):
        if len(block) < min_cells:
            continue
        idx = block["position"].to_numpy(dtype=int)
        row = {
            "dataset": dataset,
            "replicate": str(replicate),
            "age": float(age),
            "broad_type": str(broad_type),
            "n_cells": len(idx),
        }
        row.update({g: float(np.mean(expression[idx, j])) for j, g in enumerate(genes)})
        row.update({program: float(meta.iloc[idx][program].mean()) for program in score_columns})
        rows.append(row)
    return pd.DataFrame(rows)


zf_pb = make_pseudobulk(zf_meta, zf_human_log, zf_human_genes, zf_keep, 5, "zebrafish")
pop_pb = make_pseudobulk(pop_meta, pop_log, pop_found, pop_keep, 20, "Popescu")
gse_pb = make_pseudobulk(gse_meta, gse_log, human_genes, gse_fetal, 20, "GSE189161_fetal")

for frame, filename in [
    (zf_pb, "zebrafish_replicate_lineage_pseudobulk.csv"),
    (pop_pb, "popescu_replicate_lineage_pseudobulk.csv"),
    (gse_pb, "gse189161_fetal_replicate_lineage_pseudobulk.csv"),
]:
    frame.to_csv(RESULTS_DIR / filename, index=False)

pseudobulk_audit = pd.concat(
    [
        x.groupby("dataset").agg(
            pseudobulks=("replicate", "size"),
            replicates=("replicate", "nunique"),
            lineages=("broad_type", "nunique"),
            cells=("n_cells", "sum"),
        )
        for x in (zf_pb, pop_pb, gse_pb)
    ]
).reset_index()
display(pseudobulk_audit)

# %% [markdown]
# ## 9. Primary test: conserved profiles across homologous blood/immune classes
#
# For a common set of broad lineages, donor/fish pseudobulks are averaged and each gene
# is standardized across lineages. Concordance is computed on the resulting
# gene-by-lineage profiles. The null independently permutes human gene identities
# within proliferation and apoptosis classes:
#
# $$
# p_{\mathrm{perm}}=
# \frac{1+\sum_{b=1}^{B}I\!\left(r_b\ge r_{\mathrm{obs}}\right)}{B+1}.
# $$
#
# Donor/fish bootstrap intervals quantify replicate sensitivity. The permutation null
# asks whether the exact ortholog pairing contributes more than generic module-level
# co-expression; it is not a causal test.

# %%
def standardized_type_profile(pb: pd.DataFrame, types: list[str], genes: list[str]) -> np.ndarray:
    profile = pb[pb.broad_type.isin(types)].groupby("broad_type")[genes].mean().loc[types].to_numpy()
    mean = profile.mean(axis=0, keepdims=True)
    std = profile.std(axis=0, keepdims=True)
    std[std < 1e-8] = 1.0
    return (profile - mean) / std


def resampled_type_profile(
    pb: pd.DataFrame,
    types: list[str],
    genes: list[str],
    rng: np.random.Generator,
) -> np.ndarray:
    # Resample biological-replicate pseudobulks *within each lineage*.  A global
    # replicate bootstrap can omit a rare lineage completely when that lineage is
    # observed in only a subset of fish/donors, making the requested profile
    # undefined.  Stratification keeps the estimand fixed while still propagating
    # replicate-level uncertainty.
    pieces = []
    for broad_type in types:
        block = pb[pb.broad_type.eq(broad_type)]
        assert len(block) > 0, f"No pseudobulks available for {broad_type}"
        sampled_rows = rng.choice(block.index.to_numpy(), size=len(block), replace=True)
        sampled = pb.loc[sampled_rows].copy()
        sampled["bootstrap_draw"] = np.arange(len(sampled))
        pieces.append(sampled)
    sampled_pb = pd.concat(pieces, ignore_index=True)
    return standardized_type_profile(sampled_pb, types, genes)


def profile_concordance(
    left: pd.DataFrame,
    right: pd.DataFrame,
    comparison: str,
    programs=("proliferation", "apoptosis"),
    permutations: int = 1_000,
    bootstraps: int = 500,
) -> pd.DataFrame:
    common_types = sorted(set(left.broad_type).intersection(right.broad_type))
    results = []
    for program in programs:
        genes = [g for g in PROGRAMS[program] if g in left.columns and g in right.columns]
        left_profile = standardized_type_profile(left, common_types, genes)
        right_profile = standardized_type_profile(right, common_types, genes)
        x = left_profile.ravel()
        y = right_profile.ravel()
        observed_spearman = spearmanr(x, y).statistic
        observed_pearson = pearsonr(x, y).statistic

        null = np.empty(permutations)
        for i in range(permutations):
            null[i] = spearmanr(x, right_profile[:, RNG.permutation(len(genes))].ravel()).statistic
        p_perm = (1 + np.sum(null >= observed_spearman)) / (permutations + 1)

        bootstrap = np.empty(bootstraps)
        for i in range(bootstraps):
            a = resampled_type_profile(left, common_types, genes, RNG)
            b = resampled_type_profile(right, common_types, genes, RNG)
            bootstrap[i] = spearmanr(a.ravel(), b.ravel()).statistic
        ci_low, ci_high = np.nanpercentile(bootstrap, [2.5, 97.5])
        results.append(
            {
                "comparison": comparison,
                "program": program,
                "n_genes": len(genes),
                "n_common_types": len(common_types),
                "common_types": ";".join(common_types),
                "spearman_rho": observed_spearman,
                "pearson_r": observed_pearson,
                "bootstrap_ci_low": ci_low,
                "bootstrap_ci_high": ci_high,
                "ortholog_permutation_p": p_perm,
            }
        )
    return pd.DataFrame(results)


primary_concordance = profile_concordance(zf_pb, pop_pb, "zebrafish_vs_Popescu")

# GSE contains progenitor biases rather than mature classes. Collapse Popescu and
# zebrafish myeloid classes only for this explicitly secondary, exploratory comparison.
def collapse_for_gse(pb: pd.DataFrame) -> pd.DataFrame:
    result = pb.copy()
    result["broad_type"] = result["broad_type"].replace(
        {"macrophage": "myeloid", "neutrophil": "myeloid"}
    )
    numeric = ["n_cells", *human_genes, *PROGRAMS.keys()]
    result = (
        result.groupby(["dataset", "replicate", "age", "broad_type"], as_index=False)[numeric]
        .mean()
    )
    return result


secondary_concordance = profile_concordance(
    collapse_for_gse(zf_pb),
    gse_pb[gse_pb.broad_type.isin(["stem_progenitor", "erythroid", "myeloid"])],
    "zebrafish_vs_GSE189161_exploratory",
    bootstraps=300,
)
concordance = pd.concat([primary_concordance, secondary_concordance], ignore_index=True)
concordance["qvalue"] = multipletests(concordance.ortholog_permutation_p, method="fdr_bh")[1]
concordance.to_csv(RESULTS_DIR / "cross_species_lineage_profile_concordance.csv", index=False)
display(concordance)

# %% [markdown]
# ## 10. Age-adjusted program tests at the biological-replicate level
#
# Ages are scaled to $[0,1]$ separately in each cohort; coefficients therefore encode
# direction and standardized developmental-window effect, not equivalent absolute
# time between fish and humans. The GSE test below is limited to fetal-liver donors.

# %%
def fit_temporal_programs(pb: pd.DataFrame) -> pd.DataFrame:
    frame = pb.copy()
    minimum, maximum = frame.age.min(), frame.age.max()
    frame["age_scaled"] = (frame.age - minimum) / (maximum - minimum)
    dummies = pd.get_dummies(frame["broad_type"], prefix="type", drop_first=True, dtype=float)
    design = pd.concat(
        [pd.Series(1.0, index=frame.index, name="intercept"), frame[["age_scaled"]], dummies],
        axis=1,
    )
    rows = []
    for program in PROGRAMS:
        model = sm.OLS(frame[program].astype(float), design.astype(float)).fit(
            cov_type="cluster", cov_kwds={"groups": frame["replicate"]}
        )
        ci = model.conf_int().loc["age_scaled"]
        rows.append(
            {
                "dataset": frame.dataset.iloc[0],
                "program": program,
                "n_pseudobulks": len(frame),
                "n_replicates": frame.replicate.nunique(),
                "age_min": minimum,
                "age_max": maximum,
                "beta_age_scaled": model.params["age_scaled"],
                "se_clustered": model.bse["age_scaled"],
                "ci_low": ci.iloc[0],
                "ci_high": ci.iloc[1],
                "pvalue": model.pvalues["age_scaled"],
            }
        )
    return pd.DataFrame(rows)


temporal_models = pd.concat(
    [fit_temporal_programs(frame) for frame in (zf_pb, pop_pb, gse_pb)],
    ignore_index=True,
)
temporal_models["qvalue_within_dataset"] = temporal_models.groupby("dataset")["pvalue"].transform(
    lambda p: multipletests(p, method="fdr_bh")[1]
)
temporal_models.to_csv(RESULTS_DIR / "replicate_level_temporal_program_models.csv", index=False)
display(temporal_models)

# %% [markdown]
# ## 11. Independent human HSC/MPP replication using gene-level temporal slopes
#
# Popescu HSC/MPP cells and the GSE189161 HSC/MPP clusters are analyzed separately.
# For each mapped human gene, its donor-pseudobulk temporal slope is estimated. We then
# ask whether the vector of gene slopes agrees between the two independent human fetal
# datasets. This is a stronger use of GSE189161 than pretending it is a broad mature
# immune atlas.

# %%
def donor_gene_means(
    meta: pd.DataFrame,
    expression: np.ndarray,
    genes: list[str],
    keep_mask: np.ndarray,
    dataset: str,
    min_cells: int = 20,
) -> pd.DataFrame:
    positions = np.flatnonzero(keep_mask)
    groups = meta.iloc[positions][["replicate", "age"]].copy()
    groups["position"] = positions
    rows = []
    for (replicate, age), block in groups.groupby(["replicate", "age"]):
        if len(block) < min_cells:
            continue
        idx = block.position.to_numpy(dtype=int)
        row = {"dataset": dataset, "replicate": str(replicate), "age": float(age), "n_cells": len(idx)}
        row.update({g: float(np.mean(expression[idx, j])) for j, g in enumerate(genes)})
        rows.append(row)
    return pd.DataFrame(rows)


pop_hspc_mask = pop_meta["cell.labels"].eq("HSC_MPP").to_numpy()
gse_hspc_mask = gse_fetal & gse_meta["cluster_name"].isin(["HSC", "MPP-1", "MPP-2", "MPP-3"]).to_numpy()
pop_hspc = donor_gene_means(pop_meta, pop_log, pop_found, pop_hspc_mask, "Popescu_HSC_MPP")
gse_hspc = donor_gene_means(gse_meta, gse_log, human_genes, gse_hspc_mask, "GSE189161_HSC_MPP")

def gene_slopes(frame: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    age = (frame.age - frame.age.min()) / (frame.age.max() - frame.age.min())
    design = sm.add_constant(age.rename("age_scaled"))
    rows = []
    for gene in genes:
        response = frame[gene].astype(float)
        response_std = response.std(ddof=0)
        if response_std > 1e-8:
            response = (response - response.mean()) / response_std
        model = sm.OLS(response, design).fit()
        rows.append(
            {
                "dataset": frame.dataset.iloc[0],
                "human_gene": gene,
                "beta_age_scaled": model.params["age_scaled"],
                "pvalue": model.pvalues["age_scaled"],
            }
        )
    out = pd.DataFrame(rows)
    out["qvalue"] = multipletests(out.pvalue, method="fdr_bh")[1]
    return out


pop_slopes = gene_slopes(pop_hspc, human_genes)
gse_slopes = gene_slopes(gse_hspc, human_genes)
hspc_slopes = pop_slopes.merge(gse_slopes, on="human_gene", suffixes=("_Popescu", "_GSE189161"))
hspc_slopes = hspc_slopes.merge(
    ortholog_map[["human_gene", "marker_class"]].drop_duplicates(), on="human_gene", how="left"
)
hspc_slopes.to_csv(RESULTS_DIR / "human_hspc_gene_temporal_slopes.csv", index=False)

hspc_replication_rows = []
for marker_class, block in hspc_slopes.groupby("marker_class"):
    rho, pvalue = spearmanr(block.beta_age_scaled_Popescu, block.beta_age_scaled_GSE189161)
    same_direction = (
        np.sign(block.beta_age_scaled_Popescu) == np.sign(block.beta_age_scaled_GSE189161)
    )
    direction_fraction = np.mean(same_direction)
    direction_pvalue = binomtest(
        int(same_direction.sum()), len(block), p=0.5, alternative="greater"
    ).pvalue
    hspc_replication_rows.append(
        {
            "marker_class": marker_class,
            "n_genes": len(block),
            "spearman_rho": rho,
            "pvalue": pvalue,
            "same_direction_fraction": direction_fraction,
            "direction_pvalue": direction_pvalue,
        }
    )
hspc_replication = pd.DataFrame(hspc_replication_rows)
hspc_replication["qvalue"] = multipletests(hspc_replication.pvalue, method="fdr_bh")[1]
hspc_replication["direction_qvalue"] = multipletests(
    hspc_replication.direction_pvalue, method="fdr_bh"
)[1]
hspc_replication.to_csv(RESULTS_DIR / "popescu_gse189161_hspc_replication.csv", index=False)
display(hspc_replication)

# %% [markdown]
# ## 12. Black-and-white figures

# %%
dataset_counts = pd.DataFrame(
    {
        "dataset": ["Zebrafish saved", "Popescu", "GSE fetal", "GSE all"],
        "cells": [len(zf_meta), len(pop_meta), int(gse_fetal.sum()), len(gse_meta)],
        "replicates": [zf_meta.replicate.nunique(), pop_meta.replicate.nunique(),
                       gse_meta.loc[gse_fetal, "replicate"].nunique(), gse_meta.replicate.nunique()],
    }
)
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].bar(dataset_counts.dataset, dataset_counts.cells, color="white", edgecolor="black", hatch="///")
axes[0].set_ylabel("Cells")
axes[0].tick_params(axis="x", rotation=35)
axes[1].bar(dataset_counts.dataset, dataset_counts.replicates, color="0.75", edgecolor="black")
axes[1].set_ylabel("Biological replicates")
axes[1].tick_params(axis="x", rotation=35)
fig.suptitle("Datasets used in cross-species validation")
fig.tight_layout()
fig.savefig(FIGURES_DIR / "01_dataset_audit.pdf", bbox_inches="tight")
plt.show()

# %%
def plot_program_pseudobulks(program: str):
    frames = []
    for label, frame in [("Zebrafish", zf_pb), ("Popescu", pop_pb)]:
        temp = frame[["broad_type", program]].copy()
        temp["species"] = label
        frames.append(temp)
    plot_data = pd.concat(frames, ignore_index=True)
    types = sorted(set(plot_data.broad_type))
    fig, axes = plt.subplots(1, 2, figsize=(13, 4), sharey=True)
    for ax, species in zip(axes, ["Zebrafish", "Popescu"]):
        block = plot_data[plot_data.species.eq(species)]
        values = [block.loc[block.broad_type.eq(t), program].to_numpy() for t in types]
        bp = ax.boxplot(values, tick_labels=types, showfliers=False, patch_artist=True)
        for patch in bp["boxes"]:
            patch.set(facecolor="white", edgecolor="black", hatch="//")
        for element in ("whiskers", "caps", "medians"):
            for artist in bp[element]:
                artist.set(color="black")
        ax.axhline(0, color="0.5", lw=0.8, ls="--")
        ax.set_title(species)
        ax.tick_params(axis="x", rotation=35)
        ax.set_ylabel(f"Within-species {program} score")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / f"02_{program}_pseudobulks.pdf", bbox_inches="tight")
    plt.show()


plot_program_pseudobulks("proliferation")
plot_program_pseudobulks("apoptosis")

# %%
common_types = sorted(set(zf_pb.broad_type).intersection(pop_pb.broad_type))
genes = PROGRAMS["proliferation"]
x_profile = standardized_type_profile(zf_pb, common_types, genes)
y_profile = standardized_type_profile(pop_pb, common_types, genes)
markers = ["o", "s", "^", "D", "x", "+"]
fig, ax = plt.subplots(figsize=(6, 6))
for index, broad_type in enumerate(common_types):
    marker = markers[index % len(markers)]
    style = {"color": "black"} if marker in {"x", "+"} else {
        "facecolors": "none", "edgecolors": "black"
    }
    ax.scatter(x_profile[index], y_profile[index], s=28, marker=marker,
               label=broad_type, alpha=0.8, **style)
limits = [min(x_profile.min(), y_profile.min()), max(x_profile.max(), y_profile.max())]
ax.plot(limits, limits, color="black", ls="--", lw=0.8)
ax.set(xlabel="Zebrafish standardized lineage profile", ylabel="Popescu standardized lineage profile")
ax.legend(frameon=False, fontsize=8)
ax.set_title("Ortholog-mapped proliferation profile")
fig.tight_layout()
fig.savefig(FIGURES_DIR / "03_zebrafish_popescu_proliferation_profile.pdf", bbox_inches="tight")
plt.show()

# %%
fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
for ax, frame, title in zip(
    axes,
    [zf_pb, pop_pb, gse_pb],
    ["Zebrafish", "Popescu fetal liver", "GSE189161 fetal HSPCs"],
):
    donor_means = frame.groupby(["replicate", "age"], as_index=False).proliferation.mean()
    ax.scatter(donor_means.age, donor_means.proliferation, facecolors="white", edgecolors="black")
    if len(donor_means) > 1:
        slope, intercept = np.polyfit(donor_means.age, donor_means.proliferation, 1)
        grid = np.linspace(donor_means.age.min(), donor_means.age.max(), 100)
        ax.plot(grid, intercept + slope * grid, color="black")
    ax.set_title(title)
    ax.set_xlabel("hpf" if title == "Zebrafish" else "PCW")
axes[0].set_ylabel("Mean proliferation score")
fig.tight_layout()
fig.savefig(FIGURES_DIR / "04_replicate_proliferation_time.pdf", bbox_inches="tight")
plt.show()

# %%
fig, ax = plt.subplots(figsize=(6, 6))
for marker_class, block in hspc_slopes.groupby("marker_class"):
    marker = "o" if marker_class == "proliferation" else "x"
    style = {"color": "black"} if marker == "x" else {
        "facecolors": "none", "edgecolors": "black"
    }
    ax.scatter(
        block.beta_age_scaled_Popescu,
        block.beta_age_scaled_GSE189161,
        marker=marker,
        label=marker_class,
        **style,
    )
ax.axhline(0, color="0.6", lw=0.8)
ax.axvline(0, color="0.6", lw=0.8)
ax.set(
    xlabel="Popescu HSC/MPP standardized temporal slope",
    ylabel="GSE189161 HSC/MPP standardized temporal slope",
    title="Independent human HSPC replication",
)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(FIGURES_DIR / "05_popescu_gse_hspc_gene_slopes.pdf", bbox_inches="tight")
plt.show()

# %% [markdown]
# ## 13. Evidence scorecard and interpretation
#
# “Strong” below means strong **transcriptomic conservation within these datasets**,
# never causal or transport validation. A causal claim would additionally need human
# lineage tracing, perturbation, or an independently measured proliferation/death rate.

# %%
def safe_bool(value) -> bool:
    return bool(value) if pd.notna(value) else False


primary_pro = concordance[
    concordance.comparison.eq("zebrafish_vs_Popescu") & concordance.program.eq("proliferation")
].iloc[0]
primary_apo = concordance[
    concordance.comparison.eq("zebrafish_vs_Popescu") & concordance.program.eq("apoptosis")
].iloc[0]
secondary_pro = concordance[
    concordance.comparison.eq("zebrafish_vs_GSE189161_exploratory")
    & concordance.program.eq("proliferation")
].iloc[0]
secondary_apo = concordance[
    concordance.comparison.eq("zebrafish_vs_GSE189161_exploratory")
    & concordance.program.eq("apoptosis")
].iloc[0]
human_pro = hspc_replication[hspc_replication.marker_class.eq("proliferation")].iloc[0]
human_apo = hspc_replication[hspc_replication.marker_class.eq("apoptosis")].iloc[0]

scorecard = pd.DataFrame(
    [
        {
            "criterion": "Reviewed markers detected in all datasets",
            "pass": coverage.fraction.min() >= 0.90,
            "metric": f"minimum coverage={coverage.fraction.min():.1%}",
        },
        {
            "criterion": "Primary Popescu proliferation program profile",
            "pass": safe_bool(primary_pro.spearman_rho > 0 and primary_pro.bootstrap_ci_low > 0),
            "metric": (
                f"rho={primary_pro.spearman_rho:.3f}, bootstrap CI="
                f"[{primary_pro.bootstrap_ci_low:.3f}, {primary_pro.bootstrap_ci_high:.3f}]"
            ),
        },
        {
            "criterion": "Primary Popescu exact proliferation ortholog specificity",
            "pass": safe_bool(primary_pro.spearman_rho > 0 and primary_pro.ortholog_permutation_p < 0.05),
            "metric": f"ortholog-label permutation p={primary_pro.ortholog_permutation_p:.4g}",
        },
        {
            "criterion": "Primary Popescu apoptosis program profile",
            "pass": safe_bool(primary_apo.spearman_rho > 0 and primary_apo.bootstrap_ci_low > 0),
            "metric": (
                f"rho={primary_apo.spearman_rho:.3f}, bootstrap CI="
                f"[{primary_apo.bootstrap_ci_low:.3f}, {primary_apo.bootstrap_ci_high:.3f}]"
            ),
        },
        {
            "criterion": "Exploratory GSE189161 proliferation progenitor profile",
            "pass": safe_bool(secondary_pro.spearman_rho > 0 and secondary_pro.qvalue < 0.05),
            "metric": f"rho={secondary_pro.spearman_rho:.3f}, permutation FDR={secondary_pro.qvalue:.4g}",
        },
        {
            "criterion": "Independent human HSPC proliferation slope directions",
            "pass": safe_bool(human_pro.same_direction_fraction > 0.5 and human_pro.direction_qvalue < 0.05),
            "metric": (
                f"same direction={human_pro.same_direction_fraction:.1%}, "
                f"binomial FDR={human_pro.direction_qvalue:.4g}"
            ),
        },
        {
            "criterion": "Independent human HSPC proliferation slope magnitudes",
            "pass": safe_bool(human_pro.spearman_rho > 0 and human_pro.qvalue < 0.05),
            "metric": f"rho={human_pro.spearman_rho:.3f}, rank-correlation FDR={human_pro.qvalue:.4g}",
        },
        {
            "criterion": "Independent human HSPC apoptosis slope directions",
            "pass": safe_bool(human_apo.same_direction_fraction > 0.5 and human_apo.direction_qvalue < 0.05),
            "metric": (
                f"same direction={human_apo.same_direction_fraction:.1%}, "
                f"binomial FDR={human_apo.direction_qvalue:.4g}"
            ),
        },
    ]
)

pro_primary = safe_bool(primary_pro.spearman_rho > 0 and primary_pro.bootstrap_ci_low > 0)
pro_exact = safe_bool(primary_pro.ortholog_permutation_p < 0.05)
pro_secondary = safe_bool(
    (secondary_pro.spearman_rho > 0 and secondary_pro.qvalue < 0.05)
    or (human_pro.same_direction_fraction > 0.5 and human_pro.direction_qvalue < 0.05)
)
pro_magnitude = safe_bool(human_pro.spearman_rho > 0 and human_pro.qvalue < 0.05)
if pro_primary and pro_exact and pro_secondary and pro_magnitude:
    proliferation_grade = "strong transcriptomic support"
elif pro_primary and pro_secondary:
    proliferation_grade = "moderate program-level transcriptomic support"
else:
    proliferation_grade = "limited transcriptomic support"

apo_primary = safe_bool(primary_apo.spearman_rho > 0 and primary_apo.bootstrap_ci_low > 0)
apo_secondary = safe_bool(
    (secondary_apo.spearman_rho > 0 and secondary_apo.qvalue < 0.05)
    or (human_apo.same_direction_fraction > 0.5 and human_apo.direction_qvalue < 0.05)
)
apoptosis_grade = (
    "moderate program-level transcriptomic support"
    if apo_primary and apo_secondary
    else "limited transcriptomic support"
)
evidence_grade = (
    f"Proliferation: {proliferation_grade}; apoptosis: {apoptosis_grade}; "
    "not causal or dynamic validation"
)

scorecard["evidence_grade"] = evidence_grade
scorecard.to_csv(RESULTS_DIR / "validation_scorecard.csv", index=False)
display(scorecard)
display(Markdown(f"### Overall grade\n\n**{evidence_grade}**"))

# %%
report = f"""# Zebrafish-to-human cross-species validation report

## Scope

This analysis reused the saved zebrafish checkpoint and did not rerun GraphVelo,
moscot, or the earlier transport workflow. The full zebrafish count H5AD was used only
to recover marker genes absent from the checkpoint's 2,986-gene selected feature space.

## Cohorts

- Zebrafish: {len(zf_meta):,} saved cells, {zf_meta.replicate.nunique()} fish.
- Popescu E-MTAB-7407: {len(pop_meta):,} annotated fetal-liver cells,
  {pop_meta.replicate.nunique()} donors, {pop_meta.age.min():.2f}-{pop_meta.age.max():.2f} PCW.
- GSE189161: {len(gse_meta):,} CD34-enriched cells from {gse_meta.replicate.nunique()} samples;
  the replication analysis used {gse_meta.loc[gse_fetal, 'replicate'].nunique()} fetal-liver donors.

## Overall evidence grade

**{evidence_grade}**

## Primary zebrafish-Popescu lineage result

- Proliferation: Spearman rho={primary_pro.spearman_rho:.3f},
  95% donor/fish bootstrap CI [{primary_pro.bootstrap_ci_low:.3f}, {primary_pro.bootstrap_ci_high:.3f}],
  ortholog-label permutation p={primary_pro.ortholog_permutation_p:.4g}.
- Apoptosis: Spearman rho={primary_apo.spearman_rho:.3f},
  95% donor/fish bootstrap CI [{primary_apo.bootstrap_ci_low:.3f}, {primary_apo.bootstrap_ci_high:.3f}],
  ortholog-label permutation p={primary_apo.ortholog_permutation_p:.4g}.

## Independent human HSC/MPP replication

- Proliferation gene temporal slopes: rho={human_pro.spearman_rho:.3f},
  rank-correlation FDR={human_pro.qvalue:.4g}, same-direction
  fraction={human_pro.same_direction_fraction:.1%}, direction FDR={human_pro.direction_qvalue:.4g}.
- Apoptosis gene temporal slopes: rho={human_apo.spearman_rho:.3f},
  rank-correlation FDR={human_apo.qvalue:.4g}, same-direction
  fraction={human_apo.same_direction_fraction:.1%}, direction FDR={human_apo.direction_qvalue:.4g}.

Slope-direction agreement and slope-magnitude rank correlation answer different
questions. Significant direction agreement means the two human cohorts tend to move
the same way with fetal age; a nonsignificant rank correlation means they do not
reliably reproduce the ordering of gene-specific effect sizes.

## Interpretation boundary

Agreement supports conservation of the mapped expression program. It does not validate
an absolute growth rate, a cell-to-cell transport coupling, or an RNA-velocity vector in
humans. Strong causal validation would require independent human lineage tracing,
perturbation, or measured proliferation/death outcomes.
"""
(RESULTS_DIR / "validation_report.md").write_text(report)
display(Markdown(report))

print("Saved results:", RESULTS_DIR)
print("Saved figures:", FIGURES_DIR)
