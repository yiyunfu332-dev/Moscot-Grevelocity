"""Build a frozen, provenance-rich zebrafish growth-marker mapping table.

This script does not select genes from the validation outcomes.  It combines an
existing core panel with external cell-cycle/pathway resources, then requires a
reviewed ZFIN human-zebrafish orthology record and presence in all_genes.txt.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
SOURCE_DIR = ROOT / "growth_marker_reference"
CORE_CSV = ROOT / "reviewed_zebrafish_to_human_growth_orthologs.csv"
ZFIN_TSV = ROOT / "human_orthos.txt"
ALL_GENES = ROOT / "all_genes.txt"
OUTPUT_CSV = ROOT / "reviewed_zebrafish_to_human_growth_orthologs_expanded.csv"
HUMAN_PANEL_CSV = SOURCE_DIR / "expanded_growth_marker_human_panel.csv"
MANIFEST_CSV = SOURCE_DIR / "source_manifest.csv"

SEURAT_URL = "https://satijalab.org/seurat/reference/cc.genes.updated.2019.html"
TIROSH_DOI = "https://doi.org/10.1126/science.aad0501"
KEGG_CELL_CYCLE_URL = "https://www.kegg.jp/pathway/hsa04110"
KEGG_APOPTOSIS_URL = "https://www.kegg.jp/pathway/hsa04210"
REACTOME_MITOTIC_URL = "https://reactome.org/content/detail/R-HSA-69278"
REACTOME_DNA_REPLICATION_URL = "https://reactome.org/content/detail/R-HSA-69306"
REACTOME_APOPTOSIS_URL = "https://reactome.org/content/detail/R-HSA-109581"
REACTOME_INTRINSIC_URL = "https://reactome.org/content/detail/R-HSA-109606"
ZFIN_URL = "https://zfin.org/downloads/human_orthos.txt"

# Core genes retained from the previous reviewed panel even when absent from the
# phase-marker resource.  They are assigned to a mechanistically coherent phase.
CORE_S_PHASE = {
    "PCNA", "CDC6", "ORC1", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6",
    "CCNE1", "MCM10", "TYMS", "GINS1", "GINS2", "GINS3", "POLA2",
    "RFC2", "RFC3", "GMNN",
}
CORE_G2M_PHASE = {
    "MKI67", "TOP2A", "UBE2C", "CDK1", "CDC20", "CCNB1", "CCNB2",
    "PLK1", "AURKB", "AURKA", "BUB1", "BUB3", "MAD2L1", "CENPF",
    "CENPK", "SMC4", "KIF11", "TPX2", "NDC80", "BIRC5",
}

# Directionally curated pro-death genes.  The full KEGG/Reactome apoptosis sets
# contain anti-apoptotic and context-dependent survival genes, so the entire pathway
# must not be passed to moscot as a same-direction death score.
PRO_APOPTOTIC_SUBPROGRAMS = {
    "intrinsic_apoptosis": {
        "BAX", "BAD", "BID", "BBC3", "BCL2L11", "APAF1", "DIABLO",
        "AIFM1", "ENDOG", "CASP9", "BNIP3L", "OMA1",
    },
    "extrinsic_apoptosis": {
        "FAS", "FASLG", "FADD", "TNFSF10", "TNFRSF1A", "TRADD",
        "CASP8", "DAB2IP",
    },
    "execution_apoptosis": {
        "CASP2", "CASP3", "CASP7", "DFFA", "DFFB", "ACIN1", "ROCK1",
    },
    "stress_pro_apoptotic": {
        "TP53", "DAP3", "DAPK1", "DAPK2", "DAPK3", "DDIT3", "EIF2AK3",
        "TP53BP2", "UNC5B", "MAP3K5", "PRKCD",
    },
}


def read_gmt(path: Path) -> dict[str, set[str]]:
    result = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            result[fields[0]] = {gene.upper() for gene in fields[2:] if gene}
    return result


def compact_join(values) -> str:
    return ";".join(sorted({str(value) for value in values if pd.notna(value) and str(value)}))


core = pd.read_csv(CORE_CSV)
core["human_gene"] = core["human_gene"].astype(str).str.upper()
core_human = {
    marker_class: set(block.human_gene)
    for marker_class, block in core.groupby("marker_class")
}

seurat = pd.read_csv(SOURCE_DIR / "seurat_cc_genes_updated_2019.tsv", sep="\t")
seurat["human_gene"] = seurat["human_gene"].astype(str).str.upper()
seurat_phase = dict(zip(seurat.human_gene, seurat.subprogram))

zfin = pd.read_csv(ZFIN_TSV, sep="\t", header=None, dtype=str)
zfin.columns = [
    "zfin_id", "zebrafish_gene", "zebrafish_name", "human_gene", "human_name",
    "omim_id", "human_entrez", "zebrafish_entrez", "evidence_code",
    "zfin_publication", "evidence_method", "eco_code", "evidence_description",
]
zfin["zebrafish_gene"] = zfin["zebrafish_gene"].astype(str).str.casefold()
zfin["human_gene"] = zfin["human_gene"].astype(str).str.upper()

all_genes_lookup = {
    line.strip().casefold(): line.strip()
    for line in ALL_GENES.read_text().splitlines()
    if line.strip()
}

reactome = read_gmt(SOURCE_DIR / "reactome_selected_pathways.gmt")
reactome_cell_cycle = reactome["Cell Cycle, Mitotic"] | reactome["DNA Replication"]
reactome_apoptosis = reactome["Apoptosis"] | reactome["Intrinsic Pathway for Apoptosis"]

kegg_cell_entrez = {
    line.rstrip().split("\t")[1].split(":")[1]
    for line in (SOURCE_DIR / "kegg_hsa04110_cell_cycle.tsv").read_text().splitlines()
}
kegg_apoptosis_entrez = {
    line.rstrip().split("\t")[1].split(":")[1]
    for line in (SOURCE_DIR / "kegg_hsa04210_apoptosis.tsv").read_text().splitlines()
}
kegg_cell_cycle = set(zfin.loc[zfin.human_entrez.isin(kegg_cell_entrez), "human_gene"])
kegg_apoptosis = set(zfin.loc[zfin.human_entrez.isin(kegg_apoptosis_entrez), "human_gene"])

proliferation_human = set(seurat.human_gene) | core_human["proliferation"]
apoptosis_subclass = {
    gene: subclass
    for subclass, genes in PRO_APOPTOTIC_SUBPROGRAMS.items()
    for gene in genes
}
apoptosis_human = set(apoptosis_subclass) | core_human["apoptosis"]


def sources_for(gene: str, marker_class: str) -> list[str]:
    sources = []
    if marker_class == "proliferation":
        if gene in seurat_phase:
            sources.append("Seurat_cc.genes.updated.2019/Tirosh2016")
        if gene in kegg_cell_cycle:
            sources.append("KEGG_hsa04110_Cell_cycle")
        if gene in reactome_cell_cycle:
            sources.append("Reactome_Cell_Cycle_Mitotic_or_DNA_Replication")
    else:
        if gene in kegg_apoptosis:
            sources.append("KEGG_hsa04210_Apoptosis")
        if gene in reactome_apoptosis:
            sources.append("Reactome_Apoptosis_or_Intrinsic_Apoptosis")
    if gene in core_human.get(marker_class, set()):
        sources.append("previous_reviewed_core_panel")
    return sources


def subclass_for(gene: str, marker_class: str) -> str:
    if marker_class == "apoptosis":
        return apoptosis_subclass.get(gene, "previous_core_pro_apoptotic")
    if gene in seurat_phase:
        return seurat_phase[gene]
    if gene in CORE_S_PHASE:
        return "S_phase"
    if gene in CORE_G2M_PHASE:
        return "G2M_phase"
    raise ValueError(f"No directional cell-cycle subclass for {gene}")


panel_rows = []
mapping_rows = []
for marker_class, genes in [
    ("proliferation", proliferation_human),
    ("apoptosis", apoptosis_human),
]:
    for human_gene in sorted(genes):
        sources = sources_for(human_gene, marker_class)
        if marker_class == "apoptosis" and not sources:
            # Directional seed genes still require pathway membership, unless they
            # were already part of the reviewed core panel.
            continue
        subclass = subclass_for(human_gene, marker_class)
        block = zfin[
            zfin.human_gene.eq(human_gene)
            & zfin.zebrafish_gene.isin(all_genes_lookup)
        ].copy()
        if block.empty:
            continue
        unique_pairs = block[["zebrafish_gene", "human_gene"]].drop_duplicates()
        panel_rows.append(
            {
                "human_gene": human_gene,
                "marker_class": marker_class,
                "marker_subclass": subclass,
                "pathway_sources": ";".join(sources),
                "n_external_sources": sum("previous_reviewed" not in value for value in sources),
                "n_zebrafish_orthologs_in_all_genes": len(unique_pairs),
            }
        )
        for _, pair in unique_pairs.iterrows():
            evidence = block[
                block.zebrafish_gene.eq(pair.zebrafish_gene)
                & block.human_gene.eq(pair.human_gene)
            ]
            urls = []
            if human_gene in seurat_phase:
                urls.extend([SEURAT_URL, TIROSH_DOI])
            if human_gene in kegg_cell_cycle:
                urls.append(KEGG_CELL_CYCLE_URL)
            if human_gene in kegg_apoptosis:
                urls.append(KEGG_APOPTOSIS_URL)
            if human_gene in reactome_cell_cycle:
                urls.extend([REACTOME_MITOTIC_URL, REACTOME_DNA_REPLICATION_URL])
            if human_gene in reactome_apoptosis:
                urls.extend([REACTOME_APOPTOSIS_URL, REACTOME_INTRINSIC_URL])
            urls.append(ZFIN_URL)
            mapping_rows.append(
                {
                    "zebrafish_gene": all_genes_lookup[pair.zebrafish_gene],
                    "human_gene": human_gene,
                    "source_gene": human_gene,
                    "marker_class": marker_class,
                    "marker_subclass": subclass,
                    "mapping_direction": "zebrafish_to_human",
                    "mapping_database": "ZFIN human orthology download",
                    "evidence_status": "reviewed_external_consensus",
                    "pathway_sources": ";".join(sources),
                    "n_external_sources": sum("previous_reviewed" not in value for value in sources),
                    "selection_rule": (
                        "core_or_Seurat_phase_marker+ZFIN_ortholog+present_in_all_genes"
                        if marker_class == "proliferation"
                        else "directionally_curated_pro_death+KEGG_or_Reactome_or_core+ZFIN_ortholog+present_in_all_genes"
                    ),
                    "zfin_ids": compact_join(evidence.zfin_id),
                    "zfin_evidence_codes": compact_join(evidence.evidence_code),
                    "zfin_publications": compact_join(evidence.zfin_publication),
                    "source_urls": ";".join(sorted(set(urls))),
                }
            )

panel = pd.DataFrame(panel_rows).sort_values(["marker_class", "marker_subclass", "human_gene"])
mapping = pd.DataFrame(mapping_rows).sort_values(
    ["marker_class", "marker_subclass", "human_gene", "zebrafish_gene"]
)

assert not mapping.duplicated(["zebrafish_gene", "human_gene", "marker_class"]).any()
assert mapping.mapping_direction.eq("zebrafish_to_human").all()
assert set(core.zebrafish_gene.str.casefold()).issubset(set(mapping.zebrafish_gene.str.casefold()))
assert set(mapping.zebrafish_gene.str.casefold()).issubset(all_genes_lookup)
assert mapping.groupby("marker_class").human_gene.nunique()["proliferation"] >= 75
assert mapping.groupby("marker_class").human_gene.nunique()["apoptosis"] >= 25

mapping.to_csv(OUTPUT_CSV, index=False)
panel.to_csv(HUMAN_PANEL_CSV, index=False)

manifest = pd.DataFrame(
    [
        ["Seurat cc.genes.updated.2019", "S/G2M phase markers", SEURAT_URL, TIROSH_DOI, "2019 list"],
        ["KEGG hsa04110", "Cell cycle pathway membership", KEGG_CELL_CYCLE_URL, "https://rest.kegg.jp/link/hsa/hsa04110", "retrieved 2026-09-08"],
        ["KEGG hsa04210", "Apoptosis pathway membership", KEGG_APOPTOSIS_URL, "https://rest.kegg.jp/link/hsa/hsa04210", "retrieved 2026-09-08"],
        ["Reactome R-HSA-69278/R-HSA-69306", "Mitotic cell cycle and DNA replication", REACTOME_MITOTIC_URL, REACTOME_DNA_REPLICATION_URL, "release downloaded 2026-09-08"],
        ["Reactome R-HSA-109581/R-HSA-109606", "Apoptosis and intrinsic apoptosis", REACTOME_APOPTOSIS_URL, REACTOME_INTRINSIC_URL, "release downloaded 2026-09-08"],
        ["ZFIN human orthology", "Reviewed human-zebrafish orthology", ZFIN_URL, str(ZFIN_TSV), "local snapshot"],
    ],
    columns=["source", "use", "url", "secondary_url_or_local_file", "version_or_retrieval"],
)
manifest.to_csv(MANIFEST_CSV, index=False)

summary = mapping.groupby("marker_class").agg(
    mapping_rows=("zebrafish_gene", "size"),
    zebrafish_genes=("zebrafish_gene", "nunique"),
    human_genes=("human_gene", "nunique"),
    subclasses=("marker_subclass", "nunique"),
)
print(summary.to_string())
print(f"\nWrote {OUTPUT_CSV}")
print(f"Wrote {HUMAN_PANEL_CSV}")
print(f"Wrote {MANIFEST_CSV}")
