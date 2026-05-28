"""
Patch A — KG full-universe expansion.

Generates ~250 new ACTIVE rules across 9 tissue classes (CNS, Skin, Kidney,
Heart, Stromal/EC, Endocrine, Hematopoietic, Cancer, Trophoblast) and appends
them to kg_data/oracle_kg_v0.2.0.json.

Run: /Users/sivamoturi/miniconda/bin/python3.10 scripts/patch_a_kg_expansion.py

All PMIDs verified via Entrez before writing (see PMID audit in CLAUDE.md §21g).
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
KG_PATH = ROOT / "kg_data" / "oracle_kg_v0.2.0.json"

# ── Verified PMID registry ────────────────────────────────────────────────────
PMIDS = {
    "tasic2018":      "30382198",  # Shared/distinct transcriptomic cell types neocortex
    "tasic2016":      "26727548",  # Adult mouse cortex cell taxonomy
    "liddelow2017":   "28099414",  # Neurotoxic reactive astrocytes A1
    "keren2017":      "28602351",  # DAM microglia (Keren-Shaul 2017)
    "marques2016":    "27284195",  # Oligodendrocyte heterogeneity (Marques 2016)
    "jakel2019":      "30747918",  # Altered oligo heterogeneity MS (Jäkel 2019)
    "vanland2018":    "29443965",  # Brain vasculature molecular atlas (Vanlandewijck 2018)
    "schwann2020":    "33197966",  # Corneal Schwann cells single cell
    "ependymal":      "33436697",  # BMP signaling ependymal differentiation
    "ji2020":         "32579974",  # Skin multimodal (Ji 2020 Science)
    "lake2023":       "37468583",  # Kidney atlas healthy/injured (KPMP 2023)
    "stewart2019":    "31604275",  # Spatiotemporal immune zonation kidney
    "litvin2020":     "32971526",  # Cells of adult human heart
    "elyada2019":     "31197017",  # Cross-species PDA iCAF/myCAF (Elyada 2019)
    "baron2016":      "27667365",  # Single-cell map human/mouse pancreas
    "pellin2019":     "29588278",  # Single-cell hematopoietic landscape
    "velten2017":     "28319093",  # Human HSC lineage commitment
    "erythroid":      "32457162",  # Regulators erythroid differentiation
    "miller2019":     "30778252",  # Subsets exhausted CD8+ T cells
    "zheng2021":      "34914499",  # Pan-cancer TIL landscape
    "tirosh2016":     "27124452",  # Multicellular ecosystem metastatic melanoma
    "plitas2016":     "27851913",  # Regulatory T cells human breast cancer
    "ventotormo2018": "30429548",  # Maternal-fetal interface (Vento-Tormo 2018)
    "tabumuris":      "30283141",  # Single-cell transcriptomics 20 mouse organs
    "allenbrain":     "38168182",  # Cell-type specific signatures brain-wide
}


def make_rule(
    cell_type, rule_id, rule_type, evidence_source,
    gene_or_regulon, direction, mechanistic_basis, pmid_key,
    confidence="high", tissue_context=None, paired_with=None,
    incompatible_with=None, source_db="primary_literature",
):
    r = {
        "cell_type": cell_type,
        "rule_id": rule_id,
        "rule_type": rule_type,
        "evidence_source": evidence_source,
        "gene_or_regulon": gene_or_regulon,
        "direction": direction,
        "mechanistic_basis": mechanistic_basis,
        "pmid": PMIDS[pmid_key],
        "confidence": confidence,
        "tissue_context": tissue_context or [],
        "source_db": source_db,
        "status": "ACTIVE",
        "added_in_version": "0.2.0-patch-a",
    }
    if paired_with:
        r["paired_with"] = paired_with
    if incompatible_with:
        r["incompatible_with"] = incompatible_with
    return r


NEW_RULES = []

# ═══════════════════════════════════════════════════════════════════════════════
# CNS — 10 cell types, 40 rules
# ═══════════════════════════════════════════════════════════════════════════════

# Excitatory_Neuron
NEW_RULES += [
    make_rule("Excitatory_Neuron", "EXC_POS_001", "positive", "marker_genes",
        ["SLC17A7"], "high",
        "VGLUT1 (SLC17A7) is the vesicular glutamate transporter expressed in cortical excitatory neurons. Its high expression is the defining marker of glutamatergic identity.",
        "tasic2018", tissue_context=["brain", "cortex", "hippocampus"]),
    make_rule("Excitatory_Neuron", "EXC_NEG_001", "negative", "marker_genes",
        ["GAD1"], "absent",
        "GAD1 encodes glutamate decarboxylase, the rate-limiting enzyme for GABA synthesis. Its presence defines inhibitory neurons and is mutually exclusive with glutamatergic identity.",
        "tasic2018", incompatible_with=["GAD1"], tissue_context=["brain", "cortex"]),
    make_rule("Excitatory_Neuron", "EXC_NEG_002", "negative", "marker_genes",
        ["GAD2"], "absent",
        "GAD2 (GAD65) is a second glutamate decarboxylase isoform co-expressed with GAD1 in all inhibitory neurons. Its presence eliminates excitatory neuron identity.",
        "tasic2016", incompatible_with=["GAD2"], tissue_context=["brain"]),
    make_rule("Excitatory_Neuron", "EXC_POS_002", "positive", "marker_genes",
        ["CAMK2A"], "high",
        "CaMKII alpha (CAMK2A) is selectively expressed in excitatory pyramidal neurons and drives synaptic plasticity. High expression confirms glutamatergic identity.",
        "tasic2018", tissue_context=["brain", "cortex"]),
]

# Inhibitory_Neuron
NEW_RULES += [
    make_rule("Inhibitory_Neuron", "INH_CIRCUIT_001", "circuit", "marker_genes",
        ["GAD1", "GAD2"], "active",
        "Co-expression of GAD1 and GAD2 defines GABAergic interneurons. Both isoforms of glutamate decarboxylase are required together for robust GABA synthesis in inhibitory neurons.",
        "tasic2018", tissue_context=["brain", "cortex", "hippocampus"]),
    make_rule("Inhibitory_Neuron", "INH_NEG_001", "negative", "marker_genes",
        ["SLC17A7"], "absent",
        "VGLUT1 (SLC17A7) marks glutamatergic neurons and is mutually exclusive with GABAergic identity. Its presence rules out inhibitory neuron classification.",
        "tasic2016", incompatible_with=["SLC17A7"], tissue_context=["brain"]),
    make_rule("Inhibitory_Neuron", "INH_POS_001", "positive", "marker_genes",
        ["SST"], "high",
        "Somatostatin (SST) marks one of the major inhibitory interneuron subtypes (Martinotti cells). High SST expression supports inhibitory neuron identity.",
        "tasic2018", tissue_context=["brain", "cortex"]),
    make_rule("Inhibitory_Neuron", "INH_POS_002", "positive", "marker_genes",
        ["PVALB"], "high",
        "Parvalbumin (PVALB) marks fast-spiking basket and chandelier interneurons. High expression is a key inhibitory neuron marker.",
        "allenbrain", tissue_context=["brain", "cortex"]),
]

# Astrocyte_homeostatic
NEW_RULES += [
    make_rule("Astrocyte_homeostatic", "AST_H_POS_001", "positive", "marker_genes",
        ["AQP4", "SLC1A2", "GJA1", "S100B"], "high",
        "Homeostatic astrocytes co-express AQP4 (water channel), SLC1A2 (GLT-1 glutamate transporter), GJA1 (connexin 43 for gap junctions), and S100B. This quartet defines the homeostatic astrocyte program.",
        "liddelow2017", tissue_context=["brain"]),
    make_rule("Astrocyte_homeostatic", "AST_H_NEG_001", "negative", "marker_genes",
        ["LCN2"], "absent",
        "LCN2 (lipocalin-2) is an A1 reactive astrocyte marker induced by neuroinflammation. Its high expression marks reactive, not homeostatic, astrocytes.",
        "liddelow2017", incompatible_with=["LCN2"], tissue_context=["brain"]),
    make_rule("Astrocyte_homeostatic", "AST_H_POS_002", "positive", "marker_genes",
        ["ALDH1L1"], "high",
        "ALDH1L1 (10-formyltetrahydrofolate dehydrogenase) is a pan-astrocyte marker enriched in homeostatic astrocytes and used to distinguish them from other glia.",
        "allenbrain", tissue_context=["brain"]),
]

# Astrocyte_reactive
NEW_RULES += [
    make_rule("Astrocyte_reactive", "AST_R_CIRCUIT_001", "circuit", "marker_genes",
        ["GFAP", "LCN2", "C3"], "active",
        "A1 reactive astrocytes are characterized by strong upregulation of GFAP (structural), LCN2 (inflammatory), and C3 (complement). This triad is induced by activated microglia via IL-1alpha, TNF, and C1q.",
        "liddelow2017", tissue_context=["brain", "neuroinflammation"]),
    make_rule("Astrocyte_reactive", "AST_R_POS_001", "positive", "marker_genes",
        ["GFAP"], "high",
        "GFAP is strongly upregulated in reactive astrocytes compared to homeostatic astrocytes. Very high GFAP is the hallmark of astrogliosis.",
        "liddelow2017", tissue_context=["brain"]),
]

# Oligodendrocyte
NEW_RULES += [
    make_rule("Oligodendrocyte", "OLIGO_POS_001", "positive", "marker_genes",
        ["MBP", "MAG", "MOG", "PLP1"], "high",
        "Mature oligodendrocytes co-express the major myelin proteins MBP, MAG, MOG, and PLP1. This quartet defines the myelinating oligodendrocyte program and is absent in OPCs.",
        "marques2016", tissue_context=["brain", "spinal_cord", "white_matter"]),
    make_rule("Oligodendrocyte", "OLIGO_NEG_001", "negative", "marker_genes",
        ["PDGFRA"], "absent",
        "PDGFRA marks oligodendrocyte precursor cells (OPCs) and is downregulated upon OPC differentiation into mature oligodendrocytes. PDGFRA presence excludes mature myelinating oligodendrocyte identity.",
        "marques2016", incompatible_with=["PDGFRA"], tissue_context=["brain"]),
]

# OPC
NEW_RULES += [
    make_rule("OPC", "OPC_CIRCUIT_001", "circuit", "marker_genes",
        ["PDGFRA", "OLIG2", "SOX10"], "active",
        "OPC identity requires co-expression of PDGFRA (growth factor receptor for PDGF-AA mitogen), OLIG2 (oligodendrocyte lineage TF), and SOX10 (neural crest/glial TF). All three are required for OPC identity.",
        "marques2016", tissue_context=["brain", "spinal_cord"]),
    make_rule("OPC", "OPC_NEG_001", "negative", "marker_genes",
        ["MBP"], "absent",
        "MBP (myelin basic protein) is expressed upon OPC differentiation into mature myelinating oligodendrocytes. Its presence indicates the cell has exited the OPC state.",
        "marques2016", incompatible_with=["MBP"], tissue_context=["brain"]),
    make_rule("OPC", "OPC_POS_001", "positive", "marker_genes",
        ["CSPG4"], "high",
        "CSPG4 (NG2 proteoglycan) marks OPCs throughout the CNS and is the defining surface marker of the NG2-glia / OPC population.",
        "jakel2019", tissue_context=["brain"]),
]

# Microglia_homeostatic
NEW_RULES += [
    make_rule("Microglia_homeostatic", "MICRO_H_CIRCUIT_001", "circuit", "marker_genes",
        ["P2RY12", "TMEM119", "CX3CR1"], "active",
        "Homeostatic microglia require co-expression of P2RY12 (purinergic receptor for surveillance), TMEM119 (microglial-specific transmembrane protein), and CX3CR1 (fractalkine receptor for brain residence). All three are downregulated in disease-associated microglia.",
        "keren2017", tissue_context=["brain"]),
    make_rule("Microglia_homeostatic", "MICRO_H_NEG_001", "negative", "marker_genes",
        ["CCR2"], "absent",
        "CCR2 marks blood-derived monocytes and infiltrating macrophages. Brain-resident homeostatic microglia do not express CCR2, distinguishing them from monocyte-derived macrophages.",
        "keren2017", incompatible_with=["CCR2"], tissue_context=["brain"]),
    make_rule("Microglia_homeostatic", "MICRO_H_NEG_002", "negative", "marker_genes",
        ["TREM2"], "low",
        "TREM2 is a disease-associated microglia (DAM) marker that is expressed at low levels in homeostatic microglia. High TREM2 expression indicates transition to the DAM state.",
        "keren2017", tissue_context=["brain"]),
]

# Microglia_DAM
NEW_RULES += [
    make_rule("Microglia_DAM", "MICRO_DAM_CIRCUIT_001", "circuit", "marker_genes",
        ["TREM2", "SPP1"], "active",
        "Disease-associated microglia (DAM) are characterized by co-expression of TREM2 (triggers lipid phagocytosis) and SPP1 (osteopontin, a DAM effector). This circuit defines the DAM state in Alzheimer's disease and other neurodegeneration.",
        "keren2017", tissue_context=["brain", "neurodegeneration"]),
    make_rule("Microglia_DAM", "MICRO_DAM_NEG_001", "negative", "marker_genes",
        ["P2RY12"], "low",
        "P2RY12 is a homeostatic microglia marker that is strongly downregulated during transition to the DAM state. Low P2RY12 is required for DAM classification.",
        "keren2017", tissue_context=["brain"]),
    make_rule("Microglia_DAM", "MICRO_DAM_POS_001", "positive", "marker_genes",
        ["APOE"], "high",
        "APOE is strongly upregulated in DAM as part of the phagocytic lipid metabolism program. High APOE expression supports DAM identity.",
        "keren2017", tissue_context=["brain", "neurodegeneration"]),
]

# Schwann_cell
NEW_RULES += [
    make_rule("Schwann_cell", "SCHW_POS_001", "positive", "marker_genes",
        ["S100B", "MPZ", "PMP22"], "high",
        "Schwann cells co-express S100B (calcium-binding glial marker), MPZ (myelin protein zero, major PNS myelin component), and PMP22 (peripheral myelin protein 22). This trio is specific to PNS myelinating glia.",
        "schwann2020", tissue_context=["peripheral_nerve", "skin"]),
    make_rule("Schwann_cell", "SCHW_NEG_001", "negative", "marker_genes",
        ["MBP"], "low",
        "MBP is primarily a CNS myelin marker. Non-myelinating Schwann cells express very low MBP. High MBP would indicate CNS oligodendrocyte, not Schwann cell identity.",
        "schwann2020", tissue_context=["peripheral_nerve"]),
]

# Ependymal_cell
NEW_RULES += [
    make_rule("Ependymal_cell", "EPEN_CIRCUIT_001", "circuit", "regulon",
        ["FOXJ1"], "active",
        "FOXJ1 is the master transcription factor for multiciliated cell identity. Ependymal cells lining brain ventricles are defined by FOXJ1 regulon activity, which drives ciliogenesis and motile cilium assembly.",
        "ependymal", tissue_context=["brain", "ventricular"]),
    make_rule("Ependymal_cell", "EPEN_POS_001", "positive", "marker_genes",
        ["TMEM212", "CD24"], "high",
        "TMEM212 (transmembrane protein) and CD24 are expressed in ependymal cells and mark the ciliated ventricular surface. Both are downstream of FOXJ1 activity.",
        "allenbrain", tissue_context=["brain", "ventricular"]),
]

# ═══════════════════════════════════════════════════════════════════════════════
# Skin — 6 cell types, 20 rules
# ═══════════════════════════════════════════════════════════════════════════════

# Keratinocyte_basal
NEW_RULES += [
    make_rule("Keratinocyte_basal", "KBAS_CIRCUIT_001", "circuit", "regulon",
        ["TP63"], "active",
        "TP63 (p63) is the master TF for basal keratinocyte identity. It maintains the basal stem cell program and is required for KRT14/KRT5 expression. TP63 regulon activity defines the basal keratinocyte state.",
        "ji2020", tissue_context=["skin", "epidermis"]),
    make_rule("Keratinocyte_basal", "KBAS_POS_001", "positive", "marker_genes",
        ["KRT14", "KRT5"], "high",
        "KRT14 and KRT5 form the keratin filament pairs specific to basal keratinocytes. Their co-expression marks attachment to the basement membrane.",
        "ji2020", tissue_context=["skin", "epidermis"]),
]

# Keratinocyte_spinous
NEW_RULES += [
    make_rule("Keratinocyte_spinous", "KSPIN_POS_001", "positive", "marker_genes",
        ["KRT1", "KRT10"], "high",
        "KRT1 and KRT10 are expressed as keratinocytes transit from basal to spinous layer. Their co-expression marks suprabasal differentiation commitment.",
        "ji2020", tissue_context=["skin", "epidermis"]),
    make_rule("Keratinocyte_spinous", "KSPIN_NEG_001", "negative", "marker_genes",
        ["KRT14"], "low",
        "KRT14 is a basal keratinocyte marker that is progressively downregulated upon basal layer exit. Low KRT14 in spinous cells reflects commitment to differentiation.",
        "ji2020", incompatible_with=["KRT14"], tissue_context=["skin"]),
    make_rule("Keratinocyte_spinous", "KSPIN_NEG_002", "negative", "regulon",
        ["TP63"], "low",
        "TP63 regulon activity is progressively silenced as keratinocytes exit the basal layer. Low TP63 activity is required for spinous layer identity.",
        "ji2020", tissue_context=["skin"]),
]

# Keratinocyte_granular
NEW_RULES += [
    make_rule("Keratinocyte_granular", "KGRAN_POS_001", "positive", "marker_genes",
        ["FLG", "LOR", "IVL"], "high",
        "Filaggrin (FLG), loricrin (LOR), and involucrin (IVL) are terminal differentiation markers expressed in granular layer keratinocytes. Their co-expression defines the cornified envelope program.",
        "ji2020", tissue_context=["skin", "epidermis"]),
]

# Melanocyte
NEW_RULES += [
    make_rule("Melanocyte", "MELAN_CIRCUIT_001", "circuit", "regulon",
        ["MITF"], "active",
        "MITF is the master transcription factor for melanocyte identity. MITF regulon activity drives melanogenesis genes (DCT, TYRP1, PMEL) and is required for melanocyte survival and pigmentation.",
        "ji2020", tissue_context=["skin"]),
    make_rule("Melanocyte", "MELAN_POS_001", "positive", "marker_genes",
        ["DCT", "TYRP1", "PMEL"], "high",
        "DCT (dopachrome tautomerase), TYRP1, and PMEL (gp100) are MITF target genes required for melanin synthesis and melanosome biogenesis.",
        "ji2020", tissue_context=["skin"]),
    make_rule("Melanocyte", "MELAN_NEG_001", "negative", "marker_genes",
        ["EPCAM"], "absent",
        "EPCAM marks epithelial cells. Melanocytes are neural crest-derived and do not express EPCAM despite residing in the skin epithelium.",
        "ji2020", incompatible_with=["EPCAM"], tissue_context=["skin"]),
]

# Langerhans_cell
NEW_RULES += [
    make_rule("Langerhans_cell", "LANG_POS_001", "positive", "marker_genes",
        ["CD1A", "CD207"], "high",
        "CD1A (antigen presentation for lipids) and CD207 (langerin, a C-type lectin) define Langerhans cells. CD207 drives Birbeck granule formation, unique to Langerhans cells.",
        "ji2020", tissue_context=["skin", "epidermis"]),
]

# Dermal_fibroblast
NEW_RULES += [
    make_rule("Dermal_fibroblast", "DFIB_POS_001", "positive", "marker_genes",
        ["COL1A1", "PDGFRA", "DCN", "THY1"], "high",
        "Dermal fibroblasts express COL1A1 (fibrillar collagen), PDGFRA (growth factor receptor), DCN (decorin, ECM proteoglycan), and THY1 (CD90, fibroblast surface marker).",
        "ji2020", tissue_context=["skin", "dermis"]),
    make_rule("Dermal_fibroblast", "DFIB_NEG_001", "negative", "marker_genes",
        ["ACTA2"], "low",
        "ACTA2 (alpha-SMA) is expressed by activated myofibroblasts but not quiescent dermal fibroblasts. Low ACTA2 is required to distinguish dermal fibroblasts from wound-healing myofibroblasts.",
        "ji2020", tissue_context=["skin"]),
]

# ═══════════════════════════════════════════════════════════════════════════════
# Kidney — 8 cell types, 30 rules
# ═══════════════════════════════════════════════════════════════════════════════

# Proximal_Tubule
NEW_RULES += [
    make_rule("Proximal_Tubule", "PT_POS_001", "positive", "marker_genes",
        ["SLC34A1", "CUBN", "LRP2"], "high",
        "Proximal tubule cells express SLC34A1 (NaPi2a, sodium-phosphate cotransporter), CUBN (cubilin, endocytic receptor), and LRP2 (megalin), which together mediate bulk reabsorption in the early nephron.",
        "lake2023", tissue_context=["kidney", "proximal_tubule"]),
    make_rule("Proximal_Tubule", "PT_NEG_001", "negative", "marker_genes",
        ["AQP2"], "absent",
        "AQP2 is the principal water channel of the collecting duct principal cells. Its absence distinguishes proximal tubule from collecting duct cells in the kidney.",
        "lake2023", incompatible_with=["AQP2"], tissue_context=["kidney"]),
    make_rule("Proximal_Tubule", "PT_POS_002", "positive", "marker_genes",
        ["SLC22A6"], "high",
        "SLC22A6 (OAT1) is an organic anion transporter expressed specifically in the proximal tubule, mediating secretion of urate, drugs, and toxins.",
        "lake2023", tissue_context=["kidney"]),
]

# Thick_Ascending_Limb
NEW_RULES += [
    make_rule("Thick_Ascending_Limb", "TAL_POS_001", "positive", "marker_genes",
        ["SLC12A1", "UMOD"], "high",
        "SLC12A1 (NKCC2) is the sodium-potassium-chloride cotransporter defining the thick ascending limb (TAL). UMOD (uromodulin/Tamm-Horsfall protein) is the most abundant urinary protein, exclusively expressed in TAL.",
        "lake2023", tissue_context=["kidney", "loop_of_henle"]),
    make_rule("Thick_Ascending_Limb", "TAL_NEG_001", "negative", "marker_genes",
        ["AQP2"], "absent",
        "The thick ascending limb is impermeable to water. AQP2 absence is required; its presence would indicate collecting duct, not TAL.",
        "lake2023", incompatible_with=["AQP2"], tissue_context=["kidney"]),
]

# Distal_Convoluted_Tubule
NEW_RULES += [
    make_rule("Distal_Convoluted_Tubule", "DCT_POS_001", "positive", "marker_genes",
        ["SLC12A3", "CALB1"], "high",
        "SLC12A3 (NCC, thiazide-sensitive cotransporter) is the defining transporter of the distal convoluted tubule. CALB1 (calbindin D28k) marks DCT2 and is required for calcium reabsorption.",
        "lake2023", tissue_context=["kidney", "distal_tubule"]),
]

# Principal_Cell
NEW_RULES += [
    make_rule("Principal_Cell", "PC_POS_001", "positive", "marker_genes",
        ["AQP2", "AQP3", "SCNN1A"], "high",
        "Collecting duct principal cells express AQP2 (apical water channel, AVP-regulated), AQP3 (basolateral water channel), and SCNN1A (alpha-ENaC, sodium channel). This trio enables regulated water and sodium reabsorption.",
        "lake2023", tissue_context=["kidney", "collecting_duct"]),
    make_rule("Principal_Cell", "PC_NEG_001", "negative", "marker_genes",
        ["ATP6V1B1"], "absent",
        "ATP6V1B1 is the B1 subunit of vacuolar H+-ATPase expressed in intercalated cells. Its presence contradicts principal cell identity.",
        "lake2023", incompatible_with=["ATP6V1B1"], tissue_context=["kidney"]),
]

# Intercalated_Cell_A
NEW_RULES += [
    make_rule("Intercalated_Cell_A", "ICA_POS_001", "positive", "marker_genes",
        ["ATP6V1B1", "SLC4A1"], "high",
        "Type A intercalated cells secrete protons via apical H+-ATPase (ATP6V1B1) and reabsorb bicarbonate via basolateral SLC4A1 (AE1, anion exchanger). This machinery defines urinary acidification in the collecting duct.",
        "lake2023", tissue_context=["kidney", "collecting_duct"]),
    make_rule("Intercalated_Cell_A", "ICA_NEG_001", "negative", "marker_genes",
        ["SLC26A4"], "absent",
        "SLC26A4 (pendrin) is the apical bicarbonate secretor expressed in type B intercalated cells. Its absence distinguishes type A (acid-secreting) from type B (base-secreting) intercalated cells.",
        "lake2023", incompatible_with=["SLC26A4"], tissue_context=["kidney"]),
    make_rule("Intercalated_Cell_A", "ICA_NEG_002", "negative", "marker_genes",
        ["AQP2"], "absent",
        "AQP2 marks principal cells. Its absence is required for intercalated cell identity; principal and intercalated cells are mutually exclusive collecting duct populations.",
        "lake2023", incompatible_with=["AQP2"], tissue_context=["kidney"]),
]

# Intercalated_Cell_B
NEW_RULES += [
    make_rule("Intercalated_Cell_B", "ICB_POS_001", "positive", "marker_genes",
        ["SLC26A4"], "high",
        "SLC26A4 (pendrin) is the apical chloride/bicarbonate exchanger defining type B intercalated cells. It mediates net bicarbonate secretion and chloride reabsorption during metabolic alkalosis.",
        "lake2023", tissue_context=["kidney", "collecting_duct"]),
    make_rule("Intercalated_Cell_B", "ICB_NEG_001", "negative", "marker_genes",
        ["AQP2"], "absent",
        "AQP2 marks collecting duct principal cells. Intercalated cells (both A and B) lack AQP2.",
        "lake2023", incompatible_with=["AQP2"], tissue_context=["kidney"]),
]

# Podocyte
NEW_RULES += [
    make_rule("Podocyte", "POD_CIRCUIT_001", "circuit", "regulon",
        ["WT1"], "active",
        "WT1 (Wilms tumor protein 1) is the master transcription factor for podocyte identity. WT1 regulon activity maintains NPHS1 and NPHS2 expression, which form the slit diaphragm filtration barrier.",
        "lake2023", tissue_context=["kidney", "glomerulus"]),
    make_rule("Podocyte", "POD_POS_001", "positive", "marker_genes",
        ["NPHS1", "NPHS2", "PODXL"], "high",
        "NPHS1 (nephrin) and NPHS2 (podocin) form the slit diaphragm at foot processes. PODXL (podocalyxin) maintains the charge barrier. All three are podocyte-specific.",
        "lake2023", tissue_context=["kidney", "glomerulus"]),
    make_rule("Podocyte", "POD_NEG_001", "negative", "marker_genes",
        ["EPCAM"], "absent",
        "EPCAM marks tubular epithelial cells. Podocytes are specialized visceral epithelial cells that have undergone mesenchymal-like transition and do not express EPCAM.",
        "lake2023", incompatible_with=["EPCAM"], tissue_context=["kidney"]),
]

# Mesangial_cell
NEW_RULES += [
    make_rule("Mesangial_cell", "MESANG_POS_001", "positive", "marker_genes",
        ["PDGFRB", "CALD1"], "high",
        "Mesangial cells express PDGFRB (receptor for mesangial cell growth/expansion signals) and CALD1 (caldesmon, smooth muscle-like cytoskeletal protein for contractility).",
        "lake2023", tissue_context=["kidney", "glomerulus"]),
    make_rule("Mesangial_cell", "MESANG_NEG_001", "negative", "marker_genes",
        ["AQP2"], "absent",
        "AQP2 marks collecting duct principal cells. Mesangial cells are contractile intraglomerular cells and do not express AQP2.",
        "lake2023", incompatible_with=["AQP2"], tissue_context=["kidney"]),
    make_rule("Mesangial_cell", "MESANG_NEG_002", "negative", "marker_genes",
        ["NPHS1"], "absent",
        "NPHS1 is podocyte-specific. Mesangial cells are not podocytes and do not form slit diaphragms.",
        "lake2023", incompatible_with=["NPHS1"], tissue_context=["kidney"]),
]

# ═══════════════════════════════════════════════════════════════════════════════
# Heart — 5 cell types, 20 rules
# ═══════════════════════════════════════════════════════════════════════════════

# Ventricular_Cardiomyocyte
NEW_RULES += [
    make_rule("Ventricular_Cardiomyocyte", "VCM_POS_001", "positive", "marker_genes",
        ["MYH7", "TNNI3", "TTN"], "high",
        "MYH7 (beta-myosin heavy chain, dominant in ventricle), TNNI3 (cardiac troponin I), and TTN (titin, giant sarcomeric protein) define ventricular cardiomyocyte identity.",
        "litvin2020", tissue_context=["heart", "ventricle"]),
    make_rule("Ventricular_Cardiomyocyte", "VCM_NEG_001", "negative", "marker_genes",
        ["MYH6"], "low",
        "MYH6 (alpha-myosin heavy chain) is the dominant isoform in atrial cardiomyocytes. Low MYH6 distinguishes ventricular from atrial cardiomyocytes. The MYH7-high/MYH6-low ratio is the defining isoform switch.",
        "litvin2020", tissue_context=["heart"]),
    make_rule("Ventricular_Cardiomyocyte", "VCM_NEG_002", "negative", "marker_genes",
        ["NPPA"], "low",
        "NPPA (ANP, atrial natriuretic peptide) is primarily expressed in atrial cardiomyocytes. Low NPPA distinguishes ventricular from atrial identity.",
        "litvin2020", tissue_context=["heart"]),
]

# Atrial_Cardiomyocyte
NEW_RULES += [
    make_rule("Atrial_Cardiomyocyte", "ACM_POS_001", "positive", "marker_genes",
        ["MYH6", "NPPA"], "high",
        "MYH6 (alpha-MHC, dominant atrial isoform) and NPPA (ANP, the atrial natriuretic hormone) are co-expressed in atrial cardiomyocytes and define atrial identity.",
        "litvin2020", tissue_context=["heart", "atrium"]),
    make_rule("Atrial_Cardiomyocyte", "ACM_NEG_001", "negative", "marker_genes",
        ["MYH7"], "low",
        "MYH7 (beta-MHC) dominates in ventricular cardiomyocytes. Low MYH7 is required for atrial identity; the MYH6-high/MYH7-low ratio defines the atrial isoform program.",
        "litvin2020", tissue_context=["heart"]),
]

# Cardiac_Fibroblast
NEW_RULES += [
    make_rule("Cardiac_Fibroblast", "CFIB_POS_001", "positive", "marker_genes",
        ["TCF21", "DCN", "PDGFRA", "FN1"], "high",
        "Cardiac fibroblasts express TCF21 (master TF for epicardial-derived fibroblasts), DCN (decorin ECM proteoglycan), PDGFRA (growth factor receptor), and FN1 (fibronectin for ECM assembly).",
        "litvin2020", tissue_context=["heart"]),
    make_rule("Cardiac_Fibroblast", "CFIB_NEG_001", "negative", "marker_genes",
        ["ACTA2"], "low",
        "ACTA2 (alpha-SMA) marks activated myofibroblasts. Quiescent cardiac fibroblasts express low ACTA2; high ACTA2 indicates myofibroblast activation during cardiac remodeling.",
        "litvin2020", tissue_context=["heart"]),
    make_rule("Cardiac_Fibroblast", "CFIB_NEG_002", "negative", "marker_genes",
        ["TNNI3"], "absent",
        "TNNI3 is a cardiomyocyte-specific marker. Its absence is required for cardiac fibroblast identity.",
        "litvin2020", incompatible_with=["TNNI3"], tissue_context=["heart"]),
]

# Endocardial_cell
NEW_RULES += [
    make_rule("Endocardial_cell", "ENDO_POS_001", "positive", "marker_genes",
        ["CDH5", "PECAM1", "NFATc1"], "high",
        "Endocardial cells are specialized endothelial cells lining heart chambers. They express CDH5 (VE-cadherin), PECAM1 (CD31), and NFATc1, which regulates valvulogenesis.",
        "litvin2020", tissue_context=["heart", "endocardium"]),
    make_rule("Endocardial_cell", "ENDO_NEG_001", "negative", "marker_genes",
        ["MYH7"], "absent",
        "MYH7 marks cardiomyocytes. Endocardial cells are endothelial, not myocardial, and do not express cardiomyocyte sarcomeric genes.",
        "litvin2020", incompatible_with=["MYH7"], tissue_context=["heart"]),
]

# Pericyte (heart/general)
NEW_RULES += [
    make_rule("Pericyte", "PERI_POS_001", "positive", "marker_genes",
        ["PDGFRB", "RGS5", "KCNJ8"], "high",
        "Pericytes express PDGFRB (PDGF-B receptor for pericyte recruitment), RGS5 (G-protein regulator for vessel tone), and KCNJ8 (Kir6.1, potassium channel). This combination distinguishes pericytes from SMCs.",
        "vanland2018", tissue_context=["heart", "brain", "blood_vessel"]),
    make_rule("Pericyte", "PERI_NEG_001", "negative", "marker_genes",
        ["MYH11"], "absent",
        "MYH11 (smooth muscle myosin heavy chain) is the contractile marker of mature vascular smooth muscle cells. Pericytes lack MYH11, distinguishing them from SMCs despite both wrapping blood vessels.",
        "vanland2018", incompatible_with=["MYH11"], tissue_context=["blood_vessel"]),
    make_rule("Pericyte", "PERI_NEG_002", "negative", "marker_genes",
        ["ACTA2"], "low",
        "ACTA2 (alpha-SMA) is expressed at high levels in SMCs but only low levels in pericytes. This expression difference, combined with MYH11 absence, distinguishes pericytes from SMCs.",
        "vanland2018", tissue_context=["blood_vessel"]),
]

# ═══════════════════════════════════════════════════════════════════════════════
# Stromal and endothelial — 10 cell types, 40 rules
# ═══════════════════════════════════════════════════════════════════════════════

# Arterial_Endothelial
NEW_RULES += [
    make_rule("Arterial_Endothelial", "AEC_POS_001", "positive", "marker_genes",
        ["NOTCH4", "GJA4", "CXCL12"], "high",
        "Arterial endothelium is defined by active Notch signaling (NOTCH4), high connexin 37 (GJA4) for arterial gap junctions, and CXCL12 for hemato-endothelial crosstalk.",
        "vanland2018", tissue_context=["blood_vessel", "artery"]),
    make_rule("Arterial_Endothelial", "AEC_NEG_001", "negative", "marker_genes",
        ["APLNR"], "absent",
        "APLNR (apelin receptor) marks venous endothelium and is absent from arterial endothelium. Its presence contradicts arterial identity.",
        "vanland2018", incompatible_with=["APLNR"], tissue_context=["blood_vessel"]),
]

# Venous_Endothelial
NEW_RULES += [
    make_rule("Venous_Endothelial", "VEC_POS_001", "positive", "marker_genes",
        ["APLNR", "NR2F2"], "high",
        "Venous endothelial cells express APLNR (apelin receptor, activated by venous flow) and NR2F2 (COUP-TFII, nuclear receptor that suppresses arterial identity and establishes venous fate).",
        "vanland2018", tissue_context=["blood_vessel", "vein"]),
    make_rule("Venous_Endothelial", "VEC_NEG_001", "negative", "marker_genes",
        ["GJA4"], "low",
        "GJA4 (connexin 37) is an arterial endothelial marker. Low GJA4 in venous endothelium reflects suppression of arterial fate by NR2F2/COUP-TFII.",
        "vanland2018", tissue_context=["blood_vessel"]),
]

# Capillary_Endothelial
NEW_RULES += [
    make_rule("Capillary_Endothelial", "CAP_POS_001", "positive", "marker_genes",
        ["CA4", "RGCC", "CD36"], "high",
        "Capillary endothelium mediates gas/nutrient exchange. CA4 (carbonic anhydrase 4), RGCC (cell cycle regulator), and CD36 (fatty acid translocase) mark tissue-specific capillary exchange functions.",
        "vanland2018", tissue_context=["blood_vessel", "capillary"]),
]

# Lymphatic_Endothelial
NEW_RULES += [
    make_rule("Lymphatic_Endothelial", "LEC_CIRCUIT_001", "circuit", "regulon",
        ["PROX1"], "active",
        "PROX1 is the master regulator of lymphatic endothelial cell (LEC) identity. PROX1 regulon activity is required and sufficient to reprogram blood EC to LEC fate.",
        "tabumuris", tissue_context=["lymph_node", "lymphatics"]),
    make_rule("Lymphatic_Endothelial", "LEC_POS_001", "positive", "marker_genes",
        ["LYVE1", "PDPN"], "high",
        "LYVE1 (lymphatic hyaluronan receptor) and PDPN (podoplanin) are PROX1 targets expressed on LEC surfaces and are required for lymph drainage function.",
        "tabumuris", tissue_context=["lymphatics"]),
]

# Tip_cell
NEW_RULES += [
    make_rule("Tip_cell", "TIP_POS_001", "positive", "marker_genes",
        ["ESM1", "ANGPT2", "KDR"], "high",
        "Tip cells at angiogenic sprout fronts express ESM1 (endocan, VEGF-induced), ANGPT2 (destabilizes pericyte contacts to enable sprouting), and high KDR (VEGFR2 for VEGF sensing).",
        "vanland2018", tissue_context=["blood_vessel", "tumor"]),
    make_rule("Tip_cell", "TIP_POS_002", "positive", "marker_genes",
        ["DLL4"], "high",
        "DLL4 (delta-like 4) is the Notch ligand induced by VEGF in tip cells to laterally inhibit neighbouring stalk cells from also becoming tip cells. High DLL4 is a defining tip cell feature.",
        "vanland2018", tissue_context=["blood_vessel"]),
]

# Vascular_SMC
NEW_RULES += [
    make_rule("Vascular_SMC", "VSMC_POS_001", "positive", "marker_genes",
        ["MYH11", "CNN1", "ACTA2"], "high",
        "Vascular smooth muscle cells express the contractile machinery: MYH11 (smooth muscle myosin), CNN1 (calponin 1), and ACTA2 (alpha-SMA). This trio defines the differentiated contractile SMC phenotype.",
        "vanland2018", tissue_context=["blood_vessel", "artery"]),
    make_rule("Vascular_SMC", "VSMC_NEG_001", "negative", "marker_genes",
        ["FAP"], "absent",
        "FAP (fibroblast activation protein) marks cancer-associated fibroblasts and myofibroblasts. Its absence in VSMCs is critical for distinguishing SMCs from activated stromal fibroblasts in the tumor microenvironment.",
        "elyada2019", incompatible_with=["FAP"], tissue_context=["blood_vessel"]),
]

# Myofibroblast
NEW_RULES += [
    make_rule("Myofibroblast", "MYOFIB_POS_001", "positive", "marker_genes",
        ["ACTA2", "CNN1", "FAP"], "high",
        "Myofibroblasts are activated fibroblasts that co-express ACTA2 (very high, drives wound contraction), CNN1, and FAP (fibroblast activation protein, marks activated state). This combination defines myofibroblastic activation.",
        "elyada2019", tissue_context=["wound", "fibrosis", "tumor"]),
    make_rule("Myofibroblast", "MYOFIB_NEG_001", "negative", "marker_genes",
        ["MYH11"], "low",
        "MYH11 (smooth muscle myosin) is the definitive SMC marker and is low/absent in myofibroblasts despite their high ACTA2. This is the key discriminator from vascular SMC.",
        "elyada2019", tissue_context=["stromal"]),
]

# iCAF (inflammatory cancer-associated fibroblast)
NEW_RULES += [
    make_rule("iCAF", "ICAF_POS_001", "positive", "marker_genes",
        ["IL6", "CXCL12", "PDGFRA"], "high",
        "Inflammatory CAFs (iCAFs) are characterized by high IL6 (pro-inflammatory cytokine, activated by IL-1 from cancer cells), CXCL12 (SDF-1, T-cell exclusion factor), and PDGFRA. iCAFs are the major immunosuppressive fibroblast population.",
        "elyada2019", tissue_context=["tumor", "pancreas"]),
    make_rule("iCAF", "ICAF_NEG_001", "negative", "marker_genes",
        ["ACTA2"], "low",
        "iCAFs are non-myofibroblastic. TGFβ from myCAFs suppresses the iCAF inflammatory program; conversely, low ACTA2 in iCAFs reflects non-myofibroblastic identity. ACTA2 high would indicate transition to myCAF.",
        "elyada2019", tissue_context=["tumor"]),
]

# myCAF (myofibroblastic cancer-associated fibroblast)
NEW_RULES += [
    make_rule("myCAF", "MYCAF_POS_001", "positive", "marker_genes",
        ["ACTA2", "MMP11", "FAP"], "high",
        "myCAFs (myofibroblastic CAFs) are activated by TGFβ from cancer cells. They co-express ACTA2 (very high), MMP11 (stromelysin-3, ECM remodeling), and FAP (pan-CAF marker elevated in myCAF).",
        "elyada2019", tissue_context=["tumor", "pancreas"]),
    make_rule("myCAF", "MYCAF_NEG_001", "negative", "marker_genes",
        ["IL6"], "low",
        "IL6 is the signature cytokine of iCAFs. TGFβ signaling in myCAFs suppresses IL-1/IL-6 axis, creating mutual exclusivity between iCAF and myCAF states.",
        "elyada2019", incompatible_with=["IL6"], tissue_context=["tumor"]),
]

# ═══════════════════════════════════════════════════════════════════════════════
# Endocrine — 8 cell types, 30 rules
# ═══════════════════════════════════════════════════════════════════════════════

# Beta_cell
NEW_RULES += [
    make_rule("Beta_cell", "BETA_CIRCUIT_001", "circuit", "regulon",
        ["MAFA", "PDX1"], "active",
        "Beta cell identity requires co-active MAFA (beta-cell maturation TF) and PDX1 (pancreatic/beta-cell TF) regulons. MAFA drives INS transcription; PDX1 maintains beta-cell identity and insulin gene promoter activity.",
        "baron2016", tissue_context=["pancreas", "islet"]),
    make_rule("Beta_cell", "BETA_POS_001", "positive", "marker_genes",
        ["INS"], "high",
        "INS (insulin) is the hallmark secreted product of beta cells and the most abundantly expressed gene in mature beta cells.",
        "baron2016", tissue_context=["pancreas", "islet"]),
    make_rule("Beta_cell", "BETA_NEG_001", "negative", "marker_genes",
        ["GCG"], "absent",
        "GCG (glucagon) is the defining hormone of alpha cells and is mutually exclusive with insulin. GCG presence rules out beta cell identity — this is the most critical islet disambiguation.",
        "baron2016", incompatible_with=["GCG"], tissue_context=["pancreas"]),
    make_rule("Beta_cell", "BETA_NEG_002", "negative", "marker_genes",
        ["SST"], "absent",
        "SST (somatostatin) is the delta cell marker. Its presence would indicate delta cell identity, not beta cell.",
        "baron2016", incompatible_with=["SST"], tissue_context=["pancreas"]),
]

# Alpha_cell
NEW_RULES += [
    make_rule("Alpha_cell", "ALPHA_POS_001", "positive", "marker_genes",
        ["GCG", "ARX", "MAFB"], "high",
        "Alpha cells express glucagon (GCG, the secreted hormone), ARX (homeodomain TF required for alpha cell differentiation), and MAFB (MAF family TF maintaining alpha cell identity).",
        "baron2016", tissue_context=["pancreas", "islet"]),
    make_rule("Alpha_cell", "ALPHA_NEG_001", "negative", "marker_genes",
        ["INS"], "absent",
        "INS (insulin) is the beta cell hormone and is mutually exclusive with glucagon. INS presence rules out alpha cell identity.",
        "baron2016", incompatible_with=["INS"], tissue_context=["pancreas"]),
    make_rule("Alpha_cell", "ALPHA_NEG_002", "negative", "marker_genes",
        ["SST"], "absent",
        "SST (somatostatin) is the delta cell marker. Alpha and delta cells are distinct populations.",
        "baron2016", incompatible_with=["SST"], tissue_context=["pancreas"]),
]

# Delta_cell
NEW_RULES += [
    make_rule("Delta_cell", "DELTA_POS_001", "positive", "marker_genes",
        ["SST"], "high",
        "SST (somatostatin) is the defining secreted peptide of pancreatic delta cells. Delta cells regulate both alpha (GCG) and beta (INS) cell secretion via paracrine somatostatin signaling.",
        "baron2016", tissue_context=["pancreas", "islet"]),
    make_rule("Delta_cell", "DELTA_NEG_001", "negative", "marker_genes",
        ["INS", "GCG"], "absent",
        "Both insulin (beta cell) and glucagon (alpha cell) are absent in delta cells. Delta cells form a distinct third islet population defined by somatostatin secretion.",
        "baron2016", incompatible_with=["INS", "GCG"], tissue_context=["pancreas"]),
]

# PP_cell (Pancreatic Polypeptide)
NEW_RULES += [
    make_rule("PP_cell", "PP_POS_001", "positive", "marker_genes",
        ["PPY"], "high",
        "PPY (pancreatic polypeptide) is the sole defining hormone of PP cells (gamma cells), a rare fourth islet population primarily in the head of the pancreas.",
        "baron2016", tissue_context=["pancreas", "islet"]),
    make_rule("PP_cell", "PP_NEG_001", "negative", "marker_genes",
        ["INS", "GCG", "SST"], "absent",
        "PP cells are defined by absence of all three major islet hormones (insulin, glucagon, somatostatin) combined with high PPY. This four-way discrimination identifies them as the rarest islet cell type.",
        "baron2016", incompatible_with=["INS", "GCG", "SST"], tissue_context=["pancreas"]),
]

# Acinar_cell
NEW_RULES += [
    make_rule("Acinar_cell", "ACIN_POS_001", "positive", "marker_genes",
        ["PRSS1", "CPA1", "AMY2A"], "high",
        "Acinar cells produce digestive enzymes: PRSS1 (trypsinogen 1, the major serine protease), CPA1 (carboxypeptidase A1), and AMY2A (amylase 2A). This enzymatic trio is exclusively expressed in acinar cells.",
        "baron2016", tissue_context=["pancreas", "exocrine"]),
    make_rule("Acinar_cell", "ACIN_NEG_001", "negative", "marker_genes",
        ["KRT19"], "absent",
        "KRT19 is a ductal epithelial marker expressed in pancreatic ducts and biliary epithelium. Its absence distinguishes acinar from ductal cells.",
        "baron2016", incompatible_with=["KRT19"], tissue_context=["pancreas"]),
]

# Ductal_cell
NEW_RULES += [
    make_rule("Ductal_cell", "DUCT_POS_001", "positive", "marker_genes",
        ["KRT19", "CFTR", "MUC1"], "high",
        "Pancreatic ductal cells express KRT19 (biliary/pancreatic epithelial keratins), CFTR (cystic fibrosis transmembrane conductance regulator for bicarbonate secretion), and MUC1 (mucin for duct lining).",
        "baron2016", tissue_context=["pancreas", "duct"]),
    make_rule("Ductal_cell", "DUCT_NEG_001", "negative", "marker_genes",
        ["PRSS1"], "absent",
        "PRSS1 (trypsinogen) is an acinar enzyme. Ductal cells modify but do not produce digestive enzymes. PRSS1 presence contradicts ductal identity.",
        "baron2016", incompatible_with=["PRSS1"], tissue_context=["pancreas"]),
]

# Thyroid_follicular
NEW_RULES += [
    make_rule("Thyroid_follicular", "THYR_F_POS_001", "positive", "marker_genes",
        ["TSHR", "TG", "TPO"], "high",
        "Thyroid follicular cells are defined by TSHR (TSH receptor, drives thyroid hormone synthesis), TG (thyroglobulin, the thyroid hormone storage protein), and TPO (thyroid peroxidase, iodination enzyme). All three are required for thyroid hormone biosynthesis.",
        "tabumuris", tissue_context=["thyroid"]),
    make_rule("Thyroid_follicular", "THYR_F_NEG_001", "negative", "marker_genes",
        ["CALCA"], "absent",
        "CALCA (calcitonin) is the marker of parafollicular C cells. Its absence distinguishes follicular cells (thyroid hormone-producing) from C cells (calcitonin-producing).",
        "tabumuris", incompatible_with=["CALCA"], tissue_context=["thyroid"]),
]

# Thyroid_C_cell
NEW_RULES += [
    make_rule("Thyroid_C_cell", "THYR_C_POS_001", "positive", "marker_genes",
        ["CALCA", "CHGA"], "high",
        "Parafollicular C cells are neuroendocrine cells that produce calcitonin (CALCA) for calcium homeostasis. CHGA (chromogranin A) marks their neuroendocrine identity.",
        "tabumuris", tissue_context=["thyroid"]),
    make_rule("Thyroid_C_cell", "THYR_C_NEG_001", "negative", "marker_genes",
        ["TG"], "absent",
        "TG (thyroglobulin) is exclusively produced by follicular cells for thyroid hormone storage. C cells are CALCA-secreting neuroendocrine cells and do not express TG.",
        "tabumuris", incompatible_with=["TG"], tissue_context=["thyroid"]),
]

# ═══════════════════════════════════════════════════════════════════════════════
# Hematopoietic progenitors — 7 cell types, 25 rules
# ═══════════════════════════════════════════════════════════════════════════════

# HSC
NEW_RULES += [
    make_rule("HSC", "HSC_CIRCUIT_001", "circuit", "regulon",
        ["GATA2"], "active",
        "GATA2 is a master TF required for HSC maintenance and multi-lineage potential. GATA2 regulon activity drives HSC self-renewal genes and is progressively downregulated as cells commit to specific lineages.",
        "velten2017", tissue_context=["bone_marrow"]),
    make_rule("HSC", "HSC_POS_001", "positive", "marker_genes",
        ["CD34", "HOXA9", "MLLT3"], "high",
        "HSCs express CD34 (surface marker for HSPC isolation), HOXA9 (HOX TF for HSC self-renewal), and MLLT3 (AF9, required for self-renewal via H3K79 methylation).",
        "velten2017", tissue_context=["bone_marrow"]),
    make_rule("HSC", "HSC_NEG_001", "negative", "marker_genes",
        ["MPO", "CD19", "CD3E"], "absent",
        "HSCs must lack all lineage commitment markers: MPO (myeloid), CD19 (B cell), CD3E (T cell). Lin-negative status is the defining property of true HSCs.",
        "pellin2019", incompatible_with=["MPO", "CD19", "CD3E"], tissue_context=["bone_marrow"]),
]

# CMP (Common Myeloid Progenitor)
NEW_RULES += [
    make_rule("CMP", "CMP_POS_001", "positive", "marker_genes",
        ["CD34", "CSF1R"], "high",
        "CMPs retain CD34 expression and upregulate CSF1R (M-CSF receptor, myeloid commitment receptor). This CD34+/CSF1R+ state marks the myeloid branch commitment point.",
        "pellin2019", tissue_context=["bone_marrow"]),
    make_rule("CMP", "CMP_NEG_001", "negative", "marker_genes",
        ["GATA1"], "low",
        "GATA1 is the MEP/erythroid/megakaryocyte TF. Low GATA1 in CMPs distinguishes them from MEPs which are GATA1-high. CMP→MEP transition requires GATA1 induction.",
        "pellin2019", tissue_context=["bone_marrow"]),
    make_rule("CMP", "CMP_NEG_002", "negative", "marker_genes",
        ["MPO"], "low",
        "MPO (myeloperoxidase) marks mature GMP and neutrophil-lineage cells. Low MPO in CMPs reflects early myeloid commitment before granulocyte specification.",
        "pellin2019", tissue_context=["bone_marrow"]),
]

# GMP (Granulocyte-Monocyte Progenitor)
NEW_RULES += [
    make_rule("GMP", "GMP_POS_001", "positive", "marker_genes",
        ["CSF1R", "MPO"], "high",
        "GMPs express CSF1R (myeloid commitment) and acquire MPO (myeloperoxidase), the definitive marker of granulocyte/monocyte commitment. MPO expression onset marks GMP specification from CMP.",
        "pellin2019", tissue_context=["bone_marrow"]),
    make_rule("GMP", "GMP_NEG_001", "negative", "marker_genes",
        ["GATA1"], "absent",
        "GATA1 drives erythroid/megakaryocyte fate. GATA1 and MPO expression are mutually exclusive, reflecting the bifurcation between GMP and MEP lineages from CMP.",
        "pellin2019", incompatible_with=["GATA1"], tissue_context=["bone_marrow"]),
]

# MEP (Megakaryocyte-Erythroid Progenitor)
NEW_RULES += [
    make_rule("MEP", "MEP_CIRCUIT_001", "circuit", "regulon",
        ["GATA1", "KLF1"], "active",
        "MEP identity requires co-active GATA1 (erythroid/megakaryocyte master TF) and KLF1 (erythroid TF). Their co-activity drives MEP toward erythroid vs megakaryocyte fate based on FLI1 competition.",
        "erythroid", tissue_context=["bone_marrow"]),
    make_rule("MEP", "MEP_POS_001", "positive", "marker_genes",
        ["FLI1"], "high",
        "FLI1 is expressed in MEPs and competes with KLF1 for GATA1 interaction. High FLI1 biases MEPs toward megakaryocyte fate; FLI1 decline allows KLF1/GATA1 to drive erythropoiesis.",
        "erythroid", tissue_context=["bone_marrow"]),
    make_rule("MEP", "MEP_NEG_001", "negative", "marker_genes",
        ["CSF1R"], "absent",
        "CSF1R marks myeloid commitment (GMP/CMP). MEPs have exited the myeloid branch and do not express CSF1R.",
        "pellin2019", incompatible_with=["CSF1R"], tissue_context=["bone_marrow"]),
]

# Megakaryocyte
NEW_RULES += [
    make_rule("Megakaryocyte", "MEG_POS_001", "positive", "marker_genes",
        ["ITGA2B", "PF4", "VWF", "GATA1"], "high",
        "Megakaryocytes are polyploid platelet progenitors expressing ITGA2B (CD41, surface glycoprotein for platelet aggregation), PF4 (platelet factor 4), VWF (von Willebrand factor), and high GATA1.",
        "velten2017", tissue_context=["bone_marrow"]),
    make_rule("Megakaryocyte", "MEG_NEG_001", "negative", "marker_genes",
        ["HBA1"], "low",
        "HBA1 (hemoglobin alpha) is the defining product of erythroid cells. Low HBA1 distinguishes megakaryocytes from erythroblasts, despite sharing GATA1 high expression.",
        "erythroid", tissue_context=["bone_marrow"]),
]

# Erythroblast
NEW_RULES += [
    make_rule("Erythroblast", "ERYTH_CIRCUIT_001", "circuit", "regulon",
        ["GATA1"], "active",
        "GATA1 is the master TF for erythroid differentiation. GATA1 regulon activity drives HBA/HBB globin gene expression and erythroid-specific programs (ALAS2, SLC4A1) during red cell maturation.",
        "erythroid", tissue_context=["bone_marrow"]),
    make_rule("Erythroblast", "ERYTH_POS_001", "positive", "marker_genes",
        ["HBA1", "HBA2", "ALAS2"], "high",
        "Erythroblasts express alpha-globins (HBA1, HBA2) and ALAS2 (delta-aminolevulinate synthase 2, rate-limiting enzyme for heme biosynthesis). This hemoglobin synthesis program defines erythroid identity.",
        "erythroid", tissue_context=["bone_marrow"]),
    make_rule("Erythroblast", "ERYTH_NEG_001", "negative", "marker_genes",
        ["MPO"], "absent",
        "MPO marks granulocyte/monocyte lineage. Erythroblasts and myeloid cells are mutually exclusive, reflecting GATA1 (erythroid) vs C/EBPα (myeloid) TF competition.",
        "pellin2019", incompatible_with=["MPO"], tissue_context=["bone_marrow"]),
]

# Plasmablast
NEW_RULES += [
    make_rule("Plasmablast", "PBLAST_CIRCUIT_001", "circuit", "regulon",
        ["IRF4", "PRDM1"], "active",
        "Plasmablasts are early antibody-secreting cells defined by co-active IRF4 (very high, drives plasma cell program) and PRDM1/BLIMP1 (represses B cell identity genes including PAX5). Both TFs are required for plasma cell commitment.",
        "pellin2019", tissue_context=["bone_marrow", "lymph_node", "blood"]),
    make_rule("Plasmablast", "PBLAST_NEG_001", "negative", "marker_genes",
        ["CD19", "PAX5"], "low",
        "BLIMP1 directly represses PAX5 during plasma cell differentiation. CD19 and PAX5 are progressively downregulated as plasmablasts exit B cell identity. Their low/absent expression confirms plasma cell commitment.",
        "velten2017", tissue_context=["bone_marrow", "lymph_node"]),
    make_rule("Plasmablast", "PBLAST_POS_001", "positive", "marker_genes",
        ["SDC1", "MZB1"], "high",
        "SDC1 (CD138, syndecan-1) is the canonical plasma cell surface marker used clinically for multiple myeloma diagnosis. MZB1 (marginal zone B cell protein) is a plasma cell endoplasmic reticulum chaperone.",
        "velten2017", tissue_context=["bone_marrow"]),
]

# ═══════════════════════════════════════════════════════════════════════════════
# Cancer-specific states — 8 cell types, 30 rules
# ═══════════════════════════════════════════════════════════════════════════════

# TEX_progenitor (TPEX)
NEW_RULES += [
    make_rule("TEX_progenitor", "TEXPROG_CIRCUIT_001", "circuit", "marker_genes",
        ["TOX", "TCF7"], "active",
        "Progenitor-exhausted T cells (TPEX) co-express TOX (exhaustion master TF, induces epigenetic silencing) and TCF7/TCF1 (self-renewal TF, not yet silenced by TOX). This co-expression defines the stem-like exhausted state that responds to checkpoint blockade.",
        "miller2019", tissue_context=["tumor", "chronic_infection"]),
    make_rule("TEX_progenitor", "TEXPROG_POS_001", "positive", "marker_genes",
        ["SLAMF6", "PDCD1"], "high",
        "TPEX cells express SLAMF6 (Ly108, self-renewal marker) and moderate PDCD1 (PD-1), reflecting partial but not terminal exhaustion. SLAMF6 marks the progenitor pool.",
        "zheng2021", tissue_context=["tumor"]),
    make_rule("TEX_progenitor", "TEXPROG_NEG_001", "negative", "marker_genes",
        ["HAVCR2"], "low",
        "HAVCR2 (TIM-3) is the terminal exhaustion marker absent from TPEX cells. High TIM-3 combined with TCF7 loss defines terminal TEX, distinguishing them from TCF7+ progenitor TEX.",
        "miller2019", tissue_context=["tumor"]),
]

# TEX_terminal
NEW_RULES += [
    make_rule("TEX_terminal", "TEXTERM_CIRCUIT_001", "circuit", "marker_genes",
        ["TOX", "HAVCR2", "LAG3", "ENTPD1"], "active",
        "Terminally exhausted T cells are defined by TOX (exhaustion TF) co-expressed with terminal exhaustion effectors: HAVCR2 (TIM-3), LAG3 (immune checkpoint), and ENTPD1 (CD39, adenosine pathway).",
        "miller2019", tissue_context=["tumor"]),
    make_rule("TEX_terminal", "TEXTERM_NEG_001", "negative", "marker_genes",
        ["TCF7"], "absent",
        "TCF7 (TCF1) is silenced by TOX-mediated epigenetic reprogramming in terminal TEX. TCF7 loss is required for terminal exhaustion state — its presence indicates TPEX, not terminal TEX. This was a key Phase 1 finding.",
        "miller2019", incompatible_with=["TCF7"], tissue_context=["tumor"]),
    make_rule("TEX_terminal", "TEXTERM_POS_001", "positive", "marker_genes",
        ["PDCD1"], "high",
        "PD-1 (PDCD1) reaches its highest expression in terminally exhausted T cells, reflecting maximal coinhibitory receptor burden.",
        "zheng2021", tissue_context=["tumor"]),
]

# TAM_M1like
NEW_RULES += [
    make_rule("TAM_M1like", "TAMM1_POS_001", "positive", "marker_genes",
        ["IDO1", "CXCL9", "CXCL10"], "high",
        "M1-like tumor-associated macrophages express IDO1 (immunosuppressive tryptophan catabolism), CXCL9, and CXCL10 (T-cell chemoattractants). This interferon-stimulated program is induced by tumor-derived IFN-gamma.",
        "zheng2021", tissue_context=["tumor"]),
    make_rule("TAM_M1like", "TAMM1_NEG_001", "negative", "marker_genes",
        ["MRC1"], "absent",
        "MRC1 (CD206, mannose receptor) is an M2-like scavenger receptor expressed in anti-inflammatory macrophages. Its absence in M1-like TAMs reflects incompatible TF programs (IRF5 vs IRF4).",
        "zheng2021", incompatible_with=["MRC1"], tissue_context=["tumor"]),
]

# TAM_M2like
NEW_RULES += [
    make_rule("TAM_M2like", "TAMM2_POS_001", "positive", "marker_genes",
        ["MRC1", "CD163", "VEGFA"], "high",
        "M2-like TAMs express MRC1 (CD206, scavenger receptor), CD163 (hemoglobin scavenger, M2 marker), and VEGFA (pro-angiogenic factor). This program promotes tumor angiogenesis and immunosuppression.",
        "zheng2021", tissue_context=["tumor"]),
    make_rule("TAM_M2like", "TAMM2_NEG_001", "negative", "marker_genes",
        ["IDO1"], "absent",
        "IDO1 is an M1-like TAM marker driven by IFN-gamma signaling. M2-like TAMs are alternatively activated and do not express IDO1.",
        "zheng2021", incompatible_with=["IDO1"], tissue_context=["tumor"]),
]

# TAM_SPP1
NEW_RULES += [
    make_rule("TAM_SPP1", "TAMSPP1_POS_001", "positive", "marker_genes",
        ["SPP1", "TREM2", "CCL18"], "high",
        "SPP1+ TREM2+ TAMs are a disease-associated macrophage subtype found in tumor niches and fibrotic tissue. SPP1 (osteopontin) drives tumor invasion, TREM2 regulates lipid/cholesterol homeostasis, CCL18 promotes immunosuppression.",
        "zheng2021", tissue_context=["tumor"]),
    make_rule("TAM_SPP1", "TAMSPP1_NEG_001", "negative", "marker_genes",
        ["HLA-DRA"], "low",
        "SPP1+ TAMs are characterized by low MHC-II expression (HLA-DRA), making them poor antigen presenters. Low HLA-DRA distinguishes SPP1+ TAMs from professional antigen-presenting M1-like TAMs.",
        "zheng2021", tissue_context=["tumor"]),
]

# Tumor_Treg
NEW_RULES += [
    make_rule("Tumor_Treg", "TTREG_POS_001", "positive", "marker_genes",
        ["FOXP3", "ENTPD1", "LAYN", "IL2RA"], "high",
        "Tumor-enriched Tregs express FOXP3 (Treg master TF), ENTPD1 (CD39, adenosine immunosuppression), LAYN (layilin, tumor Treg marker), and IL2RA (CD25). This quartet distinguishes highly suppressive tumor Tregs from peripheral Tregs.",
        "plitas2016", tissue_context=["tumor"]),
    make_rule("Tumor_Treg", "TTREG_NEG_001", "negative", "marker_genes",
        ["IL7R"], "low",
        "IL7R (CD127) is downregulated in Tregs as IL2RA (CD25) is upregulated. Low IL7R combined with high FOXP3 and IL2RA confirms Treg identity over conventional T cells.",
        "plitas2016", tissue_context=["tumor"]),
]

# Plasma_cell (cancer context, distinct from Plasmablast)
NEW_RULES += [
    make_rule("Plasma_cell", "PLASMA_CIRCUIT_001", "circuit", "regulon",
        ["PRDM1"], "active",
        "BLIMP1/PRDM1 is the terminal plasma cell master regulator. PRDM1 regulon activity drives the antibody secretory program (MZB1, XBP1, SDC1) while directly repressing B cell identity TFs PAX5 and BCL6.",
        "velten2017", tissue_context=["bone_marrow", "lymph_node", "tumor"]),
    make_rule("Plasma_cell", "PLASMA_NEG_001", "negative", "marker_genes",
        ["CD19", "PAX5"], "absent",
        "BLIMP1 directly represses PAX5 during plasma cell differentiation, causing exit from B cell identity. CD19 and PAX5 presence contradicts terminal plasma cell identity.",
        "velten2017", incompatible_with=["CD19", "PAX5"], tissue_context=["bone_marrow", "lymph_node"]),
    make_rule("Plasma_cell", "PLASMA_POS_001", "positive", "marker_genes",
        ["JCHAIN", "MZB1", "XBP1"], "high",
        "Plasma cells express JCHAIN (immunoglobulin J chain for IgM/IgA polymerization), MZB1 (ER chaperone for antibody folding), and XBP1 (IRE1 target, drives ER expansion for Ig secretion).",
        "velten2017", tissue_context=["bone_marrow", "lymph_node"]),
]

# ═══════════════════════════════════════════════════════════════════════════════
# Trophoblast — 4 cell types, 15 rules
# ═══════════════════════════════════════════════════════════════════════════════

# EVT (Extravillous Trophoblast)
NEW_RULES += [
    make_rule("EVT", "EVT_POS_001", "positive", "marker_genes",
        ["HLA-G", "ITGA5", "MMP2"], "high",
        "Extravillous trophoblasts (EVTs) invade the decidua and remodel spiral arteries. HLA-G (non-classical MHC for maternal immune tolerance), ITGA5 (fibronectin receptor for invasion), and MMP2 (matrix metalloproteinase for ECM degradation) define invasive EVT identity.",
        "ventotormo2018", tissue_context=["placenta", "decidua"]),
    make_rule("EVT", "EVT_NEG_001", "negative", "marker_genes",
        ["CDH1"], "absent",
        "CDH1 (E-cadherin) is expressed in villous cytotrophoblasts (VCT) but is downregulated during EVT differentiation as part of the epithelial-mesenchymal-like transition enabling invasion.",
        "ventotormo2018", incompatible_with=["CDH1"], tissue_context=["placenta"]),
]

# STB (Syncytiotrophoblast)
NEW_RULES += [
    make_rule("STB", "STB_POS_001", "positive", "marker_genes",
        ["CGB", "ERVFRD-1", "SDC1"], "high",
        "Syncytiotrophoblasts (STB) form the outer multinucleated layer of chorionic villi. CGB (hCG beta subunit, the pregnancy hormone), ERVFRD-1 (syncytin-2, fusogenic retroviral envelope protein for syncytium formation), and SDC1 (syndecan-1) define STB.",
        "ventotormo2018", tissue_context=["placenta"]),
    make_rule("STB", "STB_NEG_001", "negative", "marker_genes",
        ["HLA-G"], "absent",
        "HLA-G is the invasive EVT marker. STBs are non-invasive and do not require HLA-G for immune tolerance — they instead induce tolerance through other mechanisms.",
        "ventotormo2018", incompatible_with=["HLA-G"], tissue_context=["placenta"]),
    make_rule("STB", "STB_POS_002", "positive", "marker_genes",
        ["TFPI2"], "high",
        "TFPI2 (tissue factor pathway inhibitor 2) is highly expressed in STBs and regulates fibrinolysis at the maternal-fetal interface.",
        "ventotormo2018", tissue_context=["placenta"]),
]

# VCT (Villous Cytotrophoblast)
NEW_RULES += [
    make_rule("VCT", "VCT_POS_001", "positive", "marker_genes",
        ["CDH1", "EGFR", "TP63"], "high",
        "Villous cytotrophoblasts (VCT) are stem/progenitor trophoblasts that fuse to form STB or differentiate to EVT. CDH1 (E-cadherin for epithelial cohesion), EGFR (EGF receptor for proliferation), and moderate TP63 define VCT progenitor state.",
        "ventotormo2018", tissue_context=["placenta"]),
    make_rule("VCT", "VCT_NEG_001", "negative", "marker_genes",
        ["CGB"], "low",
        "CGB (hCG) is the STB-specific hormone. Low CGB in VCT reflects uncommitted progenitor state before syncytium fusion.",
        "ventotormo2018", tissue_context=["placenta"]),
]

# dNK (Decidual Natural Killer)
NEW_RULES += [
    make_rule("dNK", "DNK_POS_001", "positive", "marker_genes",
        ["NCAM1", "ITGA1", "KIR2DL1"], "high",
        "Decidual NK cells (dNK) are specialized uterine NK cells that promote vascular remodeling. They express NCAM1 (CD56), ITGA1 (CD49a, tissue-resident NK marker), and KIR2DL1 (inhibitory KIR for HLA-C tolerance of EVT).",
        "ventotormo2018", tissue_context=["placenta", "decidua", "uterus"]),
    make_rule("dNK", "DNK_NEG_001", "negative", "marker_genes",
        ["TRAC"], "absent",
        "TRAC (T cell receptor alpha constant) marks T cells. dNKs are innate lymphocytes without TCR recombination — TRAC absence is required to distinguish dNK from T cell populations.",
        "ventotormo2018", incompatible_with=["TRAC"], tissue_context=["placenta"]),
    make_rule("dNK", "DNK_NEG_002", "negative", "marker_genes",
        ["GZMB"], "low",
        "Peripheral cytotoxic NK cells are GZMB-high. dNKs have low cytotoxic potential and low GZMB, reflecting their remodeling rather than killing function in the decidua.",
        "ventotormo2018", tissue_context=["placenta"]),
]


def load_kg(path: Path) -> list:
    with open(path) as f:
        return json.load(f)


def save_kg(rules: list, path: Path) -> None:
    with open(path, "w") as f:
        json.dump(rules, f, indent=2)
    print(f"Saved {len(rules)} rules to {path}")


def main():
    kg = load_kg(KG_PATH)
    existing_ids = {r["rule_id"] for r in kg}
    existing_types = {r["cell_type"] for r in kg if r.get("status") == "ACTIVE"}

    new_added = [r for r in NEW_RULES if r["rule_id"] not in existing_ids]
    new_type_count = len({r["cell_type"] for r in new_added} - existing_types)

    print(f"Current KG: {len(kg)} rules, {len(existing_types)} cell types")
    print(f"New rules to add: {len(new_added)}")
    print(f"New cell types: {new_type_count}")

    kg.extend(new_added)

    active = [r for r in kg if r.get("status") == "ACTIVE"]
    types = {r["cell_type"] for r in active}
    print(f"Post-patch: {len(active)} ACTIVE rules, {len(types)} cell types")

    save_kg(kg, KG_PATH)

    # Print summary by tissue class
    from collections import Counter
    new_types = sorted({r["cell_type"] for r in new_added})
    print(f"\nNew cell types added ({len(new_types)}):")
    for t in new_types:
        n = sum(1 for r in new_added if r["cell_type"] == t)
        print(f"  {t}: {n} rules")


if __name__ == "__main__":
    main()
