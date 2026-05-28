"""
Patch A supplement — adds ~30 more cell types (2-3 rules each) to reach
≥150 cell types / ≥450 ACTIVE rules.

Cell types covered: adipose, skeletal muscle, additional CNS subtypes,
CLP, adrenal cortex/medulla, hepatic stellate, additional endothelial.

Run: /Users/sivamoturi/miniconda/bin/python3.10 scripts/patch_a_supplement.py
"""
import json
from pathlib import Path

ROOT = Path(__file__).parent.parent
KG_PATH = ROOT / "kg_data" / "oracle_kg_v0.2.0.json"

PMIDS = {
    "tabumuris":   "30283141",  # Tabula Muris
    "pellin2019":  "29588278",  # Single-cell hematopoietic landscape
    "velten2017":  "28319093",  # Human HSC lineage commitment
    "vanland2018": "29443965",  # Brain vasculature molecular atlas
    "litvin2020":  "32971526",  # Cells of adult human heart
    "baron2016":   "27667365",  # Pancreas single cell
    "lake2023":    "37468583",  # Kidney KPMP atlas
    "ji2020":      "32579974",  # Skin multimodal
    "marques2016": "27284195",  # Oligodendrocyte heterogeneity
    "allenbrain":  "38168182",  # Allen Brain 2023
    "elyada2019":  "31197017",  # iCAF/myCAF Elyada 2019
    "zheng2021":   "34914499",  # Pan-cancer TIL
    "ventotormo":  "30429548",  # Trophoblast
    "erythroid":   "32457162",  # Erythroid regulators
    "litvin2020":  "32971526",  # Heart
}


def r(ct, rid, rt, es, genes, direction, basis, pmid_key,
       confidence="high", tissue=None, incompatible=None):
    rule = {
        "cell_type": ct, "rule_id": rid, "rule_type": rt,
        "evidence_source": es, "gene_or_regulon": genes,
        "direction": direction, "mechanistic_basis": basis,
        "pmid": PMIDS[pmid_key], "confidence": confidence,
        "tissue_context": tissue or [], "source_db": "primary_literature",
        "status": "ACTIVE", "added_in_version": "0.2.0-patch-a",
    }
    if incompatible:
        rule["incompatible_with"] = incompatible
    return rule


NEW = [
    # ── Adipose ──────────────────────────────────────────────────────────────
    r("Adipocyte", "ADIPO_POS_001", "positive", "marker_genes",
      ["PPARG", "ADIPOQ", "FABP4"], "high",
      "Mature adipocytes express PPARG (master adipogenic TF), ADIPOQ (adiponectin, adipokine secreted by adipocytes), and FABP4 (fatty acid binding protein 4, lipid transport). This trio defines the lipid-storing adipocyte program.",
      "tabumuris", tissue=["adipose"]),
    r("Adipocyte", "ADIPO_POS_002", "positive", "marker_genes",
      ["LPL"], "high",
      "LPL (lipoprotein lipase) is expressed at high levels in adipocytes for triglyceride hydrolysis and fatty acid uptake from circulating lipoproteins.",
      "tabumuris", tissue=["adipose"]),

    r("Pre_Adipocyte", "PREADIPO_POS_001", "positive", "marker_genes",
      ["PDGFRA", "PREF1"], "high",
      "Pre-adipocytes are adipogenic progenitors expressing PDGFRA (fibroblast growth factor receptor) and PREF1/DLK1 (preadipocyte factor-1, inhibitor of adipogenesis that is lost upon terminal differentiation).",
      "tabumuris", tissue=["adipose"]),
    r("Pre_Adipocyte", "PREADIPO_NEG_001", "negative", "marker_genes",
      ["ADIPOQ"], "absent",
      "ADIPOQ (adiponectin) is a terminal adipocyte marker absent in pre-adipocytes. Pre-adipocyte identity requires PREF1 high and ADIPOQ absent.",
      "tabumuris", incompatible=["ADIPOQ"], tissue=["adipose"]),

    r("Adipose_Macrophage", "ADIPMAC_POS_001", "positive", "marker_genes",
      ["CSF1R", "MRC1", "TREM2"], "high",
      "Adipose tissue macrophages express CSF1R (M-CSF receptor for tissue-resident macrophage maintenance), MRC1 (CD206, alternatively activated), and TREM2 (lipid-associated macrophage marker in metabolically stressed adipose).",
      "tabumuris", tissue=["adipose"]),

    # ── Skeletal muscle ───────────────────────────────────────────────────────
    r("Skeletal_Muscle_Fiber", "SKEL_POS_001", "positive", "marker_genes",
      ["MYH1", "ACTN2", "TNNT3"], "high",
      "Skeletal muscle fibers express fast-twitch myosin heavy chain MYH1, alpha-actinin ACTN2 (sarcomere cross-linking), and fast-type troponin T (TNNT3). This defines the mature skeletal myofiber program.",
      "tabumuris", tissue=["skeletal_muscle"]),
    r("Skeletal_Muscle_Fiber", "SKEL_NEG_001", "negative", "marker_genes",
      ["MYH7"], "low",
      "MYH7 (beta-myosin) dominates in cardiac muscle and slow-twitch skeletal fibers. Most fast-twitch skeletal fibers have low MYH7, distinguishing them from cardiac cardiomyocytes.",
      "litvin2020", tissue=["skeletal_muscle"]),

    r("Satellite_cell", "SAT_POS_001", "positive", "marker_genes",
      ["PAX7", "CD34", "NCAM1"], "high",
      "Muscle satellite cells are muscle stem cells expressing PAX7 (muscle stem cell master TF), CD34 (progenitor marker), and NCAM1 (CD56, marks quiescent satellite cells).",
      "tabumuris", tissue=["skeletal_muscle"]),
    r("Satellite_cell", "SAT_NEG_001", "negative", "marker_genes",
      ["MYH1"], "absent",
      "MYH1 marks terminally differentiated myofibers. Satellite cells are undifferentiated progenitors and do not express structural myosin.",
      "tabumuris", incompatible=["MYH1"], tissue=["skeletal_muscle"]),

    r("Tenocyte", "TENO_POS_001", "positive", "marker_genes",
      ["COL1A1", "TNMD", "SCX"], "high",
      "Tenocytes are the primary cells of tendons and ligaments, expressing high COL1A1 (type I collagen), TNMD (tenomodulin, tenocyte-specific glycoprotein), and SCX (scleraxis, tendon TF).",
      "tabumuris", tissue=["tendon"]),

    # ── Additional CNS subtypes ───────────────────────────────────────────────
    r("Dopaminergic_Neuron", "DA_POS_001", "positive", "marker_genes",
      ["TH", "SLC6A3", "NR4A2"], "high",
      "Dopaminergic neurons express TH (tyrosine hydroxylase, rate-limiting enzyme for dopamine synthesis), SLC6A3 (DAT, dopamine transporter for reuptake), and NR4A2/Nurr1 (TF required for DA neuron differentiation and maintenance).",
      "allenbrain", tissue=["brain", "substantia_nigra", "VTA"]),
    r("Dopaminergic_Neuron", "DA_NEG_001", "negative", "marker_genes",
      ["SLC17A7"], "absent",
      "SLC17A7 (VGLUT1) marks glutamatergic neurons. Dopaminergic neurons are not glutamatergic, though some co-release glutamate — VGLUT1 high is incompatible with canonical DA neuron identity.",
      "allenbrain", incompatible=["SLC17A7"], tissue=["brain"]),

    r("Cholinergic_Neuron", "CHOL_POS_001", "positive", "marker_genes",
      ["CHAT", "SLC18A3"], "high",
      "Cholinergic neurons express ChAT (choline acetyltransferase, synthesizes acetylcholine) and SLC18A3 (VAChT, vesicular acetylcholine transporter). This pair defines all cholinergic identity.",
      "allenbrain", tissue=["brain", "spinal_cord", "basal_forebrain"]),

    r("Purkinje_cell", "PRKJ_POS_001", "positive", "marker_genes",
      ["CALB1", "PCP2", "GRID2"], "high",
      "Cerebellar Purkinje cells express CALB1 (calbindin D28k, calcium buffer), PCP2 (L7, Purkinje cell-specific protein), and GRID2 (GluRdelta2, cerebellar synapse receptor). This trio is specific to Purkinje cells.",
      "allenbrain", tissue=["brain", "cerebellum"]),
    r("Purkinje_cell", "PRKJ_NEG_001", "negative", "marker_genes",
      ["GAD1"], "high",
      "Purkinje cells are GABAergic (GAD1+) — this confirms inhibitory identity — but their defining markers (CALB1, PCP2) must also be present. GAD1 alone does not distinguish Purkinje from other inhibitory neurons.",
      "allenbrain", tissue=["brain", "cerebellum"]),

    r("Bergmann_glia", "BERG_POS_001", "positive", "marker_genes",
      ["GFAP", "BLBP", "PTPRZ1"], "high",
      "Bergmann glia are specialized radial glia of the cerebellar cortex. They express GFAP (astrocytic marker), BLBP/FABP7 (brain lipid binding protein), and PTPRZ1 (phosphacan, radial glia marker).",
      "allenbrain", tissue=["brain", "cerebellum"]),

    # ── Adrenal cortex/medulla ────────────────────────────────────────────────
    r("Adrenal_Cortex_ZG", "ADCORTEX_ZG_POS_001", "positive", "marker_genes",
      ["CYP11B2", "NR5A1"], "high",
      "Zona glomerulosa (ZG) cells produce aldosterone via CYP11B2 (aldosterone synthase, only expressed in ZG). NR5A1 (SF-1/Ad4BP) is the adrenal cortex master TF required in all cortical zones.",
      "tabumuris", tissue=["adrenal_gland"]),
    r("Adrenal_Cortex_ZG", "ADCORTEX_ZG_NEG_001", "negative", "marker_genes",
      ["CYP11B1"], "low",
      "CYP11B1 (11-beta-hydroxylase) is expressed in the zona fasciculata for cortisol synthesis. Low CYP11B1 distinguishes ZG (aldosterone) from ZF (cortisol) cells.",
      "tabumuris", tissue=["adrenal_gland"]),

    r("Adrenal_Cortex_ZF", "ADCORTEX_ZF_POS_001", "positive", "marker_genes",
      ["CYP11B1", "CYP17A1", "NR5A1"], "high",
      "Zona fasciculata (ZF) cells produce cortisol via CYP11B1 (11-beta-hydroxylase) and CYP17A1 (17-hydroxylase). Both enzymes are required for glucocorticoid synthesis absent from ZG.",
      "tabumuris", tissue=["adrenal_gland"]),

    r("Adrenal_Chromaffin", "ADCHROM_POS_001", "positive", "marker_genes",
      ["CHGA", "TH", "PNMT"], "high",
      "Adrenal chromaffin cells are neuroendocrine cells that produce epinephrine (adrenaline). CHGA (chromogranin A, neuroendocrine vesicle marker), TH (dopamine/norepinephrine synthesis), and PNMT (phenylethanolamine N-methyltransferase, epinephrine synthesis) define them.",
      "tabumuris", tissue=["adrenal_gland", "adrenal_medulla"]),
    r("Adrenal_Chromaffin", "ADCHROM_NEG_001", "negative", "marker_genes",
      ["NR5A1"], "absent",
      "NR5A1 (SF-1) is the adrenal cortex master TF. Chromaffin cells are neural crest-derived (medulla) and do not express the cortical TF NR5A1.",
      "tabumuris", incompatible=["NR5A1"], tissue=["adrenal_gland"]),

    # ── Liver additions ───────────────────────────────────────────────────────
    r("Hepatic_Stellate_cell", "HSC_HPST_POS_001", "positive", "marker_genes",
      ["ACTA2", "PDGFRB", "COL1A1", "LRAT"], "high",
      "Hepatic stellate cells (HSCs) are vitamin A-storing perisinusoidal cells that become activated myofibroblasts in fibrosis. ACTA2 (alpha-SMA), PDGFRB (PDGF receptor), COL1A1, and LRAT (lecithin-retinol acyltransferase, retinoid storage) define activated HSCs.",
      "tabumuris", tissue=["liver"]),
    r("Hepatic_Stellate_cell", "HSC_HPST_NEG_001", "negative", "marker_genes",
      ["ALB"], "absent",
      "ALB (albumin) is the hallmark hepatocyte marker. Hepatic stellate cells are perisinusoidal non-parenchymal cells and do not express albumin.",
      "tabumuris", incompatible=["ALB"], tissue=["liver"]),

    r("Cholangiocyte", "CHOL_BILE_POS_001", "positive", "marker_genes",
      ["KRT7", "KRT19", "CFTR", "SOX9"], "high",
      "Cholangiocytes are the biliary epithelial cells lining bile ducts. KRT7 and KRT19 (biliary keratins), CFTR (bicarbonate secretion into bile), and SOX9 (biliary/ductal TF) define cholangiocyte identity.",
      "tabumuris", tissue=["liver", "biliary"]),
    r("Cholangiocyte", "CHOL_BILE_NEG_001", "negative", "marker_genes",
      ["ALB"], "absent",
      "ALB marks hepatocytes. Cholangiocytes are biliary epithelial cells and are mutually exclusive with hepatocyte identity.",
      "tabumuris", incompatible=["ALB"], tissue=["liver"]),

    # ── Additional immune progenitors ─────────────────────────────────────────
    r("CLP", "CLP_POS_001", "positive", "marker_genes",
      ["IL7R", "CD34", "LTB"], "high",
      "Common lymphoid progenitors (CLPs) express IL7R (IL-7 receptor, lymphoid lineage signal), CD34 (progenitor marker), and LTB (lymphotoxin beta). IL7R acquisition marks lymphoid commitment from multipotent HSPCs.",
      "velten2017", tissue=["bone_marrow"]),
    r("CLP", "CLP_NEG_001", "negative", "marker_genes",
      ["MPO"], "absent",
      "MPO marks granulocyte/myeloid commitment. CLPs have exited the myeloid path and do not express MPO.",
      "pellin2019", incompatible=["MPO"], tissue=["bone_marrow"]),

    r("ILC_precursor", "ILCP_POS_001", "positive", "marker_genes",
      ["PLZF", "RORC", "IL7R"], "high",
      "ILC precursors (ILCPs) express PLZF/ZBTB16 (ILC commitment TF), RORC (lineage priming for ILC3), and IL7R. ILCPs are committed lymphoid progenitors that have lost NK/T potential but retain ILC1-3 potential.",
      "pellin2019", tissue=["bone_marrow", "tonsil"]),

    # ── Additional epithelial ─────────────────────────────────────────────────
    r("Goblet_cell", "GOBLET_POS_001", "positive", "marker_genes",
      ["MUC2", "CLCA1", "FCGBP"], "high",
      "Goblet cells are mucin-secreting cells of intestinal and airway epithelium. MUC2 (secreted mucin, the major intestinal gel), CLCA1 (calcium-activated chloride channel for mucus secretion), and FCGBP (Fc-gamma binding protein, mucus-associated) define them.",
      "tabumuris", tissue=["intestine", "colon", "lung"]),
    r("Goblet_cell", "GOBLET_NEG_001", "negative", "marker_genes",
      ["CHGA"], "absent",
      "CHGA marks enteroendocrine cells. Goblet cells are secretory non-endocrine cells; CHGA presence would indicate enteroendocrine, not goblet, identity.",
      "tabumuris", incompatible=["CHGA"], tissue=["intestine"]),

    r("Enteroendocrine_cell", "ENTERO_POS_001", "positive", "marker_genes",
      ["CHGA", "CHGB", "NEUROD1"], "high",
      "Enteroendocrine cells are hormone-secreting epithelial cells. CHGA and CHGB (chromogranins, neuroendocrine vesicle proteins) and NEUROD1 (neuroendocrine TF) define their neuroendocrine identity.",
      "tabumuris", tissue=["intestine", "colon"]),

    r("Paneth_cell", "PANETH_POS_001", "positive", "marker_genes",
      ["DEFA5", "DEFA6", "LYZ"], "high",
      "Paneth cells are specialized secretory cells at intestinal crypt bases. DEFA5, DEFA6 (defensins, antimicrobial peptides), and LYZ (lysozyme) define the Paneth cell antimicrobial secretory program.",
      "tabumuris", tissue=["intestine", "small_intestine"]),

    r("Tuft_cell", "TUFT_POS_001", "positive", "marker_genes",
      ["DCLK1", "POU2F3", "IL25"], "high",
      "Tuft cells (brush cells) are chemosensory epithelial cells defined by DCLK1 (doublecortin kinase 1), POU2F3 (master TF for tuft cell differentiation), and IL25 (IL-17E, activates ILC2 for type 2 immunity initiation).",
      "tabumuris", tissue=["intestine", "lung"]),

    # ── Additional stromal ────────────────────────────────────────────────────
    r("Lymph_node_FDC", "FDC_POS_001", "positive", "marker_genes",
      ["CR1", "CR2", "CXCL13"], "high",
      "Follicular dendritic cells (FDCs) are non-hematopoietic stromal cells that display antigen in germinal centers. CR1 (CD35, complement receptor for antigen trapping), CR2 (CD21), and CXCL13 (B cell chemoattractant maintaining germinal center) define FDCs.",
      "tabumuris", tissue=["lymph_node", "germinal_center"]),
    r("Lymph_node_FDC", "FDC_NEG_001", "negative", "marker_genes",
      ["CD45"], "absent",
      "CD45 (PTPRC) marks all hematopoietic cells. FDCs are stromal cells of mesenchymal origin and do not express CD45, distinguishing them from macrophages, DCs, and other immune cells in the germinal center.",
      "tabumuris", incompatible=["PTPRC"], tissue=["lymph_node"]),

    r("Peritoneal_Meso", "MESO_POS_001", "positive", "marker_genes",
      ["MSLN", "CALB2", "WT1"], "high",
      "Mesothelial cells line the peritoneal, pleural, and pericardial cavities. MSLN (mesothelin, serosal surface glycoprotein), CALB2 (calretinin), and WT1 (Wilms tumor 1, mesothelial TF) define mesothelial identity.",
      "tabumuris", tissue=["peritoneum", "pleura", "pericardium"]),

    # ── Cornea/eye ────────────────────────────────────────────────────────────
    r("Corneal_Endothelial", "CORNEND_POS_001", "positive", "marker_genes",
      ["SLC4A11", "AQP1", "COL8A1"], "high",
      "Corneal endothelial cells pump fluid to maintain corneal transparency. SLC4A11 (bicarbonate transporter, mutated in corneal dystrophy), AQP1 (aquaporin water channel), and COL8A1 (collagen VIII for Descemet membrane) define them.",
      "tabumuris", tissue=["eye", "cornea"]),

    # ── Reproductive ──────────────────────────────────────────────────────────
    r("Sertoli_cell", "SERT_POS_001", "positive", "marker_genes",
      ["SOX9", "AMH", "FSHR"], "high",
      "Sertoli cells are the nursing cells of the testis that support spermatogenesis. SOX9 (Sertoli master TF), AMH (anti-Müllerian hormone, produced by Sertoli cells), and FSHR (FSH receptor) define testicular Sertoli cells.",
      "tabumuris", tissue=["testis"]),
    r("Leydig_cell", "LEYDIG_POS_001", "positive", "marker_genes",
      ["CYP17A1", "STAR", "HSD3B1"], "high",
      "Leydig cells produce testosterone in the testicular interstitium. CYP17A1 (steroidogenic enzyme for androgen synthesis), STAR (steroidogenic acute regulatory protein), and HSD3B1 (3-beta-HSD for progesterone/androgen synthesis) define them.",
      "tabumuris", tissue=["testis"]),

    # ── Kidney additional ─────────────────────────────────────────────────────
    r("Macula_Densa", "MD_POS_001", "positive", "marker_genes",
      ["PTGS2", "NOS1"], "high",
      "Macula densa cells are specialized tubular cells at the thick ascending limb / DCT junction that sense luminal [NaCl]. PTGS2 (COX-2, prostaglandin synthesis for renin regulation) and NOS1 (neuronal NOS, nitric oxide for afferent arteriole tone) are induced by low NaCl.",
      "lake2023", tissue=["kidney", "macula_densa"]),
]


def main():
    with open(KG_PATH) as f:
        kg = json.load(f)

    existing_ids = {r["rule_id"] for r in kg}
    to_add = [rule for rule in NEW if rule["rule_id"] not in existing_ids]

    print(f"Before: {len(kg)} rules")
    kg.extend(to_add)
    active = [r for r in kg if r.get("status") == "ACTIVE"]
    types = {r["cell_type"] for r in active}
    print(f"Added {len(to_add)} rules")
    print(f"After: {len(active)} ACTIVE rules, {len(types)} cell types")

    with open(KG_PATH, "w") as f:
        json.dump(kg, f, indent=2)
    print(f"Saved to {KG_PATH}")


if __name__ == "__main__":
    main()
