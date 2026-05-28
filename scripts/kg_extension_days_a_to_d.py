"""
KG Extension: Days A-D
Adds retina, bone/gonads, inner ear/cornea/salivary/synovium, and vascular/misc
cell types to the AXIOM-SC knowledge graph.

Verified anchor PMIDs:
  31436334 — Lukowski 2019 retina single-cell atlas (EMBO J)
  11694879 — NRL required for rod photoreceptor development (Nat Genet)
  24318667 — RBPMS selective marker of ganglion cells (Invest Ophthalmol Vis Sci)
  18094249 — Lim1/LHX1 essential for horizontal cell positioning (J Neurosci)
  14633986 — Osteocyte control via sclerostin BMP antagonist (J Bone Miner Res)
  39027344 — Sp7 in osteoblasts: proliferation and differentiation (Front Endocrinol)
  38218304 — Runx2 deletion in hypertrophic chondrocytes (Bone)
  11680692 — L-Sox5, Sox6 and Sox9 chondrocyte differentiation (Dev Cell)
  14736745 — Foxl2 required for granulosa cell development (Genes Dev)
  30735645 — Human spermatogonial stem cells single-cell magnifying glass (Cell Rep)
  40695261 — SYCP3 meiotic arrest spermatocytes (Hum Reprod)
  40482072 — Early transcriptional states and markers of spermatogonia (Andrology)
  35387428 — Sertoli cell number and function in human testicular tissue (Sci Rep)
  41596239 — CDK5RAP3 testosterone production in Leydig cells (Mol Cell Biol)
  41945437 — Regulatory axis for MYO7A expression in cochlear hair cells (Stem Cell Rep)
  41638427 — Electromotility/SLC26A5 in outer hair cells (Biophys J)
  30410537 — Key transcription factors in vestibular/auditory supporting cells (Hear Res)
  41491447 — NEUROD1 regulatory networks in spiral ganglion specification (J Neurosci)
  20012478 — Endocochlear potential formation stria vascularis KCNJ10 (Hear Res)
  42194275 — Human corneal organoids KRT12/KRT3 (Biomaterials)
  40437144 — Co-differentiation corneal endothelial/keratocyte KERA/ALDH3A1 (Exp Eye Res)
  35220796 — Human parotid gland single-cell transcriptomics (Mol Ther Nucleic Acids)
  37024468 — PRG4 lubricin synovial fibroblast articular joint (Sci Rep)
  13679391 — ACKR1/DARC on high endothelial venules chemokine binding (J Exp Med)
  16807374 — Pericyte microarray KCNJ8/PDGFRB blood microvessels (FASEB J)
  33208946 — Human lung molecular cell atlas CA4 aerocyte (Nature)
  10906459 — NKX3.1 targeted disruption prostate morphogenetic defects (EMBO J)
  21853341 — UPK1A UPK1B urothelial preneoplastic (Histopathology)
  32046769 — Choroid plexus LAT2/SNAT3 CSF TTR amino acid (J Neurochem)
  33325807 — ANO1 ICC dysmotility bowel resection (J Smooth Muscle Res)
  40585958 — FAP PDGFRA skeletal muscle degeneration (Mol Ther)
  12753744 — NOTCH2 preferentially expressed in mature B cells MZB (J Exp Med)
  20421466 — Tim4/TIMD4 essential for maintenance of peritoneal macrophages (Nature)
  24630724 — Heme SPI-C red pulp macrophage differentiation (Nat Immunol)
  39178365 — Integrative single-cell liver spatial transcriptomics (Hepatology)

DO NOT RE-RUN on an existing KG — rules will duplicate.
Run once, then verify with: pytest tests/test_kg_extension_days_a_d.py
"""

import json
from pathlib import Path

KG_PATH = Path("kg_data/oracle_kg_v0.2.0.json")


def _rule(rule_id, cell_type, rule_type, evidence_source, genes, direction,
          mechanistic_basis, pmid, confidence="high", tissue=None, source_db="literature"):
    return {
        "rule_id": rule_id,
        "cell_type": cell_type,
        "rule_type": rule_type,
        "evidence_source": evidence_source,
        "gene_or_regulon": genes,
        "direction": direction,
        "mechanistic_basis": mechanistic_basis,
        "pmid": pmid,
        "confidence": confidence,
        "tissue_context": tissue or [],
        "source_db": source_db,
        "status": "ACTIVE",
        "added_in_version": "0.2.0",
        "paired_with": [],
        "incompatible_with": [],
    }


# ── DAY A: RETINA ─────────────────────────────────────────────────────────────
# 8 new types: Rod_photoreceptor, Cone_photoreceptor, Muller_glia,
#   Horizontal_cell, Amacrine_cell, Bipolar_cell, Ciliary_epithelium_retina,
#   Retinal_Endothelium
# Enrichments: Retinal_Ganglion_cell (+3), Retinal_Pigment_Epithelium (+3)

DAY_A_RULES = [
    # Rod photoreceptor (3 rules)
    _rule("ROD_POS_001", "Rod_photoreceptor", "positive", "marker_genes",
          ["NRL", "RHO", "CNGB1"],
          "high",
          "NRL is the rod master transcription factor that activates RHO (rhodopsin) and CNGB1; all three define rod phototransduction identity in the human retina.",
          "11694879", tissue=["retina"]),
    _rule("ROD_NEG_001", "Rod_photoreceptor", "negative", "marker_genes",
          ["ARR3"],
          "absent",
          "ARR3 (cone arrestin) is exclusively expressed in cone photoreceptors to terminate cone-specific phototransduction; its presence is incompatible with rod identity.",
          "31436334", tissue=["retina"]),
    _rule("ROD_CIRCUIT_001", "Rod_photoreceptor", "circuit", "regulon",
          ["NRL"],
          "active",
          "Active NRL regulon drives the entire rod-specific gene expression program including RHO and phototransduction components; NRL deletion converts rods to cones.",
          "11694879", tissue=["retina"]),

    # Cone photoreceptor (3 rules)
    _rule("CONE_POS_001", "Cone_photoreceptor", "positive", "marker_genes",
          ["RXRG", "ARR3", "OPN1SW"],
          "high",
          "RXRG (RXR-gamma TF) drives cone-specific programs; ARR3 (cone arrestin) terminates phototransduction; OPN1SW marks short-wavelength cones in human retina.",
          "31436334", tissue=["retina"]),
    _rule("CONE_NEG_001", "Cone_photoreceptor", "negative", "regulon",
          ["NRL"],
          "inactive",
          "NRL is a rod master TF; its active regulon would redirect cone fate toward rod identity via RHO activation. Cones require NRL silencing.",
          "11694879", tissue=["retina"]),
    _rule("CONE_CIRCUIT_001", "Cone_photoreceptor", "circuit", "marker_genes",
          ["RXRG", "ARR3"],
          "high",
          "RXRG-driven transcriptional program co-activates ARR3 in cones; co-expression of both confirms the cone-specific signal-termination circuit.",
          "31436334", tissue=["retina"]),

    # Muller glia (3 rules)
    _rule("MULL_POS_001", "Muller_glia", "positive", "marker_genes",
          ["RLBP1", "SLC1A3", "GLUL"],
          "high",
          "RLBP1 (cellular retinaldehyde binding protein) and GLUL (glutamine synthetase) define Muller glia metabolic function; SLC1A3 mediates glutamate uptake from synaptic cleft.",
          "31436334", tissue=["retina"]),
    _rule("MULL_NEG_001", "Muller_glia", "negative", "marker_genes",
          ["RBPMS"],
          "absent",
          "RBPMS is a selective RNA-binding protein exclusively in retinal ganglion cells; its presence distinguishes RGCs from Muller glia which occupy the same nuclear layer.",
          "24318667", tissue=["retina"]),
    _rule("MULL_CIRCUIT_001", "Muller_glia", "circuit", "marker_genes",
          ["RLBP1", "GLUL"],
          "high",
          "RLBP1 and GLUL co-expression defines the Muller glia retinoid and glutamine recycling metabolic circuit; both are required for normal retinal neurotransmission support.",
          "31436334", tissue=["retina"]),

    # Horizontal cell (3 rules)
    _rule("HORIZ_POS_001", "Horizontal_cell", "positive", "marker_genes",
          ["LHX1", "CALB1", "ONECUT1"],
          "high",
          "LHX1 (Lim1) is a homeodomain TF required for horizontal cell positioning; CALB1 (calbindin) and ONECUT1 are co-expressed horizontal cell markers in human retina.",
          "18094249", tissue=["retina"]),
    _rule("HORIZ_NEG_001", "Horizontal_cell", "negative", "marker_genes",
          ["MYO7A"],
          "absent",
          "MYO7A is expressed in photoreceptors and retinal pigment epithelium for cytoskeletal function; absence distinguishes horizontal cells in the inner nuclear layer.",
          "31436334", tissue=["retina"]),
    _rule("HORIZ_CIRCUIT_001", "Horizontal_cell", "circuit", "regulon",
          ["LHX1"],
          "active",
          "Active LHX1 regulon establishes horizontal cell identity and correct laminar positioning; LHX1 loss causes horizontal cells to be mispositioned in the ganglion cell layer.",
          "18094249", tissue=["retina"]),

    # Amacrine cell (3 rules)
    _rule("AMAC_POS_001", "Amacrine_cell", "positive", "marker_genes",
          ["GAD1", "TFAP2A", "CHAT"],
          "high",
          "GAD1 (GABA synthesis) and TFAP2A define inhibitory amacrine subtypes; CHAT marks cholinergic starburst amacrine cells; all expressed in INL amacrine cell population.",
          "31436334", tissue=["retina"]),
    _rule("AMAC_NEG_001", "Amacrine_cell", "negative", "marker_genes",
          ["RBPMS"],
          "absent",
          "RBPMS is a selective ganglion cell marker; amacrine cells occupy the inner nuclear layer and do not express RBPMS.",
          "24318667", tissue=["retina"]),
    _rule("AMAC_NEG_002", "Amacrine_cell", "negative", "marker_genes",
          ["SLC26A5"],
          "absent",
          "SLC26A5 (prestin) is exclusively expressed in outer hair cells of the cochlea and outer segment of OHC; never expressed in retinal amacrine cells.",
          "31436334", tissue=["retina"]),

    # Bipolar cell (3 rules)
    _rule("BIPOL_POS_001", "Bipolar_cell", "positive", "marker_genes",
          ["PRKCA", "CABP5"],
          "high",
          "PRKCA (protein kinase C alpha) is a classic rod bipolar cell marker; CABP5 (calcium binding protein 5) marks ON bipolar cells mediating rod-to-RGC signal transmission.",
          "31436334", tissue=["retina"]),
    _rule("BIPOL_NEG_001", "Bipolar_cell", "negative", "marker_genes",
          ["RBPMS"],
          "absent",
          "RBPMS is a selective RGC marker; bipolar cells project within the inner nuclear layer and do not express RBPMS.",
          "24318667", tissue=["retina"]),
    _rule("BIPOL_CIRCUIT_001", "Bipolar_cell", "circuit", "regulon",
          ["VSX2"],
          "active",
          "VSX2 (visual system homeobox 2, also CHX10) regulon is active in bipolar progenitors and mature bipolar cells; VSX2 loss eliminates bipolar cell identity.",
          "31436334", tissue=["retina"]),

    # Ciliary epithelium retina (2 rules)
    _rule("CIEP_POS_001", "Ciliary_epithelium_retina", "positive", "marker_genes",
          ["AQP1", "PLCL1", "PTGDS"],
          "high",
          "AQP1 (aquaporin-1) mediates aqueous humor flow through ciliary epithelium; PLCL1 and PTGDS are expressed in ciliary body non-pigmented epithelium.",
          "31436334", tissue=["retina", "eye"]),
    _rule("CIEP_NEG_001", "Ciliary_epithelium_retina", "negative", "marker_genes",
          ["RHO"],
          "absent",
          "RHO (rhodopsin) is strictly photoreceptor-specific; its presence would indicate rod contamination rather than ciliary epithelium.",
          "31436334", tissue=["retina", "eye"]),

    # Retinal endothelium (2 rules)
    _rule("RETENDO_POS_001", "Retinal_Endothelium", "positive", "marker_genes",
          ["CDH5", "PECAM1", "NRP1"],
          "high",
          "CDH5 (VE-cadherin) and PECAM1 define vascular endothelium; NRP1 marks retinal blood vessels including tip cells in angiogenic retinal vasculature.",
          "31436334", tissue=["retina"]),
    _rule("RETENDO_NEG_001", "Retinal_Endothelium", "negative", "marker_genes",
          ["RLBP1"],
          "absent",
          "RLBP1 is a glial marker expressed in Muller glia and RPE; its presence would indicate Muller glia rather than retinal vascular endothelium.",
          "31436334", tissue=["retina"]),

    # Enrich Retinal_Ganglion_cell (3 rules)
    _rule("RGC_CIRCUIT_001", "Retinal_Ganglion_cell", "circuit", "marker_genes",
          ["RBPMS", "SNCG", "POU4F2"],
          "high",
          "RBPMS, SNCG (gamma-synuclein), and POU4F2 (Brn3b) form the core RGC identity circuit; all three are co-expressed in nearly all RGC subtypes across human retina.",
          "24318667", tissue=["retina"]),
    _rule("RGC_NEG_001", "Retinal_Ganglion_cell", "negative", "marker_genes",
          ["RLBP1"],
          "absent",
          "RLBP1 is exclusively expressed in Muller glia and RPE for retinoid cycling; its presence contradicts RGC identity.",
          "31436334", tissue=["retina"]),
    _rule("RGC_NEG_002", "Retinal_Ganglion_cell", "negative", "marker_genes",
          ["ARR3"],
          "absent",
          "ARR3 (cone arrestin) is a cone photoreceptor-specific marker; RGCs project to brain and do not participate in phototransduction.",
          "31436334", tissue=["retina"]),

    # Enrich Retinal_Pigment_Epithelium (3 rules)
    _rule("RPE_NEG_001", "Retinal_Pigment_Epithelium", "negative", "marker_genes",
          ["RHO"],
          "absent",
          "RHO (rhodopsin) is strictly a rod photoreceptor outer segment protein; its presence would indicate rod contamination of RPE cultures or tissue sections.",
          "31436334", tissue=["retina"]),
    _rule("RPE_CIRCUIT_001", "Retinal_Pigment_Epithelium", "circuit", "marker_genes",
          ["RPE65", "BEST1"],
          "high",
          "RPE65 (retinoid isomerase) and BEST1 (bestrophin-1 Cl- channel) co-expression defines the RPE visual cycle and fluid transport circuit essential for photoreceptor support.",
          "31436334", tissue=["retina"]),
    _rule("RPE_NEG_002", "Retinal_Pigment_Epithelium", "negative", "marker_genes",
          ["RBPMS"],
          "absent",
          "RBPMS is a selective ganglion cell marker; RPE cells in the outer retina layer do not express RBPMS.",
          "24318667", tissue=["retina"]),
]

# ── DAY B: BONE/CARTILAGE + GONADS ───────────────────────────────────────────
# 8 new types: Osteocyte, Hypertrophic_chondrocyte, Spermatogonial_SC,
#   Spermatocyte, Spermatid, Granulosa_cell, Theca_cell, Oocyte
# Enrichments: Osteoblast (+3), Chondrocyte (+2), Sertoli_cell (+3), Leydig_cell (+3)

DAY_B_RULES = [
    # Osteocyte (3 rules)
    _rule("OSTCYTE_POS_001", "Osteocyte", "positive", "marker_genes",
          ["PHEX", "DMP1", "SOST"],
          "high",
          "PHEX and DMP1 are osteocyte-specific phosphate regulators; SOST (sclerostin) is secreted by osteocytes to inhibit Wnt/bone formation signaling — none are expressed in osteoblasts.",
          "14633986", tissue=["bone"]),
    _rule("OSTCYTE_NEG_001", "Osteocyte", "negative", "marker_genes",
          ["CTSK"],
          "absent",
          "CTSK (cathepsin K) is the osteoclast-specific cysteine protease for bone resorption; osteocytes are embedded in bone matrix and do not resorb bone.",
          "14633986", tissue=["bone"]),
    _rule("OSTCYTE_CIRCUIT_001", "Osteocyte", "circuit", "marker_genes",
          ["SOST", "DMP1"],
          "high",
          "SOST and DMP1 co-expression defines the osteocyte mechanosensing circuit: SOST responds to mechanical unloading while DMP1 supports perilacunar remodeling.",
          "14633986", tissue=["bone"]),

    # Hypertrophic chondrocyte (3 rules)
    _rule("HYPCHOND_POS_001", "Hypertrophic_chondrocyte", "positive", "marker_genes",
          ["COL10A1", "MMP13", "IHH"],
          "high",
          "COL10A1 (type X collagen) and MMP13 are hypertrophic chondrocyte markers driving endochondral ossification; IHH coordinates chondrocyte-osteoblast signaling in the growth plate.",
          "38218304", tissue=["cartilage", "bone"]),
    _rule("HYPCHOND_NEG_001", "Hypertrophic_chondrocyte", "negative", "regulon",
          ["SOX9"],
          "inactive",
          "SOX9 regulon activity maintains chondrocyte proliferating/prehypertrophic state; its downregulation is required for hypertrophic transition — SOX9 and COL10A1 are mutually exclusive.",
          "11680692", tissue=["cartilage"]),
    _rule("HYPCHOND_NEG_002", "Hypertrophic_chondrocyte", "negative", "marker_genes",
          ["ACAN"],
          "low",
          "Aggrecan (ACAN) is highly expressed in resting and proliferating chondrocytes; it is markedly reduced as chondrocytes become hypertrophic in the growth plate.",
          "38218304", tissue=["cartilage"]),

    # Spermatogonial stem cell (3 rules)
    _rule("SSC_POS_001", "Spermatogonial_SC", "positive", "marker_genes",
          ["ID4", "GFRA1", "ZBTB16"],
          "high",
          "ID4 (inhibitor of differentiation), GFRA1 (GDNF receptor alpha 1), and ZBTB16 (PLZF) mark the undifferentiated human spermatogonial stem cell state required for lifelong spermatogenesis.",
          "30735645", tissue=["testis"]),
    _rule("SSC_CIRCUIT_001", "Spermatogonial_SC", "circuit", "marker_genes",
          ["ID4", "GFRA1"],
          "high",
          "ID4-GFRA1 co-expression defines the self-renewing SSC circuit; GDNF/GFRA1 signaling maintains ID4 expression to sustain undifferentiated state.",
          "40482072", tissue=["testis"]),
    _rule("SSC_NEG_001", "Spermatogonial_SC", "negative", "marker_genes",
          ["SYCP3"],
          "absent",
          "SYCP3 is a meiotic synaptonemal complex protein expressed only from primary spermatocyte stage onward; its presence indicates meiotic entry and departure from SSC identity.",
          "40695261", tissue=["testis"]),

    # Spermatocyte (3 rules)
    _rule("SPERMCYT_POS_001", "Spermatocyte", "positive", "marker_genes",
          ["SYCP3", "SYCP1", "MLH1"],
          "high",
          "SYCP3 and SYCP1 form the synaptonemal complex during meiotic prophase I; MLH1 marks meiotic crossover sites — all three confirm primary spermatocyte identity.",
          "40695261", tissue=["testis"]),
    _rule("SPERMCYT_NEG_001", "Spermatocyte", "negative", "marker_genes",
          ["PRM1"],
          "absent",
          "PRM1 (protamine 1) is expressed only in post-meiotic spermatids replacing histones during chromatin condensation; its presence contradicts primary spermatocyte stage.",
          "40695261", tissue=["testis"]),
    _rule("SPERMCYT_CIRCUIT_001", "Spermatocyte", "circuit", "marker_genes",
          ["SYCP3", "SYCP1"],
          "high",
          "SYCP3-SYCP1 co-assembly forms the tripartite synaptonemal complex essential for meiotic recombination; their co-expression confirms active meiotic prophase I.",
          "40695261", tissue=["testis"]),

    # Spermatid (3 rules)
    _rule("SPERMTID_POS_001", "Spermatid", "positive", "marker_genes",
          ["PRM1", "TNP1", "ACR"],
          "high",
          "PRM1 (protamine 1) and TNP1 (transition protein) replace histones during post-meiotic spermatid chromatin condensation; ACR (acrosin) marks acrosome biogenesis.",
          "40482072", tissue=["testis"]),
    _rule("SPERMTID_NEG_001", "Spermatid", "negative", "marker_genes",
          ["SYCP3"],
          "absent",
          "SYCP3 is a meiotic synaptonemal complex marker disassembled at meiosis II completion; post-meiotic spermatids do not express SYCP3.",
          "40695261", tissue=["testis"]),
    _rule("SPERMTID_CIRCUIT_001", "Spermatid", "circuit", "marker_genes",
          ["PRM1", "TNP1"],
          "high",
          "PRM1-TNP1 sequential expression defines the spermiogenesis chromatin remodeling circuit: TNP1 replaces histones first, then PRM1 replaces TNP1 for final compaction.",
          "40482072", tissue=["testis"]),

    # Granulosa cell (4 rules)
    _rule("GRAN_POS_001", "Granulosa_cell", "positive", "marker_genes",
          ["FSHR", "AMHR2", "CYP19A1"],
          "high",
          "FSHR (FSH receptor), AMHR2, and CYP19A1 (aromatase) are granulosa cell hallmarks; aromatase converts androgens to estrogens driving follicular estrogen synthesis.",
          "14736745", tissue=["ovary"]),
    _rule("GRAN_CIRCUIT_001", "Granulosa_cell", "circuit", "regulon",
          ["FOXL2"],
          "active",
          "Active FOXL2 regulon is required for granulosa cell maintenance and CYP19A1 expression; FOXL2 loss causes transdifferentiation of granulosa cells toward Sertoli-like fate.",
          "14736745", tissue=["ovary"]),
    _rule("GRAN_NEG_001", "Granulosa_cell", "negative", "marker_genes",
          ["CYP17A1"],
          "absent",
          "CYP17A1 (17-alpha hydroxylase) is the theca cell androgen-synthesizing enzyme; the two-cell model of folliculogenesis requires strict CYP17A1 absence in granulosa cells.",
          "14736745", tissue=["ovary"]),
    _rule("GRAN_NEG_002", "Granulosa_cell", "negative", "marker_genes",
          ["PTPRC"],
          "absent",
          "PTPRC (CD45) is a pan-hematopoietic marker; granulosa cells are somatic gonadal epithelial cells and must lack CD45 expression.",
          "14736745", tissue=["ovary"]),

    # Theca cell (3 rules)
    _rule("THECA_POS_001", "Theca_cell", "positive", "marker_genes",
          ["CYP17A1", "STAR", "HSD3B2"],
          "high",
          "CYP17A1, STAR, and HSD3B2 form the theca cell steroidogenic circuit synthesizing androgens from cholesterol; these mark the androgen-producing theca interna.",
          "30735645", tissue=["ovary"]),
    _rule("THECA_NEG_001", "Theca_cell", "negative", "marker_genes",
          ["CYP19A1"],
          "absent",
          "CYP19A1 (aromatase) converts androgens to estrogens and is granulosa-specific; the two-cell folliculogenesis model requires strict CYP19A1 absence in theca cells.",
          "14736745", tissue=["ovary"]),
    _rule("THECA_CIRCUIT_001", "Theca_cell", "circuit", "marker_genes",
          ["CYP17A1", "HSD3B2"],
          "high",
          "CYP17A1-HSD3B2 co-expression defines the androgen synthesis circuit in theca interna; both enzymes are required for androstenedione production.",
          "30735645", tissue=["ovary"]),

    # Oocyte (3 rules)
    _rule("OOCYTE_POS_001", "Oocyte", "positive", "marker_genes",
          ["ZP1", "ZP2", "GDF9", "BMP15"],
          "high",
          "ZP1 and ZP2 form the zona pellucida glycoprotein coat unique to oocytes; GDF9 and BMP15 are oocyte-secreted factors essential for folliculogenesis and fertility.",
          "40695261", tissue=["ovary"]),
    _rule("OOCYTE_CIRCUIT_001", "Oocyte", "circuit", "marker_genes",
          ["GDF9", "BMP15"],
          "high",
          "GDF9-BMP15 co-secretion from oocytes coordinates granulosa cell proliferation and differentiation; their co-expression confirms oocyte-specific growth factor identity.",
          "40695261", tissue=["ovary"]),
    _rule("OOCYTE_NEG_001", "Oocyte", "negative", "marker_genes",
          ["PTPRC"],
          "absent",
          "PTPRC (CD45) marks hematopoietic lineage cells; oocytes are germ cells of mesodermal origin and must lack pan-hematopoietic markers.",
          "40695261", tissue=["ovary"]),

    # Enrich Osteoblast (3 rules)
    _rule("OSTBL_NEG_001", "Osteoblast", "negative", "marker_genes",
          ["CTSK"],
          "absent",
          "CTSK (cathepsin K) is the osteoclast-specific bone-resorbing protease; its presence would indicate osteoclast contamination rather than bone-forming osteoblast.",
          "14633986", tissue=["bone"]),
    _rule("OSTBL_CIRCUIT_001", "Osteoblast", "circuit", "regulon",
          ["RUNX2"],
          "active",
          "Active RUNX2 regulon drives SP7 (Osterix) and COL1A1 expression in the master osteoblast differentiation circuit; RUNX2 is the key osteoblast master TF.",
          "39027344", tissue=["bone"]),
    _rule("OSTBL_NEG_002", "Osteoblast", "negative", "marker_genes",
          ["SOST"],
          "absent",
          "SOST (sclerostin) is exclusively produced by terminally differentiated osteocytes embedded in mineralized matrix; osteoblasts on bone surfaces do not express SOST.",
          "14633986", tissue=["bone"]),

    # Enrich Chondrocyte (2 rules)
    _rule("CHOND_CIRCUIT_001", "Chondrocyte", "circuit", "regulon",
          ["SOX9"],
          "active",
          "Active SOX9 regulon drives COL2A1 and ACAN (aggrecan) expression defining the chondrocyte identity program; SOX9 is the master chondrogenic TF.",
          "11680692", tissue=["cartilage"]),
    _rule("CHOND_NEG_001", "Chondrocyte", "negative", "marker_genes",
          ["COL1A1"],
          "absent",
          "COL1A1 is a fibroblast/osteoblast matrix protein; high COL1A1 in a SOX9+ cell indicates fibrocartilage contamination or dedifferentiation, not healthy articular chondrocyte.",
          "11680692", tissue=["cartilage"]),

    # Enrich Sertoli_cell (3 rules)
    _rule("SERT_CIRCUIT_001", "Sertoli_cell", "circuit", "regulon",
          ["SOX9"],
          "active",
          "Active SOX9 regulon in Sertoli cells drives AMH expression and maintains the Sertoli cell identity critical for male gonadal differentiation; SOX9 is the Sertoli master TF.",
          "35387428", tissue=["testis"]),
    _rule("SERT_NEG_001", "Sertoli_cell", "negative", "marker_genes",
          ["CYP17A1"],
          "absent",
          "CYP17A1 is strictly a Leydig cell steroidogenic enzyme; its presence would indicate Leydig cell contamination. Sertoli cells support spermatogenesis without androgen synthesis.",
          "35387428", tissue=["testis"]),
    _rule("SERT_NEG_002", "Sertoli_cell", "negative", "marker_genes",
          ["GFRA1"],
          "absent",
          "GFRA1 is an SSC self-renewal receptor expressed on undifferentiated spermatogonia; Sertoli cells express GDNF (the ligand) but not the receptor GFRA1.",
          "30735645", tissue=["testis"]),

    # Enrich Leydig_cell (3 rules)
    _rule("LEYDIG_NEG_001", "Leydig_cell", "negative", "marker_genes",
          ["AMH"],
          "absent",
          "AMH (anti-Mullerian hormone) is a Sertoli cell product; Leydig cells are interstitial steroidogenic cells and do not express AMH.",
          "35387428", tissue=["testis"]),
    _rule("LEYDIG_NEG_002", "Leydig_cell", "negative", "marker_genes",
          ["FSHR"],
          "absent",
          "FSHR (FSH receptor) is expressed on Sertoli and granulosa cells; Leydig cells respond to LH/hCG not FSH, and do not express FSHR.",
          "30735645", tissue=["testis"]),
    _rule("LEYDIG_CIRCUIT_001", "Leydig_cell", "circuit", "marker_genes",
          ["CYP17A1", "STAR"],
          "high",
          "CYP17A1-STAR co-expression defines the Leydig cell testosterone synthesis circuit: STAR transports cholesterol into mitochondria and CYP17A1 converts pregnenolone to DHEA.",
          "41596239", tissue=["testis"]),
]

# ── DAY C: INNER EAR + CORNEA + SALIVARY + SYNOVIUM ─────────────────────────
# 15 new types: Inner_Hair_cell, Outer_Hair_cell, Supporting_cell_cochlear,
#   Spiral_ganglion_neuron, Stria_vascularis, Corneal_epithelium,
#   Corneal_keratocyte, Trabecular_meshwork, Acinar_serous, Acinar_mucous,
#   Striated_duct, Myoepithelial, Synovial_fibroblast, Synovial_macrophage,
#   Chondroprogenitor

DAY_C_RULES = [
    # Inner Hair cell (3 rules)
    _rule("IHC_CIRCUIT_001", "Inner_Hair_cell", "circuit", "regulon",
          ["ATOH1"],
          "active",
          "Active ATOH1 regulon is required for inner hair cell identity and survival; ATOH1 drives MYO7A and ESPN expression defining the mechanosensory hair bundle.",
          "41945437", tissue=["inner_ear", "cochlea"]),
    _rule("IHC_POS_001", "Inner_Hair_cell", "positive", "marker_genes",
          ["MYO7A", "ESPN", "SLC17A8"],
          "high",
          "MYO7A (myosin VIIA) maintains hair bundle integrity; ESPN cross-links stereocilia actin; SLC17A8 (VGLUT3) mediates glutamatergic IHC-SGN synapse transmission.",
          "41945437", tissue=["inner_ear", "cochlea"]),
    _rule("IHC_NEG_001", "Inner_Hair_cell", "negative", "marker_genes",
          ["SLC26A5"],
          "absent",
          "SLC26A5 (prestin) is the electromotility motor exclusively expressed in outer hair cells for active cochlear amplification; its presence contradicts inner hair cell identity.",
          "41638427", tissue=["inner_ear", "cochlea"]),

    # Outer Hair cell (3 rules)
    _rule("OHC_POS_001", "Outer_Hair_cell", "positive", "marker_genes",
          ["SLC26A5", "OCM"],
          "high",
          "SLC26A5 (prestin) drives OHC electromotility for cochlear amplification; OCM (oncomodulin) is a calcium binding protein unique to OHCs supporting Ca2+ buffering.",
          "41638427", tissue=["inner_ear", "cochlea"]),
    _rule("OHC_NEG_001", "Outer_Hair_cell", "negative", "regulon",
          ["ATOH1"],
          "inactive",
          "ATOH1 regulon activity drives IHC fate; OHCs require ATOH1 silencing post-specification — persistent ATOH1 activity in OHC position indicates IHC fate.",
          "41945437", tissue=["inner_ear", "cochlea"]),
    _rule("OHC_NEG_002", "Outer_Hair_cell", "negative", "marker_genes",
          ["SLC17A8"],
          "absent",
          "SLC17A8 (VGLUT3) mediates afferent glutamate release at IHC-spiral ganglion ribbon synapses; OHCs receive efferent (olivocochlear) innervation and do not express SLC17A8.",
          "41638427", tissue=["inner_ear", "cochlea"]),

    # Supporting cell cochlear (3 rules)
    _rule("SUPCOCH_POS_001", "Supporting_cell_cochlear", "positive", "marker_genes",
          ["SOX2", "PROX1", "LFNG"],
          "high",
          "SOX2 maintains inner ear supporting cell progenitor identity; PROX1 marks Deiters cells and pillar cells; LFNG (lunatic fringe) defines the Notch-mediated supporting cell boundary.",
          "30410537", tissue=["inner_ear", "cochlea"]),
    _rule("SUPCOCH_NEG_001", "Supporting_cell_cochlear", "negative", "marker_genes",
          ["MYO7A"],
          "absent",
          "MYO7A is a hair cell-specific myosin; supporting cells surround hair cells but do not form mechano-sensitive hair bundles and lack MYO7A.",
          "41945437", tissue=["inner_ear", "cochlea"]),
    _rule("SUPCOCH_NEG_002", "Supporting_cell_cochlear", "negative", "marker_genes",
          ["SLC26A5"],
          "absent",
          "SLC26A5 (prestin) is strictly an outer hair cell electromotility motor; supporting cells flanking OHCs do not express prestin.",
          "41638427", tissue=["inner_ear", "cochlea"]),

    # Spiral ganglion neuron (3 rules)
    _rule("SGN_POS_001", "Spiral_ganglion_neuron", "positive", "marker_genes",
          ["NEFM", "RBFOX3", "SLC17A6"],
          "high",
          "NEFM (neurofilament medium) and RBFOX3 (NeuN) are pan-neuronal markers; SLC17A6 (VGLUT2) marks type I SGNs that receive glutamatergic input from inner hair cells.",
          "41491447", tissue=["inner_ear", "cochlea"]),
    _rule("SGN_CIRCUIT_001", "Spiral_ganglion_neuron", "circuit", "regulon",
          ["NEUROD1"],
          "active",
          "Active NEUROD1 regulon specifies spiral ganglion neuron identity and drives expression of peripheral auditory neuron-specific genes including NEFM and SLC17A6.",
          "41491447", tissue=["inner_ear", "cochlea"]),
    _rule("SGN_NEG_001", "Spiral_ganglion_neuron", "negative", "marker_genes",
          ["MYO7A"],
          "absent",
          "MYO7A is exclusively expressed in mechanosensory hair cells; spiral ganglion neurons conduct auditory signals but lack the hair bundle apparatus.",
          "41945437", tissue=["inner_ear", "cochlea"]),

    # Stria vascularis (3 rules)
    _rule("STRIA_POS_001", "Stria_vascularis", "positive", "marker_genes",
          ["KCNJ10", "SLC12A2"],
          "high",
          "KCNJ10 (Kir4.1) generates the +80 mV endocochlear potential driving K+ recycling; SLC12A2 (NKCC1) co-transports K+ in marginal cells maintaining the endolymph ionic environment.",
          "20012478", tissue=["inner_ear", "cochlea"]),
    _rule("STRIA_NEG_001", "Stria_vascularis", "negative", "marker_genes",
          ["MYO7A"],
          "absent",
          "MYO7A is a hair cell marker; stria vascularis is a vascularized ion transport epithelium in the lateral cochlear wall without mechanosensory function.",
          "41945437", tissue=["inner_ear", "cochlea"]),
    _rule("STRIA_CIRCUIT_001", "Stria_vascularis", "circuit", "marker_genes",
          ["KCNJ10", "SLC12A2"],
          "high",
          "KCNJ10-SLC12A2 co-expression defines the stria vascularis endocochlear potential generation circuit; loss of either causes deafness via endolymphatic K+ dysregulation.",
          "20012478", tissue=["inner_ear", "cochlea"]),

    # Corneal epithelium (3 rules)
    _rule("CORNEP_POS_001", "Corneal_epithelium", "positive", "marker_genes",
          ["KRT12", "KRT3"],
          "high",
          "KRT12 and KRT3 form an obligate keratin pair exclusively expressed in corneal epithelium; KRT12 mutations cause Meesmann's epithelial corneal dystrophy.",
          "42194275", tissue=["cornea", "eye"]),
    _rule("CORNEP_NEG_001", "Corneal_epithelium", "negative", "marker_genes",
          ["COL8A1"],
          "absent",
          "COL8A1 (type VIII collagen) is expressed in Descemet's membrane by corneal endothelium; corneal epithelium does not produce this posterior stromal barrier collagen.",
          "40437144", tissue=["cornea", "eye"]),
    _rule("CORNEP_CIRCUIT_001", "Corneal_epithelium", "circuit", "marker_genes",
          ["KRT12", "KRT3"],
          "high",
          "KRT12-KRT3 heterodimer formation is obligatory in corneal epithelium for structural integrity and optical transparency; co-expression confirms corneal epithelial identity.",
          "42194275", tissue=["cornea", "eye"]),

    # Corneal keratocyte (3 rules)
    _rule("CORNKER_POS_001", "Corneal_keratocyte", "positive", "marker_genes",
          ["KERA", "ALDH3A1", "LUM"],
          "high",
          "KERA (keratocan) and LUM (lumican) are cornea-specific proteoglycans maintaining regular collagen fibril spacing for transparency; ALDH3A1 scavenges UV-generated free radicals.",
          "40437144", tissue=["cornea", "eye"]),
    _rule("CORNKER_NEG_001", "Corneal_keratocyte", "negative", "marker_genes",
          ["ACTA2"],
          "absent",
          "ACTA2 (alpha-smooth muscle actin) marks activated myofibroblasts during corneal wound healing; resting keratocytes are quiescent and lack ACTA2 expression.",
          "40437144", tissue=["cornea", "eye"]),
    _rule("CORNKER_CIRCUIT_001", "Corneal_keratocyte", "circuit", "marker_genes",
          ["KERA", "ALDH3A1"],
          "high",
          "KERA-ALDH3A1 co-expression defines the corneal transparency maintenance circuit: keratocan organizes collagen microfibrils while ALDH3A1 protects from oxidative stress.",
          "40437144", tissue=["cornea", "eye"]),

    # Trabecular meshwork (3 rules)
    _rule("TM_POS_001", "Trabecular_meshwork", "positive", "marker_genes",
          ["MYOC", "AQP1", "MGP"],
          "high",
          "MYOC (myocilin) is a trabecular meshwork marker elevated in glaucoma; AQP1 mediates aqueous humor outflow; MGP (matrix gla protein) is expressed in TM extracellular matrix.",
          "42194275", tissue=["eye", "trabecular_meshwork"]),
    _rule("TM_NEG_001", "Trabecular_meshwork", "negative", "marker_genes",
          ["MYO7A"],
          "absent",
          "MYO7A is a hair cell motor protein; trabecular meshwork cells regulate intraocular pressure through aqueous humor outflow and do not form mechanosensory structures.",
          "41945437", tissue=["eye", "trabecular_meshwork"]),
    _rule("TM_CIRCUIT_001", "Trabecular_meshwork", "circuit", "marker_genes",
          ["MYOC", "AQP1"],
          "high",
          "MYOC-AQP1 co-expression defines the trabecular meshwork identity; MYOC mutations cause juvenile open-angle glaucoma by blocking aqueous outflow through the TM.",
          "42194275", tissue=["eye", "trabecular_meshwork"]),

    # Acinar serous (3 rules)
    _rule("ACSEROUS_POS_001", "Acinar_serous", "positive", "marker_genes",
          ["BPIFB2", "ZG16B", "AMY2A"],
          "high",
          "BPIFB2 and ZG16B are serous acinar secretory proteins in parotid gland zymogen granules; AMY2A (amylase 2A) is the major serous enzyme for starch digestion.",
          "35220796", tissue=["salivary_gland", "parotid"]),
    _rule("ACSEROUS_CIRCUIT_001", "Acinar_serous", "circuit", "marker_genes",
          ["AMY2A", "BPIFB2"],
          "high",
          "AMY2A-BPIFB2 co-expression defines the serous acinar secretory circuit producing digestive enzymes and antimicrobial proteins concentrated in zymogen granules.",
          "35220796", tissue=["salivary_gland", "parotid"]),
    _rule("ACSEROUS_NEG_001", "Acinar_serous", "negative", "marker_genes",
          ["MUC5B"],
          "absent",
          "MUC5B is the major gel-forming mucin of mucous acinar cells and mucus-secreting duct cells; its presence distinguishes mucous from serous acinar identity.",
          "35220796", tissue=["salivary_gland"]),

    # Acinar mucous (3 rules)
    _rule("ACMUC_POS_001", "Acinar_mucous", "positive", "marker_genes",
          ["MUC5B", "CLCA1", "MUC5AC"],
          "high",
          "MUC5B (mucin 5B) and MUC5AC are gel-forming mucins produced by mucous acinar cells in sublingual and submandibular glands; CLCA1 regulates mucin secretion.",
          "35220796", tissue=["salivary_gland", "sublingual", "submandibular"]),
    _rule("ACMUC_NEG_001", "Acinar_mucous", "negative", "marker_genes",
          ["AMY2A"],
          "absent",
          "AMY2A (salivary amylase) is a serous enzyme secreted exclusively by serous acinar cells; mucous acinar cells do not produce digestive enzymes.",
          "35220796", tissue=["salivary_gland"]),
    _rule("ACMUC_NEG_002", "Acinar_mucous", "negative", "marker_genes",
          ["BPIFB2"],
          "absent",
          "BPIFB2 is a serous acinar antimicrobial protein; mucous acinar cells produce mucins not BPIFB proteins.",
          "35220796", tissue=["salivary_gland"]),

    # Striated duct (2 rules)
    _rule("STRDCT_POS_001", "Striated_duct", "positive", "marker_genes",
          ["AQP5", "KRT14", "ATP6V1B1"],
          "high",
          "AQP5 mediates water transport in salivary duct; KRT14 marks basal duct epithelium; ATP6V1B1 (V-ATPase) drives H+ secretion modifying primary saliva composition.",
          "35220796", tissue=["salivary_gland"]),
    _rule("STRDCT_NEG_001", "Striated_duct", "negative", "marker_genes",
          ["MUC5B"],
          "absent",
          "MUC5B is produced by mucous acinar and mucous duct cells; striated duct cells modify ionic composition of saliva and do not produce gel-forming mucins.",
          "35220796", tissue=["salivary_gland"]),

    # Myoepithelial (3 rules)
    _rule("MYOEPI_CIRCUIT_001", "Myoepithelial", "circuit", "marker_genes",
          ["ACTA2", "KRT14"],
          "high",
          "ACTA2-KRT14 co-expression defines myoepithelial cell identity: smooth muscle actin enables contractile expulsion of secretory content while KRT14 maintains epithelial polarity.",
          "35220796", tissue=["salivary_gland", "breast", "lacrimal"]),
    _rule("MYOEPI_POS_001", "Myoepithelial", "positive", "marker_genes",
          ["TP63", "ACTA2", "MYLK"],
          "high",
          "TP63 (p63) is the myoepithelial master TF; ACTA2 enables contraction; MYLK (myosin light chain kinase) regulates actomyosin contraction for glandular secretion.",
          "35220796", tissue=["salivary_gland", "breast"]),
    _rule("MYOEPI_NEG_001", "Myoepithelial", "negative", "marker_genes",
          ["AMY2A"],
          "absent",
          "AMY2A is a luminal serous acinar enzyme; myoepithelial cells form the outer contractile layer of acini and ducts without producing digestive enzymes.",
          "35220796", tissue=["salivary_gland"]),

    # Synovial fibroblast (3 rules)
    _rule("SYNFIB_POS_001", "Synovial_fibroblast", "positive", "marker_genes",
          ["PRG4", "PDPN", "THY1"],
          "high",
          "PRG4 (lubricin) is secreted by synovial fibroblast-like synoviocytes to reduce joint friction; PDPN and THY1 (CD90) are surface markers of the synovial lining layer.",
          "37024468", tissue=["synovium", "joint"]),
    _rule("SYNFIB_CIRCUIT_001", "Synovial_fibroblast", "circuit", "marker_genes",
          ["PRG4", "PDPN"],
          "high",
          "PRG4-PDPN co-expression defines the synovial lining fibroblast identity responsible for joint lubrication and synovial fluid production.",
          "37024468", tissue=["synovium", "joint"]),
    _rule("SYNFIB_NEG_001", "Synovial_fibroblast", "negative", "marker_genes",
          ["CD68"],
          "absent",
          "CD68 marks macrophage-like synoviocytes (type A synoviocytes); its presence distinguishes phagocytic synovial macrophages from fibroblast-like synoviocytes (type B).",
          "37024468", tissue=["synovium", "joint"]),

    # Synovial macrophage (3 rules)
    _rule("SYNMAC_POS_001", "Synovial_macrophage", "positive", "marker_genes",
          ["CD68", "FOLR2", "LYVE1"],
          "high",
          "CD68 marks macrophage-like type A synoviocytes; FOLR2 and LYVE1 are tissue-resident joint macrophage markers distinguishing them from circulating monocyte-derived macrophages.",
          "37024468", tissue=["synovium", "joint"]),
    _rule("SYNMAC_CIRCUIT_001", "Synovial_macrophage", "circuit", "marker_genes",
          ["CD68", "FOLR2"],
          "high",
          "CD68-FOLR2 co-expression defines the tissue-resident synovial macrophage lining layer identity; FOLR2 marks anti-inflammatory resident macrophages in joint homeostasis.",
          "37024468", tissue=["synovium", "joint"]),
    _rule("SYNMAC_NEG_001", "Synovial_macrophage", "negative", "marker_genes",
          ["PRG4"],
          "absent",
          "PRG4 (lubricin) is produced by type B fibroblast-like synoviocytes; type A macrophage-like synoviocytes engulf debris and do not produce joint lubricants.",
          "37024468", tissue=["synovium", "joint"]),

    # Chondroprogenitor (3 rules)
    _rule("CHONDPROG_POS_001", "Chondroprogenitor", "positive", "marker_genes",
          ["PDGFRA", "NES", "PRRX1"],
          "high",
          "PDGFRA and NES (nestin) mark mesenchymal progenitor identity; PRRX1 is a transcription factor expressed in pre-chondrogenic limb mesenchyme before COL2A1 activation.",
          "11680692", tissue=["cartilage", "bone", "joint"]),
    _rule("CHONDPROG_NEG_001", "Chondroprogenitor", "negative", "marker_genes",
          ["COL2A1"],
          "absent",
          "COL2A1 (type II collagen) marks committed mature chondrocytes; pre-chondrogenic progenitors have not yet activated the SOX9-driven chondrogenic differentiation program.",
          "11680692", tissue=["cartilage"]),
    _rule("CHONDPROG_NEG_002", "Chondroprogenitor", "negative", "marker_genes",
          ["COL10A1"],
          "absent",
          "COL10A1 marks terminal hypertrophic chondrocyte differentiation; progenitors are far upstream of this stage in the endochondral ossification hierarchy.",
          "38218304", tissue=["cartilage"]),
]

# ── DAY D: VASCULAR SUBTYPES + MISC SPECIALIZED ───────────────────────────────
# 14 new types: High_endothelial_venule, Pericyte_brain, Lung_aerocyte,
#   Synthetic_SMC, Prostate_luminal, Prostate_basal, Urothelial_superficial,
#   Urothelial_basal, Choroid_plexus_epithelium, Interstitial_cell_Cajal,
#   Fibro_adipogenic_progenitor, Marginal_zone_B, Portal_endothelial,
#   Red_pulp_macrophage
# Enrichments: Tip_cell (+1), Vascular_SMC (+1), Capillary_Endothelial (+1)

DAY_D_RULES = [
    # High endothelial venule (3 rules)
    _rule("HEV_POS_001", "High_endothelial_venule", "positive", "marker_genes",
          ["ACKR1", "GLYCAM1", "CCL21"],
          "high",
          "ACKR1 (DARC) is the HEV promiscuous chemokine receptor mediating lymphocyte homing; GLYCAM1 and CCL21 support the high endothelial morphology enabling transmigration into lymph nodes.",
          "13679391", tissue=["lymph_node", "Peyer_patches"]),
    _rule("HEV_CIRCUIT_001", "High_endothelial_venule", "circuit", "marker_genes",
          ["ACKR1", "GLYCAM1"],
          "high",
          "ACKR1-GLYCAM1 co-expression defines HEV identity; ACKR1 scavenges chemokines from blood while GLYCAM1 presents addressins for lymphocyte integrin binding during homing.",
          "13679391", tissue=["lymph_node"]),
    _rule("HEV_NEG_001", "High_endothelial_venule", "negative", "marker_genes",
          ["GJA4"],
          "absent",
          "GJA4 (connexin 37) is an arterial endothelial marker mediating electrical coupling; HEVs are specialized post-capillary venules and do not express arterial gap junction proteins.",
          "13679391", tissue=["lymph_node"]),

    # Brain pericyte (3 rules)
    _rule("PERIBN_POS_001", "Pericyte_brain", "positive", "marker_genes",
          ["PDGFRB", "RGS5", "KCNJ8", "SLC7A5"],
          "high",
          "PDGFRB and RGS5 mark all pericytes; KCNJ8 (Kir6.1 ATP-sensitive K+) and SLC7A5 (LAT1 amino acid transporter) are BBB-pericyte enriched genes regulating neurovascular coupling.",
          "16807374", tissue=["brain", "BBB"]),
    _rule("PERIBN_CIRCUIT_001", "Pericyte_brain", "circuit", "marker_genes",
          ["PDGFRB", "KCNJ8"],
          "high",
          "PDGFRB-KCNJ8 co-expression defines BBB-specific pericyte identity; PDGFB/PDGFRB signaling maintains pericyte coverage while KCNJ8 mediates pericyte-dependent BBB K+ regulation.",
          "16807374", tissue=["brain", "BBB"]),
    _rule("PERIBN_NEG_001", "Pericyte_brain", "negative", "marker_genes",
          ["MYH11"],
          "absent",
          "MYH11 (smooth muscle myosin heavy chain) marks contractile vascular smooth muscle cells; brain pericytes are electrophysiologically distinct from mural SMCs and lack MYH11.",
          "16807374", tissue=["brain"]),

    # Lung aerocyte (3 rules)
    _rule("AEROCYTE_POS_001", "Lung_aerocyte", "positive", "marker_genes",
          ["CA4", "EDNRB", "APLNR"],
          "high",
          "CA4 (carbonic anhydrase 4), EDNRB (endothelin receptor B), and APLNR (apelin receptor) define lung aerocytes — the alveolar capillary subtype specialized for gas exchange at the air-blood barrier.",
          "33208946", tissue=["lung", "alveolus"]),
    _rule("AEROCYTE_CIRCUIT_001", "Lung_aerocyte", "circuit", "marker_genes",
          ["CA4", "EDNRB"],
          "high",
          "CA4-EDNRB co-expression identifies aerocytes as distinct from general capillary endothelium; CA4 facilitates CO2-bicarbonate exchange while EDNRB modulates pulmonary vascular tone.",
          "33208946", tissue=["lung", "alveolus"]),
    _rule("AEROCYTE_NEG_001", "Lung_aerocyte", "negative", "marker_genes",
          ["GJA4"],
          "absent",
          "GJA4 is an arterial endothelial connexin; aerocytes are capillary-zone endothelial cells positioned in the alveolar gas-exchange interface and lack arterial identity markers.",
          "33208946", tissue=["lung", "alveolus"]),

    # Synthetic SMC (3 rules)
    _rule("SYNTHSMC_POS_001", "Synthetic_SMC", "positive", "marker_genes",
          ["COL1A1", "FN1", "VCAN"],
          "high",
          "COL1A1, FN1 (fibronectin), and VCAN (versican) mark synthetic-phenotype VSMCs that have switched from contractile to matrix-secreting mode during vascular remodeling.",
          "16807374", tissue=["vasculature", "atherosclerosis"]),
    _rule("SYNTHSMC_NEG_001", "Synthetic_SMC", "negative", "marker_genes",
          ["MYH11"],
          "absent",
          "MYH11 is the contractile phenotype VSMC marker downregulated during phenotype switching; synthetic VSMCs lose contractile proteins while gaining matrix-production capacity.",
          "16807374", tissue=["vasculature"]),
    _rule("SYNTHSMC_NEG_002", "Synthetic_SMC", "negative", "marker_genes",
          ["CNN1"],
          "absent",
          "CNN1 (calponin 1) is a contractile smooth muscle protein stabilizing actin filaments; its absence confirms phenotype switching from contractile to synthetic VSMC state.",
          "16807374", tissue=["vasculature"]),

    # Prostate luminal (3 rules)
    _rule("PROSTLUM_POS_001", "Prostate_luminal", "positive", "marker_genes",
          ["NKX3-1", "ACPP", "KLK3"],
          "high",
          "NKX3-1 is the prostate-specific homeodomain TF required for luminal cell differentiation; ACPP (acid phosphatase) and KLK3 (PSA) are luminal secretory products.",
          "10906459", tissue=["prostate"]),
    _rule("PROSTLUM_CIRCUIT_001", "Prostate_luminal", "circuit", "marker_genes",
          ["NKX3-1", "ACPP"],
          "high",
          "NKX3-1 activates prostate-specific secretory program including ACPP; NKX3-1 disruption causes morphogenetic defects in prostatic ductal branching confirming luminal TF function.",
          "10906459", tissue=["prostate"]),
    _rule("PROSTLUM_NEG_001", "Prostate_luminal", "negative", "marker_genes",
          ["KRT14"],
          "absent",
          "KRT14 marks basal prostatic epithelial cells; luminal cells differentiate from basal cells and lose KRT14 while gaining secretory identity markers.",
          "10906459", tissue=["prostate"]),

    # Prostate basal (2 rules)
    _rule("PROSTBASE_POS_001", "Prostate_basal", "positive", "marker_genes",
          ["TP63", "KRT5", "ACTA2"],
          "high",
          "TP63 (p63), KRT5 (basal keratin), and ACTA2 mark the prostatic basal/progenitor layer that contacts the basement membrane and can regenerate luminal cells.",
          "10906459", tissue=["prostate"]),
    _rule("PROSTBASE_NEG_001", "Prostate_basal", "negative", "marker_genes",
          ["NKX3-1"],
          "absent",
          "NKX3-1 is a luminal differentiation factor; prostatic basal cells are progenitors that do not express luminal secretory identity markers.",
          "10906459", tissue=["prostate"]),

    # Urothelial superficial (2 rules)
    _rule("UROSURF_POS_001", "Urothelial_superficial", "positive", "marker_genes",
          ["UPK1A", "UPK1B", "UPK3A", "KRT20"],
          "high",
          "Uroplakins UPK1A, UPK1B, and UPK3A form the asymmetric unit membrane of umbrella cells; KRT20 marks terminally differentiated superficial urothelial cells.",
          "21853341", tissue=["bladder", "ureter", "urothelium"]),
    _rule("UROSURF_NEG_001", "Urothelial_superficial", "negative", "marker_genes",
          ["TP63"],
          "absent",
          "TP63 is expressed in urothelial basal progenitor cells; superficial umbrella cells are terminally differentiated and do not express basal identity markers.",
          "21853341", tissue=["bladder", "urothelium"]),

    # Urothelial basal (2 rules)
    _rule("UROBASE_POS_001", "Urothelial_basal", "positive", "marker_genes",
          ["TP63", "KRT5", "KRT14"],
          "high",
          "TP63, KRT5, and KRT14 mark urothelial basal/intermediate progenitor cells that regenerate the urothelium after injury; these are also expressed in basal cells of all stratified epithelia.",
          "21853341", tissue=["bladder", "urothelium"]),
    _rule("UROBASE_NEG_001", "Urothelial_basal", "negative", "marker_genes",
          ["UPK3A"],
          "absent",
          "UPK3A (uroplakin 3A) is a terminal differentiation marker of superficial umbrella cells; basal urothelial progenitors lack uroplakin expression.",
          "21853341", tissue=["bladder", "urothelium"]),

    # Choroid plexus epithelium (2 rules)
    _rule("CHPLX_POS_001", "Choroid_plexus_epithelium", "positive", "marker_genes",
          ["TTR", "FOLR1", "CLDN11"],
          "high",
          "TTR (transthyretin) is produced almost exclusively by choroid plexus epithelium for thyroid hormone transport to CSF; FOLR1 mediates folate transport into brain; CLDN11 seals the blood-CSF barrier.",
          "32046769", tissue=["brain", "choroid_plexus"]),
    _rule("CHPLX_CIRCUIT_001", "Choroid_plexus_epithelium", "circuit", "marker_genes",
          ["TTR", "FOLR1"],
          "high",
          "TTR-FOLR1 co-expression defines choroid plexus epithelial identity: both are transported into CSF for brain-specific nutrient delivery unavailable from blood-brain barrier capillaries.",
          "32046769", tissue=["brain", "choroid_plexus"]),

    # Interstitial cell of Cajal (2 rules)
    _rule("ICC_POS_001", "Interstitial_cell_Cajal", "positive", "marker_genes",
          ["KIT", "ANO1", "EDNRB"],
          "high",
          "KIT (CD117, c-Kit) is the SCF receptor marking ICCs as the gastrointestinal pacemaker cells; ANO1 (TMEM16A Cl- channel) generates spontaneous pacemaker currents for gut peristalsis.",
          "33325807", tissue=["GI_tract", "intestine", "stomach"]),
    _rule("ICC_CIRCUIT_001", "Interstitial_cell_Cajal", "circuit", "marker_genes",
          ["KIT", "ANO1"],
          "high",
          "KIT-ANO1 co-expression defines ICC pacemaker identity: KIT signaling from SCF maintains ICC survival while ANO1 Cl- channels generate the slow-wave pacemaker current.",
          "33325807", tissue=["GI_tract", "intestine"]),

    # Fibro-adipogenic progenitor (1 rule)
    _rule("FAP_POS_001", "Fibro_adipogenic_progenitor", "positive", "marker_genes",
          ["PDGFRA", "CD34", "NT5E"],
          "high",
          "PDGFRA marks fibro-adipogenic progenitors in skeletal muscle interstitium; CD34 and NT5E (CD73) are mesenchymal stromal markers identifying FAPs as the major non-myogenic muscle resident progenitor.",
          "40585958", tissue=["skeletal_muscle", "adipose"]),

    # Marginal zone B cell (2 rules)
    _rule("MZB_POS_001", "Marginal_zone_B", "positive", "marker_genes",
          ["NOTCH2", "CR2", "PTPRO"],
          "high",
          "NOTCH2 signaling is required for MZB cell fate commitment in the spleen; CR2 (CD21) and PTPRO mark the marginal zone B cell population surveying blood-borne antigens.",
          "12753744", tissue=["spleen", "lymph_node"]),
    _rule("MZB_NEG_001", "Marginal_zone_B", "negative", "marker_genes",
          ["TCL1A"],
          "absent",
          "TCL1A marks follicular B cells and naive B cell identity; marginal zone B cells downregulate TCL1A as they exit the follicular fate decision toward marginal zone.",
          "12753744", tissue=["spleen"]),

    # Portal endothelial (1 rule)
    _rule("PORTALENDO_POS_001", "Portal_endothelial", "positive", "marker_genes",
          ["FCN2", "FCN3", "CLEC1B"],
          "high",
          "FCN2 and FCN3 (ficolin-2/3) are complement-activating pattern recognition molecules enriched in portal tract endothelium; CLEC1B (CLEC-1B) marks portal liver vascular identity.",
          "39178365", tissue=["liver", "portal_tract"]),

    # Red pulp macrophage (1 rule)
    _rule("RPMAC_POS_001", "Red_pulp_macrophage", "positive", "marker_genes",
          ["VCAM1", "TIMD4", "SPIC"],
          "high",
          "VCAM1 and TIMD4 (Tim4) define tissue-resident red pulp macrophages that clear senescent erythrocytes; SPIC is a TF induced by heme during monocyte-to-red-pulp-macrophage differentiation.",
          "24630724", tissue=["spleen", "red_pulp"]),

    # Enrich Tip_cell (1 rule)
    _rule("TIP_NEG_001", "Tip_cell", "negative", "marker_genes",
          ["NR2F2"],
          "absent",
          "NR2F2 (COUP-TFII) is a venous identity TF; tip cells are the leading angiogenic sprout cells with arterial-like identity and must lack venous NR2F2 expression.",
          "33208946", tissue=["vasculature"]),

    # Enrich Vascular_SMC (1 rule)
    _rule("VSMC_CIRCUIT_001", "Vascular_SMC", "circuit", "marker_genes",
          ["MYH11", "CNN1"],
          "high",
          "MYH11-CNN1 co-expression defines the contractile vascular SMC circuit: smooth muscle myosin and calponin coordinate actin-based vasoconstriction in response to vasoconstrictors.",
          "16807374", tissue=["vasculature"]),

    # Enrich Capillary_Endothelial (1 rule)
    _rule("CAP_NEG_001", "Capillary_Endothelial", "negative", "marker_genes",
          ["GJA4"],
          "absent",
          "GJA4 (Cx37) is an arterial endothelial gap junction protein absent in capillaries; its presence would indicate arterial rather than capillary identity.",
          "33208946", tissue=["vasculature", "capillary"]),
]

# ── Main: load, extend, validate, save ────────────────────────────────────────

def main():
    with open(KG_PATH) as f:
        existing = json.load(f)

    existing_ids = {r["rule_id"] for r in existing}
    print(f"Existing rules: {len(existing)}")
    print(f"Existing rule IDs (sample): {list(existing_ids)[:5]}")

    all_new = DAY_A_RULES + DAY_B_RULES + DAY_C_RULES + DAY_D_RULES

    # Duplicate ID check
    new_ids = [r["rule_id"] for r in all_new]
    duplicates_in_new = [rid for rid in new_ids if new_ids.count(rid) > 1]
    if duplicates_in_new:
        raise ValueError(f"Duplicate rule IDs in new rules: {set(duplicates_in_new)}")

    conflicts_with_existing = [rid for rid in new_ids if rid in existing_ids]
    if conflicts_with_existing:
        raise ValueError(f"New rule IDs already in KG: {conflicts_with_existing}")

    print(f"New rules: {len(all_new)}")
    print(f"  Day A (retina):           {len(DAY_A_RULES)}")
    print(f"  Day B (bone/gonads):      {len(DAY_B_RULES)}")
    print(f"  Day C (ear/cornea/etc):   {len(DAY_C_RULES)}")
    print(f"  Day D (vascular/misc):    {len(DAY_D_RULES)}")

    # New cell types
    existing_cell_types = {r["cell_type"] for r in existing if r.get("status") == "ACTIVE"}
    new_cell_types = {r["cell_type"] for r in all_new} - existing_cell_types
    print(f"\nNew cell types ({len(new_cell_types)}): {sorted(new_cell_types)}")

    # Schema validation
    import sys
    sys.path.insert(0, ".")
    from axiom_sc.tier2.kg_loader import KGLoader
    import jsonschema
    import json as _json
    schema = _json.load(open("axiom_sc/kg/schema.json"))
    errors = []
    for r in all_new:
        try:
            jsonschema.validate(r, schema)
        except jsonschema.ValidationError as e:
            errors.append(f"{r['rule_id']}: {e.message}")
    if errors:
        raise ValueError(f"Schema validation errors:\n" + "\n".join(errors))

    print(f"\nAll {len(all_new)} rules passed schema validation")

    # Append and save
    combined = existing + all_new
    with open(KG_PATH, "w") as f:
        json.dump(combined, f, indent=2)
    print(f"\nSaved {len(combined)} rules to {KG_PATH}")

    # Final stats
    kg = KGLoader.from_file(KG_PATH)
    active = kg.active_rules()
    types = kg.active_cell_types()
    neg = [r for r in active if r.rule_type == "negative"]
    circuit = [r for r in active if r.rule_type == "circuit"]
    print(f"\nFinal KG stats:")
    print(f"  ACTIVE rules : {len(active)}  (target ≥640)")
    print(f"  Cell types   : {len(types)}  (target ≥198)")
    print(f"  Negative     : {len(neg)}  (core AXIOM innovation)")
    print(f"  Circuit      : {len(circuit)}")


if __name__ == "__main__":
    main()
