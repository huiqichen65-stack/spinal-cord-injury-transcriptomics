# Analysis notes

Interpretation notes originally emitted in per-run prose reports. Numerical results remain in the analysis tables.

## bulk_wgcna/44_bulk_input_qc.R

- Bulk RNA-seq Step 44 | Input audit + QC
- Expected main bulk time points:
- NEXT DECISION:
- Output files:

## bulk_wgcna/45_bulk_analysis_object.R

- Bulk RNA-seq Step 45
- Batch structure assessment + formal analysis object
- Operational batch definition:
- Batch_A = jerry.* + normal.*
- Batch_B = D*
- phase definition:
- Early = Uninjured + 15 min + 1 dpi + 3 dpi
- Transition = 7 dpi
- Established late = 14 dpi + 28 dpi
- Extended late = 42 dpi
- Primary formal comparison:
- Established late vs Early
- Formal model 1:
- ~ batch + timepoint
- Formal model 2:
- ~ batch + phase_primary
- IMPORTANT:
- No samples were removed automatically.
- Batch-adjusted expression is visualization/QC only.
- The original transformed matrix must be used in formal limma models.

## bulk_wgcna/46_bulk_temporal_differential_expression.R

- Bulk RNA-seq Step 46
- Batch-aware temporal differential-expression analysis
- PRIMARY COMPARISON:
- Established late (14 + 28 dpi) vs Early (Uninjured + 15 min + 1 + 3 dpi)
- Primary phase-by-batch table:
- Primary DEG counts:
- Transition 7 dpi vs Early:
- Extended late 42 dpi vs Early:
- 28 dpi vs 14 dpi:
- Primary sensitivity analysis:
- IMPORTANT INTERPRETATION:
- 15 min vs Uninjured and 1 dpi vs Uninjured are exploratory because those two
- early time points are perfectly separated by the operational batch assignment.
- The primary phase contrast is substantially better balanced across batches.
- No manuscript figure was generated in Step46.

## bulk_wgcna/47_Clec7a_signature_GSEA.R

- Bulk RNA-seq Step 47
- scRNA-derived Clec7a+ lipid-handling signature GSEA
- PRIMARY GSEA:
- Established late (14+28 dpi) vs Early
- INTERPRETATION:
- TEMPORAL PROFILE:
- ROBUSTNESS:
- IMPORTANT:
- 1 dpi vs Uninjured is exploratory because of the previously identified operational-batch confounding.
- No manuscript figure was generated in Step47.

## bulk_wgcna/48_bulk_functional_enrichment.R

- Step48 v04 — Bulk functional enrichment report
- contrast:
- Established late (14+28 dpi) vs Early (Uninjured+15min+1dpi+3dpi)
- Step46 checks:
- Step46 gene identity handling:
- Stable unique ID = gene_id
- Annotation / core-gene matching = gene_symbol_restored
- Annotation mapping table = gene_symbol_restored (SYMBOL) -> ENTREZID
- ORA mapping QC:
- ORA complete-result counts (BH FDR < 0.05):
- GSEA setup:
- Ranking metric = limma moderated t statistic
- Duplicate Entrez rule = retain the tested row with largest absolute moderated t
- GSEA significant pathway counts (BH FDR < 0.05):
- Pre-specified biological-theme summary (BH FDR < 0.05):
- Interpretation rule:
- "enhanced lipid handling/transport and lysosomal processing, accompanied by suppression of de novo sterol biosynthesis"
- Important:
- - No manuscript figure was generated.
- - Complete enrichment outputs were retained; no pathways were manually deleted to improve appearance.
- - ORA background was the Step46 tested universe, not the whole genome.
- - GSEA ranked genes by moderated t, not logFC.
- - KEGG results depend on the KEGG resource available at run time.

## bulk_wgcna/49_WGCNA_modules.R

- Step49 v01 — Bulk WGCNA network construction and module detection
- input policy:
- - Primary WGCNA network uses Step45 expression_original.
- - Step45 batch-adjusted visualization matrix is NOT used for the primary network.
- - All 39 samples are retained.
- - Step45 review flags are retained as QC annotations only.
- - All 8 time points are retained.
- - No DEG-based gene preselection is used.
- - No biologically key module is selected in Step49.
- Input/network QC:
- Network definition:
- Module detection:
- core-gene module audit (NOT module selection):
- Interpretation boundary:
- - Step49 establishes network/module structure only.
- - Module-trait correlations, batch audit, target-state linkage,
- lipid/lysosome/sterol program linkage, key-module choice,
- and module enrichment are deferred to Step50.
- - The diagnostic PDF is QC output, not a manuscript figure.

## bulk_wgcna/50_WGCNA_module_trait_analysis.R

- Step50 v01 — Bulk WGCNA module-trait and cross-omics integration
- UPSTREAM:
- INTEGRITY:
- PRIMARY MODULE-TRAIT:
- Established late (14+28 dpi) vs Early, adjusted for operational batch.
- KEY-MODULE SELECTION:
- The module color was NOT hard-coded.
- Step48 mechanism themes were NOT used to choose the key module.
- STEP48 MECHANISM LINKAGE:
- Mechanism overlap uses significant GSEA core genes only.
- Positive themes: lipid handling/transport, lysosomal processing, myeloid activation.
- Negative theme: sterol biosynthesis.
- KEY-MODULE FUNCTIONAL ENRICHMENT:
- INTERPRETATION BOUNDARY:
- - Step50 integrates the WGCNA network with temporal and cross-omics evidence.
- - A key module is accepted only by the pre-specified statistical rule.
- - Step48 mechanism overlap is confirmatory, not part of module selection.
- - The selected module enrichment uses the WGCNA-tested gene universe.
- - No manuscript figure is generated here.
- - BayesPrism / SCIPAC remain later steps and are not run here.

## cellchat/77_CellChat_pooled.R

- Step77 | Global pooled Clec7a-centered CellChat
- Version: v05
- ANALYSIS DESIGN
- All source datasets pooled.
- All available time points pooled.
- No timepoint-specific CellChat objects.
- No dataset-specific CellChat objects.
- No temporal differential communication analysis.
- No downsampling.
- COMMUNICATION GROUPS
- Neutrophils remain a separate major lineage.
- Other cells retain finalized major-cell-type labels.
- CELLCHAT SETTINGS
- Species: mouse
- Database: CellChatDB v2
- Included: Secreted Signaling + ECM-Receptor + Cell-Cell Contact
- Excluded: Non-protein Signaling
- Average-expression method: triMean
- population.size = FALSE
- INPUT SIZE
- GLOBAL NETWORK
- CLEC7A-CENTERED NETWORK
- INTERPRETATION BOUNDARY
- These results represent a pooled/global potential communication landscape.
- They do not demonstrate time-specific changes in communication.
- Because cells from different datasets/time points are pooled, a predicted
- ligand-receptor interaction should not be interpreted as proof that the
- sender and receiver states coexist at one specific injury stage.
- Key axes should subsequently be prioritized with biological context,
- NicheNet and/or spatial plausibility where appropriate.
- FINAL QC
- Review Clec7a outgoing/incoming pathway and LR tables first.
- Then select candidate signaling axes for Figure 5 and NicheNet.

## external_cohort/115_input_check.R

- 115 - GSE281205 input audit |
- Sample-level metadata:
- Sample matrix summary:
- INTERPRETATION AND LIMITATIONS:
- 1. SOFT/Series Matrix values are original GEO metadata, not measured expression.
- 2. GSM8614665-4667 have a GEO treatment-title age conflict (Young vs Old).
- This script does not silently resolve it. Check BioSample and paper records
- before assigning these libraries to a young-animal analysis.
- 3. 10x barcodes may include filtered cells or droplets depending on upstream
- processing; raw barcode count is not necessarily a QC-passed cell count.
- 4. Raw Clec7a marker detection in all barcodes is NOT identification of
- Clec7a+ microglia. Cell-state assignment requires independent annotation
- and contamination/doublet assessment in a subsequent analysis.
- 5. Two cellular comparators (other microglia/all myeloid) CANNOT be counted
- from this technical audit. Do not claim pseudobulk feasibility yet.
- 6. GEO 'biological replicate' labels do not by themselves establish detailed
- animal pooling, independence, or matched injury-model comparability.
- 7. This step does not read the discovery 1,320-gene list and does not
- modify any signature genes or use validation-derived feature selection.
- NO ZIP, NO RDS, NO clustering, NO DEG, NO pseudobulk, NO P/FDR.
- No input source files were modified.
- Next step only if QC and metadata audit pass:
- Independently annotate myeloid/microglial states in the six libraries,
- confirm at least three eligible biological sample/library units at
- 14 dpi and 28 dpi for both comparison groups, and only then evaluate
- the unchanged frozen 1,320-gene signature.

## external_cohort/116_cell_qc.R

- 116 v02 | GSE281205 | CELL QUALITY CONTROL ONLY
- METHODS: retain >=300 detected genes and <10% mitochondrial transcripts.
- No upper nFeature/nCount limit; counts are unmodified raw UMI values.
- NO DecontX, doublet removal, clustering, integration, myeloid annotation,
- NO frozen 1,320-gene signature, pseudobulk, DE, or statistical testing.
- SAMPLE-LEVEL RESULTS:
- STAGE-LEVEL RESULTS:
- METADATA NOTE:
- GSM8614665--667 GEO title=Young-28dpi, GEO treatment=Old-28dpi.
- Preserved conflicting fields; intended age labels remain PROVISIONAL.
- Check primary BioSample and study source before definitive age inference.
- INTERPRETATION:
- QC-passed barcodes are NOT yet verified single cells: doublets and
- ambient RNA have not been removed. QC status alone is not proof of
- Clec7a+ microglia or availability of target/control groups.
- Six GSM libraries are not automatically six independent mice.
- No statement on signature validation is possible from step 116.
- NEXT STEP: ambient RNA (DecontX), then sample-wise doublets,
- then independent myeloid annotation and specimen-level analyses.

## external_cohort/118_decontx.R

- 118 v01 | GSE281205 | DecontX AMBIENT RNA CORRECTION
- METHOD:
- Six independent per-GSM corrections of 116-v02 QC-passed raw UMI matrices.
- For each GSM, DecontX z=NULL uses its own broad heuristic clustering.
- No independent empty-droplet matrix is available: background=NULL.
- No stage pooling, original signature, myeloid-state annotation, or manual
- cluster selection was used in fitting. Zero-count gene rows were skipped
- for model fitting but restored as zero rows in 118 RDS outputs.
- Original integer counts remain UNMODIFIED in 116 RDS files; 118 corrected
- counts are fractional estimates and not raw UMI.
- DecontX estimates may under/over-correct; no contamination-based cells
- were filtered or classified as singlets by this script.
- SAMPLE RESULTS:
- STAGE SUMMARY (descriptive; not animal-level inferential statistics):
- AGE PROVENANCE:
- GSM8614665-667 GEO title=Young-28dpi but treatment=Old-28dpi.
- Conflicting fields are retained, young assignment remains provisional.
- Inspect contamination/marker trends for obvious over-correction and
- extreme contamination; then do per-GSM scDblFinder (a later script).
- DO NOT interpret DecontX heuristic groups as microglial cell states.
- DO NOT use corrected fractional counts for formal pseudobulk inference.
- Use original integer UMIs from 116 and independent final cell labels.
- No frozen 1,320-gene signature has been accessed or tested.

## external_cohort/119_doublet_detection.R

- 119 v02 | GSE281205 | scDblFinder DOUBLETS ON 116 ORIGINAL UMI
- METHOD:
- Per-GSM independent scDblFinder using unmodified original 116 integer UMI counts.
- clusters=NULL; dbr=NULL; dbr.per1k=0.008; SerialParam; fixed seeds per GSM.
- A 10x-style expected rate is a model assumption, not a measured doublet rate.
- The scDblFinder model does NOT use DecontX estimated fractional expression.
- DecontX contamination is linked by exact cell IDs only for cross-tabulation.
- Neither >=20%, >=50%, nor >=80% contamination is a direct removal rule.
- Singlet-only raw integer UMI and separate DecontX corrected-count RDS are saved.
- Original 116/118 RDS inputs are never overwritten.
- SAMPLE-LEVEL RESULTS:
- STAGE-LEVEL RESULTS:
- GSM8614667 high-contamination cells require later manual marker/QC review.
- No contamination-based culling was done, including for this sample.
- GSM8614665--667: GEO young/old treatment disagreement is retained in metadata.
- 28dpi age assignment remains provisional pending final provenance check.
- GSM libraries cannot automatically be assumed to be independent mice.
- scDblFinder cannot rule out all doublets (especially homotypic doublets).
- NO cell-lineage annotation, integration, 1,320-gene signature scoring,
- pseudobulk, differential gene expression, or statistical testing performed.
- NEXT: review singlet counts and contamination cross-tabs before lineage annotation.

## external_cohort/121_normalization_PCA.R

- STEP 121 v02 | GSE281205 | LOGNORMALIZE/PCA/GSM HARMONY/AUDIT
- ANALYSIS:
- Six 119-v02 GSM-specific singlet-only RDS pairs checked by exact gene/cell IDs.
- Exploratory expression assay = DecontX-corrected estimated (possibly fractional) counts.
- Normalization: Seurat LogNormalize, scale.factor=10000.
- Harmony applied to the 30 PCs, correcting group.by.vars='gsm' only (not dpi).
- RPCA, SCTransform, cell clustering, annotation, DE, or signature scoring NOT run.
- NO regression of stage, UMI depth, contamination, or mitochondrial percent in expression scaling.
- No new cell filtering based on contamination or projected coordinates.
- SAMPLE STATISTICS:
- BATCH ASSESSMENT:
- PCA and Harmony stage/GSM R-squared are cell-level DESCRIPTIVE summaries.
- Stage is a biological effect; it must not be labeled a technical batch.
- GSM effects within each injury stage may include biological heterogeneity.
- Sample mixing against random-label expectation is descriptive, not a hypothesis test.
- Harmony was applied with 'gsm' as the correction covariate.
- WARNING: GSM is NESTED within dpi; Harmony CAN attenuate TRUE time-related biology.
- Inspect original and Harmony reductions side-by-side before any cell annotation.
- IMPORTANT PROVENANCE:
- GSM8614665-667 GEO Young/Old treatment metadata disagreement remains unresolved
- in this audit; sample titles and original Step119 labels are preserved.
- Six GSM libraries do not by themselves prove six independent mice.
- Original INTEGER UMI counts are stored in the per-GSM Step119 raw-count RDS files.
- Never run raw-count pseudobulk differential expression using DecontX estimated counts.
- The original 1,320-gene signature has not been used in Step 121.
- The present results do not establish a Clec7a+ cell state or external replication.

## external_cohort/123_major_lineage_clustering.R

- STEP123 v03 | GSE281205 | WHOLE-CELL CLUSTERING RECOVERY
- v01 generated all cluster and marker tables; failure occurred only at a plot.
- Verified complete per-cell ID/GSM/stage/replicate/contamination mappings.
- Recovered all four v01 cluster assignments without reclustering.
- Rebuilt Harmony/PCA neighbor graphs from unchanged Step121 coordinates.
- Copied v01 tables and figures into parallel v03-named files (excluding regenerated stage plot).
- Regenerated the stage composition plot using a discrete x-axis.
- Main whole-atlas solution: Harmony r0.4 (cluster identities from v01).
- The nine-state myeloid-only re-clustering has NOT been performed.
- Broad-lineage suggestions are exploratory, not manually confirmed.
- Original frozen 1320-gene signature was NOT loaded.
- No pseudobulk or differential expression validation was performed.
- Original integer UMI counts remain in Step119 per-GSM RDS files.
- High-contamination cells were marked but NOT removed.
- GEO age metadata discrepancy at 28 dpi remains unresolved.
- GSM is nested within dpi; consider Harmony overcorrection.
- CLUSTER STABILITY (DESCRIPTIVE):
- QC HIGH CONTAMINATION:
- NEXT ACTION: manually review marker evidence, candidate labels, stage/GSM
- and high-contamination patterns before extracting myeloid cells.
- STEP123 v03 | GSE281205 | WHOLE-CELL CLUSTERING + PROVISIONAL LINEAGE
- ANALYSIS:
- Whole-atlas Harmony SNN clustering at resolutions 0.4, 0.6, 0.8.
- Main whole-atlas solution: resolution 0.4 (original manuscript reference).
- The myeloid-only nine-state re-clustering has NOT been performed.
- Marker scoring uses independent canonical broad-lineage genes, NOT the frozen 1,320-gene signature.
- All broad-lineage suggestions are exploratory and require manual review.
- Neither the Clec7a+ state nor any pseudobulk contrasts are defined here.
- RNA for clustering/markers: DecontX-estimated normalized expression.
- Original integer UMI counts remain stored in Step119 per-GSM RDS objects.
- The GSM/sample biological replicate is not replaced by individual cells.
- High-contamination cells were flagged, NOT removed in this step.
- WARNING: 28dpi GSM title vs treatment age metadata disagreement remains.
- WARNING: Harmony GSM variable nested within dpi; check biological signal preservation.
- CLUSTER COUNTS:
- NEXT ACTION:
- Manually review cluster marker heatmap, candidate lineage table, stage/GSM composition,
- and high-contamination flags before extracting myeloid cells or labeling target state.
- Do not use cell-level Wilcoxon P-values as biological replicate evidence.

## external_cohort/124_myeloid_candidate_extraction.R

- STEP124 v01 | GSE281205 | PROVISIONAL LINEAGE REVIEW + CANDIDATE MYELOID EXTRACTION
- PROVENANCE & GUARDRAILS:
- Validated all 42,361 cell IDs, GSMs, stages, and parent Harmony cluster assignments.
- Reused fixed 24 whole-atlas clusters: Harmony resolution 0.4 (Step123 v03).
- Cluster identity decisions are provisional until explicitly manually confirmed.
- The frozen 1,320-gene signature was NOT loaded or used to select cells.
- Clec7a expression alone NEVER establishes microglial lineage.
- Cluster 0/3 remain ambiguously myeloid; cluster 18 dendritic/APC-like.
- Clusters 8/19 are excluded from proposed main extraction but remain available for quality review.
- No high-contamination cell was removed solely on its DecontX contamination estimate.
- The Step123 whole atlas is not edited by this script.
- The candidate Seurat assay contains DecontX-estimated values, NOT original integer UMI.
- Original integer UMI files from Step119 are not read or modified; retain for pseudobulk.
- GSMs are biological LIBRARIES; confirm animal/pool independence before inferential replicates.
- GSM identity is nested in dpi, and Harmony can weaken time-dependent biology.
- 28-dpi GEO title/treatment age conflict must remain disclosed.
- COUNTS:
- PER SAMPLE:
- PER STAGE:
- ACTION REQUIRED BEFORE FORMAL STATE-LEVEL PSEUDOBULK:
- 1) Independently review broad lineage and suspected contaminants; record approved cluster decisions.
- 2) Myeloid-only normalization, HVGs, PCA, Harmony sensitivity, subclustering and marker review.
- 3) Define target vs other microglia vs all other myeloid without using frozen 1,320 signature.
- 4) Verify adequate target/comparator cells in each of 3 GSM replicates per stage.
- 5) Aggregate Step119 ORIGINAL integer UMI counts by GSM x reviewed cell state.
- No p-values from exploratory single-cell markers are treated as specimen-level inference.

## external_cohort/126_targeted_lineage_validation.R

- STEP126 v01 | GSE281205 | TARGETED LINEAGE REVIEW
- NEXT ACTION: Review 07_126_cluster_review_dossier_v01.csv,
- 04_126_pairwise_codetection_by_cluster_GSM_v01.csv and figures.
- Record genuinely confirmed labels in the ORIGINAL Step125:
- Use Step125 allowed lineage/gate/status values; do not overwrite
- the Step125 decisions file from Step126 automatically.
- INTERPRETATION LIMITS:
- - Marker co-detection is descriptive, not proof of ancestry or function.
- - DecontX-corrected expression and contamination estimates are modeled.
- - Cluster 13 mixed endothelial/mural signal may represent true heterogeneity
- rather than same-cell coexpression; assess the pair table by GSM.
- - A high estimated contamination score does not itself warrant cell removal.
- - Nonsignificant/unmeasured gene signal does not exclude marker-negative states.
- - No frozen 1320-gene signature, no cell-level DE, no pseudobulk.
- - Original integer counts for future pseudobulk remain at Step119.
- - GSM nested in 14/28 dpi and 28dpi GEO title/treatment age discordance persists.
- - Harmony may attenuate genuine stage-related biological differences.

## external_cohort/127_lineage_evidence_check.R

- STEP127 v01 | GSE281205 | cluster13 raw UMI and manual review aid
- ANALYSIS:
- - Raw UMI count >=1/2/5 compared for vascular lineage marker pairs.
- - Compare all retained singlets with contamination<20% sensitivity subset.
- - scDblFinder scores are from Step119; Step119 already removed predicted doublets.
- - Neither co-detection nor low doublet scores proves or excludes heterotypic doublets.
- - No frozen 1,320-gene signature used; no pseudobulk or DE performed.
- MANUAL ANNOTATION POLICY:
- - Step125 original manual decision CSV is left unchanged.
- - Step127 proposed labels are suggestions, not review_status=CONFIRMED.
- - Review ALL 24 clusters; preserve 8/19 as unresolved until warranted.
- - Cluster13 remains Vascular_ambiguous / HOLD until evidence adjudicated.
- - Copy *human-edited* decisions to Step125 ORIGINAL 04_125_manual_decisions_v01.csv
- and rerun Step125 (NOT Step127) for final cell-level lineage assignment.
- LIMITATIONS:
- - Mouse lineage cannot be causally assigned from RNA markers alone.
- - Six GSM libraries are nested in 14 vs 28 dpi; the 28 dpi GEO age metadata conflict persists.
- - Harmony may attenuate biological stage differences.
- - This targeted audit only studies clusters 5, 10 and 13 at original-UMI level.
- - No automatic cell removal; all excluded/unresolved candidates remain traceable.

## external_cohort/128_lineage_finalization.R

- STEP128 v02 | GSE281205 | evidence-based CURATED WORKING LINEAGE ANNOTATION
- RESULT: existing Steps 125-127 recommendations adopted directly.
- There is no further 24-row manual-signoff gate or PREPARE-only pass.
- The user authorized moving ahead with the evidence-supported annotations.
- Do NOT call these 24 labels individually human-confirmed/signed-off.
- All 24 cluster gate decisions are now assigned as curated working decisions.
- INCLUDE clusters: 0,1,2,3,6,18,22
- HOLD clusters: 8,13,19 (explicitly reviewed as unresolved; not biological cell identities)
- CHECK: Step124 expanded-myeloid cell list matches Step128 INCLUDE exactly.
- The Step124 expression RDS was NOT loaded, altered or duplicated by this script.
- The Step125 manual_decisions file is not modified.
- NEXT: myeloid-only renormalization, PCA and reclustering.
- Use Step124 existing candidate Seurat RDS and Step128 curated cell IDs.
- For pseudobulk, use Step119 original integer UMI counts per GSM, NOT DecontX/Harmony.
- This script did NOT annotate the target Clec7a+ microglial substate.
- This script did NOT load the 1320-gene discovery signature or perform DE.
- WARNING: GSM is nested within dpi; Harmony may remove biological stage signal.
- WARNING: 28dpi GEO title/treatment age metadata conflict persists.
- WARNING: curating broad lineages does NOT prove microglial ontogeny.

## external_cohort/130_myeloid_reclustering.R

- STEP130 v02 | GSE281205 | myeloid-only renormalization + PCA + Harmony + reclustering
- INPUT AND CONSISTENCY:
- 14dpi 3 GSMs, 28dpi 3 GSMs; labels preserved as separate specimen/library IDs.
- Step128 WORKING broad lineages and exact INCLUDE cell set checked.
- GEO Young/Old treatment discrepancy in the three 28dpi GSMs remains noted.
- METHOD:
- DecontX estimated counts -> LogNormalize(scale.factor=10000), 3000 HVGs,
- ScaleData HVGs -> PCA 30 PCs -> Harmony by GSM (theta=2) -> UMAP.
- Harmony SNN FindClusters resolutions 0.10, 0.15, 0.20; primary=0.10.
- Uncorrected PCA SNN 0.10 clustering retained for sensitivity comparison.
- IMPORTANT SCIENTIFIC LIMITATIONS:
- No automatic target Clec7a+ lipid-handling microglia annotation was made.
- Canonical marker expression and exploratory HVG rankings are descriptive.
- The frozen original 1,320-gene signature was NOT loaded, scored or reselected.
- GSMs are nested within injury stage; Harmony by GSM can reduce true 14/28dpi biology.
- Any stage composition differences are cell-level descriptive, not animal-level tests.
- DecontX counts can be fractional and are NOT suitable as raw UMI pseudobulk input.
- Future specimen-by-state pseudobulk MUST use Step119 integer singlet raw counts.
- Next step: annotation of myeloid states from INDEPENDENT canonical markers and
- per-GSM feasibility checks, then test 1320 genes with replicate-aware models.
- PER GSM:
- CLUSTER ARI:

## external_cohort/131_myeloid_state_annotation.R

- STEP131 v01 | GSE281205 WORKING transcriptomic state annotation
- Primary state counts (QC/marker-gated):
- Cross-resolution caution flag is DESCRIPTIVE; not automatically removed from primary.
- Contrast A: target vs P2ry12-high other microglia, matched stage, 3 GSM/dpi.
- Contrast B: target vs all other HIGH-CONFIDENCE myeloid (including non-target microglia).
- Contrast C: target vs non-microglial myeloid only (secondary).
- Additional broad and cross-resolution partition sensitivities in manifest.
- These working states DO NOT use the 1,320 gene list for annotation.
- CAUTION: unsupervised variable gene selection could share genes with signature.
- No specimen-level expression inference is performed in Step131.
- Step132 pseudobulk MUST use original integer counts from Step119; not DecontX.
- GSM is nested within stage and was harmonized in Step130; evaluate uncorrected PCA sensitivity.
- 28dpi GEO age metadata discrepancy is not resolved by this computational analysis.
- Step131 does not establish transcriptomic groups as discrete developmental lineages.
- Contrast/stage readiness:
- SUCCESS | Step131 v01

## external_cohort/132a_pseudobulk_signature_validation.R

- Step132 v01 | GSE281205 ORIGINAL UMI, specimen-paired signature validation
- GSM: GSM8614662-4664 (14 dpi), GSM8614665-4667 (28 dpi)
- Primary model per comparison x stage: edgeR TMM QL, design ~ gsm + group.
- groupOTHER = baseline; groupTARGET coefficient positive means TARGET higher.
- Each stage: 3 paired GSM, 6 pseudobulk group libraries, residual df=2.
- GSM denotes specimen/library unit, not proven individual animals.
- 14 and 28 dpi analyzed separately. Never pool stages in a GSM-adjusted model.
- GENES: edgeR FDR corrected genome-wide among genes passing filterByExpr.
- Unmapped/untested frozen signature genes remain in per-gene output with NA.
- SIGNATURE: per-GSM mean/median frozen gene target-minus-other logCPM and
- direction, with pre-specified marker-excluded signature sensitivity.
- CAMERA gene-set p values are exploratory; small n and gene correlation limit
- statistical strength. Do not present as definitive biological validation.
- Note: 28 dpi metadata age discrepancies need source-level provenance review.
- Identity comparisons are state-defined, not Clec7a-positive vs negative.
- Sensitivity: Step131 broad and partition-consistent cell groupings.
- NONE of the 1,320 signature genes were newly selected or optimized here.
- No ZIP inputs; no DecontX counts; no cell-level pseudoreplication.
- STAGE-SPECIFIC VALIDATION SUMMARY:
- Interpret stage-specific directional consistency alongside independent
- specimen count, mapping coverage, stage confounding, and group definition.

## external_cohort/132b_pseudobulk_signature_comparison.R

- Step132 v02 | GSE281205 original UMI pseudobulk | MGI alias + CAMERA repair
- Raw original integer UMI sums reused from Step132 v01: 12 pseudobulk RDS.
- Gene symbols were reconciled AFTER original UMIs had been aggregated.
- The discovery gene membership, ordering and discovered effects did not change.
- Paired edgeR design ~gsm+group, six pseudo libraries/stage and residual df=2.
- 14 and 28 dpi analyzed separately; biological replicate is GSM specimen.
- Do not present GSM as independent mice without checking source provenance.
- CAMERA: 12 x 2 pre-specified sets; PValue competitive test vs expressed genome.
- CAMERA BH FDR is computed across the 12 stage/comparison tests within each
- set, NOT across 1,320 genes. Exploratory due to n=3 GSM per stage.
- Expression markers used to classify cells are excluded from a sensitivity
- analysis; remaining gene-set membership is not reoptimized.
- Any unmatched or ambiguous gene remains present with NA validation.
- Old v01 input/output files are never overwritten. NO ZIP dependencies.
- GSM and dpi are nested; age discrepancies at 28dpi remain a limitation.
- GENE COVERAGE BY DISCOVERY EFFECT RANK:
- STAGE SIGNATURE SUMMARY:
- CAMERA RESULTS:

## external_cohort/142_lipid_composition_robustness.R

- STEP142 v01 | GSE281205 GSM-level robustness
- No ZIP/archive input. Step20 gene sets, Step131 compositions,
- and the TWO locally stored Step132 contrast-A pseudobulk RDS checkpoints.
- Independent animal IDs are not established by this script; unit is GSM.
- COMPOSITION (per-GSM fractions; 3 GSM per stage):
- LIPID PROGRAMS:
- SCIENTIFIC INTERPRETATION
- Lipid score = mean paired TARGET-OTHER log2 TMM CPM of mapped genes.
- It is NOT an AUCell activity score; do not equate these metrics.
- No spot-level or cell-level statistical tests; every plotted replicate is a GSM.
- Stage comparisons are exploratory: n=3 per 14dpi / 28dpi; possible batch
- and 28dpi metadata-age issues must be acknowledged.
- Exact two-sided sample-label permutation with 3+3 has minimum p=0.10.
- Exact two-sided sign flip with 3 GSM has minimum p=0.25.
- Do not claim non-significant tests demonstrate no biology.
- Do not infer 0-28 dpi time-course or within-animal longitudinal changes.
- GSM identity does not establish a distinct animal; verify public metadata.
- Step131 states remain working transcriptomic labels; cell recovery varies.
- Marker-excluded scores are sensitivity checks, not separate discovery sets.
- 14 vs 28 dpi contrast reflects observed specimens, not causal stage effect.
- RECOMMENDED DISPLAY
- New Fig 142 lipid-program GSM effects and composition denominator figure.
- CSV 03/04: sample and stage lipid programs. CSV 06/07/08: composition.
- Existing S4D/S4E and Table S19 remain separate prior evidence.
- Assign final journal supplementary figure/table numbers at submission.

## figures/fig3/52_fig3B_bulk_stage_volcanoes.R

- Figure 3B — Bulk multi-volcano temporal DEG landscape
- source:
- Comparisons:
- DEG rule:
- Effect-significant = adj.P.Val < 0.05 AND |logFC| >= 0.5
- Plot rule:
- Only genes with |logFC| >= 0.5 are displayed, following the supplied multi-volcano visual concept.
- Up = red; Down = blue; effect-large but FDR-nonsignificant = light gray.
- DEG counts:
- Important interpretation note:
- 15 min vs Uninjured and 1 dpi vs Uninjured are exploratory because of the previously operational-batch imbalance.
- The figure is intended to summarize the temporal DEG landscape, not to establish those two earliest contrasts as independent causal effects.
- Formatting:
- Arial; no panel letter B; no internal title.

## figures/fig3/53_fig3C_bulk_primary_volcano.R

- Figure 3C — PRIMARY Established late vs Early volcano
- input:
- Comparison:
- Established late (14 + 28 dpi) vs Early (Uninjured + 15 min + 1 + 3 dpi)
- thresholds:
- adj.P.Val < 0.05 AND |logFC| >= 0.5
- Pre-specified lipid-handling labels:
- Pre-specified sterol-biosynthesis labels:
- Label rule:
- A pre-specified gene is labeled only when present and effect-significant in the Primary table.
- X = log2FC
- Y = -log10(BH-adjusted P value)
- Up = red; Down = blue; Not significant = gray.
- Dashed vertical lines = log2FC +/- 0.5.
- Dashed horizontal line = FDR 0.05.
- No panel letter C and no internal title.
- Interpretation boundary:
- This panel visualizes the Primary comparison.
- It does not rerun or redefine differential expression.

## figures/fig3/54_fig3D_Clec7a_signature_GSEA.R

- Figure 3D — scRNA-derived Clec7a+ lipid-handling signature GSEA
- source:
- Primary contrast:
- Established late (14+28 dpi) vs Early
- Ranking metric:
- limma moderated t statistic from Step46 Primary comparison.
- Plot construction:
- Top = weighted running enrichment score.
- Middle = signature hit positions; frozen leading-edge genes highlighted.
- Bottom = moderated-t ranking distribution.
- Left side corresponds to genes increased in Established late.
- Right side corresponds to genes increased in Early.
- Important:
- No fgseaMultilevel call is made in this figure script.
- NES / ES / P / FDR are read directly from the Step47 result.
- The enrichment curve is reconstructed only for visualization.
- No panel letter D and no internal plot title.

## figures/fig3/55_fig3E_bulk_GO.R

- Figure 3E — GO BP ORA barplot (manual pathway set; v04)
- source files:
- Display:
- Top panel = Down
- Bottom panel = Up
- x = -log10(pvalue)
- fill = -log10(p.adjust)
- pathway labels are placed inside bars
- bars are ordered from longest to shortest within each panel
- Manual DOWN pathways:
- Manual UP pathways:
- Important:
- This script does not rerun GO enrichment.
- It only reads and plots the selected terms from the Step48 GO BP results.

## figures/fig3/56_fig3F_bulk_GSEA.R

- Figure 3F — GO BP GSEA ridgeline
- sources:
- GSEA ranking:
- limma moderated t statistic.
- Ridgeline x-axis:
- gene-level log2FC from the Step46 Primary comparison.
- Important interpretation:
- The ridges visualize the log2FC distribution of ALL mapped tested genes
- belonging to each selected GO gene set.
- The ridges are not leading-edge-only distributions.
- NES, P and FDR are read directly from Step48 GSEA results.
- GSEA is not recomputed from logFC.
- Plot encoding:
- ridge fill = -log10(BH FDR)
- left dot color = NES
- red NES = Established late enrichment
- blue NES = Early enrichment
- vertical dashed line = log2FC 0
- v02 visual refinement:
- geom_density_ridges scale = 0.84
- rel_min_height = 0 to preserve continuous KDE tails and avoid clipped density islands
- No pathway/statistical changes relative to v01.
- Selected positive pathways:
- Selected negative pathways:
- No panel letter F and no internal title.

## figures/fig3/57_fig3G_WGCNA_module_trait_heatmap.R

- Figure 3G v05 — Reference-style cross-omics WGCNA module relationships
- Frozen Step18 state positive-marker signatures x frozen WGCNA modules.
- Two-sided Fisher enrichment/depletion; BH FDR across 11 modules within each state.
- Fill = signed log2(odds ratio) only for FDR<0.05 cells.
- FDR<0.05 & OR>1 = red enrichment; FDR<0.05 & OR<1 = blue depletion; non-significant cells = neutral light gray.
- Blue cells indicate significant depletion, not negative correlation.
- LEFT TWO-SIDED FISHER SUMMARY:
- CENTER:
- WGCNA module colors and ME labels.
- RIGHT:
- Step50 batch-adjusted temporal module associations.
- Fill = target-minus-reference ME estimate.
- V05 LAYOUT CHANGE:
- Center WGCNA-module panel widened to ensure complete rendering of long labels,
- especially MEgreenyellow. No data/statistical changes from v04.
- Interpretation:
- Left and right use different statistical quantities and therefore separate legends.
- This supports cross-omics signature linkage, not direct cell-proportion inference.
- No deconvolution was performed.
- No scRNA marker analysis, WGCNA reconstruction, temporal limma rerun, or module reselection.

## figures/fig3/58_fig3H_brown_module_GO.R

- Figure 3H v03 — Brown-module GO BP barplot
- INPUTS:
- PLOT DESIGN:
- - Horizontal barplot.
- - X-axis = Number of Gene (Count).
- - Fill color = -log10(p.adjust).
- - GO term labels start inside each horizontal bar.
- - A small plot title is included.
- For short bars, long labels are allowed to extend beyond the bar endpoint.
- A panel border box is included.
- No panel letter H.
- SELECTION PRINCIPLE:
- - The brown-module GO BP ORA table contains many significant terms.
- - Terms were manually selected to represent:
- 1) lipid handling,
- 2) lysosome / uptake / phagocytosis,
- 3) myeloid activation,
- 4) tissue repair / remodeling.
- - ORA was NOT rerun.
- - WGCNA was NOT rerun.
- SELECTED TERMS (ordered in the final plot from long to short bars):
- FOCUS CORE GENE SUPPORT:

## figures/fig4/69_fig4A_spatial_regions.R

- Figure 4A / Step 69
- Biologically concordant representative spatial transcriptional-region maps
- SELECTION PRINCIPLE
- 3 dpi:
- choose the section whose Clec7a section-level mean is closest
- to the 3 dpi timepoint median.
- 7/14/28 dpi:
- section must contain the Myeloid/immune–lipid-handling enriched
- region, must have >= minimum target-region spots, and must show
- positive direct enrichment for that region.
- among eligible sections, choose the one closest to the
- timepoint median Clec7a section-level mean.
- This avoids choosing the strongest outlier section.
- VISUAL STYLE
- Figure4A v11 intentionally returns to the v10 visual design.
- Original Step64 array_x / array_y coordinates and v10 crop windows are retained.
- No recentering, rescaling, stretching, rotation or interpolation is used.
- The v10 larger solid-square spots, darker saturated palette, compact 2 x 2
- arrangement and centered right-side legend area are retained.
- v12 retains the v11 explicit color-key + text legend, but compresses
- the vertical spacing so all 10 legend entries form a compact centered block.
- No H&E/pixel registration is claimed.
- SELECTED REPRESENTATIVE SECTIONS
- OUTPUT FILES
- FINAL QC

## figures/fig4/70_fig4B_Clec7a_RCTD_weights.R

- Figure 4B / Step 70
- Clec7a spatial abundance map
- DISPLAY MODE
- signal_plot = raw normalized RCTD weight
- One shared scale is used across all four representative sections.
- Display lower limit = 0.
- Display upper limit = q99 across all plotted spots.
- Values above q99 are squished to the maximum color.
- Color palette = custom darker inferno-like gradient
- (black -> deep purple -> dark magenta -> orange -> deep yellow)
- VISUAL ITERATION (v08)
- Biological content, selected sections, raw RCTD signal and q99 limits are unchanged.
- The continuous palette remains the darker, saturated v07 version.
- v08 increases spot size from 1.28 to 1.60 for denser visual coverage.
- Crop geometry, panel layout, q99 scale and white facet-strip backgrounds are unchanged.
- dplyr::select() remains explicitly namespaced to avoid package-function masking.
- REPRESENTATIVE SECTIONS
- OUTPUT FILES
- FINAL QC

## figures/fig4/71_fig4C_niche_burden.R

- Figure 4C / Step 71
- Clec7a spatial-niche burden heatmap
- BIOLOGICAL QUESTION
- Which transcriptionally defined spatial niches contribute most
- to the tissue-level Clec7a+ lipid-handling Microglia signal
- over time?
- SECTION-LEVEL METRIC
- spatial_niche_burden =
- n_region_spots * mean_region_target_weight
- / total_section_spots
- This integrates both region occupancy and local Clec7a RCTD signal.
- If a region is absent from a section, its burden is zero.
- TIMEPOINT SUMMARY
- Each heatmap cell is the arithmetic mean of section-level burdens
- for that timepoint x spatial region.
- Sections are not treated as confirmed independent animals.
- DISPLAY
- Columns: 3 dpi, 7 dpi, 14 dpi, 28 dpi.
- Rows: fixed biological spatial-region order.
- Fill: raw section-normalized spatial niche burden.
- Palette: pale cream -> orange -> red.
- Pale-cream anchor = 0.
- Orange anchor = q70 of raw burden values.
- Red anchor = q95 of raw burden values.
- Values above q95 are visually squished to deepest red.
- Raw numeric cell labels are unchanged.
- Colors represent lower -> higher spatial burden only.
- No row z-score.
- No clustering / dendrogram.
- No hypothesis testing.
- VISUAL ITERATION (v04)
- Only the heatmap palette is changed to match the supplied reference:
- pale cream -> orange -> red.
- Burden values, q70/q95 display anchors, row/column order, and labels are unchanged.
- TARGET NICHE
- INTERPRETATION LIMIT
- This is a descriptive tissue-level spatial burden summary.
- It does not establish independent-animal inference.
- Spatial regions are transcriptionally defined and are not
- H&E-defined anatomical lesion compartments.
- FINAL QC

## figures/fig4/72_fig4D_MoransI.R

- Figure 4D / Step 72
- Clec7a section-aware Moran's I summary
- BIOLOGICAL QUESTION
- Is the Clec7a+ lipid-handling Microglia spatial signal
- non-randomly organized across individual tissue sections?
- PLOT DESIGN
- One semi-transparent point = one coordinate-defined major tissue section.
- Large filled point = timepoint median Moran's I.
- The line connects TIMEPOINT MEDIANS only and is descriptive.
- It does not indicate paired or longitudinal tracking of sections.
- Plot background is white; no positive-zone shading is used.
- Dashed horizontal line = Moran's I = 0 reference.
- No violin / boxplot / significance stars.
- No between-timepoint hypothesis testing.
- DESCRIPTIVE SUMMARY
- INTERPRETATION LIMIT
- This panel is descriptive across injury timepoints.
- The coordinate-defined sections are not assumed to be
- independent biological replicates or independent animals.
- Do not claim that 14 dpi is statistically higher than 3/7/28 dpi
- from this panel.
- FINAL QC

## figures/fig4/75_fig4F_Clec7a_lipid_14dpi.R

- Figure4F / Step75
- 14 dpi Clec7a-lipid density scatter
- DISPLAY
- Each point = one eligible 14 dpi Visium spot.
- X = Clec7a+ lipid-handling Microglia abundance (Clec7a_RCTDnorm).
- Y = lipid-handling program score (score_Core_lipid_handling).
- Point color = local point density via densCols.
- Black line = linear trend for visualization only.
- Only Spearman rho is shown in the panel, at the top-right.
- Spot-level P and n are NOT shown in the main figure.
- Figure4F and Figure4G use the SAME X/Y axis limits.
- DESCRIPTIVE SUMMARY
- INTERPRETATION LIMIT
- This is a descriptive spot-level association plot.
- Visium spots are multicellular mixtures and not independent animals.
- The displayed rho is descriptive; spot-level P is not used as manuscript inference.
- Therefore the panel supports association, not causal inference.
- FINAL QC

## figures/fig4/76_fig4G_Clec7a_lipid_28dpi.R

- Figure4G / Step76
- 28 dpi Clec7a-lipid density scatter
- DISPLAY
- Each point = one eligible 28 dpi Visium spot.
- X = Clec7a+ lipid-handling Microglia abundance (Clec7a_RCTDnorm).
- Y = lipid-handling program score (score_Core_lipid_handling).
- Point color = local point density via densCols.
- Black line = linear trend for visualization only.
- Only Spearman rho is shown in the panel, at the top-right.
- Spot-level P and n are NOT shown in the main figure.
- Figure4F and Figure4G use the SAME X/Y axis limits.
- DESCRIPTIVE SUMMARY
- INTERPRETATION LIMIT
- This is a descriptive spot-level association plot.
- Visium spots are multicellular mixtures and not independent animals.
- The displayed rho is descriptive; spot-level P is not used as manuscript inference.
- Therefore the panel supports association, not causal inference.
- FINAL QC

## figures/fig6/100_fig6F_virtual_knockout.R

- Step 100 v03: final main-text Figure 6F (2 x 3 layout)
- ANALYSIS POLICY
- No scTenifoldKnk analysis was rerun.
- Step 91 v04 results were used directly.
- DISPLAY RULE
- 1. Keep the KO TF itself as row 1.
- 2. Remove only low-information gene-symbol categories:
- mitochondrial (^mt-), ribosomal (^Rpl/^Rps), predicted Gm,
- Rik-containing clone symbols, pseudogene-style -ps/ps#, and Rn###s rRNA-style genes.
- 3. Rank all remaining genes only by Step 91 perturbation_rank.
- 4. Keep top 19 non-KO genes.
- 5. KO TF + 19 non-KO genes = 20 displayed genes per panel.
- NO STORY-DRIVEN MANUAL SELECTION
- No Clec7a, lipid-handling, or signature gene was manually added.
- FDR was not used to decide display membership.
- BIOLOGICAL AUDIT OF DISPLAYED NON-KO GENES
- IMPORTANT INTERPRETATION
- log2_perturbation_FC is a perturbation-distance fold change.
- It is NOT an RNA-expression fold change and does NOT encode expression up/down direction.
- OUTPUTS

## figures/fig6/103_fig6G_Zfp148_GO_GSEA.R

- Step 103 v03: Zfp148 Figure 6G GO-GSEA ridgeplot
- SELECTION POLICY:
- 16 FDR<0.05 terms were selected for SCI/lipid biological relevance.
- This is a representative, theme-driven display, not the statistically top 16 terms.
- All NES/FDR are from Step102 Zfp148 PRIMARY filtered GO-GSEA.
- No scTenifoldKnk KO or GO-GSEA rerun.
- v03 changes are display-only: readable labels, four category gaps,
- vertical separators, legend spacing, and export dimensions.
- All full GO terms are preserved in Description; LDL is label abbreviation.
- RIDGE POLICY:
- Each density contains full GOALL membership intersected with PRIMARY ranked genes.
- Leading-edge/core enrichment genes alone are NOT used.
- Annotation version can affect GO memberships; differences against setSize are audited.
- INTERPRETATION:
- Perturbation Z and NES reflect virtual-KO network perturbation ranking.
- They do NOT imply actual up/downregulated RNA, increased lipid clearance,
- nor direct Zfp148 regulation of Lpl/Npc2/Abca1.
- SELECTED GO TERMS (top-to-bottom):

## figures/fig6/95_fig6A_regulon_RSS.R

- Step 95: Figure 6A RSS specificity ranking
- STATUS: PASS
- REFERENCE CODE LOGIC
- All regulons are ranked by specificity score.
- Top regulons are colored red and labeled; remaining regulons are gray.
- X-axis = Regulons Rank; Y-axis = Specificity Score.
- PROJECT ADAPTATION
- IMPORTANT
- RSS was NOT recalculated.
- The figure reads the Step89 RSS ranking directly.

## figures/fig6/96_fig6B_regulon_activity_umap.R

- Step 96: Figure 6B top10 regulon activity UMAP
- STATUS: PASS
- REFERENCE-CODE ADAPTATION
- The uploaded reference script uses FeaturePlot on the continuous Regulon assay.
- This project uses the Step88 AUCell matrix directly on the myeloid UMAP.
- LAYOUT
- 2 rows x 5 columns.
- COLOR SCALE
- Values above the pooled 99th percentile are visually squished only; underlying AUC values are unchanged.
- v03 color direction follows the reference: purple=low, blue/green=intermediate, yellow=high.
- IMPORTANT
- No RSS, regulon, or AUCell score was recalculated.
- Figure6B is a publication-only visualization of Step88/89 results.

## figures/fig6/97_fig6C_regulon_heatmap.R

- Step 97 v03: Figure 6C state-associated regulon heatmap
- REFERENCE ADAPTATION
- The reference script selects top regulons per cell subset and then displays the union in a scaled heatmap.
- This project uses state-mean AUCell values to perform the equivalent publication-stage selection.
- SELECTION
- TARGET STATE TOP10 BY RELATIVE ACTIVITY
- DISPLAY
- Rows and columns are both clustered.
- Blue-white-red shows row-Z-scored relative regulon activity.
- Figure6A RSS top10 row labels are red; other regulons are black.
- Figure-only state abbreviations:
- IMPORTANT
- No pySCENIC, RSS, or AUCell calculation was rerun.

## figures/fig6/98_fig6D_TF_target_network.R

- Step 98 v05: Figure 6D Cytoscape GraphML
- STATUS: PASS
- PURPOSE
- Avoid Cytoscape node-table key matching problems by embedding
- all node and edge attributes directly in a single GraphML network.
v03 also fixes the v02 igraph vertex-key bug by forcing gene ID
to be the first vertices-table column.
- IMPORTANT
- No TF/target selection was rerun.
- The Step98 v01 node/edge tables are used unchanged.
- CYTOSCAPE IMPORT
- File -> Import -> Network from File...
- After import, open Node Table and confirm:
- Mef2a -> node_type = TF; regulon_size = 847
- Qk    -> node_type = Target
- Nr3c1 -> node_type = TF; regulon_size = 881
- Then Style mappings can use node_type, regulon_size, log2FC,
- display_name, cy_node_size, cy_fill_hex_final, and edge attributes directly.
- FINAL FIGURE 6D FILL COLOR
- CYTOSCAPE FILL COLOR SETTING
- Style -> Node -> Fill Color
- Column = cy_fill_hex_final
- Mapping Type = Passthrough Mapping

## figures/fig6/99_fig6E_TF_gene_correlations.R

- Step 99 v01: Figure 6E reference-style correlation heatmap
- ANALYSIS POLICY
- No correlation was recalculated.
- The Step 90 v05 target-state correlation table was used directly.
- PANEL DEFINITION
- Rows: 20 key functional genes.
- Columns: top 10 RSS-specific TFs.
- Fill: Spearman rho.
- Stars: BH-adjusted FDR.
- SELECTED FIVE TFs
- OUTPUTS

## regulatory/102_Zfp148_GO_GSEA.R

- Step 102 | Zfp148 full-rank GO-GSEA v01
- Step 91 virtual knockouts are NOT recalculated.
- Rank: signed scTenifoldKnk perturbation Z; KO TF itself excluded.
- Primary: interpretable gene-filtered; sensitivity: unfiltered mapped genes.
- Annotation: org.Mm.eg.db; symbols mapped to Entrez with multiVals='first'.
- Multiple testing: BH within each ontology and within each TF/analysis.
- Z/NES DO NOT INDICATE GENE UPREGULATION, DOWNREGULATION, PATHWAY ACTIVATION OR REPRESSION.
- GO member-gene Z ridges are descriptive density plots, not additional statistical tests.
- The Step 91 5000-gene network forcibly preserved selected genes; interpret GSEA cautiously.
- The same algorithm is run anew for BOTH TFs for fair Zfp148-vs-Mef2a comparisons.
- This is not an exact replay of Step 101 plotting/selection without its source code.
- --- Counts ---
- --- Versions ---

## regulatory/84_pySCENIC_input.R

- Step 84: pySCENIC input preparation / QC
- PURPOSE
- QC-only preparation for Figure 6 pySCENIC analysis.
- No normalization, integration, clustering, annotation, GRN inference,
- motif enrichment, regulon inference or AUCell analysis was rerun.
- OBJECT
- TARGET STATE
- GENE IDENTIFIERS
- FINAL QC
- NEXT STEP

## regulatory/85_GRNBoost2_input.R

- Step 85: pySCENIC GRNBoost2 input preparation
- ANALYSIS STATUS
- MYELOID ATLAS
- GRN DISCOVERY SAMPLING
- IMPORTANT SAMPLING PRINCIPLE
- Time points were not forced to equal sample sizes.
- GENE FILTERING
- MATRIX
- EXPORTED GRNBOOST2 INPUT
- MATRIX ORIENTATION
- The exported Matrix Market matrix is genes x cells.
- Step 86 must read/transpose explicitly for pySCENIC/Arboreto as required.
- WHAT WAS NOT DONE
- No normalization was performed.
- No clustering or annotation was rerun.
- No RDS was modified.
- GRNBoost2 was NOT run.
- cisTarget was NOT run.
- Regulon AUCell was NOT run.
- NEXT STEP

## regulatory/87_cisTarget_pruning.R

- Step 87: pySCENIC cisTarget motif enrichment / pruning
- METHOD
- Command: pyscenic ctx
- OUTPUT

## regulatory/88_regulon_AUCell.R

- Step 88: pySCENIC regulon AUCell on full myeloid atlas
- METHOD
- Command: pyscenic aucell
- OUTPUT

## regulatory/89_regulon_RSS.R

- Step 89: pySCENIC RSS / target-specific regulon ranking
- METHOD
- Function: pyscenic.rss.regulon_specificity_scores
- Metric: Jensen-Shannon-distance-based RSS
- RSS OUTPUT
- PRIMARY OUTPUTS
- INTERPRETATION RULE
- RSS is the primary specificity metric.
- Mean-AUC contrasts are auxiliary descriptors and are not used to preselect biological conclusions.

## regulatory/90a_TF_target_network.R

- Step 90: top 10 TF network + target-state signature correlation
- FROZEN TARGET-STATE SIGNATURE
- Definition: p_val_adj < 0.05 and avg_log2FC >= 0.25
- Confirmed by Step 47 signature audit.
- TOP 10 RSS TFs
- FIGURE 6D
- TF node size = full regulon size.
- Target node color = target-state avg_log2FC.
- Visualization targets are preferentially frozen UP-signature genes ranked by final Step 87 regulon weight.
- FIGURE 6E
- Expression normalization: log1p(raw_count / library_size * 10000).
- All-myeloid correlation is retained only as a sensitivity analysis.
- TOP-5 POLICY
- No top 5 TFs were automatically declared.
- The evidence table integrates RSS, regulon/network structure, signature-target coverage, and TF-signature correlations for inspection.

## regulatory/90b_TF_functional_correlations.R

- Step 90 v05: final Figure 6E functional-gene panel
- Figure 6D/network remains from Step 90 v03.
- Figure 6E was rebuilt using a predeclared 20-gene functional panel.
- FINAL 20 KEY FUNCTIONAL GENES
- SELECTION SAFEGUARDS
- All 20 genes belong to the frozen 1,320-gene Step 21 UP signature.
- Each gene has support from the Step 20 lipid manifest and/or selected enriched Step 22 GO-BP terms.
- No v05 TF-gene correlation was used to select genes.
- No top-5 TFs were preselected.
- CORRELATION
- Expression normalization: log1p(raw_count / library_size * 10000).
- DESCRIPTIVE V05 CORRELATION RANKING
- TOP-5 POLICY
- Do not select the five scTenifoldKnk TFs from correlation alone.
- Review RSS specificity, Figure 6D network/regulon evidence, and final Figure 6E correlation together.

## regulatory/91_virtual_knockout.R

- Step 91: scTenifoldKnk top-5 virtual knockout
- ENVIRONMENT
- KO TFS
- INPUT / QC
- Manual >=25-cell gene filtering was applied because scTenifoldKnk v1.0.3 does not assign its qc_minCells subset.
- NETWORK PARAMETERS
- KO SUMMARY
- IMPORTANT INTERPRETATION
- scTenifoldKnk FC is a fold change of perturbation distance relative to expectation.
- It is NOT an expression fold change and does not encode up- versus down-regulation.
- DESCRIPTIVE PERTURBATION RANKING

## sensitivity/bulk_wgcna/102_WGCNA_sensitivity.R

- Step102 v01 -
- WGCNA preprocessing / batch / outlier sensitivity analysis
- INPUTS:
- - No ZIP/archive input is used.
- NETWORK DESIGN:
- WGCNA SETTINGS:
- - TOM files are not saved by this sensitivity script.
- BATCH PROVENANCE:
- - Batch_A / Batch_B are operational labels assigned in Step45 from sample-ID prefixes.
- - Batch_A = jerry. / normal.; Batch_B = D-prefixed sample IDs.
- - Do not describe Batch_A / Batch_B as independently verified experimental batches without source/laboratory confirmation.
- BATCH-ADJUSTED MATRIX:
- - Step45 expression_batch_adjusted_visualization was generated by limma::removeBatchEffect with timepoint effects preserved.
- - It is used here only as a sensitivity preprocessing matrix.
- BROWN MODULE STABILITY:
- BROWN MODULE ASSOCIATIONS:
- INTERPRETATION LIMITS:
- - Module colors are not assumed to correspond across networks; modules are matched by gene overlap.
- - Jaccard, retention, eigengene correlation and association statistics must be interpreted together.
- - The batch-adjusted expression matrix was originally designated for visualization/QC in Step45; this is a sensitivity analysis, not proof that batch has been fully removed.
- - Batch_A / Batch_B are not independently confirmed experimental batches in the supplied Step44 metadata.
- - This script does not perform tissue deconvolution and cannot determine whether Brown reflects cell composition versus cell-intrinsic activation.
- - A stable Brown module does not by itself prove Clec7a+-specific regulation.
- OUTPUTS:
- No upstream input files are modified.

## sensitivity/bulk_wgcna/109_bulk_timepoint_sensitivity.R

- Step 109: Bulk RNA-seq temporal-stage sensitivity analyses
- group definitions:
- Five pre-specified sensitivity analyses:
- S1. Established late vs Early, excluding Uninjured.
- S2. Established late vs Early, excluding 15 min.
- S3. Established late vs Early, excluding Uninjured and 15 min.
- S4. 14 dpi vs Early.
- S5. 28 dpi vs Early.
- Model results and estimability:
- DEG count summary:
- Concordance against the full-sample Step46 primary contrast:
- Interpretation notes:
- 1. The reference primary model is re-fitted only to provide a consistent
- comparator; it is not counted as one of the five sensitivity analyses.
- 2. Excluding Uninjured and/or 15 min changes the composition and sample
- size of the Early reference group. Interpret these as robustness checks.
- 3. The 14 dpi and 28 dpi contrasts are each compared with the same Early
- group. Review the batch-by-group sample distribution and full-rank
- diagnostics before interpreting these contrasts.
- 4. NOT_ESTIMABLE models were not forced. Their reasons are recorded in
- tables/03 design_diagnostics_and_estimability.csv.
- 5. The expression scale, filtering, and normalization were not changed;
- the Step45 expression_original matrix was used throughout.
- No manuscript figure was generated in Step109.

## sensitivity/bulk_wgcna/111_myeloid_signal_adjustment.R

- Step111: Bulk myeloid-marker-score-adjusted signature enrichment
- Marker-score phase tests:
- Limma model diagnostics:
- Clec7a-UP signature GSEA:
- Gene-effect concordance:
- DEG count summary:
- Interpretation safeguards:
- 1. The Step46 continuous expression_original matrix was used unchanged.
- 2. The original Bulk CSV was not interpreted as raw counts.
- 3. No rounding, offset/shift, inverse transformation, additional log transform
- or DESeq2 analysis was performed.
- 4. The Myeloid9 and Myeloid10 covariates are marker enrichment scores,
- not RNA fractions, cell proportions or absolute cell counts.
- 5. Marker genes were derived from Step09/Step20 single-cell references
- and the frozen 1,320-gene Clec7a-UP signature was excluded from selection.
- 6. The adjusted limma models contain phase and marker score only; no batch
- variable was inferred from sample names or used.
- 7. Persistence of signature enrichment after marker-score adjustment supports
- robustness to this independent transcriptional covariate, but does not
- prove cell-intrinsic regulation or eliminate all tissue-composition effects.
- 8. Inspect marker gene rankings, cross-state prevalence, marker-score/phase
- association, VIF and design condition number before manuscript use.
- No Step09/20/21/46/49 source objects or previous outputs were modified.

## sensitivity/cellchat/139_stage_specific_CellChat.R

- Step 139 v03 | full late-stage CellChat / spatial evidence
- REPRODUCIBLE INPUTS
- DESIGN
- Step09 v04 whole atlas + Step19 v02 refined myeloid annotations.
- Exact cell-ID joins; all nine original states; Clec7a renamed only for display.
- 14dpi and 28dpi NEVER pooled into one CellChat probability model.
- No deliberate cell downsampling; groups < min.cells excluded only within each stage.
- CellChat link tables from subsetCommunication (usually significant links only).
- Cross-stage comparison of significant links is DESCRIPTIVE ONLY.
- CellChat outputs are in RDS/, figure PDFs in figures/, tables and warnings in tables/.
- Resume checkpoints in checkpoints_v03/; using FORCE_RESTART affects only this v03 step.
- FOUR PREDECLARED CANONICAL AXES
- Apoe -> Trem2 + Tyrobp; Psap -> Gpr37l1; Spp1 -> Cd44; App -> Trem2 + Tyrobp.
- All additional ligand==Apoe/Psap/Spp1/App LR are also output, but the spatial
- pair-specific test is restricted to these four canonical receptor definitions.
- Same-stage RNA: source group ligand detected fraction; target group ALL receptor
- subunits co-detected in same single cells; expression is NOT proof of interaction.
- SPATIAL TEST
- Step64 Spatial log-normalized expression + Step65 v06 section-aware RCTD weights.
- Exact (GSM,array_x,array_y) coordinate joins. Major eligible sections only.
- For each confirmed significant canonical LR, count source spots with ligand RNA
- plus source RCTD group; receiver spots with ALL receptor genes plus receiver group.
- Report same-spot and adjacent-spot compatibility by capture area, not as binding.
- All spot-level values are descriptive. No per-spot inferential p-values.
- RCTD weights and spot gene expression derive from the same transcriptomic assay;
- their agreement does not constitute independent molecular validation.
- GSE195783 has two capture areas per late stage; areas may not be independent animals.
- EXTERNAL COHORT
- GSE281205 Step131 working myeloid groups are NOT identical to Step19 states.
- The optional Step130 DecontX assay supplies corrected-count detection, not raw UMI.
- GSE281205 GSM-level expression is NOT replicated communication probabilities.
- UNRESOLVED BIOLOGICAL LIMITATIONS
- 14dpi CellChat only source Figshare 17702045; 28dpi only source GSE196928.
- Stage effect cannot be separated from source. Do not claim same-stage cross-source reproducibility.
- Step77 sample identifiers do not prove individual biological-animal independence.
- CellChat bootstraps do not measure independent-animal variability.
- A model-predicted Clec7a-to-Clec7a APOE edge does NOT validate autocrine activity.
- Expression and RCTD-based spatial coexistence do NOT demonstrate receptor binding
- or functional signal transduction and do not resolve anatomical ground-truth.
- No updated Figure 5 should be finalized before inspecting all new QC tables.

## sensitivity/regulatory/135_TF_candidate_comparison.R

- Step 135 v02 optimized
- INPUT: direct original files under the configured project root; NEVER reads a ZIP archive.
- upstream: Step20,84,88,89,90,91,92, and 135 v01 historical output.
- RSS / sample-label-level activity / target sets / lipid co-activity: recomputed identically to 135 v01.
- INDEPENDENT ANIMALS NOT VERIFIED: label bootstrap is descriptive ONLY.
- NO late-stage independent-animal replication established.
- Network: original 5000-gene WT imported read-only; new 3000-gene networks have nNet=10.
- Seed and cell80 scenarios have 3 independently initialized runs each.
- Paired sensitivity: remove forced 20 gene panel OR change q OR change K; alter ONE factor at a time.
- 3k B00 uses seed9101, matching Step91 original 5k to isolate gene-universe sensitivity.
- 3000 genes are selected from CP01 5000-gene input universe, NOT de novo transcriptome.
- All 10 TFs forced into each network so their KO can be evaluated fairly.
- Source-specific networks remain confounded with source-specific stage/cell recovery.
- Lipid fgsea uses Step20 pathway manifest, ranked scTenifoldKnk Z;
- Positive/negative Z and NES describe network distance ranking position, NOT RNA regulation direction.
- Cross-network comparisons use the common NON-KO gene universe for each TF.
- MA randomness: fixed TF-specific seeds matched to original 135 v01 S00.
- Lipid enrichment is a secondary descriptive computational association, NOT causal function.
- Network stability does NOT establish regulon/motif-pruning robustness.
- Multiple data-dependent decisions are exploratory; no pre-registered composite TF score.
- Do not claim Zfp148 dominance unless supported by all candidate comparisons.

## sensitivity/regulatory/136_myelin_term_annotation_check.R

- STEP 136 | Myelin-sheath GO annotation audit | v01
- DecontX raw-versus-corrected pairing status:
- Interpretation:
- 1. Enriched GO:0043209 is based on NETWORK PERTURBATION Z rankings, not KO-RNA fold change.
- 2. GO terms can share leading-edge genes; these are not independent functional assays.
- 3. Canonical myelin/oligodendrocyte-marker and phagocytosis-marker detection is descriptive.
- 4. Any paired DecontX comparison is a computational contamination-sensitivity check only.
- 5. No scRNA-seq count-only analysis can prove myelin ingestion, ambient origin, or endogenous myelin synthesis.
- 6. Tissue repair/remyelination is not established by GO CC enrichment or network virtual knockout.
- 7. Unknown source/stage/sample metadata or failed exact cell-ID matching MUST be reported, not imputed.
- Input files:
- Analysis status:
- Created outputs:

## sensitivity/regulatory/137_decontx_cell_id_check.R

- STEP 137 | Cell-ID provenance and DecontX | v01
- Raw/corrected pairing QC:
- Atlas-to-DecontX crosswalk QC:
- CRITICAL LINKAGE RULES:
- - Only exact identifiers, explicit original cell IDs, or SAMPLE-QUALIFIED 16nt barcodes may form accepted links.
- - Any conflicting candidate, sample mismatch, or Control/PLX uncertainty invalidates the corresponding join.
- - Source accession is never guessed; raw/corrected sources are defined by existing project input file paths.
- - Barcode identity alone, approximate string matching, suffix stripping, matching by ordering, and cell-number proximity are prohibited.
- - If an exact original-ID pair cannot be reconstructed, report failure instead of manufacturing pairings.
- ASSAY CAUTION:
- - Verify whether Step03 recorded RNA/counts are genuinely DecontX-corrected. Comparing nominal RDS stages alone is insufficient.
- - Crosswalk reports count-level concordance with Step20 metadata and whether individual gene counts actually changed.
- - Step20 inherited decontX_contamination is descriptive and does not estimate the source of each myelin transcript.
- - No RNA count-only result establishes phagocytosed origin, genuine endogenous transcription, myelin synthesis or repair.
- Source-level warnings:
- Output: all tables in tables/; this report and sessionInfo in parent directory; runtime log in logs/.

## sensitivity/regulatory/138_decontx_myelin_sensitivity.R

- STEP 138 | Multilayer count/ID DecontX provenance | v01
- Raw/corrected pairing QC:
- Atlas-to-DecontX crosswalk QC:
- Count-layer inventory:
- Split-layer / independent UMI validation:
- CRITICAL LINKAGE RULES:
- - Only exact identifiers, explicit original cell IDs, or SAMPLE-QUALIFIED 16nt barcodes may form accepted links.
- - Any conflicting candidate, sample mismatch, or Control/PLX uncertainty invalidates the corresponding join.
- - Source accession is never guessed; raw/corrected sources are defined by existing project input file paths.
- - Barcode identity alone, approximate string matching, suffix stripping, matching by ordering, and cell-number proximity are prohibited.
- - If an exact original-ID pair cannot be reconstructed, report failure instead of manufacturing pairings.
- ASSAY CAUTION:
- - RNA counts.* layers are processed separately; any overlapping cell IDs across layers abort that input.
- - Gene absent in a given counts.* layer is unassayed (not set to zero).
- - Verified ID pairs failing Step20 UMI references or corrected>raw checks are excluded from expression comparison.
- - At least one observed count decrease in an input branch is required before reporting corrected-count sensitivity.
- - Crosswalk reports count-level concordance with Step20 metadata and whether individual gene counts actually changed.
- - Step20 inherited decontX_contamination is descriptive and does not estimate the source of each myelin transcript.
- - No RNA count-only result establishes phagocytosed origin, genuine endogenous transcription, myelin synthesis or repair.
- Source-level warnings:
- Output: all tables in tables/; this report and sessionInfo in parent directory; runtime log in logs/.

## sensitivity/regulatory/141_TF_correlation_sensitivity.R

- Step 141 v03 -- / Figure 6E technical-complexity sensitivity
- QUESTION:
- Do the pooled TF-expression versus functional-gene-expression correlations in Figure 6E
- persist after adjustment for per-cell RNA technical complexity?
- SCOPE:
- v03 performs technical-complexity adjustment only.
- No sample_id, source, stage, treatment, library/pool, specimen, or unit effects are modeled.
- No source-stratified, sample-level, or animal-level replication analysis is performed.
- INPUT:
- Original Step90 v05 10 TFs x 20 functional genes, Step84 per-cell state/QC metadata,
- Step20 RNA counts Seurat RDS, and the Step90 raw pooled Figure 6E correlation reference.
- No ZIP inputs. No SCENIC re-inference, no network reconstruction, no virtual KO, and no GSEA.
- METHOD:
- Expression is reproduced exactly as Step90: log1p(10000 * count / total RNA UMI).
- Original Figure 6E Spearman correlation is reproduced by correlating expression ranks.
- Technical-adjusted correlations are Pearson correlations between ranked TF and target-gene
- expression residuals after conditioning on ranked log1p(total RNA UMI) and ranked log1p
- number of detected genes; this is a partial Spearman correlation.
- Cell-based p-values and BH-adjusted values are approximate descriptive indicators only.
- INTERPRETATION:
- Figure 6E measures TF RNA expression versus functional-gene RNA expression.
- It does not measure regulon AUCell activity, motif-supported direct targeting, or virtual-KO direction.
- Attenuation after technical adjustment indicates that part of the pooled correlation structure can be
- explained by per-cell RNA detection complexity. Residual positive correlations remain associative and
- do not establish direct TF regulation or specimen-level reproducibility.
- OUTPUTS:
- 00 input audit; 01 panel; 02 raw reproduction; 03 pair-level raw vs technical-adjusted;
- 04 TF-level summary; 05 technical-complexity QC; 06 heatmap; 07 final QC; 08 report; 09 sessionInfo.

## sensitivity/regulatory/143a_network_reconstruction_sensitivity.R

- Step 143 v01: whole-transcriptome de novo alternative gene-selection sensitivity
- genome-wide gene-selection sensitivity ONLY (other seed/subsample/q/K completed by Step135 v02)
- Upstream work is read-only. No ZIP used. One new 3000-gene network is reconstructed.
- Algorithm: log1p(CPM/10k) variance; 20 mean-expression bins; median/MAD standardized dispersion; top-ranked genes across full RNA; ten KO TFs always retained.
- Paired comparison with B00 fixes cells, seed, nNet=10, nCells=500, q=0.90 and K=3.
- CAVEATS: mean-bin ranking and absence of functional20 forcing change jointly; no isolated estimate of forced-panel effect.
- CAVEATS: gene-network KO measures are computational associations, NOT gene expression up/down or proven causality.
- CAVEATS: gene-selecting from full transcriptome is an alternative reconstruction, NOT motif-pruning robustness.
- Output status: COMPLETE

## sensitivity/regulatory/143b_network_sensitivity_summary.R

- Step 143 v02: continuation of whole-transcriptome gene-selection sensitivity
- This script resumes post-processing from saved Step 143 v01 KO CSV.GZ files.
- No network reconstruction, virtual KO, expression matrix, or RDS loading was performed.
- FGSEA uses ranked scTenifoldKnk network Z; positive/negative NES is NOT RNA expression direction.
- The A01-vs-B00 comparison changes the whole-transcriptome gene-selection strategy jointly; it does not isolate the forced functional-panel effect.
- Network robustness is computational reproducibility, not experimental causality.

## sensitivity/regulatory/144_full_universe_enrichment.R

- Step 144 v01: full-universe lipid-program enrichment sensitivity analysis
- Each network strategy was tested in its own complete non-KO network gene universe.
- The 842-gene paired intersection was NOT used for the primary FGSEA in this script.
- FDR was BH-adjusted across all evaluable TF x pathways separately within B00 and A01.
- NES and leading-edge comparisons across strategies are descriptive because the ranked universes differ.
- scTenifoldKnk Z and FGSEA NES are network-distance ranking statistics, not expression direction.
- No network reconstruction, virtual KO rerun, expression matrix loading, or RDS reading was performed.

## sensitivity/single_cell/104a_source_harmony_sensitivity.R

- 104 | SOURCE/INTEGRATION SENSITIVITY
- READ BEFORE USING IN RESPONSE LETTER:
- (1) Original myeloid cell cohort and reference states are from Step 18.
- (1b) -label denominators include all 58,533 original retained cells.
- (1c) A small number of zero-count cells, if present, are EXCLUDED ONLY from PCA/Harmony refits.
- (1d) See 00_104_RNA_count_QC_by_source_v03.csv and problem-cell audit for exact counts.
- (2) Leave-out proportions on labels are DESCRIPTIVE reweighting only.
- (3) Branch PCA/Harmony clustering is rerun; label overlap remains non-independent.
- (3b) 05 includes a recalculated marker-rich CLUSTER proxy; it is not a cell identity call.
- (3c) marker co-detection is invariant to PCA/Harmony on identical cells.
- (3d) Reclustered marker-rich proxy thresholds are fixed from the full RNA
- cell set (mean + 0.5SD score; co-detection >= full-set mean).
- (4) Fixed-marker co-detection is an expression proxy, NOT lineage or function.
- (5) Harmony source mixing alone cannot prove biological signal preservation.
- (6) Cell counts cannot replace specimen-level replication or provide P values.
- (7) 14 dpi vs 28 dpi is fully confounded with dataset source.
- (8) GSE196928 Control corresponds to SHAM; do not call uninjured.
- (9) Only GSE196928 Vehicle samples are present; PLX5622 cannot be analyzed.
- (10) 4 h and 38 dpi are atlas-only (outside the main six-stage comparison).
- RECOMPUTATION STATUS:
- KEY OUTPUTS:
- 00 RNA count integrity audit by source/stage; zero-count problem cells
- 01 source/stage/treatment coverage and denominators
- 02 within-source labels and fixed-marker co-detection
- 03 leave-one-out -label proportions including absent stages
- 04 recomputed cluster/marker ranking (PCA/Harmony branches)
- 05 marker score/positivity by stage and dataset (same fixed cells)
- 06 pre/post Harmony adjusted Rand index for newly computed clusters
- 07 stage-matched neighborhood source mixing; requires RANN
- 08 analysis successes/failures (do not cite failed branches)
- 09 non-identifiable comparisons, suitable for revision summary
- 10 PDF: descriptive trends and UMAP where enabled
- NOTE: The stage-by-source contrasts are observational and partly confounded.
- RUN FAILURE = NO EVIDENCE. Examine the .log and 08 status table.

## sensitivity/single_cell/104b_marker_diagnostics.R

- 104 | MARKER DIAGNOSTICS v04
- STRICT MARKER CO-DETECTION DIAGNOSIS:
- HOW TO INTERPRET:
- 1. A strict same-cell multi-marker rule is sensitive to dropout, feature coverage and UMI depth.
- 2. Within-source percentile is calculated on the full retained cell distribution of each source,
- with no use of labels when calculating scores. However, the original label comparison
- is NOT an independent validation, since the upstream atlas and annotation were used as input.
- 3. Expression scores represent RNA programs only, not proven microglial lineage or lipid function.
- 4. Within-source top-decile status is source-relative and NOT a validated Clec7a+ cell identity.
- 5. The leave-one-out table uses original cell annotations or fixed RNA scores; no PCA/Harmony
- or clustering was rerun in v04. Interpret as descriptive composition sensitivity only.
- 6. v03 did not save per-cell PCA/Harmony/LODO cluster IDs. The v03 cluster proxy cannot be
- recalculated from its aggregate CSV files. For cluster-level reassessment, rerun clustering
- with saved per-cell memberships, or use explicitly stored branch cell IDs if later available.
- 7. Zero raw RNA counts for a gene may reflect feature incompatibility or assay construction,
- not biological absence. Investigate source-specific matrices if a core marker is always zero.
- 8. Biological independent specimens remain uncertain (Figshare) or limited (GSE196928 pooled).
- Cell counts are not independent replication, and no biological-replicate P values are generated.
- 9. 14 dpi occurs only in Figshare; 28 dpi only in GSE196928. Stage and source are not separable.
- 10. GSE196928 Control refers to sham surgery under Vehicle; PLX5622 cells are excluded.
- OUTPUT INDEX:
- 00 RNA count audit
- 01-04 Per-gene marker detection/expression by source/stage/original identity
- 05-07 Strict/pairwise/broad detection cascade and fixed scores
- 08 Sequencing-depth stratum detection (descriptive)
- 09 Source-stage scores and annotation fractions
- 10 Leave-one-out descriptive scores / stage coverage
- 11 Diagnoses of why strict co-detection was zero
- 12 Optional summary of v03 cluster-proxy caveat
- 13 Cell-level marker score audit (.csv.gz) for future branch cluster matching
- 14 Diagnostic PDF plots
- 15 Session info
- IMPORTANT: DO NOT claim v04 proves state stability across PCA/Harmony or leave-one-out re-clustering.

## sensitivity/single_cell/104c_state_stability.R

- 104 v08 | common-marker STATE STABILITY
- ANALYSIS DESIGN:
- - Refit HVG/PCA/neighbors/clusters per branch, plus Harmony for >=2 sources.
- - Harmony corrects source_dataset ONLY; not injury stage.
- - PCA-only and Harmony models start from the SAME positive-count cells.
- - Two predetermined exploratory proxy thresholds; labels NOT used to pick clusters.
- - Each gene rank-normalized within SOURCE, averaged; top 20% is SOURCE-relative.
- - Source-relative module score is not an absolute cross-dataset marker metric.
- - Cell-level module score invariant to PCA/Harmony; selected CLUSTERS can change.
- - Results are conditional on Step18 selected myeloid cells and upstream annotations.
- - All cell-level comparisons are descriptive, NOT animal-level inferential statistics.
- - Clec7a transcript cannot be independently assayed in two of three sources.
- - Original 14 dpi and 28 dpi are from different sources; stage confounding remains.
- - GSE196928 sham baseline is NOT uninjured; PLX5622 groups were excluded.
- BRANCH STATUS:
- MAIN RESULTS:
- 00 RNA count integrity by source
- 01 common marker feature availability and RNA detection by source
- 02 thresholds fixed before comparing cluster assignments
- 03 source-stage state plus expression-score trends
- 04 refitted cluster proxy audit with mean score and original-state overlap
- 05 source-stage reclustered proxy proportion, PCA/Harmony/LODO
- 06 within-branch PCA-versus-Harmony ARI and cell-level proxy agreement
- 07 success/failure per refitting branch
- 08 LODO-versus-all same-cell ARI / marker proxy overlap
- 09 non-identifiable comparisons
- 10 compressed per-cell cluster and proxy assignments
- 11 PDF with descriptive evidence and sensitivity plots
- IMPORTANT: A branch failure supplies NO evidence; examine 07 output.
- Proxy thresholds are operational exploratory definitions, NOT independently validated states.

## sensitivity/single_cell/140_doublet_rate_sensitivity.R

- fixed whole-dataset doublet-rate sensitivity
- Dataset: Figshare 17702045 ONLY (no GSE196928 PLX5622)
- Direct input: 03_DecontX_17702045_v01.rds; expected 19,956 QC-passing cells.
- Reason for fixed 7%: source contains no verified 10x capture/library ID;
- 7% was an operational expected-rate assumption, NOT a measured empirical rate.
- Treating all 19,956 cells as one capture for automatic dbr would be unjustified.
- Timepoint-derived sample_id labels were NOT treated as true capture batches.
- Method: same scDblFinder call with set.seed(1234) for each of 5%, 7%, 10%;
- only the dbr parameter changes; scDblFinder is actually rerun each time.
- Alternative scenarios:
- -label descriptive phase percentages (NOT specimen-level inference):
- Scope / limitations:
- 1. scDblFinder ran once per entire 17702045 dataset for each dbr.
- 2. annotation comes from historical 7% atlas; only original annotated
- cells that remain singlets are summarized. Newly rescued cells are NOT
- silently assigned identities. Conservative bounds are supplied separately.
- 3. This does NOT establish the response of a fully re-integrated/re-annotated
- atlas to alternative filtering; do not claim that without further analysis.
- 4. Lack of real capture/library labels and independent biological replicates
- prevents per-capture expected-rate calibration or specimen-level inference.
- 5. 28 dpi interpretation is across the three-source atlas, not a claim
- that Figshare 17702045 by itself contains a 28 dpi sample.
- 6. No original results, RDS, figures, or manuscript files were altered.
- Review QC flags in tables/09_140_FINAL_QC_v01.csv before writing a reply.

## sensitivity/spatial/114_uninjured_section_validation.R

- Step 114: GSE195783 candidate uninjured-section validation
- Targeted gene-set scoring from the eight already-extracted 10x matrices.
- No ZIP/TAR archive was used.
- No RDS/RData file was read or loaded.
- No global clustering, differential expression, or RCTD was rerun.
- Scores are descriptive and do not by themselves establish histological status.
- Gene-wise z-scores are standardized within each capture area; use them primarily
- for comparisons among the four coordinate-defined pieces in the same capture area.
- logCP10K summaries are also exported to retain normalized expression-scale context.
- SECTION IDENTITY LIMIT
- Step 68 section_uid_internal identifies coordinate-defined disconnected pieces.
- Without original animal/tissue-block records or H&E registration, these remain
- coordinate-defined section candidates, not verified original animal IDs.
- CANDIDATE SECTIONS
- GENE-SET COVERAGE
- INTERPRETATION
- A candidate is directionally supported when it has higher normal-like scores,
- lower injury-associated scores, and a lower fraction of injury-high spots
- than the other three coordinate-defined pieces in the same capture area.
- This is a within-capture-area descriptive comparison, not an inferential test.
- The Q75/Q80 high-score fractions are relative thresholds within each capture area.
- Moran's I is reported descriptively without permutation-based p-values.
- IMPORTANT BIOLOGICAL LIMIT
- Normal neuronal/myelin signals can persist in spared tissue within an injured
- spinal cord section. Therefore, high normal scores alone do not prove that an
- entire section is uninjured. Injury-myeloid and fibrotic/ECM signatures and
- their spatial distributions must also be considered.
- No candidate was automatically reclassified or removed.

## sensitivity/spatial/117_capture_area_heterogeneity.R

- DATA AND CLASSIFICATION
- - Raw GSE195783 labels: 8 GSM/capture areas, 4 coordinate pieces each (32 pieces).
- - Primary sensitivity cohort: 27 putative injured pieces (3/7/14/28dpi: 8/7/6/6).
- - Five inferred-uninjured pieces, using Step114-derived correction:
- - Classification is inferred from the original paper's totals + Step 114 phenotype.
- Original animal, OCT block, and physical section IDs remain unavailable.
- PRIMARY SUMMARIES
- - Each of 10 regions: proportion pooled by REGION SPOTS / QC-passed RCTD SPOTS.
- - Secondary: equal mean of 3 or 4 within-area piece proportions.
- - Clec7a primary: equal mean of piece-specific RCTD means; pooled-spot mean sensitivity.
- - Injury-associated region classes are the four expression-defined categories
- beginning Myeloid/immune, Fibrotic/ECM, Acute proliferative, IFN-responsive.
- - These classes are not histologically validated lesion-core labels.
- - Niche burden = sum of Clec7a predicted weights in region / all spots of piece.
- At capture level we report a spot-weighted pooled version and equal-piece mean.
- - JSD is a composition distance only, not a hypothesis-test p-value.
- - Analyses retain 8 separate capture areas rather than pooling all spots as replicates.
- STAGE-MATCHED CAPTURE-AREA COMPARISON (before vs inferred controls excluded)
- (JS is base-2 Jensen-Shannon distance: lower = more similar)
- REPRESENTATIVE SECTIONS
- - Original Figure4A representative IDs are read from Step69, not selected retrospectively.
- - Rank against all corrected-stage putative injured pieces, using an equal-capture-area
- mean regional composition as the stage center; these are descriptive ranks.
- LIMITATIONS / INTERPRETATION (MANDATORY IN REPLY)
- - Normal-like transcriptomes can also indicate spared injured tissue.
- - The five corrected labels are PUTATIVE, not confirmed original experimental IDs.
- - Per-stage capture-area count is only two, so do not claim formal replicated
- stage-level inference, capture-area confidence intervals or batch-adjusted significance.
- - Section identity is coordinate-defined, not histologically registered.
- - Spatial anatomy, exact distance from epicenter, OCT block and animal IDs remain unknown.
- - Capture-area and specimen effects are not independently identifiable from these metadata.
- - Step62/63 batch context predates inferred-control relabeling; not causal evidence.
- - Five inferred-control inclusion/exclusion comparison documents sensitivity;
- changes do not by themselves prove that control admixture caused the difference.
- - Preexisting Step71 annotated regions and Step65 RCTD weights share expression data.
- This is not independent validation of niche biology.
- - Figure4A representative sections are illustrations; whether representative should be
- stated with descriptive distances and full-section supplemental display.

## sensitivity/spatial/122_matched_control_spatial_sensitivity.R

- Step 122 | Matched controls and spatial dependence
- NO RCTD refitting. Original 37-gene score reproduced from Step64 log-expression.
- Primary null: median of within-portion Fisher z; descriptive comparison to matched 37-gene panels.
- Independent biological specimens are INFERRED from original sample design; animal IDs are unavailable.
- 8 capture areas are technical capture groups, not the biological replicate count.
- Spatial dependence: graph radius 1.55 x median within-portion nearest-neighbor distance;
- per-portion Moran's I (descriptive), neighborhood Spearman, and quadratic spatial rank residualization.
- Combined spatial+QC+total-myeloid partial correlation is exploratory; target is part of myeloid total.
- Sign-flip sensitivity tests are capture-area level (8 groups); they are not spot-level tests.
- Region labels are transcriptionally annotated, NOT histopathologically registered lesion zones.
- Spot-level tests and claims of intracellular lipid metabolism are NOT warranted.
- Matched negative-control comparison (effect units are Fisher z):
- Capture-area sign-flip sensitivity:
- CAVEATS:
- 1) Random matched control programs can share general injury-associated biology; the test is for relative specificity only.
- 2) Step120 feature membership was reconstructed, not directly archived from Step65.
- 3) Spatial detrending controls smooth quadratic gradients but not every form of spatial autocorrelation.
- 4) Using same-spot expression for RCTD and program scoring leaves other shared technical factors.
- 5) Step117 exclusion of five inferred uninjured portions is a sensitivity, not verified ground truth.
- 6) Empirical control P reflects control-panel sampling, not independent animal or causal evidence.
- 7) Full 37-gene program is NOT a held-out or feature-disjoint program.

## sensitivity/trajectory/129_microglia_trajectory_sensitivity.R

- 129 | high-confidence MG-only Monocle3 sensitivity | v01
- No ZIP, prior Script129, or separately processed expression matrix used.
- High-confidence means transcriptionally microglia-CONSISTENT, NOT demonstrated resident origin.
- Step106 scores are gene-Z standardized WITHIN source across all myeloid cells;
- fixed threshold is operational, not experimentally validated. Injury can reduce MG markers.
- Monocle3 primary: fresh microglia-only RNA-count PCA/UMAP, clusters and graph, seed1234.
- Primary root in each partition ONLY if >=10 Control+Homeostatic cells; otherwise UNROOTED.
- Alternatives: all-stage Homeostatic, earliest-available (diagnostic only), farthest endpoint (stress test).
- No source regression: source and injury stage are confounded; omission/within-source refits assess dependency.
- Cycle sensitivity: source-wise S/G2M program Z scores; remove annotated proliferating MG and
- top 10% scored cells per source. This is not formal cell-cycle phase inference or causal adjustment.
- Reported correlations compare common cells in matched partition(s); source-refit partitions are
- matched by maximal cell Jaccard overlap. Even a high Spearman is not evidence of lineage progression.
- Never combine numerical pseudotime across different partitions, root rules, or new graph fits.
- Do not treat 17702045 stage-derived sample IDs as biological replicates; GSE196928 pooled libraries
- and sham controls are not independent uninjured animals. 14/28 dpi source confounding persists.
- No sample-level hypothesis tests or lineage tracing; conclusions remain descriptive.

## spatial/65c_RCTD_section_summary.R

- GSE195783 Step 65 v06
- Section-aware finalize-only update of normalized RCTD results
- STATUS
- v06 DOES NOT rerun RCTD.
- v06 DOES NOT renormalize RCTD weights.
- v06 joins Step65 v04 spot weights to the Step68
- coordinate-defined tissue-section map.
- INPUT RCTD
- INPUT SECTION MAP
- FORMAL JOIN
- GSM + array_x + array_y
- SECTION STRUCTURE
- PRIMARY TARGET
- SECTION-AWARE DESCRIPTIVE TEMPORAL SUMMARY
- INTERPRETATION
- Section-level summaries are descriptive.
- The 32 coordinate-defined tissue sections are not automatically
- interpreted as 32 independent animals because mouse identity is
- not established by the Step68 coordinate map.
- No section-level or per-spot hypothesis tests are performed in Step65 v06.
- Formal spatial inference will be handled in Step66 v02 with
- section-aware analysis plus capture-area-blocked sensitivity checks.
- FINAL QC

## spatial/66_spatial_autocorrelation.R

- GSE195783 Step 66
- Section-aware corrected Clec7a+ spatial autocorrelation and neighborhood analysis
- CORRECTED SPATIAL UNIT
- Each section is analyzed independently for graph construction and Moran's I.
- Section identity does NOT establish animal identity.
- SPATIAL GRAPH
- Graphs are rebuilt separately within each tissue section.
- No edges can connect different tissue sections.
- SECTION-LEVEL MORAN'S I
- DIRECT REGION ENRICHMENT - SECTION-LEVEL EFFECT SIZES
- Effect size = mean Clec7a normalized RCTD weight in region /
- mean weight in all other spots of the SAME tissue section.
- Section-level pooled P values below are EXPLORATORY because animal
- independence among the 32 coordinate-defined sections is not established.
- DIRECT REGION ENRICHMENT - CAPTURE-AREA-BLOCKED SENSITIVITY
- The median section-level log2 effect within each capture area is used
- as one blocked effect before one-sided Wilcoxon testing.
- NEIGHBORHOOD ENRICHMENT - SECTION-LEVEL EFFECT SIZES
- Neighbor means are calculated only from immediate neighbors within
- the SAME tissue section.
- Effect size = mean neighbor signal in region / mean neighbor signal
- outside that region within the SAME tissue section.
- Section-level pooled P values below are EXPLORATORY.
- NEIGHBORHOOD ENRICHMENT - CAPTURE-AREA-BLOCKED SENSITIVITY
- The median section-level log2 neighborhood effect within each capture
- area is used as one blocked effect before one-sided Wilcoxon testing.
- STATISTICAL INTERPRETATION
- Moran's I permutation P values test within-section spatial autocorrelation.
- Section-level direct/neighborhood effects describe spatial consistency
- at the coordinate-defined tissue-section scale.
- Because mouse identity is not established, section-level pooled Wilcoxon
- P values are exploratory and must not be presented as independent-animal
- biological-replicate inference.
- Capture-area-blocked tests are the conservative sensitivity analysis.
- No timepoint-specific hypothesis tests are performed.
- No per-spot differential P values are used.
- SPATIAL LIMITATION
- Coordinates are reconstructed Visium array coordinates.
- No H&E registration or H&E-defined lesion-core/border claim is made.
- FINAL QC

## spatial/67_section_identity.R

- GSE195783 Step 67 spatial section identity audit
- PURPOSE
- Audit disconnected tissue pieces within each Visium capture area.
- This audit is formally recorded as Step 67.
- METHOD
- PRIMARY RESULT
- THRESHOLD-STABILITY RESULT
- INTERPRETATION
- Major disconnected components are treated as candidate tissue sections.
- Tiny isolated components are reported separately and are not automatically
- interpreted as independent tissue sections.
- IMPORTANT LIMITS
- This is a coordinate-topology audit only.
- It does not establish animal identity by itself.
- It does not alter the Step64-66 analyses.
- It does not treat spots as biological replicates.
- No H&E-defined lesion center/border is inferred.
- FINAL QC

## spatial/68_section_mapping.R

- GSE195783 finalized coordinate-defined tissue-section map
- STATUS
- Step 68 is a lightweight finalize-only mapping operation.
- It is formally recorded as Step 68.
- Step60-66 were NOT rerun.
- PRIMARY INPUT
- FINALIZED STRUCTURE
- INTERNAL SECTION ASSIGNMENT
- Four major components in each capture area were assigned deterministic
- internal positions from their array-coordinate centroids:
- upper_left
- upper_right
- lower_left
- lower_right
- These IDs are computational grouping labels only.
- They are NOT printed on manuscript figures.
- THRESHOLD SENSITIVITY
- JOIN STRATEGY FOR STEP65/66 CORRECTION
- Future joins should use:
- GSM + array_x + array_y
- rather than barcode alone.
- IMPORTANT INTERPRETATION LIMIT
- A disconnected major component is treated as a coordinate-defined
- candidate tissue section.
- This map does NOT establish mouse/animal identity.
- Therefore section-level results must not automatically be described
- as independent animal-level biological replicates unless external
- metadata establish that correspondence.
- Tiny disconnected components remain unassigned.
- No H&E/pixel registration is claimed.
- NEXT USE
- This lightweight map can be joined to the Step65 v04 normalized
- RCTD spot-weight table without rerunning RCTD.
- FINAL QC

## spatial/73a_Clec7a_lipid_spatial_coupling.R

- GSE195783 Step 73
- Clec7a - lipid-handling spatial coupling
- INPUTS
- LIPID-HANDLING PROGRAM
- Prespecified from the same curated modules used in scRNA AUCell:
- FA uptake / transport
- cholesterol efflux / transport
- lipoprotein processing
- lysosomal lipid catabolism
- Core score = mean within-sample gene z-score across the union of
- detected program genes.
- Clec7a is explicitly excluded from the program.
- FORMAL SPOT JOIN
- Step64 lipid-program scores -> Step65 v06 RCTD weights:
- GSM + array_x + array_y.
- Barcode strings are not used as the formal join key.
- SPATIAL UNIT
- Each section is analyzed independently.
- Section identity does not establish animal identity.
- COUPLING METRICS
- Direct same-spot coupling:
- Spearman rho(Clec7a normalized RCTD, lipid-handling program score).
- Local spatial coupling:
- Spearman rho(Clec7a normalized RCTD at a spot, mean lipid-handling
- program score among immediate neighbors within the same section).
- TIMEPOINT DESCRIPTIVE SUMMARY
- No timepoint-specific hypothesis tests were performed.
- CAPTURE-AREA-BLOCKED FORMAL SENSITIVITY
- Section correlations are Fisher-z transformed, then median-aggregated
- within each capture area before one-sided Wilcoxon testing across
- capture areas.
- COMPONENT MODULES - DESCRIPTIVE ONLY
- INTERPRETATION
- Positive direct coupling supports spatial coincidence between
- Clec7a+ lipid-handling Microglia abundance and a lipid-handling
- transcriptional program.
- Positive local-neighbor coupling additionally supports organization
- of that metabolic program in the immediate spatial environment.
- INTERPRETATION LIMIT
- This analysis supports spatial coupling, not causal regulation.
- Visium spots are multicellular mixtures.
- Coordinate-defined sections are not assumed to be independent animals.
- No H&E-defined lesion core/border is inferred.
- FINAL QC

## spatial/73b_local_bivariate_coupling.R

- GSE195783 Step73 v05
- Local bivariate Clec7a-lipid spatial coupling
- METHOD
- The Step73 v03 core lipid-handling score is reused without rescoring.
- Within each coordinate-defined major section:
- z_Clec7a = standardized Clec7a normalized RCTD weight.
- z_lipid  = standardized core lipid-handling program score.
- lag(z_lipid) = mean z_lipid among immediate spatial neighbors.
- local bivariate I = z_Clec7a * lag(z_lipid).
- QUADRANT INTERPRETATION
- High-High: high Clec7a at a spot + high lipid program in its local neighborhood.
- High-Low : high Clec7a + low neighboring lipid program.
- Low-High : low Clec7a + high neighboring lipid program.
- Low-Low  : low Clec7a + low neighboring lipid program.
- IMPORTANT
- This is a descriptive local spatial map.
- No spot-level permutation P values are calculated.
- Formal inferential support remains the capture-area-blocked Step73 v03 analysis.
- Coordinate-defined sections are not assumed to be independent animals.
- TIMEPOINT DESCRIPTIVE SUMMARY
- REPRESENTATIVE SECTIONS
- Same sections as Figure4A/B:
- MAIN-FIGURE CANDIDATE
- Grey = tissue scaffold.
- Red intensity = High-High local Clec7a-lipid coupling strength.
- DIAGNOSTIC / SUPPLEMENTARY MAPS
- INTERPRETATION LIMIT
- High-High spots are descriptive spatial coupling hotspots,
- not independently replicated biological units.
- This analysis supports spatial co-organization, not causal regulation.
- Visium spots remain multicellular mixtures.
- No H&E-defined lesion core/border is inferred.
- FINAL QC

## spatial/74_same_spot_coupling.R

- GSE195783 Step74
- Same-spot Clec7a-lipid spatial coupling
- INPUTS
- METHOD
- Step73 v03 core lipid-handling score is reused without rescoring.
- Within each coordinate-defined major tissue section:
- z_Clec7a = standardized Clec7a normalized RCTD weight.
- z_lipid  = standardized core lipid-handling program score.
- Same-spot state is assigned from the signs of z_Clec7a and z_lipid.
- High–High = z_Clec7a > 0 AND z_lipid > 0 in the SAME Visium spot.
- HH strength = z_Clec7a * z_lipid for HH spots.
- BIOLOGICAL MEANING
- High–High spots are Clec7a-enriched locations that simultaneously
- show above-section-average lipid-handling transcriptional activity.
- This supports spatial coincidence / co-organization.
- IMPORTANT LIMIT
- Visium spots are multicellular mixtures.
- This analysis does NOT prove that every lipid-program transcript
- originates from the Clec7a+ cells themselves.
- No spot-level P values are calculated.
- No H&E-defined lesion compartment is inferred.
- TIMEPOINT DESCRIPTIVE SUMMARY
- No timepoint-specific hypothesis tests were performed.
- REPRESENTATIVE SECTIONS
- Same sections as Figure4A/B:
- MAIN-FIGURE CANDIDATE
- Grey = all tissue spots.
- Red intensity = same-spot High–High Clec7a–lipid coupling strength.
- DIAGNOSTIC / SUPPLEMENTARY
- Full High–High / High–Low / Low–High / Low–Low map.
- FINAL QC

## external_cohort/125_major_lineage_annotation.R

- STEP125 v01 | GSE281205 | REVIEWED MAJOR-LINEAGE ANNOTATION
- NOTE: Myeloid INCLUDE does not imply Clec7a+ microglia identity.
- NOTE: Labels were applied from the manually edited Step125 decision file.
- NOTE: No frozen 1320-gene signature and no pseudobulk/DE in this script.
- NOTE: Data are DecontX-corrected; original integer counts stay at Step119.
- WARNING: GEO Young/Old treatment inconsistency for 28dpi remains.
- WARNING: Harmony may attenuate 14 vs 28 dpi differences.

## regulatory/86_GRNBoost2.R

- Step 86 v03: pySCENIC GRNBoost2
- METHOD
- pySCENIC 0.12.1
- GRNBoost2
- OUTPUT

## spatial/61a_spatial_qc_normalization.R

- GSE195783 Step 61
- Spatial spot QC + conservative filtering + normalization
- QC METHOD
- No upper nCount/nFeature cutoff was used.
- PER-SAMPLE THRESHOLDS
- RETENTION
- NORMALIZATION
- Samples normalized independently.
- No cross-sample integration performed.
- SPATIAL INTERPRETATION LIMIT
- Coordinates remain reconstructed Visium V1 ARRAY coordinates.
- They are NOT H&E-registered pixel coordinates.
- OUTPUT CANDIDATES
- NEXT STEP
- FINAL QC

## spatial/62_spatial_clustering.R

- GSE195783 Step 62
- Joint unintegrated embedding + batch/sample assessment + clustering
- ANALYSIS DESIGN
- Joint assay: Spatial
- Step61 LogNormalize data reused.
- No new NormalizeData call was performed.
- No cross-sample integration/batch correction was performed.
- Rationale: sample is nested within biological timepoint, so correction
- should not be applied before checking whether sample-specific structure
- materially exceeds true timepoint structure.
- JOINT DIMENSION REDUCTION
- PRELIMINARY CLUSTERING
- PC SAMPLE-STRUCTURE AUDIT
- R2_timepoint = variance in each PC associated with biological timepoint.
- R2_sample_label = variance associated with the 8 individual samples.
- Extra sample R2 = max(R2_sample_label - R2_timepoint, 0).
- Because sample is nested in timepoint, this extra term is used only as
- a descriptive screen for replicate/sample-specific structure.
- Top PCs by variance-weighted extra sample structure:
- SPATIAL INTERPRETATION LIMIT
- Coordinates remain reconstructed Visium V1 ARRAY coordinates.
- They are NOT H&E-registered pixel coordinates.
- NEXT DECISION
- Review files 07, 08, 11 and 12 together.
- If replicate/sample-specific separation is weak relative to timepoint
- structure, retain the unintegrated representation.
- If strong sample-specific separation is present, perform a separate,
- conservative batch-correction sensitivity analysis rather than blindly
- replacing this unintegrated result.
- FINAL QC
- FINAL OBJECT

## spatial/64a_spatial_region_markers.R

- GSE195783 Step 64
- Spatial cluster marker characterization + transcriptional program annotation
- PRIMARY REPRESENTATION
- Step62 UNINTEGRATED representation retained.
- Step63 Harmony results are sensitivity analysis only.
- MARKER ANALYSIS
- SAMPLE-AWARE CONSISTENCY
- robust_directional_support requires:
- >=2 evaluable samples
- >=75% samples with positive within-sample expression difference
- positive median within-sample difference
- ROBUST MARKER COUNTS
- PROVISIONAL TRANSCRIPTIONAL PROGRAMS
- IMPORTANT INTERPRETATION LIMIT
- These are transcriptionally defined programs/regions.
- They are NOT histology-defined lesion-core / scar compartments.
- No H&E-registered pixel coordinates are available in the processed GEO files.
- NEXT STEP
- Review the marker table, program heatmap, and spatial maps.
- After freezing transcriptional-region interpretation, proceed to
- scRNA-to-spatial mapping/deconvolution focused on the Clec7a+
- lipid-handling Microglia state.
- FINAL QC

## spatial/64b_spatial_region_annotation.R

- GSE195783 Step 64 v03
- Curated spatial-cluster annotation refinement
- IMPORTANT
- FindAllMarkers was NOT rerun.
- Sample-aware marker consistency was NOT rerun.
- Program scoring was NOT rerun.
- PCA / UMAP / clustering were NOT rerun.
- Step62 unintegrated representation remains primary.
- CURATED LABELS
- KEY REFINEMENTS
- Cluster 2 -> Myeloid/immune–lipid-handling enriched.
- Cluster 5 -> Acute proliferative injury-associated myeloid.
- Cluster 6 -> Peripheral-myelin / Schwann-like region.
- Cluster 8 -> Fibrotic/ECM–myeloid mixed injury niche.
- Cluster 9 -> IFN-responsive inflammatory/astroglial region.
- INTERPRETATION LIMIT
- These are transcriptionally defined spatial regions/programs.
- They are NOT H&E-defined anatomical lesion compartments.
- Visium spots remain multicellular mixtures.
- No pure scRNA-defined cell subtype is assigned from Step64 alone.
- DOWNSTREAM USE
- This object can be used as spatial-region context for Step65
- scRNA-to-spatial mapping/deconvolution.
- FINAL QC

## spatial/65a_RCTD_mapping.R

- GSE195783 Step 65
- Hierarchical RCTD scRNA-to-spatial mapping
- REFERENCE DESIGN
- 10 broad non-myeloid classes from the global atlas.
- 9 finalized myeloid states replace the 3 broad myeloid classes.
- spacexr-safe internal labels use only letters/numbers/underscores;
- all outputs are mapped back to manuscript display labels.
- RCTD DESIGN
- Run independently for each of the 8 spatial samples.
- Checkpointing enabled per sample.
- SAMPLE QC
- PRIMARY TARGET
- TARGET TEMPORAL SUMMARY
- SANITY CHECK
- Correlation is descriptive only and not independent validation.
- IMPORTANT INTERPRETATION
- RCTD weights are reference-informed estimated cell-type proportions.
- Visium spots remain multicellular mixtures.
- Region/time summaries are descriptive at this step.
- No per-spot P values were used.
- Coordinates are Visium array coordinates, not H&E-registered pixels.
- NEXT DECISION
- Review target temporal pattern, region enrichment and spatial maps.
- If technically stable, freeze Step65 and proceed to spatial enrichment/
- neighborhood statistics using biological sample as the replication unit.
- FINAL QC
- FINAL ANNOTATED OBJECT

## spatial/65b_RCTD_weight_normalization.R

- GSE195783 Step 65 FINAL
- Hierarchical RCTD scRNA-to-spatial mapping
- RCTD computation version: v02
- IMPORTANT
- All 8 expensive RCTD computations were completed in v02.
- v04 loaded those checkpoints and did NOT rerun RCTD.
- WEIGHT NORMALIZATION
- For every Visium spot:
- normalized weight = raw weight / sum(raw weights across all 19 types).
- All downstream summaries and plots in v04 use normalized weights.
- RCTD DESIGN
- Mode: full.
- Reference: 10 broad non-myeloid classes + 9 finalized myeloid states.
- Eight spatial samples analyzed independently.
- SAMPLE QC
- PRIMARY TARGET
- NORMALIZED TARGET TEMPORAL SUMMARY
- INTERPRETATION
- Normalized RCTD weights are reference-informed estimated proportions.
- Visium spots remain multicellular mixtures.
- Region/time summaries remain descriptive at Step65.
- No per-spot P values were used.
- Coordinates are reconstructed Visium array coordinates, not H&E pixels.
- FINAL QC
- FINAL ANNOTATED OBJECT

