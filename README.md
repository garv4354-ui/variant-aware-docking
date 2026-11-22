# Variant-aware drug screening: Gene → Mutation → Drug (EGFR, lung cancer)

**Author:** Garv Gupta, B.Tech Biotechnology
**Status:** in progress (docking and re-docking done; decoy validation running)

## Question
Does a mutation in a cancer-driver gene change how well approved drugs bind its protein, and can a docking pipeline reproduce known resistance mutations?

## Pipeline
1. **RNA-seq (GSE81089, 19 normal vs 199 tumor lung samples):** DESeq2 differential expression, GO enrichment, target selection.
2. **Mutations:** EGFR L858R, T790M, C797S, G719S and the L858R/T790M double mutant (table in `02_variants/`).
3. **Structures:** wild-type from PDB 1M17 (erlotinib complex), all mutants built from the same scaffold with PDBFixer; AlphaFold model P00533 kept for comparison.
4. **Docking:** AutoDock Vina, 10 compounds (7 EGFR inhibitors + 3 unrelated controls) against 6 EGFR variants.
5. **Validation:** re-docking of erlotinib, DUD-E decoy test (ROC-AUC).

## Results so far
- **RNA-seq:** 5,540 genes up and 2,303 down in tumor (padj < 0.05, |log2FC| > 1). Control genes behaved as expected (SPP1 up; AGER, SFTPC down). EGFR is modestly but significantly up (log2FC 0.64, padj 0.0068). It was chosen as target because it is a clinically validated kinase whose relevance is driven by mutations, not expression.
- **Docking:** EGFR inhibitors score about -8 to -9 kcal/mol, unrelated drugs about -5 to -7. Mutant-minus-wild-type differences are mostly under 0.5 kcal/mol, within Vina's error. The setup does not reproduce T790M resistance to gefitinib/erlotinib.
- **Re-docking (erlotinib, 1M17):** the top-ranked pose is 4.16 Å from the crystal pose; the crystal-like pose (1.31 Å) is ranked 2nd, 0.02 kcal/mol lower. A tighter box with higher exhaustiveness did not improve this (best 3.5 Å).
- **Decoy test (DUD-E EGFR subset):** TODO (ROC-AUC and 95% CI, to be filled in when the run finishes).

## Limitations
- Vina's score surface is flat for erlotinib: pose ranking is unreliable and differences of a few tenths of kcal/mol are not meaningful.
- Afatinib, dacomitinib and osimertinib bind covalently to Cys797; Vina cannot model this, so C797S results for these drugs are not interpretable.
- Mutants are built by side-chain replacement without relaxation, so shifts in the pocket may be artefacts.
- T790M resistance to first-generation drugs partly reflects tighter ATP binding, which this setup does not measure.
- The RMSD used is an approximate, element-matched value that ignores bonds.
- The DUD-E test uses a small random subset (20 actives, 80 decoys).

## Reproduce
Scripts are in numbered folders (`01_rnaseq` to `05_validation`). Environments: conda (`rnaseq`, `docking`).

## Data sources
GEO GSE81089; RCSB PDB 1M17; AlphaFold DB P00533; PubChem (ligands); DUD-E (EGFR actives and decoys).
