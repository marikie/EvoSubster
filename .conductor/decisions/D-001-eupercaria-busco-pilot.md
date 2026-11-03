# D-001 - Eupercaria BUSCO pilot Reference-override judgment

Date: 2026-08-07
Tier: 0 (conductor re-derivation plus independent adversarial verification)
Conductor: GPT-5 Codex
Inputs judged: `.conductor/handoffs/T003-handoff.md`, `.conductor/handoffs/T004-handoff.md`, `results/eupercaria_20260807/busco_experiment/external_busco_qc.tsv`, six specific BUSCO JSON summaries, and `results/eupercaria_20260807/busco_experiment/stage0_rerank/baseline_vs_busco.tsv`

## Question

Under the current EvoSubster Stage 0 ranking rule, can comparable external BUSCO genome-mode evidence replace a Reference assembly with a distinct non-reference assembly for species represented by the supplied Eupercaria tree?

## Verdict

| Species | Verdict | Reason |
|---|---|---|
| `Plectropomus leopardus` | VALID clear override | The distinct non-reference `GCA_026936395.1` has 99.9% Complete, 0 Fragmented, and 5 Missing of 3,640 markers, versus 96.2%, 77 Fragmented, and 61 Missing for Reference `GCF_008729295.1`. Duplicated also improves from 46 to 11. |
| `Larimichthys crocea` | CONDITIONAL mechanical override only | Current code selects `GCA_059630365.1` at 99.8% Complete over Reference `GCF_000972845.2` at 99.7%, but Duplicated worsens from 10 to 215 and Missing from 3 to 6. Do not promote this as clearly better. |
| `Takifugu rubripes` | INVALID as non-reference override evidence | The apparent winner change is `GCF_901000725.3` to its paired `GCA_901000725.3`, with identical BUSCO scores. It is the same Reference assembly under an accession tie-break. |

## Outcome level + paired non-claim

- Ledger transition: draft -> reviewed -> accepted for the bounded Plectropomus claim; Larimichthys downgraded; Takifugu rejected as override evidence.
- O-level: O3. This is a reproducible, evidence-bounded three-species pilot using six real genome FASTAs, not an exhaustive or externally replicated conclusion.
- Paired non-claim: The result does not establish that the non-reference is globally superior for every downstream use, that Larimichthys is improved, or that any of the other 142 retained Eupercaria species would be overridden.
- Independent judgment panel: The independent JSON re-derivation/adversarial check agreed on Plectropomus, exposed the Larimichthys duplication confound, and rejected the Takifugu paired-accession false positive.

## Load-bearing reasoning

The pilot scope was fixed before any BUSCO score was observed: the previously discussed Takifugu example plus the two largest non-reference/reference contig-N50 ratios in the saved top-three shortlist. All six measured inputs used BUSCO 6.1.0 genome mode, Miniprot 0.18-r281, HMMER 3.4, three CPUs, and the same `actinopterygii_odb10` dataset dated 2024-01-08 with 3,640 markers. Six raw JSON summaries independently reproduce all nine exact-accession TSV rows, including declared paired-accession propagation. The current R ranking and an independent rerank agree.

Plectropomus is separated by 133 additional Complete BUSCOs, 77 fewer Fragmented, 56 fewer Missing, and 35 fewer Duplicated markers, so the conclusion does not depend on one-decimal rounding or a Complete-only tie. Larimichthys differs by only four Complete markers while gaining 205 Duplicated markers; because the current external-QC schema and ranking omit Duplicated, its mechanical selection is not sufficient evidence of superior assembly completeness. Takifugu never selects a biologically distinct non-reference assembly.

## Conditions attached

- Report the lineage as `actinopterygii_odb10`, not odb12.2; the switch was recorded before any candidate result because the 1.48 GB odb12.2 download was impractical locally.
- Preserve the Miniprot internal-stop-codon warning (6.0-9.2% of Complete matches) as a method caveat.
- Treat paired GCA/GCF accessions as one assembly when interpreting winner changes.
- Before making BUSCO quality-first selection a default policy, add Duplicated BUSCO to the external schema/ranking or define a guard against duplication-confounded marginal improvements.

## Applied

The external TSV is `results/eupercaria_20260807/busco_experiment/external_busco_qc.tsv`; the reranking comparison is `results/eupercaria_20260807/busco_experiment/stage0_rerank/baseline_vs_busco.tsv`; the adversarial comparison is `results/eupercaria_20260807/busco_experiment/adversarial_rerank_comparison.tsv`. Six FASTA SHA-256 values, six raw JSON summaries, and all commands/logs remain under the ignored experiment directory. Current Stage 0 gate, quality-selection, and trio-helper R tests passed. No application source was changed.

## Residual / follow-up

This pilot sampled 3 of 145 retained species and used odb10. A tree-wide result would require the other species' complete shortlist FASTAs, substantially more storage/compute, and the same duplication-aware interpretation. The immediate design follow-up is to decide whether BUSCO Duplicated and a minimum-improvement threshold should be added before expanding the experiment.

delegation_effect: difficulty=hard; agents=4 delegated workers; useful=parallel inventory/runtime audit, controlled execution, independent JSON re-derivation, and detection of one false and one confounded override; misses=no full-tree or odb12.2 replication; next_mix=one executor plus one adversarial verifier is sufficient for another bounded lineage pilot.
