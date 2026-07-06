# Canonical methods & definitions (mirror of config/definitions.yaml)

This is the single human-readable statement of the frozen conventions. If it and
`definitions.yaml` ever disagree, `definitions.yaml` wins (it drives the code).

- **S′** = asinh( (E_max / EC50) × C_ref ), C_ref = 1 µM. E_max = upper − lower 4PL asymptote.
  ⚠ E_max convention (percent vs fraction) TBD — see concerns.csv C1.
- **pS′** = mean S′ over a cohort; **ΔpS′** = pS′_WT − pS′_mutant.
- **Synthetic-lethal rule:** WT pS′ > 0, mutant pS′ > 0, ΔpS′ ≤ −2 (≤ −4 = high stringency).
- **Genotype:** DepMap damaging-allele-frequency sum > 0.95 → mutant (value 2); 0 = WT; 1 = excluded.
  Single-gene unit of analysis; genes PTEN/CDKN2A/RB1/TP53.
- **Data:** PRISM 19Q4 SSDRC (screen HTS002 + MTS010 overlay, passed_str_profiling); DepMap 24Q2
  damaging mutations; DEMETER2 RNAi (validation).
- **Statistics:** bootstrap 95% CI on ΔpS′ (1000 iters). BH-FDR computed but not used to gate hits.
