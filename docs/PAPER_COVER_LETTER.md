# Cover Letter — Nature Genetics Submission

**Manuscript title:** Relational biological structure improves fine-mapping of causal GWAS variants under weak signal

**Corresponding author:** Jian-Feng Mao (Umeå University; jianfeng.mao@umu.se)

**Submission type:** Article (or Technical Report — per editor discretion)

---

Dear Editor,

We submit for consideration at *Nature Genetics* a methodology paper describing **two new fine-mapping algorithms — Hierarchical Belief Propagation (HBP) and L1 Bayesian fine-mapping — that carry multi-omics biological structure through the inference as a relational factor graph rather than as a flat per-variant annotation vector.** The contribution type is a new algorithmic class with a theoretical guarantee (Banach-contraction convergence for HBP, causal-variant ranking theorem for L1) and a decisive empirical result in a scientifically important regime: head-to-head against SuSiE at weak signal with tissue-specific eQTL priors, L1 wins 27–2 on 79 independent replicates (sign-test p < 10⁻⁵). HBP additionally matches the accuracy of all six methods benchmarked by Wu *et al.* 2026 — SuSiE, FINEMAP, SuSiE-inf, FINEMAP-inf, Polyfun-proxy and SBayesRC — at 6–60× the speed. All methods are PIP-calibrated with 0% null false-positive rate across 100 null replicates, and a sumstats-only entry path consumes real Pan-UK Biobank summary statistics across four ancestries with no algorithmic modification.

The relevant literature context is tight. Cui *et al.* (*Nat Genet* 2024) showed that adding an infinitesimal random-effect background to SuSiE and FINEMAP improves fine-mapping under non-sparse architectures; Zheng *et al.* (*Nat Genet* 2024) released SBayesRC as a genome-wide Bayesian mixture with HapMap-3 LD reference; Wu *et al.* (*Nat Genet* 2026) used SBayesRC as the substrate for genome-wide fine-mapping at biobank scale. All of these methods integrate annotations as flat per-variant priors. The methodological gap our paper addresses is whether the *relational* structure connecting variants to tissue-specific eQTLs, pathways and protein–protein interactions — information that flat priors discard by construction — can be carried through the inference. We show that it can, and that it changes the accuracy–speed frontier in a regime-dependent way: graph and flat priors converge on accuracy when the prior information is identical (Polyfun-proxy tie), but when the graph captures tissue-specific eQTL structure that a scalar prior averages out, the graph wins decisively on rank-#1 rate, mean rank and calibrated PIP. As multi-omics resources continue to densify with tissue-resolved annotations, this regime widens, not narrows.

The paper bounds its claims honestly. L1's 27–2 headline is specific to a weak-signal regime with tissue-specific eQTL causal enrichment; in a 100-replicate 1KG replication without that enrichment, L1 and SuSiE are at parity on rank and L1 retains the speed advantage only. We state this explicitly (main text §2.3). The Pan-UKB real-data demonstration uses 1KG-matched LD as a stand-in for Pan-UKB in-sample LD; this is a software-plumbing limitation (hadoop-aws Spark classpath) we are transparent about. SBayesRC is not compared on matched simulations because it targets a different problem (genome-wide, N ≳ 10⁵); we position it as complementary to HBP/L1, not competitive. GraphGWAS is a broader platform of which fine-mapping is the first method class rigorously benchmarked; Supplementary Note S3 provides an honest benchmark-status table for the other method classes (epistasis, heritability, PRS, MR, multivariate, GNN) that are implemented in the released code but either previewed for a companion manuscript in preparation (epistasis) or awaiting dedicated benchmarks.

We believe this work fits *Nature Genetics* because it reframes an active methodological debate (flat-prior versus relational-prior fine-mapping) and delivers quantitative evidence that the reframing matters in the regime most relevant to post-GWAS discovery. The Pan-UK Biobank cross-ancestry demonstration directly engages the journal's emphasis on ancestry-inclusive genetics, and the software-paradigm framing (graph as a first-class structure for statistical genetics, not a post-hoc enrichment lookup) is one we believe will be of interest across fine-mapping, polygenic score, and Mendelian-randomisation communities.

We have no prior submissions of this work to other journals. A preprint is being deposited on bioRxiv in parallel with submission. None of the authors have competing interests relevant to this manuscript.

Sincerely,

Jian-Feng Mao
on behalf of the GraphGWAS team (Estaji, Zhao, Chen, Nie, Mao)
Umeå Plant Science Centre, Umeå University
jianfeng.mao@umu.se
