# linearSNPfoldR: estimate the effect of single nucleotide polymorphisms in ribonucleic acid secondary structure.
This R package implements the algorithm described by Halvorsen et al[1] to estimate the effect of single nucleotide polymorphisms (SNPs) in the secondary structure of ribonucleic acid (RNA). The algorithm, originally called SNPfold, uses a cubic time RNA folding routine implemented in ViennaRNA's RNAfold[2]. In addition to leveraging the algorithm in R, in this package the folding routine is powered by LinearPartition[3]. This is a recent algorithm that implements heuristics to fold RNAs in approximately linear time.

## References
- `1` Halvorsen M, Martin JS, Broadaway S, Laederach A (2010) Disease-Associated Mutations That Alter the RNA Structural Ensemble. PLOS Genetics 6(8): e1001074. https://doi.org/10.1371/journal.pgen.1001074
- `2` Lorenz, R., Bernhart, S.H., Höner zu Siederdissen, C. et al. ViennaRNA Package 2.0. Algorithms Mol Biol 6, 26 (2011). https://doi.org/10.1186/1748-7188-6-26
- `3` He Zhang, Liang Zhang, David H Mathews, Liang Huang, LinearPartition: linear-time approximation of RNA folding partition function and base-pairing probabilities, Bioinformatics, Volume 36, Issue Supplement_1, July 2020, Pages i258–i267, https://doi.org/10.1093/bioinformatics/btaa460
