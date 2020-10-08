# haplotype simulations  
**tool to perform the haplotype simulation study described in "An exact, unifying framework for rare variant association testing in family-based designs, including higher criticism approaches, SKATs, multivariate and burden tests" (Hecker et al. 2020)**

- haps= haplotype file containing the haplotypes for the simulation (nvar=1000, nhaps=1006). Haplotypes were extracted from the EUR set of the 1000 Genomes phase3 dataset.
- assumes file "haps" is in this same folder. Performs the FBAT tests, as well as gTDT and RV-TDT BRV (see below for references)

how to run simulations: (number of replicates set to 100 in main.cpp)

*./main seed number_trios number_causal_variants index_1_causal_variant ... index_last_causal_variant effect_size_1 ... effect_size_last*

**CAUTION: no index checking, relies heavily on reasonable input data. No input checking!**

example scenario 1 in the manuscript:

*./main 1 1000 2 494 924 1.8 1.8*

**gTDT:** Chen,R. et al. (2015) A haplotype-based framework for group-wise transmission/disequilibrium tests for rare variant association analysis. Bioinforma. Oxf. Engl., 31, 1452–1459.  
**RV-TDT:** He,Z. et al. (2014) Rare-variant extensions of the transmission disequilibrium test: application to autism exome sequence data. Am. J. Hum. Genet., 94, 33–46.  
**1000 Genomes:** A global reference for human genetic variation, The 1000 Genomes Project Consortium, Nature 526, 68-74 (01 October 2015) doi:10.1038/nature15393  
