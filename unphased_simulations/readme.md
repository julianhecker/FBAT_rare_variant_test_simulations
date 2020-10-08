# genotype-based simulations

**tool to simulate family data according to the haplotype data and specifications described in the scenario file. Corresponds to the data simulated in the unphased simulation study described in "An exact, unifying framework for rare variant association testing in family-based designs, including higher criticism approaches, SKATs, multivariate and burden tests" (Hecker et al. 2020)**

- tool simulates and outputs data. Analysis can be performed with the corresponding tools. Data is generated for FBAT, gTDT, and RV-GDT (see references below).  

general usage:

*./main output_folder scenario_name_with_nvars seed scenario_file haplotype_file number_variants number_families number_replications scenario_name (for RV-GDT) number_haplotypes*

example scenario 1, p=30 variants

*./main output/ sig1_30 1 data/s1_sig1_30 data/hap_99ceu_91gbr_30 30 1000 1000 sig1 380*



**example scenario 1 file**

- number founder: 2 
- number offspring: 1
- nhaps: 380 
- nvars: 30 
- ncausal: 3 
- affection status of the members (founders first, then offspring): 0 0 1 
- offspring-founder relationship (member 3=offspring is offspring of founder 1 and 2): 1 2 3
- haplotype distribution founder 1: prob_hap_1 ... prob_hap_nvar
- haplotype distribution founder 2: prob_hap_1 ... prob_hap_nvar 
- indices causal variants (range 1 to nvar): 5 12 26
- effect sizes causal variants: 1.0319134 0.6319134 0.9115014
- which members to output (here: all): 1 1 1 
- sex information (for some tools required): 1 2 1


**RV-GDT:** He,Z. et al. (2014) Rare-variant extensions of the transmission disequilibrium test: application to autism exome sequence data. Am. J. Hum. Genet., 94, 33–46.  
**gTDT:** Chen,R. et al. (2015) A haplotype-based framework for group-wise transmission/disequilibrium tests for rare variant association analysis. Bioinforma. Oxf. Engl., 31, 1452–1459.  
