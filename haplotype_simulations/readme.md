## tool to perform the simulations with haplotype data

#assumes file "haps" is in this same folder

#run simulations (number of replicates set to 100 in main.cpp)

./main seed number_trios number_causal_variants index_1_causal_variant ... index_last_causal_variant effect_size_1 ... effect_size_last

#CAUTION: no index checking, relies heavily on reasonable input data. No input checking!

#example scenario 1 manuscript

./main 1 1000 2 494 924 1.8 1.8
