# 

# for BZD and LIC
# store allele frequencies in a variable:
allele_frequency_bzd <- df$BZD
allele_frequency_lic <- df$LIC

# do normality test <- not normally distributed
ad.test(allele_frequency_bzd)
ad.test(allele_frequency_lic)

# load in your allele frequencies produced from the poly freq tool for the populations you want to compare
# and store the allele frequencies in a variable. then create a data frame with the vector of allele frequencies so we filter it for intermediate allele frequencies
# for bzd 
intermediate_freq_bzd <- df$BZD[df$BZD>=0.45 & df$BZD <= 0.55]
intermediate_df_bzd <- data.frame(x=intermediate_freq_bzd)

# for lic
intermediate_freq_lic <- df$LIC[df$LIC>=0.45 & df$LIC <= 0.55]
intermediate_df_lic <- data.frame(x=intermediate_freq_lic)

# compare the distributions using the wilcoxon/mann whitney U test
# bzd vs lic <- significant difference
wilcox.test(intermediate_freq_bzd,intermediate_freq_lic)

# load in average AF values at sites with fixed differences in the in-silico allopolyploid. Fixed difference is >0.8 and <0.2 here
# to get only AF in range 0.4-0.6 which is the intermediate allele frequency range we are using
intermediate_freq_allopolyploid <- new_df[,2][new_df[,2]>=0.4 & new_df[,2] <= 0.6]

# get allele frequencies for allopolyploid
allele_frequency_allopolyploid <- new_df[,2]

# test for significance between our populations
# bzd vs allopolyploid <- significant difference
wilcox.test(intermediate_freq_bzd,intermediate_freq_allopolyploid)

# allopolyploid vs lic <- significant difference
wilcox.test(intermediate_freq_lic,intermediate_freq_allopolyploid)
