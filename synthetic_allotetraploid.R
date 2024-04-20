#This script creates a synthetic allotetraploid
#Read the data in:
arenosa_632 <- read.table(file ='arenosa_632.txt' ,header = TRUE,sep = '\t')
lyrata_272 <- read.table(file ='lyrata_272_with_some_hybrids.txt' ,header = TRUE,sep = '\t')

#Merge the two data frames based on the 'POS' column:
merged_df <- merge(arenosa_632, lyrata_272, by = "POS", suffixes = c("_arenosa", "_lyrata"))

#Calculate the mean of the 'AF' column for matching rows:
#merged_df$Mean_AF <- (merged_df$AF_arenosa + merged_df$AF_lyrata) / 2

#Create a new data frame with 'POS' and 'Mean_AF':
new_df <- merged_df[, c("POS", "Mean_AF")]

#This new data frame has the mean allele frequencies at the same sites for lyrata and arenosa.

#Plot:
ggplot(data = new_df, aes(Mean_AF)) +
  geom_histogram(color='black',fill='white', bins = 100)

#Filter the data as the uncommon SNPs conceal potential trends in the data (you can see a slight peak in the middle however its masked because theres such a high count of in-frequent alleles), so filter the mean allele frequency to be greater than 0.1:

filtered_df <- new_df[new_df$Mean_AF > 0.1, ]

#Now plot:
ggplot(data = filtered_df, aes(Mean_AF)) +
  geom_histogram(color='black',fill='white', bins = 100)
