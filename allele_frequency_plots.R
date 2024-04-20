#This script creates a histogram for allele frequencies greater than 0.1 for a population you specify

#Create a populations file:
individual_names <- indNames(aa.genlight)
populations <- as.character(pop(aa.genlight))
data <- data.frame(individual_names, populations)
write.table(data, "pops.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

#In a Unix environment:
#Remove the quotation marks:
sed 's/"//g' pops.txt > populations.txt

#Execute the poly_freq.c script created by Tuomas:
##Compile the scipt into an executable environment called poly_freq:
gcc poly_freq.c -o poly_freq -lm

##Execute the script:
./poly_freq -vcf filtered_tetraploids.vcf -pops populations.txt > info.tsv

#This will create a file with population-specific allele frequencies

#Now read the tsv file into R:
df <- read.table(file ='info.tsv', header = TRUE, sep = '\t')

#Store allele frequencies of the population you want to plot in a variable:
allele_frequencies <- df$HAB[df$HAB > 0.1]

#Plot a histogram of the results

ggplot(data = data.frame(allele_frequencies), aes(x = allele_frequencies)) +
  geom_histogram(color = 'black', fill = 'white', bins = 15) +
  labs(x = "Allele Frequencies", y = "Count", title = "AFS of HAB")

#Repeat the previous two steps for each population by changing the 3 instances of 'HAB'. 
