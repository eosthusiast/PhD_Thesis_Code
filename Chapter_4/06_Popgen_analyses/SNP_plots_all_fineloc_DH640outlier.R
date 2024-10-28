#install.packages("tidyverse", repos = "https://cloud.r-project.org/")
#library(tidyverse)
library(ggplot2)
library(readr)
library(tibble)
library(dplyr)
library(vcfR)
library(adegenet)
library(poppr)
library(viridis)
library(pheatmap)
library(geosphere)
library(vegan)
library(grid)
library(gridExtra)

setwd("/nesi/nobackup/ga03488/David_H/Ppo_Popgen/14_whole_genome_snpcalling/03_SNP_subset_80/")

# Run script in command line: Rscript SNP_plots_ArthursPass.R *.vcf

# script.R
args <- commandArgs(trailingOnly = TRUE)
my_vcf <- args[1]  # First argument passed to the script
my_vcf <- "C:/Users/herad/OneDrive - MWLR/Fungi, Food and Fodder/03 - Molecular work/PopGen Bioinformatics/04_WG_SNPs/variants_snpeff_1k_WG_80.vcf"
basename_vcf <- tools::file_path_sans_ext(basename(my_vcf))  # Get the basename without extension

# Print or use the variable
print(paste("The variable passed is:", basename_vcf))

## Import data
WG_vcf_80 <- read.vcfR(my_vcf, verbose = FALSE)
#sample_valleys <- as.data.frame(read_tsv("my.strata_valleys.tsv"))
#sample_valleys_only <- sample_valleys$POPULATION
#Ppo_coords <- read_csv("Ppo_Coordinates.csv")
# Remove DH641 from geodist
#Ppo_coords <- Ppo_coords[-13,]
#sample_ids <- Ppo_coords$ID
# Subset for AP
#Ppo_coords <- Ppo_coords[c(-1, -5, -6, -23, -24, -25, -28, -27, -12),]
# Generate a vector with sample IDs from column names
#sample_ids_AP <- sample_ids
sample_valleys_only[28] <- "Klondyke"
sample_valleys28 <- sample_valleys[-13,]
sample_valleys28_only <- sample_valleys_only[-13]

# Create genome data objects
WG_vcf_80.gl <- vcfR2genlight(WG_vcf_80)
WG_vcf_80.gid <- vcfR2genind(WG_vcf_80)

# Load strata
strata(WG_vcf_80.gid) <- sample_valleys28
pop(WG_vcf_80.gl) <- sample_valleys28_only
pop(WG_vcf_80.gid) <- sample_valleys28_only

# Calculate PCA
WG_vcf_80.pca <- glPca(WG_vcf_80.gl, nf = 10, scale = TRUE)
WG_vcf_80.pca

# convert scores of WG_vcf_80.pca into a tibble
WG_vcf_80.pca.scores <- as_tibble(WG_vcf_80.pca$scores)

# add the country data into a column of WG_vcf_80.pca.scores tibble
WG_vcf_80.pca.scores$site <- sample_valleys28_only

# We will also determine the variance each PC contributes the data, which will help us understand potential drivers of patterns in our dataset. Lets plot the eigenvectors to try an understand this a bit more.

png("C:/Users/herad/OneDrive - MWLR/Fungi, Food and Fodder/03 - Molecular work/PopGen Bioinformatics/04_WG_SNPs/eig_80.png")
eigenplot <- barplot(100 * WG_vcf_80.pca$eig / sum(WG_vcf_80.pca$eig), col="darkgrey")
eigenplot <- title(ylab = "Percent of variance explained (%)") 
eigenplot <- title(xlab = "Eigenvalues")
eigenplot
dev.off()

# Lets extract the variance associated with the top 4 PCs, so we can use them in our plots.

# first we sum all the eigenvalues
eig.total <- sum(WG_vcf_80.pca$eig)

# sum the variance
PC1.variance <- formatC(head(WG_vcf_80.pca$eig)[1]/eig.total * 100)
PC2.variance <- formatC(head(WG_vcf_80.pca$eig)[2]/eig.total * 100)
PC3.variance <- formatC(head(WG_vcf_80.pca$eig)[3]/eig.total * 100)
PC4.variance <- formatC(head(WG_vcf_80.pca$eig)[4]/eig.total * 100)

xlabel <- paste0("PC1 variance = ",PC1.variance,"%")
ylabel <- paste0("PC2 variance = ", PC2.variance, "%")

plot12 <- ggplot(WG_vcf_80.pca.scores, aes(PC1, PC2, color = site)) + 
  geom_point() + 
  geom_text(aes(label = sample_ids), vjust = -0.5, hjust = 0.5) +  
  scale_color_viridis(discrete = TRUE, end = 0.9, option = "B") +
  labs(x = xlabel, y = ylabel, title = paste(basename_vcf, "- PC1 vs PC2"))

plot12

# Lets quickly look at PC3/PC4, and compare to the first plot.

plot34 <- ggplot(WG_vcf_80.pca.scores, aes(PC3, PC4, color = site)) + 
    geom_point() + geom_text(aes(label = sample_ids), vjust = -0.5, hjust = 0.5) +  
  scale_color_viridis(discrete = TRUE, end = 0.9, option = "B") +
    labs(x = paste0("PC3 variance = ", PC3.variance,"%"), y = paste0("PC4 variance = ", PC4.variance, "%"), 
       title = paste(basename_vcf, "- PC3 vs PC4"))
plot34

# Calculate the mean value of the principal components for each country. 
# We can use this to make some labels for our plots
means <- WG_vcf_80.pca.scores %>% group_by(site) %>% 
  dplyr:::summarize(meanPC1 = mean(PC1), meanPC2 = mean(PC2),meanPC3 = mean(PC3), meanPC4 = mean(PC4))

plot12ellipse <- ggplot(WG_vcf_80.pca.scores, aes(PC1, PC2, color = site)) + 
  geom_point() + 
  geom_text(aes(label = sample_ids), vjust = -0.5, hjust = 0.5) +  
  scale_color_viridis(discrete = TRUE, end = 0.9, option = "B") +
  labs(x = xlabel, y = ylabel, title = paste(basename_vcf, "- PC1 vs PC2 (Ellipse)")) +
  stat_ellipse(level = 0.95, size = 1)

plot12ellipse

plot34ellipse <- ggplot(WG_vcf_80.pca.scores, aes(PC3, PC4, color = site)) + 
    geom_point() + geom_text(aes(label = sample_ids), vjust = -0.5, hjust = 0.5) +  
  scale_color_viridis(discrete = TRUE, end = 0.9, option = "B") +
    labs(x = paste0("PC3 variance = ", PC3.variance,"%"), y = paste0("PC4 variance = ", PC4.variance, "%"), 
       title = paste(basename_vcf, "- PC3 vs PC4 (Ellipse)")) +
  stat_ellipse(level = 0.95, size = 1)
plot34ellipse

## Investigate outlier in PC1: DH640
which.max(-1*(WG_vcf_80.pca$loadings[,1]))
plot(WG_vcf_80.pca$loadings[,1:2], xlim = c(-0.02, -0.01))
plot(jitter(WG_vcf_80.pca$loadings[, 1]), jitter(WG_vcf_80.pca$loadings[, 2]), 
     xlim = c(-0.02, -0.01))

# Get the loadings for the first principal component
loadings_pc1 <- WG_vcf_80.pca$loadings[, 1]

# Get the indices of the top 10 highest absolute values
top10_indices <- order(abs(loadings_pc1), decreasing = TRUE)[1:50]

# Extract the top 10 highest loadings
top10_loadings <- loadings_pc1[top10_indices]

# Extract corresponding locus names
locus_names <- locNames(WG_vcf_80.gl)[top10_indices]

# Create a barplot to visualize the top 10 highest loadings
barplot(top10_loadings, names.arg = locus_names, las = 2, col = "blue", main = "Top 50 Highest Loadings for PC1",
        ylab = "Loading Value", xlab = "", cex.names = 0.5)


# Define the loci of interest
loci_of_interest <- c("k141_10516_1935", "k141_28181_4523")

# Extract genotype information for the specific loci
genotype_info <- WG_vcf_80.gid@tab[, loci_of_interest]

# Display the genotype information
genotype_info

# Separate loci into a list of genind objects
loci_list <- seploc(WG_vcf_80.gid)

# Extract specific loci of interest using their names
locus_1 <- tab(loci_list[["k141_10516_1935"]])
locus_2 <- tab(loci_list[["k141_28181_4523"]])

# Combine the loci information into a data frame
genotype_data <- data.frame(locus_1, locus_2)

# Assign column names for easier interpretation
colnames(genotype_data) <- c("k141_10516_1935_allele_1", "k141_10516_1935_allele_2", "k141_28181_4523_allele_1", "k141_28181_4523_allele_2")

# Display the genotype data
print(genotype_data)

# Optional: Calculate the average allele frequencies per population
freq_locus_1 <- apply(locus_1, 2, function(e) tapply(e, (WG_vcf_80.gid@strata$INDIVIDUALS), mean, na.rm = TRUE))
freq_locus_2 <- apply(locus_2, 2, function(e) tapply(e, (WG_vcf_80.gid@strata$INDIVIDUALS), mean, na.rm = TRUE))

# Display allele frequencies for both loci
print(freq_locus_1)
print(freq_locus_2)



# DAPC
dapc.coarse<- dapc(WG_vcf_80.gid, var.contrib = TRUE, scale = FALSE, n.pca = nPop(WG_vcf_80.gid) - 1, n.da = nPop(WG_vcf_80.gid) - 1)
scatter(dapc.coarse, scree.pca=TRUE, scree.da=TRUE, cell = 0, pch = 19, cstar = 1, mstree = TRUE, lwd = 1, lty = 1, clabel = 1, legend=F)
contrib <- loadingplot(dapc.coarse$var.contr, axis = 2, thres = 0.00015, lab.jitter = 1)

summary(dapc.coarse)
predict(dapc.coarse)
a.score(dapc.coarse)
optim.a.score(dapc.coarse)

crossval <- xvalDapc(WG_vcf_80.gid, grp = pop(WG_vcf_80.gid), n.pca.max = 20)


# cross validation identified 3 axes to best explain variation
scatter(crossval$DAPC, scree.pca=TRUE, scree.da=TRUE, cell = 0, pch = 19, cstar = 1, mstree = TRUE, lwd = 1, lty = 1, clabel = 1, legend=F)
contrib <- loadingplot(crossval$DAPC$var.contr, axis = 2, thres = 0.00013, lab.jitter = 1)
summary(crossval$DAPC)
predict(crossval$DAPC) # posterior assignment into 5 groups, no more Otago, EastCoast, Nina or Fiordland samples. some strains switched between the other groups
a.score(crossval$DAPC, n.sim= 999) # very weak, indicates either weak discrimination or instability of the results

# check heatmap with outlier removed
# Get the list of locus names
all_loci <- locNames(WG_vcf_80.gid)

homo_alt_outliers <- c("k141_10516_1935",  "k141_21593_50810", "k141_17821_19887", "k141_4121_920", "k141_9751_2147")

# Exclude the specific locus to be removed
loci_to_keep <- all_loci[all_loci != homo_alt_outliers]

# Subset the genind object to keep only the desired loci
WG_vcf_80.gid_filtered <- WG_vcf_80.gid[loc = loci_to_keep]

# Subset the genlight object to keep only the desired loci
WG_vcf_80.gl_filtered <- WG_vcf_80.gl[, loci_to_keep]

# Display information about the filtered genlight object
WG_vcf_80.gl_filtered

# Display information about the filtered genind object
WG_vcf_80.gid_filtered

# Separate loci into a list of genind objects
#loci_list_filtered <- seploc(WG_vcf_80.gid_filtered)

# Extract specific loci of interest using their names
#filtertest <- tab(loci_list_filtered[["k141_10516_1934"]])

nei_dist_filtered <- nei.dist(WG_vcf_80.gid_filtered)

# Check PCA with k141_10516_1935 removed
# Calculate PCA
WG_vcf_80.pca <- glPca(WG_vcf_80.gl_filtered, nf = 10, scale = TRUE)
WG_vcf_80.pca

# convert scores of WG_vcf_80.pca into a tibble
WG_vcf_80.pca.scores <- as_tibble(WG_vcf_80.pca$scores)

# add the country data into a column of WG_vcf_80.pca.scores tibble
WG_vcf_80.pca.scores$site <- sample_valleys28_only

# We will also determine the variance each PC contributes the data, which will help us understand potential drivers of patterns in our dataset. Lets plot the eigenvectors to try an understand this a bit more.

eigenplot <- barplot(100 * WG_vcf_80.pca$eig / sum(WG_vcf_80.pca$eig), col="green")
eigenplot <- title(main = paste(basename_vcf, " - Eigenvalues"))  # Add VCF file basename in the title
eigenplot <- title(ylab = "Percent of variance explained") 
eigenplot <- title(xlab = "Eigenvalues")
eigenplot

# Lets extract the variance associated with the top 4 PCs, so we can use them in our plots.

# first we sum all the eigenvalues
eig.total <- sum(WG_vcf_80.pca$eig)

# sum the variance
PC1.variance <- formatC(head(WG_vcf_80.pca$eig)[1]/eig.total * 100)
PC2.variance <- formatC(head(WG_vcf_80.pca$eig)[2]/eig.total * 100)
PC3.variance <- formatC(head(WG_vcf_80.pca$eig)[3]/eig.total * 100)
PC4.variance <- formatC(head(WG_vcf_80.pca$eig)[4]/eig.total * 100)

xlabel <- paste0("PC1 variance = ",PC1.variance,"%")
ylabel <- paste0("PC2 variance = ", PC2.variance, "%")

plot12 <- ggplot(WG_vcf_80.pca.scores, aes(PC1, PC2, color = site)) + 
  geom_point() + 
  geom_text(aes(label = sample_ids), vjust = -0.5, hjust = 0.5) +  
  scale_color_viridis(discrete = TRUE, end = 0.9, option = "B") +
  labs(x = xlabel, y = ylabel, title = paste(basename_vcf, "- PC1 vs PC2"))

plot12

### Heatmap 
nei_dist_filtered_matrix <- as.matrix(nei_dist_filtered)

nei_dist_filtered_heatmap <- pheatmap(nei_dist_filtered_matrix, 
                             cluster_rows = TRUE, 
                             cluster_cols = TRUE, 
                             color = rev(viridis(256)),  # Use viridis color palette
                             labels_row = sample_ids,  # Label the rows
                             labels_col = sample_ids,  # Label the rows
                             show_rownames = TRUE, 
                             show_colnames = TRUE, 
                             main = paste(basename_vcf, "- Nei's distance"),
                             legend = TRUE,  # Add legend
                             border_color = NA)  # Remove border color


# Create Heatmap, Mantel test and geography vs genomic distance matrix
# Calculate Nei's genetic distance
nei_dist <- nei.dist(WG_vcf_80.gid)

### Heatmap 
nei_dist_matrix <- as.matrix(nei_dist)

nei_dist_heatmap <- pheatmap(nei_dist_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         color = rev(viridis(256)),  # Use viridis color palette
         labels_row = sample_ids,  # Label the rows
         labels_col = sample_ids,  # Label the rows
         show_rownames = TRUE, 
         show_colnames = TRUE, 
        main = paste(basename_vcf, "- Nei's distance"),
         legend = TRUE,  # Add legend
         border_color = NA)  # Remove border color


###
# Calculate proportion of homozygous sites for each sample

# Extract genotype matrix from the genlight object
geno_matrix <- as.matrix(WG_vcf_80.gl)

# Function to calculate homozygosity for a sample
calculate_homozygosity <- function(genotypes) {
  # Count loci that are homozygous (0 or 2) and divide by total number of loci
  homozygous_count <- sum(genotypes == 0 | genotypes == 2, na.rm = TRUE)
  total_loci <- length(genotypes) - sum(is.na(genotypes))
  homozygosity <- homozygous_count / total_loci
  return(homozygosity)
}

# Apply the function to each sample (row)
homozygosity_per_sample <- apply(geno_matrix, 1, calculate_homozygosity)

# Display the results
homozygosity_per_sample

# Create the plot of homozygous alternate proportion per sample
png("C:/Users/herad/OneDrive - MWLR/Fungi, Food and Fodder/03 - Molecular work/PopGen Bioinformatics/04_WG_SNPs/homozygous_80.png", res = 100)
plot(homozygosity_per_sample, 
     main = "", 
     xlab = "Samples", 
     ylab = "Proportion of Homozygous Sites", 
     pch = 16, 
     col = "blue")

# Add rownames as labels to each point
text(x = 1:length(homozygosity_per_sample), 
     y = homozygosity_per_sample, 
     labels = WG_vcf_80.gid@strata$prefix, 
     pos = 3, cex = 0.6, col = "black")
dev.off()

# Function to calculate the proportion of homozygous alternate sites for a sample
calculate_homozygous_alt <- function(genotypes) {
  # Count loci that are homozygous alternate (value == 2)
  homozygous_alt_count <- sum(genotypes == 2, na.rm = TRUE)
  
  # Calculate the total number of loci that are non-missing
  total_loci <- length(genotypes) - sum(is.na(genotypes))
  
  # Calculate the proportion of homozygous alternate sites
  homozygous_alt_proportion <- homozygous_alt_count / total_loci
  
  return(homozygous_alt_proportion)
}

# Apply the function to each sample (row)
homozygous_alt_per_sample <- apply(geno_matrix, 1, calculate_homozygous_alt)

# Display the results
homozygous_alt_per_sample




# Geographic distance
geo_dist <- distm(Ppo_coords[, c("longitude", "latitude")], fun = distHaversine)

#  Mantel test
# Perform the Mantel test
mantel_result <- vegan:::mantel(nei_dist_matrix, geo_dist, method = "pearson", permutations = 9999)

# Add Mantel test result to the PDF as text
mantel_text <- capture.output(print(mantel_result))
mantelGrob <- grid.text(paste(mantel_text, collapse = "\n"), x = 0.1, y = 0.9, just = "left", gp = gpar(fontsize = 8), check.overlap=TRUE)

# Geographic vs genomic distance
# Convert matrices to vectors and remove diagonal
geo_vector <- geo_dist[lower.tri(geo_dist)]
geno_vector <- nei_dist_matrix[lower.tri(nei_dist_matrix)]

# Plot: see below in PDF script

# Calculate average distances for each sample
avg_geo_dist <- rowMeans(geo_dist)
avg_geno_dist <- rowMeans(nei_dist_matrix)

# Create a data frame for plotting
data <- data.frame(Sample = sample_ids_AP,
                   Avg_Geo_Distance = avg_geo_dist,
                   Avg_Geno_Similarity = avg_geno_dist)

# Plot
avg_geogenoplot <- ggplot(data, aes(x = Avg_Geo_Distance, y = Avg_Geno_Similarity, label = Sample)) +
  geom_point(color = "blue", size = 3) +
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
  labs(x = "Average Geographic Distance", y = "Average Genomic Similarity", 
       title = paste(basename_vcf, "- Avg. Geographic vs. Genomic Distances")) #+ xlim(0, 300000)
  theme_minimal()

## AMOVA
print(paste("Start AMOVA for", basename_vcf))

#load("data.RData")
amova_result <- poppr.amova(WG_vcf_80.gid, ~POPULATION, within = FALSE)
amova_text <- capture.output(print(amova_result))
amovaGrob <- grid.text(paste(amova_text, collapse = "\n"), x = 0.1, y = 0.9, just = "left", gp = gpar(fontsize = 8), check.overlap=TRUE)
amova.test <- randtest(amova_result, nrepet = 999) # Test for significance
plot(amova.test)

## Poppr stats
print(paste("Start calulating poppr stats for", basename_vcf))
# Poppr summary table
#load("poppr_p.RData")
#print(poppr_p)

poppr_p <- poppr(WG_vcf_80.gid, sample = 9, plot = FALSE, quiet = FALSE) 

# Save data
save(amova_result, poppr_p, file = "data.RData")

# Saving to PDF
pdf(paste("Plots_all_fine_", basename_vcf, ".pdf", sep=""))
eigenplot
print(plot12)
print(plot34)
print(plot12ellipse)
print(plot34ellipse)
grid.newpage()  # Force a new page
nei_dist_heatmap
grid.newpage()  # Force a new page
tt <- ttheme_default(base_size = 6)
grid.text(paste(amova_text, collapse = "\n"), x = 0.1, y = 0.7, just = "left", gp = gpar(fontsize = 8), check.overlap=TRUE)
plot(amova.test)
# Add Poppr summary table to the PDF
grobPoppr <- tableGrob(poppr_p[,-15], theme=tt)
grid.arrange(grobPoppr, mantelGrob, ncol=1, top="Poppr summary stats", bottom="Mantel test results")
plot(geo_vector, geno_vector,
     xlab = "Geographic Distance",
     ylab = "Genomic Distance",
     main = paste(basename_vcf, "- Geographic vs. Genomic Distances"),
     pch = 19, col = "blue")
     #ylim = c(0.16,0.2))
print(avg_geogenoplot)
dev.off() 
