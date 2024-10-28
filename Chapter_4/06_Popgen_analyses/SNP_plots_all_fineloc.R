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
basename_vcf <- tools::file_path_sans_ext(basename(my_vcf))  # Get the basename without extension

# Print or use the variable
print(paste("The variable passed is:", basename_vcf))

## Import data
ArthursPass_vcf_80 <- read.vcfR(my_vcf, verbose = FALSE)
sample_valleys <- as.data.frame(read_tsv("my.strata_valleys.tsv"))
sample_valleys_only <- sample_valleys$POPULATION
Ppo_coords <- read_csv("Ppo_Coordinates.csv")
# Remove DH641 from geodist
Ppo_coords <- Ppo_coords[-13,]
sample_ids <- Ppo_coords$ID
# Subset for AP
#Ppo_coords <- Ppo_coords[c(-1, -5, -6, -23, -24, -25, -28, -27, -12),]
# Generate a vector with sample IDs from column names
sample_ids_AP <- sample_ids

# Create genome data objects
ArthursPass_vcf_80.gl <- vcfR2genlight(ArthursPass_vcf_80)
ArthursPass_vcf_80.gid <- vcfR2genind(ArthursPass_vcf_80)

# Load strata
strata(ArthursPass_vcf_80.gid) <- sample_valleys
pop(ArthursPass_vcf_80.gl) <- sample_valleys_only
pop(ArthursPass_vcf_80.gid) <- sample_valleys_only

# Calculate PCA
ArthursPass_vcf_80.pca <- glPca(ArthursPass_vcf_80.gl, nf = 10, scale = TRUE)
ArthursPass_vcf_80.pca

# convert scores of ArthursPass_vcf_80.pca into a tibble
ArthursPass_vcf_80.pca.scores <- as_tibble(ArthursPass_vcf_80.pca$scores)

# add the country data into a column of ArthursPass_vcf_80.pca.scores tibble
ArthursPass_vcf_80.pca.scores$site <- sample_valleys_only

# We will also determine the variance each PC contributes the data, which will help us understand potential drivers of patterns in our dataset. Lets plot the eigenvectors to try an understand this a bit more.

eigenplot <- barplot(100 * ArthursPass_vcf_80.pca$eig / sum(ArthursPass_vcf_80.pca$eig), col="green")
eigenplot <- title(main = paste(basename_vcf, " - Eigenvalues"))  # Add VCF file basename in the title
eigenplot <- title(ylab = "Percent of variance explained") 
eigenplot <- title(xlab = "Eigenvalues")
eigenplot

# Lets extract the variance associated with the top 4 PCs, so we can use them in our plots.

# first we sum all the eigenvalues
eig.total <- sum(ArthursPass_vcf_80.pca$eig)

# sum the variance
PC1.variance <- formatC(head(ArthursPass_vcf_80.pca$eig)[1]/eig.total * 100)
PC2.variance <- formatC(head(ArthursPass_vcf_80.pca$eig)[2]/eig.total * 100)
PC3.variance <- formatC(head(ArthursPass_vcf_80.pca$eig)[3]/eig.total * 100)
PC4.variance <- formatC(head(ArthursPass_vcf_80.pca$eig)[4]/eig.total * 100)

xlabel <- paste0("PC1 variance = ",PC1.variance,"%")
ylabel <- paste0("PC2 variance = ", PC2.variance, "%")

plot12 <- ggplot(ArthursPass_vcf_80.pca.scores, aes(PC1, PC2, color = site)) + 
  geom_point() + 
  geom_text(aes(label = sample_ids_AP), vjust = -0.5, hjust = 0.5) +  
  scale_color_viridis(discrete = TRUE, end = 0.9, option = "B") +
  labs(x = xlabel, y = ylabel, title = paste(basename_vcf, "- PC1 vs PC2"))

plot12

# Lets quickly look at PC3/PC4, and compare to the first plot.

plot34 <- ggplot(ArthursPass_vcf_80.pca.scores, aes(PC3, PC4, color = site)) + 
    geom_point() + geom_text(aes(label = sample_ids_AP), vjust = -0.5, hjust = 0.5) +  
  scale_color_viridis(discrete = TRUE, end = 0.9, option = "B") +
    labs(x = paste0("PC3 variance = ", PC3.variance,"%"), y = paste0("PC4 variance = ", PC4.variance, "%"), 
       title = paste(basename_vcf, "- PC3 vs PC4"))
plot34

# Calculate the mean value of the principal components for each country. 
# We can use this to make some labels for our plots
means <- ArthursPass_vcf_80.pca.scores %>% group_by(site) %>% 
  summarize(meanPC1 = mean(PC1), meanPC2 = mean(PC2),meanPC3 = mean(PC3), meanPC4 = mean(PC4))

plot12ellipse <- ggplot(ArthursPass_vcf_80.pca.scores, aes(PC1, PC2, color = site)) + 
  geom_point() + 
  geom_text(aes(label = sample_ids_AP), vjust = -0.5, hjust = 0.5) +  
  scale_color_viridis(discrete = TRUE, end = 0.9, option = "B") +
  labs(x = xlabel, y = ylabel, title = paste(basename_vcf, "- PC1 vs PC2 (Ellipse)")) +
  stat_ellipse(level = 0.95, size = 1)

plot12ellipse

plot34ellipse <- ggplot(ArthursPass_vcf_80.pca.scores, aes(PC3, PC4, color = site)) + 
    geom_point() + geom_text(aes(label = sample_ids_AP), vjust = -0.5, hjust = 0.5) +  
  scale_color_viridis(discrete = TRUE, end = 0.9, option = "B") +
    labs(x = paste0("PC3 variance = ", PC3.variance,"%"), y = paste0("PC4 variance = ", PC4.variance, "%"), 
       title = paste(basename_vcf, "- PC3 vs PC4 (Ellipse)")) +
  stat_ellipse(level = 0.95, size = 1)
plot34ellipse

# Create Heatmap, Mantel test and geography vs genomic distance matrix
# Calculate Nei's genetic distance
nei_dist <- nei.dist(ArthursPass_vcf_80.gid)

### Heatmap 
nei_dist_matrix <- as.matrix(nei_dist)

nei_dist_heatmap <- pheatmap(nei_dist_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         color = rev(viridis(256)),  # Use viridis color palette
         labels_row = sample_ids_AP,  # Label the rows
         labels_col = sample_ids_AP,  # Label the rows
         show_rownames = TRUE, 
         show_colnames = TRUE, 
        main = paste(basename_vcf, "- Nei's distance"),
         legend = TRUE,  # Add legend
         border_color = NA)  # Remove border color

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
amova_result <- poppr.amova(ArthursPass_vcf_80.gid, ~POPULATION, within = FALSE)
amova_text <- capture.output(print(amova_result))
amovaGrob <- grid.text(paste(amova_text, collapse = "\n"), x = 0.1, y = 0.9, just = "left", gp = gpar(fontsize = 8), check.overlap=TRUE)
amova.test <- randtest(amova_result, nrepet = 999) # Test for significance
plot(amova.test)

## Poppr stats
print(paste("Start calulating poppr stats for", basename_vcf))
# Poppr summary table
#load("poppr_p.RData")
#print(poppr_p)

poppr_p <- poppr(ArthursPass_vcf_80.gid, sample = 9, plot = FALSE, quiet = FALSE) 

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
