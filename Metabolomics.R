# Load required libraries
library(magrittr)  # For piping operations (%>%)
library(ggplot2)   # For advanced data visualization
library(ggrepel)   # For avoiding overlapping labels in plots
library(tidyverse) # For data manipulation and cleaning
library(limma)     # For differential expression analysis
library(vegan)     # For ecological data analysis (e.g., PCA)
library(cluster)   # For clustering algorithms
library(factoextra)# For visualizing clustering results
library(gridExtra) # For arranging multiple plots in a grid
library(PerformanceAnalytics) # For performance metrics
library(corrplot)  # For correlation matrix visualization
library(Hmisc)     # For statistical functions
library(RColorBrewer) # For color palettes
library(impute)    # For imputing missing data
library(glmnet)    # For regularized regression models

# Load metabolomic data
data <- read.csv('liver_cancer_metabolomic_data.csv',
                 row.names = 1, stringsAsFactors = FALSE, 
                 na.strings = c("NA", "NaN", ""), check.names = FALSE)
dim(data)  # Check dimensions of the data
head(data) # Display first few rows
is.numeric(data) # Check if data is numeric
str(data)  # Display structure of the data

# QC Data: Convert data into numeric format
int.mat <- data
num.mat <- matrix(as.numeric(gsub(",", "", as.matrix(int.mat))),
                  nrow = nrow(int.mat),
                  dimnames = dimnames(int.mat))
head(num.mat) # Display first few rows of numeric matrix
mode(num.mat) # Check the data type of the matrix
dim(data)     # Check dimensions of the original data

# Remove zero-variance features
dim(num.mat) # Check dimensions of numeric matrix
varCol <- apply(num.mat, 2, var, na.rm = TRUE) # Calculate variance for each column
constCol <- (varCol == 0 | is.na(varCol)) # Identify columns with zero or NA variance
sum(constCol) # Count number of constant columns
num.mat <- num.mat[, !constCol] # Remove constant columns
dim(num.mat) # Check dimensions after removal

# Compute missingness rate in the metabolites
nrow(num.mat) * ncol(num.mat) # Total number of data points
sum(is.na(num.mat)) # Total number of missing values
sum(is.na(num.mat)) / (nrow(num.mat) * ncol(num.mat)) * 100 # Percentage of missing values
round(sum(is.na(num.mat)) / (nrow(num.mat) * ncol(num.mat)) * 100, 1) # Rounded percentage

# Calculate missingness rate per metabolite
missingRatePerMetabolite <- (apply(is.na(num.mat), 2, sum) / nrow(num.mat)) * 100
min(missingRatePerMetabolite) # Minimum missingness rate
max(missingRatePerMetabolite) # Maximum missingness rate

# Histogram for missingness rate
options(repr.plot.width = 10, repr.plot.height = 8)
h <- hist((apply(is.na(num.mat), 2, sum) / nrow(num.mat)) * 100, breaks = 10,
          main = "Histogram for Missingness",
          xlab = "Percentage of missingness")
text(h$mids, h$counts, labels = h$counts, adj = c(0.5, -0.5)) # Add count labels
abline(v = 50, col = "red", lwd = 3, lty = 2) # Add vertical line at 50%

# Filter data to keep metabolites with <50% missingness
good.inx <- apply(is.na(num.mat), 2, sum) / nrow(num.mat) < 0.5
num.mat <- num.mat[, good.inx]
dim(num.mat) # Check dimensions after filtering

# Impute missing values using k-nearest neighbors (k=10)
num.mat.imputed <- impute.knn(num.mat, k = 10)$data
class(num.mat.imputed) # Check class of imputed data
head(num.mat.imputed)  # Display first few rows
dim(num.mat.imputed)   # Check dimensions after imputation

# Log-transform the imputed data
num.mat.imputed.logged <- log2(num.mat.imputed + 1)
head(num.mat.imputed.logged) # Display first few rows of log-transformed data
write.csv(num.mat.imputed.logged, 'Normalized_data.csv') # Save normalized data

# Density plots to compare distributions before and after log transformation
options(repr.plot.width = 15, repr.plot.height = 8)
par(mfrow = c(1, 2))
plot(density(num.mat.imputed[, 1]), main = 'Before log2')
plot(density(num.mat.imputed.logged[, 1]), main = 'After log2')

# Boxplots to visualize distributions before and after log transformation
options(repr.plot.width = 10, repr.plot.height = 8)
par(mar = c(10, 5, 2, 5), mfrow = c(1, 2))
boxplot(num.mat.imputed[, 1:10], main = "Before log2", horizontal = TRUE, 
        names = colnames(num.mat.imputed)[1:10], las = 2, col = "lightgreen")
boxplot(num.mat.imputed.logged[, 1:10], main = "After log2", horizontal = TRUE, 
        names = colnames(num.mat.imputed.logged)[1:10], las = 2, col = "lightgreen")

# Scale the log-transformed data
num.mat.imputed.logged.scaled <- scale(num.mat.imputed.logged, center = TRUE, scale = TRUE)
head(num.mat.imputed.logged.scaled) # Display first few rows of scaled data

# PCA for dimensionality reduction and visualization
df_pca <- prcomp(Normalized_data)
df_out <- as.data.frame(df_pca$x)
ggplot(df_out, aes(x = PC1, y = PC2, color = samples_meta$Stage, shape = samples_meta$Gender)) +
  geom_point() + ggtitle("") + labs(color = '') +
  geom_point(size = 8, alpha = 0.5) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 0.5, vjust = 0.5, face = "plain"),
        axis.text.y = element_text(size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0, face = "bold"),
        axis.title.y = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, face = "plain"),
        strip.text = element_text(size = 15, face = "plain"),
        legend.position = "right",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        axis.line = element_line(colour = "black")) +
  xlab(paste0("PC 1 (", round(df_pca$sdev[1], 1), "%)")) +
  ylab(paste0("PC 2 (", round(df_pca$sdev[2], 1), "%)"))

# K-means clustering
kmeans2 <- kmeans(Normalized_data, centers = 2, nstart = 25)
kmeans3 <- kmeans(Normalized_data, centers = 3, nstart = 25)
kmeans4 <- kmeans(Normalized_data, centers = 4, nstart = 25)
kmeans5 <- kmeans(Normalized_data, centers = 5, nstart = 25)

# Visualize clustering results
plot1 <- fviz_cluster(kmeans2, geom = "point", data = Normalized_data) + ggtitle("k = 2")
plot2 <- fviz_cluster(kmeans3, geom = "point", data = Normalized_data) + ggtitle("k = 3")
plot3 <- fviz_cluster(kmeans4, geom = "point", data = Normalized_data) + ggtitle("k = 4")
plot4 <- fviz_cluster(kmeans5, geom = "point", data = Normalized_data) + ggtitle("k = 5")
grid.arrange(plot1, plot2, plot3, plot4, nrow = 2)

# Hierarchical clustering
d <- dist(Normalized_data, method = "euclidean") # Distance matrix
fit <- hclust(d, method = "ward") # Hierarchical clustering
plot(fit) # Display dendrogram
groups <- cutree(fit, k = 5) # Cut tree into 5 clusters
rect.hclust(fit, k = 5, border = "red") # Draw borders around clusters

# Differential expression analysis using limma
design <- model.matrix(~0 + factor(type) + samples_meta$Gender)
colnames(design) <- c(levels(factor(type)), 'Group')
contrast <- makeContrasts(Cancer - Healthy, levels = design)
fit <- lmFit(as.matrix(t(Normalized_data)), design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# Extract differentially expressed metabolites (DEGs)
DEGs <- topTable(fit2, adjust.method = 'fdr', number = Inf, p.value = 1, coef = 1)
dim(DEGs) # Check dimensions of DEGs
head(DEGs) # Display first few rows
tail(DEGs) # Display last few rows

# Filter DEGs based on p-value and log fold change
DEGs %>% filter(P.Value < 0.05 & logFC > 1) %>% dim  # Upregulated in Cancer
DEGs %>% filter(P.Value < 0.05 & logFC < -1) %>% dim # Downregulated in Cancer

# Save DEGs to a CSV file
write.csv(DEGs, "DEGS.csv")

# Pathway analysis
overlapping <- read.delim("C:/My pc/Egcombio/MODA/Metabolomics/ORA_results.tab")
overlapping_filtered <- overlapping

# Visualize pathway analysis results
# Pathway analysis visualization
p6 <- ggplot() +
  geom_point(data = overlapping_filtered, 
             mapping = aes(x = overlapping_filtered$size,
                           y = -log(overlapping_filtered$q.value),
                           color = overlapping_filtered$source,
                           size = overlapping_filtered$size)) +
  scale_size(range = c(10, 30), guide = 'none') +
  labs(x = "Size of genes pathway", y = "-Log(q-value)", 
       color = "Pathway source", size = "# of overlaped metabolites") +
  ggtitle("Pathways analysis") +
  geom_label_repel(aes(x = overlapping_filtered$size,
                       y = -log(overlapping_filtered$q.value), 
                       color = overlapping_filtered$source,
                       label = str_wrap(overlapping_filtered$pathway, width = 20)),
                   min.segment.length = unit(2, 'lines'),
                   size = 3.5, force = 1, 
                   arrow = arrow(length = unit(0.02, "npc")),
                   segment.color = 'red',
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.5, "lines"), show.legend = FALSE) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"))

print(p6)

# Save pathway analysis results
write.csv(overlapping, "liver_cancer_normal_LIMMA_adj_0.05_1.csv")