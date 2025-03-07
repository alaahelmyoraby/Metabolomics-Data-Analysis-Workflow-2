This R script is designed to perform a comprehensive analysis of **liver cancer metabolomic data**. The goal is to preprocess, analyze, and visualize the data to identify key metabolites and pathways associated with liver cancer. The analysis includes the following steps:

1. **Data Loading and Inspection**: The metabolomic data is loaded from a CSV file, and basic checks are performed to understand its structure and content.
2. **Data Preprocessing**: The data is cleaned, missing values are imputed, and zero-variance features are removed to ensure high-quality input for downstream analysis.
3. **Data Transformation**: The data is log-transformed to stabilize variance and normalize the distribution.
4. **Exploratory Data Analysis (EDA)**: Density plots and boxplots are used to visualize the distribution of metabolites before and after transformation.
5. **Dimensionality Reduction**: Principal Component Analysis (PCA) is performed to reduce the dimensionality of the data and identify patterns.
6. **Clustering Analysis**: Both k-means and hierarchical clustering are applied to group similar metabolites and samples.
7. **Differential Expression Analysis**: The `limma` package is used to identify differentially expressed metabolites (DEGs) between cancer and healthy samples.
8. **Pathway Analysis**: Enriched pathways are identified and visualized to understand the biological significance of the results.

This script leverages a variety of R packages, including `tidyverse` for data manipulation, `ggplot2` for visualization, `limma` for differential expression analysis, and `factoextra` for clustering. The results are saved as CSV files for further exploration and interpretation.

---

### **Key Objectives**
- **Preprocess raw metabolomic data** to handle missing values, remove noise, and normalize the data.
- **Identify differentially expressed metabolites** associated with liver cancer.
- **Visualize patterns** in the data using clustering and dimensionality reduction techniques.
- **Perform pathway analysis** to uncover biological insights.

---

### **Input Files**
1. **`liver_cancer_metabolomic_data.csv`**: Contains the raw metabolomic data.
2. **`liver_cancer_sample_metadata.csv`**: Contains metadata for the samples (e.g., cancer stage, gender).
3. **`liver_cancer_metabolite_metadata.csv`**: Contains metadata for the metabolites (e.g., biochemical names).
4. **`ORA_results.tab`**: Contains pathway analysis results.

---

### **Output Files**
1. **`Normalized_data.csv`**: Log-transformed and normalized metabolomic data.
2. **`DEGS.csv`**: Differentially expressed metabolites (DEGs).
3. **`liver_cancer_normal_LIMMA_adj_0.05_1.csv`**: Pathway analysis results.

### **1. Loading Required Libraries**
The script begins by loading the necessary R libraries, which provide functions for data manipulation, visualization, statistical analysis, and more.

```R
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
```

---

### **2. Loading and Inspecting the Data**
The metabolomic data is loaded from a CSV file, and basic checks are performed to understand its structure and content.

```R
# Load metabolomic data
data <- read.csv('liver_cancer_metabolomic_data.csv',
                 row.names = 1, stringsAsFactors = FALSE, 
                 na.strings = c("NA", "NaN", ""), check.names = FALSE)
dim(data)  # Check dimensions of the data
head(data) # Display first few rows
is.numeric(data) # Check if data is numeric
str(data)  # Display structure of the data
```

---

### **3. Data Preprocessing**
#### **3.1. Converting Data to Numeric**
The data is converted to a numeric matrix to ensure compatibility with downstream analysis.

```R
# QC Data: Convert data into numeric format
int.mat <- data
num.mat <- matrix(as.numeric(gsub(",", "", as.matrix(int.mat))),
                  nrow = nrow(int.mat),
                  dimnames = dimnames(int.mat))
head(num.mat) # Display first few rows of numeric matrix
mode(num.mat) # Check the data type of the matrix
dim(data)     # Check dimensions of the original data
```

#### **3.2. Removing Zero-Variance Features**
Columns with zero variance are removed to reduce noise in the data.

```R
# Remove zero-variance features
dim(num.mat) # Check dimensions of numeric matrix
varCol <- apply(num.mat, 2, var, na.rm = TRUE) # Calculate variance for each column
constCol <- (varCol == 0 | is.na(varCol)) # Identify columns with zero or NA variance
sum(constCol) # Count number of constant columns
num.mat <- num.mat[, !constCol] # Remove constant columns
dim(num.mat) # Check dimensions after removal
```

#### **3.3. Handling Missing Data**
The script calculates the percentage of missing values and visualizes the missingness rate.

```R
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
```

#### **3.4. Filtering Data**
Metabolites with more than 50% missing values are removed.

```R
# Filter data to keep metabolites with <50% missingness
good.inx <- apply(is.na(num.mat), 2, sum) / nrow(num.mat) < 0.5
num.mat <- num.mat[, good.inx]
dim(num.mat) # Check dimensions after filtering
```

#### **3.5. Imputing Missing Values**
Missing values are imputed using the k-nearest neighbors (k-NN) algorithm.

```R
# Impute missing values using k-nearest neighbors (k=10)
num.mat.imputed <- impute.knn(num.mat, k = 10)$data
class(num.mat.imputed) # Check class of imputed data
head(num.mat.imputed)  # Display first few rows
dim(num.mat.imputed)   # Check dimensions after imputation
```

#### **3.6. Log Transformation**
The data is log-transformed to stabilize variance and normalize the distribution.

```R
# Log-transform the imputed data
num.mat.imputed.logged <- log2(num.mat.imputed + 1)
head(num.mat.imputed.logged) # Display first few rows of log-transformed data
write.csv(num.mat.imputed.logged, 'Normalized_data.csv') # Save normalized data
```

---

### **4. Data Visualization**
#### **4.1. Density Plots**
Density plots are used to compare the distribution of data before and after log transformation.

```R
# Density plots to compare distributions before and after log transformation
options(repr.plot.width = 15, repr.plot.height = 8)
par(mfrow = c(1, 2))
plot(density(num.mat.imputed[, 1]), main = 'Before log2')
plot(density(num.mat.imputed.logged[, 1]), main = 'After log2')
```

#### **4.2. Boxplots**
Boxplots are generated to visualize the distribution of metabolites.

```R
# Boxplots to visualize distributions before and after log transformation
options(repr.plot.width = 10, repr.plot.height = 8)
par(mar = c(10, 5, 2, 5), mfrow = c(1, 2))
boxplot(num.mat.imputed[, 1:10], main = "Before log2", horizontal = TRUE, 
        names = colnames(num.mat.imputed)[1:10], las = 2, col = "lightgreen")
boxplot(num.mat.imputed.logged[, 1:10], main = "After log2", horizontal = TRUE, 
        names = colnames(num.mat.imputed.logged)[1:10], las = 2, col = "lightgreen")
```

#### **4.3. Scaling Data**
The data is scaled to have zero mean and unit variance.

```R
# Scale the log-transformed data
num.mat.imputed.logged.scaled <- scale(num.mat.imputed.logged, center = TRUE, scale = TRUE)
head(num.mat.imputed.logged.scaled) # Display first few rows of scaled data
```

---

### **5. Principal Component Analysis (PCA)**
PCA is performed to reduce dimensionality and visualize patterns in the data.

```R
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
```

---

### **6. Clustering Analysis**
#### **6.1. K-Means Clustering**
K-means clustering is performed for different numbers of clusters.

```R
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
```

#### **6.2. Hierarchical Clustering**
Hierarchical clustering is performed using the Ward method.

```R
# Hierarchical clustering
d <- dist(Normalized_data, method = "euclidean") # Distance matrix
fit <- hclust(d, method = "ward") # Hierarchical clustering
plot(fit) # Display dendrogram
groups <- cutree(fit, k = 5) # Cut tree into 5 clusters
rect.hclust(fit, k = 5, border = "red") # Draw borders around clusters
```

---

### **7. Differential Expression Analysis**
Differential expression analysis is performed using the `limma` package.

```R
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
```

---

### **8. Pathway Analysis**
Pathway analysis is performed to identify enriched pathways.

```R
# Pathway analysis
overlapping <- read.delim("C:/My pc/Egcombio/MODA/Metabolomics/ORA_results.tab")
overlapping_filtered <- overlapping

# Visualize pathway analysis results
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
```

---

### **Key Notes**
1. **No changes were made to the original code.**
2. **Comments were added to explain each step and its purpose.**
3. Ensure all required input files are in the correct paths before running the script.
