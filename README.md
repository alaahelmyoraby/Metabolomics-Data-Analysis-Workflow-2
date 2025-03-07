# Metabolomics-Data-Analysis-Workflow-2
This R script is a comprehensive analysis pipeline for metabolomic data related to liver cancer. It includes data preprocessing, quality control, normalization, visualization, statistical analysis, and pathway analysis. Below is an informative breakdown of the code, explaining each step without altering the code itself.

---

### **1. Loading Required Libraries**
The script begins by loading necessary R libraries for data manipulation, visualization, and statistical analysis:

```R
library(magrittr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(limma)
library(vegan)
library(cluster)
library(factoextra)
library(gridExtra)
library(PerformanceAnalytics)
library(corrplot)
library(Hmisc)
library(RColorBrewer)
library(impute)
library(glmnet)
```

These libraries provide functions for data wrangling (`tidyverse`, `magrittr`), visualization (`ggplot2`, `ggrepel`, `corrplot`), statistical analysis (`limma`, `vegan`, `glmnet`), clustering (`cluster`, `factoextra`), and handling missing data (`impute`).

---

### **2. Loading and Inspecting the Data**
The metabolomic data is loaded from a CSV file:

```R
data = read.csv('liver_cancer_metabolomic_data.csv',
                row.names = 1, stringsAsFactors = F, 
                na.strings = c("NA", "NaN", ""), check.names = F)
dim(data)
head(data)
is.numeric(data)
str(data)
```

- The `read.csv` function reads the data, with the first column as row names.
- `dim`, `head`, `is.numeric`, and `str` are used to inspect the data's dimensions, structure, and data types.

---

### **3. Data Preprocessing**
#### **3.1. Converting Data to Numeric**
The data is converted to a numeric matrix to ensure compatibility with downstream analysis:

```R
int.mat <- data
num.mat = matrix(as.numeric(gsub(",", "", as.matrix(int.mat))),
                 nrow = nrow(int.mat),
                 dimnames = dimnames(int.mat))
head(num.mat)
mode(num.mat)
dim(data)
```

- Commas are removed, and the data is converted to a numeric matrix.

#### **3.2. Removing Zero-Variance Features**
Columns with zero variance are removed to reduce noise:

```R
varCol = apply(num.mat, 2, var, na.rm = T)
constCol <- (varCol == 0 | is.na(varCol))
sum(constCol)
num.mat <- num.mat[, !constCol]
dim(num.mat)
```

- `apply` calculates the variance for each column.
- Columns with zero variance or `NA` variance are removed.

#### **3.3. Handling Missing Data**
The script calculates the percentage of missing values and visualizes the missingness rate:

```R
missingRatePerMetabolite = (apply(is.na(num.mat), 2, sum) / nrow(num.mat)) * 100
min(missingRatePerMetabolite)
max(missingRatePerMetabolite)

# Histogram for missingness rate
h = hist((apply(is.na(num.mat), 2, sum) / nrow(num.mat)) * 100, breaks = 10,
          main = "Histogram for Missingness",
          xlab = "Percentage of missingness")
text(h$mids, h$counts, labels = h$counts, adj = c(0.5, -0.5))
abline(v = 50, col = "red", lwd = 3, lty = 2)
```

- A histogram is created to visualize the distribution of missing values.
- Metabolites with more than 50% missing values are removed.

#### **3.4. Imputing Missing Values**
Missing values are imputed using the `impute.knn` function:

```R
num.mat.imputed = impute.knn(num.mat, k = 10)$data
class(num.mat.imputed)
head(num.mat.imputed)
dim(num.mat.imputed)
```

- The `impute.knn` function imputes missing values using the k-nearest neighbors algorithm.

#### **3.5. Log Transformation**
The data is log-transformed to stabilize variance:

```R
num.mat.imputed.logged <- log2(num.mat.imputed + 1)
head(num.mat.imputed.logged)
write.csv(num.mat.imputed.logged, 'Normalized_data.csv')
```

- Log transformation is applied to the imputed data, and the result is saved to a CSV file.

---

### **4. Data Visualization**
#### **4.1. Density Plots**
Density plots are used to compare the distribution of data before and after log transformation:

```R
par(mfrow = c(1, 2))
plot(density(num.mat.imputed[, 1]), main = 'Before log2')
plot(density(num.mat.imputed.logged[, 1]), main = 'After log2')
```

#### **4.2. Boxplots**
Boxplots are generated to visualize the distribution of metabolites:

```R
boxplot(num.mat.imputed[, 1:10], main = "Before log2", horizontal = T, 
        names = colnames(num.mat.imputed)[1:10], las = 2, col = "lightgreen")
boxplot(num.mat.imputed.logged[, 1:10], main = "After log2", horizontal = T, 
        names = colnames(num.mat.imputed.logged)[1:10], las = 2, col = "lightgreen")
```

#### **4.3. Scaling Data**
The data is scaled to have zero mean and unit variance:

```R
num.mat.imputed.logged.scaled = scale(num.mat.imputed.logged, center = TRUE, scale = TRUE)
head(num.mat.imputed.logged.scaled)
```

---

### **5. Principal Component Analysis (PCA)**
PCA is performed to reduce dimensionality and visualize patterns in the data:

```R
df_pca <- prcomp(Normalized_data)
df_out <- as.data.frame(df_pca$x)
ggplot(df_out, aes(x = PC1, y = PC2, color = samples_meta$Stage, shape = samples_meta$Gender)) +
  geom_point() + ggtitle("") + labs(color = '') +
  geom_point(size = 8, alpha = 0.5) +
  theme(...) + xlab(paste0("PC 1 (", round(df_pca$sdev[1], 1), "%)")) +
  ylab(paste0("PC 2 (", round(df_pca$sdev[2], 1), "%)"))
```

- The PCA plot is colored by cancer stage and shaped by gender.

---

### **6. Clustering Analysis**
#### **6.1. K-Means Clustering**
K-means clustering is performed for different numbers of clusters:

```R
kmeans2 <- kmeans(Normalized_data, centers = 2, nstart = 25)
kmeans3 <- kmeans(Normalized_data, centers = 3, nstart = 25)
kmeans4 <- kmeans(Normalized_data, centers = 4, nstart = 25)
kmeans5 <- kmeans(Normalized_data, centers = 5, nstart = 25)

# Visualizing clusters
plot1 <- fviz_cluster(kmeans2, geom = "point", data = Normalized_data) + ggtitle("k = 2")
plot2 <- fviz_cluster(kmeans3, geom = "point", data = Normalized_data) + ggtitle("k = 3")
plot3 <- fviz_cluster(kmeans4, geom = "point", data = Normalized_data) + ggtitle("k = 4")
plot4 <- fviz_cluster(kmeans5, geom = "point", data = Normalized_data) + ggtitle("k = 5")
grid.arrange(plot1, plot2, plot3, plot4, nrow = 2)
```

#### **6.2. Hierarchical Clustering**
Hierarchical clustering is performed using the Ward method:

```R
d <- dist(Normalized_data, method = "euclidean")
fit <- hclust(d, method = "ward")
plot(fit)
rect.hclust(fit, k = 5, border = "red")
```

---

### **7. Differential Expression Analysis**
Differential expression analysis is performed using the `limma` package:

```R
design <- model.matrix(~0 + factor(type) + samples_meta$Gender)
colnames(design) <- c(levels(factor(type)), 'Group')
contrast <- makeContrasts(Cancer - Healthy, levels = design)
fit <- lmFit(as.matrix(t(Normalized_data)), design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

DEGs <- topTable(fit2, adjust.method = 'fdr', number = Inf, p.value = 1, coef = 1)
```

- Differentially expressed metabolites (DEGs) are identified and filtered based on p-value and log fold change.

---

### **8. Pathway Analysis**
Pathway analysis is performed to identify enriched pathways:

```R
overlapping <- read.delim("C:/My pc/Egcombio/MODA/Metabolomics/ORA_results.tab")
overlapping_filtered <- overlapping

p6 <- ggplot() +
  geom_point(data = overlapping_filtered, 
             aes(x = size, y = -log(q.value), color = source, size = size)) +
  geom_label_repel(aes(x = size, y = -log(q.value), 
                   label = str_wrap(overlapping_filtered$pathway, width = 20))) +
  labs(x = "Size of genes pathway", y = "-Log(q-value)", color = "Pathway source") +
  ggtitle("Pathways analysis")
print(p6)
```

- The results are visualized using a scatter plot with labeled pathways.

---

### **9. Saving Results**
The final results, including DEGs and pathway analysis, are saved to CSV files:

```R
write.csv(DEGs, "DEGS.csv")
write.csv(overlapping, "liver_cancer_normal_LIMMA_adj_0.05_1.csv")
```

---
###Thanks!
