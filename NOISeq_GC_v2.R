library(readxl)
library(tidyverse)
library(NOISeq) # version 2.14.1 last updated feb 11 2016

# read in data
counts <- read_excel('/Users/gretacook/Desktop/Thesis2/c.vio-gene-counts.GC.xlsx')
  #4391 observations means 4391 mapped genes -(?) 

# filter 0-count genes - necessary, results change with filtering, 
  #filtering out genes with zero counts under both conditions, 
  #keeps genes that has at least one count for one condition
counts_fil <- counts %>%
  filter(rowSums(select(., -1)) > 0)

nrow(counts) #4391 genes 
nrow(counts_fil) #2755 genes with at least one count for at least one sample

# Convert to a matrix, excluding the first column (gene IDs)
count_matrix <- as.matrix(counts_fil[, -1])  # Exclude the first column
rownames(count_matrix) <- counts_fil[[1]]     # Set row names to gene IDs

# samples
meta_data <- data.frame(Sample = c("263-light", "260-dark"))
meta_data

# create noiseq object
mydata <- readData(data = count_matrix, factors = meta_data)
mydata

# quick dge
myresults <- noiseq(mydata, factor = "Sample", k = NULL, norm = "tmm", pnr = 0.2, 
                    nss = 5, v = 0.02, lc = 1, replicates = "no")

  #my data is filtered counts
  #factor is sample (light vs dark)
  #k is null - zero is set to midpoint between zero and next nonzero value 
  #norm is tmm - trimmed Mean of M values method of normalization
  #pnr is proportion of nonresponders - sets the proportion of genes that do not change between conditions
    #a value of 0.2 means 20% of the genes will not be differentially expressed
  #nss is number of simulated replicates
  #v is variance - sets the variability of simulated data. 0.02 means variablilty in simulated expression levels is low.  
  #lc is lower count threshold - discards data below set threshold

myresults.deg <- degenes(myresults, q = 0.8, M = NULL) 
  #759 differentially expressed genes with filtering out genes with zero counts for both conditions
  #761 differentially expressed features without removing any zeros
  

  # q is the quantile threshold for choosing DEGs. 
  #setting 0.8 considers the top 20% of genes based off their significance/expression differences
  # M specifies a log2 fold change threshold. 
  #setting to null sets no specific threshold and all genes meeting the quantile threshold will be included regardless of fold cahnge


myresults.deg


  # mean expression leves of the genes for each condition
  # M is the log 2-ratio of two conditions, 
  #positive value indicates higher expression in the dark sample
  #negative value indicates higher expression in the light sample, downregulation in dark sample 
  # D is the absolute difference between two conditions
  # prob is the probability that the gene is differentially expressed based on the analysis. 
  #values closer to 1 indicate higher likelihood of being differeentially expressed 
  # ranking the gene based on its expression change, significance, 
  #higher values may indicate higher significance   

head(myresults.deg)

            #X260.dark_mean X263.light_mean         M        D      prob   ranking
#CV_RS10275      28.247313       0.6361237  5.472663 27.61119 0.9938374  28.14832

  
 #highly expressed in the dark condition (28.25) compared to the light condition (0.64).
  #The large positive M value indicates strong upregulation in the dark.
  #The high probability value suggests that this gene is differentially expressed based on the analysis.
  #The high ranking value may help prioritize this gene in further analyses or studies.

##visualizations

#expression plot 
#plot of the average expression values of each condition, highlighting the features declared as differentially expressed 
DE.plot(myresults, q = 0.8, graphic = "expr", log.scale = TRUE)
#dots in red consider the threshold of 0.8
##[1] "734 differentially expressed features"


#MD plot 
#Plot the log fold change (M) values and the absolute value of the difference in expression between light and dark contioions
DE.plot(myresults, q = 0.8, graphic = "MD")

# Define breaks for consistent coloring
breaks<- seq(-5.5, 5.5, length.out = 101)  # 100 breaks from -3 to 3

# Extract all rows and the column "M" - log2FC values for heatmap
heatmap_data <- myresults.deg %>%
  select(M)

# Generate heatmap #cant read the genes and see the entire heatmap at the same time. 
pheatmap(heatmap_data,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = breaks,
         main = "Heatmap of Differentially Expressed Genes",
         cellwidth = 100,  # Adjust width of cells
         cellheight = 6,
         fontsize = 8,            # General font size
         fontsize_row = 6,       # Font size for row names
         fontsize_col = 6) # Adjust height of cells
#red means upregulated in dark condition, downregulated in light condition

#heatmap of top 10 DEG

# Filter for top 10 values in column M
top_ten_DEG <- myresults.deg[order(myresults.deg$M, decreasing = TRUE), ][1:10, ]
print(top_ten)

#extract M values for all genes

top_ten_DEG <- top_ten_DEG %>%
  select(M)


#generate heatmap for top 10 DEG 

pheatmap(top_ten_DEG,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = breaks,
         main = "Heatmap of Upregulated Genes in Dark Conditions",
         cellwidth = 100,  # Adjust width of cells
         cellheight = 6,
         fontsize = 8,            # General font size
         fontsize_row = 6,       # Font size for row names
         fontsize_col = 6) # Adjust height of cells

# Filter for bottom 10 values in column M
bottom_ten_DEG <- myresults.deg[order(myresults.deg$M), ][1:10, ]
print(bottom_ten_DEG)

bottom_ten_DEG <- bottom_ten_DEG %>%
  select(M)

pheatmap(bottom_ten_DEG,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = breaks,
         main = "Heatmap of Upregulated Genes in Light Conditions",
         cellwidth = 100,  # Adjust width of cells
         cellheight = 6,
         fontsize = 8,            # General font size
         fontsize_row = 6,       # Font size for row names
         fontsize_col = 6) # Adjust height of cells


#look at violacein genes 

#table of violacein genes
violacein.genes <- c("CV_RS16100", "CV_RS16105", "CV_RS16110", "CV_RS16115", "CV_RS16125", "CV_RS16130", "CV_RS16140", 
                     "CV_RS16150", "CV_RS16155", "CV_RS16160", "CV_RS16165", "CV_RS16170", "CV_RS16175", 
                     "sph", "vioB", "vioE")

vio.gene.data <- myresults.deg
vio.gene.data <- vio.gene.data[rownames(vio.gene.data) %in% violacein.genes, ]

vio.gene.data.h <- vio.gene.data %>%
  select(M)

#generate heatmap of violacein genes
pheatmap(vio.gene.data.h,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = breaks,
         main = "Heatmap of Differentially Expressed Violacein Genes",
         cellwidth = 100,  # Adjust width of cells
         cellheight = 6,
         fontsize = 8,            # General font size
         fontsize_row = 6,       # Font size for row names
         fontsize_col = 6) # Adjust height of cells


# Extract normalized counts - use for heatmap or use log2FC?
#normalized_counts <- exprs(mydata)



