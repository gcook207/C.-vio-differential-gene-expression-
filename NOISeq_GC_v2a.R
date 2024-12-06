library(readxl)
library(tidyverse)
library(dplyr)
library(pheatmap)
library(NOISeq) # version 2.14.1 last updated feb 11 2016

# read in data
counts <- read_excel(file.choose())
 counts <- read_excel('/Users/gretacook/Desktop/Thesis2/c.vio-gene-counts.GC.xlsx')
  #4391 observations means 4391 identified genes 

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

myresults.deg <- degenes(myresults, q = 0.9, M = NULL) 
  #368 differentially expressed genes with filtering out genes with zero counts for both conditions

myresults.deg.up <- degenes(myresults, q = 0.9, M = "up") 
  #218 differentially expressed features (up in dark condition)

myresults.deg.down <- degenes(myresults, q = 0.9, M = "down") 
  #150 differentially expressed features (down in light condition)

  ## A feature is considered to be differentially expressed if its corresponding M and D values are likely to be higher than in noise
  #Noise distribution is obtained by comparing all pairs of replicates within the same condition. 
  #The corresponding M and D values are pooled together to generate the distribution
  #by comparing the (M, D) values of a given feature against the noise distribution, NOISeq obtains the
    #“probability of differential expression” for this feature
  #simulates technical replicates in order to estimate the differential expression probability
  #simulates technical replicates from a multinomial distribution

  # q is the quantile threshold for choosing DEGs. 
  #setting 0.8 considers the top 20% of genes based off their significance/expression differences 
  ###when no replicates are available, then it is preferable to use a higher threshold such as q = 0.9
  # M specifies a log2 fold change threshold. 
  #setting to null sets no specific threshold and all genes meeting the quantile threshold will be included regardless of fold cahnge


myresults.deg


  # mean expression leves of the genes for each condition
  # M is the log 2-ratio of two conditions, 
  #positive value indicates higher expression in the dark sample
  #negative value indicates higher expression in the light sample, downregulation in dark sample 
  # D is the absolute difference between two conditions
  # prob is the probability that the gene is differentially expressed based on the analysis. 
      #(probabilities are not equivalent to p-values), the probability of differential expression is not equivalent to
        #1 − pvalue. 
      # The higher the probability, the more likely that the difference in expression is due to the 
      #change in the experimental condition and not to chance
  #values closer to 1 indicate higher likelihood of being differeentially expressed 
  # ranking the gene based on its expression change, significance, 
  #higher values may indicate higher significance   
  #when no replicates are available, then it is preferable to use a higher threshold such as q = 0.9

head(myresults.deg)

            #X260.dark_mean X263.light_mean         M        D      prob   ranking
#CV_RS10275      28.247313       0.6361237  5.472663 27.61119 0.9938374  28.14832

  
 #highly expressed in the dark condition (28.25) compared to the light condition (0.64).
  #The large positive M value indicates strong upregulation in the dark.
  #The high probability value suggests that this gene is differentially expressed based on the analysis.
  #The high ranking value may help prioritize this gene in further analyses or studies.

##table includidng all relevant data to include in appendix of thesis document
all.data <- myresults.deg

all.data <- merge(all.data, count_matrix, by = "row.names", all.x = TRUE)

all.data <-  all.data %>%
  rename(
    "normalized light counts" = X263.light_mean,
    "normalized dark counts" = X260.dark_mean,
    "log2FC" = M, 
    "difference" = D,
    "probability" = prob,
    "ranking" = ranking, 
    "raw light counts" = `263-light`,
    "raw dark counts" = `260-dark`,
    "gene id" = Row.names
  )

all.data <- all.data[c("gene id", "raw light counts", "raw dark counts", "normalized light counts", "normalized dark counts", "log2FC", "difference", "probability", "ranking" )]

library(openxlsx)
write.xlsx(all.data, "all.data.xlsx")

library(officer)

doc <- read_docx()

# Add the table to the Word document
doc <- doc %>%
  body_add_table(value = all.data, style = "table_template")

# Save the document
print(doc, target = "my_report.docx")

##visualizations

#expression plot 
#trouble with x and y labels

#plot of the average expression values of each condition, highlighting the features declared as differentially expressed 

colnames(myresults@results[[1]])[1] <- "Normalized Light Counts"

colnames(myresults@results[[1]])[2] <- "Normalized Dark Counts"

DE.plot(myresults, q = 0.9, graphic = "expr", log.scale = TRUE, main = expression("Differentially Expressed Genes in " ~ italic("C. violaceum") ~ " between Light and Dark Conditions"))
#dots in red consider the threshold of 0.9
##331 differentially expressed features

library(ggplot2)
filtered_results <- myresults[myresults$prob <= 0.9, ]

# edits by TM
#heatmap of top 25 DEG

# Filter for top & down-regulated genes
top_up_DEG <- myresults.deg[order(myresults.deg$M, decreasing = TRUE), ][1:25,]
print(top_up_DEG)
top_down_DEG <- myresults.deg[order(myresults.deg$M, decreasing = FALSE), ][1:25,]
print(top_down_DEG)

# combine dataframes
top_genes <- bind_rows(top_up_DEG, top_down_DEG)
top_genes


# extract counts for all top genes
top_genes_counts <- top_genes %>% select(X260.dark_mean, X263.light_mean)
top_genes_counts

colnames(top_genes_counts) <- c("Dark", "Light")

#generate heatmap for top 25 DEG 

pheatmap(top_genes_counts,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Top 50 Differentially Expressed Genes",
         cellwidth = 100,  # Adjust width of cells
         cellheight = 6,
         fontsize = 8,            # General font size
         fontsize_row = 6,       # Font size for row names
         fontsize_col = 6, 
         ) # Adjust height of cells


#analysis of special metabolite gene clusters

#violacein genes 

#table of violacein genes
violacein.genes <- c("CV_RS16100", "CV_RS16105", "CV_RS16110", "CV_RS16115", "CV_RS16125", "CV_RS16130", "CV_RS16140", 
                     "CV_RS16150", "CV_RS16155", "CV_RS16160", "CV_RS16165", "CV_RS16170", "CV_RS16175", 
                     "sph", "vioB", "vioE")

vio.gene.data <- myresults.deg
vio.gene.data <- vio.gene.data[rownames(vio.gene.data) %in% violacein.genes, ]

vio.gene.data.h <- vio.gene.data %>% select(X260.dark_mean, X263.light_mean)

#doesnt change the sample column labels on heatmap
#vio.gene.data.h <-  vio.gene.data.h %>%
#  rename(
#    "normalized light counts" = X263.light_mean,
#    "normalized dark counts" = X260.dark_mean
#  )


#generate heatmap of violacein genes
pheatmap(vio.gene.data.h,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Heatmap of Differentially Expressed Violacein Genes",
         cellwidth = 100,  # Adjust width of cells
         cellheight = 6,
         fontsize = 8,            # General font size
         fontsize_row = 6,       # Font size for row names
         fontsize_col = 6,
         ) # Adjust height of cells

#chromobactin genes

chromobactin.genes <- c("CV_RS07170","CV_RS07180", "CV_RS07190","CV_RS07195","CV_RS07200","CV_RS07205","CV_RS07210", "CV_RS07215", "CV_RS07220",
                      "CV_RS07225", "CV_RS07235", "CV_RS07240", "CV_RS07245", "CV_RS07250", "CV_RS07255",
                     "CV_RS07260", "CV_RS07265", "CV_RS07270", "CV_RS07275", "CV_RS07280", "CV_RS07285", "CV_RS07290", "CV_RS07295", "CV_RS07320", 
                     "CV_RS23430", "CV_RS23440", "CV_RS24215", "aruF", "astA", "astB", "astD", "dhbA" )

chromobactin.gene.data <- myresults.deg
chromobactin.gene.data <- chromobactin.gene.data[rownames(chromobactin.gene.data) %in% chromobactin.genes, ]

chromobactin.gene.data.h <- chromobactin.gene.data %>% select(X260.dark_mean, X263.light_mean)

#generate heatmap of Chromobactin genes
pheatmap(chromobactin.gene.data.h,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Heatmap of Differentially Expressed Chromobactin Genes",
         cellwidth = 100,  # Adjust width of cells
         cellheight = 6,
         fontsize = 8,            # General font size
         fontsize_row = 6,       # Font size for row names
         fontsize_col = 6,
) # Adjust height of cells

#viobactin genes
viobactin.genes <- c("CV_RS10870","CV_RS10875", "CV_RS10880","CV_RS10890","CV_RS10895","CV_RS10900","CV_RS10905", "CV_RS10915", "CV_RS10920",
                        "CV_RS10925", "CV_RS10930", "CV_RS10935", "CV_RS10945", "CV_RS10950", "CV_RS10960",
                        "CV_RS10965", "CV_RS10970", "CV_RS10990", "CV_RS11005", "CV_RS11010", "CV_RS11015", "CV_RS11020", "CV_RS11025", "CV_RS11030", 
                        "CV_RS11035", "CV_RS11040", "CV_RS11045", "CV_RS11050", "CV_RS11055", "CV_RS24250", "entS", "fepB", "fepD", "fepG", "fes" )

viobactin.gene.data <- myresults.deg
viobactin.gene.data <- viobactin.gene.data[rownames(viobactin.gene.data) %in% viobactin.genes, ]

viobactin.gene.data.h <- viobactin.gene.data %>% select(X260.dark_mean, X263.light_mean)

#generate heatmap of Chromobactin genes
pheatmap(viobactin.gene.data.h,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Heatmap of Differentially Expressed Viobactin Genes",
         cellwidth = 100,  # Adjust width of cells
         cellheight = 6,
         fontsize = 8,            # General font size
         fontsize_row = 6,       # Font size for row names
         fontsize_col = 6,
) # Adjust height of cells



#hydrogen-cyanide genes
hydrogen.cyanide.genes <- c("CV_RS08230","CV_RS08235", "CV_RS08240","CV_RS08245","CV_RS08255","CV_RS08265","CV_RS08275", "CV_RS08280", "CV_RS08285",
                        "CV_RS08290", "CV_RS22120", "hcnA", "hcnC")

hydrogen.cyanide.gene.data <- myresults.deg
hydrogen.cyanide.gene.data <- hydrogen.cyanide.gene.data[rownames(hydrogen.cyanide.gene.data) %in% hydrogen.cyanide.genes, ]

hydrogen.cyanide.gene.data.h <- hydrogen.cyanide.gene.data %>% select(X260.dark_mean, X263.light_mean)

#generate heatmap of Chromobactin genes
pheatmap(hydrogen.cyanide.gene.data.h,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Heatmap of Differentially Expressed Hydrogen-cyanide Genes",
         cellwidth = 100,  # Adjust width of cells
         cellheight = 6,
         fontsize = 8,            # General font size
         fontsize_row = 6,       # Font size for row names
         fontsize_col = 6,
) # Adjust height of cells


#normalized_counts <- exprs(mydata)



