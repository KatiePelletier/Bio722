##################################IAN DWORKIN##########################################
library(DESeq2)
library(tximport)
library(readr)
library(RColorBrewer)
library(gplots)

# This will differ for you!!!
setwd("~/University/Masters/Courses/Bio 722/Data/salmon_counts/")
# Setting it up for the import (this is not the import itself)
quant_files <- list.files(pattern = "quant.sf", recursive = TRUE) # Easier than what Ian did (for me anyways)

# Let's take a look at this variable we have created
quant_files

samples <- dir() # Easier than what Ian did
samples <- head(samples, length(quant_files))
names(quant_files) <- samples

# Prepare the translation from transcript to gene ID
tx2gene <- read.table("txp_to_gene.tsv", col.names=c("TXNAME", "GENEID")) # In a different location.  Will want to look at the Reference folder

# basic exploratin
head(tx2gene)
dim(tx2gene)

length(table(tx2gene$TXNAME))
length(table(tx2gene$GENEID))

# Importing! Type is the program used here
txi <- tximport(quant_files,
                type = "salmon",
                tx2gene = tx2gene)

#####################################
# Preparing the experimental design #
#####################################

# Tissue used
tissue <- c(rep("genital", 5),
            rep("wing", 6),
            rep("genital", 3),
            rep("wing", 6))

tissue <- as.factor(tissue)
length(tissue)


# What temperature?
temperature <- c(rep(17, 11), rep(24,9))
temperature <- as.factor(temperature)
length(temperature)

food <- c(rep("fed", 3),
          rep("starved", 2),
          rep("fed", 3),
          rep("starved", 3),
          rep("fed", 2),
          rep("starved", 1),
          rep("fed", 3),
          rep("starved", 3))
food <- as.factor(food)
length(food)

# What lane?
lane <- c(2,4,5,5,2,3,2,4,4,3,2,4,2,3,2,4,3,3,2,4)
lane <- factor(lane) # we will want to treat this as a factor
length(lane)

# Creates our RNA design plot
rna.design <- data.frame(sample=samples,
                         file=quant_files,
                         tissue=tissue,
                         food=food,
                         temperature=temperature,
                         lane = lane)

rna.design

# and we can start with a simple model (back to this later)
load.model <- formula(~ tissue)

# Creating the model.  Theses objects are in S4
all.data <- DESeqDataSetFromTximport(txi, 
                                     rna.design,
                                     design=load.model)

# Controling for the lanes, do they have an important effect on the results?
load.model <- formula(~ lane)

test_lane_effects <- DESeqDataSetFromTximport(txi, rna.design, design=load.model)

test_lane_effects2 <- DESeq(test_lane_effects)

# Getting the results
test_lane_effects2_results <- results(test_lane_effects2, alpha = 0.05) # alpha = 0.05 is the FWER (Family Wise Error Rate)

summary(test_lane_effects2_results)
# 2 genes which may show  evidence of lane effects, but this is a bit incomplete for the full data set.
# Low counts are genes that were thrown out prior to the anlaysis.  Not entirely clear on the definition of low here....

head(test_lane_effects2_results)

# let's re-order the data to look at genes.
test_lane_effects2_results <- test_lane_effects2_results[order(test_lane_effects2_results$padj),]


load.model <- formula(~ lane + tissue) #DESeq2 will assume the last variable is what we're interested in
test_tissue_effects <- DESeqDataSetFromTximport(txi, rna.design, design=load.model)

test_tissue_effects2 <- DESeq(test_tissue_effects)
#####################
# The real analysis #
#####################
load.model <- formula(~ lane + tissue) #DESeq2 will assume the last variable is what we're interested in
test_tissue_effects <- DESeqDataSetFromTximport(txi, rna.design, design=load.model)

test_tissue_effects2 <- DESeq(test_tissue_effects)
tissue_results <- results(test_tissue_effects2, contrast = c("tissue", "genital", "wing"), alpha = 0.05)
##################################IAN DWORKIN##########################################
##################################GEORGE LONG##########################################
# Libraries
library(topGO)
library(dplyr)
# Let's do some cleaning
rm(list = setdiff(ls(), "tissue_results"))

# Let's get our Data ready
gene_associations <- read.delim("C:/Users/getsl/Desktop/gene_association.fb", comment.char = "!", header = FALSE, as.is = TRUE) # Don't want factors
colnames(gene_associations) <- c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID",
                                 "DB:Reference", "Evidence", "With_From", "Aspect", "DB_Object_Name",
                                 "DB_Object_Synonym", "DB_Object_Type", "Taxon", "Date", "Assigned_by")
gene_associations <- gene_associations[,c(2,3,5,7,9,14)] # Trimming the dataframe so that it's only what we're interested in



gene_GO <- lapply(unique(gene_associations$DB_Object_ID), function(x){tmp <- gene_associations %>% filter(DB_Object_ID == x)
  return(tmp$GO_ID)}) # Could easily convert so that it only looks at one of the molecular functions
names(gene_GO) <- unique(gene_associations$DB_Object_ID)


#################Only to save your work##############################
# Making a mapping file
gene_GO <- sapply(gene_GO, function(x){paste(x, collapse = ", ")}) # Collapsing the GOIDS for the mapping file
gene_GO <- data.frame("GOID" = gene_GO) # Making a dataframe
write.table(gene_GO, file = "C:/Users/getsl/Desktop/fly_to_GO.delim", sep = "\t", quote = FALSE, col.names = FALSE) # Writing the dataframe in the correct format.  If you want to save te results

#####################################################################

# Gene of interest
genes_of_interest <- tissue_results@rownames[tissue_results@listData$padj < 0.05] # What genes had a padj < 0.05
genes_of_interest <- genes_of_interest[complete.cases(genes_of_interest)] # Removing the NAs

geneList <- factor(as.integer(names(gene_GO) %in% genes_of_interest)) # What genes are of interest?
names(geneList) <- names(gene_GO) # Adding the names

GOdata <- new("topGOdata",
              description = "BIO722 Tutorial on Gene Enrichment",
              ontology = "MF",
              allGenes = geneList,
              nodeSize = 10,
              annotationFun = annFUN.gene2GO,
              gene2GO = gene_GO) # Creating the topGO dataset

##################################GEORGE LONG##########################################
