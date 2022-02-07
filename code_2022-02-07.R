# ----------------------------------------------------------------------------------
# Set functions
# ----------------------------------------------------------------------------------

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95) {
    library(doBy)
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))
    datac <- summaryBy(formula, data=data, FUN=c(length2,mean,sd), na.rm=na.rm)

    names(datac)[ names(datac) == paste(measurevar, ".mean",    sep="") ] <- measurevar
    names(datac)[ names(datac) == paste(measurevar, ".sd",      sep="") ] <- "sd"
    names(datac)[ names(datac) == paste(measurevar, ".length2", sep="") ] <- "N"
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
}

# select significant ASVs between groups
gm_mean = function(x, na.rm=TRUE){
	exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
select_sigtab <- function(physeq, group) {
	diagdds <- phyloseq_to_deseq2(physeq, ~ group)
	geoMeans <- apply(counts(diagdds), 1, gm_mean)
	diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans)
	diagdds <- DESeq(diagdds, test="Wald", fitType="parametric")
	res <- results(diagdds, cooksCutoff = FALSE)
	alpha <- 0.05
	sigtab <- res[which(res$pvalue < alpha), ]		
	sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(physeq_sf_1)[rownames(sigtab), ], "matrix"))
	otu_ID <- rownames(otu_table(physeq)[which(res$pvalue < alpha)])
	return(sigtab)
}

# select gene count data from the groups to be compared
select_countData <- function(countData, group_name, comparison){
	library(dplyr)
	Gene1 <- colnames(countData)[1]
	m <- as.data.frame(t(as.matrix(countData[,-1])))
	m$group <- group_name
	m <- filter(m, group==comparison[1] | group==comparison[2])
	m$group <- NULL
	m <- as.data.frame(t(m))
	m <- cbind(countData[1], m)
	colnames(m)[1] <- Gene1
	rownames(m) <- c(1:nrow(m))
	return(m)
}

# select metadata from the groups to be compared
select_metaData <- function(metaData, group_name, comparison){
	m <- metaData
	m$group <- group_name
	m <- filter(metaData, group==comparison[1] | group==comparison[2])
	return(m)
}

# construct DESeqDataSet Object
dds_object <- function(countData, metaData){
	dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~group, tidy = TRUE)
	dds$group <- factor(dds$group, levels = comparison)
	keep <- rowSums(counts(dds)) >= 10	# Pre-filtering
	dds <- dds[keep,]
	dds <- DESeq(dds)	
	return(dds)
}

# select significant genes between groups
select_siggene <- function(dds){
	res <- results(dds) 
	alpha <- 0.05
	siggene <- res[which(res$padj < alpha), ]
	return(siggene)
}

# bar plot for CCA
barplot_cca <- function(u1,is_desc,title2,xlab="",ylab="",topK=50,fillcol="Steelblue") {
	u1$order <- c(1:dim(u1)[1])
	u1 <- u1[order(cca_coeff,decreasing = is_desc),][1:topK]
	u1 <- u1[!is.na(name),]
	u1 <- u1[order(order,decreasing = T),]
	p <- ggplot(u1,aes(x=reorder(name, -order),y=cca_coeff)) +
		geom_bar(stat = "identity", fill=fillcol) +
		ggtitle(title2) + labs(x=ylab, y=xlab) +
		coord_flip()  
	p
}

# get colors for heatmap
get_red_white_blue_for_heatmap <- function(a,paletteLength=50,debug2=0) {
	if (debug2==1) {browser()}
	def_min_neg <- -1e-3
	def_max_pos <- 1e-3
	myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
	myBreaks <- c(seq(min(min(a[!is.na(a)]),def_min_neg), 0, length.out=ceiling(paletteLength/2) + 1),
			seq(max(a[!is.na(a)])/paletteLength, max(max(a[!is.na(a)]),def_max_pos), length.out=floor(paletteLength/2)))
	dupindc <- duplicated(myBreaks)
	myBreaks <- myBreaks[!dupindc]
	return(list(colors=myColor,breaks=myBreaks))
}

# scaling for heatmap
scale_and_clip <- function(expriv,scale.center=TRUE) {
	expriv.scaled = t(apply(expriv, 1, function(x) {
		q10 = quantile(x, 0.1)
		q90 = quantile(x, 0.9)
		x[x < q10] = q10
		x[x > q90] = q90
		scale(x,center = scale.center)
	}))
	dim(expriv.scaled)
	colnames(expriv.scaled) = colnames(expriv)
	j <- rowSums(is.nan(expriv.scaled))==0
	expriv.scaled[j,]
}

library(phyloseq)
library(ggplot2)
library(vegan)
library(DESeq2)
library(pheatmap)

################################################################################################
## Gastric microbiome (16S rRNA gene sequencing)
################################################################################################

# ----------------------------------------------------------------------------------
# Reading 16S samples 
# (* Run this code if you wish to analyze the BIOPSY cohort)
# ----------------------------------------------------------------------------------
m <- read.csv("./biopsy, 16S/ASV_table_biopsy.txt", header=TRUE, sep='\t', stringsAsFactor=FALSE)
ASV_matrix <- as.matrix(m[,-1])
rownames(ASV_matrix) <- m[,1]

m <- read.csv("./biopsy, 16S/Taxonomy_table_biopsy.txt", header=TRUE, sep='\t', stringsAsFactor=FALSE)
TAX_matrix <- as.matrix(m[,-1])
rownames(TAX_matrix) <- m[,1]

metaData <- read.csv("metadata_biopsy.txt", header=TRUE, sep='\t', na.strings=9999, stringsAsFactor=TRUE)
metaData[ ,"group"] <- factor(metaData[ ,"group"], level=c("Healthy", "HP_eradication", "HP_gastritis", "Gastric_cancer"))

treefile <- "./biopsy, 16S/tree_biopsy.nwk"

# ----------------------------------------------------------------------------------
# Reading 16S samples from the surgery cohort
# (* Run this code if you wish to analyze the SURGERY cohort)
# ----------------------------------------------------------------------------------
m <- read.csv("./surgery, 16S/ASV_table_surgery.txt", header=TRUE, sep='\t', stringsAsFactor=FALSE)
ASV_matrix <- as.matrix(m[,-1])
rownames(ASV_matrix) <- m[,1]

m <- read.csv("./surgery, 16S/Taxonomy_table_surgery.txt", header=TRUE, sep='\t', stringsAsFactor=FALSE)
TAX_matrix <- as.matrix(m[,-1])
rownames(TAX_matrix) <- m[,1]

metaData <- read.csv("metadata_surgery.txt", header=TRUE, sep='\t', na.strings=9999, stringsAsFactor=TRUE)
metaData[ ,"group"] <- factor(metaData[ ,"group"], level=c("Severe_gastritis", "Gastric_cancer"))

treefile <- "./surgery, 16S/tree_surgery.nwk"

# ----------------------------------------------------------------------------------
# Making a phyloseq object
# ----------------------------------------------------------------------------------
ASV <- otu_table(ASV_matrix, taxa_are_rows = TRUE)
TAX <- tax_table(TAX_matrix)
Sample <- sample_data(metaData); rownames(Sample) <- metaData$sample_ID
Tree <- read_tree(treefile)

physeq_data <- phyloseq(ASV, TAX, Sample, Tree)
physeq_prune <- prune_taxa(taxa_sums(physeq_data) > 0, physeq_data)	# prefiltering for diversity analysis
physeq <- filter_taxa(physeq_prune, function(x) sum(x > 10) > (0.05*length(x)), TRUE)	# filtering for abundance analysis
physeq

# Standardize abundances to the median sequencing depth
total <- median(sample_sums(physeq))
standf <- function(x, t=total) round(t * (x / sum(x)))
physeq_sf <- transform_sample_counts(physeq, standf)
physeq_sfr <- transform_sample_counts(physeq_sf, function(x) x  / sum(x) )

# ----------------------------------------------------------------------------------
# Alpha diversity
# ----------------------------------------------------------------------------------
physeq_prune <- prune_taxa(taxa_sums(physeq_data) > 0, physeq_data)
physeq_prune <- transform_sample_counts(physeq_prune, function(x) round(x))
p <- plot_richness(physeq_prune, x="group", color="group", measures=c("Observed", "Chao1", "Shannon", "InvSimpson")) +
		geom_boxplot(size=0.1, outlier.shape=NA) + 
		geom_jitter(size=0.2, alpha=1) + 
		labs(x = "Group", y = "Measures") + 
		theme(aspect.ratio=2, axis.title=element_text(size=8), axis.text=element_text(size=5), legend.title=element_blank(), legend.text=element_text(size=6), axis.ticks=element_line(size=0.1))
p$layers <- p$layers[c(-1, -2)]
p

rich <- estimate_richness(physeq_prune, measures=c("Observed", "Chao1", "Shannon", "InvSimpson"))
head(rich)

# Mann-Whitney U test (Wilcoxonn rank-sum test)
pairwise.wilcox.test(rich$Observed, sample_data(physeq_prune)$group)
pairwise.wilcox.test(rich$Chao1, sample_data(physeq_prune)$group)
pairwise.wilcox.test(rich$Shannon, sample_data(physeq_prune)$group)
pairwise.wilcox.test(rich$InvSimpson, sample_data(physeq_prune)$group)

# ----------------------------------------------------------------------------------
# Beta diversity: MDS (*PCoA*) on Bray-Curtis Dissimilarity (by group)
# 	(for both the biopsy and surgery cohorts)
# ----------------------------------------------------------------------------------
dist = phyloseq::distance(physeq_prune, method="bray")
ordu <- ordinate(physeq_prune, "PCoA", dist)
p <- plot_ordination(physeq_prune, ordu, color="group") +
	geom_point(size=3, alpha=0.7) + theme(aspect.ratio=1, axis.text=element_text(size=8), legend.title=element_blank(), legend.text=element_text(size=6), axis.ticks=element_line(size=0.1))
p$layers <- p$layers[-1]
p

adonis(dist ~ sample_data(physeq_prune)$group)	# permutational ANOVA (PERMANOVA) analysis
ano <- anosim(t(otu_table(physeq_prune)), sample_data(physeq_prune)$group, distance = "bray", permutations = 9999)	# ANOSIM analysis
summary(ano)

# ----------------------------------------------------------------------------------
# Beta diversity: MDS (*PCoA*) on Bray-Curtis Dissimilarity (by gastric cancer histology)
#	(* Only for the surgery cohort)
# ----------------------------------------------------------------------------------
dist = phyloseq::distance(physeq_prune, method="bray")
ordu <- ordinate(physeq_prune, "PCoA", dist)
p <- plot_ordination(physeq_prune, ordu, color="tumor_histology") +
	stat_ellipse(type="norm", level=0.95)+
	geom_point(size=3, alpha=0.7) + theme(aspect.ratio=1, axis.text=element_text(size=8), legend.title=element_blank(), legend.text=element_text(size=6), axis.ticks=element_line(size=0.1)) + theme_bw()
p$layers <- p$layers[-1]
p

adonis(dist ~ sample_data(physeq_min)$group)	# permutational ANOVA (PERMANOVA) analysis
ano <- anosim(t(otu_table(physeq_min)), sample_data(physeq_min)$group, distance = "bray", permutations = 9999)	# ANOSIM analysis
summary(ano)

# ----------------------------------------------------------------------------------
# Relative abundance of bacterial taxa
# ----------------------------------------------------------------------------------
# Select Palette
library(scales)
library(RColorBrewer)
library(doBy)
taxa_level <- 'Phylum'
psdat_gen_r <- tax_glom(physeq_sfr, taxrank = "Phylum")
ps_melt_r <- psmelt(psdat_gen_r)
summary_r_sample <- summarySE(ps_melt_r, measurevar="Abundance", groupvars=c("sample_ID", "group", taxa_level))
taxa <- levels(as.factor(summary_r_sample[,3]))
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
col_Phylum <- getPalette(length(taxa))
col_group <- hue_pal()(4)

# plot
psdat_gen_r <- tax_glom(physeq_sfr, taxrank = "Phylum")
ps_melt_r <- psmelt(psdat_gen_r)
summary_r <- summarySE(ps_melt_r, measurevar="Abundance", groupvars=c("group", "Phylum"))
ggplot(data=summary_r, aes(x=Phylum, y=Abundance * 100, fill=Phylum))	+
	geom_col(data=summary_r) + facet_grid("group~.", scales="free_x") +
	scale_y_continuous(limits=c(0,100), breaks=seq(0,100,25) ) +
	scale_fill_manual(values = col_Phylum) + scale_color_manual(values = col_Phylum) +
      labs(x = "", y = "Relative abundance (%)") +
	theme(aspect.ratio=0.1, axis.title=element_text(size=14), axis.text.x = element_text(size=12, angle = 45, hjust = 1), legend.title=element_blank(), legend.text=element_text(size=12), axis.ticks=element_line(size=0.1))


# ----------------------------------------------------------------------------------
# Differential abundance analysis
# ----------------------------------------------------------------------------------
# Identifying significant ASVs between the groups (* Only for the biopsy cohort)
# 	Compared groups:
#		(1) HP_eradication vs. Healthy
#		(2) HP_gastritis vs. Healthy
#		(3) Gastric_cancer vs. Healthy
# ----------------------------------------------------------------------------------
ASV_ID_sig <- NULL

# (1) HP_eradication vs. Healthy)
comparison <- c("Healthy", "HP_eradication")
physeq_sf_1 <- merge_phyloseq(subset_samples(physeq_sf, group == comparison[1]), subset_samples(physeq_sf, group == comparison[2]))
sigtab <- select_sigtab(physeq_sf_1, group)
ASV_ID_sig <- unique(c(ASV_ID_sig, rownames(sigtab)))
# ploting
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, FALSE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))
ggplot(sigtab, aes(x=log2FoldChange, y=Family, color=Phylum)) + geom_jitter(size=2, height=0.2, alpha=0.7) + 
  scale_color_manual(values = col_Phylum) +
  labs(x = "Bacterial taxa", y = "Log2 fold change") +
  theme(aspect.ratio=2, axis.title=element_blank(), axis.text.y = element_text (size=7), axis.text.x = element_text(size=7), legend.title=element_text(size=10), legend.text=element_text(size=7), axis.ticks=element_line(size=0.1))

# (2) HP_gastritis vs. Healthy
comparison <- c("Healthy", "HP_gastritis")
physeq_sf_1 <- merge_phyloseq(subset_samples(physeq_sf, group == comparison[1]), subset_samples(physeq_sf, group == comparison[2]))
sigtab <- select_sigtab(physeq_sf_1, group)
ASV_ID_sig <- unique(c(ASV_ID_sig, rownames(sigtab)))
# ploting
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, FALSE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))
ggplot(sigtab, aes(x=log2FoldChange, y=Family, color=Phylum)) + geom_jitter(size=2, height=0.2, alpha=0.7) + 
  scale_color_manual(values = col_Phylum) +
  labs(x = "Bacterial taxa", y = "Log2 fold change") +
  theme(aspect.ratio=2, axis.title=element_blank(), axis.text.y = element_text (size=7), axis.text.x = element_text(size=7), legend.title=element_text(size=10), legend.text=element_text(size=7), axis.ticks=element_line(size=0.1))

# (3) Gastric_cancer vs. Healthy
comparison <- c("Healthy", "Gastric_cancer")
physeq_sf_1 <- merge_phyloseq(subset_samples(physeq_sf, group == comparison[1]), subset_samples(physeq_sf, group == comparison[2]))
sigtab <- select_sigtab(physeq_sf_1, group)
ASV_ID_sig <- unique(c(ASV_ID_sig, rownames(sigtab)))
# ploting
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, FALSE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))
ggplot(sigtab, aes(x=log2FoldChange, y=Family, color=Phylum)) + geom_jitter(size=2, height=0.2, alpha=0.7) + 
  scale_color_manual(values = col_Phylum) +
  labs(x = "Bacterial taxa", y = "Log2 fold change") +
  theme(aspect.ratio=2, axis.title=element_blank(), axis.text.y = element_text (size=7), axis.text.x = element_text(size=7), legend.title=element_text(size=10), legend.text=element_text(size=7), axis.ticks=element_line(size=0.1))

ASV_ID_sig	# significant ASVs between any diseased group (HP_eradication, HP_gastritis, or Gastric_cancer) and Healthy group

# ----------------------------------------------------------------------------------
# Identifying significant ASVs between the groups (* Only for the surgery cohort)
#   (Do not run this code if you are analyzing the biopsy cohort)
# 	Compared group:
#		(1) Gastric_cancer vs. Severe_gastritis
# ----------------------------------------------------------------------------------
ASV_ID_sig <- NULL

# (1) Gastric_cancer vs. Severe_gastritis
comparison <- c("Severe_gastritis", "Gastric_cancer")
physeq_sf_1 <- merge_phyloseq(subset_samples(physeq_sf, group == comparison[1]), subset_samples(physeq_sf, group == comparison[2]))
sigtab <- select_sigtab(physeq_sf_1, group)
ASV_ID_sig <- unique(c(ASV_ID_sig, rownames(sigtab)))
# ploting
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, FALSE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))
ggplot(sigtab, aes(x=log2FoldChange, y=Family, color=Phylum)) + geom_jitter(size=2, height=0.2, alpha=0.7) + 
  scale_color_manual(values = col_Phylum) +
  labs(x = "Bacterial taxa", y = "Log2 fold change") +
  theme(aspect.ratio=2, axis.title=element_blank(), axis.text.y = element_text (size=7), axis.text.x = element_text(size=7), legend.title=element_text(size=10), legend.text=element_text(size=7), axis.ticks=element_line(size=0.1))

ASV_ID_sig	# significant ASVs between Gastric_cancer and Severe_gastritis groups

# ----------------------------------------------------------------------------------
# Make a "significant ASVs" table
# ----------------------------------------------------------------------------------
ASV_ID_sig <- data.frame(ASV = ASV_ID_sig)
ASV_tab <- as.data.frame(otu_table(physeq_sf))
ASV_tab$ASV <- rownames(ASV_tab)
ASV_tab <- merge(ASV_tab, ASV_ID_sig, by = "ASV", sort=FALSE)
tax_tab <- as.data.frame(tax_table(physeq_sf))
tax_tab$ASV <- rownames(tax_tab)
tax_tab2 <- data.frame(tax_tab$ASV, tax_tab$Family); colnames(tax_tab2) <- c("ASV", "Taxonomy")
ASV_tab2 <- merge(tax_tab2, ASV_tab, by = "ASV", sort=FALSE)

# ----------------------------------------------------------------------------------
# Replacing ASV unique identifier with ID number
# ----------------------------------------------------------------------------------
ASV_ID <- read.csv("unique_identifier.txt", header=TRUE, sep='\t', stringsAsFactor=FALSE)
colnames(ASV_ID)[2] <- "ASV"; ASV_ID <- ASV_ID[,-3]
ASV_tab2 <- merge(ASV_tab2, ASV_ID, by = "ASV", sort=FALSE)
rownames(ASV_tab2) <- paste(ASV_tab2$ASV_ID, "; ", ASV_tab2$Taxonomy, sep='')
ASV_tab2 <- ASV_tab2[,c(-1,-2)]
ASV_tab2 <- ASV_tab2[,-dim(ASV_tab2)[2]]
ASV_tab2 <- t(ASV_tab2)

################################################################################################
## Human gene expression (RNA-seq)
################################################################################################
library(dplyr)
library(pcaExplorer)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(gage)
library(gageData)

# ----------------------------------------------------------------------------------
# Reading RNA-seq samples (* Only for the biopsy cohort)
# ----------------------------------------------------------------------------------
m <- read.csv("./biopsy, RNA-seq/gene_table_biopsy.txt", header=TRUE, sep='\t', stringsAsFactor=FALSE)
countData_original <- m
metaData_original <- read.csv("metadata_biopsy.txt", header=TRUE, sep='\t', na.strings=9999, stringsAsFactor=TRUE)
metaData_original[ ,"group"] <- factor(metaData_original[ ,"group"], level=c("Healthy", "HP_eradication", "HP_gastritis", "Gastric_cancer"))
comparison <- c("Healthy", "HP_eradication", "HP_gastritis", "Gastric_cancer")

# ----------------------------------------------------------------------------------
# Reading RNA-seq samples (* Only for the the surgery cohort)
#   (Do not run this code if you wish to analyze the biopsy cohort)
# ----------------------------------------------------------------------------------
m <- read.csv("./surgery, RNA-seq/gene_table_surgery.txt", header=TRUE, sep='\t', stringsAsFactor=FALSE)
countData_original <- m
metaData_original <- read.csv("metadata_surgery.txt", header=TRUE, sep='\t', na.strings=9999, stringsAsFactor=TRUE)
metaData_original[ ,"group"] <- factor(metaData[ ,"group"], level=c("Severe_gastritis", "Gastric_cancer"))
comparison <- c("Severe_gastritis", "Gastric_cancer")

# ----------------------------------------------------------------------------------
# PCA plot by group
#	(for both the biopsy and surgery cohorts)
# ----------------------------------------------------------------------------------
countData <- countData_original
metaData <- metaData_original
dds <- dds_object(countData, metaData)
vsd <- vst(dds, blind=FALSE)
pcaplot(vsd, intgroup="group", ntop=1000, pcX=1, pcY=2, text_labels=0, ellipse=0) 

# ----------------------------------------------------------------------------------
# PCA plot by gastric cancer histology
#	(* Only for the surgery cohort)
# ----------------------------------------------------------------------------------
countData <- countData_original
metaData <- metaData_original
dds <- dds_object(countData, metaData)
vsd <- vst(dds, blind=FALSE)
pcaplot(vsd, intgroup="tumor_histology", ntop=1000, pcX=1, pcY=2, text_labels=0, ellipse=0) + stat_ellipse(type="norm", level=0.95)

# ----------------------------------------------------------------------------------
# Differential gene exprsesion analysis (* Only for the biopsy cohort)
# ----------------------------------------------------------------------------------
# Identifying significant genes between the groups
# 	Compared groups:
#		(1) HP_eradication vs. Healthy
#		(2) HP_gastritis vs. Healthy
#		(3) Gastric_cancer vs. Healthy
# ----------------------------------------------------------------------------------
gene_ID_sig <- NULL

# (1) HP_eradication vs. Healthy
comparison <- c("Healthy", "HP_eradication")
countData <- select_countData(countData_original, group_name=metaData_original$group, comparison=comparison)
metaData <- select_metaData(metaData_original, group_name=metaData_original$group, comparison)
dds <- dds_object(countData, metaData)
res <- results(dds) 
res[order(res$pvalue),]
siggene <- select_siggene(dds)
gene_ID_sig <- unique(c(gene_ID_sig, rownames(siggene)))
# KEGG pathway
res$entrez = mapIds(org.Hs.eg.db,	keys=row.names(res),	column="ENTREZID",	keytype="SYMBOL",	multiVals="first")
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
lapply(keggres, head, 20)

# (2) HP_gastritis vs. Healthy
comparison <- c("Healthy", "HP_gastritis")
countData <- select_countData(countData_original, group_name=metaData_original$group, comparison=comparison)
metaData <- select_metaData(metaData_original, group_name=metaData_original$group, comparison)
dds <- dds_object(countData, metaData)
res <- results(dds) 
res[order(res$pvalue),]
siggene <- select_siggene(dds)
gene_ID_sig <- unique(c(gene_ID_sig, rownames(siggene)))
# KEGG pathway
res$entrez = mapIds(org.Hs.eg.db,	keys=row.names(res),	column="ENTREZID",	keytype="SYMBOL",	multiVals="first")
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
lapply(keggres, head, 20)

# (3) Gastric_cancer vs. Healthy
comparison <- c("Healthy", "Gastric_cancer")
countData <- select_countData(countData_original, group_name=metaData_original$group, comparison=comparison)
metaData <- select_metaData(metaData_original, group_name=metaData_original$group, comparison)
dds <- dds_object(countData, metaData)
res <- results(dds) 
res[order(res$pvalue),]
siggene <- select_siggene(dds)
gene_ID_sig <- unique(c(gene_ID_sig, rownames(siggene)))
# KEGG pathway
res$entrez = mapIds(org.Hs.eg.db,	keys=row.names(res),	column="ENTREZID",	keytype="SYMBOL",	multiVals="first")
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
lapply(keggres, head, 20)

gene_ID_sig		# significant genes between any diseased group (HP_eradication, HP_gastritis, or Gastric_cancer) and Healthy group

# ----------------------------------------------------------------------------------
# Differential gene exprsesion analysis (* Only for the surgery cohort)
#   (Do not run this code if you are analyzing the biopsy cohort)
# ----------------------------------------------------------------------------------
# Identifying significant genes between the groups
# 	Compared group:
#		(1) Gastric_cancer vs. Severe_gastritis
# ----------------------------------------------------------------------------------
gene_ID_sig <- NULL

# (1) Gastric_cancer vs. Severe_gastritis
comparison <- c("Severe_gastritis", "Gastric_cancer")
countData <- select_countData(countData_original, group_name=metaData_original$group, comparison=comparison)
metaData <- select_metaData(metaData_original, group_name=metaData_original$group, comparison)
dds <- dds_object(countData, metaData)
res <- results(dds) 
res[order(res$pvalue),]
siggene <- select_siggene(dds)
gene_ID_sig <- unique(c(gene_ID_sig, rownames(siggene)))
# KEGG pathway
res$entrez = mapIds(org.Hs.eg.db,	keys=row.names(res),	column="ENTREZID",	keytype="SYMBOL",	multiVals="first")
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
lapply(keggres, head, 20)

gene_ID_sig		# significant genes between Gastric_cancer and Severe_gastritis groups

# ----------------------------------------------------------------------------------
# Make a "significant genes" table
# ----------------------------------------------------------------------------------
gene_ID_sig <- data.frame(Gene = gene_ID_sig)
countData <- countData_original
countData2 <- merge(countData, gene_ID_sig, by = "Gene")
rownames(countData2) <- countData2$Gene
countData2 <- countData2[,-1]
countData2 <- t(countData2)

################################################################################################
## CCA (assocation between ASVs and genes)
################################################################################################
library(PMA)
library(data.table)
library(patchwork)

countData2.scaled <- t(apply(countData2,1,scale))
colnames(countData2.scaled) <- colnames(countData2)
ASV_tab2.scaled <- t(apply(ASV_tab2,1,scale))
colnames(ASV_tab2.scaled) <- colnames(ASV_tab2)
# ----------------
hfs <- countData2
mfs <- ASV_tab2
o1s <- apply(scale(hfs),2,sd)
hfs <- hfs[,names(o1s)[!is.na(o1s)|o1s!=0.]]
o1s <- apply(scale(mfs),2,sd)
mfs <- mfs[,names(o1s)[!is.na(o1s)|o1s!=0.]]
in1 <- list(host=hfs,otu=mfs)
# ----------------
perm.out <- CCA.permute(hfs, mfs, typex="standard", typez="standard", nperms=7)
out <- CCA(hfs, mfs, typex="standard", typez="standard", K=2, niter=1000, trace=TRUE, xnames=colnames(hfs), znames=colnames(mfs), 
		penaltyx=perm.out$bestpenaltyx, penaltyz=perm.out$bestpenaltyz, v=perm.out$v.init)
sdx <- apply(hfs, 2, sd)
hfsz <- scale(hfs,TRUE,sdx)
sdx <- apply(mfs, 2, sd)
mfsz <- scale(mfs,TRUE,sdx)
cca_r2 <- cor(hfsz %*% out$u, mfsz %*% out$v)
message(sprintf("cca[U1,V1]=%g,cca[U2,V2]=%g",cca_r2[1,1],cca_r2[2,2]))
# ----------------
cca_variates <- as.data.frame(cbind(hfsz %*% out$u,mfsz %*% out$v))
colnames(cca_variates) <- c("U1","U2","V1","V2")
# ----------------
# barplot of CCA coefficient
cca_coeff <- list()
cca_coeff$U <- data.table(name=out$xnames,out$u)
colnames(cca_coeff$U) <-c("name","U1","U2")
cca_coeff$V <- data.table(name=out$znames,out$v)
colnames(cca_coeff$V) <-c("name","V1","V2")
# ----------------
cca_coeff_export <- cca_coeff; colnames(cca_coeff_export$U) <- c("Gene", "U1", "U2"); colnames(cca_coeff_export$V) <- c("OTU_ID", "V1", "V2")
cca_coeff_export$U <- cca_coeff_export$U[order(-abs(cca_coeff_export$U$U2)),]
cca_coeff_export$U <- cca_coeff_export$U[order(-abs(cca_coeff_export$U$U1)),]
cca_coeff_export$V <- cca_coeff_export$V[order(-abs(cca_coeff_export$V$V2)),]
cca_coeff_export$V <- cca_coeff_export$V[order(-abs(cca_coeff_export$V$V1)),]
cca_coeff_export$U
cca_coeff_export$V
# ----------------
# plotting (U1, V1)
plist <- list()
# CCA V1 coefficient >0 (ASV)
dt2 <- cca_coeff$V[V1>0,]
setnames(dt2,"V1","cca_coeff")
plist[[1]] <- barplot_cca(dt2,TRUE,"CCA V1 coefficient > 0",xlab="Coefficient",ylab="OTU",fillcol="#FF4848") + theme(plot.title=element_blank(), axis.title=element_blank(), axis.text=element_text(size=5), axis.ticks=element_line(size=0.1))
# CCA U1 coefficient >0 (Gene)
dt2 <- cca_coeff$U[U1>0,]
setnames(dt2,"U1","cca_coeff")
plist[[2]] <- barplot_cca(dt2,TRUE,"CCA U1 coefficient > 0",xlab="Coefficient",ylab="Gene",fillcol="#FF4848") +	theme(plot.title=element_blank(), axis.title=element_blank(), axis.text=element_text(size=5), axis.ticks=element_line(size=0.1))
# CCA V1 coefficient <0 (ASV)
dt2 <- cca_coeff$V[V1<0,]
setnames(dt2,"V1","cca_coeff")
plist[[3]] <- barplot_cca(dt2,FALSE,"CCA V1 coefficient < 0",xlab="Coefficient",ylab="OTU",fillcol="#368AFF") + theme(plot.title=element_blank(), axis.title=element_blank(), axis.text=element_text(size=5), axis.ticks=element_line(size=0.1))
# CCA U1 coefficient <0 (Gene)
dt2 <- cca_coeff$U[U1<0,]
setnames(dt2,"U1","cca_coeff")
plist[[4]] <- barplot_cca(dt2,FALSE,"CCA U1 coefficient < 0",xlab="Coefficient",ylab="Gene",fillcol="#368AFF") + theme(plot.title=element_blank(), axis.title=element_blank(), axis.text=element_text(size=5), axis.ticks=element_line(size=0.1))
p <- wrap_plots(plist,ncol=2,nrow=2)
plot(p)
# ----------------
# plotting (U2, V2)
plist <- list()
# CCA V1 coefficient >0 (ASV)
dt2 <- cca_coeff$V[V2>0,]
setnames(dt2,"V2","cca_coeff")
plist[[1]] <- barplot_cca(dt2,TRUE,"CCA V2 coefficient > 0",xlab="Coefficient",ylab="OTU",fillcol="#FF4848") + theme(plot.title=element_blank(), axis.title=element_blank(), axis.text=element_text(size=5), axis.ticks=element_line(size=0.1))
# CCA U1 coefficient >0 (Gene)
dt2 <- cca_coeff$U[U2>0,]
setnames(dt2,"U2","cca_coeff")
plist[[2]] <- barplot_cca(dt2,TRUE,"CCA U2 coefficient > 0",xlab="Coefficient",ylab="Gene",fillcol="#FF4848") +	theme(plot.title=element_blank(), axis.title=element_blank(), axis.text=element_text(size=5), axis.ticks=element_line(size=0.1))
# CCA V1 coefficient <0 (ASV)
dt2 <- cca_coeff$V[V2<0,]
setnames(dt2,"V2","cca_coeff")
plist[[3]] <- barplot_cca(dt2,FALSE,"CCA V2 coefficient < 0",xlab="Coefficient",ylab="OTU",fillcol="#368AFF") + theme(plot.title=element_blank(), axis.title=element_blank(), axis.text=element_text(size=5), axis.ticks=element_line(size=0.1))
# CCA U1 coefficient <0 (Gene)
dt2 <- cca_coeff$U[U2<0,]
setnames(dt2,"U2","cca_coeff")
plist[[4]] <- barplot_cca(dt2,FALSE,"CCA U2 coefficient < 0",xlab="Coefficient",ylab="Gene",fillcol="#368AFF") + theme(plot.title=element_blank(), axis.title=element_blank(), axis.text=element_text(size=5), axis.ticks=element_line(size=0.1))
p <- wrap_plots(plist,ncol=2,nrow=2)
plot(p)
# ----------------
# plotting (scatter plot for canonical coefficients)
title_pref <- "human GI host vs. 16S microbiome (U1,V1)"
title2 <- sprintf("CCA [%s],r2=%g",title_pref,cca_r2[1,1])
p <- ggplot(data = cca_variates,aes(U1,V1)) + 
	geom_point(size=1) +
	geom_text(
		label=rownames(cca_variates), 
		size=2, color="#5D5D5D", nudge_x = 0.25, nudge_y = -0.25, 
		check_overlap = T
	) +
	geom_smooth(method=lm, colour="#353535", size=0.5) +
	ggtitle(title2)
plist <- list()
plist[["cca1"]] <- p
title_pref <- "human GI host vs. 16S microbiome (U2,V2)"
title2 <- sprintf("CCA [%s],r2=%g",title_pref,cca_r2[2,2])
p <- ggplot(data = cca_variates,aes(U2,V2)) + 
	geom_point(size=1) +
	geom_text(
		label=rownames(cca_variates), 
		size=2, color="#5D5D5D", nudge_x = 0.25, nudge_y = -0.25, 
		check_overlap = T
	) +
	geom_smooth(method=lm, colour="#353535", size=0.5) +
	ggtitle(title2)

plist[["cca2"]] <- p
p <- wrap_plots(plist,ncol=2)
plot(p)
# ---------------
summary (lm(formula = cca_variates$V1 ~ cca_variates$U1))
summary (lm(formula = cca_variates$V2 ~ cca_variates$U2))

################################################################################################
## Cell enrichment analysis (xCell)
################################################################################################
library(devtools)
library(xCell)

metaData <- metaData_original
metaData$group <- factor(metaData$group, levels=c("Healthy", "HP_eradication", "HP_gastritis", "Severe_gastritis", "Gastric_cancer"))
metaData_xCell <- data.frame(sample_ID=metaData$sample_ID, group=metaData$group)

tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}
countData_TPM <- tpm(countData_original[,-1], dim(countData_original)[1])
rownames(countData_TPM) <- countData_original$Gene

raw.scores = rawEnrichmentAnalysis(as.matrix(countData_TPM),
	xCell.data$signatures,
	xCell.data$genes)

non.spill.over <- c("aDC", "B-cells", "Basophils", "CD4+ memory T-cells", "CD4+ naive T-cells", "CD4+ T-cells", "CD4+ Tcm", "CD4+ Tem", "CD8+ naive T-cells", "CD8+ T-cells",
				"CD8+ Tcm", "CD8+ Tem", "cDC", "Class-switched memory B-cells", "DC", "Endothelial cells", "Eosinophils", "Epithelial cells", "Erythrocytes",
				"Fibroblasts", "iDC", "ly Endothelial cells", "Macrophages", "Macrophages M1", "Macrophages M2", "Mast cells", "Memory B-cells", "Monocytes",
				"mv Endothelial cells", "naive B-cells", "Neutrophils", "NK cells", "NKT", "pDC", "Pericytes", "Plasma cells", "Tgd cells", "Th1 cells", "Th2 cells", "Tregs")
non.spill.over.rename <- c("Myeloid dendritic cell activated", "B cell", "Basophil", "T cell CD4+ memory", "T cell CD4+ naive", "T cell CD4+ (non-regulatory)", 
					"T cell CD4+ central memory", "T cell CD4+ effector memory", "T cell CD8+ naive", "T cell CD8+", "T cell CD8+ central memory", 
					"T cell CD8+ effector memory", "Myeloid dendritic cell conventional", "Class-switched memory B cell", "Myeloid dendritic cell", 
					"Endothelial cell", "Eosinophil", "Epithelial cell", "Erythrocyte", "Fibroblast", "Myeloid dendritic cell immature", "Lymphatic endothelial cell", 
					"Macrophage", "Macrophage M1", "Macrophage M2", "Mast cell", "B cell memory", "Monocyte", "Microvascular endothelial cell", "B cell naive", 
					"Neutrophil", "NK cell", "T cell NK", "Plasmacytoid dendritic cell", "Pericyte", "B cell plasma", "T cell gamma delta", "T cell CD4+ Th1",
					"T cell CD4+ Th2", "T cell regulatory")
cell.types.use = intersect(colnames(xCell.data$spill$K), non.spill.over)
transformed.scores = transformScores(raw.scores[cell.types.use,], xCell.data$spill$fv)
scores = spillOver(transformed.scores, xCell.data$spill$K)
# ---------------
xCell <- t(scores)
colnames(xCell) <- non.spill.over.rename
xCell.scaled <- t(apply(xCell,1,scale))
colnames(xCell.scaled) <- colnames(xCell)
# ----------------
ASV_tab2.scaled <- t(apply(ASV_tab2,1,scale))
colnames(ASV_tab2.scaled) <- colnames(ASV_tab2)

# ------------------------------------------------------------------------------------------
# PCA for xCell scores	(* Only for the biopsy cohort)
# ------------------------------------------------------------------------------------------
library(grid)
library(gridExtra)

xCell.scaled_pca <- prcomp(xCell.scaled)
xCell.scaled_out <- as.data.frame(xCell.scaled_pca$x)
xCell.scaled_out$group <- metaData$group
xCell.scaled_out$group <- ifelse(xCell.scaled_out$group=="Healthy", "Healthy stomach", "Non-healthy stomach")
percentage <- round(xCell.scaled_pca$sdev / sum(xCell.scaled_pca$sdev) * 100, 1)
percentage <- paste(colnames(xCell.scaled_out), "(", paste(as.character(percentage), "%", ")", sep="") )

p <-  ggplot(xCell.scaled_out,aes(x=PC1,y=PC2,color=group )) +
	geom_point(size=4, alpha=0.7) + xlab(percentage[1]) + ylab(percentage[2]) +
	theme(aspect.ratio=1, axis.title=element_text(size=10), axis.text=element_text(size=10), legend.title=element_blank(), legend.text=element_text(size=10), axis.ticks=element_line(size=0.1))
p

# ------------------------------------------------------------------------------------------
# PCA for xCell scores	(* Only for the surgery cohort)
# ------------------------------------------------------------------------------------------
library(grid)
library(gridExtra)

xCell.scaled_pca <- prcomp(xCell.scaled)
xCell.scaled_out <- as.data.frame(xCell.scaled_pca$x)
xCell.scaled_out$group <- metaData$group
percentage <- round(xCell.scaled_pca$sdev / sum(xCell.scaled_pca$sdev) * 100, 1)
percentage <- paste(colnames(xCell.scaled_out), "(", paste(as.character(percentage), "%", ")", sep="") )

p <-  ggplot(xCell.scaled_out,aes(x=PC1,y=PC2,color=group )) +
	geom_point(size=4, alpha=0.7) + xlab(percentage[1]) + ylab(percentage[2]) +
	theme(aspect.ratio=1, axis.title=element_text(size=10), axis.text=element_text(size=10), legend.title=element_blank(), legend.text=element_text(size=10), axis.ticks=element_line(size=0.1))
p
# ----------------
# PCA plot, clustering by gastric cancer histology
xCell.scaled_out$group <- metaData$tumor_histology
p <-  ggplot(xCell.scaled_out,aes(x=PC1,y=PC2,color=group )) +
	geom_point(size=4, alpha=0.7) + xlab(percentage[1]) + ylab(percentage[2]) +
	theme(aspect.ratio=1, axis.title=element_text(size=10), axis.text=element_text(size=10), legend.title=element_blank(), legend.text=element_text(size=10), axis.ticks=element_line(size=0.1)) +
	  stat_ellipse(type="norm", level=0.95) + theme_bw()
p

# ------------------------------------------------------------------------------------------
# Bar graph for xCell scores
# ------------------------------------------------------------------------------------------
# plotting: xCell score
library(reshape)
library(RColorBrewer)

xCell.melted <- melt(xCell)
colnames(xCell.melted) <- c("sample_ID", "immune_cell", "relative_abundance")
xCell.melted <- merge(xCell.melted, metaData_xCell, by = "sample_ID")

immune_cell_type <- levels(as.factor(xCell.melted$immune_cell))
getPalette <- colorRampPalette(brewer.pal(9, "Set1")) 
col_immune_cell <- getPalette(length(immune_cell_type))
col_group <- hue_pal()(4)

sample_ID2 <- gsub("N", "SG", xCell.melted$sample_ID)	# for the surgery cohort (not affected for the biopsy cohort)
sample_ID2 <- gsub("T", "GC", sample_ID2)			# for the surgery cohort (not affected for the biopsy cohort)
xCell.melted$sample_ID <- sample_ID2
ggplot(xCell.melted, aes(x=sample_ID, y=relative_abundance, color=immune_cell, fill=immune_cell)) + 
	geom_bar(stat = "identity") + facet_wrap("group~.", scales="free_x", nrow=1)  +
	theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
	scale_fill_manual(values = col_immune_cell) + scale_color_manual(values = col_immune_cell) +
	theme(legend.text=element_text(size=15), axis.ticks=element_line(size=0.1))

# ------------------------------------------------------------------------------------------
# 2 way hierachical clustering for xCell data (* Only for the biopsy cohort)
# ------------------------------------------------------------------------------------------
xCell.scaled.clip <- scale_and_clip(t(xCell))
ann_colors = list(group = c("Healthy stomach"="#F8766D", "Non-healthy stomach"="#00BFC4"))
annotation <- data.frame(group=metaData$group)
annotation <- data.frame(group=ifelse(annotation$group=="Healthy", "Healthy stomach", "Non-healthy stomach"))
rownames(annotation) <- metaData$sample_ID

myCol <- get_red_white_blue_for_heatmap(xCell.scaled.clip, paletteLength=49)	
pheatmap(xCell.scaled.clip, fontsize=9, fontsize_row=9, fontsize_col=9, color=myCol$colors, breaks=myCol$breaks,
			cutree_cols=3, cutree_rows=2,
			annotation_col=annotation, annotation_colors = ann_colors, annotation_legend=TRUE)	

# ------------------------------------------------------------------------------------------
# 2 way hierachical clustering for xCell data (* Only for the surgery cohort)
# ------------------------------------------------------------------------------------------
xCell.scaled.clip <- scale_and_clip(t(xCell))
ann_colors = list(group = c(Severe_gastritis="#F8766D", Gastric_cancer="#00BFC4"))
sample_ID2 <- gsub("N", "SG", rownames(xCell))
sample_ID2 <- gsub("T", "GC", sample_ID2)
colnames(xCell.scaled.clip) <- sample_ID2
metaData$sample_ID2 <- sample_ID2
annotation <- data.frame(group=metaData$group)
rownames(annotation) <- metaData$sample_ID2

myCol <- get_red_white_blue_for_heatmap(xCell.scaled.clip, paletteLength=49)	
pheatmap(xCell.scaled.clip, fontsize=9, fontsize_row=9, fontsize_col=9, color=myCol$colors, breaks=myCol$breaks,
			cutree_cols=3, cutree_rows=2,
			annotation_col=annotation, annotation_colors = ann_colors, annotation_legend=TRUE)	

################################################################################################
## Correlation between ASVs and Immune cells
################################################################################################
library(ggcorrplot)

# ------------------------------------------------------------------------------------------
# Heatmap for corrleation
# ------------------------------------------------------------------------------------------
ASV_tab2.scaled <- t(apply(ASV_tab2,1,scale))
colnames(ASV_tab2.scaled) <- colnames(ASV_tab2)

m <- cbind(ASV_tab2.scaled, xCell.scaled)
corr_original  <- cor(m, method = "spearman")
pvalue_original <- cor_pmat(m)

corr_data <- corr_original[1:dim(ASV_tab2.scaled)[2], (dim(ASV_tab2.scaled)[2]+1):(dim(ASV_tab2.scaled)[2]+dim(xCell.scaled)[2])]
pvalue_data <- pvalue_original[1:dim(ASV_tab2.scaled)[2], (dim(ASV_tab2.scaled)[2]+1):(dim(ASV_tab2.scaled)[2]+dim(xCell.scaled)[2])]

myCol <- get_red_white_blue_for_heatmap(corr_data, paletteLength=49)	
pheatmap(t(corr_data), fontsize=16, fontsize_row=10, fontsize_col=6,
		color=myCol$colors, breaks=myCol$breaks, cutree_cols = 2, cutree_rows = 2,
		display_numbers=matrix(ifelse(t(pvalue_data) < 0.05, "*", ""), nrow(t(pvalue_data))) 	)

# ------------------------------------------------------------------------------------------
# Top 20 high correlations according to clinical phenotype
# ------------------------------------------------------------------------------------------
n_control <- NULL
n_control <- nrow(metaData[metaData$group=="Healthy",])		# * for the biopsy cohort
n_control <- nrow(metaData[metaData$group=="Severe_gastritis",])	# * for the surgery 

xCell.scaled.sub1 <- xCell.scaled[1:n_control,]
xCell.scaled.sub2 <- xCell.scaled[(n_control+1):nrow(metaData),]
ASV_tab2.scaled.sub1 <- ASV_tab2.scaled[1:n_control,]
ASV_tab2.scaled.sub2 <- ASV_tab2.scaled[(n_control+1):nrow(metaData),]

m.sub1 <- cbind(ASV_tab2.scaled.sub1, xCell.scaled.sub1)
m.sub2 <- cbind(ASV_tab2.scaled.sub2, xCell.scaled.sub2)

corr_original.sub1  <- cor(m.sub1, method = "spearman")
pvalue_original.sub1 <- cor_pmat(m.sub1)
corr_original.sub2  <- cor(m.sub2, method = "spearman")
pvalue_original.sub2 <- cor_pmat(m.sub2)

corr_data.sub1 <- corr_original.sub1[1:dim(ASV_tab2.scaled)[2], (dim(ASV_tab2.scaled)[2]+1):(dim(ASV_tab2.scaled)[2]+dim(xCell.scaled)[2])]
pvalue_data.sub1 <- pvalue_original.sub1[1:dim(ASV_tab2.scaled)[2], (dim(ASV_tab2.scaled)[2]+1):(dim(ASV_tab2.scaled)[2]+dim(xCell.scaled)[2])]
corr_data.sub2 <- corr_original.sub2[1:dim(ASV_tab2.scaled)[2], (dim(ASV_tab2.scaled)[2]+1):(dim(ASV_tab2.scaled)[2]+dim(xCell.scaled)[2])]
pvalue_data.sub2 <- pvalue_original.sub2[1:dim(ASV_tab2.scaled)[2], (dim(ASV_tab2.scaled)[2]+1):(dim(ASV_tab2.scaled)[2]+dim(xCell.scaled)[2])]
# ----------------
corr_data.sub1.melted <- melt(corr_data.sub1)
corr_data.sub1.melted <- corr_data.sub1.melted[order(-abs(corr_data.sub1.melted$value)),]
corr_data.sub1.melted <- corr_data.sub1.melted[1:20,]
corr_data.sub1.melted <- corr_data.sub1.melted[order(-corr_data.sub1.melted$value),]
corr_data.sub1.melted$comparison <- paste(corr_data.sub1.melted$X2,"  ::  ", corr_data.sub1.melted$X1, sep='')
corr_data.sub1.melted$color <- ifelse(corr_data.sub1.melted$value>=0, "positive", "negative")

order <- corr_data.sub1.melted$comparison[order(corr_data.sub1.melted$value)]
corr_data.sub1.melted$comparison <- factor(corr_data.sub1.melted$comparison, levels = order)

ggplot(corr_data.sub1.melted, aes(x=value, y=comparison, fill=color)) +
	geom_bar(stat='identity') +
	scale_fill_manual(values=c("#368AFF", "#FF4848")) + theme(axis.ticks=element_line(size=0.1), legend.position="none") + 
	xlab("Correlation coefficient") + ylab("Cell type :: Bacterial ASV")
# ----------------
corr_data.sub2.melted <- melt(corr_data.sub2)
corr_data.sub2.melted <- corr_data.sub2.melted[order(-abs(corr_data.sub2.melted$value)),]
corr_data.sub2.melted <- corr_data.sub2.melted[1:20,]
corr_data.sub2.melted <- corr_data.sub2.melted[order(-corr_data.sub2.melted$value),]
corr_data.sub2.melted$comparison <- paste(corr_data.sub2.melted$X2,"  ::  ", corr_data.sub2.melted$X1, sep='')
corr_data.sub2.melted$color <- ifelse(corr_data.sub2.melted$value>=0, "positive", "negative")

order <- corr_data.sub2.melted$comparison[order(corr_data.sub2.melted$value)]
corr_data.sub2.melted$comparison <- factor(corr_data.sub2.melted$comparison, levels = order)

ggplot(corr_data.sub2.melted, aes(x=value, y=comparison, fill=color)) +
	geom_bar(stat='identity') +
	scale_fill_manual(values=c("#368AFF", "#FF4848")) + theme(axis.ticks=element_line(size=0.1), legend.position="none") + 
	xlab("Correlation coefficient") + ylab("Cell type :: Bacterial ASV")

# ------------------------------------------------------------------------------------------