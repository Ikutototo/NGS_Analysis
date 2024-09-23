

# save.RData --------------------------------


# load.RData --------------------------------

setwd("/Users/saitoikuto/Documents/015.RStudio/NGS_Practice/The phyllosphere microbiome shifts toward combating melanose pathogen")

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(tidyverse); packageVersion("tidyverse")


load(file = "start_phyloseq.RData")


# Accessors ---------------------------------

ntaxa(EN_EP_phyloseq)
nsamples(EN_EP_phyloseq)
sample_names(EN_EP_phyloseq)
rank_names(EN_EP_phyloseq)
sample_variables(EN_EP_phyloseq)
otu_table(EN_EP_phyloseq)
tax_table(EN_EP_phyloseq)
phy_tree(EN_EP_phyloseq)
taxa_names(EN_EP_phyloseq)
taxa_sums(EN_EP_phyloseq)
prune_taxa(names(sort(taxa_sums(EN_EP_phyloseq),
                      decreasing = TRUE)[1:10])
           ,EN_EP_phyloseq)

names(sort(taxa_sums(EN_EP_phyloseq), decreasing = TRUE)[1:10])


# PreProcess --------------------------------

EN_EP_phyloseq_otu_table_.Data <- EN_EP_phyloseq@otu_table@.Data
EN_EP_phyloseq_tax_table_.Data <- EN_EP_phyloseq@tax_table@.Data

(EN_EP_phyloseq_subset = subset_taxa(EN_EP_phyloseq, Order!="Chloroplast"))
(EN_EP_phyloseq_subset = subset_taxa(EN_EP_phyloseq_subset, Family!="Mitochondria"))

EN_EP_phyloseq_RA = transform_sample_counts(EN_EP_phyloseq_subset, function(x) x / sum(x) )
EN_EP_phyloseq_RA_filter = filter_taxa(EN_EP_phyloseq_RA, function(x) mean(x) > 0.01, TRUE)


a <- EN_EP_phyloseq_RA@otu_table@.Data
b <- EN_EP_phyloseq_RA@tax_table@.Data
c <- EN_EP_phyloseq_RA_filter@otu_table@.Data

total = median(sample_sums(EN_EP_phyloseq_subset))
standf = function(x, t=total) round(t * (x / sum(x)))
EN_EP_phyloseq_z = transform_sample_counts(EN_EP_phyloseq_subset, standf)
a <- EN_EP_phyloseq_z@otu_table@.Data

EN_EP_phyloseq_cv = filter_taxa(EN_EP_phyloseq_z, function(x) sd(x)/mean(x) > 3.0, TRUE)

a <- EP_phyloseq@otu_table@.Data
b <- EP_phyloseq@tax_table@.Data



c <- EP_subset@otu_table@.Data
d <- EP_subset@tax_table@.Data

e <- EP_subset_RA@otu_table@.Data
f <- EP_subset_RA@tax_table@.Data


samEP_df <- sample_data(EP_phyloseq)
samEP_df$Disease
rownames(subset(samEP_df, samEP_df$Disease == "H"))

filtered_phyloseq_object <- prune_samples(rownames(subset(samEP_df, samEP_df$Disease == "H")),
                                          samEP_df)
prune_samples(sample_data(EP_phyloseq)$Disease == "H", EP_phyloseq)


#------

EP_subset = subset_taxa(EP_phyloseq, Order!="Chloroplast")
EP_subset = subset_taxa(EP_subset, Family!="Mitochondria")

EP_Hl_subset <- prune_samples(sample_data(EP_subset)$Disease == "H", EP_subset)
EP_Di_subset <- prune_samples(sample_data(EP_subset)$Disease == "D", EP_subset)

EP_Hl_RA = transform_sample_counts(EP_Hl_subset, function(x) x / sum(x) )
EP_Di_RA = transform_sample_counts(EP_Di_subset, function(x) x / sum(x) )
EP_RA = transform_sample_counts(EP_subset, function(x) x / sum(x) )

EP_Hl_RA_filt <- filter_taxa(EP_Hl_RA, function(x) mean(x) > 0.001, TRUE)
EP_Di_RA_filt <- filter_taxa(EP_Di_RA, function(x) mean(x) > 0.001, TRUE)
EP_RA_filt <- filter_taxa(EP_RA, function(x) mean(x) > 0.001, TRUE)
#EP_RA_filt <- filter_taxa(EP_RA, function(x) max(x) > 0.001, TRUE)


i <- t(EP_Hl_RA_filt@otu_table@.Data)
j <- t(EP_Di_RA_filt@otu_table@.Data)
k <- t(EP_RA_filt@otu_table@.Data)


l <- EP_Hl_RA_filt@tax_table@.Data
m <- EP_Di_RA_filt@tax_table@.Data
n <- EP_RA_filt@tax_table@.Data

l <- as.data.frame(l)
l$ASV <- rownames(l)
colnames(l)

m <- as.data.frame(m)
m$ASV <- rownames(m)

df_ml <- inner_join(l,m, by = "ASV")
colnames(df_ml)
df_ml <- df_ml %>% select("ASV","Kingdom.x","Phylum.x","Class.x",
                          "Order.x","Family.x","Genus.x","Species.x")
colnames(df_ml) <- c("ASV","Kingdom","Phylum","Class","Order","Family","Genus","Species")


mean = colSums(t(j)*100/10)
data = bind_cols(j*100,Mean_infected=mean)
data = as.data.frame(data)
data$ASV <- rownames(j)

mean2 = colSums(t(i)*100/10)
data2 = bind_cols(i*100,Mean_uninfected=mean2)
data2 = as.data.frame(data2)
data2$ASV <- rownames(i)


merge_df <- inner_join(data,data2, by = "ASV")
Data_ENEP <- merge_df %>% select(c("ASV", "Mean_infected",  "Mean_uninfected", 
                                   "EP_16S_D_1a","EP_16S_D_1j","EP_16S_D_1i","EP_16S_D_1h","EP_16S_D_1g",
                                   "EP_16S_D_1f","EP_16S_D_1e","EP_16S_D_1d","EP_16S_D_1c","EP_16S_D_1b",   
                                   "EP_16S_H_1j","EP_16S_H_1i","EP_16S_H_1h","EP_16S_H_1g","EP_16S_H_1f",
                                   "EP_16S_H_1e","EP_16S_H_1d","EP_16S_H_1c","EP_16S_H_1b","EP_16S_H_1a"))


Data <- Data %>% mutate(log2_Fold_Change = log2(Mean_infected/Mean_uninfected))


write.csv(x = Data, file = "data_wilcox_test.csv", )

rownames(merge_df) <- merge_df$ASV
merge_df$ASV <- NULL
merge_df <- t(merge_df)








library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library("vegan"); packageVersion("vegan")


dist_methods <- unlist(distanceMethodList)
print(dist_methods)


# These require tree
dist_methods[(1:3)]
# Remove them from the vector
dist_methods <- dist_methods[-(1:3)]
# This is the user-defined method:
dist_methods["designdist"]
# Remove the user-defined distance
dist_methods = dist_methods[-which(dist_methods=="ANY")]

plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- phyloseq::distance(enterotype, method=i)
  # Calculate ordination
  iMDS  <- phyloseq::ordinate(enterotype, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- phyloseq::plot_ordination(enterotype, iMDS, color="SeqTech", shape="Enterotype")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}


df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=SeqTech, shape=Enterotype))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for Enterotype dataset")
p


