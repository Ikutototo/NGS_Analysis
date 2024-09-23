setwd("/Users/saitoikuto/Documents/015.RStudio/NGS_Practice/The phyllosphere microbiome shifts toward combating melanose pathogen")

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(tidyverse); packageVersion("tidyverse")
library(ggtext); packageVersion("ggtext")

load(file = "start_phyloseq.RData")

EP_subset = subset_taxa(EP_phyloseq, Order!="Chloroplast")
EP_subset = subset_taxa(EP_subset, Family!="Mitochondria")

EP_RA = transform_sample_counts(EP_subset, function(x) x / sum(x) )

EP_RA_filt <- filter_taxa(EP_RA, function(x) mean(x) > 0.001, TRUE)
EP_RA_filt = transform_sample_counts(EP_RA_filt, function(x) x / sum(x) )

taxa <- as.data.frame(EP_RA_filt@tax_table@.Data)

EP_RA_filt_data <- as.data.frame(t(EP_RA_filt@otu_table@.Data)) %>% 
  select(contains("EP_16S_D"), everything())

Mean_infected <- EP_RA_filt_data %>% 
  select(contains("EP_16S_D")) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate_all(~ . * 100) %>% 
  apply(., 2, mean)

EP_RA_filt_data <- EP_RA_filt_data %>% 
  mutate(Mean_infected = Mean_infected)

Mean_uninfected <- EP_RA_filt_data %>% 
  select(contains("EP_16S_H")) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate_all(~ . * 100) %>% 
  apply(., 2, mean)

EP_RA_filt_data <- EP_RA_filt_data %>% 
  mutate(Mean_uninfected = Mean_uninfected)

Mean_all <- EP_RA_filt_data %>% 
  select(contains("EP_16S")) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate_all(~ . * 100) %>% 
  apply(., 2, mean)

EP_RA_filt_data <- EP_RA_filt_data %>% 
  mutate(Mean_all_samples = Mean_all)

EP_RA_filt_data <- bind_cols(EP_RA_filt_data, taxa)

rm(list = setdiff(ls(), "EP_RA_filt_data"))
              
EP_RA_filt_data[1,][grep("EP_16S_D", colnames(EP_RA_filt_data), value = TRUE)]

EP_RA_filt_data$p.value <- EP_RA_filt_data %>% 
  dplyr::select(!c(Mean_infected, Mean_uninfected)) %>% 
  apply(., 1, \(x){
  EP_D <- as.numeric(x[grep("EP_16S_D", colnames(EP_RA_filt_data), value = TRUE)])
  EP_H <- as.numeric(x[grep("EP_16S_H", colnames(EP_RA_filt_data), value = TRUE)])
  result <- wilcox.test(EP_D,EP_H)
  result$p.value
})

EP_RA_filt_data <- EP_RA_filt_data %>% 
  mutate(log10_p.value = -log10(p.value)) %>% 
  mutate(log2_FoldChange = log2(Mean_infected/Mean_uninfected)) 
  
EP_RA_filt_data <- EP_RA_filt_data %>% 
  mutate(ratio = Mean_infected / Mean_uninfected) 


EP_RA_filt_data$ASV_ID <- rownames(EP_RA_filt_data)


EP_RA_filt_data$diff <- EP_RA_filt_data %>% 
  dplyr::select(c(p.value, log2_FoldChange)) %>% 
  apply(., 1, \(x){
    if (x["p.value"] > 0.05) {
      "NotSig"
    } else if(x["log2_FoldChange"] < 0){
      "Depleted"
    } else{
      "Enriched"
    }
  }) 

EP_RA_filt_data <- EP_RA_filt_data %>% 
  mutate(Mean_all = Mean_infected + Mean_uninfected)

EP_RA_filt_data <- EP_RA_filt_data %>% 
  select(ASV_ID, log2_FoldChange, p.value, log10_p.value,
         diff, ratio , Kingdom, Phylum, Class, Order, Family, Genus, Species,
         Mean_infected, Mean_uninfected, Mean_all, Mean_all_samples, everything())

rownames(EP_RA_filt_data) <- NULL
colnames(EP_RA_filt_data)

write.csv(x = EP_RA_filt_data, file = "final_data.csv", row.names = FALSE)



ggplot(EP_RA_filt_data, aes(
  x = Genus, y = log10_p.value, 
  size = log10_p.value,
  shape = as.factor(case_when(
    diff == "NotSig" ~ 16, 
    diff == "Depleted" ~ if_else(ratio < 0.8, 25, 6), # FALSEになるASVがない
    diff == "Enriched" ~ if_else(ratio > 1.5, 24, 17) # FALSEになるASVがない
  )))) + 
  geom_hline(
    yintercept = -log10(0.05), color = "grey40",
    linetype = "dashed"
  ) +
  geom_point(alpha = 0.8, position = position_jitter(width = 1),
             aes(color = as_factor(Genus))) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

