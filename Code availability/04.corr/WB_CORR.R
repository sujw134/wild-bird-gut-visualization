rm(list = ls())

options(stringsAsFactors = F)
pacman::p_load(tidyverse,ggrepel,ggforce,vegan,dplyr,tidyr, tibble, purrr)
setwd("05/")

##################################################################################gene_arrow################################################################
# ==========================================
# 1.
# ==========================================

#
gene_info_args <- read.delim("ARGs_info_rename.txt") %>% 
  select(gene_id = Gene, gene_name = ARG) %>%
  distinct()

# 
gene_info_mges <- read.delim("all.info", sep = "\t", header = FALSE) %>% 
  rename(gene_id = V1, mge_name = V2) %>% 
  select(gene_id, gene_name = mge_name) %>%
  distinct()

# ==========================================
# 2
# ==========================================

# 
colinearity <- read.delim("arg_mge_find.tsv", header = FALSE, 
                          col.names = c("neighborhood_id", "gene_id", "type", "start", "end", "strand")) %>%  
  filter(type != "Other") %>% 
  select(gene_id, type) %>%
  distinct() %>%
  left_join(gene_info_args %>% rename(name_arg = gene_name), by = "gene_id") %>%
  left_join(gene_info_mges %>% rename(name_mge = gene_name), by = "gene_id") %>%
  mutate(final_gene_name = coalesce(name_arg, name_mge))

#
gene_map_cleaned <- colinearity %>%
  group_by(gene_id) %>%
  summarize(
    label_name = first(na.omit(final_gene_name)),
    final_type = {
      all_types <- unlist(strsplit(as.character(type), "\\|"))
      paste(sort(unique(all_types)), collapse = "|")
    },
    .groups = "drop"
  )

# ==========================================
# 3
# ==========================================
dat <- read.delim("arg_mge_find.tsv", header = FALSE, 
                  col.names = c("neighborhood_id", "gene_id", "type", "start", "end", "strand"))

plot_data <- dat %>%
  left_join(gene_map_cleaned, by = "gene_id") %>%
  mutate(
    is_forward = ifelse(strand == "+", TRUE, FALSE),
    final_type = ifelse(is.na(final_type), "Other", final_type),
    label_text = ifelse(final_type == "Other", "", label_name)
  )

# ==========================================
# 4
# ==========================================
library(gggenes)

p <- ggplot(plot_data, aes(xmin = start, xmax = end, y = neighborhood_id, 
                           fill = final_type, forward = is_forward)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), 
                  arrowhead_width = unit(1.5, "mm"),
                  arrow_body_height = unit(2, "mm"),
                  color = "black", size = 0.2) +
  geom_text(aes(x = (start + end)/2, label = label_text), 
            size = 2, vjust = -1.5, fontface = "italic") +
  facet_wrap(~ neighborhood_id, scales = "free", ncol = 1) +
  
  scale_fill_manual(values = c(
    "ARGs" = "#F59B7B",       
    "MGEs" = "#FCC41E",       
    "ARGs|MGEs" = "#c2bdde",  
    "Other" = "#FFFFFF"       
  )) +
  
  theme_genes() +
  theme(
    axis.text.y = element_blank(),
    strip.text = element_text(size = 6, face = "bold", hjust = 0),
    legend.position = "top",
    panel.spacing = unit(0.3, "lines")
  ) +
  labs(title = "ARG-MGE Genetic Neighborhood",
       fill = "Type", x = "Genomic Position (bp)")

# 保存
ggsave("ARG_MGE_Map_Fixed2.pdf", p, width = 5, height = 10, limitsize = FALSE)
#########################################################################################################################################################################
###################################### co-abundance ##############################################################################
#
source("profile_process.R")
metadata_args <- read.delim("../ARGs_info_rename.txt") 
metadata_mges <- read.delim("../all.info",sep = "\t",header = F) %>% 
  rename(gene = V1, mag = V5, mge = V2,mge_type = V3, mge_type2 = V4 )

#
tem <- read.delim("../02/ARGs_prevalence.txt") %>% 
  .[order(.$avg_tpm, decreasing = TRUE),] %>% head(20) 

profile_args <- read.delim("../ARO_abundance.txt", row.names = 1) %>%
  mutate(ARG=metadata_args$ARG[match(rownames(.),metadata_args$ARO)])%>% 
  group_by(ARG) %>% summarise_all(sum) %>% column_to_rownames(var = "ARG") %>% 
  mutate(tem_value=tem$ARG[match(rownames(.),tem$ARG)]) %>% filter(!is.na(tem_value)) %>% 
  select(-tem_value) %>% t()%>% data.frame(check.names = F)

tem2 <- read.delim("../MGEs_prevalence.txt") %>% 
  .[order(.$avg_tpm, decreasing = TRUE),] %>% head(20) 


profile_mges <- read.delim("../MGEs_abundance.txt",sep = "\t",row.names = 1) %>% 
  t()%>% data.frame(check.names = F) %>% 
  select(all_of(tem2$name))

####install.packages("psych")###########
source("corr_process.R")
corr <- psych::corr.test(x = profile_args, y = profile_mges, method = "spearman", adjust = "BH")

corr_res <- corr_process(corr, rho = 0, padj = 0.05, out_df = T) %>%
  mutate(direct = case_when(
    r >= 0.6  ~ "pos",
    r <= -0.6 ~ "neg"
  ))%>%
  drop_na(direct)  

#write.table(corr_res, "top20prevalence_mge_arg_corr_res.tsv", sep = "\t", row.names = F, quote = F)

####
pacman::p_load(tidygraph,ggraph,igraph)
plotdat<- corr_res 
edges <- plotdat %>%rename(from = name_x, to = name_y) %>% 
  mutate(abs=abs(r),dir = factor(direct, levels = c("pos", "neg")))

nodes <-rbind(data.frame(name = edges$from) %>% 
                add_column(type = "ARG"),
              data.frame(name = edges$to) %>% 
                add_column(type = "MGE")) %>% 
  mutate(type = factor(type, labels = c("ARG", "MGE")))

graph <- tbl_graph(nodes = nodes, edges = edges)
nodes$degree <- degree(graph)    
nodes2 <-nodes[nodes$degree!=0,]
graph <- tbl_graph(nodes = nodes2, edges = edges)

class_color <- 
  set.seed(2025)


##################################################################
ggraph(graph, layout = "nicely") + 
  geom_edge_link(aes(edge_width = abs, edge_linetype = dir), 
                 color = "#F59B7B", show.legend = T) + 
  scale_edge_width(range = c(.1, 1)) +
  scale_edge_alpha_manual(values = c(0.2, 0.4)) +
  geom_node_point(aes(fill = type, size = degree), 
                  shape = 21, stroke = .1) +
  scale_size(range = c(6, 15)) +
  #scale_color_manual(values = c("#ED8828","#FFE4B5" )) +
  scale_fill_manual(values = c("ARG" = "#ED8828", "MGE" = "#FFE4B5")) +
  geom_node_text(aes(label = name),size = 1.5, fontface = "plain") +
  theme(panel.background = element_rect(fill = NA, color = "black"),
        legend.key = element_rect(fill = NA),
        aspect.ratio = 1)
#ggsave("corr_cooccurence_50.pdf", width = 10, height = 6)

###############################################################arg_mge heatmap###############################################################################
metadata_args <- read.delim("../ARGs_info_rename.txt") 
metadata_mges <- read.delim("../all.info",sep = "\t",header = F) %>% 
  rename(gene = V1, mag = V5, mge = V2,mge_type = V3, mge_type2 = V4 )

#选取耐药流行率前xxx
tem <- read.delim("../ARGs_prevalence.txt") %>% 
  .[order(.$avg_tpm, decreasing = TRUE),] %>% head(100) 

profile_args <- read.delim("../ARO_abundance.txt", row.names = 1) %>%
  mutate(ARG=metadata_args$ARG[match(rownames(.),metadata_args$ARO)])%>% 
  group_by(ARG) %>% summarise_all(sum) %>% column_to_rownames(var = "ARG") %>% 
  mutate(tem_value=tem$ARG[match(rownames(.),tem$ARG)]) %>% filter(!is.na(tem_value)) %>% 
  select(-tem_value) %>% t()%>% data.frame(check.names = F)

#
tem2 <- read.delim("../MGEs_prevalence.txt") %>% 
  .[order(.$avg_tpm, decreasing = TRUE),] %>% head(700) 


profile_mges <- read.delim("../MGEs_abundance.txt",sep = "\t",row.names = 1) %>% 
  #profile_filter(min_abundance = 1, min_prevalence = 0.05, by_group = F) %>% 
  t()%>% data.frame(check.names = F) %>% 
  select(all_of(tem2$name))

# 
####install.packages("psych")###########
source("corr_process.R")
corr <- psych::corr.test(x = profile_args, y = profile_mges, method = "spearman", adjust = "BH")

corr_res <- corr_process(corr, rho = 0, padj = 0.05, out_df = T) %>%
  #mutate(direct = ifelse(r > 0.5, "pos", "neg"))
  mutate(direct = case_when(
    r >= 0.6  ~ "pos",
    r <= -0.6 ~ "neg"
  ))%>%
  drop_na(direct)  #
plotdat<-corr_res

ggplot(plotdat, aes(x = name_x, y = name_y, fill = r)) +
  geom_tile(color = "white") +
  # geom_text(aes(label = round(r, 2)), size = 2, color = "black") +
  scale_fill_gradient2(low = "#8AB1D2", high = "#ED8828", mid = "white", name = "Correlation") +
  theme_minimal() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(size = 10),strip.text = element_text(size = 10, color = "black", face = "italic"),
        panel.background = element_rect(linewidth = .4, color = "black", fill = "transparent"),
        panel.grid = element_line(linewidth = .5, color = "grey85")) +
  labs(title = "Correlation Heatmap", x = "Variables", y = "Variables") +
  coord_fixed(ratio = 0.8)
#ggsave("corr_abundance_heatmap.pdf", width = 10, height = 6)


#### div cor ####
pacman::p_load(tidyr, dplyr, tibble, purrr)
library(ggpubr)

source("diversity.R")
profile_args <- read.delim("../ARO_abundance.txt", row.names = 1)
profile_vfgs <- read.delim("../vfdb_tpm.profile", row.names = 1)
profile_mges <- read.delim("../MGEs_abundance.txt", row.names = 1)

dat <- list(
  args = calu_alpha(profile_args, method = "shannon", out_colnames = "args"),
  vfgs = calu_alpha(profile_vfgs, method = "shannon", out_colnames = "vfgs"),
  mges = calu_alpha(profile_mges, method = "shannon", out_colnames = "mges")
) %>% reduce(\(x, y) merge(x, y, by = "sample"))

grp_pair <- list(c("args", "mges"),c("args", "vfgs"), c("vfgs", "mges"))
p_list <- list()

# 
for (x in 1:length(grp_pair)) {
  p_list[[x]] <- ggscatter(dat, x = grp_pair[[x]][1], y = grp_pair[[x]][2], rug = T, add = "reg.line",
                           conf.int = T, conf.int.level = .95, title = "shannon",
                           xlab = grp_pair[[x]][1], ylab = grp_pair[[x]][2], size = .8,
                           add.params = list(color = "#F59B7B", fill = "lightgray", size = 1)) +
    stat_cor(label.sep = "\n", color = "black", method = "spearman", size = 3) +
    theme(aspect.ratio = 1)
}
cowplot::plot_grid(plotlist = p_list, nrow  = 1)
#ggsave("div_arg_vfg_mge_richness.pdf", width = 6, height = 3)


###########
library(vegan)
profile_args <- read.delim("../ARO_abundance.txt", row.names = 1)
profile_mges <- read.delim("../MGEs_abundance.txt", row.names = 1)

#
profile_mges_clean <- profile_mges[, colSums(profile_mges != 0) > 0]
profile_args_clean <- profile_args[, colSums(profile_args != 0) > 0]

#
dist_1 <- vegdist(t(profile_args_clean))
dist_2 <- vegdist(t(profile_mges_clean))

#
PCoA_1 <- cmdscale(dist_1)
PCoA_2 <- cmdscale(dist_2)

###
common_samples <- intersect(rownames(PCoA_1), rownames(PCoA_2))
PCoA_1 <- PCoA_1[common_samples, ]
PCoA_2 <- PCoA_2[common_samples, ]

#
proc <- procrustes(PCoA_1, PCoA_2, symmetric = T)
summary(proc)

plot(proc, kind = 2)
residuals(proc)

set.seed(2025)
proc_test <- protest(PCoA_1, PCoA_2, permutations = 999)
proc_test

proc_test$ss

proc_test$signif


proc_point <- cbind(
  data.frame(proc$Yrot) %>% rename(X1_rotated = X1, X2_rotated = X2),
  data.frame(proc$X) %>% rename(X1_target = Dim1, X2_target = Dim2) 
)
proc_coord <- data.frame(proc$rotation) 

# 
ggplot(proc_point) + #
  geom_segment(aes(x = X1_rotated, y = X2_rotated, xend = (X1_rotated + X1_target)/2, yend = (X2_rotated + X2_target)/2), 
               arrow = arrow(length = unit(0, 'cm')), color = "#8AB1D2", size = .4) +
  geom_segment(aes(x = (X1_rotated + X1_target)/2, y = (X2_rotated + X2_target)/2, xend = X1_target, yend = X2_target),
               arrow = arrow(length = unit(0.2, 'cm')), color = "#F59B7B", size = .4) +
  geom_point(aes(X1_rotated, X2_rotated), color = "#8AB1D2", size = 1.6, shape = 16) + 
  geom_point(aes(X1_target, X2_target), color = "#F59B7B", size = 1.6, shape = 16) + 
  labs(x = 'Dim 1', y = 'Dim 2',
       subtitle = paste0("coefficients: M2 = ", round(proc_test$ss, 4), ", p = ", proc_test$signif)) +
  labs(title = "Correlation analysis by Procrustes analysis") +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.4) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.4) +
  geom_abline(intercept = 0, slope = proc_coord[1,2]/proc_coord[1,1], size = 0.4) +
  geom_abline(intercept = 0, slope = proc_coord[2,2]/proc_coord[2,1], size = 0.4) +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.line = element_blank(),
        plot.title = element_text(size = 10, color = "black"),
        plot.subtitle = element_text(size = 10, color = "black"),
        panel.border = element_rect(linewidth = .4, color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        aspect.ratio = 3/4)


################################################## vfg##############################################################################################

library(vegan)
profile_args <- read.delim("../ARO_abundance.txt", row.names = 1)

profile_vfgs <- read.delim("../vfdb_tpm.profile",sep = "\t",row.names = 1) 

#
profile_vfgs_clean <- profile_vfgs[, colSums(profile_vfgs != 0) > 0]
profile_args_clean <- profile_args[, colSums(profile_args != 0) > 0]


dist_1 <- vegdist(t(profile_vfgs_clean), na.rm = TRUE)
dist_2 <- vegdist(t(profile_args_clean), na.rm = TRUE)


PCoA_1 <- cmdscale(dist_1)
PCoA_2 <- cmdscale(dist_2)


common_samples <- intersect(rownames(PCoA_1), rownames(PCoA_2))
PCoA_1 <- PCoA_1[common_samples, ]
PCoA_2 <- PCoA_2[common_samples, ]


proc <- procrustes(PCoA_1, PCoA_2, symmetric = T)
summary(proc)


plot(proc, kind = 2)
residuals(proc)

set.seed(2025)
proc_test <- protest(PCoA_1, PCoA_2, permutations = 999)
proc_test

proc_test$ss

proc_test$signif

proc_point <- cbind(
  data.frame(proc$Yrot) %>% rename(X1_rotated = X1, X2_rotated = X2),
  data.frame(proc$X) %>% rename(X1_target = Dim1, X2_target = Dim2) 
proc_coord <- data.frame(proc$rotation) #

# 
ggplot(proc_point) + 
  geom_segment(aes(x = X1_rotated, y = X2_rotated, xend = (X1_rotated + X1_target)/2, yend = (X2_rotated + X2_target)/2), 
               arrow = arrow(length = unit(0, 'cm')), color = "#8AB1D2", size = .4) +
  geom_segment(aes(x = (X1_rotated + X1_target)/2, y = (X2_rotated + X2_target)/2, xend = X1_target, yend = X2_target),
               arrow = arrow(length = unit(0.2, 'cm')), color = "#F59B7B", size = .4) +
  geom_point(aes(X1_rotated, X2_rotated), color = "#8AB1D2", size = 1.6, shape = 16) + 
  geom_point(aes(X1_target, X2_target), color = "#F59B7B", size = 1.6, shape = 16) + 
  labs(x = 'Dim 1', y = 'Dim 2',
       subtitle = paste0("coefficients: M2 = ", round(proc_test$ss, 4), ", p = ", proc_test$signif)) +
  labs(title = "Correlation analysis by Procrustes analysis") +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.4) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.4) +
  geom_abline(intercept = 0, slope = proc_coord[1,2]/proc_coord[1,1], size = 0.4) +
  geom_abline(intercept = 0, slope = proc_coord[2,2]/proc_coord[2,1], size = 0.4) +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.line = element_blank(),
        plot.title = element_text(size = 10, color = "black"),
        plot.subtitle = element_text(size = 10, color = "black"),
        panel.border = element_rect(linewidth = .4, color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        aspect.ratio = 3/4)
#ggsave("procrustes_vfgs_taxo.pdf", width = 6, height = 4.5)

#
source("profile_process.R")
metadata_args <- read.delim("../ARGs_info_rename.txt") 
metadata_vfgs <- read.delim("../all_info",sep = "\t",header = T)

#选取耐药流行率前xxx
tem <- read.delim("../ARGs_prevalence.txt") %>% 
  .[order(.$avg_tpm, decreasing = TRUE),] %>% head(20) 

profile_args <- read.delim("../ARO_abundance.txt", row.names = 1) %>%
  mutate(ARG=metadata_args$ARG[match(rownames(.),metadata_args$ARO)])%>% 
  group_by(ARG) %>% summarise_all(sum) %>% column_to_rownames(var = "ARG") %>% 
  mutate(tem_value=tem$ARG[match(rownames(.),tem$ARG)]) %>% filter(!is.na(tem_value)) %>% 
  select(-tem_value) %>% t()%>% data.frame(check.names = F)

tem2 <- read.delim("../VFDB_prevalence.txt") %>% 
  .[order(.$avg_tpm, decreasing = TRUE),] %>% head(20) 


profile_vfgs <- read.delim("../VFs_abundance.txt",sep = "\t",row.names = 1) %>% 
  t()%>% data.frame(check.names = F) %>% 
  select(all_of(tem2$name))

####install.packages("psych")###########
source("corr_process.R")
corr <- psych::corr.test(x = profile_args, y = profile_vfgs, method = "spearman", adjust = "BH")

corr_res <- corr_process(corr, rho = 0, padj = 0.05, out_df = T) %>%
  mutate(direct = case_when(
    r >= 0.6  ~ "pos",
    r <= -0.6 ~ "neg"
  ))%>%
  drop_na(direct)  

#write.table(corr_res, "corr_res.tsv", sep = "\t", row.names = F, quote = F)

#####
pacman::p_load(tidygraph,ggraph,igraph)
plotdat<- corr_res 
edges <- plotdat %>%rename(from = name_x, to = name_y) %>% 
  mutate(abs=abs(r),dir = factor(direct, levels = c("pos", "neg")))

nodes <-rbind(data.frame(name = edges$from) %>% 
                add_column(type = "ARG"),
              data.frame(name = edges$to) %>% 
                add_column(type = "VF")) %>% 
  mutate(type = factor(type, labels = c("ARG", "VF")))

graph <- tbl_graph(nodes = nodes, edges = edges)
nodes$degree <- degree(graph)   
nodes2 <-nodes[nodes$degree!=0,]
graph <- tbl_graph(nodes = nodes2, edges = edges)

class_color <- set.seed(2025)

ggraph(graph, layout = "circle") + 
  geom_edge_link(aes(edge_width = abs, edge_linetype = dir), 
                 color = "#F59B7B", show.legend = T) + 
  scale_edge_width(range = c(.1, 1)) +
  scale_edge_alpha_manual(values = c(0.2, 0.4)) +
  geom_node_point(aes(fill = type, size = degree), 
                  shape = 21, stroke = .1) +
  scale_size(range = c(2, 6)) +
  scale_fill_manual(values = c("ARG" = "#ED8828", "VF" = "#FFE4B5")) + 
  geom_node_text(aes(label = name),size = 1.5, fontface = "plain") +
  theme(panel.background = element_rect(fill = NA, color = "black"),
        legend.key = element_rect(fill = NA),
        aspect.ratio = 1)
#ggsave("corr_cooccurence_50.pdf", width = 10, height = 6)

######################## mantel test for taxa and ARGs phenotype ################################################################
#remotes::install_github("Hy4m/linkET")
library(linkET)
source("taxa.R")

metadata_args <- read.delim("../ARGs_info_rename.txt") 

profile_args <- read.delim("../ARO_abundance.txt", row.names = 1) %>%
  mutate(ARG=metadata_args$Drug_Rename[match(rownames(.),metadata_args$ARO)])%>% 
  group_by(ARG) %>% summarise_all(sum) %>% column_to_rownames(var = "ARG") %>% 
  t()%>% data.frame(check.names = F)


taxonomy <- read.delim("../all_tax_split",header = F) %>% 
  mutate_all(~replace(., . == "", "Unknown"))
colnames(taxonomy) <- c("feature", "kingdom", "phylum", "class", "order", "family", "genus", "species")


profile_taxa <- read.delim("mag_catalog_tpm.profile", row.names = 1) %>% 
  t() %>% 
  as.data.frame()


tmp_df <- data.frame(taxa = taxonomy$family[match(colnames(profile_taxa), taxonomy$feature)],
                     ncol = 1:ncol(profile_taxa)) %>% arrange(taxa)

tmp_list <- list()
for (i in unique(tmp_df$taxa)) {
  tmp_list[[i]] = c(tmp_df$ncol[tmp_df$taxa %in% i])
}

#
mantel <- mantel_test(profile_taxa, profile_args, spec_select = tmp_list, permutations = 999,
                      method = "spearman", 
                      parallel =  16)
#write.table(mantel, "mantel_test_args_family", sep = "\t", row.names = F, quote = F)
#######################################
#
mag_means <- colMeans(profile_taxa, na.rm = TRUE)

#
mag_to_family <- taxonomy$family[match(colnames(profile_taxa), taxonomy$feature)]

#
family_means <- tapply(mag_means, mag_to_family, mean, na.rm = TRUE)

#
top_families <- names(sort(family_means, decreasing = TRUE))
top_families <- top_families[top_families != "Unknown"]

#
top20_families <- top_families[1:20]

#
selected_mags <- colnames(profile_taxa)[mag_to_family %in% top20_families]

# 
profile_taxa_filtered <- profile_taxa[, selected_mags, drop = FALSE]

# 
filtered_mag_to_family <- mag_to_family[match(selected_mags, colnames(profile_taxa))]

# 
tmp_df_filtered <- data.frame(taxa = filtered_mag_to_family, ncol = 1:length(selected_mags))

#
tmp_list <- list()
for (i in unique(tmp_df_filtered$taxa)) {
  tmp_list[[i]] <- tmp_df_filtered$ncol[tmp_df_filtered$taxa == i]
}

#
mantel <- mantel_test(profile_taxa_filtered, profile_args, spec_select = tmp_list, permutations = 999, 
                      method = "spearman",  parallel = 16)

#write.table(mantel, "mantel_test_args_family_top20.tsv", sep = "\t", row.names = F, quote = F)


mantel_df <- mutate(mantel, 
                    rd = cut(r, breaks = c(-Inf, 0.3, 0.5, 0.7, Inf), labels = c("< 0.3", "0.3 - 0.6", "0.6 - 0.9", ">= 0.9")),
                    pd = cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("< 0.001", "0.001 - 0.01", "0.01 - 0.05", ">= 0.05"))) %>% 
  filter(r > 0.3 & p < 0.05)  

qcorrplot(correlate(profile_args), type = "upper", diag = FALSE) +
  geom_square() +
  geom_couple(data = mantel_df, aes(colour = pd, size = rd), curvature = nice_curvature()) +
  scale_fill_gradient2(
    low = "#FFE4B5",  
    high = "#E74C3C", 
    limits = c(0, 1) 
  ) +
  scale_size_manual(values = c(0.6, 0.9, 1.2)) +
  scale_colour_manual(values = c( "#ED8828","#F59B7B", "#FFE4B5")) +
  guides(
    size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
    colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
    fill = guide_colorbar(title = "Pearson's r", order = 3)
  )

###############################################################################################################################################################
target_families <- c("Enterobacteriaceae")

#
taxonomy_filtered <- taxonomy %>%
  filter(family %in% target_families & species != "Unknown")

#
selected_mags <- taxonomy_filtered$feature
profile_taxa_species <- profile_taxa[, colnames(profile_taxa) %in% selected_mags]

#
species_to_mag <- taxonomy_filtered$species
names(species_to_mag) <- taxonomy_filtered$feature

#
tmp_df_species <- data.frame(
  taxa = species_to_mag[match(colnames(profile_taxa_species), names(species_to_mag))],
  ncol = 1:ncol(profile_taxa_species)
) %>% filter(!is.na(taxa)) %>% arrange(taxa)

tmp_list_species <- list()
for (i in unique(tmp_df_species$taxa)) {
  tmp_list_species[[i]] = tmp_df_species$ncol[tmp_df_species$taxa %in% i]
}

##
#mantel_species <- readRDS("mantel_species_results.rds")

mantel_species <- mantel_test(
  profile_taxa_species, profile_args, 
  spec_select = tmp_list_species, 
  permutations = 999, 
  method = "spearman",
  parallel = 16
)

#write.table(mantel_species, "mantel_specise_famliy_spearman.tsv", sep = "\t", row.names = F, quote = F)

mantel_df <- mutate(mantel_species, 
                    rd = cut(r, breaks = c(-Inf, 0.3, 0.5, 0.7, Inf), labels = c("< 0.3", "0.3 - 0.6", "0.6 - 0.9", ">= 0.9")),
                    pd = cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("< 0.001", "0.001 - 0.01", "0.01 - 0.05", ">= 0.05"))) %>% 
  filter(r > 0.6 & p < 0.05)   %>% 
  mutate(spec = forcats::fct_lump_min(spec, 2)) %>% filter(spec != "Other") 

qcorrplot(correlate(profile_args), type = "upper", diag = FALSE) +
  geom_square() +
  geom_couple(data = mantel_df, aes(colour = pd, size = rd), curvature = nice_curvature()) +
  scale_fill_gradient2(
    low = "#FFE4B5",  
    high = "#E74C3C",
    limits = c(0, 1) 
  ) +
  scale_size_manual(values = c(0.6, 0.9, 1.2)) +
  scale_colour_manual(values = c( "#ED8828","#F59B7B", "#FFE4B5")) +
  guides(
    size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
    colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
    fill = guide_colorbar(title = "Pearson's r", order = 3)
  )
