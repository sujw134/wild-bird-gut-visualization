
rm(list = ls())
library("tidyverse")
setwd("06/")
source("rare_curve.R")
profile <- read.delim("../ARO_abundance.txt") %>% 
  column_to_rownames(var = "name") %>% 
  {sweep(., 2, colSums(.), FUN = "/")} 

sample_group <- read.delim("../group_species_8.txt")   
otu <- profile[, sample_group$sample]%>% 
  data.frame()

## ###########################################α-diversity############################################################################################################################
#
source("diversity.R")
method = "shannon" 

dat <-  alpha.index(otu, method = method)
group <- read.delim("../group_species_8.txt",header = T)  
plotdat <- merge(dat, group, by = "sample", all.x = T) %>% 
  select(-sample)

#install.packages("ggrain")
library(ggrain)
library(patchwork)

str(plotdat)

comparison <- diff.test(dat, group,method = "wilcox") %>% 
  filter(pval < 0.05) %>% 
  select(group_pair) %>% 
  unlist() %>% 
  as.character() %>% 
  strsplit(x = ., split = "_vs_")

####
colors <- c("#F59B7B" ,"#ED8828","#FCC41E","#FFD700","#FFE4B5","#F2EFBB", "#8AB1D2","#6BB7CA",
            "#33ABC1","#A4DDD3","#ABD7EC","#b2df8a", "#8D73BA","#C6B3D3", "#33a02c","#80BA8A","#9CD1CB", "#fccde5", "#7f7f7f",  "#ED9F9B","#81B21F" )
ggplot(plotdat, aes(x = group, y = index)) +
  geom_boxplot(width = .2, aes(fill = group), color = "black", outlier.shape = NA, lwd = .3, show.legend = FALSE) +
  geom_jitter(aes(fill = group), width = 0.1, size = 1, alpha = 0.6, shape = 21, stroke = 0.5, color = "black", show.legend = FALSE) +  # 添加散点并设置边框为黑色
  scale_fill_manual(values = colors) +
  labs(x = "", y = paste0(stringr::str_to_title(method), " Index")) +
  stat_compare_means(comparisons = comparison, method = "wilcox", method.args = list(exact = FALSE), label = "p.signif", 
                     tip.length = 0.02, step.increase = 0.05, vjust = 0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  
        axis.text = element_text(color = "black", size = 10),
        axis.ticks = element_line(linewidth = 0.5),
        axis.title = element_text(size = 12),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 1))

##
comparison_df <- diff.test(dat, group, method = "wilcox") %>% 
  mutate( p.adj = p.adjust(pval, method = "BH"), 
          p.signif = case_when( p.adj <= 0.001 ~ "***", p.adj <= 0.01 ~ "**", p.adj <= 0.05 ~ "*", TRUE ~ "ns" ) ) #%>% filter(p.adj < 0.05) # 仅保留显著组间比较

###########################################
summary_data <- plotdat %>%
  group_by(group) %>%
  summarise(
    mean_index = mean(index, na.rm = TRUE),
    se = sd(index, na.rm = TRUE) / sqrt(n())
  ) %>%
  arrange(desc(mean_index)) %>%
  mutate(group = fct_reorder(group, mean_index, .desc = TRUE))

#
gradient_colors <- colorRampPalette(c("#F59B7B", "#FFE9E0"))(nrow(summary_data))

#
ggplot(summary_data, aes(x = group, y = mean_index, fill = group)) +
  geom_col(position = "dodge", color = "black") + 
  geom_errorbar(aes(ymin = mean_index - se, ymax = mean_index + se), 
                position = position_dodge(0.9), width = 0.2) +  
  scale_fill_manual(values = gradient_colors) +  # 
  labs(x = "", y = paste0(stringr::str_to_title(method), " Index")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  
    axis.text = element_text(color = "black", size = 10),
    axis.ticks = element_line(linewidth = 0.5),
    axis.title = element_text(size = 12),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 1),
    legend.position = "none"  # 
  )

### ggsave(paste0("family_", method, ".pdf"), width = 4, height = 4.5)
######################################beta diversity####################################################################################################
source("diff_test.R")
library(vegan)
# PCoA

otu_clean <- otu[, colSums(is.na(otu)) == 0]  

distance <- vegdist(t(otu_clean), method = "bray")
####
#
PCoA <- cmdscale(distance, k = 2, eig = T)

PCoA_points <- data.frame(PCoA$points) %>% # 
  rownames_to_column(var = "sample")
PCoA_eig <- round(PCoA$eig/sum(PCoA$eig) * 100, digits = 2)

plotdat <- merge(PCoA_points, group, by = "sample", all.x = T) # %>% 
 
group <- group[match(rownames(as.matrix(distance)), group$sample),] %>% as.data.frame(row.names = NULL)  # 整体水平比较,保证样本在两个数据集中对应
adonis <- adonis2(distance ~ group, group, permutations = 999)  
label <- paste0("R2 = ",round(adonis[1,3], digits = 4),"  p = ", adonis[1,5])


color_mapping <- c(
  "Cygnus cygnus" =  "#A4DDD3",
  "Tadorna tadorna" = "#ED8828",
  "Anser cygnoides" = "#FCC41E",
  "Anser fabalis" = "#FFD700",
  "Anser indicus" = "#FFE4B5",
  "Cygnus atratus" = "#F2EFBB",
  "Anser albifrons" = "#8AB1D2",
  "Emberiza spodocephala" = "#6BB7CA",
  "Tarsiger cyanurus" = "#33ABC1",
  "Gyps himalayensis" = "#fccde5",
  "Branta canadensis" = "#ABD7EC",
  "Corvus dauuricus" = "#b2df8a",
  "Corvus frugilegus" = "#8D73BA",
  "Lagopus muta" = "#C6B3D3",
  "Larus ridibundus" = "#33a02c",
  "Pyrrhocorax pyrrhocorax" = "#80BA8A",
  "Cygnus olor" = "#F59B7B"
)


plot <- ggplot(plotdat, aes(x = X1, y = X2, fill = group, colour = group)) +
  geom_vline(xintercept = 0, lty = 2, lwd = .5) +
  geom_hline(yintercept = 0, lty = 2, lwd = .5) +
  geom_point(shape = 21, size = 3) +
  stat_ellipse(aes(color = group), geom = 'path', level = .95, show.legend = F, lty = 2) +
  scale_fill_manual(values = color_mapping) +
  scale_color_manual(values = color_mapping, guide = "none") +
  labs(x = paste("PCoA1 (", PCoA_eig[1], "%)", sep = ""), 
       y = paste("PCoA2 (", PCoA_eig[2], "%)", sep = ""),
       title = "bray_curtis PCoA", fill = "Group",
       subtitle = label) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = "black", linewidth = .5),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black", linewidth  = .5),
        axis.ticks = element_line(color = "black", linewidth  = .5),
        aspect.ratio = 3/4,
        legend.justification = c(0, 1),
        legend.position= c(0, 1),
        legend.background = element_rect(color = "black", linewidth  = .5))

###############
library("ggExtra")
ggMarginal(
  plot,
  type = 'histogram',
  margins = 'both',
  size = 8,
  groupFill = T
)
