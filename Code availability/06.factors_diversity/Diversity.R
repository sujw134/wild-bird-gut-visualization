
rm(list = ls())
library("tidyverse")
setwd("06/")
source("rare_curve.R")
profile <- read.delim("../ARO_abundance.txt") %>% 
  column_to_rownames(var = "name") %>% 
  {sweep(., 2, colSums(.), FUN = "/")}  

sample_group <- read.delim("group_residency_new.txt")
###
otu <- profile[, sample_group$sample]%>% 
  data.frame()

##################################################################
##
source("diversity.R")
method = "simpson"

dat <-  alpha.index(otu, method = method)
group <- read.delim("../group_residency_new.txt",header = T)     
plotdat <- merge(dat, group, by = "sample", all.x = T) %>% 
  #mutate(group = factor(group, levels = c("16-17","18-19","20-21", "22-23"))) %>% 
  select(-sample)

#install.packages("ggrain")
library(ggrain)
library(patchwork)

# 1
stat.test <- plotdat %>%
  wilcox_test(index ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_y_position(data = plotdat, formula = index ~ group, step.increase = 0.6)


# 2
ggplot(plotdat, aes(x = group, y = index, fill = group, color = group)) + 
  geom_rain(
    alpha = 0.6,
    violin.args = list(color = NA, alpha = 0.3), 
    point.args = list(size = 1.2, shape = 16), 
    point.args.pos = list(position = position_jitter(width = 0.1)),
    boxplot.args = list(color = "black", outlier.shape = NA, alpha = 0.7) 
  ) +
  stat_pvalue_manual(
    stat.test, 
    label = "p.adj.signif",             
    tip.length = 0.01,                 
    hide.ns = TRUE,                    
    coord.flip = TRUE,                 
    inherit.aes = FALSE,               
    vjust = -0.001 
  ) +
  scale_fill_manual(values = c("#F59B7B" ,"#ED8828","#FCC41E","#FFD700","#FFE4B5","#F2EFBB", "#DAA520")) +
  scale_color_manual(values = c("#F59B7B" ,"#ED8828","#FCC41E","#FFD700","#FFE4B5","#F2EFBB", "#DAA520")) +
  labs(
    x = "", 
    y = paste0(stringr::str_to_title(method), " Index"),
    caption = ""
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(size = 12),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 1),
    legend.position = "none"
  ) +
  expand_limits(y = max(stat.test$y.position) * 1.1) +
  coord_flip()

  select(-groups, -.y.) %>%
  write.table(file = "xxx_wilcox_BH_simpson.txt", sep = "\t", row.names = FALSE, quote = FALSE)

######################################beta diversity####################################################################################################
source("diff_test.R")
library(vegan)
# PCoA
otu_clean <- otu[, colSums(is.na(otu)) == 0]  

distance <- vegdist(t(otu_clean), method = "bray")
####
#distance[is.na(distance)] <- 0
PCoA <- cmdscale(distance, k = 2, eig = T)

PCoA_points <- data.frame(PCoA$points) %>% 
  rownames_to_column(var = "sample")
PCoA_eig <- round(PCoA$eig/sum(PCoA$eig) * 100, digits = 2)

plotdat <- merge(PCoA_points, group, by = "sample", all.x = T) 
#mutate(group = factor(group, levels = c("RB","MB")))

group <- group[match(rownames(as.matrix(distance)), group$sample),] %>% as.data.frame(row.names = NULL)  # 整体水平比较,保证样本在两个数据集中对应
adonis <- adonis2(distance ~ group, group, permutations = 999)  
label <- paste0("R2 = ",round(adonis[1,3], digits = 4),"  p = ", adonis[1,5])

color_mapping <- c(   "Migratory birds" = "#F59B7B",  "Resident birds" = "#8AB1D2")

ggplot(plotdat, aes(x = X1, y = X2, fill = group, colour = group)) +
  geom_vline(xintercept = 0, lty = 2, lwd = .5) +
  geom_hline(yintercept = 0, lty = 2, lwd = .5) +
  geom_point(shape = 21, size = 4) +
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

