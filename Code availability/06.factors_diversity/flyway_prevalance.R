
rm(list = ls())
library("tidyverse")
setwd("06/")
profile <- read.delim("../ARO_abundance.txt") %>% 
  column_to_rownames(var = "name") %>% 
  {sweep(., 2, colSums(.), FUN = "/")}  

sample_group <- read.delim("group_flyway.txt")  

###
otu <- profile[, sample_group$sample]%>% 
  data.frame()

##
otu <- otu[, colSums(is.na(otu)) == 0]

arg_metadat <- read.delim("../ARGs_info_rename.txt")

data <- otu %>%
  mutate(ARG=arg_metadat$ARG[match( rownames(.), arg_metadat$ARO)]) %>% 
rownames_to_column(var = "aro") %>% 
  select(-aro) %>% 
  column_to_rownames(var = "ARG")


#
otu_long <- data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "abundance") %>%
  left_join(sample_group, by = "sample")

#
summary_data <- otu_long %>%
  group_by(group, gene) %>%
  summarise(
    prevalence = sum(abundance > 0) / n(),  
    mean_abundance = mean(abundance, na.rm = TRUE),  
    .groups = "drop"
  )

#
top10_genes <- summary_data %>%
  group_by(gene) %>%
  summarise(total_abundance = mean(mean_abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 10) %>%
  pull(gene)

#
summary_data <- summary_data %>%
  mutate(
    gene_category = if_else(gene %in% top10_genes, gene, "other")
  )

#
custom_colors <- c(
  "#F59B7B", "#ED8828", "#FCC41E", "#FFD700", "#FFE4B5",
  "#F2EFBB", "#8AB1D2", "#6BB7CA", "#33ABC1", "#A4DDD3"
)

# 
colors <- c(
  setNames(custom_colors, top10_genes), 
  "other" = "grey50"
)

library(ggrepel)

ggplot(summary_data, aes(x = prevalence, y = mean_abundance, color = gene_category, label = ifelse(prevalence > 0.75& mean_abundance > 0.2, gene_category, NA))) +
  geom_point(size = 5, alpha = 0.7) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0.75, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_text_repel(
    na.rm = TRUE,
    size = 4,
    fontface = "bold",
    box.padding = 0.3,
    max.overlaps = Inf
  ) +
  facet_wrap(~ group, scales = "free_y") +
  scale_color_manual(values = colors) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Prevalence",
    y = "Mean Abundance",
    color = "ARG Category"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )


write.csv(summary_data,file = "prevalance.csv")
