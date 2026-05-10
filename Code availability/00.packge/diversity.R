#### info ####
# encoding: utf-8
# author: Jinxin Meng
# e-mail: mengjx855@163.com
# created data：2021-10-29
# modified data: 2023-12-04
# version: 0.2

# 2023-01-01: update function: calu_alpha.
# 2023-12-04: update function: check_file_name was deprecated.

library(dplyr)
library(tibble)
library(tidyr)
library(vegan)
library(ggplot2)
source("D:/桌面/R plot/diff_test.R")

#### calu_alpha ####
# export data.frame with colnames (sample|val)
calu_alpha <- function(profile, method = 'richness', out_colnames = NULL, tree = NULL) {
  otu_table = t(profile)
  if (method == 'richness') {
    result <- rowSums(otu_table > 0)
  } else if (method == 'chao1') {
    result <- estimateR(ceiling(otu_table))[2,]
  } else if (method == 'observed') {
    result <- estimateR(ceiling(otu_table))[1,]
  } else if (method == 'ace') {
    result <- estimateR(ceiling(otu_table))[4,]
  } else if (method == 'shannon') {
    result <- vegan::diversity(otu_table, index = 'shannon')
  } else if (method == 'simpson') {
    result <- vegan::diversity(otu_table, index = 'simpson')
  } else if (method == 'pielou') {
    result <- vegan::diversity(otu_table, index = 'shannon') / log(estimateR(otu_table)[1, ])
  } else if (method == 'gc') {
    result <- 1 - rowSums(otu_table == 1) / rowSums(otu_table)
  } else if (method == 'pd' & !is.null(tree)) {
    pd <- pd(otu_table, tree, include.root = FALSE)
    result <- pd[, 1]
    names(result) <- rownames(pd)
  }
  index <- data.frame(sample = names(result), val = result, row.names = NULL)
  if (!is.null(out_colnames) & is.character(out_colnames)) colnames(index)[2] <- out_colnames
  return(index)
}

#### plot_alpha ####
# alpha diversity index;
# dat: a data.frame with colnames (sample|val)
# group: a data.frame with colnames (sample|group)
# dat_colnames|group_colnames: a vector for current colnames of the data.frame with new name.
# for example: group_colnames = c(sample = "Sample", val = "Val")
# group_order: group order
# group_color: group color
# diff_test()
plot_alpha <- function(dat, group, group_order = NULL, group_color = NULL, dat_colnames = NULL, group_colnames = NULL, 
                       xlab = "", ylab = "", title = "", aspect_ratio = 1, diff_test = T, method = "wilcox") {
  if (!all(colnames(dat) %in% c("sample", "val")) & is.null(dat_colnames)) stop("dat field (sample|val)")
  if (!all(colnames(group) %in% c("sample", "group")) & is.null(group_colnames)) stop("group field (sample|group)")
  if (!is.null(dat_colnames)) dat <- data.frame(dat, check.names = F) %>% dplyr::rename(all_of(dat_colnames))
  if (!is.null(group_colnames)) group <- data.frame(group, check.names = F) %>% dplyr::rename(all_of(group_colnames))
  if (is.null(group_order)) group_order <- unique(group$group)
  if (is.null(group_color)) {
    color <- c("#1f78b4","#33a02c","#e31a1c","#ff7f00","#6a3d9a","#ffff99","#b15928","#a6cee3","#b2df8a","#fb9a99","#fdbf6f","#cab2d6")
    group_color <- rep(color, times = ceiling(length(group_order)/length(color)))[1:length(group_order)]
  }
  
  plotdat <- merge(dat, group, by = "sample", all = T) %>% 
    mutate(group = factor(group, levels = group_order))

  p <- ggplot(plotdat, aes(x = group, y = val, fill = group)) + 
    geom_boxplot(position = position_dodge(width = .9), color = "black", 
                 outlier.size = .3, outlier.color = "grey30", lwd = .4, show.legend = F) +
    geom_jitter(size = .3, width = .2, color = "grey30", show.legend = F) +
    scale_fill_manual(values = group_color) +
    labs(x = xlab, y = ylab, title = title) +
    theme_bw() +
    theme(axis.text = element_text(color = "black", size = 8),
          axis.title = element_text(color = "black", size = 10),
          axis.ticks = element_line(color = "black", linewidth = .4),
          axis.line = element_blank(),
          panel.background = element_rect(color = "black", linewidth = .4),
          panel.grid = element_line(color = "grey90", linewidth = .4),
          aspect.ratio = aspect_ratio)
  
  # 是否差异检验
  if (isTRUE(diff_test)) {
    comparisons <- diff_test(dat, group, method = method, group_order = group_order) %>% 
      filter(pval < 0.05) %>% select(group_pair) %>% unlist() %>% as.character() %>% strsplit(x = ., split = "_vs_")
    p <- p + stat_compare_means(comparisons = comparisons, method = method, method.args = list(exact = F), label = "p.signif", 
                                tip.length = .02, step.increase = .05, vjust = .8, size = 3) 
  }
  cat(paste0("  Suggestion: ggsave(file = \"xxx.pdf\", width = ", length(group_order)*0.6 ,", height = 4.5)\n"))
  return(p)
}
