#### info ###
# encoding: utf-8
# author: Jinxin Meng
# e-mail: mengjx855@163.com
# created data：2022-9-18
# modified date: 2023-12-15
# version: 0.1

library(dplyr)
library(tibble)
library(tidyr)
library(stringr)

#### taxa_split ####
# 拆分物种信息表
# taxonomy: a data.frame containning two field, feature|taxonomy..., also rename by taxonomy_name.
# sep: the separate label.
# na_fill: unclassified rename.
# rm_suffix: remove taxa name with suffix providing by GTDB.
# taxa_level = "k:s", 指定层级列表
taxa_split <- function(taxonomy, taxonomy_colnames = NULL, sep = ";", taxa_level = "k:s", 
                       na_fill = "Unclassified", rm_suffix = F) {
  if (!all(c("feature", "taxonomy") %in% colnames(taxonomy)) & is.null(taxonomy_colnames)) stop("taxonomy field (feature|taxonomy)")
  if (!is.null(taxonomy_colnames)) taxonomy <- data.frame(taxonomy, check.names = F) %>% dplyr::rename(all_of(taxonomy_colnames))
  # 拆分物种信息
  dat <- data.frame(taxonomy, check.names = F)
  # 判断层级
  taxa_level_vec <- c("k", "p", "c", "o", "f", "g", "s", "t")
  taxa_name_vec <- c("kingdom", "phylum", "class", "order", "family", "genus", "species", "strain")
  cut_off <- c(match(unlist(strsplit(taxa_level, ":"))[1], taxa_level_vec):match(unlist(strsplit(taxa_level, ":"))[2], taxa_level_vec))
  taxa_name_vec <- taxa_name_vec[cut_off]
  taxa_level_vec <- taxa_level_vec[cut_off] %>% paste0(., "__")
  # 处理
  res <- map_dfr(dat$taxonomy, \(x) 
      if(is.na(x)){ purrr::map_vec(taxa_level_vec, \(i) ifelse(grepl("__$", i), paste0(i, na_fill), i)) %>% t() %>% data.frame() } 
      else { purrr::map_vec(unlist(strsplit(x, sep)), \(i) ifelse(grepl("__$", i), paste0(i, na_fill), i))[cut_off] %>% t() %>% data.frame() })
  colnames(res) <- taxa_name_vec
  res <- add_column(res, feature = dat$feature, .before = 1)
  # 去除门的子群标记，例如，p__Firmicutes_A, 修改为p__Firmicutes
  if (isTRUE(rm_suffix)) {
    res <- apply(res, 2, \(x) gsub("_\\w$", "", x = x) %>% gsub("_\\w ", " ", x = ., perl = T)) %>% data.frame(check.names = F)
  }
  return(res)
}

#### taxa_trans ####
# 转换物种表的层级关系
# otu: input a otu table
# taxonomy: input a data.frame, colnames following: feature|phylum|family..., also rename by taxonomy_name parameter.
# taxonomy_name = a vector used to rename taxon table.
# to: 指定转变为哪个级别的
# top_n: select top n feature, including Other.
# top_list: select top n feature by specify feature.
# other_name: a name to rename other feature.
# out_all: output all result.
# na_fill: if some feature a not corresponding front level, rename it.
# smp2grp: merge sample to group.
# group: input a data.frame, colnames following: sample|group, also rename by group_name parameter.
# method: merge method, normal by mean or median.
# reture a otu table.
taxa_trans <- function(otu, taxonomy, group, to = "family", top_n = 12, top_list = NULL, 
                       other_name = "Other", out_all = F, na_fill = "Unclassified",
                       transRA = F, smp2grp = F, method = "mean", 
                       taxonomy_colnames = NULL, group_colnames = NULL) {
  otu <- data.frame(otu, check.names = F)
  taxonomy <- data.frame(taxonomy, check.names = F)
  if (!all(c("feature", to) %in% colnames(taxonomy)) & is.null(taxonomy_colnames)) {
    stop(paste0("taxonomy field must have \"feature\" and \"", to, "\", and other level optional input."))
  }
  if (!is.null(taxonomy_colnames)) {
    taxonomy <- dplyr::rename(taxonomy, all_of(taxonomy_colnames))
  }
  if(isTRUE(out_all)) {
    top_n <- 0
  }
  dat <- otu %>% 
    mutate(taxa = taxonomy[match(rownames(.), taxonomy$feature), to],
           taxa = ifelse(is.na(taxa), na_fill, taxa)) %>% 
    group_by(taxa) %>% summarise_all(sum) %>% ungroup() %>% 
    column_to_rownames(var = "taxa") %>% data.frame(check.names = F)
  
  # trans taxon
  if(top_n > 0 & is.null(top_list)) {
    taxa_vec <- data.frame(val = rowSums(dat)) %>%
      arrange(desc(val)) %>% head(n = top_n - 1) %>% 
      rownames(.) %>% unlist() %>% as.character()
    dat <- dat %>% rownames_to_column(var = "feature") %>% 
      mutate(feature = ifelse(feature %in% taxa_vec, feature, other_name)) %>% 
      group_by(feature) %>% summarise_all(sum) %>% ungroup() %>% 
      data.frame(check.names = F) %>% column_to_rownames(var = "feature")
    } else if(!is.null(top_list)) {
      dat <- dat %>% rownames_to_column(var = "feature") %>% 
        mutate(feature = ifelse(feature %in% top_list, feature, other_name)) %>% 
        group_by(feature) %>% summarise_all(sum) %>% ungroup() %>% 
        data.frame(check.names = F) %>% column_to_rownames(var = "feature")
    }
  
  # sample to group
  if(isTRUE(smp2grp) & !missing(group)) {
    group <- data.frame(group, check.names = F)
    if (!all(c("sample", "group") %in% colnames(group)) & is.null(group_colnames)) stop("group field must have \"sample\" and \"group\"")
    if (!is.null(group_colnames)) group <- dplyr::rename(group, all_of(group_colnames))
    dat <- data.frame(t(dat), check.names = F) %>% 
      mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
      group_by(group) %>% summarise_all(method) %>% ungroup() %>% 
      column_to_rownames(var = "group") %>% t() %>% data.frame(check.names = F)
  }
  
  # transRA
  if(isTRUE(transRA)) {
    dat <- apply(dat, 2, \(x) x/sum(x)*100) %>% data.frame(check.names = F)
  }
  return(dat)
}

#### plot_compos ####
# 物种组成,
# otu表输入行名为feature，列名为样本的表
# group输入列名为sample和group的表
# taxonomy输入taxa.split()函数输出的表格，或者是指定新表格，第一列为feature，随后是不同层级的分类信息
# display指定样本级别的组成还是分组级别的组成， 必须指定为sample或者group
# taxa_level显示哪个水平组成, 可选phylum, class, family, genus, species
# top_n输入显示多少个物种，不足数量的按最大数量算
# top_list指定物种列表
# group_order/sample_order，分别之指定sample和group的顺序
# taxa_order指定物种的排序
# taxa_color指定物种的颜色
# title指定图的标题内容
plot_compos <- function(otu, taxonomy, group, display = "group", taxa_level = "family", 
                        top_n = 12, top_list = NULL, group_order = NULL, sample_order = NULL,
                        taxa_order = NULL, taxa_color = NULL, plot_title = NULL,
                        taxonomy_colnames = NULL, group_colnames = NULL){
  if(missing(otu) & missing(taxonomy)) { 
    stop("missing parameter.") 
  }
  if(display == "group" & !missing(group)) {
    group <- data.frame(group, check.names = F) 
  } else if (display == "group" & missing(group)) { 
    stop("missing parameter.") 
  }
  # 读入数据
  otu <- data.frame(otu, check.names = F)
  taxonomy <- data.frame(taxonomy, check.names = F)
  colors <- c("#4E79A7FF","#A0CBE8FF","#F28E2BFF","#FFBE7DFF","#59A14FFF",
              "#8CD17DFF","#B6992DFF","#F1CE63FF","#499894FF","#86BCB6FF",
              "#E15759FF","#FF9D9AFF","#79706EFF","#BAB0ACFF","#D37295FF",
              "#FABFD2FF","#B07AA1FF","#D4A6C8FF","#9D7660FF","#D7B5A6FF")
  if(display == "group") {
    other_name <- paste0(str_to_lower(str_sub(taxa_level, 1, 1)), "__Other")
    dat <- taxa_trans(otu, taxonomy, group, to = taxa_level, top_n = top_n,
                      top_list = top_list, other_name = other_name, smp2grp = T, transRA = T,
                      group_colnames = group_colnames, taxonomy_colnames = taxonomy_colnames)
    if(is.null(taxa_order)) { 
      taxa_order <- names(sort(rowSums(dat), decreasing = T))
    }
    if(is.null(taxa_color)) {
      taxa_count <- nrow(dat)
      taxa_color <- rep(colors, time = ceiling(taxa_count/20))[1:taxa_count]
    }
    if(is.null(group_order)) {
      group_order <- names(sort(unlist(dat[taxa_order[1],])))
    }
    if(is.null(plot_title)) {
      plot_title = stringr::str_to_sentence(taxa_level)
    }
    fill_title = paste0(plot_title, " taxa")
    plotdat <- rownames_to_column(dat, var = "feature") %>%  
      gather(., key = "group", value = "val", -feature) %>% 
      mutate(group = factor(group, group_order), 
             feature = factor(feature, taxa_order))
    p <- ggplot(plotdat, aes(group, val, fill = feature)) +
      geom_bar(stat = "identity", position = position_stack(), color = "#000000", linewidth = .2, width = .8) +
      scale_fill_manual(values = taxa_color) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = "", y = "Relative Abundance (%)", title = plot_title, fill = fill_title) +
      theme_classic() +
      theme(axis.line = element_line(linewidth = .4, color = "#000000"),
            axis.ticks = element_line(linewidth = .4, color = "#000000"),
            axis.text = element_text(size = 8, color = "#000000"),
            axis.title = element_text(size = 8, color = "#000000"),
            plot.title = element_text(size = 10, color = "#000000"),
            legend.text = element_text(size = 8, color = "#000000", face = "italic"),
            legend.title = element_text(size = 10, color = "#000000"),
            panel.grid = element_blank())
    message("  ggsave(file = \"compos_stacked_barplot.pdf\", width = 6, height = 4)")
  } else if(display == "sample") {
    other_name <- paste0(str_to_lower(str_sub(taxa_level, 1, 1)), "__Other")
    dat <- taxa_trans(otu, taxonomy, group, to = taxa_level, top_n = top_n,
                      top_list = top_list, other_name = other_name, smp2grp = F, transRA = T,
                      group_colnames = group_colnames, taxonomy_colnames = taxonomy_colnames)
    if(is.null(taxa_order)) { 
      taxa_order <- names(sort(rowSums(dat), decreasing = T))
    }
    if(is.null(taxa_color)) {
      taxa_color <- rep(colors, time = ceiling(nrow(dat)/20))[1:nrow(dat)]
    }
    if(is.null(sample_order)) {
      sample_order <- names(sort(unlist(dat[taxa_order[1],])))
    }
    if(is.null(plot_title)) {
      plot_title = stringr::str_to_sentence(taxa_level)
    }
    fill_title = paste0(plot_title, " taxa")
    plotdat <- rownames_to_column(dat, var = "feature") %>%  
      gather(., key = "sample", value = "val", -feature) %>% 
      mutate(sample = factor(sample, sample_order),
             feature = factor(feature, taxa_order))
    p <- ggplot(plotdat, aes(sample, val, fill = feature)) +
      geom_bar(stat = "identity", position = position_stack(), color = "#000000", linewidth = .2, width = 1) +
      scale_fill_manual(values = taxa_color) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = "", y = "Relative Abundance (%)", title = plot_title, fill = fill_title) +
      theme_classic() +
      theme(axis.line = element_line(linewidth = .4, color = "#000000"),
            axis.ticks = element_line(linewidth = .4, color = "#000000"),
            axis.ticks.x = element_blank(),
            axis.text = element_text(size = 8, color = "#000000"),
            axis.text.x = element_blank(),
            axis.title = element_text(size = 8, color = "#000000"),
            plot.title = element_text(size = 10, color = "#000000"),
            legend.text = element_text(size = 8, color = "#000000", face = "italic"),
            legend.title = element_text(size = 10, color = "#000000"),
            panel.grid = element_blank())
    message("  ggsave(file = \"compos_stacked_barplot.pdf\", width = 8, height = 4)")
  } else {
    stop("error parameter.")
  }
  return(p)
}

#### plot_mcompos ####
# 物种组成,
# otu表输入行名为feature，列名为样本的表
# group输入列名为sample和group的表
# taxonomy输入taxa.split()函数输出的表格，或者是指定新表格，第一列为feature，随后是不同层级的分类信息
# display指定样本级别的组成还是分组级别的组成， 必须指定为sample或者group
# top_n输入显示多少个分类，不足数量的按最大数量算
plot_mcompos <- function(otu, taxonomy, group, display = "group", top_n = 12, 
                         top_list = NULL, group_order = NULL, sample_order = NULL,
                         taxonomy_colnames = NULL, group_colnames = NULL) {
  if(missing(otu) & missing(taxonomy)) {
    stop("missing parameter.")
  }
  taxa_level <- setdiff(colnames(taxonomy), c("feature", "kingdom"))
  
  p_list <- list()
  for (i in taxa_level) {
    p <- plot_compos(otu = otu, taxonomy = taxonomy, group = group, display = display,
                     top_n = top_n, top_list = top_list, taxa_level = i,
                     group_order = group_order,sample_order = sample_order, 
                     taxonomy_colnames = taxonomy_colnames, group_colnames = group_colnames) %>% 
      suppressMessages()
    p_list[[i]] <- p
  }
  p <- cowplot::plot_grid(plotlist = p_list, nrow = 2, align = "v")
  if (display == "group") {
    cat("  ggsave(file = \"compos_stacked_barplot.pdf\", width = 18, height = 8)")
  } else if (display == "sample") {
    cat("  ggsave(file = \"compos_stacked_barplot.pdf\", width = 30, height = 8)")
  }
  return(p)
}

#### plot_taxa_diff ####
# taxa的差异
# 必须对接taxa.compos.tib(display = "sample")
# group_list指定的分组的因子顺序
# group_color指定分组的填充颜色
# 调用diff.test()函数计算组间差异，可选的方法有wilcox和t检验
# 输出ggplot对象的list
plot_taxa_box <- function(otu, group, group_order = NULL, 
                           group_color = NULL, method = "wilcox"){
  otu <- data.frame(otu, check.names = F)
  source("D:/桌面/R plot/utilities.R")
  group <- data.frame(group, check.names = F) %>% 
    dplyr::select(sample, group)
  dat <- data.frame(t(otu), check.names = F) %>% 
    rownames_to_column(var = "sample") %>% 
    dplyr::filter(sample %in% group$sample) %>% 
    merge(., group, by = "sample", all.x = T)
  
  # group的order
  if(is.null(group_order)){
    group_order <- unique(dat$group)
  }
  # group的配色
  if(is.null(group_color)){
    color <- c("#4E79A7FF","#A0CBE8FF","#F28E2BFF","#FFBE7DFF","#59A14FFF",
               "#8CD17DFF","#B6992DFF","#F1CE63FF","#499894FF","#86BCB6FF",
               "#E15759FF","#FF9D9AFF","#79706EFF","#BAB0ACFF","#D37295FF",
               "#FABFD2FF","#B07AA1FF","#D4A6C8FF","#9D7660FF","#D7B5A6FF")
    group_color <- rep(color, times = ceiling(length(group_order)/12))[1:length(group_order)]
  }
  
  feature <- setdiff(colnames(dat), c("sample", "group"))
  
  p_list <- list()
  for (i in feature) {
    dat_i <- dat %>% 
      select(sample, group, all_of(i)) %>% 
      rename(index = all_of(i))
      
    source("D:/桌面/R plot/diff_test.R")
    comparison <- diff_test(select(dat_i, -group), select(dat_i, -index), 
                            group_order = group_order, method = method) %>% 
      filter(pval < 0.05) %>% 
      select(group_pair) %>% 
      unlist() %>% 
      as.character() %>% 
      strsplit(x = ., split = "_vs_")
    
    p <- ggplot(dat_i, aes(x = factor(group, group_order), y = index, 
                           color = factor(group, group_order))) + 
      geom_boxplot(fill = "transparent", outlier.size = .7, lwd = .4, show.legend = F) +
      geom_jitter(size = 1, width = .3) +
      scale_color_manual(values = group_color) +
      labs(x = "", y = "Relative Abundance(%)", color = "Group", title = i) +
      ggpubr::stat_compare_means(comparisons = comparison, method = method, size = 3,
                                 method.args = list(exact = F), label = "p.signif", 
                                 tip.length = .01, step.increase = .03, vjust = .9) +
      theme_classic() +
      theme(axis.line = element_line(linewidth = .4, color = "#000000"),
            axis.ticks = element_line(linewidth = .4, color = "#000000"),
            axis.text = element_text(size = 8, color = "#000000"),
            axis.title = element_text(size = 8, color = "#000000"),
            plot.title = element_text(size = 10, color = "#000000", hjust = .5, face = "italic"),
            legend.text = element_text(size = 8, color = "#000000"),
            legend.title = element_text(size = 10, color = "#000000"),
            panel.grid = element_blank())

    p_list[[i]] <- p
  }
  return(p_list)
}
# lapply(p, function(x) ggsave(x, filename = paste0("family_", as.character(x$labels[4]), ".pdf"), width = 4.5, height = 4.5))
