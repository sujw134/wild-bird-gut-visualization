#### info ####
# encoding: utf-8
# author: Jinxin Meng
# e-mail: mengjx855@163.com
# created data：2022-01-01
# modified date: 2023-12-13
# version: 0.4

# 2022-06-01: 可选择"wilcox rank-sum","one-way anova","student's t test"三种方法做差异分析；
# 2023-01-17: diff_test_profile函数对feature进行差异分析，输入的是标准otu表和group表
# 2023-12-04: 修改diff_test_profile函数中的for循环，使用purrr::map_dfr，速度上稍微快一丢丢。
# 2023-12-04: 修改diff_test函数中的rbind()，使用tibble::add_column()，速度上稍微快一丢丢。

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

#### diff_test ####
# difference test；
# 用于多组数据之间两两做差异分析；
# dat输入长数据，一列sample, 一列要比较的数据；
# group输入group信息，一般是第一列是sample，第二列是group；
# group_order输入一个参与差异检验的向量, 支持多个组比较，默认时所有组；
# sample  index1
# s1      20
# s2      31 
# s3      15
# s4      12
# sample  group
# s1      g1
# s2      g1 
# s3      g2
# s4      g2
diff_test <- function(dat, group, group_order = NULL, method = "wilcox", paired = F, add_plab = F) {
  if (is.null(group_order)) group_order <- unique(group$group)
  dat <- dplyr::rename(dat, sample = sample, val = all_of(setdiff(colnames(dat), "sample"))) %>% 
    mutate(group = group$group[match(sample, group$sample)])
  
  diff <- data.frame(group_pair = character(), pval = numeric(), method = character())
  if (method == "wilcox") { # wilcox test
    for (i in 1:(length(group_order) - 1)) {
      for (j in (i + 1):length(group_order)) {  
        dat_i <- dat %>% subset(group %in% group_order[i]) %>% 
          select(-sample, -group) %>% unlist() %>% as.numeric()
        dat_j <- dat %>% subset(group %in% group_order[j]) %>%
          select(-sample, -group) %>% unlist() %>% as.numeric()
        wilcox <- wilcox.test(dat_i, dat_j, paired = paired, exact = F)
        diff <- tibble::add_row(diff, group_pair = paste0(group_order[i], "_vs_", group_order[j]),
                                pval = as.numeric(wilcox$p.value), method = "wilcoxon-rank sum")

      }
    }
  } else if (method == "anova") { # anova 
    for (i in 1:(length(group_order) - 1)) {      
      for (j in (i + 1):length(group_order)) {      
        group_ij <- subset(group, group %in% c(group_order[i], group_order[j]))
        dat_ij <- subset(dat, sample %in% group_ij$sample) %>% 
          rename(index = all_of(setdiff(colnames(.), c("sample", "group"))))
        anova <- oneway.test(index ~ group, dat_ij)  # anova-test 
        diff <- tibble::add_row(diff, group_pair = paste0(group_order[i], "_vs_", group_order[j]),
                                pval = as.numeric(anova$p.value), method = "oneway anova")
      }        
    } 
  } else if (method == "t") {      
    for (i in 1:(length(group_order) - 1)) { # student's t        
      for (j in (i + 1):length(group_order)) {         
        group_ij <- subset(group, group %in% c(group_order[i], group_order[j]))
        dat_ij <- subset(dat, sample %in% group_ij$sample) %>% 
          rename(index = all_of(setdiff(colnames(.), c("sample", "group"))))
        t <- stats::t.test(index ~ group, dat_ij, paired = paired)
        diff <- tibble::add_row(diff, group_pair = paste0(group_order[i], "_vs_", group_order[j]),
                                pval = as.numeric(t$p.value), method = "student's t")
      }        
    }      
  }
  if (isTRUE(add_plab)) {
    diff <- mutate(diff, plab = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "")))) %>% 
      relocate(plab, .after = "pval")
  }
  return(diff)  
}

#### diff_test_profile ####
# 用于对OTU表中所有的feature进行差异分析；
# profile输入OTU表， 行为feature，列为sample；
# group输入group表，一般是第一列是sample，第二列是group；如果不是这种情况，请用group_colnames指定一个改名字的向量
# group_order输入一个参与差异检验的向量, 支持多个组比较，默认时所有组
#           s1  s2  s3  s4
# feature1  20  31  15  12
# feature2  21  32  16  13
# feature3  22  33  17  14
# feature4  23  34  18  15
# sample  group
# s1      g1
# s2      g1 
# s3      g2
# s4      g2
diff_test_profile <- function(profile, group, group_order = NULL, group_colnames = NULL, method = "wilcox", 
                              paired = F, add_plab = F) {
  if (!is.null(group_colnames)) group <- data.frame(group, check.names = F) %>% dplyr::rename(all_of(group_colnames))
  if (is.null(group_order)) group_order <- unique(group$group)
  dat <- data.frame(t(profile), check.names = F) %>% rownames_to_column(var = "sample")
  tmp_vec <- setdiff(colnames(dat), "sample")
  diff <- purrr::map_dfr(tmp_vec, \(x) diff_test(dat = dplyr::select(dat, sample, all_of(x)), group = group, 
                                         group_order = group_order, method = method, paired = paired, 
                                         add_plab = add_plab) %>% 
                   add_column(feature = x, .before = "group_pair"))
  return(diff)  
}
