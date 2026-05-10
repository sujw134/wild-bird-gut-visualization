#### info ####
# encoding: utf-8
# author: Jinxin Meng
# e-mail: mengjx855@163.com
# created data：2023-10-24
# modified data: 2023-12-19
# version: 0.1

# 2023-11-01: add all_group parameter in profile_filter.
# 2023-12-19: add function profile_replace.

library(dplyr)
library(tibble)
library(tidyr)
library(purrr)

#### LOG ####
# LOG transformation method in MaAsLin2;
# The default log transformation incorporated into MaAsLin does add a pseudo-count;
# As is best-known practice currently, the pseudo-count is half the minimum feature; 
# x [numeric]: a numeric vector.
LOG <- function(x) {
  y <- replace(x, x == 0, min(x[x>0]) / 2)
  return(log2(y))
}

# conduct LOG transformation for a otu_table, rows represent features, and columns represent samples.
profile_transLOG <- function(profile) {
  profile <- apply(profile, 1, LOG) %>% t() %>% data.frame(check.names = F)
  return(profile)
}
  
#### profile_filter ####
# profile: input a data.frame of relative abundance profile.
# group: mapping (sample|group), also specify by map_names parameter.
# by_group: filter feature in each group.
# all_group prevalence only meet all group are outputted.
# min_prevalence: threshold of prevalence of features in all sample.
# min_abundance: threshold of abundance in a sample is considered a feature presenting in the sample.
profile_filter <- function(profile, group, group_colnames = NULL, by_group = T, 
                           all_group = F, min_prevalence = 0.1, min_abundance = 0.0) {
  
  profile <- data.frame(profile, check.names = F)
  
  if (isTRUE(by_group) & !missing(group)) {
    if (!all(colnames(group) %in% c("sample", "group")) & is.null(group_colnames)) stop("group field (sample|group)")
    if (!is.null(group_colnames)) group <- data.frame(group, check.names = F) %>% dplyr::rename(all_of(group_colnames))
    
    prevalence <- data.frame(t(profile), check.names = F) %>% 
      mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
      group_by(group) %>% 
      group_modify(~ purrr::map_df(.x, \(x) sum(x > min_abundance)/length(x))) %>%
      ungroup() %>% select(-group)
    
    if (isTRUE(all_group)) { 
      flag_vec <- map_vec(prevalence, \(x) all(x > min_prevalence))
    } else { 
      flag_vec <- map_vec(prevalence, \(x) sum(x > min_prevalence) >= 1) 
    }
    flag_vec <- names(flag_vec[flag_vec])
    profile <- profile[flag_vec, ]
    
  } else if (!isTRUE(by_group)) {
    dat <- data.frame(t(profile), check.names = F)
    flag_vec <- apply(dat, 2, \(x) sum(x > min_abundance)/length(x) > min_prevalence)
    flag_vec <- names(flag_vec[flag_vec])
    profile <- profile[flag_vec, ]
  } else {
    stop("in parameter, if by_group = TRUE, group data.frame should be provided ..")
  }
  return(profile)
}

#### profile_smp2grp ####
# 2023-11-18
# profile column transform
# unsmp_colnames: a vector. some column are not samples
# method: merge column method. please get help in summarise_all() function.
profile_smp2grp <- function(profile, group, group_colnames = NULL, unsmp_colnames = NULL, method = "mean"){
  if (!all(colnames(group) %in% c("sample", "group")) & is.null(group_colnames)) stop("group field (sample|group)")
  if (!is.null(group_colnames)) group <- data.frame(group, check.names = F) %>% dplyr::rename(all_of(group_colnames))
  profile <- data.frame(profile, check.names = F)
  if(!is.null(unsmp_colnames)) {
    profile <- select(profile, -all_of(unsmp_colnames)) %>% t() %>% data.frame(check.names = F) %>% 
      mutate(group = group$group[match(rownames(.), group$sample)]) %>% group_by(group) %>% 
      summarise_all(method) %>% ungroup() %>% column_to_rownames(var = "group") %>%
      t() %>% data.frame(check.names = F) %>% cbind(select(profile, all_of(unsmp_colnames)), .) %>% 
      data.frame(check.names = F)
  } else {
    profile <- data.frame(t(profile), check.names = F) %>% 
      mutate(group = group$group[match(rownames(.), group$sample)]) %>% group_by(group) %>% 
      summarise_all(method) %>% ungroup() %>% column_to_rownames(var = "group") %>%
      t() %>% data.frame(check.names = F)
  }
  return(profile)
}

#### profile_transRA ####
# 2023-11-19
# replace some value meeting parameter to specified value.
profile_transRA <- function(profile){
  profile <- data.frame(profile, check.names = F)  %>% apply(profile, 2, x/sum(x)*100) %>% data.frame()
  return(profile)
}

#### profile_replace ####
# 2023-11-19
# replace some value meeting parameter to specified value.
profile_replace <- function(profile, limit = 1, replace = 0, transRA = F){
  profile <- data.frame(profile, check.names = F) 
  if (isTRUE(transRA)) {
    profile <- apply(profile, 2, x/sum(x)*100) %>% data.frame()
  }
  profile <- apply(profile, 2, \(x) ifelse(x > limit, x, replace)) %>% 
    data.frame()
  profile <- profile[rowSums(profile)!=0,]
  return(profile)
}

#### profile_adjacency ####
# 2023-11-19
# profile to adjacent matrix
# 1 0 1 1 0
# 0 1 1 1 0
# 1 0 0 0 1
profile_adjacency <- function(profile, limit = 0, pos = 1, neg = 0){
  profile <- data.frame(profile, check.names = F) %>% 
    apply(2, \(x) ifelse(x > limit, pos, neg)) %>% 
    data.frame()
  profile <- profile[rowSums(profile)!=0,]
  return(profile)
}

#### profile_prevalence ####
profile_prevalence <- function(profile, limit = 0) {
  prevalence <- data.frame(profile, check.names = F) %>% profile_adjacency(limit) %>% 
    rowSums() %>% data.frame(n = .) %>% rownames_to_column(var = "name") %>% 
    mutate(prevalence = n/(ncol(profile))*100)
  return(prevalence)
}



