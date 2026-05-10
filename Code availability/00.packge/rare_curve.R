#### Info ####
# Encoding: utf-8
# Author: Jinxin Meng
# Email: mengjx855@163.com
# Created Data：2022-5-29
# Modified Data: 2023-8-9
# Version: 1.1

library(vegan)

#### F1 ####
# Rarefaction curve analysis
# X: cumulative samples; Y: obs. index 
# Usage: rare.sample.obs(otu, step = 1, random = 10)
rare.sample.obs <- function(otu, step = 1, random = 10) {
  dat <- t(otu)
  n_sample <- dim(dat)[1]
  sampling <- seq(1, n_sample, step)
  if (max(sampling) != n_sample) { sampling <- c(sampling, n_sample) } else { sampleing <- sampling }
  
  # sampling for div index
  result <- data.frame(sample = numeric(), obs = numeric(), sd = numeric())
  
  for (i in sampling) {
    flag = 1
    vec_index <- c()
    while (flag <= random) { # random repeat
      flag = flag + 1
      dat_i <- dat[sample(n_sample, size = i), ]
      
      if (i == 1) { 
        vec_index <- c(vec_index, sum(dat_i > 0)) 
        } else {
        index <- sum(colSums(dat_i) > 0)
        vec_index <- c(vec_index, index)
        }
    }
    
    result <- result %>% 
      add_row(sample = i, obs = mean(vec_index), sd = sd(vec_index)) 
  }
  return(result)
}

#### F2 ####
# Rarefaction curve analysis
# X: cumulative samples; Y: obs. index 
# Usage: rare.sample.obs2(otu)
# Using specaccum function
rare.sample.obs2 <- function(otu){
  sampling <- specaccum(t(otu), method = "random", permutations = 10)
  result <- data.frame(sample = sampling$sites, obs = sampling$richness, sd = sampling$sd)
  return(result)
}

#### F2_1 ####
rare.sample.obs2.group <- function(otu, group){
  group_list <- unique(group$group)
  
  result <- rbind()
  for (i in group_list) {
    sample_i <- group %>% 
      filter(group == i) %>% 
      select(sample) %>% 
      unlist() %>% 
      as.character()
    
    otu_i <- select(otu, all_of(sample_i))
    rare_i <- rare.sample.obs2(otu = otu_i) %>% 
      add_column(group = i)
    
    result <- rbind(result, rare_i)
  }
  return(result)
}


#### F3 ####
# Rarefaction curve analysis
# X: cumulative samples; Y: obs. index 
# Usage: rare.group.obs2(otu, group)
# Using specaccum function
rare.group.obs <- function(otu, group,random=10){
  otu <- t(otu)
  group_list <- unique(group$group)
  
  result <- rbind()
  for (i in group_list) {
    sample_i <- group %>% 
      subset(group%in%i) %>% 
      select(sample) %>% 
      unlist() %>% 
      as.character()
    otu_i <- otu %>% 
      subset(rownames(.)%in%sample_i)
    
    rare <- specaccum(otu_i, method = "random", permutations = random)
    rare <- data.frame(sample = rare$sites, obs = rare$richness, sd = rare$sd) %>% 
      add_column(group = i)
    result <- rbind(result, rare)
  }
  return(result)
}

#### F4 ####
# Rarefaction curve analysis
# sequencing depth; Y: obs. index 
# Usage: rare.dep.obs(otu, step = 2000)
rare.dep.obs <- function(otu, step = 10000000) {
  otu <- t(otu)
  n_sample <- dim(otu)[1]
  dep <- rowSums(otu)
 
  # 按照测序量抽平数据
  sampling <- seq(0, max(dep), step)
  if (max(sampling) < max(dep)) sampling <- c(sampling, max(dep))
  
  dat <- rbind()
  for (i in sampling) {
    rarefied_otu <- rrarefy(otu, i)
    dat_i <- estimateR(rarefied_otu)[1,] %>% 
      data.frame(obs = .) %>% 
      rownames_to_column(var = "sample") %>% 
      add_column(dep = i)
    dat <- rbind(dat, dat_i)
  }
  return(dat)
}

