# Encoding: utf-8
# Author: Jinxin Meng
# Email: mengjx855@163.com
# Created Data：2022-10-31
# Modified Data: 2022-10-31
# Version: 1.0

library(vegan)
get.adjusted.r2 <- function(adonis_object) {
  n_observations <- adonis_object$Df[3]+1
  d_freedom <- adonis_object$Df[1]
  r2 <- adonis_object$R2[1]
  adjusted_r2 <- RsquareAdj(r2, n_observations, d_freedom)
  adjusted_r2
}