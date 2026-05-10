#### info ####
# encoding: utf-8
# author: Jinxin Meng
# created date: 2023-11-13
# modified date: 2023-12-09
# version: 0.1

library(dplyr)
library(tidyr)
library(ggplot2)

#### plot_pie ####
# dat: a data_frame contain two field, name|n;
# names: rename dat, colnames must be name and n;
# desc: descending count;
# order_list: manual list;
# add_n: add count to label;
# add_perc: add percentage to label;
# lab_cir: label circle layout;
# color: border color;
# fill; color fill for each part; if auto, fill multiple colors; if hue, fill hue palette. 
plot_pie <- function(dat, dat_colnames = NULL, desc = T, order_list = NULL, add_n = F, top_n = NULL,
                     add_perc = T, lab_cir = F, title = NULL, color = "white", fill = "grey35") {
  if (!all(colnames(dat) %in% c("name", "n")) & is.null(dat_colnames)) stop("dat field (name|n)")
  if (!is.null(dat_colnames)) dat <- data.frame(dat, check.names = F) %>% dplyr::rename(all_of(dat_colnames))
  plotdat <- data.frame(dat, check.names = F)
  
  if (!is.null(order_list)) { plotdat <- mutate(plotdat, name = factor(name, order_list)) %>% arrange(name)
  } else if (is.null(order_list) & isTRUE(desc)) { plotdat <- arrange(plotdat, desc(n))
  } else if (is.null(order_list) & isFALSE(desc)) { plotdat <- arrange(plotdat, n) }
  
  if (!is.null(top_n) & is.numeric(top_n)) {
    tmp_vec <- head(plotdat, top_n - 1) %>% dplyr::select(name) %>% unlist() %>% as.character()
    plotdat <- mutate(plotdat, name = ifelse(name%in%tmp_vec, name, "Other")) %>% 
      group_by(name) %>% summarise(n = sum(n)) %>% ungroup()
    if (isTRUE(desc)) {
      plotdat <- arrange(plotdat, desc(n))
    }
  }
  
  plotdat <- mutate(plotdat, perc = n/sum(n), ypos = cumsum(perc) - 0.5 * perc,
                    angle = ifelse(ypos < .5, 360 * ypos + 180, 360 * ypos))
  
  if (isFALSE(add_n) & isFALSE(add_perc)) { plotdat <- mutate(plotdat, lab = name)
  } else if (isTRUE(add_n) & isTRUE(add_perc)) { plotdat <- mutate(plotdat, lab = paste0(name, ", ", prettyNum(n, big.mark = ","), ", ", round(perc * 100, 2), "%"))
  } else if (isTRUE(add_n) & isFALSE(add_perc)) { plotdat <- mutate(plotdat, lab = paste0(name, ", ", prettyNum(n, big.mark = ",")))
  } else { plotdat <- mutate(plotdat, lab = paste0(name, ", ", round(perc * 100, 2), "%"))}
  
  if (fill == "auto") {
    fill <- c("#F59B7B" ,"#ED8828","#FCC41E","#FFD700","#FFE4B5","#F2EFBB", "#8AB1D2","#6BB7CA","#33ABC1","#A4DDD3","#ABD7EC","#b2df8a", "#8D73BA","#C6B3D3", "#33a02c","#80BA8A")
    fill <- rep(fill, time = (ceiling(nrow(plotdat)/12)))[1:nrow(plotdat)]
  } else if (fill == "hue") {
    fill <- scales::hue_pal()(nrow(plotdat))
  }

  p <- ggplot(plotdat, aes(x = "", y = perc)) +
    geom_col(width = 1 , color = color, fill = fill, lwd = .4, show.legend = F) +
    coord_polar("y", start = 0) +
    theme_void() +
    theme(aspect.ratio = 1,
          plot.title = element_text(color = "#000000", size = 10, hjust = 0.5))
  
  if (!is.null(title)) { p <- p + labs(title = as.character(title)) }
  
  if (isTRUE(lab_cir)) { 
    p <- p + geom_text(aes(x = 1.55, y = ypos, label = lab, angle = -angle), size = 2, hjust = 0.5, fontface = "italic")
  } else { 
    p <- p + geom_text(aes(x = 1.55, y = ypos, label = lab), size = 2, hjust = 0.5, fontface = "italic")
  }
  
  cat("  ggsave(file = \"pie.pdf\", width = 4, height = 4)\n")
  return(p)
}