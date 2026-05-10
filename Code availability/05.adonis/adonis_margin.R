
rm(list = ls())

library("tidyverse")
source("R2adjust.R")

setwd("05/")
#
group <- read.delim("group.txt", sep = "\t") %>% 
  rename(sample = "Sample")
 
#
otu <- read.delim("../ARO_abundance.txt", header = TRUE, row.names = "name") %>%
  mutate(across(everything(), ~ replace(., is.na(.), 0))) %>%
  select(where(~ !all(. == 0))) %>%
  mutate(across(everything(), ~ . / sum(.) * 100))


# Bray-Curtis
distance <- vegdist(t(otu), method = "bray", na.rm = FALSE) %>% 
  as.matrix(.) %>% 
  data.frame(.)

#
group <- group %>%
  filter(sample %in% colnames(otu)) %>%  
  select(sample, Year, location, Species, Residency.status) %>%
  na.omit() %>%
  data.frame()


dat <- rbind()
for (i in setdiff(colnames(group), "sample")) {
  group_i <- group %>%
    select("sample", all_of(i)) %>% 
    na.omit() %>% 
    data.frame()
  distance_i <- distance[group_i$sample, group_i$sample]
  
  adonis <- adonis2(formula = eval(parse(text = paste0("distance_i ~ ", i))),
                    group_i, by = "margin", permutations = 999, parallel = 4)
  
  dat <- rbind(dat, data.frame(term = i,
                               sample = nrow(group_i),
                               R2 = round(adonis[1,3], digits = 4),
                               R2adjust = get.adjusted.r2(adonis) %>% round(., digits = 4),
                               pval = adonis[1,5]) %>% 
                 mutate(plab = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "")))
                 ))
}

write.table(dat, "adonis_for_margin.txt", sep = "\t", quote = F, row.names = F)

ggplot(dat, aes(reorder(term, desc(R2adjust)), R2adjust * 100, fill = R2adjust * 100)) +
  geom_bar(stat = "identity", width = 0.8, size = 0.5) +
  geom_text(aes(label = plab), nudge_y = 0.05, color = "black", size = 15) +
  scale_fill_gradient(
    low = "#FFE4B5",  
    high = "#F59B7B",
    name = "Adjusted R² (%)" 
  ) +
  labs(x = "", y = "Adjusted R² (%)") +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "right" 
  )


