rm(list = ls())

options(stringsAsFactors = F)
pacman::p_load(tidyverse,ggrepel,ggforce,vegan,dplyr,ggplot2,ggsankey,gapminder)
setwd("01/")
##################################################################################################################################################################
library("openxlsx")


data <- read.xlsx("all_RGI_INFO.xlsx") %>%
  mutate(
    drug_list = str_split(Drug, ";") %>%
      map(~ str_trim(.x)) %>%
      map(~ tolower(.x)),
    has_macrolide     = map_lgl(drug_list, ~ any(str_detect(.x, "macrolide"))),
    has_lincosamide   = map_lgl(drug_list, ~ any(str_detect(.x, "lincosamide"))),
    has_streptogramin = map_lgl(drug_list, ~ any(str_detect(.x, "streptogramin"))),
    n_classes = map_int(drug_list, length),
    Drug_Rename = case_when(
      has_macrolide & has_lincosamide & has_streptogramin & n_classes == 3 ~ "M-L-S",
      n_classes > 1 ~ "Multi-type drug",
      TRUE ~ Drug
    ),
    Mechanism_Rename = ifelse(
      str_detect(Mechanism, ";"),
      "Multi-type mechanism",
      Mechanism
    )
  ) %>%
  select(-drug_list, -has_macrolide, -has_lincosamide, -has_streptogramin, -n_classes)


data <- data %>%
  mutate(across(everything(), ~ replace_na(.x, "Unknown")))

#write.table(data,"ARGs_info_rename.txt",sep = "\t", row.names = F, quote = F)
###################################ARO_abundance.txt#####################################################################################################################################
ARO_abundance <- read.delim("rgi_tpm.profile",sep = "\t",row.names = 1) %>% 
  apply(.,2,function(x) ifelse(x > 0.00, x, 0)) %>% data.frame() %>%    
  mutate(name=data$ARO[match(rownames(.),data$Gene)]) %>% 
  group_by(name)%>%
  summarise_all(sum) 

#write.table(ARO_abundance,"ARO_abundance.txt",sep = "\t", row.names = F, quote = F)
###############################################################################################################################################
dat <- read.delim("ARGs_info_rename.txt") %>% mutate(n=1) %>% 
  select(MAG,n) %>% group_by(MAG) %>% summarise(n=sum(n))

#write.table(dat,"number_MAG_arg.txt",sep = "\t", row.names = F, quote = F)

######
MAG_number <- read.delim("number_MAG_arg.txt") %>% 
  mutate(ARG_n = case_when(
    n == 1 ~ "1 ARG",
    n >= 2 & n <= 50 ~ "2-50 ARGs",
    n >= 51 & n <= 60 ~ "51-60 ARGs",
    n > 60 ~ ">60 ARGs",
    TRUE ~ NA_character_  
  )) %>% 
  mutate(a=1) %>% 
  select(ARG_n,a) %>% 
  group_by(ARG_n) %>% 
  summarise_all(sum)

#
MAG_number$ARG_n <- factor(MAG_number$ARG_n, 
                           levels = c("1 ARG", "2-50 ARGs", "51-60 ARGs", ">60 ARGs"),
                           ordered = TRUE)

ggplot(MAG_number, aes(x = ARG_n, y = a)) +
  geom_col(fill = "#F59B7B", width = 0.8) +
  geom_text(aes(label = a), vjust = -0.2, size = 4, color = "black") +  
  labs(x = "", y = "Number of Genome") +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black", face = "italic"),
        legend.title = element_text(size = 8, color = "black"),
        aspect.ratio = 1/2)


#################################################################################Number of ARGs in MAG###########################################################
ARG_number <- dat[order(dat$n, decreasing = TRUE),] %>% head(10)
ARG_number <- ARG_number %>%
  mutate(MAG = factor(MAG, levels = ARG_number$MAG[order(ARG_number$n,decreasing = TRUE)]))

ggplot(ARG_number, aes(x = MAG, y = n, fill = n)) +
  geom_col() +
  scale_fill_gradient(low = "#FFE4B5", high = "#F59B7B") +  
  labs(x = "", y = "ARG Number") +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = 0.4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(linewidth = 0.4),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black", face = "italic"),
        legend.title = element_text(size = 8, color = "black"),
        aspect.ratio = 1/2) +
  geom_text(aes(label = n), 
            position = position_stack(vjust = .9),
            size = 4,  
            color = "black")  

########################################################################################################################
ARGs_MAG<-read.delim("ARGs_info_rename.txt")%>% 
  select(MAG,ARG) %>%
  distinct(MAG, ARG) %>% 
  mutate(n=1) %>% 
  select(ARG,n) %>% 
  group_by(ARG) %>% 
  summarise(n=sum(n))

#write.table(ARGs_MAG,"ARGs_MAG_number.txt",sep = "\t", row.names = F, quote = F)

plotdata <- read.delim("ARGs_MAG_number.txt") %>% 
  arrange(desc(n)) %>% 
  mutate(mag_n = case_when(
    n == 1 ~ "1 host",
    n >= 2 & n <= 10 ~ "2-10 hosts",
    n >= 11 & n <= 50 ~ "11-50 hosts",
    n > 50 ~ ">50 hosts",
    TRUE ~ NA_character_  #
  )) %>% 
  mutate(a=1) %>% 
  select(mag_n,a) %>% 
  group_by(mag_n) %>% 
  summarise_all(sum)

###
plotdata$mag_n <- factor(plotdata$mag_n, 
                         levels = c("1 host", "2-10 hosts", "11-50 hosts", ">50 hosts"),
                         ordered = TRUE)

ggplot(plotdata, aes(mag_n, a)) +
  geom_col(fill = "#F59B7B") +
  geom_text(aes(label = a), vjust = 0, size = 4, color = "black") +  
  labs(x = "", y = "ARG Number") +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid.major = element_line(linewidth = .4),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black", face = "italic"),
        legend.title = element_text(size = 8, color = "black"),
        aspect.ratio = 1/2)


##########################################
dat2<-data %>% 
  select(MAG,Family) %>% 
  mutate(gen_n=1) %>% 
  group_by(MAG,Family) %>% 
  summarise_all(sum) %>% 
  ungroup() %>% 
  arrange(desc(gen_n)) %>% 
  head(10) %>% 
  mutate(MAG = factor(MAG, levels = unique(MAG)))

ggplot(dat2, aes(MAG, gen_n)) +
  geom_col(fill = "#F59B7B", width = 0.8) +
  geom_text(aes(label = gen_n), vjust = -0.2, size = 3.5, color = "black") +  
  labs(x = "", y = "ARG Number") +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = 0.4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major = element_line(linewidth = 0.4),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black", face = "italic"),
        legend.title = element_text(size = 8, color = "black"),
        aspect.ratio = 1/2)


#####################################################################################
source("plot_pie.R")
dat <- read.delim("ARGs_info_rename.txt", sep = "\t") %>% 
  select(ARO) %>% unique() %>%
  mutate(n = 1, name = data$Drug_Rename[match(ARO, data$ARO)]) %>% 
  group_by(name) %>% 
  summarise(n = sum(n))

#
tem <- dat[order(dat$n, decreasing = TRUE),] %>% head(10) %>% pull(name)

# 
dat2 <- mutate(dat, name = ifelse(name %in% tem, name, "Other")) %>% 
  group_by(name) %>% summarise(n = sum(n)) %>% 
  ungroup() %>% mutate(n = round(n, 2)) %>%
  arrange(desc(n)) 

#
dat2$name <- factor(dat2$name, levels = dat2$name)

#
myLabel <- paste(dat2$name, "(", round(dat2$n / sum(dat2$n) * 100, 2), "%)", sep = "")

#
color <- c("#F59B7B","#ED8828","#FCC41E","#FFD700","#FFE4B5","#F2EFBB", "#8AB1D2","#6BB7CA","#33ABC1","#A4DDD3","#ABD7EC","#b2df8a", "#8D73BA","#C6B3D3", "#33a02c","#80BA8A",
           "#9CD1CB", "#fccde5", "#7f7f7f",  "#ED9F9B","#81B21F" )



#
ggplot(dat2, aes(x = "", y = n, fill = name)) + 
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y", start = 160) + 
  scale_fill_manual(values = color) + 
  theme_minimal() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 15, face = "bold")) + 
  theme(axis.text.x = element_blank()) +
  labs(fill = "Antibiotic Type") +
  geom_text(aes(label = myLabel), position = position_stack(vjust = 0.5), color = "black", size = 4)  

#plot_pie(dat,top_n =11,fill = "auto")
#ggsave("number_drug2.pdf",width = 6, height = 6)
write.table(dat,"number_drug2_pie.txt",sep = "\t", row.names = F, quote = F)

#########################################################################################################
dat <- read.delim("ARGs_info_rename.txt",sep = "\t") %>% 
  select(ARO) %>% unique() %>%
  mutate(n=1,name=data$Mechanism_Rename[match(ARO,data$ARO)]) %>% 
  group_by(name) %>% 
  summarise(n=sum(n))

dat2 = dat[order(dat$n, decreasing = TRUE),]
dat2$name <- factor(dat2$name, levels = dat2$name)
myLabel = as.vector(dat2$name)   
myLabel = paste(myLabel, "(", round(dat2$n / sum(dat2$n) * 100, 2), "%)", sep = "")   


color <- c("#F59B7B","#ED8828","#FCC41E","#FFD700","#FFE4B5","#F2EFBB", "#8AB1D2","#6BB7CA","#33ABC1","#A4DDD3","#ABD7EC","#b2df8a", "#8D73BA","#C6B3D3", "#33a02c","#80BA8A",
           "#9CD1CB", "#fccde5", "#7f7f7f",  "#ED9F9B","#81B21F" )


ggplot(dat2, aes(x = "", y = n, fill = name)) +
  geom_bar(stat = "identity", width = 1, color = NA) +    
  coord_polar(theta = "y", start = 0) + 
  scale_fill_manual(values = color) +  
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  theme(axis.text.x = element_blank()) +
  geom_text(aes(label = myLabel), position = position_stack(vjust = 0.5), color = "black", size = 4) 

#ggsave("number_mechanism2.pdf",width = 6, height = 6)
write.table(dat,"number_mechanism2_pie.txt",sep = "\t", row.names = F, quote = F)

###########################################drug##abundance###################################################################################################################
source("plot_pie.R")

#
plotdat_drug <- read.delim("ARO_abundance.txt", row.names = 1) %>% 
  mutate(total = rowSums(.)) %>% 
  select(n = total) %>% 
  rownames_to_column(var = "name") %>% 
  mutate(name = data$Drug_Rename[match(name, data$ARO)]) %>% 
  group_by(name) %>% 
  summarise(n = sum(n))

#
tem <- plotdat_drug[order(plotdat_drug$n, decreasing = TRUE),] %>% head(10) %>% pull(name)

#
dat2 <- mutate(plotdat_drug, name = ifelse(name %in% tem, name, "Other")) %>% 
  group_by(name) %>% summarise(n = sum(n)) %>% 
  ungroup() %>% mutate(n = round(n, 2)) %>%
  arrange(desc(n)) 

dat2$name <- factor(dat2$name, levels = dat2$name)

#
myLabel <- paste(dat2$name, "(", round(dat2$n / sum(dat2$n) * 100, 2), "%)", sep = "")

#
custom_colors <- c("#F59B7B","#ED8828","#FCC41E","#FFD700","#FFE4B5","#F2EFBB", "#8AB1D2","#6BB7CA","#33ABC1","#A4DDD3","#ABD7EC","#b2df8a", "#8D73BA","#C6B3D3", "#33a02c","#80BA8A",
                   "#9CD1CB", "#fccde5", "#7f7f7f",  "#ED9F9B","#81B21F" )


#
num_colors <- length(unique(dat2$name))
if (num_colors > length(custom_colors)) {
  stop("Not enough colors in custom_colors for the number of unique categories.")
}

#
ggplot(dat2, aes(x = "", y = n, fill = name)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y", start = 60) + 
  scale_fill_manual(values = custom_colors[1:num_colors]) +  
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 15, face = "bold")) +
  theme(axis.text.x = element_blank()) +
  geom_text(aes(label = myLabel), position = position_stack(vjust = 0.5), color = "black", size = 4) 

#ggsave("abundance_drug2_pie.pdf", width = 6, height = 6)
write.table(plotdat_drug,"abundance_drug2_pie.txt",sep = "\t", row.names = F, quote = F)
#plot_pie(plotdat_drug,top_n= 11,fill = "auto")
##############################################mechanism abundance###################################################################################################################################
#
plotdat_mechanism <- read.delim("ARO_abundance.txt", row.names = 1) %>%
  mutate(total = rowSums(.)) %>% 
  select(n = total) %>% 
  rownames_to_column(var = "name") %>% 
  mutate(name = data$Mechanism_Rename[match(name, data$ARO)]) %>% 
  group_by(name) %>% 
  summarise(n = sum(n))

#
tem <- plotdat_mechanism[order(plotdat_mechanism$n, decreasing = TRUE), ] %>% head(10) %>% pull(name)

#
dat2 <- mutate(plotdat_mechanism, name = ifelse(name %in% tem, name, "Other")) %>% 
  group_by(name) %>% summarise(n = sum(n)) %>% 
  ungroup() %>% mutate(n = round(n, 2)) %>%
  arrange(desc(n))  #

#
dat2$name <- factor(dat2$name, levels = dat2$name)

#
myLabel <- paste(dat2$name, "(", round(dat2$n / sum(dat2$n) * 100, 2), "%)", sep = "")

#
custom_colors <- c("#F59B7B","#ED8828","#FCC41E","#FFD700","#FFE4B5","#F2EFBB", "#8AB1D2","#6BB7CA","#33ABC1","#A4DDD3","#ABD7EC","#b2df8a", "#8D73BA","#C6B3D3", "#33a02c","#80BA8A",
                   "#9CD1CB", "#fccde5", "#7f7f7f",  "#ED9F9B","#81B21F" )


#
num_colors <- length(unique(dat2$name))
if (num_colors > length(custom_colors)) {
  stop("Not enough colors in custom_colors for the number of unique categories.")
}

#
ggplot(dat2, aes(x = "", y = n, fill = name)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y", start = 0) + 
  scale_fill_manual(values = custom_colors[1:num_colors]) +  #
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  theme(axis.text.x = element_blank()) +
  geom_text(aes(label = myLabel), position = position_stack(vjust = 0.5), color = "black", size = 4)  

#plot_pie(plotdat_mechanism,top_n = 11,fill = "auto")
#ggsave("abundance_mechanism_pie.pdf", width = 6, height = 6)
write.table(plotdat_mechanism,"abundance_mechanism_pie.txt",sep = "\t", row.names = F, quote = F)

########################################################## prevalence and abundance of ARGs ##############################################################################################################
#
source("profile_process.R")
metadata <- read.delim("ARGs_info_rename.txt")
profile <- read.delim("ARO_abundance.txt", row.names = 1)
ARO_tpm <- read.delim("rgi_tpm.profile",sep = "\t",row.names = 1) %>% 
  mutate(name=data$ARO[match(rownames(.),data$Gene)]) %>% 
  group_by(name)%>%
  summarise_all(sum) %>% 
  column_to_rownames(var="name")

#
prevalence <- profile_replace(ARO_tpm, 0) %>% profile_prevalence()  #
abundance <- data.frame(avg_tpm = rowMeans(profile)) %>% rownames_to_column(var = "name")

dat <- merge(prevalence, abundance, by = "name") %>% 
  merge(., y = select(metadata, ARO, ARG, drug = Drug_Rename, mechanism = Mechanism_Rename) %>% unique.data.frame(), 
        by.x = "name", by.y = "ARO")

write.table(dat, "ARGs_prevalence.txt", sep = "\t", quote = F, row.names = F)

plotdat <- mutate(dat, drug = forcats::fct_lump_n(drug, 10, ties.method = "last"), avg_tpm = log10(avg_tpm+1e-20)) #

a<-plotdat %>% select(drug) %>% unique()

#
drug_counts <- plotdat %>%
  group_by(drug) %>%
  summarise(total_n = sum(n)) %>%
  ungroup()

#
ordered_drugs <- drug_counts %>%
  mutate(drug = as.character(drug)) %>%
  arrange(drug == "Other", desc(total_n)) %>%
  mutate(is_other = (drug == "Other")) %>%
  arrange(desc(is_other), desc(total_n)) %>%
  pull(drug)

plotdat$drug <- factor(plotdat$drug, levels = ordered_drugs)


colors <- c("grey", "#ED8828", "#FCC41E", "#FFD700", "#F59B7B", "#FFE4B5", 
            "#F2EFBB", "#8AB1D2", "#6BB7CA", "#33ABC1", "#8D73BA", "#7f7f7f", 
            "#ABD7EC", "#b2df8a", "#C6B3D3", "#33a02c", "#80BA8A", "green3",
            "#fccde5", "#ED9F9B", "#81B21F")

#
if (length(ordered_drugs) > length(colors)) {
  stop("颜色数量不足以覆盖所有类别")
}

#
color_mapping <- setNames(colors[1:length(ordered_drugs)], ordered_drugs)

#
ggplot(plotdat, aes(x = prevalence, y = avg_tpm)) +
  geom_vline(xintercept = 60, lty = 2, lwd = .5) +
  geom_hline(yintercept = -2, lty = 2, lwd = .5) +
  geom_point(aes(fill = drug), shape = 21, size = 4, stroke = 0.5, color = "black") +  
  scale_fill_manual(values = color_mapping) +
  labs(x = "Prevalence (%)", y = "Average TPM log(10)") +
  geom_text(data = filter(plotdat, prevalence > 60 & avg_tpm > -2.5), 
            aes(label = ARG), size = 3, vjust = -1) +  
  theme_bw() + 
  theme(axis.line = element_line(linewidth = .4, color = "black"),
        axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        panel.grid = element_blank(),
        aspect.ratio = 1)



################################################## taxa contribute number of ARGs stacked bar plot ################################################################################################################
source("profile_process.R")
source("R plot/taxa.R")
metadata <- read.delim("ARGs_info_rename.txt")

dat<- metadata %>% select(feature=MAG,Kingdom,Phylum,Class,Order,Family,Genus,Species)

#####################species
#
dat2 <- select(dat, MAG = feature, Species) %>% unique() %>% 
  merge(select(metadata, MAG, drug = Drug_Rename) %>% unique(), by = "MAG", all.x = TRUE) %>% 
  mutate(
    Species = forcats::fct_lump_n(Species, 11, ties.method = "last"), 
    drug = forcats::fct_lump_n(drug, 10, ties.method = "last")
  ) %>% 
  group_by(Species, drug) %>% 
  summarise(n = n(), .groups = "drop")


plotdat <- dat2 %>% 
  filter(!Species %in% c("Other", "Unknown", "Unknow")) %>% 
  mutate(Species = fct_drop(Species),
         drug = fct_drop(drug))


species_order <- plotdat %>% 
  group_by(Species) %>% 
  summarise(total_n = sum(n)) %>% 
  arrange(desc(total_n)) %>% 
  pull(Species) %>% as.character()

drug_order <- plotdat %>% 
  group_by(drug) %>% 
  summarise(total_n = sum(n)) %>% 
  arrange(desc(total_n)) %>% 
  # 将 "Other" 
  mutate(is_other = (drug == "Other")) %>%
  arrange(is_other, desc(total_n)) %>%
  pull(drug) %>% as.character()


plotdat <- plotdat %>% 
  mutate(Species = factor(Species, levels = species_order),
         drug = factor(drug, levels = drug_order))

#
drug_color_vec <- c("#F59B7B","#ED8828","#FCC41E","#FFD700","#FFE4B5","#F2EFBB", 
                    "#8AB1D2","#6BB7CA","#33ABC1","#A4DDD3","#ABD7EC","#b2df8a", 
                    "#8D73BA","#C6B3D3", "#33a02c","#80BA8A", "#9CD1CB", "#fccde5", 
                    "#ED9F9B","#81B21F")

#
color_mapping <- setNames(drug_color_vec[1:length(levels(plotdat$drug))], levels(plotdat$drug))
if("Other" %in% names(color_mapping)) {
  color_mapping["Other"] <- "#c1c2c2"
}

#
ggplot(plotdat, aes(Species, n, fill = drug)) +
  geom_bar(stat = "identity", position = position_stack(), color = "black", lwd = .2, width = .8) +
  geom_text(
    data = plotdat %>% group_by(Species) %>% summarize(total = sum(n)),
    aes(Species, total, label = total),
    inherit.aes = FALSE, 
    vjust = -0.5, #
    size = 3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  scale_fill_manual(values = color_mapping) +
  labs(x = "Species", y = "Number of ARGs", fill = "Drug class") +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major = element_line(linewidth = .4),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black", face = "italic"),
        legend.title = element_text(size = 8, color = "black"),
        aspect.ratio = 0.5)


#ggsave("ARGs_number_with_s_barplot.pdf", width = 8, height = 4)

select(dat, MAG = feature, Species) %>%
  merge(select(metadata, MAG, drug = Drug_Rename) %>% unique, by = "MAG", all.x = T) %>% 
  group_by(Species, drug) %>% summarise(n = n()) %>% ungroup() %>%
  write.table("ARGs_number_with_s_dat.txt", sep = "\t", row.names = F, quote = F)
#########################################################################################################################################################################
###################family
dat2 <- select(dat, MAG = feature, Family) %>% unique() %>% 
  merge(select(metadata,MAG,drug=Drug_Rename) %>% unique, by = "MAG", all.x = T) %>% 
  mutate(Family = forcats::fct_lump_n(Family, 10, ties.method = "last"), drug = forcats::fct_lump_n(drug, 10, ties.method = "last")) %>% 
  group_by(Family, drug) %>% summarise(n = n()) %>% ungroup() %>% filter(!Family %in% "Other")

Family_order <- group_by(dat2, Family) %>% summarise(n = sum(n)) %>% arrange(desc(n)) %>% select(Family) %>% unlist() %>% as.character()
drug_order <- group_by(dat2, drug) %>% summarise(n = sum(n)) %>% arrange(desc(n)) %>% select(drug) %>% unlist %>% as.character()
plotdat <- mutate(dat2, Family = factor(Family, Family_order), drug = factor(drug, drug_order))

drug_color <- c("#F59B7B","#ED8828","#FCC41E","#FFD700","#FFE4B5","#F2EFBB", "#8AB1D2","#6BB7CA","#33ABC1","#A4DDD3","#ABD7EC","#b2df8a", "#8D73BA","#C6B3D3", "#33a02c","#80BA8A",
                "#9CD1CB", "#fccde5", "#7f7f7f",  "#ED9F9B","#81B21F" )

ggplot(plotdat, aes(Family, n, fill = drug)) +
  geom_bar(stat = "identity", position = position_stack(), color = "black", lwd = .2, width = .8) +
  geom_text(
    data = plotdat %>% group_by(Family) %>% summarize(total = sum(n)),
    aes(Family, total, label = total),
    inherit.aes = FALSE, #
    vjust = 0, size = 3
  ) +
  scale_y_continuous(expand = c(.02, .02)) +
  scale_fill_manual(values = drug_color) +
  labs(x = "", y = "Number of ARGs", fill = "Drug class") +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major = element_line(linewidth = .4),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black", face = "italic"),
        legend.title = element_text(size = 8, color = "black"),
        aspect.ratio = 1/2)
#ggsave("ARGs_number_with_f_barplot.pdf", width = 8, height = 4)

select(dat, MAG = feature, Family) %>%
  merge(select(metadata, MAG, drug = Drug_Rename) %>% unique, by = "MAG", all.x = T) %>% 
  group_by(Family, drug) %>% summarise(n = n()) %>% ungroup() %>%
  write.table("ARGs_number_with_f_dat.txt", sep = "\t", row.names = F, quote = F)

##########################################################################################################################################
#####phylum

dat2 <- select(dat, MAG = feature, Phylum) %>% unique() %>% 
  merge(select(metadata,MAG,drug=Drug_Rename) %>% unique, by = "MAG", all.x = T) %>% 
  mutate(Phylum = forcats::fct_lump_n(Phylum, 10, ties.method = "last"), drug = forcats::fct_lump_n(drug, 10, ties.method = "last")) %>% 
  group_by(Phylum, drug) %>% summarise(n = n()) %>% ungroup() %>% filter(!Phylum %in% "Other")

Phylum_order <- group_by(dat2, Phylum) %>% summarise(n = sum(n)) %>% arrange(desc(n)) %>% select(Phylum) %>% unlist() %>% as.character()
drug_order <- group_by(dat2, drug) %>% summarise(n = sum(n)) %>% arrange(desc(n)) %>% select(drug) %>% unlist %>% as.character()
plotdat <- mutate(dat2, Phylum = factor(Phylum, Phylum_order), drug = factor(drug, drug_order))

drug_color <- c("#F59B7B","#ED8828","#FCC41E","#FFD700","#FFE4B5","#F2EFBB", "#8AB1D2","#6BB7CA","#33ABC1","#A4DDD3","#ABD7EC","#b2df8a", "#8D73BA","#C6B3D3", "#33a02c","#80BA8A",
                "#9CD1CB", "#fccde5", "#7f7f7f",  "#ED9F9B","#81B21F" )


ggplot(plotdat, aes(Phylum, n, fill = drug)) +
  geom_bar(stat = "identity", position = position_stack(), color = "black", lwd = .2, width = .8) +
  geom_text(
    data = plotdat %>% group_by(Phylum) %>% summarize(total = sum(n)),
    aes(Phylum, total, label = total),
    inherit.aes = FALSE, # 禁用继承全局 aes 映射
    vjust = 0, size = 3
  ) +
  scale_y_continuous(expand = c(.02, .02)) +
  scale_fill_manual(values = drug_color) +
  labs(x = "", y = "Number of ARGs", fill = "Drug class") +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major = element_line(linewidth = .4),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black", face = "italic"),
        legend.title = element_text(size = 8, color = "black"),
        aspect.ratio = 1/2)
#ggsave("ARGs_number_with_f_barplot.pdf", width = 8, height = 4)

select(dat, MAG = feature, Phylum) %>%
  merge(select(metadata, MAG, drug = Drug_Rename) %>% unique, by = "MAG", all.x = T) %>% 
  group_by(Phylum, drug) %>% summarise(n = n()) %>% ungroup() %>%
  write.table("ARGs_number_with_p_dat.txt", sep = "\t", row.names = F, quote = F)

