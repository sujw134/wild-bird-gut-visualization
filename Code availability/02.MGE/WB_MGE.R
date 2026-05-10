
rm(list = ls())

options(stringsAsFactors = F)
pacman::p_load(tidyverse,ggrepel,gggenes,ggforce,vegan,dplyr)
setwd("02/")

################################################################################################################
metadata <- read.delim("all.info",sep = "\t",header = F) %>% 
  rename( gene = V1, mag = V5, mge = V2,mge_type = V3, mge_type2 = V4 )

dat <-select(metadata,mag) %>% mutate(n=1) %>% group_by(mag) %>% 
  summarise(n=sum(n))

plotdata<-dat[order(dat$n, decreasing = TRUE),] %>% head(10) %>% 
  mutate(mag = factor(mag, levels = .$mag[order(.$n,decreasing = TRUE)]))

ggplot(plotdata, aes(mag, n)) +
  geom_col(fill = "#F59B7B", width = 0.8) +
  geom_text(aes(label = n), 
            vjust = -0.5,  
            size = 5,      
            color = "black") +  
  labs(x = "", y = "MGEs Number") +
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



###############################################
dat3 <- metadata %>% 
  select(gene) %>% unique() %>%
  mutate(n = 1, name = metadata$mge[match(gene, metadata$gene)]) %>% 
  group_by(name) %>% 
  summarise(n = sum(n))

# 
tem <- dat3[order(dat3$n, decreasing = TRUE),] %>% head(10) %>% pull(name)

#
dat3 <- mutate(dat3, name = ifelse(name %in% tem, name, "Other")) %>% 
  group_by(name) %>% summarise(n = sum(n)) %>% 
  ungroup() %>% mutate(n = round(n, 2)) %>%
  arrange(desc(n))  

#
dat3$name <- factor(dat3$name, levels = dat3$name)

#
myLabel <- paste(dat3$name, "(", round(dat3$n / sum(dat3$n) * 100, 2), "%)", sep = "")

#
color <- c("#F59B7B" ,"#ED8828","#FCC41E","#FFD700","#FFE4B5","#F2EFBB", "#8AB1D2","#6BB7CA","#33ABC1","#A4DDD3","#ABD7EC","#b2df8a", "#8D73BA","#C6B3D3", "#33a02c","#80BA8A")

#
ggplot(dat3, aes(x = "", y = n, fill = name)) + 
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y", start = 160) + 
  scale_fill_manual(values = color) +  # 
  theme_minimal() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) + 
  theme(axis.text.x = element_blank()) +
  labs(fill = "MGE") +  #
  geom_text(aes(label = myLabel), position = position_stack(vjust = 0.5), color = "black", size = 2)  # 添加标签

############################
dat2 <- metadata %>% 
  select(gene) %>% unique() %>%
  mutate(n = 1, name = metadata$mge_type[match(gene, metadata$gene)]) %>% 
  group_by(name) %>% 
  summarise(n = sum(n))

#
tem <- dat2[order(dat2$n, decreasing = TRUE),] %>% head(10) %>% pull(name)

#
dat2 <- mutate(dat2, name = ifelse(name %in% tem, name, "Other")) %>% 
  group_by(name) %>% summarise(n = sum(n)) %>% 
  ungroup() %>% mutate(n = round(n, 2)) %>%
  arrange(desc(n))  #

#
dat2$name <- factor(dat2$name, levels = dat2$name)

#
myLabel <- paste(dat2$name, "(", round(dat2$n / sum(dat2$n) * 100, 2), "%)", sep = "")

#
color <- c("#F59B7B" ,"#ED8828","#FCC41E","#FFD700","#FFE4B5","#F2EFBB", "#8AB1D2","#6BB7CA","#33ABC1","#A4DDD3","#ABD7EC","#b2df8a", "#8D73BA","#C6B3D3", "#33a02c","#80BA8A")

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
  labs(fill = "MGE_type") +  
  geom_text(aes(label = myLabel), position = position_stack(vjust = 0.5), color = "black", size = 4)  # 添加标签

####################################################################################################################################################

MGEs_abundance <- read.delim("mge_tpm.profile",sep = "\t",row.names = 1) %>% 
  mutate(name=metadata$mge[match(rownames(.),metadata$gene)]) %>% 
  group_by(name)%>%
  summarise_all(sum)  


write.table(MGEs_abundance,"MGEs_abundance.txt",sep = "\t", row.names = F, quote = F)

##abundance
source("plot_pie.R")
profile<- read.delim("MGEs_abundance.txt",sep = "\t",row.names = 1) %>% 
  mutate(total=rowSums(.)) %>% 
  select(n=total) %>% 
  rownames_to_column(var = "name") 

plot_pie(profile, top_n = 11, fill = "auto")  
#ggsave("abundance_mges_pie.pdf", width = 5, height = 5)
###################################################################################################################################################
#mge_type
dat <- profile %>% mutate(mge=metadata$mge_type[match(.$name,metadata$mge)]) %>% 
  select(-name) %>% group_by(mge) %>% summarise(n=sum(n)) %>% select(name=mge,n)

plot_pie(dat, top_n = 11, fill = "auto")  
#ggsave("abundance_type_mges_pie.pdf", width = 5, height = 5)

#### prevalence and abundance of ARGs ############################################################################################################
#
profile <- read.delim("MGEs_abundance.txt", row.names = 1)
MGEs_tpm <- read.delim("mge_tpm.profile",sep = "\t",row.names = 1) %>% 
  mutate(name=metadata$mge[match(rownames(.),metadata$gene)]) %>% 
  group_by(name)%>%
  summarise_all(sum) %>% 
  column_to_rownames(var="name")

#
source("profile_process.R")

prevalence <- profile_replace(MGEs_tpm,limit = 0) %>% 
  profile_prevalence()

abundance <- data.frame(avg_tpm = rowMeans(profile)) %>% rownames_to_column(var = "name")

dat <- merge(prevalence, abundance, by = "name") %>% 
  merge(., y = select(metadata, mge, mge_type) %>% unique.data.frame(), 
        by.x = "name", by.y = "mge")

write.table(dat, "MGEs_prevalence.txt", sep = "\t", quote = F, row.names = F)

plotdat <- mutate(dat, MGE_Type = forcats::fct_lump_n(mge_type, 10, ties.method = "last") , avg_tpm = log10(avg_tpm+1e-20)) #取前十，其余归为other#避免输入数值为0（对数为0会为负无穷）.所以后面添加10的－20次方


levels_mge <- levels(plotdat$MGE_Type)


my_colors <- c("#ED8828","#FCC41E","#FFD700","#FFE4B5","#F2EFBB", 
               "#8AB1D2","#6BB7CA","#33ABC1","#A4DDD3","#F59B7B",
               "#ABD7EC", "#b2df8a", "#8D73BA")


color_mapping <- setNames(my_colors[1:length(levels_mge)], levels_mge)
color_mapping["Other"] <- "#7f7f7f" 

# 
ggplot(plotdat, aes(x = prevalence, y = avg_tpm)) +
  geom_vline(xintercept = 40, lty = 2, lwd = 0.5) +
  geom_hline(yintercept = -2, lty = 2, lwd = 0.5) +
  geom_point(aes(fill = MGE_Type), size = 4, shape = 21, color = "black", stroke = 0.5) +
  scale_fill_manual(values = color_mapping) +  #
  labs(x = "Prevalence (%)", y = "Average TPM log(10)", fill = "MGE Type") +
  geom_text(data = filter(plotdat, prevalence > 40 & avg_tpm > -2), 
            aes(label = name), size = 2, vjust = -1) +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.4, color = "black"),
        axis.ticks = element_line(linewidth = 0.4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        panel.grid = element_blank(),
        aspect.ratio = 1)

