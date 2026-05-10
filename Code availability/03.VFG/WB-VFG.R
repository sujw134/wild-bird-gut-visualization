rm(list = ls())

options(stringsAsFactors = F)
pacman::p_load(tidyverse,ggrepel,ggforce,vegan,dplyr)
setwd("03/")

##########################################
metadata <- read.delim("all_info",sep = "\t",header = T) 

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
  labs(x = "", y = "VFGs Number") +
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
####################################################################################################################
#
dat3 <- metadata %>% 
  select(gene) %>% unique() %>%
  mutate(n = 1, name = metadata$VFCs_Name[match(gene, metadata$gene)]) %>% 
  group_by(name) %>% 
  summarise(n = sum(n))

#
tem <- dat3[order(dat3$n, decreasing = TRUE),] %>% head(100) %>% pull(name)

#
dat3 <- mutate(dat3, name = ifelse(name %in% tem, name, "Other")) %>% 
  group_by(name) %>% summarise(n = sum(n)) %>% 
  ungroup() %>% mutate(n = round(n, 2)) %>%
  arrange(desc(n))  # 

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
  labs(fill = "VFCs_Name") +  #
  geom_text(aes(label = myLabel), position = position_stack(vjust = 0.5), color = "black", size = 2)  # 添加标签

########
dat2 <- metadata %>% 
  select(gene) %>% unique() %>%
  mutate(n = 1, name = metadata$VFs_name[match(gene, metadata$gene)]) %>% 
  group_by(name) %>% 
  summarise(n = sum(n))

#
tem <- dat2[order(dat2$n, decreasing = TRUE),] %>% head(10) %>% pull(name)

#
dat2 <- mutate(dat2, name = ifelse(name %in% tem, name, "Other")) %>% 
  group_by(name) %>% summarise(n = sum(n)) %>% 
  ungroup() %>% mutate(n = round(n, 2)) %>%
  arrange(desc(n))  

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
  labs(fill = "VFs_name") +  #
  geom_text(aes(label = myLabel), position = position_stack(vjust = 0.5), color = "black", size = 4)  # 添加标签

######################################################################################
#
VFs_abundance <- read.delim("vfdb_tpm.profile",sep = "\t",row.names = 1) %>% 
  mutate(name=metadata$VFs_name[match(rownames(.),metadata$gene)]) %>% 
  group_by(name)%>%
  summarise_all(sum)  


write.table(VFs_abundance,"VFs_abundance.txt",sep = "\t", row.names = F, quote = F)

##abundanc
source("plot_pie.R")
profile<- read.delim("VFs_abundance.txt",sep = "\t",row.names = 1) %>% 
  mutate(total=rowSums(.)) %>% 
  select(n=total) %>% 
  rownames_to_column(var = "name") 

plot_pie(profile, top_n = 11, fill = "auto")  
#ggsave("abundance_mges_pie.pdf", width = 5, height = 5)
##################################################################################################
#
VFCs_abundance <- read.delim("vfdb_tpm.profile",sep = "\t",row.names = 1) %>% 
  mutate(name=metadata$VFCs_Name[match(rownames(.),metadata$gene)]) %>% 
  group_by(name)%>%
  summarise_all(sum)  

write.table(VFCs_abundance,"VFCs_abundance.txt",sep = "\t", row.names = F, quote = F)

##abundance###
profile<- read.delim("VFCs_abundance.txt",sep = "\t",row.names = 1) %>% 
  mutate(total=rowSums(.)) %>% 
  select(n=total) %>% 
  rownames_to_column(var = "name") 

plot_pie(profile, top_n = 15, fill = "auto")  
#ggsave("abundance_mges_pie.pdf", width = 5, height = 5)

##########################################流行率######################################################################
source("profile_process.R")

tpm <-VFs_abundance %>% column_to_rownames(var = "name")

prevalence <- profile_replace(tpm,limit = 0) %>% 
  profile_prevalence()

profile<- read.delim("VFs_abundance.txt",sep = "\t",row.names = 1) 

abundance <- data.frame(avg_tpm = rowMeans(profile)) %>% rownames_to_column(var = "name")

dat <- merge(prevalence, abundance, by = "name") %>% 
  merge(., y = select(metadata, VFs_name, VFCs_Name) %>% unique.data.frame(), 
        by.x = "name", by.y = "VFs_name")
write.table(dat, "VFDB_prevalence.txt", sep = "\t", quote = F, row.names = F)

plotdat <- mutate(dat, VFDB_Type = forcats::fct_lump_n(VFCs_Name, 10, ties.method = "last"), avg_tpm = log10(avg_tpm+1e-20))#取前十，其余归为other#避免输入数值为0（对数为0会为负无穷）.所以后面添加10的－20次方


all_levels <- levels(plotdat$VFDB_Type)

#
my_colors <- c("#ED8828","#FCC41E","#F59B7B","#FFD700","#FFE4B5",
               "#F2EFBB","#8AB1D2","#6BB7CA","#33ABC1","#A4DDD3",
               "#ABD7EC","#b2df8a","#8D73BA","#C6B3D3")

# 
color_mapping <- setNames(my_colors[1:length(all_levels)], all_levels)

# 
if ("Other" %in% all_levels) {
  color_mapping["Other"] <- "#7f7f7f"
}

# 
ggplot(plotdat, aes(x = prevalence, y = avg_tpm)) +
  geom_vline(xintercept = 50, lty = 2, lwd = 0.5) +
  geom_hline(yintercept = -2, lty = 2, lwd = 0.5) +
  # color = "black" 
  geom_point(aes(fill = VFDB_Type), size = 4, 
             shape = 21, color = "black", stroke = 0.5) +  
  scale_fill_manual(values = color_mapping) +  
  labs(x = "Prevalence (%)", y = "Average TPM (log10)", fill = "VFDB Type") +
  
  geom_text(data = filter(plotdat, prevalence > 50 & avg_tpm > -2), 
            aes(label = name), size = 3, vjust = -1) +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.4, color = "black"),
        axis.ticks = element_line(linewidth = 0.4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        panel.grid = element_blank(),
        aspect.ratio = 1)
