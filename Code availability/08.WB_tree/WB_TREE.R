
rm(list = ls())

# =========================
# 0
# =========================
library(rotl)
library(ape)
library(vegan)
library(dplyr)
library(tidyverse)

# =========================
# 1
# =========================
species <- c(
  "Anser albifrons",
  "Anser cygnoides",
  "Anser fabalis",
  "Anser indicus",
  "Branta canadensis",
  "Corvus dauuricus",
  "Corvus frugilegus",
  "Cygnus atratus",
  "Cygnus cygnus",
  "Cygnus olor",
  "Emberiza spodocephala",
  "Gyps himalayensis",
  "Lagopus muta",
  "Larus ridibundus",
  "Pyrrhocorax pyrrhocorax",
  "Tadorna tadorna",
  "Tarsiger cyanurus"
)

# =========================
# 2
# =========================
resolved <- tnrs_match_names(species)
tree <- tol_induced_subtree(ott_ids = resolved$ott_id)

tree_clean <- tree

tree_clean$tip.label <- gsub("_ott.*", "", tree_clean$tip.label)
tree_clean$tip.label <- gsub("_", " ", tree_clean$tip.label)
tree_clean$tip.label <- trimws(tree_clean$tip.label)

tree_clean <- compute.brlen(tree_clean, method = "Grafen")

# =========================
# 3
# =========================
profile <- read.delim("../ARO_abundance.txt") %>%
  column_to_rownames(var = "name")

profile <- profile[, colSums(profile) > 0]

profile <- sweep(profile, 2, colSums(profile), FUN = "/")

# =========================
# 4
# =========================
sample_group <- read.delim("../group_species_8.txt")

sample_group$sample <- trimws(sample_group$sample)
sample_group$group  <- trimws(sample_group$group)

# =========================
# 5
# =========================
common_samples <- intersect(sample_group$sample, colnames(profile))

sample_group <- sample_group %>%
  filter(sample %in% common_samples)

profile <- profile[, common_samples, drop = FALSE]

# =========================
# 6
# =========================
otu <- profile[, sample_group$sample, drop = FALSE]

otu_t <- as.data.frame(t(otu))
otu_t$sample <- rownames(otu_t)

dat <- merge(sample_group, otu_t, by = "sample")

# =========================
# 7
# =========================
plotdat <- dat %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup()

# =========================
# 8
# =========================
plotdat$group <- trimws(plotdat$group)

plotdat$group <- gsub("Emberiza spodocephala", "Schoeniclus spodocephala", plotdat$group)
plotdat$group <- gsub("Tarsiger cyanurus", "Luscinia cyanura", plotdat$group)
plotdat$group <- gsub("Larus ridibundus", "Chroicocephalus ridibundus", plotdat$group)

# =========================
# 9
# =========================
arg_matrix <- plotdat %>%
  column_to_rownames("group")

#
arg_matrix[is.nan(as.matrix(arg_matrix))] <- 0
arg_matrix[is.infinite(as.matrix(arg_matrix))] <- 0

# =========================
# 10
# =========================
common_species <- intersect(tree_clean$tip.label, rownames(arg_matrix))

cat("Matched species:", length(common_species), "\n")

if(length(common_species) < 2){
  stop("❌ 物种匹配失败，请检查命名")
}

# =========================
# 11
# =========================
tree2 <- drop.tip(tree_clean,
                  setdiff(tree_clean$tip.label, common_species))

arg_matrix2 <- arg_matrix[common_species, , drop = FALSE]

# =========================
# 12
# =========================
phylo_dist2 <- cophenetic(tree2)

arg_matrix2 <- decostand(arg_matrix2, method = "total")

arg_dist2 <- vegdist(arg_matrix2, method = "bray")

# =========================
# 13
# =========================
mantel_res <- mantel(
  phylo_dist2,
  arg_dist2,
  method = "spearman",
  permutations = 9999
)

print(mantel_res)

# =========================
# 14
# =========================
df_plot <- data.frame(
  phylo = as.vector(as.dist(phylo_dist2)),
  arg   = as.vector(as.dist(arg_dist2))
)

ggplot(df_plot, aes(phylo, arg)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
  theme_bw() +
  labs(
    x = "Phylogenetic distance",
    y = "ARG Bray-Curtis distance",
    title = paste("Mantel r =", round(mantel_res$statistic, 3), 
                  "P =", mantel_res$signif)
  )


########################################################################
library(phytools)

#
arg_matrix2 <- arg_matrix2[, colSums(arg_matrix2) > 0]

#
arg_matrix2 <- arg_matrix2[, apply(arg_matrix2, 2, var) > 0]

#
cat("Remaining ARG features:", ncol(arg_matrix2), "\n")

# =========================
# PCA
# =========================
pca_res <- prcomp(arg_matrix2, scale. = TRUE)

########################Phylogenetic Tree with Trait Plot
# =========================
# 0
# =========================
library(phytools)
library(ape)
library(vegan)

# =========================
# 1 resistome trait
# =========================
trait <- pca_res$x[,1]
names(trait) <- rownames(pca_res$x)

trait <- trait[tree2$tip.label]

stopifnot(!any(is.na(trait)))

# =========================
# 2
# =========================
phylosig_res <- phylosig(
  tree2,
  trait,
  method = "lambda",
  test = TRUE
)

#
lambda_val <- phylosig_res$lambda
p_val <- phylosig_res$P

cat("Pagel's lambda =", lambda_val, "\n")
cat("P-value =", p_val, "\n")

# =========================
# 3
# =========================
obj <- contMap(tree2, trait, legend = TRUE, fsize = 0.7)

# =========================
# 4
# =========================
legend(
  "topleft",
  legend = c(
    paste0("Pagel's λ = ", round(lambda_val, 4)),
    paste0("P = ", signif(p_val, 3))
  ),
  bty = "n",
  cex = 0.9
)

# =========================
# 5
# =========================
bar_cols <- rep("#F59B7B", length(trait))
bar_cols[trait < 0] <- "#F59B7B"
names(bar_cols) <- names(trait)

plotTree.barplot(
  tree2,
  trait,
  args.barplot = list(
    col = bar_cols,
    border = NA,
    xlab = "PCA1 (Resistome variation)"
  ),
  fsize = 0.7
)

legend(
  "topleft",
  legend = paste0(
    "λ = ", round(lambda_val, 4),
    "\nP = ", signif(p_val, 3)
  ),
  bty = "n"
)
# =========================
# 6
# =========================
write.csv(
  data.frame(
    species = names(trait),
    PCA1 = trait
  ),
  "PCA1_phylo_trait.csv",
  row.names = FALSE
)



