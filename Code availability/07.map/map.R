
rm(list = ls())

setwd("07/")

# =========================
# 1
# =========================
library(readxl)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(scatterpie)
library(scales)

# =========================
# 2
# =========================
sample_info <- read_xlsx("采样省份.xlsx")

# =========================
# 3
# =========================
sample_clean <- sample_info %>%
  mutate(Region_Fixed = trimws(`Country/region`)) %>%
  mutate(Region_Fixed = case_when(
    Region_Fixed %in% c("USA", "United States") ~ "United States",
    Region_Fixed == "UK" ~ "United Kingdom",
    Region_Fixed == "InnerMogonia" ~ "Inner Mongolia",
    TRUE ~ Region_Fixed
  ))

# =========================
# 4
# =========================
pie_data <- sample_clean %>%
  group_by(Region_Fixed, Host) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Host,
              values_from = n,
              values_fill = 0)

host_cols <- colnames(pie_data)[-1]

# =========================
# 5
# =========================
host_sum <- colSums(pie_data[host_cols])
top_hosts <- names(sort(host_sum, decreasing = TRUE))[1:10]

pie_data$Other <- rowSums(pie_data[, setdiff(host_cols, top_hosts)])

pie_data <- pie_data[, c("Region_Fixed", top_hosts, "Other")]
host_cols <- c(top_hosts, "Other")

pie_data[host_cols] <- lapply(pie_data[host_cols], as.numeric)

# =========================
# 6
# =========================
world_countries <- ne_countries(scale = "medium", returnclass = "sf")

world_states <- ne_states(returnclass = "sf")

# =========================
# 7
# =========================
pts <- list()

for (i in 1:nrow(pie_data)) {
  
  target <- pie_data$Region_Fixed[i]
  
  res <- world_states %>%
    filter(name_en == target | name == target)
  
  if (nrow(res) == 0) {
    res <- world_countries %>%
      filter(name_en == target | name == target | admin == target)
  }
  
  if (nrow(res) > 0) {
    
    geom <- st_centroid(res[1, ])
    
    pt <- st_sf(
      Region_Fixed = target,
      geometry = st_geometry(geom)
    )
    
  } else {
    next
  }
  
  pts[[i]] <- pt
}

centroids <- bind_rows(pts)

# =========================
# 8
# =========================
pie_map_data <- left_join(centroids, pie_data, by = "Region_Fixed")

coords <- st_coordinates(pie_map_data)

pie_map_data_df <- st_drop_geometry(pie_map_data)
pie_map_data_df$lon <- coords[,1]
pie_map_data_df$lat <- coords[,2]

# =========================
# 9
# =========================
pie_map_data_df$total_n <- rowSums(pie_map_data_df[host_cols])

pie_map_data_df$r <- rescale(
  log10(pie_map_data_df$total_n + 1),
  to = c(1.2, 6)
)

# =========================
# 10
# =========================
host_colors <- c(
  "#F59B7B","#ED8828","#FCC41E","#FFD700","#FFE4B5",
  "#F2EFBB","#8AB1D2","#6BB7CA","#33ABC1","#A4DDD3",
  "#BDBDBD"
)
names(host_colors) <- host_cols

# =========================
# 11
# =========================
legend_sizes <- data.frame(
  lon = min(pie_map_data_df$lon, na.rm = TRUE) - 20,
  lat = min(pie_map_data_df$lat, na.rm = TRUE) - 20,
  total_n = c(
    min(pie_map_data_df$total_n, na.rm = TRUE),
    median(pie_map_data_df$total_n, na.rm = TRUE),
    max(pie_map_data_df$total_n, na.rm = TRUE)
  )
)

legend_sizes$r <- rescale(
  log10(legend_sizes$total_n + 1),
  to = range(pie_map_data_df$r)
)

# =========================
# 12
# =========================
ggplot() +
  
  geom_sf(
    data = world_countries,
    fill = "#fcfcfc",
    color = "#e0e0e0",
    size = 0.5
  ) +
  
  geom_scatterpie(
    data = pie_map_data_df,
    aes(x = lon, y = lat, r = r),
    cols = host_cols,
    color = NA
  ) +
  geom_point(
    data = legend_sizes,
    aes(x = lon, y = lat, size = total_n),
    shape = 21,
    fill = "grey70",
    color = "black"
  ) +
  
  scale_size_continuous(
    name = "Sample size",
    breaks = legend_sizes$total_n,
    labels = round(legend_sizes$total_n)
  ) +
  
  scale_fill_manual(values = host_colors) +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "")


