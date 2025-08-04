
library(tidyverse)
library(sf)
library(mapview)
library(marmap) 
library(ncdf4)
library(ggmap)
library(tmap)
library(chron)
library(vegan)
library(dplyr)
library(tidyr)
library(readxl)
library(reshape2)
library(ggplot2)
library(glmnet)
library(pROC)
library(psych)
library(ca)
library(ggrepel)
library(plyr)
library(ade4)
library(magrittr)
library(ggpubr)
library(voxel)
library(caret)
library(spdep)
library(factoextra)
library(corrplot)
library(sp)
library(heatmaply)
library(lubridate)
library(irr)
library(gridExtra)
library(Metrics)
library(biscale)
library(cowplot)
library(patchwork)
library(ggradar)
library(scales)
library(showtext)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(boot)
library(RColorBrewer)
library(PRROC)
library(car)
library(sjPlot)
library(MVN)
library(biotools)
library(rrcov)
library(ggfortify)
library(glmmTMB)


# Data---


# Load pre-processed OTB and BD data
load("Data/input_data.RData")

dates_track <- unique(as.character(track_tot$Date))
dates_OTB <- unique(as.character(OTB_points$DATE))

# Create grid with cells visited at least once 
aoi_grid_m <- data.frame()

for (i in 1:length(dates_track)) {
  daily_track <- filter(track_tot, as.character(Date) == dates_track[i])
  aoi_grid_temp <- aoi_grid
  aoi_grid_temp$track = lengths(st_intersects(aoi_grid, daily_track))
  aoi_grid_temp <- aoi_grid_temp %>% filter(track > 0)
  aoi_grid_m <- rbind(aoi_grid_m, aoi_grid_temp)
}
aoi_grid_m$track <- NULL

aoi_grid_m <- aoi_grid_m %>%
  distinct(grid_id, aoi_grid, .keep_all = TRUE)

unique_OTB <- unique(OTB_points$OTB)


# Spatial Analysis----

## ## PERMANOVA + rLDA - Dolphin distribution on trawling and non-trawling days----


dates_sight <- unique(sight_lines$Date)
distdf <- data.frame(Date=dates_sight, Sight_Centroid=st_centroid(sight_lines$geometry))

# Distance from the Tiber mouth of centroids
pp <- st_point(x = c(12.236155, 41.741359))
pp <- st_sfc(pp, crs = st_crs(distdf$geometry))
distdf$dMouth <- as.numeric(st_distance(distdf$geometry, pp))

# Depth of centroids
bathy = readRDS("Data/Environment/med_bath.sqlitebathy.rData")
centroids_coord <- st_coordinates(distdf$geometry)
distdf$depth <- -1 * get.depth(mat = bathy, x = centroids_coord, 
                               locator = FALSE)$depth
# Distance from the coast of centroids
coast9 = read_sf("Data/Environment/GSA9.shp")
distdf$dCoast <- as.numeric(st_distance(distdf$geometry, coast9))

# Distribution of distances and depth on trawling and non-trawling days
distdf$Trawling <- "NTR"
distdf$Trawling[which(distdf$Date %in% dates_OTB == T)] <- "TR"

ggplot(distdf) +
  geom_histogram(aes(x = dMouth, fill = Trawling), color = "white", bins = 30) +  
  facet_wrap(~Trawling) +
  scale_fill_manual(values = c("NTR" = "lightskyblue3", "TR" = "lightcoral")) +
  labs(
    title = "Histogram of dMouth",
    x = "Distance from the river mouth",  
    y = "Count"  
  ) +
  theme_test() +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),  
    axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12),  
    strip.background = element_rect(fill = "lightsteelblue1"),  
    strip.text = element_text(size = 10, face = "bold")  
  )

ggplot(distdf) +
  geom_histogram(aes(x = dCoast, fill = Trawling), color = "white", bins = 30) +  
  facet_wrap(~Trawling) +
  scale_fill_manual(values = c("NTR" = "lightskyblue3", "TR" = "lightcoral")) +
  labs(
    title = "Histogram of dCoast",
    x = "Distance from the coast",  
    y = "Count"  
  ) +
  theme_test() +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),  
    axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12),  
    strip.background = element_rect(fill = "lightsteelblue1"),  
    strip.text = element_text(size = 10, face = "bold")  
  )

ggplot(distdf) +
  geom_histogram(aes(x = depth, fill = Trawling), color = "white", bins = 30) +  
  facet_wrap(~Trawling) +
  scale_fill_manual(values = c("NTR" = "lightskyblue3", "TR" = "lightcoral")) +
  labs(
    title = "Histogram of depth",
    x = "Depth",  
    y = "Count"  
  ) +
  theme_test() +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),  
    axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12),  
    strip.background = element_rect(fill = "lightsteelblue1"),  
    strip.text = element_text(size = 10, face = "bold")  
  )



dds <- as.data.frame(distdf[,c(3,4,5)])
dds.trw <- as.data.frame(distdf[,6])
colnames(dds.trw) <- "Trawling"

# Test for dispersion homogeneity 
dd.disp <- vegdist(dds, method="euclidean")
dispersion <- betadisper(dd.disp, distdf$Trawling)
anova(dispersion)

# Run PERMANOVA
dds.div <- adonis2(dds ~ Trawling, data = dds.trw, permutations = 999, method="euclidean")
dds.div




dds_scaled <- as.data.frame(scale(dds))

group <- as.factor(dds.trw$Trawling)
model_robust_lda <- Linda(dds_scaled, group)
summary(model_robust_lda)

predicted_classes <- predict(model_robust_lda, newdata = dds_scaled)@classification

conf_mat <- confusionMatrix(predicted_classes, group)
conf_mat

coef_mat <- model_robust_lda@ldf
ld1_vector <- as.numeric(coef_mat["TR", ] - coef_mat["NTR", ])
lda_scores <- as.matrix(dds_scaled) %*% ld1_vector

hist(lda_scores[group == "TR"], 
     breaks = 20, 
     col = adjustcolor("lightcoral", alpha.f = 0.5), 
     xlim = range(lda_scores), 
     main = "Robust LDA LD1 Scores", 
     xlab = "LD1")

hist(lda_scores[group == "NTR"], 
     breaks = 20, 
     col = adjustcolor("lightskyblue3", alpha.f = 0.5), 
     add = TRUE)

legend("topright", 
       legend = c("NTR", "TR"), 
       fill = c(adjustcolor("lightskyblue3", 0.5), adjustcolor("lightcoral", 0.5)))


lda_df <- data.frame(LD1 = lda_scores,
                     Group = dds.trw$Trawling)

ggplot(lda_df, aes(x = LD1, fill = Group)) + 
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("NTR" = "lightskyblue3", "TR" = "lightcoral")) +
  labs(
    title = "Robust LDA: Discrimination by Trawling Activity Day",
    x = "LD1",
    y = "Density"
  ) +
  theme_test() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )


## Spatial overlap metrics----

### Total effort----

BD_track <- track_sight %>% 
  distinct(Date, geometry, .keep_all = TRUE)
aoi_BD_ALL <- aoi_grid_m  # grid only with cells visited by the research vessel
aoi_BD_ALL = st_set_crs(aoi_BD_ALL, st_crs(4326))
time_interval_hours <- 1/6
aoi_OTB_ALL <- aoi_grid
aoi_OTB_ALL = st_set_crs(aoi_OTB_ALL, st_crs(4326))
BD_track$Date <- as.Date(BD_track$Date)
BD_track$MonthYear <- format(BD_track$Date, "%Y-%m")
OTB_points$MonthYear <- format(OTB_points$DATE, "%Y-%m")
BD_month_years <- unique(BD_track$MonthYear)


# Initialize an empty list to store RO values for each month-year
RO_monthly <- list()

# Initialize an empty dataframe to store n points for each month
df.combined <- data.frame()

# Loop through each month-year
for (my in BD_month_years) {
  
  # Subset BD sightings and OTB effort data for the current month-year
  BD_track_my <- BD_track[BD_track$MonthYear == my, ]
  OTB_points_my <- OTB_points[OTB_points$MonthYear == my, ]
  
  # Calculate the number of BD sightings for each grid cell for the current month-year
  aoi_BD_ALL$n_sighting_points <- lengths(st_intersects(aoi_BD_ALL, BD_track_my))
  aoi_BD_ALL$MonthYear <- my
  
  # Calculate the number of fishing points for each grid cell for the current month-year
  aoi_OTB_ALL$n_fish <- lengths(st_intersects(aoi_OTB_ALL, OTB_points_my))
  aoi_OTB_ALL$fishing_hours <- aoi_OTB_ALL$n_fish * time_interval_hours
  aoi_OTB_ALL$MonthYear <- my
  
  # Merge the two data frames for the current month-year
  aoi_BD_avg_nogeo <- st_drop_geometry(aoi_BD_ALL)
  aoi_OTB_avg_nogeo <- st_drop_geometry(aoi_OTB_ALL)
  aoi_combined <- left_join(aoi_OTB_avg_nogeo, aoi_BD_avg_nogeo, by = "grid_id")
  aoi_combined <- st_as_sf(cbind(aoi_OTB_ALL$aoi_grid, aoi_combined))
  
  # Calculate 75th percentile threshold for BD sightings
  sighting_75th_my <- quantile(aoi_combined$n_sighting_points, 0.75, na.rm = TRUE)
  aoi_combined$high_sightings <- aoi_combined$n_sighting_points >= sighting_75th_my
  
  # Calculate 75th percentile threshold for OTB effort
  # Note: don't define hotspots when total effort is zero
  total_fishing <- sum(aoi_combined$fishing_hours, na.rm = TRUE)
  if (total_fishing == 0) {
    aoi_combined$high_fishing <- FALSE 
  } else {
    fishing_75th_my <- quantile(aoi_combined$fishing_hours, 0.75, na.rm = TRUE)
    aoi_combined$high_fishing <- aoi_combined$fishing_hours >= fishing_75th_my
  }
  
  # Identify overlap cells (high BD sightings and high OTB effort)
  aoi_combined$overlap <- aoi_combined$high_sightings & aoi_combined$high_fishing
  aoi_combined$overlap[is.na(aoi_combined$high_sightings)] <- NA
  
  # Calculate RO metric for the current month-year
  A_igt_my <- sum(aoi_combined$overlap, na.rm = TRUE)  # Overlap cells
  A_it_my <- sum(aoi_combined$high_sightings, na.rm = TRUE)  # High-sighting cells
  
  if (total_fishing == 0 || A_it_my == 0) {
    RO_my <- 0  # No fishing, or no high-sighting cells → no overlap
  } else {
    RO_my <- A_igt_my / A_it_my
  }
  
  # Store the RO value for the current month-year
  RO_monthly[[as.character(my)]] <- RO_my
  
  # Add the result to the dataframe of monthly values
  m_combined <- dplyr::select(aoi_combined, c("grid_id", "n_sighting_points",
                                              "fishing_hours"))
  df.combined <- rbind(df.combined, m_combined)
}


RO_monthly_df <- data.frame(
  MonthYear = names(RO_monthly),
  RO = unlist(RO_monthly)
)

RO_monthly_df$Year <- as.numeric(substr(RO_monthly_df$MonthYear, 1, 4))
RO_monthly_df$Month <- as.numeric(substr(RO_monthly_df$MonthYear, 6, 7))


# Calculate the mean RO value across all months
mean_RO <- mean(RO_monthly_df$RO, na.rm = TRUE)
print(paste("Mean Range Overlap (RO) across all months (hotspots with 75th percentile):", mean_RO))


fermo_df <- data.frame(
  start = as.Date(c("2019-09-09", "2020-09-12", "2021-06-12", "2022-06-13", "2023-06-01")),
  end   = as.Date(c("2019-10-08", "2020-10-13", "2021-07-11", "2022-07-12", "2023-06-30"))
)

RO_monthly_df$Date <- as.Date(paste(RO_monthly_df$Year, RO_monthly_df$Month, "01", sep = "-"))


ggplot(RO_monthly_df, aes(x = Date, y = RO)) +
  # Highlight closed seasons
  geom_rect(data = fermo_df, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "firebrick", alpha = 0.3) +
  geom_point(color = "steelblue", size = 2) +
  geom_smooth(method = "loess", se = TRUE, fill = "lightblue", color = "darkblue") +
  labs(
    title = "Monthly Range Overlap (RO) Over Time, with closed seasons",
    x = "Data",
    y = "Range Overlap (RO)"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal()


# Perform bootstrap with 10000 resamples 
boot_mean <- function(data, indices) {
  return(mean(data[indices], na.rm = TRUE))
}

set.seed(123)
boot_results <- boot(data = RO_monthly_df$RO, statistic = boot_mean, R = 10000)

# Get percentile-based 95% CI
boot_ci <- boot.ci(boot_results, type = "perc")

print(boot_ci)

hist(boot_results$t, breaks = 30, main = "Bootstrap Distribution of RO Mean", xlab = "RO Mean")
abline(v = mean(RO_monthly_df$RO, na.rm = TRUE), col = "blue", lwd = 2)



# RO - Spatial Permutation Test

set.seed(123)
n_permutations <- 10000

# Store the original mean RO
observed_mean_RO <- mean(RO_monthly_df$RO, na.rm = TRUE)

# Function to compute RO for all months given randomized BD sightings
compute_permuted_RO <- function() {
  RO_monthly_perm <- list()
  
  for (my in BD_month_years) {
    # Randomly permute BD sighting coordinates across the AOI grid
    
    BD_track_my <- BD_track[BD_track$MonthYear == my, ]
    OTB_points_my <- OTB_points[OTB_points$MonthYear == my, ]
    
    # Permute sighting locations by randomly reassigning them to other grid cells
    permuted_geoms <- st_sample(aoi_OTB_ALL, size = nrow(BD_track_my), type = "random") 
    st_geometry(BD_track_my) <- permuted_geoms
    
    # Calculate sightings per grid cell (permuted)
    aoi_OTB_ALL$n_sighting_points <- lengths(st_intersects(aoi_OTB_ALL, BD_track_my))
    aoi_OTB_ALL$MonthYear <- my
    
    # Original OTB effort
    aoi_OTB_ALL$n_fish <- lengths(st_intersects(aoi_OTB_ALL, OTB_points_my))
    aoi_OTB_ALL$fishing_hours <- aoi_OTB_ALL$n_fish * time_interval_hours
    aoi_OTB_ALL$MonthYear <- my
    
    # Merge and compute RO
    aoi_combined <- aoi_OTB_ALL 
    
    if (sum(aoi_combined$fishing_hours, na.rm = TRUE) == 0) {
      RO <- 0  
    } else {
      
      sighting_75th <- quantile(aoi_combined$n_sighting_points, 0.75, na.rm = TRUE)
      fishing_75th <- quantile(aoi_combined$fishing_hours, 0.75, na.rm = TRUE)
      
      aoi_combined$high_sightings <- aoi_combined$n_sighting_points >= sighting_75th
      aoi_combined$high_fishing <- aoi_combined$fishing_hours >= fishing_75th
      aoi_combined$overlap <- aoi_combined$high_sightings & aoi_combined$high_fishing
      aoi_combined$overlap[which(is.na(aoi_combined$high_sightings))] <- NA
      
      A_igt <- sum(aoi_combined$overlap, na.rm = TRUE)
      A_it <- sum(aoi_combined$high_sightings, na.rm = TRUE)
      
      RO <- if (A_it > 0) A_igt / A_it else 0
      
    }
    
    RO_monthly_perm[[as.character(my)]] <- RO
  }
  
  # Return mean RO across months
  mean(unlist(RO_monthly_perm), na.rm = TRUE)
}


# Run the permutation test (Note: UNCOMMMENT TO RUN, can take hours)

# perm_results <- replicate(n_permutations, compute_permuted_RO())

# Progress bar setup
# pb <- txtProgressBar(min = 0, max = n_permutations, style = 3)

# Run the permutation test with progress bar
# perm_results <- sapply(1:n_permutations, function(i) {
#   result <- compute_permuted_RO()
#   setTxtProgressBar(pb, i)
#   return(result)
# })
# 
# close(pb)

# saveRDS(perm_results, "Meta/perm_results.RData")
perm_results <- readRDS("Meta/perm_results.RData")


# Calculate p-value (proportion of permuted mean RO >= observed)
p_value <- mean(perm_results >= observed_mean_RO)

cat("Observed mean RO:", round(observed_mean_RO, 3), "\n")
cat("Permutation p-value:", round(p_value, 4), "\n")

hist(perm_results, breaks = 30,
     main = "Permutation Test for Mean RO",
     xlab = "Mean RO (Permuted)",
     col = "gray", border = "white",
     xlim = range(c(perm_results, observed_mean_RO+0.01)))
abline(v = observed_mean_RO, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Observed RO"), col = "red", lwd = 2, lty = 2)




# Bivariate cloropleth maps 

# Classify the data into bivariate categories based on averaged monthly values
average_combined <- df.combined %>%
  group_by(grid_id) %>%
  summarise(mean_sighting_points = mean(n_sighting_points, na.rm=T), mean_fishing_hours = mean(fishing_hours, na.rm=T))

colnames(average_combined)[c(2,3)] <- c("n_sighting_points", "fishing_hours")

d = 3 # n breaks


# Quantile style # 

aoi_filtered <- average_combined %>% 
  filter(n_sighting_points > 0) # (without the zeros)

data <- bi_class(aoi_filtered, 
                 x = "fishing_hours", 
                 y = "n_sighting_points", 
                 style = "quantile", 
                 dim = d)

bi_class_breaks(aoi_filtered, x = "fishing_hours", 
                y = "n_sighting_points", 
                style = "quantile", 
                dim = d)


data_transformed <- st_transform(data, crs = st_crs(4326))
data_coords <- data_transformed %>%
  mutate(lon = st_coordinates(st_centroid(data_transformed))[, 1],
         lat = st_coordinates(st_centroid(data_transformed))[, 2])

map <- ggmap(map_aoi) + # base map from get_stadiamap(), maptype = "stamen_terrain_background"
  geom_sf(data = data_coords, aes(fill = bi_class), color = "white", size = 0.1, show.legend = F) +
  bi_scale_fill(pal = "GrPink2", dim = d) +  
  bi_theme() +  
  xlab("Longitude") +
  ylab("Latitude") +
  coord_sf(crs = st_crs(4326)) +
  scale_x_continuous(
    labels = function(x) sprintf("%.1f°E", x), 
    breaks = seq(12.0, 12.5, 0.1) 
  ) +
  scale_y_continuous(
    labels = function(y) sprintf("%.1f°N", y),  
    breaks = seq(41.6, 42.0, 0.1)  
  ) +
  theme(
    axis.text.x = element_text(size = 10, color = "gray30", margin = ggplot2::margin(t = -10, r = 0, b = 0, l = 0)),  
    axis.text.y = element_text(size = 10, color = "gray30", margin = ggplot2::margin(t = 0, r = -10, b = 0, l = 0)), 
    axis.title.x = element_text(size = 12, color = "black", margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0)), 
    axis.title.y = element_text(size = 12, color = "black", margin = ggplot2::margin(t = 0, r = 1, b = 0, l = 0))  
  )

legend <- bi_legend(pal = "GrPink2",
                    dim = d,
                    xlab = "Fishing Effort ",
                    ylab = "Dolphins Presence ",
                    size = 9)

finalPlot <- ggdraw() +
  draw_plot(map, 0, 0, 0.8, 1) +
  draw_plot(legend, 0.2, .6, 1.35, 0.25)

print(finalPlot)


# Jenks style #

average_combined$n_sighting_points[which(is.na(average_combined$n_sighting_points)==T)] <- 0 

data_j <- bi_class(average_combined, x = "fishing_hours", 
                   y = "n_sighting_points", 
                   style = "jenks", 
                   dim = d) 

bi_class_breaks(average_combined, x = "fishing_hours", 
                y = "n_sighting_points", 
                style = "jenks", 
                dim = d)


data_transformed_j <- st_transform(data_j, crs = st_crs(4326))
data_coords_j <- data_transformed_j %>%
  mutate(lon = st_coordinates(st_centroid(data_transformed_j))[, 1],
         lat = st_coordinates(st_centroid(data_transformed_j))[, 2])

map_j <- ggmap(map_aoi) +
  geom_sf(data = data_coords_j, aes(fill = bi_class), color = "white", size = 0.1, show.legend = F) +
  bi_scale_fill(pal = "GrPink2", dim = d) +  
  bi_theme() +  
  xlab("Longitude") +
  ylab("Latitude") +
  coord_sf(crs = st_crs(4326)) +
  scale_x_continuous(
    labels = function(x) sprintf("%.1f°E", x), 
    breaks = seq(12.0, 12.5, 0.1) 
  ) +
  scale_y_continuous(
    labels = function(y) sprintf("%.1f°N", y),  
    breaks = seq(41.6, 42.0, 0.1)  
  ) +
  theme(
    axis.text.x = element_text(size = 10, color = "gray30", margin = ggplot2::margin(t = -10, r = 0, b = 0, l = 0)),  
    axis.text.y = element_text(size = 10, color = "gray30", margin = ggplot2::margin(t = 0, r = -10, b = 0, l = 0)),  
    axis.title.x = element_text(size = 12, color = "black", margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, color = "black", margin = ggplot2::margin(t = 0, r = 1, b = 0, l = 0)) 
  )

legend_j <- bi_legend(pal = "GrPink2",
                      dim = d,
                      xlab = "Fishing Effort ",
                      ylab = "Dolphins Presence ",
                      size = 9)

finalPlot_jenks <- ggdraw() +
  draw_plot(map_j, 0, 0, 0.8, 1) +
  draw_plot(legend_j, 0.2, .6, 1.35, 0.25)

print(finalPlot_jenks)



### Single vessels----

# Create grid cell x day dataframe

grid_sight_fish <- data.frame()

for (i in 1:length(dates_track)) {
  daily_track <- dplyr::filter(track_tot, as.character(Date) == dates_track[i])
  daily_sight <- dplyr::filter(track_sight, Date == dates_track[i])
  daily_OTB <- dplyr::filter(OTB_points, as.character(DATE) == dates_track[i])
  aoi_grid_temp <- aoi_grid
  aoi_grid_temp$track = lengths(st_intersects(aoi_grid, daily_track))
  aoi_grid_temp <- aoi_grid_temp %>% filter(track > 0) # for each day, only keep cells visited by the research vessel (with at least one GPS ping)
  aoi_grid_temp$date <- dates_track[i]
  
  # - Sighting --> 1 = at least one sight point; 0 = no sight points
  # - Fishing --> 1 = at least one OTB point; 0 = no OTB points
  aoi_grid_temp$sighting <- lengths(st_intersects(aoi_grid_temp, daily_sight))
  aoi_grid_temp$sighting <- ifelse(aoi_grid_temp$sighting > 0, 1, aoi_grid_temp$sighting)
  aoi_grid_temp$fishing <- lengths(st_intersects(aoi_grid_temp, daily_OTB))
  aoi_grid_temp$fishing <- ifelse(aoi_grid_temp$fishing > 0, 1, aoi_grid_temp$fishing)
  
  # For each vessel:
  # - Fishing.xxx --> 1 = at least one OTB point; 0 = no OTB points
  for (j in 1:length(unique_OTB)) {
    OTB_j <- dplyr::filter(daily_OTB, OTB == unique_OTB[j])
    intersect_counts <- lengths(st_intersects(aoi_grid_temp, OTB_j))
    intersect_counts <- ifelse(intersect_counts > 0, 1, intersect_counts)
    aoi_grid_temp$fishing.i <- intersect_counts
    colnames(aoi_grid_temp)[ncol(aoi_grid_temp)] <- paste0("fishing.", unique_OTB[j])
  }
  grid_sight_fish <- rbind(grid_sight_fish, aoi_grid_temp)
  
}

cdf <- grid_sight_fish[,c(1,3)]




## Daily Overlap Analysis between single fishing vessels and BD

df_for_overlap <- grid_sight_fish

# List of fishing vessel columns
fishing_vessels <- grep("^fishing\\.", names(df_for_overlap), value = TRUE)

# For each cell, check if there's overlap (sighting == 1 and fishing.x > 0)
df_for_overlap <- df_for_overlap %>%
  mutate(across(all_of(fishing_vessels), 
                ~ . > 0 & sighting == 1,            # Check if both fishing and sighting occurred in the cell
                .names = "overlap_{.col}"))

# Summarize for each date and vessel if there is at least one cell with overlap
overlap_summary <- df_for_overlap %>%
  group_by(date) %>%
  summarise(across(starts_with("overlap_"), 
                   ~ any(. == TRUE),               # Check if there is at least one cell with overlap
                   .names = "any_overlap_{.col}"),
            across(all_of(fishing_vessels), 
                   ~ any(. > 0),                   # Check if the vessel was active on that date
                   .names = "any_active_{.col}")) %>%
  ungroup()


# For each vessel, calculate the number of overlap days and total fishing days
overlap_results <- data.frame(vessel = fishing_vessels, overlap_days = NA, total_fishing_days = NA)

for (vessel in fishing_vessels) {
  
  # Calculate total fishing days (days where vessel was active)
  total_fishing_days <- sum(overlap_summary[[paste0("any_active_", vessel)]])
  
  # Calculate overlap days (days where at least one cell had both fishing and dolphin sighting)
  overlap_days <- sum(overlap_summary[[paste0("any_overlap_overlap_", vessel)]])
  
  # Store results
  overlap_results[overlap_results$vessel == vessel, "overlap_days"] <- overlap_days
  overlap_results[overlap_results$vessel == vessel, "total_fishing_days"] <- total_fishing_days
}

# Calculate overlap index
overlap_results <- overlap_results %>%
  mutate(overlap_index = ifelse(total_fishing_days == 0, 0, overlap_days / total_fishing_days))


# Calculate DCR = DOLPHIN CO-OCCURENCE RATE

# Note: the index is weighted by the total number of active fishing days: more importance to vessels with a larger number of active days, making sure that vessels with only a few days of fishing don't appear disproportionately high due to having a small denominator (total_fishing_days)

overlap_results <- overlap_results %>%
  mutate(DCR = overlap_index * (log(total_fishing_days + 1) / log(max(total_fishing_days + 1))))


# Calculate FiVEI = FISHING VESSEL EXPOSURE INDEX

total_sighting_days <- length(dates_sight)

overlap_results <- overlap_results %>%
  mutate(FiVEI = overlap_days / total_sighting_days)

overlap_results <- overlap_results %>%
  arrange(desc(FiVEI))

overlap_results <- overlap_results %>%
  mutate(OTB = gsub("fishing\\.", "", vessel))


# Visualization - DCR vs FiVEI 

dcr_melt <- melt(overlap_results[, c("OTB", "DCR")], id.vars = "OTB")
fivei_melt <- melt(overlap_results[, c("OTB", "FiVEI")], id.vars = "OTB")

merged_data <- dcr_melt %>%
  rename(DCR = value) %>%
  dplyr::select(OTB, DCR) %>%
  inner_join(fivei_melt %>% rename(FiVEI = value) %>% dplyr::select(OTB, FiVEI), by = "OTB")


ggplot() +
  geom_point(data = merged_data, aes(x = DCR, y = FiVEI), fill = "cadetblue3", 
             size = 5, shape = 21, color = "black", stroke = 0.8) +
  geom_smooth(data = merged_data, aes(x = DCR, y = FiVEI), 
              method = "loess", se = FALSE, color = "blue", linetype = "dashed", size = 1.2) +
  labs(title = "DCR vs FiVEI + LOESS Smooth",
       x = "DCR",
       y = "FiVEI") +
  theme_test() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid = element_blank()
  )

summary(merged_data$DCR)
summary(merged_data$FiVEI)



# Logistic Regression----

# Load df integrated with environmental covariates
scdf <- readRDS("Meta/scdf.rData")



scdf$date <- as.character(scdf$date)
mdf <- left_join(scdf, as.data.frame(grid_sight_fish), by=c("grid_id", "date", "aoi_grid"))
str(mdf)
mdf <- st_drop_geometry(mdf)
hist(mdf$sighting) # distribution of response variable


# Variables selection

# Remove variables with zero variance
if(is_empty(which(apply(mdf[,-c(1,2)], 2, sum) == 0))==F) {
  mdf <- mdf[,-which(colnames(mdf) %in% names(which(apply(mdf[,-c(1,2)], 2, sum) == 0)))]
} 

# Remove highly correlated variables
var <- setdiff(names(mdf[,-c(1,2)]), c("sighting", "fishing", "track"))
fglm <- as.formula(paste("sighting ~ ", paste(var, collapse= "+")))
glm1 <- glm(fglm, data = mdf, family = "binomial")
mdf_vif <- vif(glm1)
names(mdf_vif[which(mdf_vif > 5)])

mdf_cor <- cor(mdf[,which(colnames(mdf) %in% var)])
image(mdf_cor, main = "Correlation Matrix", col = colorRampPalette(c("blue", "white", "red"))(20))

var_env <- var[!grepl("^fishing", var)]
env_cor <- cor(mdf[,which(colnames(mdf) %in% var_env)])
corrplot.mixed(env_cor,
               lower = "ellipse",
               upper = "number",
               tl.pos = "lt",
               diag = "l",
               tl.col = "black")

var2 <- var[-which(var == "dCoast")]
fglm2 <- as.formula(paste("sighting ~ ", paste(var2, collapse= "+")))
glm2 <- glm(fglm2, data = mdf, family = "binomial")
mdf_vif2 <- vif(glm2)
names(mdf_vif[which(mdf_vif2 > 5)])

# Remove fishing variables with too few positive cases 
colSums(mdf[, grepl("^fishing\\.", names(mdf))])
summary(colSums(mdf[, grepl("^fishing\\.", names(mdf))]))
sparse_vars <- names(which(colSums(mdf[, grepl("^fishing\\.", names(mdf))]) < 10))
mdf <- mdf[, !names(mdf) %in% sparse_vars]


# Variables scaling
mdf_og <- mdf
mdf$depth <- scale(mdf$depth)
mdf$dMouth <- scale(mdf$dMouth)
mdf$dSPM <- scale(mdf$dSPM)
mdf$sst <- scale(mdf$sst)
mdf$chl <- scale(mdf$chl)
mdf$chl_lag <- scale(mdf$chl_lag)
mdf$sss <- scale(mdf$sss)


# Fit Unweighted Model
var_m <- setdiff(names(mdf[,-c(1,2)]), c("sighting", "fishing", "track", "dCoast"))
var_m <- c(var_m, "(1 | grid_id)", "(1 | date)")
fglm_m <- as.formula(paste("sighting ~ ", paste(var_m, collapse= "+")))
glmm1 <- glmmTMB(
  formula = fglm_m,
  family = binomial,
  data = mdf,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 2e5, eval.max = 2e5)
  )
)

summary(glmm1)


# Fit Weighted Model

# inverse class-frequency weights
table(mdf$sighting)
total <- nrow(mdf)
n_class0 <- sum(mdf$sighting == 0)
n_class1 <- sum(mdf$sighting == 1)
weights <- ifelse(mdf$sighting == 1, total / (2 * n_class1), total / (2 * n_class0))

# effort scores for non-detections (mixture of true absences and false negatives)
# based on survey effort (track points/time spent in the cell per day)

range_effort <- range(mdf$track[mdf$sighting == 0], na.rm = TRUE)
mdf$effort_score <- mdf$track
mdf$effort_score[mdf$sighting == 0] <- scales::rescale(
  mdf$effort_score[mdf$sighting == 0],
  to = c(1,2) 
)
mdf$effort_score <- ifelse(mdf$sighting == 1, 1, mdf$effort_score)

# composite weights (effort-adjusted weights reflecting both class balance and non-detection reliability)
weights_eff <- ifelse(mdf$sighting == 1, weights, weights*mdf$effort_score)
mdf$weights_eff <- ifelse(mdf$sighting == 1, weights, weights*mdf$effort_score)

# normalize weights to avoid changing the overall magnitude of the log likelihood
weights_eff <- mdf$weights_eff / mean(mdf$weights_eff)


glmm2 <- glmmTMB(
  formula = fglm_m,
  family = binomial,
  data = mdf,
  weights = weights_eff,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 2e5, eval.max = 2e5)
  )
)

summary(glmm2)
plot_model(glmm2, transform = "exp", show.values = TRUE, value.offset = 0.3)
tab_model(glmm2, transform = "exp")


# Diagnostics to compare weighted and unweighted model
roc_curve_m <- roc(response = mdf$sighting, predictor = predict(glmm1, type = "response"))
pROC::auc(roc_curve_m)
actual <- mdf$sighting
predicted_m <- ifelse(predict(glmm1, type = "response") > 0.5, 1, 0)
confusionMatrix(factor(predicted_m), factor(actual), positive = "1")
pr_m <- pr.curve(scores.class0 = predict(glmm1, type = "response")[actual == 1],
                 scores.class1 = predict(glmm1, type = "response")[actual == 0],
                 curve = TRUE)
pr_m$auc.inte

roc_curve_m2 <- roc(response = mdf$sighting, predictor = predict(glmm2, type = "response"))
pROC::auc(roc_curve_m2)

actual <- mdf$sighting
predicted_m2 <- ifelse(predict(glmm2, type = "response") > 0.5, 1, 0)
confusionMatrix(factor(predicted_m2), factor(actual), positive = "1")
pr_m2 <- pr.curve(scores.class0 = predict(glmm2, type = "response")[actual == 1],
                  scores.class1 = predict(glmm2, type = "response")[actual == 0],
                  curve = TRUE)
pr_m2$auc.inte


# Extract vessels with significant effect in the weighted model
model_summary <- summary(glmm2)
model_coefs <- model_summary$coefficients
model_coefs <- model_coefs$cond

mterms <- rownames(model_coefs)
mpvalues <- model_coefs[,4]

significant_fishing_terms <- mterms[grepl("^fishing", mterms) & mpvalues < 0.05]
signfish <- sub(".*\\.", "", significant_fishing_terms)

model_ors <- exp(model_coefs[,1])
signors <- model_ors[which(names(model_ors) %in% significant_fishing_terms)]

signfish_df <- data.frame(OTB = signfish, odds_ratio = signors)



# Analysis of Similarities of vessels' catch----


# Classify vessels based on significance and effect in the weighted glmm

ind_coef_df <- merged_data %>%
  left_join(signfish_df, by = "OTB") %>%
  mutate(category = case_when(
    odds_ratio > 1  ~ "significant, positive effect",
    odds_ratio < 1  ~ "significant, negative effect",
    is.na(odds_ratio) ~ "not significant"
  ))


ggplot() +
  geom_point(data = ind_coef_df, aes(x = DCR, y = FiVEI, fill = category), 
             size = 5, shape = 21, color = "black", stroke = 0.8) +
  geom_smooth(data = ind_coef_df, aes(x = DCR, y = FiVEI), 
              method = "loess", se = FALSE, color = "blue", linetype = "dashed", size = 1.2) +
  scale_fill_manual(values = c(
    "significant, positive effect" = brewer.pal(3, "BrBG")[3],
    "significant, negative effect" = brewer.pal(3, "BrBG")[1],
    "not significant" = "gray60"
  )) +
  labs(title = "DCR vs FiVEI + LOESS Smooth",
       x = "DCR",
       y = "FiVEI",
       fill = "Significance in the GLM") +
  theme_test() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  )


# Load logbooks data
FAR_19_23 <- readRDS("Data/Fishery/FAR_19_23.rData")

# Keep only relevant months
FAR_19_23$DATA <- as.Date(FAR_19_23$DATA, format = "%d/%m/%Y")
FAR_19_23$MONTH <- as.numeric(format(FAR_19_23$DATA, "%m"))
FAR_19_23 <- dplyr::filter(FAR_19_23, MONTH %in% c(4:11))

# Get sum of fished kgs per species per month year
FAR_19_23$MonthYear <- format(FAR_19_23$DATA, "%Y-%m")
FAR_19_23 <- FAR_19_23[,-c(1,3,6)]
aggrFAR <- aggregate(FAR_19_23, .~ OTB + MonthYear + Species, FUN=sum)

# Filter to match fishing effort data
OTB_points$MonthYear <- format(OTB_points$DATE, "%Y-%m")
otbs <- unique(OTB_points$OTB)
aggrFAR_f <- data.frame()
for (i in 1:length(otbs)) {
  OTB_i <- dplyr::filter(OTB_points, OTB==otbs[i])
  FAR_i <- dplyr::filter(aggrFAR, OTB==otbs[i])
  
  FAR_i <- dplyr::filter(FAR_i, MonthYear %in% unique(OTB_i$MonthYear))
  aggrFAR_f <- rbind(aggrFAR_f, FAR_i)
}

# Convert from long to wide format
wideFAR <- dcast(aggrFAR_f, OTB+MonthYear ~ Species, value.var = "Kg")
wideFAR[is.na(wideFAR)] <- 0
wideFAR_clean = wideFAR[which(apply(wideFAR[, -c(1,2)], 1, sum) > 0),] # removing records with sum 0

str(wideFAR_clean)

# Classify vessels based on significance and effect in the weighted model
df_for_comp <- left_join(wideFAR_clean, ind_coef_df[,c(1,5)])
df_for_comp$category <- as.factor(df_for_comp$category)


# nMDS 
df_for_comp = df_for_comp[which(apply(df_for_comp[, -c(1,2,ncol(df_for_comp))], 1, sum) > 0),] # removing records with sum 0 
com = df_for_comp[,3:(ncol(df_for_comp)-1)]
m_com = as.matrix(com)

set.seed(123)
dummy <- capture.output({
  nmds <- metaMDS(m_com, distance = "bray", k = 4)
})
nmds

data.scores = as.data.frame(scores(nmds, display = "sites"))
data.scores$category = df_for_comp$category

ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(shape = category, colour = category)) + 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Category", y = "NMDS2", shape = "Category")  + 
  scale_colour_manual(values = c("brown", "darkblue", "darkgreen")) 


ggplot(data.scores, aes(x = NMDS2, y = NMDS3)) + 
  geom_point(size = 4, aes(shape = category, colour = category)) + 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS2", colour = "Category", y = "NMDS3", shape = "Category")  + 
  scale_colour_manual(values = c("brown", "darkblue", "darkgreen")) 


# Remove outliers and repeat nMDS
center <- colMeans(data.scores[, c("NMDS1", "NMDS2")])
data.scores$ID <- rownames(data.scores)
data.scores$mahal <- mahalanobis(data.scores[, c("NMDS1", "NMDS2")],
                                 center = center,
                                 cov = cov(data.scores[, c("NMDS1", "NMDS2")]))
threshold <- quantile(data.scores$mahal, 0.995)
outliers <- data.scores[data.scores$mahal > threshold, ]
df_clean <- df_for_comp[!rownames(df_for_comp) %in% outliers$ID, ]

com <- df_clean[, 3:(ncol(df_clean)-1)]
m_com <- as.matrix(com)

set.seed(123)
dummy <- capture.output({
  nmds <- metaMDS(m_com, distance = "bray", k = 3)
})

data.scores = as.data.frame(scores(nmds, display = "sites"))
data.scores$category = df_clean$category

ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  stat_ellipse(aes(group = category, colour = category), linetype = 2) +
  geom_point(size = 3, aes(shape = category, colour = category)) + 
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 15), 
        axis.title.x = element_text(face = "bold", size = 15, colour = "black"), 
        legend.title = element_text(size = 15, colour = "black", face = "bold"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Category", y = "NMDS2", shape = "Category")  + 
  scale_colour_manual(values = c(
    "significant, positive effect" = brewer.pal(3, "BrBG")[3],
    "significant, negative effect" = brewer.pal(3, "BrBG")[1],
    "not significant" = "gray60"
  )) 


# ANOSIM
ano = anosim(m_com, df_clean$category, distance = "bray", permutations = 999)
ano 

# Test for dispersion homogeneity 
com.cat <- as.data.frame(df_clean$category)
colnames(com.cat) <- "Category"
temp <- vegdist(m_com, method="bray")
dispersion <- betadisper(temp, df_clean$category)
anova(dispersion)

# PERMANOVA
dds.div <- adonis2(m_com ~ Category, data = com.cat, permutations = 999, method="bray")
dds.div


# SIMPER
simp = simper(m_com, df_clean$category, permutations = 999)
simp

# Get average dissimilarity 
dist_ma <- as.matrix(vegdist(m_com, method = "bray"))
categories <- df_clean$category
average_between_dissimilarity <- function(cat1, cat2) {
  inds1 <- which(categories == cat1)
  inds2 <- which(categories == cat2)
  submatrix <- dist_ma[inds1, inds2]
  mean(submatrix)
}
category_levels <- levels(categories)
between_combinations <- combn(category_levels, 2, simplify = FALSE)
between_dissimilarities <- sapply(between_combinations, function(x) {
  average_between_dissimilarity(x[1], x[2])
})
names(between_dissimilarities) <- sapply(between_combinations, paste, collapse = " vs ")
between_dissimilarities


# Get species with >5% contrubution
simp_list <- list(
  "AB" = list(data = simp$`not significant_significant, positive effect`, cats = c("not significant", "significant, positive effect")),
  "AC" = list(data = simp$`significant, negative effect_not significant`, cats = c("significant, negative effect", "not significant")),
  "BC" = list(data = simp$`significant, negative effect_significant, positive effect`, cats = c("significant, negative effect", "significant, positive effect"))
)


simp_long <- purrr::map_dfr(simp_list, function(obj) {
  simp_df <- obj$data
  cats <- obj$cats
  
  data.frame(
    species = simp_df$species,
    contrib = simp_df$average,
    group1 = simp_df$ava,
    group2 = simp_df$avb
  ) %>%
    filter(contrib > 0.05) %>%
    dplyr::select(-contrib) %>%
    pivot_longer(cols = c(group1, group2),
                 names_to = "group",
                 values_to = "value") %>%
    mutate(category = ifelse(group == "group1", cats[1], cats[2])) %>%
    dplyr::select(species, category, value)
})


df_radar <- simp_long %>%
  group_by(species, category) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = species, values_from = mean_value)

# Radar chart with mean fished kg per species
df_radar <- df_radar[, !names(df_radar) %in% "AES"]

font_add_google("Roboto", "roboto")
showtext_auto()

df_radar$category <- as.factor(df_radar$category)

ggradar(df_radar,
        font.radar = "roboto",
        grid.label.size = 0,  
        axis.label.size = 11,        
        group.point.size = 4, 
        group.colours = c(
          "gray60", 
          RColorBrewer::brewer.pal(3, "BrBG")[3],
          RColorBrewer::brewer.pal(3, "BrBG")[1]
        )) + 
  theme(
    legend.position = "bottom",         
    legend.text = element_text(size = 20, family = "roboto"),
    legend.key = element_rect(fill = NA, color = NA),
    legend.background = element_blank()
  )

