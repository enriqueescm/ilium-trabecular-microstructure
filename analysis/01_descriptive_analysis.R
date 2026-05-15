library(tidyverse)
library(readxl)
library(ggplot2)
library(patchwork)
library(viridis)

setwd("~/ilium-trabecular-microstructure")

# Load the 4 regional datasets
sciatic <- read_excel("data/raw/BASE DATOS - SCIATIC NOTCH.xlsx") %>%
  mutate(region = "Sciatic Notch")

posterior <- read_excel("data/raw/BASE DATOS - POSTERIOR AS.xlsx") %>%
  mutate(region = "Posterior AS")

arcuate <- read_excel("data/raw/BASE DATOS - ARCUATE LINE.xlsx") %>%
  mutate(region = "Arcuate Line")

anterior <- read_excel("data/raw/BASE DATOS - ANTERIOR AS.xlsx") %>%
  mutate(region = "Anterior AS")

# Standardize column names across all datasets
standardize <- function(df) {
  df %>%
    rename_with(~ gsub(" ", "_", .x)) %>%
    rename_with(tolower) %>%
    rename(
      age = matches("age"),
      bvtv = matches("bv.tv|bv/tv"),
      tb_th = matches("tb.th|tb,th"),
      tb_n = matches("tb.n|tb,n")
    ) %>%
    mutate(
      sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
      age = as.numeric(age),
      bvtv = as.numeric(bvtv),
      tb_th = as.numeric(tb_th),
      tb_n = as.numeric(tb_n),
      region = region
    ) %>%
    select(specimen, sex, age, bvtv, tb_th, tb_n, region)
}

sciatic   <- standardize(sciatic)
posterior <- standardize(posterior)
arcuate   <- standardize(arcuate)
anterior  <- standardize(anterior)

# Combine all regions into a single long dataset
data <- bind_rows(sciatic, posterior, arcuate, anterior) %>%
  mutate(region = factor(region, levels = c("Sciatic Notch", "Anterior AS",
                                            "Posterior AS", "Arcuate Line")))

# Quick check
glimpse(data)

# Remove empty rows
data <- data %>% filter(!is.na(specimen) & !is.na(age))

# Verify final sample size
cat("Total observations:", nrow(data), "\n")
cat("Unique individuals:", n_distinct(data$specimen), "\n")
cat("Age range:", min(data$age), "-", max(data$age), "months\n")
cat("Sex distribution:\n")
print(table(data$sex, data$region))

# Check specimen IDs
data %>% 
  filter(region == "Sciatic Notch") %>%
  arrange(specimen) %>%
  pull(specimen)

# ---- DESCRIPTIVE VISUALIZATIONS ----

# Color palette by region
region_colors <- c(
  "Sciatic Notch" = "#E41A1C",
  "Anterior AS"   = "#377EB8",
  "Posterior AS"  = "#4DAF4A",
  "Arcuate Line"  = "#FF7F00"
)

# BV/TV by age and region
p1 <- ggplot(data, aes(x = age, y = bvtv, color = region)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_color_manual(values = region_colors, name = "Region") +
  scale_x_continuous(breaks = 0:14) +
  labs(
    title = "Bone Volume Fraction (BV/TV) by Age",
    x = "Age (months)",
    y = "BV/TV"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# Tb.Th by age and region
p2 <- ggplot(data, aes(x = age, y = tb_th, color = region)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_color_manual(values = region_colors, name = "Region") +
  scale_x_continuous(breaks = 0:14) +
  labs(
    title = "Trabecular Thickness (Tb.Th) by Age",
    x = "Age (months)",
    y = "Tb.Th (mm)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# Tb.N by age and region
p3 <- ggplot(data, aes(x = age, y = tb_n, color = region)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_color_manual(values = region_colors, name = "Region") +
  scale_x_continuous(breaks = 0:14) +
  labs(
    title = "Trabecular Number (Tb.N) by Age",
    x = "Age (months)",
    y = "Tb.N (1/mm)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# Combined panel
p_combined <- p1 / p2 / p3 +
  plot_annotation(
    title = "Trabecular Bone Microstructure of the Ilium — First Year of Life",
    subtitle = "65 individuals | 0-14 months | 4 anatomical regions | microCT analysis",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  ) +
  plot_layout(guides = "collect")

p_combined
ggsave("figures/01_descriptive_by_age.png", p_combined, 
       width = 12, height = 14, dpi = 300)
