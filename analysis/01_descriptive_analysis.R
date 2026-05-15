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

# ---- REGRESSION ANALYSIS ----

# Create age groups: pre-crawling (0-5 months) vs crawling onset (6-14 months)
data <- data %>%
  mutate(age_group = factor(
    ifelse(age <= 5, "Pre-crawling (0-5m)", "Crawling onset (6-14m)"),
    levels = c("Pre-crawling (0-5m)", "Crawling onset (6-14m)")
  ))

# Linear regression BV/TV ~ age for each region
regression_results <- data %>%
  group_by(region) %>%
  summarise(
    slope_bvtv = coef(lm(bvtv ~ age))[2],
    pval_bvtv  = summary(lm(bvtv ~ age))$coefficients[2, 4],
    slope_tbth = coef(lm(tb_th ~ age))[2],
    pval_tbth  = summary(lm(tb_th ~ age))$coefficients[2, 4],
    slope_tbn  = coef(lm(tb_n ~ age))[2],
    pval_tbn   = summary(lm(tb_n ~ age))$coefficients[2, 4]
  ) %>%
  mutate(across(starts_with("slope"), ~ round(.x, 4)),
         across(starts_with("pval"), ~ round(.x, 4)))

print(regression_results)

# Boxplot by age group and region
p4 <- data %>%
  pivot_longer(cols = c(bvtv, tb_th, tb_n),
               names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable,
                           levels = c("bvtv", "tb_th", "tb_n"),
                           labels = c("BV/TV", "Tb.Th (mm)", "Tb.N (1/mm)"))) %>%
  ggplot(aes(x = age_group, y = value, fill = region)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  facet_grid(variable ~ region, scales = "free_y") +
  scale_fill_manual(values = region_colors, guide = "none") +
  labs(
    title = "Trabecular Parameters: Pre-crawling vs Crawling Onset",
    subtitle = "Comparison across anatomical regions | 65 individuals",
    x = NULL,
    y = "Value",
    caption = "Pre-crawling: 0-5 months | Crawling onset: 6-14 months"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

p4
ggsave("figures/02_precrawling_vs_crawling.png", p4,
       width = 14, height = 10, dpi = 300)
