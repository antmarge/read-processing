
# Plot Laser PCA
# MLA 2019.08.06

require(cowplot)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(viridis)

# Set to home directory
home_dir = "./"

# Filepaths
ref_file = paste0(home_dir,"HGDP_938.RefPC.coord")
seq_file  = paste0(home_dir,"laser_hadza-all.SeqPC.coord")
coverage_file = paste0(home_dir,"coverage")
meta_file = paste0(home_dir, "metadata.rds")

colnames = c("Population", "Sample", paste0("PC", seq(1, 100, by = 1)))


# Just plot 15 people per population for reference panel
ref_dat <-
  read.table(
    ref_file,
    col.names = colnames,
    header = T,
    stringsAsFactors = F
  ) %>%
  select(Population, Sample, PC1, PC2) %>%
  group_by(Population) %>%
  filter(row_number() < 15) %>%
  mutate(Data = "Reference")

# Read in coverage data
cov_dat <-
  read.table(
    coverage_file,
    header = F,
    col.names = c("Sample", "Coverage")
  ) %>%
  mutate(Sample = as.character(Sample))
colnames_seq = c("Population", "Sample", "X1", "X2", "X3", "X4", "X5", "PC1", "PC2")

# Read in metadata
meta <- readRDS(file = meta_file) %>%
  select(Sample, BUSH_CAMP) %>%
  mutate(Sample = as.character(Sample))


study_dat <-
  read.table(
    seq_file,
    col.names = colnames_seq,
    header = T,
    stringsAsFactors = F
  ) %>%
  select(Population, Sample, PC1, PC2) %>%
  mutate(Data = "Study") %>%
  mutate(Population = "HadzaStudy") %>%
  mutate(Sample = gsub("_.*", "", Sample)) %>%
  inner_join(cov_dat, by = "Sample") %>%
  left_join(meta, by = "Sample")


ref_dat_labels <-
  ref_dat %>%
  group_by(Population) %>%
  filter(row_number() == 1) %>%
  ungroup()


data <- dplyr::bind_rows(study_dat, ref_dat)

ggplot(data = ref_dat, aes(x = PC1, y = PC2)) +
  geom_point(data = ref_dat,
             alpha = 0.8,
             color = "gray") +
  geom_text_repel(
    data = ref_dat_labels,
    aes(label = Population),
    segment.color = "transparent",
    fontface = "bold",
    point.padding = 1,
    color = "gray",
    size = 8
  ) +
  geom_text_repel(
    data = study_dat %>%
      filter(!is.na(PC1)) %>%
      mutate(BUSH_CAMP = ifelse(is.na(BUSH_CAMP), "No Info", BUSH_CAMP)) %>%
      mutate(BUSH_CAMP = ifelse(
        !BUSH_CAMP %in% c("JEFF_CAMP", "No Info"),
        "Other Camp",
        BUSH_CAMP
      )),
    fontface = "bold",
    segment.color = "gray",
    size = 4,
    aes(label = Sample, color = BUSH_CAMP)
  ) +
  scale_color_manual(name = "BUSH CAMP",
                     values = c("firebrick1", "dodgerblue", "black")) +
  theme(legend.position = "top")

study_dat_bins <- study_dat %>%
  mutate(CovBin = "NA") %>%
  mutate(CovBin = ifelse(Coverage > 0.5, "> 0.5", CovBin)) %>%
  mutate(CovBin = ifelse(Coverage <= 0.5 &
                           Coverage > 0.1, "0.5 >= x > 0.1", CovBin)) %>%
  mutate(CovBin = ifelse(Coverage <= 0.1 &
                           Coverage > 0.05, "0.1 >= x > 0.05", CovBin)) %>%
  mutate(CovBin = ifelse(Coverage <= 0.05 &
                           Coverage > 0.01, "0.05 >= x > 0.01", CovBin)) %>%
  
  mutate(CovBin = ifelse(Coverage <= 0.01, "x <= 0.01", CovBin)) %>%
  filter(!is.na(PC1))


# Create color bins by coverage
bins <- study_dat_bins %>%
  arrange(Coverage) %>%
  pull(CovBin) %>%
  unique()

col_bin = viridis_pal(end = 0.8, option = "A")(n = length(unique(study_dat_bins$CovBin))) %>% rev()
names(col_bin) <- bins

# plot with labels
ggplot(data = ref_dat,
       aes(x = PC1, y = PC2)) +
  geom_point(data = ref_dat,
             alpha = 0.8,
             color = "gray") +
  geom_text_repel(
    data = ref_dat_labels,
    aes(label = Population),
    segment.color = "transparent",
    fontface = "bold",
    point.padding = 1,
    color = "gray",
    size = 8
  ) +
  geom_text_repel(
    data = study_dat_bins,
    fontface = "bold",
    segment.color = "gray",
    size = 4
  ) +
  scale_color_manual(values = col_bin) +
  theme(legend.position = "top")




ggplot(data = ref_dat, aes(x = PC1, y = PC2)) +
  geom_point(data = ref_dat,
             alpha = 0.8,
             color = "gray") +
  geom_text_repel(
    data = ref_dat_labels,
    aes(label = Population),
    segment.color = "transparent",
    fontface = "bold",
    point.padding = 1,
    color = "gray",
    size = 8
  ) +
  geom_point(
    data = study_dat_bins %>%
      mutate(BUSH_CAMP = ifelse(is.na(BUSH_CAMP), "No Info", BUSH_CAMP)) %>%
      mutate(BUSH_CAMP = ifelse(
        !BUSH_CAMP %in% c("JEFF_CAMP", "No Info"),
        "Other Camp",
        BUSH_CAMP
      )),
    size = 3,
    stroke = 2,
    alpha = 0.8,
    aes(color = factor(CovBin, levels = bins), shape = BUSH_CAMP)
  ) +
  scale_shape_manual(values = c(4, 5, 6)) +
  scale_color_manual(name = "Coverage Bin", values = col_bin) +
  theme(legend.position = "top")




