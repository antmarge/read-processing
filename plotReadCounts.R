
# Quick one to plot read counts
# MLA

# File paths
home_dir = "./"
coverage_file = paste0(home_dir, "coverage")
read_file = paste0(home_dir, "all_reads.csv")

# Color palette
cbbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Read in coverage data
cov_dat <-
  read.table(
    coverage_file,
    header = F,
    col.names = c("Sample", "Coverage")
  ) %>%
  mutate(Sample = as.character(Sample)) %>%
  mutate(CovBin = ifelse(Coverage <= 0.01, "x <= 0.01", "x > 0.01"))



# Here split by coverage
# Depends on project
# Just plot all if not that many

dat <- read.csv(read_file, header = F, 
                col.names = c("Sample", "Step", "NumReads", "Filepath")) %>%
  left_join(cov_dat, by = "Sample")

bin = "x > 0.01"

dat %>%
  filter(CovBin == bin) %>%
ggplot(aes(x = Step, y = NumReads, fill = Step)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Sample, scales = "free") +
  theme(
    axis.text.x = element_blank()
  ) +
scale_fill_manual(values = cbbPalette) +
  ggtitle(paste0("Hadza project human reads, Coverage ", bin))

