# Install RColorBrewer if not already installed
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}

# Load necessary libraries
library(ggplot2)
library(RColorBrewer)

# Load the data from CSV
data <- read.csv('all_abundance.tsv', sep = '\t', header = TRUE)

# Reorder the species so that 'Other' comes last if it exists
if("Other" %in% data$species) {
  data$species <- factor(data$species, levels = c(unique(data$species[data$species != "Other"]), "Other"))
} else {
  data$species <- factor(data$species)
}

# Factor the Group column to ensure consistent ordering across facets
data$species <- factor(data$species)
data$Group <- factor(data$Group)

# Manually defined colors with separate RGB hexadecimal components
species_colors <- c(
  "Corynebacterium sp. GD7" = ("#A5CEE0"),
  "Cutibacterium acnes" = ("#2176B9"),
  "Haemophilus parainfluenzae" = ("#B2E08A"),
  "Lactococcus lactis" = ("#34A02C"),
  "Rothia mucilaginosa" = ("#FB9A99"),
  "Staphylococcus aureus" = ("#E41A1A"),
  "Staphylococcus capitis" = ("#FDBF70"),
  "Staphylococcus epidermidis" = ("#FE8100"),
  "Staphylococcus hominis" = ("#C9B3D6"),
  "Staphylococcus saprophyticus" = ("#6A3D9A"),
  "Staphylococcus sp. AL1" = ("#FDFF99"),
  "Streptococcus mitis" = ("#B05929"),
  "Other" = ("#CCCCCC")
)

# Create a faceted stacked barplot with specified colors
plot <- ggplot(data, aes(x = Sample, y = relative_abundance, fill = species)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Group, scales = "free_x", nrow = 1) +  # All plots in one row
  scale_fill_manual(values = species_colors) +
  labs(x = "Sample", y = "Relative Abundance (%)", fill = "Species") +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    strip.background = element_rect(colour = "black", fill = "grey80", size = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Species Abundance by Group") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# Display the plot
print(plot)

# Save the plot
ggsave("species_abundance_by_group_plot.png", plot = plot, width = 12, height = 8, dpi = 300)
