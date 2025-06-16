main_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(main_dir)

if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(
  tidyverse,
  patchwork,
  paletteer,
  forcats
)

#####################
# PART 1: BAR PLOTS #
#####################

# Define COG category descriptions
cog_levels <- c(
  'S',
  'G',
  'Q',
  'O',
  'E',
  'U',
  'L',
  'C',
  'J',
  'K',
  'I',
  'T',
  'P',
  'A',
  'H',
  'B',
  'Z',
  'F',
  'M',
  'D',
  'V',
  'Y',
  'W',
  'N'
)
cog_labels <- c(
  "S: Function unknown",
  "G: Carbohydrate transport and metabolism",
  "Q: Secondary metabolites biosynthesis, transport and catabolism",
  "O: Posttranslational modification, protein turnover, chaperones",
  "E: Amino acid transport and metabolism",
  "U: Intracellular trafficking, secretion, and vesicular transport",
  "L: Transcription",
  "C: Energy production and conversion",
  "J: Ribosomal structure and biogenesis",
  "K: Replication, recombination and repair",
  "I: Lipid transport and metabolism",
  "T: Signal transduction mechanisms",
  "P: Inorganic ion transport and metabolism",
  "A: RNA processing and modification",
  "H: Coenzyme transport and metabolism",
  "B: Chromatin structure and dynamics",
  "Z: Cell wall/membrane/envelope biogenesis",
  "F: Cytoskeleton",
  "M: Nucleotide transport and metabolism",
  "D: Defense mechanisms",
  "V: Cell cycle control, cell division, chromosome partitioning",
  "Y: Nuclear structure",
  "W: Extracellular structures",
  "N: Cell motility"
)

# Function to load and process data
load_cog_data <- function(file) {
  read.table(
    file,
    sep = "\t",
    header = FALSE,
    col.names = c("COG_Category", "Count")
  ) %>%
    arrange(desc(Count)) %>%
    filter(Count >= 10) %>%
    mutate(COG_Description = factor(COG_Category, levels = cog_levels, labels = cog_labels))
}

my_palette <- paletteer_d("ggsci::default_igv")[1:21]
cog_levels_for_palette <- cog_levels[1:(length(cog_levels) - 3)]
cog_palette <- setNames(my_palette, cog_levels_for_palette)
cog_labels_for_palette <- cog_labels[1:(length(cog_labels) - 3)]
cog_palette_barplot <- setNames(my_palette, cog_labels_for_palette)

# Load and plot data
complete <- load_cog_data("eggNOG/complete/complete_cog_category_counts_clean.tsv")
characterized <- load_cog_data("eggNOG/characterized/characterized_cog_category_counts_clean.tsv")

# Create the bar chart
complete_plot <- ggplot(complete, aes(
  x = reorder(COG_Category, -Count),
  y = Count,
  fill = COG_Description
)) +
  geom_bar(stat = "identity",
           color = "black",
           width = 0.7) +
  scale_fill_manual(values = cog_palette_barplot, guide = guide_legend(title = NULL)) +
  labs(x = "Functional Category", y = "Number of Sequences") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(
      angle = 0,
      vjust = 1,
      hjust = 0.5
    ),
    panel.grid.major.y = element_line(color = "gray80", linetype = "dashed"),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.text = element_text(size = 8),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      size = 0.5,
      linetype = "solid"
    ),
  ) +
  scale_y_continuous(breaks = seq(0, max(complete$Count), by = 500)) +
  geom_hline(yintercept = 0,
             color = "black",
             size = 1) +
  geom_vline(xintercept = 0,
             color = "black",
             size = 1)

ggsave(
  "imgs/complete_profile.png",
  complete_plot,
  width = 9,
  height = 6,
  dpi = 600
)

characterized_plot <- ggplot(characterized, aes(
  x = reorder(COG_Category, -Count),
  y = Count,
  fill = COG_Description
)) +
  geom_bar(stat = "identity",
           color = "black",
           width = 0.7) +
  scale_fill_manual(values = cog_palette_barplot, guide = guide_legend(title = NULL)) +
  labs(x = "Functional Category", y = "Number of Sequences") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(
      angle = 0,
      vjust = 1,
      hjust = 0.5
    ),
    panel.grid.major.y = element_line(color = "gray80", linetype = "dashed"),
    legend.position = "none"
  ) +
  scale_y_continuous(breaks = seq(0, max(characterized$Count), by = 50)) +
  geom_hline(yintercept = 0,
             color = "black",
             size = 1) +
  geom_vline(xintercept = 0,
             color = "black",
             size = 1)

ggsave(
  "imgs/characterized_profile.png",
  characterized_plot,
  width = 9,
  height = 4,
  dpi = 600
)

###############################
# PART 2: STACKED BAR CHARTS #
##############################

uncharacterized <- load_cog_data("eggNOG/uncharacterized/uncharacterized_cog_category_counts_clean.tsv")

uncharacterized_ra <- uncharacterized %>%
  mutate(rel_abund = (Count / sum(uncharacterized$Count)) * 100, Group = "Uncharacterized") %>%
  select(Group, category = COG_Category, rel_abund) %>%
  arrange(Group, desc(rel_abund))

characterized_ra <- characterized %>%
  mutate(rel_abund = (Count / sum(characterized$Count)) * 100, Group = "Characterized") %>%
  select(Group, category = COG_Category, rel_abund) %>%
  arrange(Group, desc(rel_abund))

complete_ra <- complete %>%
  mutate(rel_abund = (Count / sum(complete$Count)) * 100, Group = "Complete") %>%
  select(Group, category = COG_Category, rel_abund) %>%
  arrange(Group, desc(rel_abund))

complete_ra$rel_abund[complete_ra$category == "J"] <- 3.9161072
complete_ra$rel_abund[complete_ra$category == "K"] <- 3.9161067

combined <- rbind(complete_ra, characterized_ra, uncharacterized_ra)

combined$Group <- factor(combined$Group,
                         levels = c("Characterized", "Uncharacterized", "Complete"))


category_order <- combined %>%
  group_by(category) %>%
  summarise(total_abund = sum(rel_abund)) %>%
  arrange(desc(total_abund)) %>%
  pull(category)

combined$category <- factor(combined$category, levels = category_order)

combined <- combined %>%
  mutate(category = fct_reorder(category, rel_abund, .fun = sum, .desc = FALSE))

# Determine category order for Characterized and Uncharacterized
category_order_characterized <- combined %>%
  filter(Group == "Characterized") %>%
  group_by(category) %>%
  summarise(total_abund = sum(rel_abund)) %>%
  arrange(total_abund) %>%
  pull(category)

category_order_complete <- combined %>%
  filter(Group == "Complete") %>%
  group_by(category) %>%
  summarise(total_abund = sum(rel_abund)) %>%
  arrange(total_abund) %>%
  pull(category)

# Apply order for first two groups and use Uncharacterized's order for Complete
combined <- combined %>%
  mutate(
    category = case_when(
      Group == "Characterized" ~ factor(category, levels = category_order_characterized),
      Group == "Uncharacterized" ~ factor(category, levels = category_order_characterized),
      Group == "Complete" ~ factor(category, levels = category_order_characterized)
    )
  )

# Recreate the bar chart
barchart <- ggplot(combined, aes(x = Group, y = rel_abund, fill = category)) +
  geom_bar(
    stat = "identity",
    width = 0.9,
    color = "black",
    size = 0.3,
    position = "stack"
  ) +
  scale_fill_manual(values = cog_palette) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     expand = c(0, 0)) +
  scale_x_discrete(
    breaks = c("Characterized", "Uncharacterized", "Complete"),
    labels = c("Characterized", "Uncharacterized", "Complete")
  ) +
  labs(x = "Profile", y = "Relative Abundance (%)") +
  theme_classic() +
  theme(
    axis.text.x = ggtext::element_markdown(
      size = 10,
      angle = 0,
      hjust = 0.5
    ),
    legend.position = "none"
  ) +
  coord_flip() +
  theme(plot.margin = margin(5, 15, 5, 5))

# Combine the plots with specific heights
combined_plot <- complete_plot + characterized_plot + barchart +
  plot_layout(heights = c(6, 4, 2)) +
  plot_annotation(tag_levels = 'A')

# Save the combined plot
ggsave(
  "imgs/combined_functional_profiles.png",
  combined_plot,
  width = 10,
  height = 12,
  dpi = 600
)
