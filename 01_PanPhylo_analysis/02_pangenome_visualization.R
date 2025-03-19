#################################################
# CHAPTER 2. Investigating pangenome structure #
################################################

# Set the dir

main_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(main_dir)

# Load/install required packages

if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(ggplot2, dplyr, tidyr, ggnewscale, ggrepel, patchwork)

# Read the gene_countsset
pangenome_df <- read.csv("pangenome/Pangenome/PanGenome-LeMy.All.prt-clust-0.9-mode1.lst.summary.txt")

#######################
# Part 1: Donut chart #
#######################

# Add new columns to categorize genes
pangenome_df_1 <- pangenome_df %>%
  mutate(
    Core = ifelse(nb_members == 24, "Core", NA),
    # Core genes (present in all genomes)
    Shared = ifelse(nb_members > 1 &
                      nb_members < 24, "Shell", NA),
    # Shared genes (present in some genomes)
    Unique = ifelse(nb_members == 1, "Cloud", NA)  # Unique genes (present in only one genome)
  )

# Count the number of genes in each category
gene_counts <- pangenome_df_1 %>%
  gather(key = "Category", value = "Type", Core, Shared, Unique) %>%
  filter(!is.na(Type)) %>%
  count(Type)

# Compute percentages
gene_counts$fraction = gene_counts$n / sum(gene_counts$n)

# Compute the cumulative percentages (top of each rectangle)
gene_counts$ymax = cumsum(gene_counts$fraction)

# Compute the bottom of each rectangle
gene_counts$ymin = c(0, head(gene_counts$ymax, n = -1))

# Compute label position
gene_counts$labelPosition <- (gene_counts$ymax + gene_counts$ymin) / 2

# Compute a good label
gene_counts$label <- paste0(gene_counts$Type,
                            "\n value: ",
                            gene_counts$n,
                            "\n",
                            round((gene_counts$fraction * 100), 2),
                            "%")

colors <- c(
  "Core" = "#d9e7f1",
  "Shell" = "#743f96",
  "Cloud" = "#8dadd1"
)

# Make the plot
pangenome_donut <- ggplot(gene_counts,
                          aes(
                            ymax = ymax,
                            ymin = ymin,
                            xmax = 4,
                            xmin = 3,
                            fill = Type
                          )) +
  geom_rect() +
  geom_text(x = 1.5,
            aes(y = labelPosition, label = label, color = Type),
            size = 6) + # x here controls label position (inner / outer)
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  coord_polar(theta = "y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")

ggsave(
  "imgs/donut_chart.png",
  pangenome_donut,
  width = 8,
  height = 8,
  dpi = 600
)

########################
# Part 2: Scatter plot #
########################

# Create the base plot
pangenome_scatter <- ggplot(pangenome_df, aes(x = num_fam, y = nb_members)) +
  
  # Highlight cloud genome
  geom_point(
    data = pangenome_df[pangenome_df$nb_members == 1, ],
    aes(x = num_fam, y = nb_members),
    color = "#8dadd1",
    alpha = 0.75,
    size = 1,
    shape = 19,
    stroke = 2
  ) +
  
  # Highlight shell genome
  geom_point(
    data = pangenome_df[pangenome_df$nb_members > 1 & pangenome_df$nb_members < 16, ],
    aes(x = num_fam, y = nb_members),
    color = "#743f96",
    alpha = 0.75,
    size = 1,
    shape = 19,
    stroke = 2
  ) +
  
  # Highlight gene family 314
  geom_point(
    data = pangenome_df[pangenome_df$nb_members >= 16, ],
    aes(x = num_fam, y = nb_members),
    color = "green",
    alpha = 0.75,
    size = 1,
    shape = 21,
    stroke = 2
  ) +
  
  geom_text_repel(
    data = pangenome_df[pangenome_df$nb_members >= 16, ],
    aes(
      x = num_fam,
      y = nb_members,
      label = paste("Gene family:", num_fam)
    ),
    size = 4,
    fontface = "bold",
    nudge_x = c(0, 0, 0, 25, 0, 0, 0, 0),
    nudge_y = c(1, -1, 5, 2, -1 , 1, 1, -1)
  ) +
  
  # Add a horizontal dashed line at y = 70
  geom_hline(
    yintercept = 15,
    color = "black",
    linetype = "dashed",
    alpha = 0.5
  ) +
  
  # Customize axis labels and breaks
  scale_x_continuous(breaks = seq(0, 614, by = 75)) +  # Display every x-axis label
  scale_y_continuous(breaks = seq(0, 24, by = 5)) +  # Display y-axis labels every 10 units
  
  # Remove title and legend
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Gene Family ID", y = "Number of Genomes Sharing the Gene") +
  
  # Ensure everything fits into the plot area
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

ggsave(
  "imgs/scatter_plot.png",
  pangenome_scatter,
  width = 10,
  height = 6,
  dpi = 600
)

######################
###### COMBINED ######
######################

everything <- (pangenome_donut + pangenome_scatter) + plot_annotation(tag_levels = list(c("A", "B")))
ggsave(
  "imgs/combined_pangenome.png",
  plot = everything,
  width = 18,
  height = 8,
  dpi = 600
)