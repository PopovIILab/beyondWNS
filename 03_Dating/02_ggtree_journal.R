# Set the dir

main_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(main_dir)

# Load/install required packages

#install.packages('BiocManager')

if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(
  ggtree,
  ggimage,
  ggtext,
  phangorn,
  ggplot2,
  treeio,
  ggnewscale,
  viridis,
  phytools,
  patchwork,
  deeptime,
  readr,
  dplyr,
  scales,
  ggrepel
)

####################
# METADATA PARSING #
####################

# Load `metadata`

metadata <- read.table(
  '../01_PanPhylo_analysis/metadata/metadata.tsv',
  sep = '\t',
  header = T
)


##############
# DATED TREE #
##############

# Read the tree file

dated_tree <- read.beast("data/dating_ultra_tree_ready.tree")

# Draft tree

dated_tree_fig <- ggtree(dated_tree, mrsd = "0-01-01") %<+% metadata +
  xlim(-330, 80) +
  geom_range(
    range = 'height_0.95_HPD',
    color = 'red',
    alpha = .6,
    size = 2
  ) +
  geom_rect(
    data = data.frame(
      xmin = c(-350.5, -298.9, -251.902, -201.4, -145, -66, -23),
      xmax = c(-298.9, -251.902, -201.4, -145, -66, -23, -2.58),
      ymin = -Inf,
      ymax = Inf,
      fill = c(
        "#67a599",
        "#f04028",
        "#812b92",
        "#34b2c9",
        "#7fc64e",
        "#fd9a52",
        "#ffe619"
      )
    ),
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = fill
    ),
    inherit.aes = FALSE,
    alpha = 0.15
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  coord_geo(
    xlim = c(-330, 80),
    pos = list("bottom", "bottom"),
    dat = list("epochs", "periods"),
    abbrv = list(TRUE, FALSE),
    expand = TRUE,
    size = "auto",
    neg = TRUE
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x.bottom = element_text(margin = margin(t = 5)),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  
  #geom_text2(aes(subset = !is.na(height_median) &
  #(abs(height_median) > 110 |
  #(abs(height_median) > 50 & abs(height_median) < 55) |
  #(abs(height_median) > 10 & abs(height_median) < 12) |
  #abs(height_median) > 0.01 & abs(height_median) < 0.07),
  #label = round(height_median, 2)),
  #fill = "white", hjust = -0.5, vjust = 0.5, size = 4, color = "black", fontface = "bold") +
  
  scale_x_continuous(
    breaks = c(-2.58, -23, -66, -145, -201.4, -251.902, -298.9, -316.85),
    labels = c("2.58", "23", "66", "145", "201.4", "251.9", "298.9", "316.85")
  ) +
  
  # Add semi-transparent dashed lines at specific positions
  geom_vline(
    xintercept = c(-2.58, -23, -66, -145, -201.4, -251.902, -298.9),
    linetype = "dashed",
    color = "black",
    alpha = 0.15
  ) +
  geom_hilight(
    mapping = aes(subset = node %in% c(34), fill = S),
    fill = "steelblue",
    alpha = .6,
    extend = 100
  ) +
  geom_tiplab(
    aes(
      label = AN_OrganismName,
      fontface = ifelse(grepl("^NC_033907", AN_OrganismName), "bold", "plain"),
    ),
    align = TRUE,
    geom = "label",
    fill = "white",
    label.size = 0
  ) +
  # Add node ages as text labels, excluding NA and 0 values
  geom_label2(
    aes(
      subset = !is.na(height_median) & abs(height_median) > 1e-5,
      label = round(height_median, 2)
    ),
    hjust = 1.2,
    vjust = -0.5,
    size = 4,
    fill = "white",
    alpha = .5,
    fontface = "bold"
  ) +
  geom_label2(aes(subset = !is.na(posterior) & abs(posterior) < 1,
                  x = branch, label = round(posterior, 2)),
              vjust = 1.5,
              size = 3,
              fill = "white",
              alpha = .5,) +
  geom_point2(
    aes(subset = node == 34),
    shape = 21,
    size = 5,
    fill = "red",
    #alpha = 0.5,
    color = "black",
    stroke = 1
  ) +
  geom_point2(
    aes(subset = node == 37),
    shape = 21,
    size = 5,
    fill = "yellow",
    #alpha = 0.5,
    color = "black",
    stroke = 1
  )

ggsave(
  'imgs/dated_ultra_tree.png',
  dated_tree_fig,
  width = 16,
  height = 12,
  dpi = 600
)

#############################

data <- read_csv("data/PhanDA_GMSTandCO2_percentiles.csv")

cols_to_convert <- c("AverageAge", "GMST_05", "GMST_50", "GMST_95")
data[cols_to_convert] <- lapply(data[cols_to_convert], as.numeric)

data_trimmed <- data %>%
  filter(AverageAge <= abs(ggplot_build(dated_tree_fig)$layout$panel_params[[1]]$x.range[1]))

data_trimmed$MYA_neg <- -data_trimmed$AverageAge

time_plot <- ggplot(data_trimmed, aes(x = MYA_neg)) +
  geom_rect(
    data = data.frame(
      xmin = c(
        ggplot_build(dated_tree_fig)$layout$panel_params[[1]]$x.range[1],-298.9,-251.902,-201.4,-145,-66,-23
      ),
      xmax = c(-298.9, -251.902, -201.4, -145, -66, -23, -2.58),
      ymin = -Inf,
      ymax = Inf,
      fill = c(
        "#67a599",
        "#f04028",
        "#812b92",
        "#34b2c9",
        "#7fc64e",
        "#fd9a52",
        "#ffe619"
      )
    ),
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = fill
    ),
    inherit.aes = FALSE,
    alpha = 0.15
  ) +
  scale_fill_identity() +
  geom_ribbon(aes(ymin = GMST_05, ymax = GMST_95),
              fill = "grey80",
              alpha = 0.5) +
  geom_line(aes(y = GMST_50), size = 1.2) +
  scale_x_continuous(
    limits = c(
      ggplot_build(dated_tree_fig)$layout$panel_params[[1]]$x.range[1],
      ggplot_build(dated_tree_fig)$layout$panel_params[[1]]$x.range[2]
    ),
    breaks = c(-2.58, -23, -66, -145, -201.4, -251.902, -298.9, -316.85),
    labels = c("2.58", "23", "66", "145", "201.4", "251.9", "298.9", "316.85"),
    expand = c(0, 0),
    sec.axis = dup_axis(labels = NULL, name = NULL)
  ) +
  scale_y_continuous(breaks = seq(0, 40, by = 5)) +
  labs(x = "Million Years Ago", y = "Global Mean Surface Temperature (°C)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.ticks.x = element_line(linetype = "dashed"),
    axis.ticks.y = element_line(),
    # Add Y axis ticks
    axis.ticks.length = unit(5, "pt"),
    axis.line.x.bottom = element_blank(),
    axis.line.x.top = element_blank(),
    axis.line.y.left = element_line(),
    # Add Y axis line
    axis.text.x.top = element_blank(),
    axis.title.x.top = element_blank(),
    axis.text.x.bottom = element_text(margin = margin(t = 5)),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  annotate(
    "segment",
    x = ggplot_build(dated_tree_fig)$layout$panel_params[[1]]$x.range[1],
    xend = 0,
    y = Inf,
    yend = Inf,
    color = "black",
    size = 1
  ) +
  annotate(
    "segment",
    x = ggplot_build(dated_tree_fig)$layout$panel_params[[1]]$x.range[1],
    xend = 0,
    y = -Inf,
    yend = -Inf,
    color = "black",
    size = 1
  ) +
  geom_vline(
    xintercept = c(-2.58, -23, -66, -145, -201.4, -251.902),
    linetype = "dashed",
    color = "black",
    alpha = 0.15
  ) +
  geom_text_repel(
    data = data.frame(
      MYA_neg = -28.22,
      GMST = approx(data_trimmed$MYA_neg, data_trimmed$GMST_50, xout = -28.22)$y
    ),
    aes(
      x = MYA_neg,
      y = GMST,
      label = sprintf("%.2f°C", GMST)
    ),
    size = 4,
    fontface = "bold",
    nudge_y = 6,
    direction = "y",
    segment.color = "black",
    segment.size = 0.3,
    max.overlaps = Inf
  ) +
  geom_point(
    data = data.frame(
      MYA_neg = -28.22,
      GMST = approx(data_trimmed$MYA_neg, data_trimmed$GMST_50, xout = -28.22)$y
    ),
    aes(x = MYA_neg, y = GMST),
    shape = 21,
    fill = "yellow",
    color = "black",
    size = 5,
    stroke = 1
  ) +
  geom_text_repel(
    data = data.frame(
      MYA_neg = -140.93,
      GMST = approx(data_trimmed$MYA_neg, data_trimmed$GMST_50, xout = -140.93)$y
    ),
    aes(
      x = MYA_neg,
      y = GMST,
      label = sprintf("%.2f°C", GMST)
    ),
    size = 4,
    fontface = "bold",
    nudge_y = 12,
    direction = "y",
    segment.color = "black",
    segment.size = 0.3,
    max.overlaps = Inf
  ) +
  geom_point(
    data = data.frame(
      MYA_neg = -140.93,
      GMST = approx(data_trimmed$MYA_neg, data_trimmed$GMST_50, xout = -140.93)$y
    ),
    aes(x = MYA_neg, y = GMST),
    shape = 21,
    fill = "red",
    color = "black",
    size = 5,
    stroke = 1
  )

TEST <- (dated_tree_fig / time_plot)  +
  plot_layout(guides = 'collect', heights = c(12, 5))

ggsave(
  "imgs/dated_wombo_combo_tree.png",
  plot = TEST,
  width = 17,
  height = 17,
  dpi = 600
)
