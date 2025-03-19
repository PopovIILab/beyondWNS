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
  patchwork
)

####################
# METADATA PARSING #
####################

# Load `metadata`

metadata <- read.table('metadata/metadata.tsv', sep = '\t', header = T)

# `Year` column

meta.year <- as.data.frame(metadata[, 'Year'])
colnames(meta.year) <- 'Year'
rownames(meta.year) <- metadata$Name
meta.year$Year[meta.year$Year == "ND"] <- NA

# `Host` column

meta.host <- as.data.frame(metadata[, 'Host'])
colnames(meta.host) <- 'Host'
rownames(meta.host) <- metadata$Name
meta.host$Host[meta.host$Host == "ND"] <- NA

# `Country` column

meta.loc <- as.data.frame(metadata[, 'Country'])
colnames(meta.loc) <- 'Country'
rownames(meta.loc) <- metadata$Name
meta.loc$Country[meta.loc$Country == "ND"] <- NA
meta.loc$Country[meta.loc$Country == "USA: New York, Williams Hotel"] <- "USA"
meta.loc$Country[meta.loc$Country == "Switzerland: Boedmerenwald"] <- "Switzerland"
colnames(meta.loc) <- c("Location")


##############
# NAD4L TREE #
##############

# Read the tree file

nad4l_tree <- read.tree("pangenome/tree/nad4l_ufb_ready.treefile")

# Midpoint root the tree

midpoint.root(nad4l_tree)

# Draft tree

nad4l_tree_fig <- ggtree(nad4l_tree) %<+% metadata +
  xlim(0, 0.1) +
  geom_hilight(
    mapping = aes(subset = node %in% c(37), fill = S),
    fill = "steelblue",
    alpha = .6,
    extend = 1
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
  scale_color_identity() +
  geom_treescale(x = 0, y = -0.25, width = 0.05)

# Onehot encode bootstrap values (<70 = 0; >70 = 1)

nad4l_tree_boot <- nad4l_tree_fig$data
nad4l_tree_boot <- nad4l_tree_boot[!nad4l_tree_boot$isTip, ]
nad4l_tree_boot$label <- as.numeric(nad4l_tree_boot$label)
nad4l_tree_boot$bootstrap <- '0'
nad4l_tree_boot$bootstrap[nad4l_tree_boot$label >= 70] <- '1'
nad4l_tree_boot$bootstrap[is.na(nad4l_tree_boot$label)] <- '1'

# Add bootstrap values to the tree (black branches = bootstrap >70; grey branches = bootstrap <70)

nad4l_tree_fig <- nad4l_tree_fig + new_scale_color() +
  geom_tree(data = nad4l_tree_boot, aes(color = bootstrap == '1')) +
  scale_color_manual(name = 'Bootstrap',
                     values = setNames(c("black", "grey"), c(T, F)),
                     guide = "none")

ggsave(
  'imgs/nad4l.png',
  nad4l_tree_fig,
  width = 11,
  height = 8,
  dpi = 600
)

#############
# COX2 TREE #
#############

# Read the tree file

cox2_tree <- read.tree("pangenome/tree/cox2_ufb_ready.treefile")

# Midpoint root the tree

midpoint.root(cox2_tree)

# Draft tree

cox2_tree_fig <- ggtree(cox2_tree) %<+% metadata +
  xlim(0, 0.1) +
  geom_hilight(
    mapping = aes(subset = node %in% c(33), fill = S),
    fill = "steelblue",
    alpha = .6,
    extend = 1
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
  scale_color_identity() +
  geom_treescale(x = 0, y = -0.25, width = 0.05)

# Onehot encode bootstrap values (<70 = 0; >70 = 1)

cox2_tree_boot <- cox2_tree_fig$data
cox2_tree_boot <- cox2_tree_boot[!cox2_tree_boot$isTip, ]
cox2_tree_boot$label <- as.numeric(cox2_tree_boot$label)
cox2_tree_boot$bootstrap <- '0'
cox2_tree_boot$bootstrap[cox2_tree_boot$label >= 70] <- '1'
cox2_tree_boot$bootstrap[is.na(cox2_tree_boot$label)] <- '1'

# Add bootstrap values to the tree (black branches = bootstrap >70; grey branches = bootstrap <70)

cox2_tree_fig <- cox2_tree_fig + new_scale_color() +
  geom_tree(data = cox2_tree_boot, aes(color = bootstrap == '1')) +
  scale_color_manual(name = 'Bootstrap',
                     values = setNames(c("black", "grey"), c(T, F)),
                     guide = "none")

ggsave(
  'imgs/cox2.png',
  cox2_tree_fig,
  width = 11,
  height = 8,
  dpi = 600
)

############
# COB TREE #
############

# Read the tree file

cob_tree <- read.tree("pangenome/tree/cob_ufb_ready.treefile")

# Midpoint root the tree

midpoint.root(cob_tree)

# Draft tree

cob_tree_fig <- ggtree(cob_tree) %<+% metadata +
  xlim(0, 0.1) +
  geom_hilight(
    mapping = aes(subset = node %in% c(29), fill = S),
    fill = "steelblue",
    alpha = .6,
    extend = 1
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
  scale_color_identity() +
  geom_treescale(x = 0, y = -0.25, width = 0.05)

# Onehot encode bootstrap values (<70 = 0; >70 = 1)

cob_tree_boot <- cob_tree_fig$data
cob_tree_boot <- cob_tree_boot[!cob_tree_boot$isTip, ]
cob_tree_boot$label <- as.numeric(cob_tree_boot$label)
cob_tree_boot$bootstrap <- '0'
cob_tree_boot$bootstrap[cob_tree_boot$label >= 70] <- '1'
cob_tree_boot$bootstrap[is.na(cob_tree_boot$label)] <- '1'

# Add bootstrap values to the tree (black branches = bootstrap >70; grey branches = bootstrap <70)

cob_tree_fig <- cob_tree_fig + new_scale_color() +
  geom_tree(data = cob_tree_boot, aes(color = bootstrap == '1')) +
  scale_color_manual(name = 'Bootstrap',
                     values = setNames(c("black", "grey"), c(T, F)),
                     guide = "none")

ggsave(
  'imgs/cob.png',
  cob_tree_fig,
  width = 14,
  height = 8,
  dpi = 600
)

##############################################
# COX1 TREE #
##############################################

# Read the tree file

cox1_tree <- read.tree("pangenome/tree/cox1_ufb_ready.treefile")

# Midpoint root the tree

midpoint.root(cox1_tree)

# Draft tree

cox1_tree_fig <- ggtree(cox1_tree) %<+% metadata +
  xlim(0, 0.2) +
  geom_hilight(
    mapping = aes(subset = node %in% c(28), fill = S),
    fill = "steelblue",
    alpha = .6,
    extend = 1
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
  scale_color_identity() +
  geom_treescale(x = 0, y = -0.25, width = 0.05)

# Onehot encode bootstrap values (<70 = 0; >70 = 1)

cox1_tree_boot <- cox1_tree_fig$data
cox1_tree_boot <- cox1_tree_boot[!cox1_tree_boot$isTip, ]
cox1_tree_boot$label <- as.numeric(cox1_tree_boot$label)
cox1_tree_boot$bootstrap <- '0'
cox1_tree_boot$bootstrap[cox1_tree_boot$label >= 70] <- '1'
cox1_tree_boot$bootstrap[is.na(cox1_tree_boot$label)] <- '1'

# Add bootstrap values to the tree (black branches = bootstrap >70; grey branches = bootstrap <70)

cox1_tree_fig <- cox1_tree_fig + new_scale_color() +
  geom_tree(data = cox1_tree_boot, aes(color = bootstrap == '1')) +
  scale_color_manual(name = 'Bootstrap',
                     values = setNames(c("black", "grey"), c(T, F)),
                     guide = "none")

ggsave(
  'imgs/cox1.png',
  cox1_tree_fig,
  width = 11,
  height = 8,
  dpi = 600
)

##############
# ALL TREE #
##############

# Read the tree file

ALL_tree <- read.tree("phylogenomics/tree/All/All_ready.treefile")

# Midpoint root the tree

midpoint.root(ALL_tree)

# Draft tree

ALL_tree_fig <- ggtree(ALL_tree) %<+% metadata +
  xlim(0, 0.3) +
  geom_hilight(
    mapping = aes(subset = node %in% c(29), fill = S),
    fill = "steelblue",
    alpha = .6,
    extend = 1
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
  scale_color_identity() +
  geom_treescale(x = 0, y = -0.25, width = 0.05)

# Onehot encode bootstrap values (<70 = 0; >70 = 1)

ALL_tree_boot <- ALL_tree_fig$data
ALL_tree_boot <- ALL_tree_boot[!ALL_tree_boot$isTip, ]
ALL_tree_boot$label <- as.numeric(ALL_tree_boot$label)
ALL_tree_boot$bootstrap <- '0'
ALL_tree_boot$bootstrap[ALL_tree_boot$label >= 70] <- '1'
ALL_tree_boot$bootstrap[is.na(ALL_tree_boot$label)] <- '1'

# Add bootstrap values to the tree (black branches = bootstrap >70; grey branches = bootstrap <70)

ALL_tree_fig <- ALL_tree_fig + new_scale_color() +
  geom_tree(data = ALL_tree_boot, aes(color = bootstrap == '1')) +
  scale_color_manual(name = 'Bootstrap',
                     values = setNames(c("black", "grey"), c(T, F)),
                     guide = "none")

ggsave(
  'imgs/ALL.png',
  ALL_tree_fig,
  width = 11,
  height = 8,
  dpi = 600
)

######################
###### COMBINED ######
######################

conserved <- (nad4l_tree_fig + cox2_tree_fig) / (cob_tree_fig + cox1_tree_fig) + plot_annotation(tag_levels = list(c("A", "B", "C", "D")))
ggsave(
  "imgs/conserved_tree.png",
  plot = conserved,
  width = 28,
  height = 16,
  dpi = 600
)

everything <- ALL_tree_fig + (nad4l_tree_fig + cox2_tree_fig) / (cob_tree_fig + cox1_tree_fig) + plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E")))

everything <- ALL_tree_fig + 
  (nad4l_tree_fig + cox2_tree_fig) / 
  (cob_tree_fig + cox1_tree_fig) + 
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E"))) +
  plot_layout(heights = c(12, 1), widths = c(0.6, 1.6))  # Make second row wider

ggsave(
  "imgs/everything_tree.png",
  plot = everything,
  width = 38,
  height = 16,
  dpi = 600
)
