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

metadata <- read.table('../01_PanPhylo_analysis/metadata/metadata.tsv', sep = '\t', header = T)

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
# DATED TREE #
##############

# Read the tree file

dated_tree <- read.beast("data/dating_super_tree_ready.tree")

# Midpoint root the tree

#midpoint.root(dated_tree)

# Draft tree

dated_tree_fig <- ggtree(dated_tree) %<+% metadata +
  xlim(0, 150) +
  geom_hilight(
    mapping = aes(subset = node %in% c(34), fill = S),
    fill = "steelblue",
    alpha = .6,
    extend = 40
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
  geom_treescale(x = 0, y = -0.25, width = 10) +
  
  # Add node ages as text labels, excluding NA and 0 values
  geom_text2(aes(subset = !is.na(height_median) & abs(height_median) > 1e-5, 
                 label = round(height_median, 2)), 
             hjust = -0.5, vjust = 0.5, size = 4, color = "black", fontface = "bold")

# Onehot encode bootstrap values (<70 = 0; >70 = 1)

dated_tree_posterior <- dated_tree_fig$data 
dated_tree_posterior <- dated_tree_posterior[!dated_tree_posterior$isTip, ]
dated_tree_posterior$posterior <- as.numeric(dated_tree_posterior$posterior)
dated_tree_posterior$posterior_vis <- '0'
dated_tree_posterior$posterior_vis[dated_tree_posterior$posterior >= 0.9] <- '1'
dated_tree_posterior$posterior_vis[is.na(dated_tree_posterior$posterior)] <- '1'

# Add bootstrap values to the tree (black branches = bootstrap >70; grey branches = bootstrap <70)

dated_tree_fig <- dated_tree_fig + new_scale_color() +
  geom_tree(data = dated_tree_posterior, aes(color = posterior_vis == '1')) +
  scale_color_manual(name = 'Posterior',
                     values = setNames(c("black", "grey"), c(T, F)),
                     guide = "none")

ggsave(
  'imgs/dated_super_tree.png',
  dated_tree_fig,
  width = 16,
  height = 12,
  dpi = 600
)

