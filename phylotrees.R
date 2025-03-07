library(tidyverse)
library(ggtree)
library(tidytree)
library(patchwork)
# detach("package:phyloseq", unload = T)

# Function to get the MRCA node in a tree of specified taxonomy level
# Best used with ICTV confirmed species
get_mrca <- function(tree, meta, group, subset = NULL) {
  treeplot <- ggtree(tree)
  tipdata <- subset(treeplot$data, is.element(treeplot$data$label, meta$label))

  ugroup <- names(table(tipdata[[group]]))


  mat <- matrix(nrow = length(ugroup), ncol = 2)
  for (i in 1:length(ugroup)) {
    labels <- tipdata[tipdata[[group]] == ugroup[i], ]$label
    labels <- labels[!is.na(labels) & !is.element(labels, subset)]
    # filtered_nodes <- tipdata[tipdata[[group]] == ugroup[i],]$node
    # filtered_nodes <- filtered_nodes[!is.na(filtered_nodes)]

    mat[i, ] <- c(tidytree::MRCA(tree@phylo, .node1 = labels), as.character(ugroup[i]))

    # mat[i,] <- c(ape::getMRCA(tree@phylo, filtered_nodes), as.character(ugroup[i]))
  }
  df <- as.data.frame(mat)
  colnames(df) <- c("node", "id")
  # df <- df[!is.na(as.numeric(as.character(df$node))),]
  rownames(df) <- NULL
  df <- transform(df, node = as.numeric(node))
  return(df)
}

# function to reduce legend size of plots
addSmallLegend <- function(myPlot, pointSize = 1, textSize = 2, spaceLegend = 0.1) {
  myPlot +
    guides(
      shape = "none", # guide_legend(override.aes = list(size = pointSize-.5)),
      color = guide_legend(override.aes = list(size = pointSize - .5, shape = 15))
    ) +
    theme(
      legend.title = element_text(size = textSize, face = "bold"),
      legend.text = element_text(size = textSize),
      legend.key.size = unit(spaceLegend, "lines"),
      legend.text.align = 0,
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
}

# Create a function to find matches
find_matches <- function(vmr_value) {
  if (!is.na(vmr_value)) {
    any(str_detect(gsub("\\.[0-9]$", "", all_meta$Nucleotide), vmr_value))
  } else {
    FALSE
  }
}

# Define tree theme
treetheme <- function() {
  theme_tree(bgcolor = "NA", plot.margin = unit(c(0, 0, 0, 0), "mm"))
}

# Create metadata for trees
vmr <- readxl::read_excel("data/VMR_MSL39_v1.xlsx")

## all_meta <- do.call(rbind, lapply(list.files(path=here::here("data/conftrees/metadata/"), pattern = "*.csv", full.names = T), read.csv))
# all_meta <- read_csv("data/conftrees/metadata/all_meta.csv")
#
## Apply find_matches function to each element in vmr$`Virus GENBANK accession`
# matches_gb <- sapply(vmr$`Virus GENBANK accession`, find_matches)
# matches_refseq <- sapply(vmr$`Virus REFSEQ accession`, find_matches)
#
# ictv_matches_gb <- vmr[matches_gb,]
# ictv_matches_refseq <- vmr[matches_refseq,]
#
# ictv_matches <- rbind(ictv_matches_gb, ictv_matches_refseq)
#
# host_species <- read_tsv("data/conftrees/metadata/host_species2.tsv") |>
#  unique()
#
# all_meta <- left_join(all_meta, host_species, by=join_by("Host" == "Original Species"))
#
# ictv_matches_long <- ictv_matches |>
#  pivot_longer(cols=c(`Virus REFSEQ accession`, `Virus GENBANK accession`), names_to = "DB", values_to = "Accession")
#
# all_meta <- all_meta |>
#  left_join(ictv_matches_long, by=join_by("Nucleotide" == "Accession"), suffix = c(".NCBI", ".ICTV"))

# ICTV metadata
ictv_meta <- do.call(rbind, lapply(
  list.files(
    path = here::here("data/ictv_trees/"),
    pattern = "\\.csv$",
    full.names = TRUE
  )[!grepl("_host\\.csv$", list.files(
    path = here::here("data/ictv_trees/"),
    pattern = "\\.csv$",
    full.names = FALSE
  ))],
  read.csv
))

ictv_meta <- ictv_meta |>
  left_join(
    vmr |>
      separate_rows(`Virus GENBANK accession`, sep = " ") |>
      mutate(`Virus GENBANK accession` = gsub(";", "", `Virus GENBANK accession`)),
    by = join_by("Nucleotide" == "Virus GENBANK accession"),
    suffix = c(".NCBI", ".ICTV")
  )


# Read the two Host source files
host_durna <- read_csv(here::here("data/ictv_trees/blast_durnavirales_host.csv"))
host_ghabri <- read_csv(here::here("data/ictv_trees/blast_ghabrivirales_host.csv"))

# Combine the host source data
host_data <- bind_rows(host_durna, host_ghabri)

# Merge host data with ictv_meta, filling missing "Host source"
ictv_meta <- ictv_meta %>%
  left_join(host_data, by = "Accession", suffix = c("", "_new")) %>%
  mutate(`Host source` = ifelse(is.na(`Host source`), `Host source_new`, `Host source`)) %>%
  select(-`Host source_new`, -Host_new) %>% # Remove the extra column
  mutate(`Host source` = gsub(" \\(S\\)", "", `Host source`)) # Clean up

# Set colors for Host
# cols <- viridis::viridis(length(unique(all_meta$Broad_category))-1)
# cols <- viridis::viridis(length(unique(ictv_meta$`Host source`)))

# give colors for 15 categories
cols <- ggthemes::tableau_color_pal(palette = "Tableau 20")(16)

# Assuming tree@data$Category is your vector
unique_categories <- unique(ictv_meta$`Host source`)

# Remove NA and "Blackflies"
filtered_categories <- unique_categories[!is.na(unique_categories)]

# Display the result
names(cols) <- filtered_categories[order(filtered_categories)]

cols <- c(cols, "Blackflies" = "black")

# Start for each tree

naked_tree <- function(tree, clade, title) {
  p <- ggtree(tree) +
    geom_tiplab(aes(label = "", color = Source, subset = Source != "Blackflies"), linesize = .2, align = T, show.legend = F, linetype = 1) +
    geom_tiplab(aes(label = "", color = Source, subset = Source == "Blackflies"), linesize = .5, align = T, show.legend = F, linetype = 1) +
    scale_color_manual(values = c("ICTV" = "darkgrey", "Blackflies" = "red"), na.translate = F) +
    # scale_x_continuous(expand=expansion(mult=c(0, 0.2)))+
    geom_cladelab(
      data = clade, mapping = aes(node = node, label = id), color = "black",
      fontface = "italic", fontsize = 2, barsize = .25, align = T, offset = .01, parse = T
    ) +
    hexpand(0.2) +
    geom_rootedge(.05) +
    coord_cartesian(clip = "off") +
    ggtitle(title) +
    # theme_tree2("transparent")+
    theme(plot.title = element_text(face = "italic", size = 8))

  y_max <- layer_scales(p)$y$range$range[2]
  x_max <- layer_scales(p)$x$range$range[2]
  y_set <- y_max - (y_max * 0.05)

  if (x_max > 1.5) {
    w <- 0.5
  } else {
    w <- 0.1
  }

  p +
    geom_treescale(x = 0, y = y_set, width = w, offset = 1, fontsize = 2)
}

detailed_tree <- function(tree, clade, title) {
  p <- ggtree::ggtree(tree, ladderize = T) +
    geom_tiplab(aes(label = glue::glue("{label} {Organism_Name}", .na = "")), align = F, size = 1, offset = .01) +
    geom_tippoint(aes(color = `Host source`), size = 1) +
    geom_nodelab(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 70), size = 1, hjust = -.01) +
    geom_rootedge(.05) +
    hexpand(0.2) +
    scale_colour_manual(values = cols, na.value = "grey")

  y_max <- layer_scales(p)$y$range$range[2]
  x_max <- layer_scales(p)$x$range$range[2]
  y_set <- y_max - (y_max * 0.05)

  offset <- x_max * 0.15

  if (x_max > 1) {
    w <- 0.5
  } else {
    w <- 0.1
  }


  p2 <- p +
    geom_cladelab(
      data = clade, mapping = aes(node = node, label = id), color = "black",
      fontface = "italic", fontsize = 2, barsize = .25, align = F, offset = offset, parse = T
    )

  addSmallLegend(p2, pointSize = 3, textSize = 6) +
    coord_cartesian(clip = "off") +
    ggtitle(title) +
    geom_treescale(x = 0, y = y_set, width = w, offset = y_max * 0.01, fontsize = 2) +
    theme(
      plot.title = element_text(face = "italic", size = 8),
      legend.position = "bottom"
    )
}

# Picornavirales
tree <- treeio::read.newick("data/ictv_trees/conftrees/picornavirales_conftree.newick") |>
  phytools::midpoint.root()

tree$tip.label <- gsub("\\.[0-9]$", "", tree$tip.label)

combined_tree <- as_tibble(tree) |>
  left_join(ictv_meta, by = join_by("label" == "Accession")) |>
  mutate(Source = ifelse(str_detect(label, "NODE"), "Blackflies", "ICTV"))

tree <- as.treedata(combined_tree) |>
  mutate(`Host source` = ifelse(startsWith(label, "NODE"), "Blackflies", `Host source`))

clade <- get_mrca(tree = tree, meta = combined_tree, group = "Family.ICTV", subset = c("AFQ98017", "AAL06289", "AAQ64627", "AKJ23343"))

p <- ggtree::ggtree(tree, linewidth = .25, ladderize = T) +
  geom_tiplab(aes(label = glue::glue("{label} {Organism_Name}", .na = "")), align = F, size = 1, offset = .01) +
  geom_tippoint(aes(color = `Host source`), size = 1) +
  geom_nodelab(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 70), size = 1, hjust = -.01) +
  geom_rootedge(.05, linewidth = .25) +
  scale_colour_manual(values = cols, na.value = "grey")

# scale Picornaviridae
p2 <- scaleClade(p, 941, .1)

# scale Secoviridae
p2 <- scaleClade(p2, 1028, .1)

# scale Caliciviridae
p2 <- scaleClade(p2, 991, .2)

p3 <- p2 +
  geom_cladelab(
    data = clade, mapping = aes(node = node, label = id), color = "black",
    fontface = "italic", fontsize = 2, barsize = .25, align = F, offset = 1.5, parse = T
  )

p4 <- p3 |>
  collapse(941, "mixed", fill = cols["vertebrates"]) |>
  collapse(1028, "mixed", fill = cols["plants"]) |>
  collapse(991, "mixed", fill = cols["vertebrates"])

addSmallLegend(p4, pointSize = 3, textSize = 6) +
  geom_treescale(y = 50, width = .5, linesize = .25, fontsize = 2) +
  theme(legend.position = "left") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2)))

ggsave("figures/ictv_trees/picorna.pdf", width = 180, height = 215, units = "mm")

picorna_tree <- tree

picorna_clade <- clade |>
  filter(id %in% c("Dicistroviridae", "Iflaviridae", "Noraviridae", "Marnaviridae"))

picorna <- naked_tree(picorna_tree, picorna_clade, "Picornavirales")

# Flaviviridae
tree <- treeio::read.newick("data/ictv_trees/conftrees/flaviviridae_conftree.newick")

tree <- phytools::midpoint.root(tree)

tree$tip.label <- gsub("\\.[0-9]$", "", tree$tip.label)

combined_tree <- as_tibble(tree) |>
  left_join(ictv_meta, by = join_by("label" == "Accession")) |>
  mutate(Source = ifelse(str_detect(label, "NODE"), "Blackflies", "ICTV"))

tree <- as.treedata(combined_tree)

tree <- tree |>
  mutate(`Host source` = ifelse(startsWith(label, "NODE"), "Blackflies", `Host source`))

clade <- get_mrca(tree = tree, meta = combined_tree, group = "Genus.ICTV")

flavi_tree <- tree

detailed_tree(flavi_tree, clade, "Flaviviridae")
ggsave("figures/ictv_trees/flavi.pdf", dpi = 300, width = 180, height = 215, units = "mm")

flavi_clade <- clade |>
  filter(id %in% "Orthoflavivirus")

flavi <- naked_tree(flavi_tree, flavi_clade, "Flaviviridae") +
  geom_tippoint(aes(label = "", color = Source, subset = Source == "Blackflies"), show.legend = F)

# Martellivirales
tree <- treeio::read.newick("data/ictv_trees/conftrees/martellivirales_conftree.newick")

tree <- phytools::midpoint.root(tree)

tree$tip.label <- gsub("\\.[0-9]$", "", tree$tip.label)

combined_tree <- as_tibble(tree) |>
  left_join(ictv_meta, by = join_by("label" == "Accession")) |>
  mutate(Source = ifelse(str_detect(label, "NODE"), "Blackflies", "ICTV"))

tree <- as.treedata(combined_tree)

tree <- tree |>
  mutate(`Host source` = ifelse(startsWith(label, "NODE"), "Blackflies", `Host source`))

clade <- get_mrca(tree = tree, meta = combined_tree, group = "Family.ICTV")

martelli_tree <- tree

detailed_tree(martelli_tree, clade, "Martellivirales")

ggsave("figures/ictv_trees/martelli.pdf", dpi = 300, width = 180, height = 215, units = "mm")

martelli_clade <- clade |>
  filter(id %in% c("Togaviridae", "Closteroviridae", "Mayoviridae", "Virgaviridae", "Bromoviridae", "Kitaviridae", "Endornaviridae"))

martelli <- naked_tree(martelli_tree, martelli_clade, "Martellivirales") +
  geom_tippoint(aes(label = "", color = Source, subset = label == "NODE_A21_length_7442_cov_187.677936_BlackFly31_2"), show.legend = F)

# Durnavirales
tree <- treeio::read.newick("data/ictv_trees/conftrees/durnavirales_conftree.newick")

tree <- phytools::midpoint.root(tree)

tree$tip.label <- gsub("\\.[0-9]$", "", tree$tip.label)

combined_tree <- as_tibble(tree) |>
  left_join(ictv_meta, by = join_by("label" == "Accession")) |>
  mutate(Source = ifelse(str_detect(label, "NODE"), "Blackflies", "ICTV"))

tree <- as.treedata(combined_tree)

tree <- tree |>
  mutate(`Host source` = ifelse(startsWith(label, "NODE"), "Blackflies", `Host source`))

clade <- get_mrca(tree = tree, meta = combined_tree, group = "Family.ICTV")

# Split Partiti clades
get_mrca(tree = tree, meta = combined_tree, group = "Genus.ICTV")

clade <- clade |>
  filter(id != "Partitiviridae") |>
  bind_rows(data.frame(node = c(286, 336, 312, 358), id = c("Gammapartitivirus", "Alphapartitivirus", "Betapartitivirus", "Deltapartitivirus")))

durna_tree <- tree
detailed_tree(durna_tree, clade, "Durnavirales")

ggsave("figures/ictv_trees/durna.pdf", dpi = 300, width = 180, height = 215, units = "mm")

durna_clade <- clade |>
  filter(id %in% c(
    "Gammapartitivirus", "Alphapartitivirus", "Betapartitivirus", "Deltapartitivirus",
    "Curvulaviridae", "Hypoviridae", "Picobirnaviridae"
  ))

durna <- naked_tree(durna_tree, durna_clade, "Durnavirales")

# Ghabrivirales
tree <- treeio::read.newick("data/ictv_trees/conftrees/ghabrivirales_conftree.newick")

tree <- phytools::midpoint.root(tree)

tree$tip.label <- gsub("\\.[0-9]$", "", tree$tip.label)

combined_tree <- as_tibble(tree) |>
  left_join(ictv_meta, by = join_by("label" == "Accession")) |>
  mutate(Source = ifelse(str_detect(label, "NODE"), "Blackflies", "ICTV"))

tree <- as.treedata(combined_tree)

tree <- tree |>
  mutate(`Host source` = ifelse(startsWith(label, "NODE"), "Blackflies", `Host source`))

clade <- get_mrca(tree = tree, meta = combined_tree, group = "Family.ICTV")

ghabri_tree <- tree

detailed_tree(ghabri_tree, clade, "Ghabrivirales")

ggsave("figures/ictv_trees/ghabri.pdf", dpi = 300, width = 180, height = 215, units = "mm")

ghabri_clade <- clade |>
  filter(id %in% c("Chrysoviridae", "Pseudototiviridae", "Orthototiviridae"))

ghabri <- naked_tree(ghabri_tree, ghabri_clade, "Ghabrivirales")

# Rhabdoviridae
tree <- treeio::read.newick("data/ictv_trees/conftrees/rhabdoviridae_conftree.newick")

tree <- phytools::midpoint.root(tree)

tree$tip.label <- gsub("\\.[0-9]$", "", tree$tip.label)

combined_tree <- as_tibble(tree) |>
  left_join(ictv_meta, by = join_by("label" == "Accession")) |>
  mutate(Source = ifelse(str_detect(label, "NODE"), "Blackflies", "ICTV"))

tree <- as.treedata(combined_tree)

tree <- tree |>
  mutate(`Host source` = ifelse(startsWith(label, "NODE"), "Blackflies", `Host source`))

clade <- get_mrca(tree = tree, meta = combined_tree, group = "Genus.ICTV")

rhabdo_tree <- tree
detailed_tree(rhabdo_tree, clade, "Rhabdoviridae")

ggsave("figures/ictv_trees/rhabdo.pdf", dpi = 300, width = 180, height = 290, units = "mm")

rhabdo_clade <- clade |>
  filter(id %in% c("Almendravirus", "Lyssavirus", "Amplylivirus"))

rhabdo <- naked_tree(rhabdo_tree, rhabdo_clade, "Rhabdoviridae")

# Phenuiviridae
tree <- treeio::read.newick("data/ictv_trees/conftrees/phenuiviridae_conftree.newick")

tree <- phytools::midpoint.root(tree)

tree$tip.label <- gsub("\\.[0-9]$", "", tree$tip.label)

combined_tree <- as_tibble(tree) |>
  left_join(ictv_meta, by = join_by("label" == "Accession")) |>
  mutate(Source = ifelse(str_detect(label, "NODE"), "Blackflies", "ICTV"))

tree <- as.treedata(combined_tree)

tree <- tree |>
  mutate(`Host source` = ifelse(startsWith(label, "NODE"), "Blackflies", `Host source`))

clade <- get_mrca(tree = tree, meta = combined_tree, group = "Genus.ICTV")

phenui_tree <- tree
detailed_tree(phenui_tree, clade, "Phenuiviridae")

ggsave("figures/ictv_trees/phenui.pdf", dpi = 300, width = 180, height = 215, units = "mm")

phenui_clade <- clade |>
  filter(id %in% c("Laulavirus"))

phenui <- naked_tree(phenui_tree, phenui_clade, "Phenuiviridae")

# Mymonaviridae
tree <- treeio::read.newick("data/ictv_trees/conftrees/mymonaviridae_conftree.newick")

tree <- phytools::midpoint.root(tree)

tree$tip.label <- gsub("\\.[0-9]$", "", tree$tip.label)

combined_tree <- as_tibble(tree) |>
  left_join(ictv_meta, by = join_by("label" == "Accession")) |>
  mutate(Source = ifelse(str_detect(label, "NODE"), "Blackflies", "ICTV"))

tree <- as.treedata(combined_tree)

tree <- tree |>
  mutate(`Host source` = ifelse(startsWith(label, "NODE"), "Blackflies", `Host source`))

clade <- get_mrca(tree = tree, meta = combined_tree, group = "Genus.ICTV") |>
  filter(id != "Botrytimonavirus")

mymona_tree <- tree

detailed_tree(mymona_tree, clade, "Mymonaviridae")

ggsave("figures/ictv_trees/mymona.pdf", dpi = 300, width = 180, height = 215, units = "mm")

mymona_clade <- clade |>
  filter(id %in% c("Penicillimonavirus", "Sclerotimonavirus", "Plasmopamonavirus"))

mymona_clade <- mymona_clade |>
  mutate(id = str_replace(id, "monavirus", "mona-\nvirus"))

mymona <- naked_tree(mymona_tree, mymona_clade, "Mymonaviridae")

# Tymoviridae
tree <- treeio::read.newick("data/ictv_trees/conftrees/tymoviridae_conftree.newick")

tree <- phytools::midpoint.root(tree)

tree$tip.label <- gsub("\\.[0-9]$", "", tree$tip.label)

combined_tree <- as_tibble(tree) |>
  left_join(ictv_meta, by = join_by("label" == "Accession")) |>
  mutate(Source = ifelse(str_detect(label, "NODE"), "Blackflies", "ICTV"))

tree <- as.treedata(combined_tree)

tree <- tree |>
  mutate(`Host source` = ifelse(startsWith(label, "NODE"), "Blackflies", `Host source`))

clade <- get_mrca(tree = tree, meta = combined_tree, group = "Genus.ICTV")

tymo_tree <- tree
detailed_tree(tymo_tree, clade, "Tymoviridae")

ggsave("figures/ictv_trees/tymo.pdf", dpi = 300, width = 180, height = 215, units = "mm")

tymo_clade <- clade |>
  filter(id %in% c("Marafivirus", "Maculavirus"))

tymo <- naked_tree(tymo_tree, tymo_clade, "Tymoviridae") +
  geom_tippoint(aes(label = "", color = Source, subset = Source == "Blackflies"), show.legend = F)

# Orthomyxoviridae
tree <- treeio::read.newick("data/ictv_trees/conftrees/orthomyxoviridae_conftree.newick")

tree <- phytools::midpoint.root(tree)

tree$tip.label <- gsub("\\.[0-9]\\|.*$", "", tree$tip.label)

orthomyxo_meta <- ictv_meta |>
  mutate(Accession = ifelse(Nucleotide %in% tree$tip.label, Nucleotide, Accession))

ortho_hosts <- read_tsv("data/ictv_trees/orthomyxo_host_species.tsv") |>
  unique()

orthomyxo_meta <- orthomyxo_meta |>
  left_join(ortho_hosts) |>
  mutate(`Host source` = ifelse(is.na(`Host source`), Category, `Host source`))

combined_tree <- as_tibble(tree) |>
  left_join(orthomyxo_meta, by = join_by("label" == "Accession")) |>
  mutate(Source = ifelse(str_detect(label, "NODE"), "Blackflies", "ICTV"))

tree <- as.treedata(combined_tree)

tree <- tree |>
  mutate(
    `Host source` = ifelse(startsWith(label, "NODE"), "Blackflies", `Host source`),
    across(where(is.character), ~ na_if(., ""))
  )

clade <- get_mrca(tree = tree, meta = combined_tree, group = "Genus.NCBI")

orthomyxo_tree <- tree
detailed_tree(orthomyxo_tree, clade, "Orthomyxoviridae")

ggsave("figures/ictv_trees/orthomyxo.pdf", dpi = 300, width = 180, height = 215, units = "mm")

orthomyxo_clade <- clade |>
  filter(id %in% "Quaranjavirus")

orthomyxo <- naked_tree(orthomyxo_tree, orthomyxo_clade, "Orthomyxoviridae")

# Parvoviridae
tree <- treeio::read.newick("data/ictv_trees/conftrees/parvoviridae_conftree.newick")

tree <- phytools::midpoint.root(tree)

tree$tip.label <- gsub("\\.[0-9]$", "", tree$tip.label)

combined_tree <- as_tibble(tree) |>
  left_join(ictv_meta, by = join_by("label" == "Accession")) |>
  mutate(Source = ifelse(str_detect(label, "NODE"), "Blackflies", "ICTV"))

tree <- as.treedata(combined_tree)

tree <- tree |>
  mutate(`Host source` = ifelse(startsWith(label, "NODE"), "Blackflies", `Host source`))

clade <- get_mrca(tree = tree, meta = combined_tree, group = "Subfamily", subset = c("WDB01642", "CAA52899"))

parvo_tree <- tree
detailed_tree(parvo_tree, clade, "Parvoviridae")

ggsave("figures/ictv_trees/parvo.pdf", dpi = 300, width = 180, height = 215, units = "mm")

parvo_clade <- clade |>
  filter(id %in% c("Densovirinae"))

parvo <- naked_tree(parvo_tree, parvo_clade, "Parvoviridae")

# Genomoviridae
tree <- treeio::read.newick("data/ictv_trees/conftrees/genomoviridae_conftree.newick")

tree <- phytools::midpoint.root(tree)

tree$tip.label <- gsub("\\.[0-9]$", "", tree$tip.label)

combined_tree <- as_tibble(tree) |>
  left_join(ictv_meta, by = join_by("label" == "Accession")) |>
  mutate(Source = ifelse(str_detect(label, "NODE"), "Blackflies", "ICTV"))

tree <- as.treedata(combined_tree)

tree <- tree |>
  mutate(`Host source` = ifelse(startsWith(label, "NODE"), "Blackflies", `Host source`))

clade <- get_mrca(tree = tree, meta = combined_tree, group = "Genus.ICTV", subset = c("QCX35071"))

genomo_tree <- tree
detailed_tree(genomo_tree, clade, "Genomoviridae")

ggsave("figures/ictv_trees/genomo.pdf", dpi = 300, width = 180, height = 215, units = "mm")

ggtree(genomo_tree) +
  geom_text(aes(label = node), size = 2)

# gemykibi
# node 303, 309, 352

genomo_clade <- clade |>
  filter(id %in% c("Gemycircularvirus", "Gemyvongvirus")) |>
  bind_rows(data.frame(node = c(303, 309, 352), id = c("Gemykibivirus; clade 1", "Gemykibivirus; clade 2", "Gemykibivirus; clade 3")))

genomo <- naked_tree(genomo_tree, genomo_clade, "Genomoviridae")


# All trees
# dsRNA + ssDNA
((parvo | genomo) / (ghabri | durna)) +
  plot_layout(
    heights = c(.5, 1),
    tag_level = "keep"
  ) &
  plot_annotation(tag_levels = "A")
ggsave("figures/ictv_trees/dsRNA_ssDNA.pdf", dpi = 300, height = 215, width = 180, units = "mm")

# ssRNA+

(picorna | (martelli / flavi / tymo) +
  plot_layout(heights = c(1, .5, .5))) +
  plot_layout(tag_level = "keep") &
  plot_annotation(tag_levels = "A")

ggsave("figures/ictv_trees/ssRNA(+).pdf", dpi = 300, height = 215, width = 180, units = "mm")

# ssRNA-
(rhabdo | (phenui / orthomyxo / mymona) +
  plot_layout(heights = c(1, 1, .5))) +
  plot_layout(tag_level = "keep") &
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.margin = margin(0),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
  )

ggsave("figures/ictv_trees/ssRNA(-).pdf", dpi = 300, height = 215, width = 180, units = "mm")

# All
(ghabri / durna) | (picorna / tymo) + plot_layout(heights = c(1, .1)) | (martelli / rhabdo) | (flavi / phenui / orthomyxo / mymona) + plot_layout(heights = c(1, 1, 1, .7))
ggsave("figures/ictv_trees/all.pdf", dpi = 300, height = 215, width = 180, units = "mm")


# Get accessions
vmr |>
  filter(Family == "Parvoviridae" & str_detect(`Genome coverage`, "[Cc]omplete")) |>
  pull(`Virus GENBANK accession`) |>
  as.data.frame() |>
  write_tsv("data/ictv_parvoviridae.tsv", col_names = F)
