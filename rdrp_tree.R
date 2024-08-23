library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(readxl)
library(ggnewscale)
library(ggstar)
library(ggtext)
library(magick)

figwidth_mm <- 164
figheight_mm <- 164

convert_mm_to_pixels <- function(size){
  return(as.integer((size/25.4)*300))
}
  
#input_tree <- treeio::read.newick("RdRP_tree/rooted_blackflies_RdRP_conftree.ladderized.newick", node.label='support')
input_tree <- treeio::read.newick("RdRP_tree/rooted_fasttree_maxcc_stripped_conf.newick", node.label='support')
#input_tree <- treeio::read.newick("RdRP_tree/rooted_palmcore_flanked_conftree.ladderized.newick", node.label='support')

# Create a data frame with tip labels
tips_data <- data.frame(tip = input_tree@phylo$tip.label)

# Add a column to the data frame indicating whether the tip label contains "BlackFly"
tips_data$contains_blackfly <- grepl("BlackFly", tips_data$tip)

# Read RCR90 data
rcr90 <- read_excel("RdRP_tree/mmc1.xlsx")
palmscan_report <- read_tsv("RdRP_tree/palmscan_blackflies_RVMT_rt.tsv")

rcr90_tree <- left_join(tips_data, rcr90, by=c("tip" = "RCR90")) %>% 
  left_join(palmscan_report, by=c("tip" = "Label"))

rcr90_tree <- rcr90_tree %>%
  mutate(
    new_phylum = case_when(
      Order == "Durnavirales" ~ "*Pisuviricota* (dsRNA)",
      Order == "Amarillovirales" ~ "*Kitrinoviricota*",
      startsWith(tip, "rt.") ~ "Reverse transcriptases",
      !is.na(name) ~ glue::glue("*{Phylum}*"),
      startsWith(Family, "f.") ~ "Neri *et al.* Cell (2022)",
      grepl("BlackFly", tip) ~ "Blackflies",
      TRUE ~ glue::glue("*{Phylum}*")  # Default condition
    )
  )

rcr90_tree$new_order <- if_else(startsWith(rcr90_tree$Order, "o."), NA, rcr90_tree$Order)

unique_orders <- unique(rcr90_tree$new_order)

library(RColorBrewer)
n <- 27
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col27 <- sample(col_vector, n)

names(col27) <- unique_orders

unique_phyla <- sort(as.character(unique(rcr90_tree$new_phylum)))

phyla <- subset(unique_phyla, !unique_phyla %in% c("Blackflies", "Neri *et al.* Cell (2022)"))

phyl_col <- c(brewer.pal(7, "Paired"), "lightgrey", "black")
names(phyl_col) <- c(phyla, "Neri *et al.* Cell (2022)", "Blackflies")

rcr90_tree$new_phylum <- factor(rcr90_tree$new_phylum, levels = names(phyl_col))


tree <- left_join(input_tree, rcr90_tree, by=c("label" = "tip"))
tree_data <- as_tibble(tree)

# Set color vector
cols <- c("TRUE"="red", "FALSE"="black")


base_tree <- ggtree(tree, layout = "circular", branch.length ="none", size=.1)


annotated_tree <- base_tree +
  geom_tippoint(aes(shape=contains_blackfly, color=contains_blackfly), 
                size=1, alpha=0.5)+
  scale_shape_manual(values=c("TRUE"=16, "FALSE"=NA), guide="none")+
  scale_color_manual(values=cols, guide="none")+
  #geom_fruit(
  #  geom=geom_star,
  #  mapping=aes(starshape=contains_blackfly, fill=contains_blackfly),
  #  position="identity",
  #  starstroke=0.1, 
  #  show.legend=F, 
  #  offset = 0.03,
  #  size=3
  #) +
  #scale_starshape_manual(
  #  values=c("TRUE"=1, "FALSE"=NA))+
  #new_scale_fill() +
  #geom_nodepoint(aes(fill=support), shape=21, color="transparent")+
  #new_scale_fill() + 
  geom_fruit(aes(fill=new_phylum), 
             geom=geom_tile, 
             pwidth=2, 
             alpha=1, 
             offset = .03)+
  scale_fill_manual(values=phyl_col, name="Phylum",
                    guide = guide_legend(order = 1, 
                                         override.aes = c(starshape=NA)))+
  #new_scale_fill() + 
  #geom_fruit(aes(fill=new_order), 
  #           geom=geom_tile, 
  #           pwidth=2, 
  #           alpha=1, 
  #           offset = .01)+
  #scale_fill_manual(values=col27)+
  #new_scale_fill() + 
  #geom_fruit(aes(fill=ABC), 
  #           geom=geom_tile, 
  #           pwidth=2, 
  #           alpha=1, 
  #           offset = .03)+
  #scale_fill_brewer(palette = "Set2", name="Catalytic domain order",
  #                  guide = guide_legend(order = 2, ))+
  new_scale_fill() + 
  new_scale("shape")+
  geom_fruit(
    geom=geom_star,
    mapping=aes(starshape=ABC, fill=ABC),
    starstroke=0,
    color=NA,
    show.legend=T, 
    offset = 0.03,
    alpha=0.5,
    size=1.5
  ) +
  scale_starshape_manual(
    values=c("CAB"=1, "ABC"=NA), name="Catalytic domain order",
    guide = guide_legend(order = 2))+
  scale_fill_brewer(palette = "Set2", name="Catalytic domain order",
                    guide = guide_legend(order = 2, override.aes = list(alpha = 1)))


#fig <- image_graph(width = convert_mm_to_pixels(figwidth_mm), height = convert_mm_to_pixels(figheight_mm), res=300)
theme_tree <- annotated_tree +
  theme(legend.title = element_text(size = 5), 
        legend.text  = element_markdown(size = 4),
        legend.direction = "vertical",
        legend.spacing = unit(0, "mm"),
        legend.box = "vertical",
        legend.text.align = 0,
        legend.key.size = unit(2, "mm"),
        legend.background = element_blank(),
        legend.position = c(0.15, 0.17),
        plot.margin =  margin(t = 0,
                              r = 0,
                              b = 0,
                              l = 0))
#dev.off()

#ggsave("figures/RdRP_tree.pdf", width = 20, height = 20, dpi=300)
#ggsave("figures/RdRP_tree_stripped.pdf", width = 20, height = 20, dpi=300)
ggsave(plot = theme_tree, filename = "figures/RdRP_tree_palmcore_abc.pdf", 
       width = figwidth_mm, height = figheight_mm, units = "mm" , dpi=300)

flavi_mrca <- MRCA(tree_data, tree_data %>% filter(Family=="Flaviviridae" & node != 19416) %>% select(node) %>% pull()) %>% pull(node)
rhabdo_mrca <- MRCA(tree_data, tree_data %>% filter(Family=="Rhabdoviridae") %>% pull(node)) %>% pull(node)

label_tree <- theme_tree +
  geom_cladelab(geom="text", node = flavi_mrca, label = "Flaviviridae", fontface="italic", fontsize=2, 
                #fontcolour="gold", barcolor="gold", offset=5.5,
                fontcolour="black", barcolor="black", offset=5.5,
                horizontal=FALSE, angle="auto",vjust=1, hjust=0.5)+
  geom_cladelab(geom="text", node = rhabdo_mrca, label = "Rhabdoviridae", fontface="italic", fontsize=2, 
                #fontcolour="orange", barcolor="orange", offset=6,
                fontcolour="black", barcolor="black", offset=6,
                horizontal=FALSE, angle="auto",vjust=-.75, hjust=0.5)

#ggsave(plot = label_tree, filename = "figures/RdRP_tree_palmcore_abc_labels.pdf", 
#       width = figwidth_mm, height = figheight_mm, units = "mm" , dpi=300)

ggsave(plot = label_tree, filename = "figures/tiff/figure5.tiff", 
       width = figwidth_mm, height = figheight_mm, units = "mm" , dpi=300)

support_tree <- base_tree +
  geom_nodepoint(aes(color=support), shape=16, size=1)+
  scale_color_gradientn(colors = c("red", "yellow2", "green3"))

ggsave(plot = support_tree, filename = "figures/RdRP_tree_palmcore_abc_support.pdf", 
       width = figwidth_mm, height = figheight_mm, units = "mm" , dpi=300)

fig <- image_read_pdf("figures/RdRP_tree_stripped_publication.pdf", density = 300)
trimmed_fig <- image_trim(fig)
image_write(trimmed_fig, path = "figures/RdRP_tree_stripped_publication_magick.pdf", 
            format = "pdf", density=300, flatten = F)


picobirna_tree_subset <- tidytree::tree_subset(tree, "Rv4_193103", levels_back = 30)

ggtree(picobirna_tree_subset, layout = "fan", open.angle = 10, branch.length ="none", size=.1)+
  geom_tippoint(aes(shape=contains_blackfly, color=contains_blackfly), 
                size=2,  alpha=0.2)+
  scale_shape_manual(values=c("TRUE"=16, "FALSE"=NA), guide = "none")+
  scale_color_manual(values=cols, guide = "none")+
  geom_fruit(aes(fill=new_phylum), 
           geom=geom_tile, 
           pwidth=2, 
           alpha=1, 
           offset = .03)+
  scale_fill_manual(values=phyl_col, 
                    guide = guide_legend(order = 1, override.aes = c(starshape=NA)))+
  new_scale_fill() + 
  new_scale("shape")+
  new_scale_color()+
  geom_fruit(
    geom=geom_star,
    mapping=aes(starshape=ABC, fill=ABC),
    starstroke=0.01, 
    show.legend=T, 
    offset = 0.03,
    alpha=0.5,
    size=1.5
  ) +
  scale_starshape_manual(
    values=c("CAB"=1, "ABC"=NA), name="Catalytic domain order",
    guide = guide_legend(order = 2))+
  scale_fill_brewer(palette = "Set2", name="Catalytic domain order",
                    guide = guide_legend(order = 2, override.aes = list(alpha = 1)))+
  theme(legend.title = element_text(size = 5), 
        legend.text  = element_markdown(size = 4),
        legend.direction = "vertical",
        legend.spacing = unit(0, "mm"),
        legend.box = "vertical",
        legend.text.align = 0,
        legend.key.size = unit(2, "mm"),
        legend.background = element_blank(),
        legend.position = c(0.15, 0.17),
        plot.margin =  margin(t = 0,
                              r = 0,
                              b = 0,
                              l = 0))

kitrino_tree_subset <- tidytree::tree_subset(tree, "rv20_v2_2758", levels_back = 35)

ggtree(kitrino_tree_subset, layout = "fan", open.angle = 10, branch.length ="none", size=.1)+
  #geom_tippoint(aes(shape=contains_blackfly, color=contains_blackfly), size=2, show.legend = F)+
  geom_fruit(
    geom=geom_star,
    mapping=aes(starshape=contains_blackfly, fill=contains_blackfly),
    position="identity",
    starstroke=0.1, 
    show.legend=F
  ) +
  #scale_shape_manual(values=c("TRUE"=16, "FALSE"=NA))+
  scale_fill_manual(values=cols)+
  scale_starshape_manual(
    values=c("TRUE"=1, "FALSE"=NA),
    guide=guide_legend(
      keywidth=0.5,
      keyheight=0.5,
      order=1
    ))+
  #geom_nodepoint(aes(fill=support), shape=21, color="transparent")+
  new_scale_fill() + 
  geom_fruit(aes(fill=new_phylum), geom=geom_tile, pwidth=2, alpha=1, offset = .08)+
  scale_fill_brewer(palette="Paired")#+
  #new_scale_fill() + 
  #geom_fruit(aes(fill=Order), geom=geom_tile, pwidth=2, alpha=1, offset = .08)+
  #scale_fill_brewer(palette="Set1")
#geom_nodelab(geom='label', aes(label=support, subset=support > 99))
#scale_shape_manual(values=named_vector)


pisu_ds_tree_subset <- tidytree::tree_subset(tree, "Rv4_318210", levels_back = 27)

ggtree(pisu_ds_tree_subset, layout = "fan", open.angle = 10, branch.length ="none", size=.1)+
  #geom_tippoint(aes(shape=contains_blackfly, color=contains_blackfly), size=2, show.legend = F)+
  geom_fruit(
    geom=geom_star,
    mapping=aes(starshape=contains_blackfly, fill=contains_blackfly),
    position="identity",
    starstroke=0.1, 
    show.legend=F
  ) +
  #scale_shape_manual(values=c("TRUE"=16, "FALSE"=NA))+
  scale_fill_manual(values=cols)+
  scale_starshape_manual(
    values=c("TRUE"=1, "FALSE"=NA),
    guide=guide_legend(
      keywidth=0.5,
      keyheight=0.5,
      order=1
    ))+
  new_scale_fill() + 
  geom_nodepoint(aes(fill=support, subset=support<5), shape=21, color="transparent")+
  new_scale_fill() + 
  geom_fruit(aes(fill=new_phylum), geom=geom_tile, pwidth=2, alpha=1, offset = .08)+
  scale_fill_brewer(palette="Paired")+
  new_scale_fill() + 
  geom_fruit(aes(fill=Order), geom=geom_tile, pwidth=2, alpha=1, offset = .08)+
  scale_fill_brewer(palette="Set1")




p0002_tree_subset <- tidytree::tree_subset(tree, "Rv4_084634", levels_back = 7)

ggtree(p0002_tree_subset, layout = "fan", open.angle = 10, branch.length ="none", size=.1)+
  #geom_tippoint(aes(shape=contains_blackfly, color=contains_blackfly), size=2, show.legend = F)+
  geom_fruit(
    geom=geom_star,
    mapping=aes(starshape=contains_blackfly, fill=contains_blackfly),
    position="identity",
    starstroke=0.1, 
    show.legend=F
  ) +
  #scale_shape_manual(values=c("TRUE"=16, "FALSE"=NA))+
  scale_fill_manual(values=cols)+
  scale_starshape_manual(
    values=c("TRUE"=1, "FALSE"=NA),
    guide=guide_legend(
      keywidth=0.5,
      keyheight=0.5,
      order=1
    ))+
  new_scale_fill() + 
  geom_nodepoint(aes(fill=support, subset=support>80), shape=21, color="transparent")+
  new_scale_fill() + 
  geom_fruit(aes(fill=new_phylum), geom=geom_tile, pwidth=2, alpha=1, offset = .08)+
  scale_fill_brewer(palette="Paired")+
  new_scale_fill() + 
  geom_fruit(aes(fill=new_order), geom=geom_tile, pwidth=2, alpha=1, offset = .08)+
  scale_fill_manual(values=col27)


### Test env
tree <- treeio::read.newick("RdRP_tree/test_tree.newick")

# Create a data frame with tip labels
tips_data <- data.frame(tip = tree$tip.label)

# Add a column to the data frame indicating whether the tip label contains "BlackFly"
tips_data$contains_blackfly <- grepl("BlackFly", tips_data$tip)
tips_data$Order <- c("Picornavirales", "Durnavirales", "Caudovirales", NA, 
                    "Bunyavirales", "Tymovirales", "Megavirales", NA, 
                    "Polycipivirales", "Tymovirales")


tree <- left_join(tree, tips_data, by=c("label" = "tip"))

cols<- c("TRUE"="red", "FALSE"="black")

# Create a named vector with unique instances
unique_orders <- unique(tree@extraInfo$Order)

# Create a named vector with shape numbers
shape_numbers <- ifelse(unique_orders == "Durnavirales", 16, NA)

# Combine names and values to create the named vector
named_vector <- setNames(shape_numbers, unique_orders)

named_vector

ggtree(tree, layout = "circular", branch.length = "none")+
  #geom_tippoint(aes(color=contains_blackfly, shape=contains_blackfly), size=2, show.legend = F)+
  geom_tippoint(aes(shape=Order))+
  geom_tiplab()+
  scale_shape_manual(values=shapes)
