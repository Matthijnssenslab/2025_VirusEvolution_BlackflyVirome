library(tidyverse)
library(patchwork)
library(ggnewscale)
library(readxl)
library(ggtext)

vmr <- read_excel("data/VMR_MSL38_v2.xlsx")

#' ## Rhabdoviridae
rhabdo <- vmr %>% 
  filter(Family == "Rhabdoviridae") 

rhabdo %>% 
  pull(`Virus REFSEQ accession`) %>% 
  na.omit() %>% 
  str_remove(".*: ") %>% 
  str_split("; ", simplify = TRUE) %>% 
  as.vector() %>% 
  purrr::discard(. == "") %>% 
  write_lines(file = "data/dinuq/rhabdo/rhabdovirus_accessions.txt")

#' ### RDA
#' 
lytras_rhabdo <- read_excel("~/Downloads/viruses-746103-supplementary.xlsx", sheet = "Rhabdoviridae") %>% 
  select(accession, Host) %>% 
  distinct() %>% 
  rename(acc=accession)

rda_dinuq <- read_tsv("data/dinuq/rhabdo/dinuq_rhabdo_rda2.tsv")

rda_pca <- rda_dinuq %>% 
  column_to_rownames("acc") %>% 
  t() %>% 
  prcomp()

joined_rda_df <- as.data.frame(rda_pca$rotation) %>% 
  rownames_to_column("acc") %>% 
  mutate(acc=str_remove(acc, "\\.[0-9].*"),
         acc=str_remove(acc, "join\\("),
         acc=str_remove(acc, "complement\\(")) %>% 
  full_join(rhabdo, by=join_by(acc==`Virus REFSEQ accession`)) %>% 
  rename(Host=`Host source`) %>% 
  full_join(lytras_rhabdo %>% mutate(acc=str_remove(acc, "\\.[0-9].*")), by="acc") %>% 
  mutate(Host.x = coalesce(Host.x, Host.y),
         Host.x=if_else(Host.x=="vertebrates", "vertebrate", Host.x), 
         `Virus name abbreviation(s)`=case_when(startsWith(acc, "Blackfly") == T ~ acc,
                                                startsWith(acc, "OVRV1") == T ~ acc,
                                                TRUE ~ `Virus name abbreviation(s)`)) %>% 
  filter(!is.na(PC1))

rhabdo_pca_plot <- joined_rda_df %>% 
  filter(Host.x %in% c("vertebrate", "invertebrates", "invertebrates, vertebrates", NA)) %>% 
  ggplot(aes(x=PC1, y=PC2))+
  stat_density_2d(data=. %>% filter(!is.na(Host.x)), aes(fill=Host.x), 
                  color=NA, alpha=.15, geom="polygon", contour_var = "count")+
  geom_point(data=. %>% filter(!is.na(Host.x)), aes(color=Host.x), size=1)+
  scale_fill_viridis_d(name="Host", guide=guide_legend(order=1))+
  scale_color_viridis_d(name="Host", guide=guide_legend(order=1))+
  theme_bw()
rhabdo_pca_plot

rhabdo_pca_plot2 <- rhabdo_pca_plot +
  new_scale_color()+
  geom_point(data=. %>% filter(startsWith(acc, "OVRV1")), aes(x=PC1, y=PC2, color="#ff7f00"), size=3)+
  geom_point(data=. %>% filter(startsWith(acc, "Blackfly")), aes(x=PC1, y=PC2, color="#377eb8"), size=3)+
  geom_point(data=. %>% filter(startsWith(acc, "NC_077192")), aes(x=PC1, y=PC2, color="#a65628"), size=3)+
  scale_color_manual(name = 'Virus of interest',
                     values =c("#ff7f00"="#ff7f00", "#377eb8"="#377eb8", "#a65628"="#a65628"), 
                     labels = c("#ff7f00"='OVRV1', "#377eb8"='Blackfly rhabdoviruses',  "#a65628"="Mundri virus"), 
                     guide=guide_legend(order=2))+
  geom_label_repel(data=. %>% filter(str_starts(acc,"OVRV1|Blackfly|NC_077192") | str_starts(`Virus name abbreviation(s)`,"RABV|CHPV|LBV")), aes(label=gsub("_", " ", `Virus name abbreviation(s)`)), 
                   size=2, min.segment.length = 0, fill = alpha(c("white"),0.8))+
  #geom_label_repel(data=. %>% filter(str_starts(`Virus name abbreviation(s)`,"RABV|VSV|CHPV|LBV")), aes(label=`Virus name abbreviation(s)`), 
  #                 size=2, min.segment.length = 0, fill = alpha(c("white"),0.8))+
  labs(title = "*Rhabdoviridae*")+
  theme(plot.title = element_markdown())
rhabdo_pca_plot2

#rda_dinuq %>% 
#  mutate(acc=str_remove(acc, "\\.[0-9]")) %>% 
#  left_join(rhabdo, by=join_by(acc==`Virus GENBANK accession`)) %>% 
#  filter(`Host source` %in% c("vertebrates", "invertebrates")) %>% 
#  ggpairs(aes(colour =`Host source`), columns = 2:17)

#' ### SDUc
rhabdo_sduc <- read_tsv("data/dinuq/rhabdo/dinuq_rhabdo_cds_sdu.tsv")
nt_cont <- read_tsv("data/dinuq/rhabdo/rhabdo_refseq_ntcont.tsv")

sdu_pca <- rhabdo_sduc %>% 
  column_to_rownames("acc") %>% 
  t() %>% 
  prcomp()

joined_sdu_df <- as.data.frame(sdu_pca$rotation) %>% 
  rownames_to_column("acc") %>% 
  mutate(acc=str_remove(acc, "\\.[0-9]")) %>% 
  full_join(rhabdo, by=join_by(acc==`Virus REFSEQ accession`))

joined_sdu_df <- as.data.frame(sdu_pca$rotation) %>% 
  rownames_to_column("acc") %>% 
  mutate(acc=str_remove(acc, "\\.[0-9].*"),
         acc=str_remove(acc, "join\\("),
         acc=str_remove(acc, "complement\\(")) %>% 
  full_join(rhabdo, by=join_by(acc==`Virus REFSEQ accession`)) %>% 
  rename(Host=`Host source`) %>% 
  full_join(lytras_rhabdo %>% mutate(acc=str_remove(acc, "\\.[0-9].*")), by="acc") %>% 
  mutate(Host.x = coalesce(Host.x, Host.y),
         Host.x=if_else(Host.x=="vertebrates", "vertebrate", Host.x)) %>% 
  filter(!is.na(PC1))

rhabdo_sdu_pca_plot <- joined_sdu_df %>% 
  filter(Host.x %in% c("vertebrate", "invertebrates", "invertebrates, vertebrates", NA)) %>% 
  ggplot(aes(x=PC1, y=PC2))+
  stat_density_2d(data=. %>% filter(!is.na(Host.x)), aes(fill=Host.x), 
                  color=NA, alpha=.15, geom="polygon", contour_var = "count")+
  geom_point(data=. %>% filter(!is.na(Host.x)), aes(color=Host.x), size=1)+
  scale_fill_viridis_d(name="Host", guide=guide_legend(order=1))+
  scale_color_viridis_d(name="Host", guide=guide_legend(order=1))+
  theme_bw()
rhabdo_sdu_pca_plot

rhabdo_sdu_pca_plot2 <- rhabdo_sdu_pca_plot +
  new_scale_color()+
  geom_point(data=. %>% filter(startsWith(acc, "OVRV1")), aes(x=PC1, y=PC2, color="#ff7f00"), size=3)+
  geom_point(data=. %>% filter(startsWith(acc, "blackfly")), aes(x=PC1, y=PC2, color="#377eb8"), size=3)+
  geom_point(data=. %>% filter(startsWith(acc, "NC_077192")), aes(x=PC1, y=PC2, color="#a65628"), size=3)+
  scale_color_manual(name = 'Virus of interest',
                     values =c("#ff7f00"="#ff7f00", "#377eb8"="#377eb8", "#a65628"="#a65628"), 
                     labels = c("#ff7f00"='OVRV1', "#377eb8"='Blackfly rhabdoviruses',  "#a65628"="Mundri virus"), 
                     guide=guide_legend(order=2))+
  geom_label_repel(data=. %>% filter(str_starts(acc,"OVRV1|blackfly|NC_077192")), aes(label=acc), 
                   size=3, min.segment.length = 0, fill = alpha(c("white"),0.8))+
  labs(title = "*Rhabdoviridae*")+
  theme(plot.title = element_markdown())
rhabdo_sdu_pca_plot2

#' ### Distance based

rhabdo_dist <- vegan::vegdist(rhabdo_sduc %>% column_to_rownames("acc"))
rhabdo_pcoa <- ape::pcoa(rhabdo_dist)

pcoa_df <- data.frame(rhabdo_pcoa$vectors) %>% 
  rownames_to_column("acc") %>% 
  mutate(acc=str_remove(acc, "\\.[0-9].*"),
         acc=str_remove(acc, "join\\("),
         acc=str_remove(acc, "complement\\(")) %>% 
  full_join(rhabdo, by=join_by(acc==`Virus REFSEQ accession`)) %>% 
  rename(Host=`Host source`) %>% 
  full_join(lytras_rhabdo %>% mutate(acc=str_remove(acc, "\\.[0-9].*")), by="acc") %>% 
  mutate(Host.x = coalesce(Host.x, Host.y),
         Host.x=if_else(Host.x=="vertebrates", "vertebrate", Host.x)) %>% 
  filter(!is.na(Axis.1))

rhabdo_sdu_pcoa_plot <- pcoa_df %>% 
  filter(Host.x %in% c("vertebrate", "invertebrates", "invertebrates, vertebrates", NA)) %>% 
  ggplot(aes(x=Axis.1, y=Axis.2))+
  stat_density_2d(data=. %>% filter(!is.na(Host.x)), aes(fill=Host.x), 
                  color=NA, alpha=.15, geom="polygon", contour_var = "count")+
  geom_point(data=. %>% filter(!is.na(Host.x)), aes(color=Host.x), size=1)+
  scale_fill_viridis_d(name="Host", guide=guide_legend(order=1))+
  scale_color_viridis_d(name="Host", guide=guide_legend(order=1))+
  theme_bw()
rhabdo_sdu_pcoa_plot

rhabdo_sdu_pcoa_plot2 <- rhabdo_sdu_pcoa_plot +
  new_scale_color()+
  geom_point(data=. %>% filter(startsWith(acc, "OVRV1")), aes(x=Axis.1, y=Axis.2, color="#ff7f00"), size=3)+
  geom_point(data=. %>% filter(startsWith(acc, "blackfly")), aes(x=Axis.1, y=Axis.2, color="#377eb8"), size=3)+
  geom_point(data=. %>% filter(startsWith(acc, "NC_077192")), aes(x=Axis.1, y=Axis.2, color="#a65628"), size=3)+
  scale_color_manual(name = 'Virus of interest',
                     values =c("#ff7f00"="#ff7f00", "#377eb8"="#377eb8", "#a65628"="#a65628"), 
                     labels = c("#ff7f00"='OVRV1', "#377eb8"='Blackfly rhabdoviruses',  "#a65628"="Mundri virus"), 
                     guide=guide_legend(order=2))+
  geom_label_repel(data=. %>% filter(str_starts(acc,"OVRV1|blackfly|NC_077192")), aes(label=acc), 
                   size=3, min.segment.length = 0, fill = alpha(c("white"),0.8))+
  labs(title = "*Rhabdoviridae*")+
  theme(plot.title = element_markdown())
rhabdo_sdu_pcoa_plot2

#rhabdo_sduc %>% 
#  left_join(nt_cont) %>% 
#  mutate(acc=str_remove(acc, "\\..*")) %>% 
#  full_join(rhabdo, by=join_by(acc==`Virus REFSEQ accession`)) %>% 
#  rowwise() %>% 
#  mutate(GC=sum(G, C)) %>% 
#  select(acc, starts_with("CpG"), GC, `Host source`) %>% 
#  filter(`Host source` %in% c("vertebrates", "invertebrates", NA)) %>% 
#  #group_by(acc) %>% 
#  #mutate(CpGpos1=mean(CpGpos1), CpGpos2=mean(CpGpos2), CpGbridge=mean(CpGbridge), GC=mean(GC)) %>% 
#  #distinct() %>% 
#  pivot_longer(cols = starts_with("CpG"), names_to = "pos", values_to = "SDUc") %>% 
#  ggplot(aes(x=GC, y=SDUc, color=`Host source`, shape=pos))+
#  geom_point()+
#  new_scale_color()+
#  geom_point(data=. %>% filter(startsWith(acc, "G")), aes(x=GC, y=SDUc), size=5, color="red")

#' ## Flaviviridae
flavi <- vmr %>% 
  filter(Genus == "Orthoflavivirus")

lytras_flavi <- read_excel("~/Downloads/viruses-746103-supplementary.xlsx", sheet = "Flaviviridae") %>% 
  select(accession, Host) %>% 
  distinct() %>% 
  rename(acc=accession)

#' ### RDA

flavi_rda <- read_tsv("data/dinuq/flavi/dinuq_flavi_rda.tsv")

rda_pca_flavi <- flavi_rda %>% 
  column_to_rownames("acc") %>% 
  t() %>% 
  prcomp()

joined_rda_df_flavi <- as.data.frame(rda_pca_flavi$rotation) %>% 
  rownames_to_column("acc") %>% 
  mutate(acc=str_remove(acc, "\\.[0-9].*")) %>% 
  full_join(flavi, by=join_by(acc==`Virus GENBANK accession`)) %>% 
  rename(Host=`Host source`) %>% 
  full_join(lytras_flavi %>% mutate(acc=str_remove(acc, "\\.[0-9].*")), by="acc") %>% 
  mutate(Host.x = coalesce(Host.x, Host.y),
         Host.x=if_else(Host.x=="vertebrates", "vertebrate", Host.x)) %>% 
  filter(!is.na(PC1))

flavi_rda_pca_plot <- joined_rda_df_flavi %>% 
  ggplot(aes(x=PC1, y=PC2))+
  stat_density_2d(data=. %>% filter(!is.na(Host.x)), aes(fill=Host.x), 
                  color=NA, alpha=.15, geom="polygon")+
  geom_point(data=. %>% filter(!is.na(Host.x)), aes(color=Host.x), size=1)+
  scale_color_viridis_d(name="Host",
                        guide=guide_legend(order=1))+
  scale_fill_viridis_d(name="Host",
                       guide=guide_legend(order=1))+
  new_scale_color()+
  geom_point(data=. %>% filter(startsWith(acc, "Blackfly")), aes(color="#D2352B"), size=3)+
  scale_color_manual(name = 'Virus of interest', 
                     values =c("#D2352B"="#D2352B"), 
                     labels = c("#D2352B"="Blackfly flavivirus"),
                     guide=guide_legend(order=2))+
  geom_label_repel(data=. %>% filter(str_starts(`Virus name abbreviation(s)`,"WNV|ZIKV|JEV|TBEV")), aes(label=`Virus name abbreviation(s)`), 
                   size=2, min.segment.length = 0, fill=alpha("white", 0.8))+
  scale_x_continuous(limit=c(-0.155, -0.055))+
  scale_y_continuous(limit=c(-0.4, 0.15))+
  labs(title = "*Flaviviridae*")+
  theme_bw()+
  theme(plot.title = element_markdown())
flavi_rda_pca_plot


#' ### SDUc
flavi_sduc <- read_tsv("data/dinuq/flavi/dinuq_flavi_cds_sdu.tsv")
nt_cont_flavi <- read_tsv("data/dinuq/flavi/flavi_ntcont.tsv")

sdu_pca_flavi <- flavi_sduc %>% 
  column_to_rownames("acc") %>% 
  t() %>% 
  prcomp()

joined_sdu_df_flavi <- as.data.frame(sdu_pca_flavi$rotation) %>% 
  rownames_to_column("acc") %>% 
  mutate(acc=str_remove(acc, "\\.[0-9].*")) %>% 
  full_join(flavi, by=join_by(acc==`Virus GENBANK accession`)) %>% 
  rename(Host=`Host source`) %>% 
  full_join(lytras_flavi %>% mutate(acc=str_remove(acc, "\\.[0-9].*")), by="acc") %>% 
  mutate(Host.x = coalesce(Host.x, Host.y),
         Host.x=if_else(Host.x=="vertebrates", "vertebrate", Host.x))

flavi_pca_plot <- joined_sdu_df_flavi %>% 
  ggplot(aes(x=PC1, y=PC2))+
  stat_density_2d(data=. %>% filter(!is.na(Host.x)), aes(fill=Host.x), 
                  color=NA, alpha=.15, geom="polygon")+
  geom_point(data=. %>% filter(!is.na(Host.x)), aes(color=Host.x), size=1)+
  scale_color_viridis_d(name="Host",
                        guide=guide_legend(order=1))+
  scale_fill_viridis_d(name="Host",
                       guide=guide_legend(order=1))+
  new_scale_color()+
  geom_point(data=. %>% filter(startsWith(acc, "blackfly")), aes(color="#D2352B"), size=3)+
  scale_color_manual(name = 'Virus of interest', 
                     values =c("#D2352B"="#D2352B"), 
                     labels = c("#D2352B"="Blackfly flavivirus"),
                     guide=guide_legend(order=2))+
  scale_x_continuous(limit=c(-0.14, -0.0395))+ 
  #scale_x_continuous(limit=c(-0.127, -0.065))+ # prcomp scale=T
  labs(title = "*Flaviviridae*")+
  theme_bw()
flavi_pca_plot
ggsave("figures/flavi_dinuq_pca.pdf", dpi=300)


#' ## Combine plots
p <- (flavi_rda_pca_plot + rhabdo_pca_plot2)&
  theme(legend.position = "bottom",
        plot.title = element_markdown())

p2 <- p+
  plot_layout(axis_titles = "collect", guides = "collect")+
  plot_annotation(tag_levels = "A")&
  theme(plot.tag = element_text(face="bold"))
#ggsave("figures/dinuq.pdf", dpi=300, width=7, height = 4.5)
ggsave("figures/tiff/figure6.pdf", dpi=300, width=7, height = 4.5)

# Create combined legend
addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize)),
             fill="none") +
    theme(legend.title = element_text(size = textSize+1), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box = "vertical", 
          legend.background = element_blank(),
          legend.spacing.y = unit(-0.1, "cm"))
}

leg <- rbind(joined_rda_df_flavi %>% rename("Virus accession"=`Virus REFSEQ accession`), 
      joined_rda_df %>% rename("Virus accession"=`Virus GENBANK accession`)) %>% 
  filter(Host.x %in% c("vertebrate", "invertebrates", "invertebrates, vertebrates", NA)) %>% 
  ggplot(aes(x=PC1, y=PC2))+
  stat_density_2d(data=. %>% filter(!is.na(Host.x)), aes(fill=Host.x), 
                  color=NA, alpha=.15, geom="polygon")+
  geom_point(data=. %>% filter(!is.na(Host.x)), aes(color=Host.x), size=1)+
  scale_color_viridis_d(name="Host",
                        guide=guide_legend(order=1))+
  scale_fill_viridis_d(name="Host",
                       guide=guide_legend(order=1))+
  new_scale_color()+
  geom_point(data=. %>% filter(startsWith(acc, "Blackfly")), aes(color="#D2352B"), size=3)+
  geom_point(data=. %>% filter(startsWith(acc, "OVRV1")), aes(x=PC1, y=PC2, color="#ff7f00"), size=3)+
  geom_point(data=. %>% filter(startsWith(acc, "Blackfly")), aes(x=PC1, y=PC2, color="#377eb8"), size=3)+
  geom_point(data=. %>% filter(startsWith(acc, "NC_077192")), aes(x=PC1, y=PC2, color="#a65628"), size=3)+
  scale_color_manual(name = 'Virus of interest',
                     values =c("#ff7f00"="#ff7f00", "#377eb8"="#377eb8", "#a65628"="#a65628", "#D2352B"="#D2352B"), 
                     labels = c("#D2352B"="Blackfly flavivirus", "#377eb8"='Blackfly rhabdoviruses', "#ff7f00"='OVRV1', "#a65628"="Mundri virus"), 
                     limits=c("#D2352B", "#377eb8", "#ff7f00", "#a65628"),
                     guide=guide_legend(order=2))+
  theme_bw()

leg2 <- get_legend(addSmallLegend(leg, pointSize = 2, textSize = 8))
p2 <- p+
  plot_layout(axis_titles = "collect", guides = "collect")+
  plot_annotation(tag_levels = "A")&
  theme(plot.tag = element_text(face="bold"),
        legend.position = "none", 
        plot.background = element_blank())

p3 <- as_ggplot(leg2)

# Plot combined + custom legend
p2/p3+plot_layout(heights = c(1,.1))+plot_annotation(tag_levels = list(c("A", "B", "")))
ggsave("figures/pdf/figure5.pdf", dpi=300, width=7, height = 4.5)
ggsave("figures/tiff/figure5.tiff", dpi=300, width=7, height = 4.5)
