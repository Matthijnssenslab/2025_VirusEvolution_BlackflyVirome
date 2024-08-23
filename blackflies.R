necessary_packages <- c('tidyverse', 'data.table', 'scales', 'knitr', 'grid', 'magrittr', 'DT', 'here', 'gridtext', 'gtools',
                        'ggpubr', 'ggthemes', 'reshape2', 'viridis', 'pals', 'circlize', 'readxl', 'VennDiagram', 'ggtext',
                        'vegan', 'phyloseq', 'metagenomeSeq', 'ComplexHeatmap', 'decontam', 'ggrepel', 'glue', 'patchwork')
lapply(necessary_packages, library, character.only = TRUE)
here::i_am("blackflies.R")

# Read statistics
raw_reads <- readr::read_tsv("data/blackflies_raw_reads.tsv") %>% 
  separate(col=file, into=c("Sample"), sep = "\\.") %>% 
  group_by(Sample) %>% 
  summarise(num_seqs=sum(num_seqs), sum_len=sum(sum_len), .groups = "drop")

trimmed_reads <- readr::read_tsv("data/blackflies_trimmed_reads.tsv") %>% 
  separate(col=file, into=c("Sample"), sep = "\\.") %>% 
  group_by(Sample) %>% 
  summarise(num_seqs=sum(num_seqs), sum_len=sum(sum_len), .groups = "drop")

reads <- full_join(raw_reads, trimmed_reads,
          by="Sample", suffix = c("_raw", "_trimmed"))

write_tsv(reads %>% arrange(order(mixedorder(Sample))), "tables/supplementary_table_1.tsv")
xlsx::write.xlsx(x=as.data.frame(reads %>% arrange(order(mixedorder(Sample)))), 
                 file="tables/supplementary_table_1.xlsx", 
                 row.names = F)

summary(reads %>% select(-Sample))

reads %>% 
  summarise(total_raw=sum(num_seqs_raw), total_trimmed=sum(num_seqs_trimmed), loss=100-total_trimmed/total_raw*100)

# Eukaryotic virome analysis
# Prepare the data
## Load the OTU table, taxonomy file and metadata into R
OTU <- read.table("data/abundance_table", header=TRUE, row.names=1, sep="\t", dec=".")

names(OTU) <- gsub(x = names(OTU), pattern = "\\.", replacement = "-")
summary(rowSums(OTU))
summary(colSums(OTU))

# Taxonomy information
tax <- read.table("data/blackflies_viruses_1000bp_taxfile.tsv", header=TRUE, row.names=1, sep="\t", dec=".")
genomad <- read.table("data/blackflies.clustered_virus_summary.tsv", header=TRUE, row.names=1, sep="\t", dec=".")
checkv <- read.table("data/quality_summary.tsv", header=TRUE, row.names=1, sep="\t", dec=".")
blastx <- read.table("data/blackflies_virus_blastx.tab", header=TRUE, row.names=1, sep="\t", dec=".") %>% 
  select(-taxID)

genomad_1000 <- genomad %>% 
  filter(length >= 1000) 

genomad_taxonomy <- genomad_1000 %>% 
  select(taxonomy) %>% 
  separate(taxonomy, 
           into = c("Virus","Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = ";") %>% 
  mutate(
    Family = ifelse(grepl("viridae$", Order), Order, Family),
    Order = ifelse(grepl("viridae$", Order), NA, Order)
  ) %>% 
  mutate(
    Family = ifelse(grepl("viridae$", Realm), Realm, Family),
    Realm = ifelse(grepl("viridae$", Realm), NA, Realm)
  ) %>% 
  mutate(
    Family = ifelse(grepl("viridae$", Kingdom), Kingdom, Family),
    Kingdom = ifelse(grepl("viridae$", Kingdom), NA, Kingdom)
  )

# Select only scaffolds > 1000 bases
tax_1000 <- tax %>%
  rownames_to_column(var = "RowName") %>%
  mutate(length = as.double(str_split(RowName, "_", simplify = TRUE)[, 4])) %>% 
  filter(length >= 1000) %>% 
  select(-length) %>% 
  column_to_rownames(var="RowName") %>% 
  rename(Virus=Kingdom)

genomad_taxonomy <- genomad_taxonomy[order(rownames(genomad_taxonomy)), ]
tax_1000 <- tax_1000[order(rownames(tax_1000)), ]
checkv <- checkv[order(rownames(checkv)), ]

merged_taxonomy <- left_join(genomad_taxonomy %>% rownames_to_column(), 
                             tax_1000 %>% rownames_to_column(), 
                             by=c("rowname"), suffix = c(".genomad", ".blastx"))

merged_taxonomy_full <- full_join(genomad_taxonomy %>% rownames_to_column(), 
                                  tax_1000 %>% rownames_to_column(), 
                                  by=c("rowname"), suffix = c(".genomad", ".blastx"))

merged_tax_checkv <- left_join(merged_taxonomy_full,
                               checkv %>% rownames_to_column(),
                               by="rowname")

merged_tax_checkv_blastx <- left_join(merged_tax_checkv,
                                      blastx %>% rownames_to_column(),
                                      by="rowname")

# Metadata
meta <- read.table("data/blackfly_meta.csv", header=TRUE, row.names = 1, sep=";", dec=".")
meta <- cbind(Sample=rownames(meta), meta)

## Make a phyloseq object
OTU.UF <- otu_table(as.matrix(OTU), taxa_are_rows=T)
tax.UF <- tax_table(as.matrix(merged_tax_checkv_blastx %>% column_to_rownames(var = "rowname")))
meta.UF <- sample_data(meta)

BF <- phyloseq(OTU.UF, tax.UF, meta.UF)
BF

## Remove contamination
### Visualize library sizes of samples and negative controls
decontam <- as.data.frame(sample_data(BF))
decontam
decontam$LibrarySize <- sample_sums(BF)
decontam <- decontam[order(decontam$LibrarySize),]
decontam$Index <- seq(nrow(decontam))
ggplot(data=decontam, aes(x=Index, y=LibrarySize, color=Control)) + 
  geom_point()+
  ggtitle("Library sizes")

ggplot(data=decontam, aes(y="libraries", x=LibrarySize, color=Control)) + 
  geom_jitter(height = .01)+
  scale_x_log10()+
  ggtitle("Library sizes")

### Detect contaminants
sample_data(BF)$is.neg <- sample_data(BF)$Control == "Yes"
sample_data(BF)$is.neg
contamdf.prev <- decontam::isContaminant(BF, method="prevalence", neg="is.neg")


# **Number of contaminating contigs:**
table(contamdf.prev$contaminant)

contaminants <- rownames(contamdf.prev[contamdf.prev$contaminant==T,])

filtered_tax <- merged_tax_checkv_blastx |> 
  # Filter out all contaminants
  filter(!rowname %in% contaminants,
          # Filter out all sequences with lower than 50% completeness prediction, except if the Realm is NA (blastx prediciton) or they belong to the Riboviria
          completeness >= 50 | is.na(Realm) | Realm == "Riboviria",
          # Filter out all BlastX that are only "bacteriophage" species and have a completeness < 50%
          !(str_detect(Species.blastx, "(?i)bacteriophage") & completeness < 50 & Realm!="Riboviria"),
          # Filter out Blastx caudoviricetes not detected by genomad and lower than 50% completeness prediction
          !(Class.blastx=="Caudoviricetes" & is.na(Realm) & completeness < 50)
         ) |> 
  mutate(across(Family.blastx:Species.blastx, ~ifelse(perc_identity < 70, NA, .)))

#write_delim(as.data.frame(merged_tax_checkv_blastx[is.na(merged_tax_checkv_blastx$Virus.x),]$rowname), "blastx_virus.txt",
#            col_names = F)

ICTV_classification <- readxl::read_excel("data/VMR_MSL38_v2.xlsx") %>% 
  select(Realm, Kingdom, Phylum, Class, Order, Family, Genus) %>% 
  distinct() %>% 
  filter(rowSums(!is.na(.)) > 0)

columns_to_fill <- c("Phylum", "Class", "Order", "Family", "Genus")

filtered_tax2 <- filtered_tax

# Complete taxonomy based on ICTV
filtered_tax2 <- reduce(columns_to_fill, function(df, col_name) {
  join_condition <- paste0(col_name, ".blastx")
  
  df %>%
    left_join(
      ICTV_classification %>% select(Realm, Kingdom, all_of(col_name)) %>% distinct(),
      by = join_by({{ join_condition }} == {{ col_name }}),
      na_matches = "never"
    ) %>%
    mutate(
      Realm.x = coalesce(Realm.x, Realm.y),
      Kingdom.x = coalesce(Kingdom.x, Kingdom.y)
    ) %>%
    rename(Realm = Realm.x, Kingdom = Kingdom.x) %>%
    select(-Realm.y, -Kingdom.y)
}, .init = filtered_tax2)


filtered_tax2 <- filtered_tax2 %>% 
  filter(is.na(Realm) | Realm == "Riboviria" | completeness >= 50)

filtered_tax3 <- filtered_tax2 %>% 
  # If genomad had no prediction, use blastx taxonomy
  mutate(Phylum.genomad = coalesce(Phylum.genomad, Phylum.blastx),
         Class.genomad = ifelse(Phylum.genomad == Phylum.blastx,
                          coalesce(Class.genomad, Class.blastx),
                          Class.genomad),
         Order.genomad = ifelse(Class.genomad == Class.blastx,
                          coalesce(Order.genomad, Order.blastx),
                          Order.genomad),
         Family.genomad = ifelse(Order.genomad == Order.blastx,
                           coalesce(Family.genomad, Family.blastx),
                           Family.genomad),
         Genus.genomad = ifelse(Family.genomad == Family.blastx,
                          coalesce(Genus.genomad, Genus.blastx),
                          Genus.genomad),
         Species.genomad = ifelse(Genus.genomad == Genus.blastx,
                            coalesce(Species.genomad, Species.blastx),
                            Species.genomad),
         ) %>% 
  mutate(across(Phylum.genomad:Species.genomad, ~ifelse(. == "unclassified", NA, .)))

taxonomy_df <- filtered_tax3 %>% 
  select(rowname, Realm,  Kingdom, Phylum.genomad, Class.genomad, Order.genomad, Family.genomad, Genus.genomad, Species.genomad) %>% 
  rename_with(~sub("\\.genomad$", "", .), ends_with(".genomad"))

## Pie charts
# Realm pie chart
realms <- unique(taxonomy_df$Realm)

realm_col <- c(rev(RColorBrewer::brewer.pal(4, "Paired"))[1:length(realms)-1], "lightgrey")
names(realm_col) <- c(realms[order(realms)][1:3], "Unclassified")

realm_df <- taxonomy_df %>% 
  mutate_at(vars(Realm), ~replace_na(., "Unclassified")) %>% 
  group_by(Realm) %>% 
  count() %>% 
  ungroup() %>% 
  arrange(desc(Realm)) %>% 
  mutate(text_y = cumsum(n) - n/2, n_label=comma(n),
         perc=n/sum(n))

realm_plot <- realm_df %>% 
  ggplot(aes(x = 1, y = perc, fill = Realm)) +
  geom_col(show.legend = T) +
  coord_polar(theta = "y")+
  xlim(c(-1.5, 1.5)) +
  geom_label(aes(label = n_label), show.legend = F, 
             position = position_stack(vjust = .5)) +
  scale_fill_manual(values = realm_col)+
  #labs(title = "All virus contigs")+
  theme_void()+
  theme(legend.text = element_text(face="italic"))+
  annotate(geom = "text",
           label = "All virus contigs",
           #size = 6,
           x = -1.5,
           y = 0)
realm_plot

realm_df %>% 
  ggplot(aes(x=Realm, y=n, fill=Realm))+
  geom_col(show.legend = F)+
  scale_fill_manual(values = realm_col)+
  labs(y="# Contigs")+
  theme_bw()+
  theme(axis.title.x = element_blank())

realm_read_df <- taxonomy_df %>% 
  left_join(as.data.frame(rowSums(OTU)) %>% rownames_to_column()) %>% 
  rename(Abundance=`rowSums(OTU)`) %>%
  mutate_at(vars(Realm), ~replace_na(., "Unclassified")) %>% 
  group_by(Realm) %>% 
  summarise(n=sum(Abundance), .groups = "drop") %>% 
  arrange(desc(Realm)) %>% 
  mutate(text_y = cumsum(n) - n/2, n_label=comma(n),
         perc=n/sum(n))

realm_read_plot <- realm_read_df %>% 
  ggplot(aes(x = 1, y = perc, fill = Realm)) +
  geom_col(show.legend = T) +
  coord_polar(theta = "y")+
  xlim(c(-1.5, 1.5)) +
  geom_label(aes(label = n_label), show.legend = F, 
             position = position_stack(vjust = .5)) +
  scale_fill_manual(values = realm_col)+
  theme_void()+
  theme(legend.text = element_text(face="italic"))+
  annotate(geom = "text",
           label = "All virus reads",
           #size = 6,
           x = -1.5,
           y = 0)

realm_read_df %>% 
  ggplot(aes(x=Realm, y=log10(n), fill=Realm))+
  geom_col(show.legend = F)+
  geom_label(aes(label = n_label), show.legend = F) +
  scale_fill_manual(values = realm_col)+
  labs(y="# Reads")+
  theme_bw()+
  theme(axis.title.x = element_blank())

# *Riboviria* pie chart
phyla <- unique(taxonomy_df[taxonomy_df$Realm=="Riboviria",]$Phylum)

#phyla_col <- c(brewer.pal(length(phyla)-1, "Paired"), "lightgrey")
phyla_col <- c(brewer.pal(10, "Paired")[5:10], "lightgrey")
names(phyla_col) <- c(phyla[order(phyla)][-length(phyla)], "Unclassified")

riboviria <- taxonomy_df %>% 
  filter(Realm=="Riboviria") 

riboviria_clean <- riboviria %>% 
  mutate_at(vars(Phylum), ~replace_na(., "Unclassified")) %>% 
  group_by(Phylum) %>% 
  count() %>% 
  ungroup() %>% 
  arrange(desc(Phylum)) %>% 
  mutate(text_y = cumsum(n) - n/2, n_label=comma(n),
         perc=n/sum(n))

ribophylum_plot <- riboviria_clean %>% 
  ggplot(aes(x = 1, y = perc, fill =  Phylum)) +
  geom_col(show.legend = T) +
  coord_polar(theta = "y")+
  xlim(c(-1.5, 1.5)) +
  ggrepel::geom_label_repel(aes(label = n_label), show.legend = F, 
             position = position_stack(vjust = .5), ) +
  scale_fill_manual(values = phyla_col)+
  theme_void()+
  theme(legend.text = element_text(face="italic"))+
  annotate(geom = "text",
           label = expression(paste(italic("Riboviria"), " contigs")),
           #size = 6,
           x = -1.5,
           y = 0)


riboviria_read_abundance <- riboviria %>% 
  left_join(as.data.frame(rowSums(OTU)) %>% rownames_to_column()) %>% 
  rename(Abundance=`rowSums(OTU)`) %>% 
  mutate_at(vars(Phylum), ~replace_na(., "Unclassified")) %>% 
  group_by(Phylum) %>% 
  summarise(n=sum(Abundance), .groups = "drop") %>% 
  arrange(desc(Phylum)) %>% 
  mutate(text_y = cumsum(n) - n/2, n_label=comma(n),
         perc=n/sum(n))

ribophylum_read_plot <- riboviria_read_abundance %>% 
  ggplot(aes(x = 1, y = perc, fill =  Phylum)) +
  geom_col(show.legend = T) +
  coord_polar(theta = "y")+
  xlim(c(-1.5, 1.5)) +
  ggrepel::geom_label_repel(aes(label = n_label), show.legend = F, 
             position = position_stack(vjust = .5)) +
  scale_fill_manual(values = phyla_col)+
  theme_void()+
  theme(legend.text = element_text(face="italic"))+
  annotate(geom = "text",
           label = expression(paste(italic("Riboviria"), " reads")),
           #size = 6,
           x = -1.5,
           y = 0)
ribophylum_read_plot

pie1 <- ggpubr::ggarrange(realm_plot, realm_read_plot,  
          labels = c('A', 'B'), common.legend = T,
          legend = "right", align="hv",
          font.label = list(size=20))

pie2 <- ggpubr::ggarrange(ribophylum_plot, ribophylum_read_plot, 
                  labels = c('C', 'D'), common.legend = T,
                  legend = "right", align="hv",
                  font.label = list(size=20))

ggpubr::ggarrange(pie1, pie2, ncol=1, align = 'hv')
ggsave("figures/piecharts.pdf", dpi=300)

# Barplots
abundance_df <- rbind(realm_df %>% mutate(level="Contigs", class="All viruses") %>% rename(tax=Realm),
                      realm_read_df %>% mutate(level="Reads", class="All viruses") %>% rename(tax=Realm),
                      riboviria_clean %>% mutate(level="Contigs", class="*Riboviria*") %>% rename(tax=Phylum),
                      riboviria_read_abundance %>% mutate(level="Reads", class="*Riboviria*") %>% rename(tax=Phylum) 
                        )

realm_p <- abundance_df %>% 
  filter(class=="All viruses") %>% 
  ggplot(aes(x=level, y=perc, fill=tax))+
  geom_col()+
  geom_label(aes(label = n_label), show.legend = F, 
             position = position_stack(vjust=0.5), size=3) +
  scale_fill_manual(values = c(realm_col, phyla_col), name="Realm")+
  facet_wrap(~class)+
  labs(y="Relative abundance")+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        strip.text = ggtext::element_markdown())

riboviria_p <- abundance_df %>% 
  filter(class=="*Riboviria*") %>% 
  ggplot(aes(x=level, y=perc, fill=tax))+
  geom_col()+
  geom_label(aes(label = n_label), show.legend = F, 
             position = position_stack(vjust=0.5), size=3) +
  scale_fill_manual(values = c(realm_col, phyla_col), name="Phylum")+
  facet_wrap(~class)+
  labs(y="Relative abundance")+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        strip.text = element_markdown())

realm_p + riboviria_p +
  plot_annotation(tag_level="A") &
  theme(plot.tag = element_text(face="bold"))

# Rarefaction curves

virome_df <- left_join(taxonomy_df, 
                       OTU %>% rownames_to_column()) %>% 
  select(!contains("NC"))

riboviria_df <- virome_df %>% 
  filter(Realm == "Riboviria") %>% 
  mutate(sample = sub('.*_', '', rowname))

duplodnaviria_df <- virome_df %>% 
  filter(Realm == "Duplodnaviria") %>% 
  mutate(sample = sub('.*_', '', rowname))

monodnaviria_df <- virome_df %>% 
  filter(Realm == "Monodnaviria") %>% 
  mutate(sample = sub('.*_', '', rowname))

rest_df <- virome_df %>% 
  filter(is.na(Realm)) %>% 
  mutate(sample = sub('.*_', '', rowname))

# Function to merge rows and sum across columns starting with "Blackfly"
merge_and_sum <- function(data) {
  data %>%
    group_by(sample, Realm, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    summarise(rowname = first(rowname),
              across(starts_with("BlackFly"), sum), 
              .groups = 'drop')
}

#riboviria_df_rare <- merge_and_sum(riboviria_df)
rest_df_rare <- merge_and_sum(rest_df)

# !! Make sure to run the RdRP tree script
riboviria_df_rare <- virome_df %>% 
  filter(rowname %in% 
           sub("_[^_]+$", "", tips_data[tips_data$contains_blackfly==T,]$tip)) %>% 
  mutate(sample = sub('.*_', '', rowname))

set.seed(123)

collect <- function(data) {
  data %>% 
    select(rowname, contains("BlackFly")) %>% 
    pivot_longer(cols=starts_with("BlackFly"), names_to = "Sample", values_to = "Reads") %>% 
    filter(Reads != 0) %>% 
    select(rowname, Sample) %>% 
    sample_n(n()) %>% 
    mutate(observation = row_number()) %>% 
    group_by(Sample) %>%
    mutate(min_observation=min(observation)) %>% 
    arrange(min_observation) %>%   
    ungroup() %>%
    group_by(rowname) %>% 
    mutate(distinct = row_number() == 1) %>% 
    ungroup() %>%
    group_by(Sample) %>% 
    mutate(s = sum(distinct)) %>% 
    select(Sample, s) %>% 
    distinct() %>% 
    ungroup() %>% 
    mutate(observation=row_number(),
           species=cumsum(s)) %>% 
    select(observation, species)
}

create_curve_df <- function(data, iterations=1000){
  collect_curves <- map_dfr(1:iterations, ~collect(data), .id="iteration")
  
  rarefaction_curve <- collect_curves %>%
    group_by(observation) %>%
    summarize(r = mean(species),
              stdev=sd(species),
              stderr=stdev/sqrt(length(collect_curves)),
              max_s=max(species),
              min_s=min(species))
  return(rarefaction_curve)
}

all_df_rare <- rbind(riboviria_df_rare, rest_df_rare, monodnaviria_df, duplodnaviria_df)


ribo_rarefaction_curve <- create_curve_df(riboviria_df_rare)%>% 
  mutate(Group="*Riboviria*")
duplo_rarefaction_curve <- create_curve_df(duplodnaviria_df)%>% 
  mutate(Group="*Duplodnaviria*")
mono_rarefaction_curve <- create_curve_df(monodnaviria_df) %>% 
  mutate(Group="*Monodnaviria*")
rest_rarefaction_curve <- create_curve_df(rest_df_rare)%>% 
  mutate(Group="Unclassified")
all_rarefaction_curve <- create_curve_df(all_df_rare)%>% 
  mutate(Group="Total")

curve_df <- rbind(ribo_rarefaction_curve,
                  duplo_rarefaction_curve,
                  mono_rarefaction_curve,
                  rest_rarefaction_curve,
                  all_rarefaction_curve)

curve_df$Group <- factor(curve_df$Group, levels = c("Total",
                                                    "*Duplodnaviria*",
                                                    "*Monodnaviria*",
                                                    "*Riboviria*",
                                                    "Unclassified"))

# Get unique values in the "Group" column
unique_groups <- unique(curve_df$Group)

# Create a data frame with all values set to 0
new_rows_df <- data.frame(
  Group = rep(unique_groups, each = 1)
)

curve_df <- bind_rows(curve_df, new_rows_df) %>% 
  replace(is.na(.), 0)

curve_col <- c(realm_col, "Total"="black")
names(curve_col) <- sub("^(.*viria)$", "*\\1*", names(curve_col))

rarefaction_plot <- curve_df %>%
  group_by(Group) %>%
  mutate(max_observation = max(observation),
         max_r = max(r)) %>%
  ggplot(aes(x=observation, y=r, fill=Group, color=Group, group=Group)) +
  geom_ribbon(aes(ymin = r - (stderr*qnorm(0.975)), ymax = r + (stderr*qnorm(0.975))), 
              alpha = 0.2, color=NA, show.legend = F) +
  geom_ribbon(aes(ymin = min_s, ymax = max_s), alpha = 0.2, color=NA, show.legend = F) +
  geom_line(size = 1.2, alpha = 0.95, show.legend = F)+ 
  ggtext::geom_richtext(data=. %>% select(Group, max_observation, max_r) %>% distinct(),
                aes(label=Group, x=max_observation, y=max_r), fill=NA, label.color=NA,
                show.legend = F, hjust=1, nudge_y = 25)+
  scale_color_manual(values=curve_col)+
  scale_fill_manual(values=curve_col)+
  scale_x_continuous(breaks = seq(min(curve_df$observation), max(curve_df$observation), by = 10)) +
  scale_y_continuous(breaks = seq(min(curve_df$r), max(curve_df$r)+100, by = 100)) +
  theme_bw()+
  labs(y="# vOTUs", x="Samples")+
  theme(legend.text = element_markdown(),
        legend.title = element_blank(),
        axis.ticks = )

rarefaction_plot

fig1 <- ((realm_plot | realm_read_plot) + plot_layout(guides = "collect")) /
  ((ribophylum_plot | ribophylum_read_plot) + plot_layout(guides = "collect")) /
  rarefaction_plot

fig1 & patchwork::plot_annotation(tag_levels = "A")&
  theme(plot.tag = element_text(face = 'bold'))

ggsave("figures/figure1.pdf", dpi=300, height = 8.2, width=7)

fig2 <- ((realm_p + riboviria_p)+ plot_layout(guides = "collect"))/
  rarefaction_plot+
  plot_layout(heights = c(1, .5))

fig2 & 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = 'bold'))

#ggsave("figures/figure_2_bis.pdf", dpi=300, height = 7.5, width=7)
ggsave("figures/tiff/figure3.tiff", dpi=300, height = 7.5, width=7)

## Venn diagram
genomad <- filtered_tax3[!is.na(filtered_tax3$Virus.genomad),]$rowname
blastx <- filtered_tax3[!is.na(filtered_tax3$Virus.blastx),]$rowname
hmm <- scan("data/HMM_RdRP_matches.txt", character(), quote = "")
#palm_annot <- scan("data/palm_annot_matches.txt", character(), quote = "")

single_blastx <- setdiff(setdiff(blastx, genomad), hmm)
single_hmm <- setdiff(setdiff(hmm, genomad), blastx)
single_genomad <- setdiff(setdiff(genomad, hmm), blastx)
#single_blastx <- setdiff(setdiff(blastx, genomad), palm_annot)
#single_palm_annot <- setdiff(setdiff(palm_annot, genomad), blastx)
#single_genomad <- setdiff(setdiff(genomad, palm_annot), blastx)

single_blastx_df  <- taxonomy_df %>% 
  filter(rowname %in% single_blastx)

single_blastx_df %>% 
  filter(str_detect(Species, "(?i)bacteriophage"))

myCol <- RColorBrewer::brewer.pal(3, "Dark2")


venn_result <- VennDiagram::venn.diagram(
  x = list(genomad, blastx, hmm),
  category.names = c("genomad" , "DIAMOND blastx", "NeoRdRP hmmsearch"),
  #x = list(genomad, blastx, palm_annot),
  #category.names = c("genomad" , "DIAMOND blastx", "palm_annot"),
  filename = NULL, #"figures/venn_diagram.tiff",
  #imagetype = "tiff",
  disable.logging=T,
  height = 480,
  width = 480,
  resolution = 300,
  #compression = "lzw",
  lwd = 1,
  col=myCol,
  fill = alpha(myCol, 0.3),
  # Numbers
  cex = .8,
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.8,
  cat.default.pos = "outer",
  cat.pos = c(-170, 153, 0),
  cat.dist = c(0.045, 0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = myCol,
  rotation = 1
)

grid.newpage()
grid.draw(venn_result)

# CheckV analysis

checkv_col <- rev(brewer.pal(6, "Blues"))[2:6]
names(checkv_col) <- c("Complete", "High-quality", "Medium-quality", 
                       "Low-quality", "Not-determined")

checkv_pie <- filtered_tax3 %>% 
  mutate(checkv_quality=factor(checkv_quality, levels=c("Complete", "High-quality", "Medium-quality", 
                                                        "Low-quality", "Not-determined"))) %>% 
  group_by(checkv_quality) %>% 
  count() %>% 
  ungroup() %>%
  mutate(text_y = sum(n) - cumsum(n) + n /2, 
         n_label=comma(n)) %>%
  ggplot(aes(x = "", y = n, fill = checkv_quality)) +
  geom_bar(stat="identity", width = 1) +
  coord_polar(theta = "y")+
  ggrepel::geom_label_repel(aes(label = n_label, y = text_y), nudge_x = 0.57, #nudge_y=0.5,
                   min.segment.length = 100, show.legend = F) +
  scale_fill_manual(values = checkv_col,  limits=names(checkv_col))+
  labs(fill="CheckV quality")+
  theme_void()+
  theme(legend.position = "right",
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=7))

ggpubr::ggarrange(venn_result, checkv_pie, labels = "AUTO", nrow=1)
#ggsave("figures/venn_checkv.pdf", dpi=300, width=160, height=80, units = "mm")
ggsave("figures/tiff/figure2.tiff", dpi=300, width=160, height=80, units = "mm")


# Correlation analysis
virus_cor <- read_tsv("data/virus_correlation.tsv")

ribo_virus_cor <- virus_cor %>% 
  filter(is.element(Contig1, riboviria$rowname) | is.element(Contig2, riboviria$rowname))

cor_1 <- ribo_virus_cor %>% 
  filter(Correlation > .7)

cor_1_tax <- left_join(cor_1,
          riboviria,
          by=join_by("Contig1" == "rowname")) %>% 
  left_join(riboviria, 
            by=join_by("Contig2" == "rowname"),
            suffix=c("_c1", "_c2")) %>% 
  select(Contig1:Species.x_c1, Virus.x_c2:Species.x_c2)
