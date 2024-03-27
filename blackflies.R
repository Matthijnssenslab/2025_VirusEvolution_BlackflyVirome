necessary_packages <- c('tidyverse', 'data.table', 'scales', 'knitr', 'grid', 'magrittr', 'DT', 'here', 'gridtext', 'gtools',
                        'ggpubr', 'ggthemes', 'reshape2', 'viridis', 'pals', 'circlize', 'readxl', 'VennDiagram',
                        'vegan', 'phyloseq', 'metagenomeSeq', 'ComplexHeatmap', 'decontam', 'ggrepel', 'glue', 'patchwork')
lapply(necessary_packages, library, character.only = TRUE)
i_am("blackflies.R")

#' ***
#' ## Read statistics
raw_reads <- read_tsv("data/blackflies_raw_reads.tsv") %>% 
  separate(col=file, into=c("Sample"), sep = "\\.") %>% 
  group_by(Sample) %>% 
  summarise(num_seqs=sum(num_seqs), sum_len=sum(sum_len), .groups = "drop")

trimmed_reads <- read_tsv("data/blackflies_trimmed_reads.tsv") %>% 
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

#' ***
#' # Eukaryotic virome analysis
#' ## Prepare the data
#' ### Load the OTU table, taxonomy file and metadata into R
OTU <- read.table("data/abundance_table", header=TRUE, row.names=1, sep="\t", dec=".")

names(OTU) <- gsub(x = names(OTU), pattern = "\\.", replacement = "-")
summary(rowSums(OTU))
summary(colSums(OTU))

#' ## Taxonomy information
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

#' Select only scaffolds > 1000 bases
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
                             by=c("rowname"))

merged_taxonomy_full <- full_join(genomad_taxonomy %>% rownames_to_column(), 
                                  tax_1000 %>% rownames_to_column(), 
                                  by=c("rowname"))

merged_tax_checkv <- left_join(merged_taxonomy_full,
                               checkv %>% rownames_to_column(),
                               by="rowname")

merged_tax_checkv_blastx <- left_join(merged_tax_checkv,
                                      blastx %>% rownames_to_column(),
                                      by="rowname")

#' Metadata
meta <- read.table("data/blackfly_meta.csv", header=TRUE, row.names = 1, sep=";", dec=".")
meta <- cbind(Sample=rownames(meta), meta)

#' ### Make a phyloseq object
OTU.UF <- otu_table(as.matrix(OTU), taxa_are_rows=T)
tax.UF <- tax_table(as.matrix(merged_tax_checkv_blastx %>% column_to_rownames(var = "rowname")))
meta.UF <- sample_data(meta)

BF <- phyloseq(OTU.UF, tax.UF, meta.UF)
BF

#' ### Remove contamination
#' #### Visualize library sizes of samples and negative controls
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

#' #### Detect contaminants
sample_data(BF)$is.neg <- sample_data(BF)$Control == "Yes"
sample_data(BF)$is.neg
contamdf.prev <- isContaminant(BF, method="prevalence", neg="is.neg")


#' **Number of contaminating contigs:**
table(contamdf.prev$contaminant)

contaminants <- rownames(contamdf.prev[contamdf.prev$contaminant==T,])

filtered_tax <- merged_tax_checkv_blastx %>% 
  # Filter out all contaminants
  filter(!rowname %in% contaminants,
          # Filter out all sequences with lower than 50% completeness prediction, except if the Realm is NA (blastx prediciton) or they belong to the Riboviria
          completeness >= 50 | is.na(Realm) | Realm == "Riboviria",
          # Filter out all BlastX that are only "bacteriophage" species and have a completeness < 50%
          !(str_detect(Species.y, "(?i)bacteriophage") & completeness < 50 & Realm!="Riboviria"),
          # Filter out Blastx caudoviricetes not detected by genomad and lower than 50% completeness prediction
          !(Class.y=="Caudoviricetes" & is.na(Realm) & completeness < 50)
         )

#write_delim(as.data.frame(merged_tax_checkv_blastx[is.na(merged_tax_checkv_blastx$Virus.x),]$rowname), "blastx_virus.txt",
#            col_names = F)

ICTV_classification <- read_excel("data/VMR_MSL38_v2.xlsx") %>% 
  select(Realm, Kingdom, Phylum, Class, Order, Family, Genus) %>% 
  distinct() %>% 
  filter(rowSums(!is.na(.)) > 0)

#filtered_tax2 <- filtered_tax %>% 
#  mutate(across(Phylum.y:Species.y, ~ifelse(. == "unclassified", NA, .))) %>% 
#  left_join(ICTV_classification %>% select(Realm, Kingdom, Phylum) %>% distinct(),
#            by=join_by(Phylum.y == Phylum),
#                     #Class.y == Class,
#                     #Order.y == Order),
#                     #Family.y == Family,
#                     #Genus.y == Genus),
#            na_matches = "never") %>% 
#  mutate(Realm.x = coalesce(Realm.x, Realm.y),
#         Kingdom.x = coalesce(Kingdom.x, Kingdom.y)) %>% 
#  rename(Realm = Realm.x, Kingdom = Kingdom.x) %>% 
#  select(-Realm.y, -Kingdom.y)
#
#filtered_tax2 <- filtered_tax2 %>% 
#  left_join(ICTV_classification %>% select(Realm, Kingdom, Class) %>% distinct(),
#            by=join_by(Class.y == Class),
#            na_matches = "never") %>% 
#  mutate(Realm.x = coalesce(Realm.x, Realm.y),
#         Kingdom.x = coalesce(Kingdom.x, Kingdom.y)) %>% 
#  rename(Realm = Realm.x, Kingdom = Kingdom.x) %>% 
#  select(-Realm.y, -Kingdom.y)
#
#filtered_tax2 <- filtered_tax2 %>% 
#  left_join(ICTV_classification %>% select(Realm, Kingdom, Order) %>% distinct(),
#            by=join_by(Order.y == Order),
#            na_matches = "never") %>% 
#  mutate(Realm.x = coalesce(Realm.x, Realm.y),
#         Kingdom.x = coalesce(Kingdom.x, Kingdom.y)) %>% 
#  rename(Realm = Realm.x, Kingdom = Kingdom.x) %>% 
#  select(-Realm.y, -Kingdom.y)
#
#filtered_tax2 <- filtered_tax2 %>% 
#  left_join(ICTV_classification %>% select(Realm, Kingdom, Family) %>% distinct(),
#            by=join_by(Family.y == Family),
#            na_matches = "never") %>% 
#  mutate(Realm.x = coalesce(Realm.x, Realm.y),
#         Kingdom.x = coalesce(Kingdom.x, Kingdom.y)) %>% 
#  rename(Realm = Realm.x, Kingdom = Kingdom.x) %>% 
#  select(-Realm.y, -Kingdom.y)
#
#filtered_tax2 <- filtered_tax2 %>% 
#  left_join(ICTV_classification %>% select(Realm, Kingdom, Genus) %>% distinct(),
#            by=join_by(Genus.y == Genus),
#            na_matches = "never") %>% 
#  mutate(Realm.x = coalesce(Realm.x, Realm.y),
#         Kingdom.x = coalesce(Kingdom.x, Kingdom.y)) %>% 
#  rename(Realm = Realm.x, Kingdom = Kingdom.x) %>% 
#  select(-Realm.y, -Kingdom.y)
#
#filtered_tax2 <- filtered_tax2 %>% 
#  filter(completeness >= 50 | is.na(Realm) | Realm == "Riboviria")


columns_to_fill <- c("Phylum", "Class", "Order", "Family", "Genus")

filtered_tax2 <- filtered_tax


filtered_tax2 <- reduce(columns_to_fill, function(df, col_name) {
  join_condition <- paste0(col_name, ".y")
  
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


#ribo_phyla <- c("Duplornaviricota", 
#                "Kitrinoviricota", 
#                "Lenarviricota", 
#                "Negarnaviricota",
#                "Pisuviricota", 
#                "Artverviricota")
#
#ribo_fam <- c("Birnaviridae", 
#              "Permutotetraviridae", 
#              "Polymycoviridae", 
#              "Sarthroviridae")
#
#ribo_gen <- c("Botybirnavirus", 
#              "Albetovirus", 
#              "Aumaivirus", 
#              "Papanivirus", 
#              "Virtovirus")
#
#duplo_phyla <- c("Uroviricota", 
#                 "Peploviricota")
#
#mono_phyla <- c("Hofneiviricota", 
#                "Phixviricota", 
#                "Cossaviricota",
#                "Cressdnaviricota", 
#                "Saleviricota")
#
#vari_phyla <- c("Preplasmiviricota",
#                "Nucleocytoviricota",
#                "Dividoviricota")
#
#adna_phyla <- c("Taleaviricota")
#
## Columns to be filled when Phylum.x is NA
#columns_to_fill <- c("Class", "Order", "Family", "Genus", "Species", "Phylum")
#
#filtered_tax2 <- filtered_tax %>%
#  # Give single blastx predictions a Realm assignment for Riboviria
#  mutate(Realm.y = ifelse(Phylum.y %in% ribo_phyla | Family.y %in% ribo_fam | Genus.y %in% ribo_gen, "Riboviria", NA)) %>%
#  # Fill in missing Realms from the blastx predictions
#  mutate(Realm = coalesce(Realm, Realm.y))
#
#filtered_tax2 <- filtered_tax2 %>%
#  # Give single blastx predictions a Realm assignment for Duplodnaviria/Monodnaviria,
#  mutate(
#    Realm.y = case_when(
#      Phylum.y %in% duplo_phyla ~ "Duplodnaviria",
#      Phylum.y %in% mono_phyla ~ "Monodnaviria",
#      Phylum.y %in% vari_phyla ~ "Varidnaviria",
#      Phylum.y %in% adna_phyla ~ "Adnaviria",
#      TRUE ~ NA
#    )
#  ) %>% 
#  mutate(Realm = coalesce(Realm, Realm.y)) %>% 
#  # Filter again on completeness
#  filter(completeness >= 50 | is.na(Realm) | Realm == "Riboviria")

filtered_tax3 <- filtered_tax2 %>% 
  # If genomad had no prediction, use blastx taxonomy
  mutate(Phylum.x = coalesce(Phylum.x, Phylum.y),
         Class.x = ifelse(Phylum.x == Phylum.y,
                          coalesce(Class.x, Class.y),
                          Class.x),
         Order.x = ifelse(Class.x == Class.y,
                          coalesce(Order.x, Order.y),
                          Order.x),
         Family.x = ifelse(Order.x == Order.y,
                           coalesce(Family.x, Family.y),
                           Family.x),
         Genus.x = ifelse(Family.x == Family.y,
                          coalesce(Genus.x, Genus.y),
                          Genus.x),
         Species.x = ifelse(Genus.x == Genus.y,
                            coalesce(Species.x, Species.y),
                            Species.x),
         ) %>% 
  mutate(across(Phylum.x:Species.x, ~ifelse(. == "unclassified", NA, .)))

taxonomy_df <- filtered_tax3 %>% 
  select(rowname, Realm,  Kingdom, Phylum.x, Class.x, Order.x, Family.x, Genus.x, Species.x) %>% 
  rename_with(~sub("\\.x$", "", .), ends_with(".x"))

#' ### Pie charts
#' Realm pie chart
realms <- unique(taxonomy_df$Realm)

realm_col <- c(rev(brewer.pal(4, "Paired"))[1:length(realms)-1], "lightgrey")
names(realm_col) <- c(realms[order(realms)][1:3], "Unclassified")

#realm_plot <- taxonomy_df %>% 
#  mutate_at(vars(Realm), ~replace_na(., "Unclassified")) %>% 
#  group_by(Realm) %>% 
#  count() %>% 
#  ungroup() %>% 
#  arrange(desc(Realm)) %>% 
#  mutate(text_y = cumsum(n) - n/2, n_label=comma(n),
#         perc=n/sum(n)) %>% 
#  ggplot(aes(x = "", y = n, fill = Realm)) +
#  geom_bar(stat = "identity", width = 1) +
#  coord_polar(theta = "y")+
#  geom_label_repel(aes(label = n_label, y = text_y), nudge_x = .57, 
#                   min.segment.length = 100, show.legend = F) +
#  scale_fill_manual(values = realm_col)+
#  labs(title = "All virus contigs")+
#  theme_void()+
#  theme(plot.title = element_text(hjust = 0.5),
#        legend.text = element_text(face="italic"))

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


#realm_read_plot <-taxonomy_df %>% 
#  left_join(as.data.frame(rowSums(OTU)) %>% rownames_to_column()) %>% 
#  rename(Abundance=`rowSums(OTU)`) %>%
#  mutate_at(vars(Realm), ~replace_na(., "Unclassified")) %>% 
#  group_by(Realm) %>% 
#  summarise(n=sum(Abundance), .groups = "drop") %>% 
#  arrange(desc(Realm)) %>% 
#  mutate(text_y = cumsum(n) - n/2, n_label=comma(n)) %>% 
#  ggplot(aes(x = "", y = n, fill = Realm)) +
#  geom_bar(stat = "identity", width = 1) +
#  coord_polar(theta = "y")+
#  geom_label_repel(aes(label = n_label, y = text_y), nudge_x = .57, 
#                   min.segment.length = 100, show.legend = F) +
#  scale_fill_manual(values = realm_col)+
#  labs(title = "All virus reads")+
#  theme_void()+
#  theme(plot.title = element_text(hjust = 0.5),
#        legend.text = element_text(face="italic"))

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

#' *Riboviria* pie chart
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

#ribophylum_plot <- riboviria_clean %>% 
#  ggplot(aes(x = "", y = n, fill = Phylum)) +
#  geom_bar(stat = "identity", width = 1) +
#  coord_polar(theta = "y")+
#  geom_label_repel(aes(label = n_label, y = text_y), nudge_x = .57, 
#                   min.segment.length = 100, show.legend = F) +
#  scale_fill_manual(values = phyla_col)+
#  labs(fill="Phylum", title = expression(paste(italic("Riboviria"), " contigs")))+
#  theme_void()+
#  theme(plot.title = element_text(hjust = 0.5),
#        legend.text = element_text(face="italic"))

ribophylum_plot <- riboviria_clean %>% 
  ggplot(aes(x = 1, y = perc, fill =  Phylum)) +
  geom_col(show.legend = T) +
  coord_polar(theta = "y")+
  xlim(c(-1.5, 1.5)) +
  geom_label_repel(aes(label = n_label), show.legend = F, 
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

#ribophylum_read_plot <- riboviria_read_abundance %>% 
#  ggplot(aes(x = "", y = n, fill = Phylum)) +
#  geom_bar(stat = "identity", width = 1) +
#  coord_polar(theta = "y")+
#  geom_label_repel(aes(label = n_label, y = text_y), nudge_x = 0.57, nudge_y=0.5,
#                   min.segment.length = 100, show.legend = F) +
#  scale_fill_manual(values = phyla_col)+
#  labs(fill="Phylum", title = expression(paste(italic("Riboviria"), " reads")))+
#  theme_void()+
#  theme(plot.title = element_text(hjust = 0.5),
#        legend.text = element_text(face="italic"))

ribophylum_read_plot <- riboviria_read_abundance %>% 
  ggplot(aes(x = 1, y = perc, fill =  Phylum)) +
  geom_col(show.legend = T) +
  coord_polar(theta = "y")+
  xlim(c(-1.5, 1.5)) +
  geom_label_repel(aes(label = n_label), show.legend = F, 
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

pie1 <- ggarrange(realm_plot, realm_read_plot,  
          labels = c('A', 'B'), common.legend = T,
          legend = "right", align="hv",
          font.label = list(size=20))

pie2 <- ggarrange(ribophylum_plot, ribophylum_read_plot, 
                  labels = c('C', 'D'), common.legend = T,
                  legend = "right", align="hv",
                  font.label = list(size=20))

ggarrange(pie1, pie2, ncol=1, align = 'hv')
ggsave("figures/piecharts.pdf", dpi=300)

#' ## Barplots
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
        strip.text = element_markdown())

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

#' ## Rarefaction curves

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

#' !! Make sure to run the RdRP tree script
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
  geom_richtext(data=. %>% select(Group, max_observation, max_r) %>% distinct(),
                aes(label=Group, x=max_observation, y=max_r), fill=NA, label.color=NA,
                show.legend = F, hjust=1, nudge_y = 25)+
  scale_color_manual(values=curve_col)+
  scale_fill_manual(values=curve_col)+
  scale_x_continuous(breaks = seq(min(curve_df$observation), max(curve_df$observation), by = 10)) +
  scale_y_continuous(breaks = seq(min(curve_df$r), max(curve_df$r)+100, by = 100)) +
  theme_bw()+
  labs(y="Observed species", x="Samples")+
  theme(legend.text = element_markdown(),
        legend.title = element_blank(),
        axis.ticks = )

fig1 <- ((realm_plot | realm_read_plot) + plot_layout(guides = "collect")) /
  ((ribophylum_plot | ribophylum_read_plot) + plot_layout(guides = "collect")) /
  rarefaction_plot

fig1 & patchwork::plot_annotation(tag_levels = "A")&
  theme(plot.tag = element_text(face = 'bold'))

ggsave("figures/figure2.pdf", dpi=300, height = 8.2, width=7)

fig2 <- ((realm_p + riboviria_p)+ plot_layout(guides = "collect"))/
  rarefaction_plot+
  plot_layout(heights = c(1, .5))

fig2 & 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = 'bold'))

#ggsave("figures/figure_2_bis.pdf", dpi=300, height = 7.5, width=7)
ggsave("figures/tiff/figure3.tiff", dpi=300, height = 7.5, width=7)

#' ### Venn diagram
genomad <- filtered_tax3[!is.na(filtered_tax3$Virus.x),]$rowname
blastx <- filtered_tax3[!is.na(filtered_tax3$Virus.y),]$rowname
hmm <- scan("data/HMM_RdRP_matches.txt", character(), quote = "")

single_blastx <- setdiff(setdiff(blastx, genomad), hmm)
single_hmm <- setdiff(setdiff(hmm, genomad), blastx)
single_genomad <- setdiff(setdiff(genomad, hmm), blastx)

single_blastx_df  <- taxonomy_df %>% 
  filter(rowname %in% single_blastx)

single_blastx_df %>% 
  filter(str_detect(Species, "(?i)bacteriophage"))

myCol <- brewer.pal(3, "Dark2")


#ggVennDiagram::ggVennDiagram(list("genomad"=genomad, "DIAMOND blastx"=blastx, "NeoRdRP hmmsearch"=hmm))+
#  scale_color_manual(values=myCol)

venn_result <- venn.diagram(
  x = list(genomad, blastx, hmm),
  category.names = c("genomad" , "DIAMOND blastx", "NeoRdRP hmmsearch"),
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

#' ## CheckV analysis

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
  geom_label_repel(aes(label = n_label, y = text_y), nudge_x = 0.57, #nudge_y=0.5,
                   min.segment.length = 100, show.legend = F) +
  scale_fill_manual(values = checkv_col,  limits=names(checkv_col))+
  labs(fill="CheckV quality")+
  theme_void()+
  theme(legend.position = "right",
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=7))

ggarrange(venn_result, checkv_pie, labels = "AUTO", nrow=1)
#ggsave("figures/venn_checkv.pdf", dpi=300, width=160, height=80, units = "mm")
ggsave("figures/tiff/figure2.tiff", dpi=300, width=160, height=80, units = "mm")

#' ### Heatmap
clean_viruses <- taxonomy_df$rowname
clean_OTU <- OTU %>% 
  filter(rownames(.) %in% clean_viruses) %>% 
  select(-contains("NC"))

joined_df <- left_join(clean_OTU %>% rownames_to_column(),
          taxonomy_df) 

grouped_order_df <- joined_df %>% 
  select(starts_with("Blackfly"), Realm,  Kingdom, Phylum, Class, Order, Family) %>% 
  group_by(Realm,  Kingdom, Phylum, Class, Order, Family) %>% 
  summarize(across(all_of(starts_with("Blackfly")), ~sum(.))) %>% 
  ungroup() %>% 
  mutate(across(!starts_with("Blackfly"), ~replace_na(., "unclassified")))


hm_df <- grouped_order_df %>% 
  mutate(rowname = paste(Realm, Kingdom, Phylum, Class, Order, Family, sep = "_")) %>%
  select(-Realm, -Kingdom, -Phylum, -Class, -Order, -Family) %>% 
  filter(rowSums(across(where(is.numeric)))!=0) %>% 
  column_to_rownames()

hm_anno_df <- hm_df %>% 
  rownames_to_column() %>%
  select(-starts_with("Blackfly")) %>% 
  separate(rowname, into = c("Realm", "Kingdom", "Phylum", "Class", "Order", "Family"), 
           sep = "_", fill = "right", remove = F) %>% 
  column_to_rownames()

heatmapCols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(200)

#heatmapCols <- viridisLite::magma(200, direction = -1)

king <- unique(hm_anno_df$Kingdom)
ordered_king <- king[order(king)]
king_col <- c(brewer.pal(length(king)-1, "Set2"), "lightgrey")
names(king_col) <- c(ordered_king[ordered_king!="unclassified"], "unclassified")


phyla <- unique(taxonomy_df[taxonomy_df$Realm=="Riboviria",]$Phylum)

phyla_col <- c(brewer.pal(10, "Paired")[5:10], "lightgrey")
names(phyla_col) <- c(phyla[order(phyla)][-length(phyla)], "unclassified")

phyla <- unique(hm_anno_df$Phylum)
ordered_phyla <- phyla[order(phyla)]

nonribo_phyla <- setdiff(phyla, names(phyla_col))
ordered_nonribo_phyla <- nonribo_phyla[order(nonribo_phyla)]
phyla_col <- c(phyla_col, c(brewer.pal(4, "Paired")))
names(phyla_col)[8:11] <- ordered_nonribo_phyla

phyla_col_order <- phyla_col[names(phyla_col) != "unclassified"]
phyla_col <- c(phyla_col_order[order(names(phyla_col_order))], "unclassified"="lightgrey")


class <- unique(hm_anno_df$Class)
ordered_class <- class[order(class)]
#class_col <- c(viridis(length(class)-1), "lightgrey")
class_col <- c(rep(brewer.pal(12, "Paired"), length.out=length(class)-1), "lightgrey")
names(class_col) <- c(ordered_class[ordered_class!="unclassified"], "unclassified")

class_pch <- rep(1:10, length.out=length(class)-1)
names(class_pch) <- c(ordered_class[ordered_class!="unclassified"])

order <- unique(hm_anno_df$Order)
ordered_order <- order[order(order)]
#order_col <- c(plasma(length(order)-1), "lightgrey")
order_col <- c(rep(brewer.pal(12, "Paired"), length.out=length(order)-1), "lightgrey")
names(order_col) <- c(ordered_order[ordered_order!="unclassified"], "unclassified")

order_pch <- rep(12:20, length.out=length(order)-1)
names(order_pch) <- c(ordered_order[ordered_order!="unclassified"])

hm_anno_df <- hm_anno_df %>% 
  mutate(across(everything(), ~ifelse(.=="unclassified", NA, .)),
         ClassNumber = class_pch[Class],  OrderNumber=order_pch[Order],
         Family_all = ifelse(is.na(Family),
                             glue("unclassified <i>{coalesce(Order, Class, Phylum, Kingdom, Realm)}</i>"),
                             glue("<i>{Family}</i>")),
         Family_all=ifelse(Family_all == "unclassified <i>NA</i>", "unclassified", Family_all),
         across(where(is.character), ~replace_na(., "unclassified"))
  )

ht_opt$message = FALSE

left_ra <- rowAnnotation(
                      'Kingdom'=anno_simple(hm_anno_df$Kingdom, col=king_col),
                      'Phylum'=anno_simple(hm_anno_df$Phylum, col=phyla_col),
                      'Class'=anno_simple(hm_anno_df$Class, col=class_col, pch=hm_anno_df$ClassNumber, 
                                          pt_size = unit(2, "mm"), pt_gp=gpar(col="black")),
                      'Order'=anno_simple(hm_anno_df$Order, col=order_col, pch=hm_anno_df$OrderNumber, 
                                          pt_size = unit(2, "mm"), pt_gp=gpar(col="black")),
                      show_annotation_name=T,
                      annotation_name_side = "top",
                      annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
                      annotation_name_rot = 45,
                      annotation_legend_param=list(labels_gp = gpar(fontface="italic"))
                      )

hm <- Heatmap(
        log2(as.matrix(hm_df+1)), 
        cluster_columns = T,
        clustering_distance_columns = vegdist(t(hm_df), method = "bray"),
        cluster_rows = F, 
        col = heatmapCols, 
        name = "Log2 Read Counts", 
        heatmap_legend_param = list(direction = "horizontal"),
        show_column_names = F,
        left_annotation = left_ra,
        row_labels = gt_render(hm_anno_df$Family_all),
        row_names_gp = gpar(fontsize=8),
        split=hm_anno_df$Realm,
        row_title_gp = gpar(fontsize=10),
        row_title_rot=0,
        row_title = gt_render(ifelse(unique(hm_anno_df$Realm) != "unclassified",
                                     glue("<i>{unique(hm_anno_df$Realm)}</i>"),
                                     "unclassified")),
        rect_gp = gpar(col = "white", lwd = .7),
        border=T
        )

lgd_king <- Legend(labels = gt_render(ifelse(names(king_col)!="unclassified", glue("<i>{names(king_col)}</i>"), "unclassified")), 
                   title = "Kingdom", legend_gp=gpar(fill=king_col))
lgd_phylum <- Legend(labels = gt_render(ifelse(names(phyla_col)!="unclassified", glue("<i>{names(phyla_col)}</i>"), "unclassified")), 
                     title="Phylum", legend_gp=gpar(fill=phyla_col))
lgd_class <- Legend(labels = gt_render(ifelse(names(class_col)!="unclassified", glue("<i>{names(class_col)}</i>"), "unclassified")), 
                    title = "Class", legend_gp=gpar(col="black"), background = class_col, pch=class_pch, type = "points")
lgd_order <- Legend(gt_render(ifelse(names(order_col)!="unclassified", glue("<i>{names(order_col)}</i>"), "unclassified")),
                    title="Order", legend_gp=gpar(col="black"), background = order_col, pch=order_pch, type = "points")


pd <- packLegend(lgd_king, lgd_phylum, lgd_class, lgd_order, max_height = unit(10, "cm"))

draw(hm, annotation_legend_list = list(lgd_king, lgd_phylum, lgd_class, lgd_order), 
     annotation_legend_side = "left", heatmap_legend_side = "left",
     merge_legend = T, legend_grouping="original")

pdf("figures/blackflies_heatmap_top_legend.pdf", width = 16.6, height =17)
draw(hm, annotation_legend_list = pd, 
     annotation_legend_side = "bottom", heatmap_legend_side = "bottom",
     merge_legend = T, legend_grouping="original")
dev.off()

tiff("figures/tiff/figure4.tiff", width = 16.6, height =17, units = "in", res=300)
draw(hm, annotation_legend_list = pd, 
     annotation_legend_side = "bottom", heatmap_legend_side = "bottom",
     merge_legend = T, legend_grouping="original")
dev.off()

# Eukaryotic viral families
hm_anno_df %>% 
  filter(Family != "unclassified", Realm!="Duplodnaviria", Class!="Leviviricetes", !Family %in% c("Picobirnaviridae", "Cystoviridae")) %>% 
  count()

# Eukaryotic viruses not classified onto Family level
hm_anno_df %>% 
  filter(Family == "unclassified", Realm!="Duplodnaviria", Class!="Leviviricetes", !Family %in% c("Picobirnaviridae", "Cystoviridae")) %>% 
  count()

# Prokaryotic viruses
hm_anno_df %>% 
  filter(Realm=="Duplodnaviria" | Class=="Leviviricetes" | Family %in% c("Picobirnaviridae", "Cystoviridae")) %>% 
  count()

# Classified eukaryotic RNA viruses
hm_anno_df %>% 
  filter(Realm == "Riboviria", Family != "unclassified", Class!="Leviviricetes", !Family %in% c("Picobirnaviridae", "Cystoviridae")) %>% 
  count()

#' ## Correlation analysis
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


# #' ### Make a phyloseq object
# OTU.UF <- otu_table(as.matrix(OTU), taxa_are_rows=T)
# tax.UF <- tax_table(as.matrix(filtered_tax2 %>% column_to_rownames(var = "rowname")))
# meta.UF <- sample_data(meta)
# 
# BF <- phyloseq(OTU.UF, tax.UF, meta.UF)
# BF
# 
# #' ### Remove contamination
# #' #### Visualize library sizes of samples and negative controls
# decontam <- as.data.frame(sample_data(BF))
# decontam   
# decontam$LibrarySize <- sample_sums(BF)
# decontam <- decontam[order(decontam$LibrarySize),]
# decontam$Index <- seq(nrow(decontam))
# ggplot(data=decontam, aes(x=Index, y=LibrarySize, color=Control)) + 
#   geom_point()+
#   ggtitle("Library sizes")
# 
# #' #### Detect contaminants
# sample_data(BF)$is.neg <- sample_data(BF)$Control == "Yes"
# sample_data(BF)$is.neg
# contamdf.prev <- isContaminant(BF, method="prevalence", neg="is.neg")
# 
# 
# #' **Number of contaminating contigs:**
# table(contamdf.prev$contaminant)
# 
# contamdf.prev[contamdf.prev$contaminant==T,]
# 
# #' #### Visualize prevalence of contaminants in samples and negative controls
# ps.pa <- transform_sample_counts(BF, function(abund) 1*(abund>0))
# ps.pa.neg <- prune_samples(sample_data(ps.pa)$Control == "Yes", ps.pa)
# ps.pa.pos <- prune_samples(sample_data(ps.pa)$Control == "No", ps.pa)
# 
# df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
#                     contaminant=contamdf.prev$contaminant)
# df.pa
# ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
#   xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")+
#   ggtitle("Prevalence of contaminants in samples vs NCs")
# 
# #' #### Remove negative controls and contaminants from phyloseq object
# ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, BF)
# ps.noncontam
# ps.noncontam <- prune_samples(sample_data(BF)$Control!='Yes', ps.noncontam)
# BF <- ps.noncontam
# BF
# 
# #' Remove negative controls from metadata table
# meta <- meta[meta$Control!="Yes",]
# 
# #' #### Subset only viruses and remove prokaryotic viruses and EVEs
# #BF.V <- subset_taxa(BF, Kingdom=="Viruses")
# BF.V <- BF
# 
# #EVE_phage <- c("Atrato Retro-like virus", "Gurupi chuvirus-like 1", "Aedes aegypti To virus 1",
# #               "Aedes aegypti To virus 2", "Guato virus", "Kaiowa virus", "Atrato Chu-like virus 1",
# #               "Chuvirus Mos8Chu0", "Chibugado virus", "Prokaryotic dsDNA virus sp.", "Bacteriophage sp.")
# 
# #BF.V2 <- subset_taxa(BF.V, !is.element(Species, EVE_phage))
# #BF.V2 <- subset_taxa(BF.V2, Order!="Caudovirales")
# #BF.V2 <- subset_taxa(BF.V2, Family!="Microviridae")
# #BF.V2 <- subset_taxa(BF.V2, Phylum!="Phage")
# #BF.V2
# BF.V2 <- BF.V
# BF.V2
# 
# #' **Info on sample variables and level of taxonomic ranks:**
# sample_variables(BF.V2)
# rank_names(BF.V2)
# 
# #' #### Agglomerate taxa on viral species level
# BF_species <- BF.V2 %>%
#   tax_glom(taxrank = "Realm")
# BF_species
# 
# #' **Only keep taxa with more than 0 reads:**
# BF_species <- prune_taxa(taxa_sums(BF_species) > 0, BF_species)
# BF_species
# 
# #' **Only keep samples with more than 0 viral reads:**
# BF_final <- prune_samples(sample_sums(BF_species) > 0, BF_species)
# BF_final
# 
# #vtax<-as.data.frame(tax_table(BF_final))
# 
# #' **Convert phyloseq object to metagenomeseq object to make heatmap:**
# BF_metaseq <- phyloseq_to_metagenomeSeq(BF_final)
# 
# #' **Aggregate by species:**
# BF_metaseq_species <- aggregateByTaxonomy(BF_metaseq, lvl = "Realm", norm = F, 
#                                              aggfun = colSums, out = "MRexperiment", alternate = T)
# 
# #' **Count number of unique viral species:**
# n_species <- length(unique(featureData(BF_metaseq_species)$Realm))
# 
# #' **Set heatmap colors:**
# heatmapCols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(20)
# 
# #' #' **Assign mosquito species for each sample:**
# #' mosquito_species <- pData(BF_metaseq_species)$Species
# 
# #' #' **Assign colors to location:**
# #' location <- pData(BF_metaseq_species)$Location
# 
# # #Draw heatmap
# 
# #n depends on the amount of taxa in 'BF_metaseq_species'
# hm <- plot_abundance(BF_metaseq_species, n = n_species, log = T, norm = F, colclust = "bray",
#                      col = heatmapCols,
#                      name = "Log2 Read Counts", 
#                      # top_annotation=column_ha,
#                      row_names_gp = gpar(fontsize = 4),
#                      show_column_names = T,
#                      column_names_rot = 45,
#                      column_names_gp=gpar(fontface=1, fontsize=10),
#                      #heatmap_legend_param = list(direction = "horizontal"),
#                      cluster_rows = T,
#                      row_title_gp = gpar(fontsize=10),
#                      row_title_rot=0,
#                      border=F)
# 
# #+ echo=TRUE, fig.width=15, fig.height=9
# draw(hm, heatmap_legend_side = "left", #annotation_legend_side = "left", 
#      merge_legend = T, legend_grouping="original")
# 
# pdf("blackflies_heatmap.pdf", width = 15, height =8.27)
# draw(hm, heatmap_legend_side = "left", annotation_legend_side = "left", 
#      merge_legend = T, legend_grouping="original")
# dev.off()
# 
# #' ## Heatmap on Family level
# #' **Aggregate by order:**
# BF_metaseq_family <- aggregateByTaxonomy(BF_metaseq, lvl = "Family", norm = F, 
#                                         aggfun = colSums, out = "MRexperiment", alternate = T)
# 
# #' **Count number of unique viral family:**
# n_family <- length(unique(featureData(BF_metaseq_family)$Family))
# 
# #' **Set heatmap colors:**
# heatmapCols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(20)
# 
# # #Draw heatmap
# 
# #n depends on the amount of taxa in 'BF_metaseq_family'
# hm <- plot_abundance(BF_metaseq_family, n = n_family, log = T, norm = F, colclust = "bray",
#                      col = heatmapCols,
#                      name = "Log2 Read Counts", 
#                      # top_annotation=column_ha,
#                      row_names_gp = gpar(fontsize = 8),
#                      show_column_names = T,
#                      column_names_rot = 45,
#                      column_names_gp=gpar(fontface=1, fontsize=10),
#                      heatmap_legend_param = list(direction = "horizontal"),
#                      cluster_rows = F,
#                      row_title_gp = gpar(fontsize=10),
#                      row_title_rot=0,
#                      border=F)
# 
# #+ echo=TRUE, fig.width=15, fig.height=9
# draw(hm, heatmap_legend_side = "left", #annotation_legend_side = "left", 
#      merge_legend = T, legend_grouping="original")
# 
# 
# 
# #' ## Heatmap on Order level
# #' **Aggregate by order:**
# BF_metaseq_order <- aggregateByTaxonomy(BF_metaseq, lvl = "Order", norm = F, 
#                                           aggfun = colSums, out = "MRexperiment", alternate = T)
# 
# #' **Count number of unique viral family:**
# n_order <- length(unique(featureData(BF_metaseq_order)$Order))
# 
# #' **Set heatmap colors:**
# heatmapCols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(20)
# 
# #' #' **Assign mosquito family for each sample:**
# #' mosquito_family <- pData(BF_metaseq_family)$family
# 
# #' #' **Assign colors to location:**
# #' location <- pData(BF_metaseq_family)$Location
# 
# # #Draw heatmap
# 
# #n depends on the amount of taxa in 'BF_metaseq_family'
# hm <- plot_abundance(BF_metaseq_order, n = n_order, log = T, norm = F, colclust = "bray",
#                      col = heatmapCols,
#                      name = "Log2 Read Counts", 
#                      # top_annotation=column_ha,
#                      row_names_gp = gpar(fontsize = 8),
#                      show_column_names = T,
#                      column_names_rot = 45,
#                      column_names_gp=gpar(fontface=1, fontsize=10),
#                      heatmap_legend_param = list(direction = "horizontal"),
#                      cluster_rows = F,
#                      row_title_gp = gpar(fontsize=10),
#                      row_title_rot=0,
#                      border=F)
# 
# #+ echo=TRUE, fig.width=15, fig.height=9
# draw(hm, heatmap_legend_side = "left", #annotation_legend_side = "left", 
#      merge_legend = T, legend_grouping="original")
# 
# 
# #' ## Onchocercidae
# BF_oncho <- subset_taxa(BF, Family=="Onchocercidae")
# BF_oncho <- tax_glom(BF_oncho, taxrank = "Family")
# BF_oncho <- prune_samples(sample_data(BF_oncho)$Control!='Yes', BF_oncho)
# BF_oncho <- prune_samples(sample_sums(BF_oncho)>0, BF_oncho)
# BF_oncho
# 
# oncho<-as.data.frame(otu_table(BF_oncho))
# oncho <- oncho %>% 
#   pivot_longer(everything(), names_to = "Sample", values_to = "Abundance")
# 
# oncho <- oncho[order(oncho$Abundance),]
# 
# p1<-ggplot(oncho, aes(x = Sample, y=Abundance))+
#   geom_col()+
#   xlab('')+
#   ggtitle("Onchocercidae")+
#   theme(axis.text.x = element_text(angle=90))
# 
# #' ## Fungi
# BF_fungi <- subset_taxa(BF,Phylum %in% c("Ascomycota", "Basidiomycota", 
#                                                      "Chytridiomycota", "Glomeromycota",
#                                                      "Zygomycota","Neocallimastigomycota",
#                                                      "Microsporidia"))
# BF_fungi <- tax_glom(BF_fungi, taxrank = "Phylum")
# BF_fungi <- prune_samples(sample_data(BF_fungi)$Control!='Yes', BF_fungi)
# BF_fungi <- prune_samples(sample_sums(BF_fungi)>0, BF_fungi)
# BF_fungi
# 
# fungi<-as.data.frame(otu_table(BF_fungi))
# fungi <- fungi %>% 
#   mutate(Phylum=rownames(.), .before=1) %>% 
#   pivot_longer(cols = 2:55, names_to="Sample", values_to = "Abundance") %>% 
#   pivot_wider(names_from = Phylum, values_from = Abundance)
# 
# colnames(fungi)<-c("Sample", as.data.frame(tax_table(BF_fungi))$Phylum)
# 
# fungi<-fungi %>% 
#   pivot_longer(cols = 2:5, names_to = "Fungi", values_to = "Abundance")
# 
# p2<-ggplot()+
#   geom_col(data = fungi, aes(x = Sample, y=Abundance, fill=Fungi))+
#   xlab('')+
#   ggtitle('Fungi')+
#   theme(axis.text.x = element_text(angle=90))
# p2
# 
# BF_fungi2 <- tax_glom(BF_fungi, taxrank = "Kingdom")
# fungi2<-as.data.frame(otu_table(BF_fungi2))
# fungi2 <- fungi2 %>% 
#   pivot_longer(everything(), names_to = "Sample", values_to = "Abundance")
# 
# #' ## Viruses
# BF_virus <- subset_taxa(BF, Kingdom=='Viruses')
# BF_virus <- tax_glom(BF_virus, taxrank = "Kingdom")
# BF_virus <- prune_samples(sample_data(BF_virus)$Control!='Yes', BF_virus)
# BF_virus <- prune_samples(sample_sums(BF_virus)>0, BF_virus)
# BF_virus
# 
# virus<-as.data.frame(otu_table(BF_virus))
# virus <- virus %>% 
#   pivot_longer(everything(), names_to = "Sample", values_to = "Abundance")
# 
# p3<-ggplot(virus, aes(x = Sample, y=Abundance))+
#   geom_col()+
#   xlab('')+
#   ggtitle("Virus")+
#   theme(axis.text.x = element_text(angle=90))
# p3
# 
# final<-left_join(virus, oncho, by="Sample") %>% 
#   left_join(., fungi2, by="Sample")%>% 
#   pivot_longer(cols=2:4, names_to = "Type", values_to = "Abundance") %>% 
#   mutate(Type=case_when(Type=="Abundance.x" ~ "Virus",
#                         Type=="Abundance.y" ~ "Onchocercidae",
#                         Type=="Abundance" ~ "Fungi"))
# 
# 
# final %>% 
#   mutate(Sample = factor(as.character(final$Sample), 
#                        levels = gtools::mixedsort(as.character(final$Sample))))
# 
# 
# ggplot(final, aes(x=Sample, y=log2(Abundance)))+
#   geom_col()+
#   facet_wrap(~Type, ncol=1)+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle=90))
# ggsave("read-overview.pdf")
