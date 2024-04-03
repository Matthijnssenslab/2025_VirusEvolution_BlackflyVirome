source("blackflies.R")

unclassified_seq <- taxonomy_df %>% 
  #mutate_at(vars(Realm), ~replace_na(., "Unclassified")) %>%
  filter(is.na(Realm))

unclassified_seq <- unclassified_seq %>% 
  mutate(length=as.integer(str_split_i(rowname, "_", 4)),
         name=sub("^(?:[^_]*_){6}(.*)$", "\\1", rowname))

#unclassified_seq %>%
#  select(rowname) %>% 
#  write_delim("data/unclassified.txt")

clusters <- read_delim("data/unclassified_clusters.tsv", 
                       col_names = c("representative", "cluster"))

unclassified_joined <- left_join(unclassified_seq, clusters, by = join_by(rowname==representative))

unclassified_joined %>% 
  filter(!is.na(Family)) %>% 
  group_by(Realm, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  summarise(rownames=paste0(rowname, collapse = ","))

unclassified_joined %>% 
  ggplot(aes(x=as.integer(length), fill=Family))+
  geom_bar()+
  scale_x_binned(limits = c(1000, 11000))+
  labs(x="length (bp)")+
  theme_bw()

family_joined <- unclassified_joined %>% 
  filter(!is.na(Family)) %>% 
  tidyr::expand(nesting(select(., -cluster)), cluster) %>% 
  left_join(unclassified_joined2 %>% select(-name), by="cluster") %>% 
  mutate(name_in_cluster = str_detect(cluster, name)) %>% 
  filter(name_in_cluster==T,
         Family.x == Family.y) %>% 
  group_by(rowname.x) %>%
  mutate(keep_row = !rowname.x %in% lead(rowname.y)) %>%
  ungroup() %>%
  filter(keep_row) %>%
  select(-keep_row, -ends_with(".y"), -starts_with("name"), -cluster) %>% 
  distinct() %>% 
  rename_with(~str_remove(., "\\.x$"), ends_with(".x"))

species_joined <- unclassified_joined %>% 
  filter(!is.na(Species) & is.na(Family)) %>% 
  tidyr::expand(nesting(select(., -cluster)), cluster) %>% 
  left_join(unclassified_joined2 %>% select(-name), by="cluster") %>% 
  mutate(name_in_cluster = str_detect(cluster, name)) %>% 
  filter(name_in_cluster==T,
         Species.x == Species.y) %>% 
  group_by(rowname.x) %>%
  mutate(keep_row = !rowname.x %in% lead(rowname.y)) %>%
  ungroup() %>%
  filter(keep_row) %>%
  select(-keep_row, -ends_with(".y"), -starts_with("name"), -cluster) %>% 
  distinct() %>% 
  rename_with(~str_remove(., "\\.x$"), ends_with(".x"))

all_unclassified <- full_join(family_joined, species_joined) %>% 
  full_join(unclassified_joined %>% 
              filter(is.na(Species) & is.na(Family)) %>% 
              select(-name, -cluster)) %>% 
  select(-length)

taxonomy_df_unique_species <- taxonomy_df %>% 
  filter(!is.na(Realm)) %>% 
  bind_rows(all_unclassified) 
