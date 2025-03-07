source("blackflies.R")
here::i_am("heatmap.R")
load("blackflies.RData")
# Prep
clean_viruses <- taxonomy_df$rowname
clean_OTU <- OTU %>%
  filter(rownames(.) %in% clean_viruses) %>%
  select(-contains("NC"))

joined_df <- left_join(
  clean_OTU %>% rownames_to_column(),
  taxonomy_df
)

tpm_scaled <- joined_df %>%
  mutate(length = as.numeric(str_extract(rowname, "(?<=length_)\\d+"))) %>%
  mutate(
    across(starts_with("BlackFly"), ~ .x / (length / 1000)),
    across(starts_with("BlackFly"), ~ .x / sum(.x))
  )
# Raw counts
ra_p <- joined_df %>%
  mutate(across(starts_with("BlackFly"), ~ .x / sum(.x))) %>%
  select(starts_with("BlackFly"), Realm, rowname) %>%
  pivot_longer(cols = starts_with("Blackfly"), names_to = "Sample", values_to = "Read_Count") %>%
  filter(Read_Count > 0) %>%
  mutate(Realm = if_else(is.na(Realm), "Unclassified", Realm)) %>%
  ggplot(aes(x = Realm, y = Read_Count, color = Realm, fill = Realm)) +
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = .5,
    ## adjust height
    width = .2,
    ## move geom to the right
    justification = -.6,
    ## remove slab interval
    .width = 0,
    point_colour = NA,
    show.legend = F
  ) +
  geom_boxplot(
    width = .15,
    ## remove outliers
    outlier.shape = NA, ## `outlier.shape = NA` or `outlier.alpha = 0` works as well
  ) +
  ## add dot plots from {ggdist} package
  ggdist::stat_dots(
    ## orientation to the left
    side = "left",
    ## move geom to the left
    justification = 1.12,
    ## adjust grouping (binning) of observations
    binwidth = 0.03,
    show.legend = F
  ) +
  scale_y_log10(
    breaks = c(0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1),
    labels = c(0, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100)
  ) +
  labs(y = "Relative abundance (%)") +
  # ggpubr::stat_pwc(hide.ns = T)+
  scale_fill_manual(values = alpha(realm_col, .3)) +
  scale_color_manual(values = realm_col) +
  theme_classic() +
  coord_cartesian(xlim = c(1.2, 3.9), clip = "off")

# TPM scaled
tpm_p <- tpm_scaled %>%
  select(starts_with("BlackFly"), Realm, rowname) %>%
  pivot_longer(cols = starts_with("Blackfly"), names_to = "Sample", values_to = "TPM") %>%
  # group_by(rowname) %>%
  # filter(TPM == max(TPM)) %>%
  # ungroup() %>%
  filter(TPM > 0) %>%
  mutate(Realm = if_else(is.na(Realm), "Unclassified", Realm)) %>%
  ggplot(aes(x = Realm, y = TPM, color = Realm, fill = Realm)) +
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = .5,
    ## adjust height
    width = .2,
    ## move geom to the right
    justification = -.6,
    ## remove slab interval
    .width = 0,
    point_colour = NA,
    show.legend = F
  ) +
  geom_boxplot(
    width = .15,
    ## remove outliers
    outlier.shape = NA, ## `outlier.shape = NA` or `outlier.alpha = 0` works as well
  ) +
  ## add dot plots from {ggdist} package
  ggdist::stat_dots(
    ## orientation to the left
    side = "left",
    ## move geom to the left
    justification = 1.12,
    ## adjust grouping (binning) of observations
    binwidth = 0.025,
    show.legend = F
  ) +
  scale_y_log10(
    breaks = c(0, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1),
    labels = c(0, 0.001, 0.01, 0.1, 1, 10, 100)
  ) +
  labs(y = "Length normalized\nrelative abundance (%)") +
  # ggpubr::stat_pwc(hide.ns = T, method = "wilcox.test", p.adjust.method = "BH") +
  scale_fill_manual(values = alpha(realm_col, .3)) +
  scale_color_manual(values = realm_col) +
  theme_classic() +
  coord_cartesian(xlim = c(1.2, 3.9), clip = "off")

tpm_p / ra_p
ggsave("figures/raincloud.pdf", dpi = 300, height = 5, width = 10)

tpm_p
ggsave(plot = tpm_p, filename = "figures/tpm_realm.pdf", dpi = 300, height = 5, width = 10)

# Heatmap
source("draw_blackfly_heatmap.R")

# Raw counts
grouped_order_df <- joined_df %>%
  select(starts_with("BlackFly"), Realm, Kingdom, Phylum, Class, Order, Family) %>%
  group_by(Realm, Kingdom, Phylum, Class, Order, Family) %>%
  summarize(across(all_of(starts_with("BlackFly")), ~ sum(.)),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(across(!starts_with("BlackFly"), ~ replace_na(., "unclassified")))

draw_blackfly_heatmap(grouped_order_df)

# TPM
tpm_order_df <- tpm_scaled %>%
  select(starts_with("BlackFly"), Realm, Kingdom, Phylum, Class, Order, Family) %>%
  group_by(Realm, Kingdom, Phylum, Class, Order, Family) %>%
  summarize(across(all_of(starts_with("BlackFly")), ~ sum(.)),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(across(!starts_with("BlackFly"), ~ replace_na(., "unclassified")))

draw_blackfly_heatmap(tpm_order_df, log2 = F, title = "Length normalized RA", legend_side = "bottom")

pdf("figures/blackflies_heatmap_left_legend.pdf", width = 16.6, height = 10)
draw_blackfly_heatmap(grouped_order_df, legend_side = "left")
dev.off()

tiff("figures/tiff/figure3.tiff", width = 16.6, height = 17, units = "in", res = 300)
draw(hm,
  annotation_legend_list = pd,
  annotation_legend_side = "bottom", heatmap_legend_side = "bottom",
  merge_legend = T, legend_grouping = "original"
)
dev.off()

# Eukaryotic viral families
hm_anno_df %>%
  filter(Family != "unclassified", Realm != "Duplodnaviria", Class != "Leviviricetes", !Family %in% c("Picobirnaviridae", "Cystoviridae")) %>%
  count()

# Eukaryotic viruses not classified onto Family level
hm_anno_df %>%
  filter(Family == "unclassified", Realm != "Duplodnaviria", Class != "Leviviricetes", !Family %in% c("Picobirnaviridae", "Cystoviridae")) %>%
  count()

# Prokaryotic viruses
hm_anno_df %>%
  filter(Realm == "Duplodnaviria" | Class == "Leviviricetes" | Family %in% c("Picobirnaviridae", "Cystoviridae")) %>%
  count()

# Classified eukaryotic RNA viruses
hm_anno_df %>%
  filter(Realm == "Riboviria", Family != "unclassified", Class != "Leviviricetes", !Family %in% c("Picobirnaviridae", "Cystoviridae")) %>%
  count()

hm_anno_df %>%
  filter(Realm == "Riboviria", Family != "unclassified", Class != "Leviviricetes", !Family %in% c("Picobirnaviridae", "Cystoviridae")) %>%
  count()
