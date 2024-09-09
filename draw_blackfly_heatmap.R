library(ComplexHeatmap)
draw_blackfly_heatmap <- function(grouped_order_df, log2 = T, title = "Log2 Read counts", legend_side="left") {
  hm_df <- grouped_order_df %>%
    mutate(rowname = paste(Realm, Kingdom, Phylum, Class, Order, Family, n, sep = "_")) %>%
    select(-Realm, -Kingdom, -Phylum, -Class, -Order, -Family, -n) %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    column_to_rownames()

  hm_anno_df <- hm_df %>%
    rownames_to_column() %>%
    select(-starts_with("BlackFly")) %>%
    separate(rowname,
      into = c("Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "n"),
      sep = "_", fill = "right", remove = F
    ) %>%
    column_to_rownames()

  heatmapCols <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(200)

  # heatmapCols <- viridisLite::magma(200, direction = -1)

  king <- unique(hm_anno_df$Kingdom)
  ordered_king <- king[order(king)]
  king_col <- c(RColorBrewer::brewer.pal(length(king) - 1, "Set2"), "lightgrey")
  names(king_col) <- c(ordered_king[ordered_king != "unclassified"], "unclassified")


  phyla <- unique(taxonomy_df[taxonomy_df$Realm == "Riboviria", ]$Phylum)

  phyla_col <- c(RColorBrewer::brewer.pal(10, "Paired")[5:10], "lightgrey")
  names(phyla_col) <- c(phyla[order(phyla)][-length(phyla)], "unclassified")

  phyla <- unique(hm_anno_df$Phylum)
  ordered_phyla <- phyla[order(phyla)]

  nonribo_phyla <- setdiff(phyla, names(phyla_col))
  ordered_nonribo_phyla <- nonribo_phyla[order(nonribo_phyla)]
  phyla_col <- c(phyla_col, c(RColorBrewer::brewer.pal(12, "Paired")[c(1:4, 11:12)]))
  names(phyla_col)[8:13] <- ordered_nonribo_phyla

  phyla_col_order <- phyla_col[names(phyla_col) != "unclassified"]
  phyla_col <- c(phyla_col_order[order(names(phyla_col_order))], "unclassified" = "lightgrey")


  class <- unique(hm_anno_df$Class)
  ordered_class <- class[order(class)]
  # class_col <- c(viridis(length(class)-1), "lightgrey")
  class_col <- c(rep(RColorBrewer::brewer.pal(12, "Paired"), length.out = length(class) - 1), "lightgrey")
  names(class_col) <- c(ordered_class[ordered_class != "unclassified"], "unclassified")

  class_pch <- rep(1:10, length.out = length(class) - 1)
  names(class_pch) <- c(ordered_class[ordered_class != "unclassified"])

  order <- unique(hm_anno_df$Order)
  ordered_order <- order[order(order)]
  # order_col <- c(plasma(length(order)-1), "lightgrey")
  order_col <- c(rep(RColorBrewer::brewer.pal(12, "Paired"), length.out = length(order) - 1), "lightgrey")
  names(order_col) <- c(ordered_order[ordered_order != "unclassified"], "unclassified")

  order_pch <- rep(12:20, length.out = length(order) - 1)
  names(order_pch) <- c(ordered_order[ordered_order != "unclassified"])

  type_col <- c(
    "Eukaryotic virus" = "black",
    "Bacteriophage" = "grey60",
    "unknown" = "grey90"
  )

  hm_anno_df <- hm_anno_df %>%
    mutate(across(everything(), ~ ifelse(. == "unclassified", NA, .)),
      ClassNumber = class_pch[Class], OrderNumber = order_pch[Order],
      Family_all = ifelse(is.na(Family),
        glue::glue("unclassified <i>{coalesce(Order, Class, Phylum, Kingdom, Realm)}</i> (n={n})"),
        glue::glue("<i>{Family}</i> (n={n})")
      ),
      Family_all = ifelse(grepl("unclassified <i>NA</i>", Family_all), gsub(" <i>NA</i>", "", x = Family_all), Family_all),
      across(where(is.character), ~ replace_na(., "unclassified"))
    )

  phages <- hm_anno_df %>%
    filter(Realm == "Duplodnaviria" | Class == "Leviviricetes" | Family %in% c("Picobirnaviridae", "Cystoviridae") | Phylum == "Hofneiviricota") %>%
    pull(Family_all)

  hm_anno_df <- hm_anno_df %>%
    mutate(type = case_when(
      Family_all %in% phages ~ "Bacteriophage",
      Kingdom == "unclassified" & Family == "unclassified" ~ "unknown",
      T ~ "Eukaryotic virus"
    ))

  ht_opt$message <- FALSE

  left_ra <- ComplexHeatmap::rowAnnotation(
    "Kingdom" = anno_simple(hm_anno_df$Kingdom, col = king_col),
    "Phylum" = anno_simple(hm_anno_df$Phylum, col = phyla_col),
    "Class" = anno_simple(hm_anno_df$Class,
      col = class_col, pch = hm_anno_df$ClassNumber,
      pt_size = unit(2, "mm"), pt_gp = gpar(col = "black")
    ),
    "Order" = anno_simple(hm_anno_df$Order,
      col = order_col, pch = hm_anno_df$OrderNumber,
      pt_size = unit(2, "mm"), pt_gp = gpar(col = "black")
    ),
    show_annotation_name = T,
    annotation_name_side = "top",
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    annotation_name_rot = 45,
    annotation_legend_param = list(labels_gp = gpar(fontface = "italic"))
  )

  right_ra <- ComplexHeatmap::rowAnnotation(
    "Virus type" = anno_simple(hm_anno_df$type,
      col = type_col,
      width = unit(1, "mm")
    ),
    show_annotation_name = F,
    annotation_name_side = "top",
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    annotation_name_rot = 45,
    annotation_legend_param = list(labels_gp = gpar(fontface = "italic"))
  )

  if (log2 == T) {
    matrix <- log2(as.matrix(hm_df+1))
  }
  
  else {
    matrix <- as.matrix(hm_df)
  }
  
  hm <- ComplexHeatmap::Heatmap(
    matrix,
    cluster_columns = T,
    clustering_distance_columns = vegan::vegdist(t(hm_df), method = "bray"),
    cluster_rows = F,
    col = heatmapCols,
    # name = "Log2 Read Counts",
    name = title,
    heatmap_legend_param = list(direction = "horizontal"),
    show_column_names = F,
    left_annotation = left_ra,
    right_annotation = right_ra,
    row_labels = gt_render(hm_anno_df$Family_all),
    row_names_gp = gpar(fontsize = 8),
    split = hm_anno_df$Realm,
    row_title_gp = gpar(fontsize = 10),
    row_title_rot = 0,
    row_title = gt_render(ifelse(unique(hm_anno_df$Realm) != "unclassified",
      glue::glue("<i>{unique(hm_anno_df$Realm)}</i>"),
      "unclassified"
    )),
    rect_gp = gpar(col = "white", lwd = .7),
    border = T
  )

  lgd_king <- ComplexHeatmap::Legend(
    labels = gt_render(ifelse(names(king_col) != "unclassified", glue::glue("<i>{names(king_col)}</i>"), "unclassified")),
    title = "Kingdom", legend_gp = gpar(fill = king_col)
  )
  lgd_phylum <- ComplexHeatmap::Legend(
    labels = gt_render(ifelse(names(phyla_col) != "unclassified", glue::glue("<i>{names(phyla_col)}</i>"), "unclassified")),
    title = "Phylum", legend_gp = gpar(fill = phyla_col)
  )
  lgd_class <- ComplexHeatmap::Legend(
    labels = gt_render(ifelse(names(class_col) != "unclassified", glue::glue("<i>{names(class_col)}</i>"), "unclassified")),
    title = "Class", legend_gp = gpar(col = "black"), background = class_col, pch = class_pch, type = "points"
  )
  lgd_order <- ComplexHeatmap::Legend(gt_render(ifelse(names(order_col) != "unclassified", glue::glue("<i>{names(order_col)}</i>"), "unclassified")),
    title = "Order", legend_gp = gpar(col = "black"), background = order_col, pch = order_pch, type = "points"
  )
  lgd_type <- ComplexHeatmap::Legend(
    labels = gt_render(names(type_col)),
    title = "Virus type", legend_gp = gpar(fill = type_col)
  )


  pd <- ComplexHeatmap::packLegend(lgd_type, lgd_king, lgd_phylum, lgd_class, lgd_order, max_height = unit(10, "cm"))

  ComplexHeatmap::draw(hm,
    annotation_legend_list = pd,
    annotation_legend_side = legend_side, heatmap_legend_side = legend_side,
    merge_legend = T, legend_grouping = "original"
  )
}