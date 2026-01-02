rm(list = ls())
graphics.off()

DATA_FILE <- "/パスを挿入/analysis_code/gut_microbiome_16S_mock_gravity_culture.csv"
OUTPUT_DIR <- "/パスを挿入/analysis_code/figures"
#これで/analysis_code/figuresにfigureが保存される

pkgs <- c("vegan", "ggplot2", "dplyr", "tidyr", "RColorBrewer", 
          "pheatmap", "gridExtra", "scales")

for (p in pkgs) {
  if (!require(p, character.only = TRUE, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org/")
    library(p, character.only = TRUE)
  }
}

gravity_colors <- c(
  "0g" = "#E64B35", "1_6g" = "#4DBBD5", "1g" = "#00A087",
  "1g_s" = "#3C5488", "5g" = "#F39B7F"
)

gravity_labels <- c(
  "0g" = "Microgravity (0g)", "1_6g" = "Lunar (1/6g)", "1g" = "Earth (1g)",
  "1g_s" = "Static Control", "5g" = "Hypergravity (5g)"
)

raw_data <- read.csv(DATA_FILE, stringsAsFactors = FALSE)

meta_cols <- c("SampleID", "SampleType", "Donor", "Gravity", "Time", "Replicate", "TotalReads")
taxa_cols <- setdiff(colnames(raw_data), meta_cols)

metadata <- raw_data[, meta_cols]
abundance <- raw_data[, taxa_cols]

rel_abundance <- abundance / rowSums(abundance)
rel_abundance[is.na(rel_abundance)] <- 0

culture_idx <- metadata$Gravity != "baseline"
culture_rel <- rel_abundance[culture_idx, ]
culture_meta <- metadata[culture_idx, ]

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Figure 1
dist_matrix <- vegdist(culture_rel, method = "bray")
pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)

eig_pos <- pcoa_result$eig[pcoa_result$eig > 0]
var_explained <- eig_pos / sum(eig_pos) * 100

pcoa_df <- data.frame(
  PC1 = pcoa_result$points[, 1],
  PC2 = pcoa_result$points[, 2],
  Gravity = culture_meta$Gravity,
  Time = culture_meta$Time,
  Donor = as.factor(culture_meta$Donor)
)

p1 <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Gravity, shape = Time)) +
  stat_ellipse(aes(fill = Gravity), geom = "polygon", alpha = 0.15, 
               level = 0.95, show.legend = FALSE) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(values = gravity_colors, labels = gravity_labels) +
  scale_fill_manual(values = gravity_colors) +
  scale_shape_manual(values = c("8h" = 16, "16h" = 17, "24h" = 15)) +
  labs(x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2]),
       title = "Figure 1: Principal Coordinates Analysis (PCoA)",
       subtitle = "Bray-Curtis distance with 95% confidence ellipses") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "gray40"))

ggsave(file.path(OUTPUT_DIR, "Figure1_PCoA.png"), p1, width = 11, height = 8, dpi = 300)

# Figure 2: Stacked Barplot
mean_abundance <- colMeans(culture_rel)
top15_taxa <- names(sort(mean_abundance, decreasing = TRUE))[1:15]

plot_data2 <- culture_rel %>%
  as.data.frame() %>%
  mutate(Gravity = culture_meta$Gravity, Time = culture_meta$Time) %>%
  group_by(Gravity, Time) %>%
  summarise(across(all_of(taxa_cols), mean), .groups = "drop") %>%
  pivot_longer(cols = all_of(taxa_cols), names_to = "Taxon", values_to = "Abundance") %>%
  mutate(Taxon = ifelse(Taxon %in% top15_taxa, Taxon, "Others")) %>%
  group_by(Gravity, Time, Taxon) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

plot_data2$Taxon <- factor(plot_data2$Taxon, levels = c(top15_taxa, "Others"))
plot_data2$Time <- factor(plot_data2$Time, levels = c("8h", "16h", "24h"))

taxa_colors <- c(brewer.pal(12, "Set3"), brewer.pal(4, "Pastel1")[1:3], "gray70")
names(taxa_colors) <- c(top15_taxa, "Others")

p2 <- ggplot(plot_data2, aes(x = Gravity, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack", width = 0.85) +
  facet_wrap(~ Time, ncol = 3) +
  scale_fill_manual(values = taxa_colors) +
  scale_x_discrete(labels = gravity_labels) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Gravity Condition", y = "Relative Abundance",
       title = "Figure 2: Microbial Community Composition",
       subtitle = "Top 15 taxa shown, others grouped") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.text = element_text(size = 7),
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(OUTPUT_DIR, "Figure2_StackedBar.png"), p2, width = 15, height = 9, dpi = 300)

# Figure 3
pcoa_df$Time_num <- as.numeric(gsub("h", "", pcoa_df$Time))
pcoa_df$Group <- paste(pcoa_df$Donor, pcoa_df$Gravity, sep = "_")

trajectory_data <- pcoa_df %>%
  arrange(Group, Time_num) %>%
  group_by(Group, Gravity, Donor) %>%
  mutate(PC1_end = lead(PC1), PC2_end = lead(PC2)) %>%
  filter(!is.na(PC1_end))

p3 <- ggplot() +
  geom_segment(data = trajectory_data,
               aes(x = PC1, y = PC2, xend = PC1_end, yend = PC2_end, color = Gravity),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               linewidth = 1, alpha = 0.7) +
  geom_point(data = pcoa_df,
             aes(x = PC1, y = PC2, color = Gravity, shape = Time), size = 3.5) +
  scale_color_manual(values = gravity_colors, labels = gravity_labels) +
  scale_shape_manual(values = c("8h" = 16, "16h" = 17, "24h" = 15)) +
  facet_wrap(~ Donor, ncol = 3) +
  labs(x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2]),
       title = "Figure 3: Temporal Trajectories in PCoA Space",
       subtitle = "Arrows indicate time progression (8h to 16h to 24h)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold", size = 12))

ggsave(file.path(OUTPUT_DIR, "Figure3_Trajectory.png"), p3, width = 14, height = 6, dpi = 300)

# Figure 4
betadisper_result <- betadisper(dist_matrix, culture_meta$Gravity)
permutest_result <- permutest(betadisper_result, permutations = 999)
p_value <- permutest_result$tab$`Pr(>F)`[1]
f_stat <- permutest_result$tab$F[1]

distances_df <- data.frame(
  Distance = betadisper_result$distances,
  Gravity = culture_meta$Gravity
)

p4 <- ggplot(distances_df, aes(x = Gravity, y = Distance, fill = Gravity)) +
  geom_violin(alpha = 0.4, show.legend = FALSE) +
  geom_boxplot(width = 0.15, alpha = 0.9, outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(width = 0.08, alpha = 0.6, size = 2.5, show.legend = FALSE) +
  scale_fill_manual(values = gravity_colors) +
  scale_x_discrete(labels = gravity_labels) +
  annotate("text", x = Inf, y = Inf, 
           label = sprintf("F = %.2f, p = %.4f", f_stat, p_value),
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
  labs(x = "Gravity Condition", y = "Distance to Centroid",
       title = "Figure 4: Beta Dispersion Analysis (BETADISPER)",
       subtitle = "Distance from samples to group centroid in PCoA space") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "Figure4_Betadisper.png"), p4, width = 10, height = 8, dpi = 300)

# Figure 5
top10_taxa <- names(sort(mean_abundance, decreasing = TRUE))[1:10]

taxa_stats <- culture_rel[, top10_taxa] %>%
  as.data.frame() %>%
  mutate(Gravity = culture_meta$Gravity) %>%
  pivot_longer(cols = all_of(top10_taxa), names_to = "Taxon", values_to = "Abundance") %>%
  group_by(Gravity, Taxon) %>%
  summarise(Mean = mean(Abundance), SE = sd(Abundance) / sqrt(n()), .groups = "drop")

taxa_stats$Taxon <- factor(taxa_stats$Taxon, levels = top10_taxa)

p5 <- ggplot(taxa_stats, aes(x = Taxon, y = Mean, fill = Gravity)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.85), width = 0.8) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                position = position_dodge(width = 0.85), width = 0.3) +
  scale_fill_manual(values = gravity_colors, labels = gravity_labels) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Taxon", y = "Mean Relative Abundance (± SE)",
       title = "Figure 5: Dominant Taxa Abundance by Gravity Condition",
       subtitle = "Top 10 taxa shown") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 50, hjust = 1, size = 9),
        legend.position = "top")

ggsave(file.path(OUTPUT_DIR, "Figure5_TaxaComparison.png"), p5, width = 14, height = 9, dpi = 300)

# Figure 6
top30_taxa <- names(sort(mean_abundance, decreasing = TRUE))[1:30]

heatmap_data <- as.matrix(culture_rel[, top30_taxa])
heatmap_zscore <- scale(heatmap_data)

rownames(heatmap_zscore) <- paste(
  culture_meta$Gravity, culture_meta$Time, 
  culture_meta$Donor, culture_meta$Replicate, sep = "_"
)

annotation_col <- data.frame(
  Gravity = culture_meta$Gravity,
  Time = culture_meta$Time
)
rownames(annotation_col) <- rownames(heatmap_zscore)

annotation_colors <- list(
  Gravity = gravity_colors,
  Time = c("8h" = "#FEE08B", "16h" = "#FDAE61", "24h" = "#F46D43")
)

png(file.path(OUTPUT_DIR, "Figure6_Heatmap.png"), 
    width = 14, height = 12, units = "in", res = 300)

pheatmap(
  t(heatmap_zscore),
  color = colorRampPalette(c("#313695", "#4575B4", "#74ADD1", "#E0E0E0", 
                             "#FDAE61", "#F46D43", "#A50026"))(100),
  breaks = seq(-3, 3, length.out = 101),
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  clustering_method = "complete",
  show_colnames = FALSE,
  fontsize_row = 8,
  main = "Figure 6: Microbial Community Heatmap (Z-score normalized)",
  border_color = NA
)

dev.off()

# Figure 7
permanova_result <- adonis2(dist_matrix ~ Gravity + Time, data = culture_meta, permutations = 999)

simper_result <- simper(culture_rel, culture_meta$Gravity)
simper_summary <- summary(simper_result)

permanova_df <- data.frame(
  Factor = c("Gravity", "Time", "Residual"),
  R2 = c(permanova_result$R2[1], permanova_result$R2[2], permanova_result$R2[3]),
  p_value = c(permanova_result$`Pr(>F)`[1], permanova_result$`Pr(>F)`[2], NA)
)
permanova_df$Significance <- ifelse(is.na(permanova_df$p_value), "",
                                    ifelse(permanova_df$p_value < 0.001, "***",
                                           ifelse(permanova_df$p_value < 0.01, "**",
                                                  ifelse(permanova_df$p_value < 0.05, "*", "ns"))))
permanova_df$Factor <- factor(permanova_df$Factor, levels = c("Gravity", "Time", "Residual"))

p7a <- ggplot(permanova_df, aes(x = Factor, y = R2 * 100, fill = Factor)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Significance), vjust = -0.5, size = 6) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "", y = "Variance Explained (%)",
       title = "A) PERMANOVA Results",
       subtitle = "* p<0.05, ** p<0.01, *** p<0.001") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  ylim(0, max(permanova_df$R2 * 100, na.rm = TRUE) * 1.3)

first_comp <- names(simper_summary)[1]
simper_data <- simper_summary[[first_comp]]

simper_df <- data.frame(
  Taxon = rownames(simper_data)[1:10],
  Contribution = simper_data$average[1:10] * 100,
  Cumulative = simper_data$cumsum[1:10] * 100
)
simper_df$Taxon <- factor(simper_df$Taxon, levels = rev(simper_df$Taxon))

p7b <- ggplot(simper_df, aes(x = Taxon, y = Contribution)) +
  geom_bar(stat = "identity", fill = "#3C5488", width = 0.7) +
  geom_line(aes(y = Cumulative, group = 1), color = "#E64B35", linewidth = 1.2) +
  geom_point(aes(y = Cumulative), color = "#E64B35", size = 3) +
  coord_flip() +
  scale_y_continuous(name = "Individual Contribution (%)",
                     sec.axis = sec_axis(~ ., name = "Cumulative (%)")) +
  labs(x = "", title = sprintf("B) SIMPER Analysis (%s)", first_comp),
       subtitle = "Taxa contribution to between-group dissimilarity") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 9))

combined_plot <- grid.arrange(p7a, p7b, ncol = 2, widths = c(1, 1.3))
ggsave(file.path(OUTPUT_DIR, "Figure7_Statistics.png"), combined_plot, 
       width = 15, height = 7, dpi = 300)

cat("\n完了: ", OUTPUT_DIR, "\n")
for (f in list.files(OUTPUT_DIR, pattern = "\\.png$")) cat(" ", f, "\n")

if (Sys.info()["sysname"] == "Darwin") system(paste("open", shQuote(OUTPUT_DIR)))
