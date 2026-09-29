# TIMES sex differences plasma proteomics analysis scripts
#
# Expected project structure:
#   data/raw/
#   data/source/
#   data/processed/
#   figures/
#

library(tidyverse)
library(arrow)
library(diann)
library(effsize)
library(FSA)
library(ggh4x)
library(ggplot2)
library(ggpmisc)
library(ggstatsplot)
library(ggtext)
library(glue)
library(patchwork)
library(rstatix)
library(svglite)

data_dir <- "data"
raw_dir <- file.path(data_dir, "raw")
source_dir <- file.path(data_dir, "source")
processed_dir <- file.path(data_dir, "processed")
figure_dir <- "figures"

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# Process DIA-NN output and generate protein-level LFQ and concentration tables

report <- read_parquet(
  file.path(raw_dir, "diann_report.parquet")
)

reference_proteins <- read.csv(
  file.path(source_dir, "reference_proteins.csv")
)

sampling_contaminants <- read.csv(
  file.path(source_dir, "sampling_contaminants.csv")
)

sample_metadata <- read.csv(
  file.path(source_dir, "sample_metadata.csv")
)

contaminant_genes <- sampling_contaminants$Genes

report <- report %>%
  mutate(
    File.Name = sample_metadata$File.Name[
      match(Run, sample_metadata$Raw.Name)
    ],
    Group = sample_metadata$Group[
      match(Run, sample_metadata$Raw.Name)
    ]
  ) %>%
  filter(!is.na(Group), Group != "") %>%
  group_by(File.Name) %>%
  mutate(Precursor.Count = n_distinct(Precursor.Id)) %>%
  ungroup()

threshold <- report %>%
  filter(Group == "sample") %>%
  summarise(threshold = median(Precursor.Count, na.rm = TRUE) * 0.6) %>%
  pull(threshold)

report <- report %>%
  filter(Precursor.Count >= threshold) %>%
  filter(
    Q.Value <= 0.01,
    Lib.Q.Value <= 0.01,
    Lib.PG.Q.Value <= 0.01
  ) %>%
  group_by(Genes) %>%
  mutate(Total.Peptide.Count = n_distinct(Stripped.Sequence)) %>%
  ungroup() %>%
  filter(Total.Peptide.Count >= 2) %>%
  filter(!grepl("\\(UniMod:35\\)", Precursor.Id))

lfq <- diann_maxlfq(
  report,
  group.header = "Genes",
  id.header = "Precursor.Id",
  quantity.header = "Precursor.Normalised"
) %>%
  as.data.frame() %>%
  rownames_to_column("Genes")

lfq_long <- lfq %>%
  pivot_longer(
    cols = -Genes,
    names_to = "File.Name",
    values_to = "LFQ"
  ) %>%
  filter(!is.na(LFQ)) %>%
  left_join(
    sample_metadata %>%
      dplyr::select(
        File.Name,
        Subject = Prob.Nr.,
        Time.point,
        Cohort,
        Plate,
        Sex,
        BMI,
        Group
      ),
    by = "File.Name"
  ) %>%
  mutate(logI = log10(LFQ))

amsbio <- lfq_long %>%
  filter(Group == "AMSBIO") %>%
  left_join(
    reference_proteins %>% dplyr::select(Genes, LogC),
    by = "Genes"
  ) %>%
  filter(!is.na(LogC))

fit <- lm(logI ~ LogC, data = amsbio)

slope <- coef(fit)[["LogC"]]
intercept <- coef(fit)[["(Intercept)"]]

lfq_long <- lfq_long %>%
  mutate(
    Concentration = 10^((logI - intercept) / slope)
  )

lfq_long <- lfq_long %>%
  filter(
    !grepl("IGLV|IGKV|IGHV", Genes)
  )

lfq_no_contaminants <- lfq_long %>%
  separate_rows(Genes, sep = ";") %>%
  filter(!Genes %in% contaminant_genes) %>%
  group_by(File.Name, Subject, Time.point, Cohort, Plate,BMI, Group, Sex, LFQ, logI,
           Concentration) %>%
  summarise(
    Genes = paste(unique(Genes), collapse = ";"),
    .groups = "drop"
  )


lfq_wide <- lfq_no_contaminants %>%
  dplyr::select(File.Name, Subject,Time.point, Cohort, Plate, Group, Sex,BMI, Genes, LFQ) %>%
  pivot_wider(names_from = Genes, values_from = LFQ)

concentration_wide <- lfq_no_contaminants %>%
  dplyr::select(File.Name, Subject,Time.point, Cohort, Plate, Group, Sex,BMI, Genes, Concentration) %>%
  pivot_wider(names_from = Genes, values_from = Concentration)

write.csv(
  filter(lfq_wide, Cohort == "TIMES"),
  file.path(processed_dir, "LFQ_TIMES.csv"),
  row.names = FALSE
)

write.csv(
  filter(concentration_wide, Cohort == "TIMES"),
  file.path(processed_dir, "Conc_TIMES.csv"),
  row.names = FALSE
)

write.csv(
  filter(lfq_wide, Cohort == "AICOVI"),
  file.path(processed_dir, "LFQ_AICOVI.csv"),
  row.names = FALSE
)

write.csv(
  filter(concentration_wide, Cohort == "AICOVI"),
  file.path(processed_dir, "Conc_AICOVI.csv"),
  row.names = FALSE
)

write.csv(
  filter(lfq_wide, Cohort == "AMSBIO"),
  file.path(processed_dir, "LFQ_AMSBIO.csv"),
  row.names = FALSE
)



# Compare protein concentrations between female and male participants

full_report <- read.csv(
  file.path(processed_dir, "Conc_TIMES.csv")
) %>%
  dplyr::select(-1) %>%
  filter(!grepl("P09_TrP12_27", File.Name, ignore.case = TRUE)) %>%
  pivot_longer(
    cols = -c(1:7),
    names_to = "Genes",
    values_to = "nM"
  ) %>%
  filter(!is.na(nM))

selected_genes <- c(
  "PZP", "SHBG", "FETUB", "AGT", "SERPINA6",
  "SERPINA7", "CP", "APOL1", "KNG1"
)

fit_gene_kruskal <- function(data) {
  gene <- unique(data$Genes)

  dunnTest(
    nM ~ Sex,
    data = data,
    method = "none"
  )$res %>%
    filter(Comparison %in% c("female - male", "male - female")) %>%
    mutate(Genes = gene)
}

stats_all <- full_report %>%
  group_split(Genes) %>%
  map_dfr(fit_gene_kruskal) %>%
  mutate(padj = p.adjust(P.unadj, method = "BH"))

fc_table <- full_report %>%
  group_by(Genes, Sex) %>%
  summarise(
    mean_nM = mean(nM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Sex, values_from = mean_nM) %>%
  mutate(fold_change = female / male)

calculate_effect_size <- function(data) {
  female <- data %>% filter(Sex == "female") %>% pull(nM)
  male <- data %>% filter(Sex == "male") %>% pull(nM)

  result <- effsize::cliff.delta(
    female,
    male,
    conf.level = 0.95
  )

  tibble(
    effect_size = unname(result$estimate),
    ci_low = result$conf.int[1],
    ci_high = result$conf.int[2]
  )
}

make_sex_plot <- function(gene_name) {
  data <- full_report %>% filter(Genes == gene_name)

  n_table <- data %>%
    group_by(Sex) %>%
    summarise(n = sum(!is.na(nM)), .groups = "drop")

  effect <- calculate_effect_size(data)

  padj <- stats_all %>%
    filter(Genes == gene_name) %>%
    pull(padj)

  fold_change <- fc_table %>%
    filter(Genes == gene_name) %>%
    pull(fold_change)

  y_max <- max(data$nM, na.rm = TRUE)

  title <- paste0(
    "FC = ", round(fold_change, 2),
    ", δ = ", round(effect$effect_size, 2),
    ", 95% CI = [",
    round(effect$ci_low, 2), ", ",
    round(effect$ci_high, 2), "]"
  )

  ggstatsplot::ggbetweenstats(
    data = data,
    x = Sex,
    y = nM,
    type = "nonparametric",
    pairwise.display = "none",
    bf.message = FALSE,
    results.subtitle = FALSE,
    xlab = "",
    ylab = "concentration [nM]",
    ggtheme = theme_light(),
    violin.args = list(width = 0.6, alpha = 0.2, linewidth = 0.2),
    boxplot.args = list(width = 0.3, alpha = 0.04, linewidth = 0.15),
    point.args = list(
      position = position_jitterdodge(dodge.width = 0.5),
      alpha = 0.4,
      size = 0.8,
      stroke = 0.5
    ),
    centrality.label.args = list(size = 1.8, nudge_y = y_max*0.25, direction = "y", segment.linetype = 4,
                                 min.segment.length = 0,label.size = 1.8,                      
                                 label.padding = unit(0.1, "lines")),
    centrality.point.args = list(
      size  = 1,
      colour = "darkred"
    )
  ) +
    scale_fill_manual(
      values = c(female = "#e5b822", male = "#04afbb")
    ) +
    scale_color_manual(
      values = c(female = "#e5b822", male = "#04afbb")
    ) +
    scale_y_continuous(
      labels = scales::scientific,
      expand = expansion(mult = c(0.05, 0.09))
    ) +
    scale_x_discrete(
      labels = c(
        female = paste0("n = ", n_table$n[n_table$Sex == "female"]),
        male = paste0("n = ", n_table$n[n_table$Sex == "male"])
      )
    ) +
    annotate(
      "text",
      x = 1.5,
      y = y_max * 1.12,
      label = paste0(
        "p[FDR-adj] == ",
        formatC(padj, format = "e", digits = 2)
      ),
      parse = TRUE,
      size = 1.9
    ) +
    labs(title = title) +
    theme(
      axis.text.x = element_text(size = 6, margin = margin(t = 4)),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 5.5),
      axis.title.y = element_text(size = 5.5),
      legend.position = "none",
      plot.title = element_text(size = 5.5, hjust = 0)
    )
}

sex_plot_dir <- file.path(figure_dir, "sex_comparisons")
dir.create(sex_plot_dir, recursive = TRUE, showWarnings = FALSE)

walk(selected_genes, function(gene) {
  ggsave(
    file.path(sex_plot_dir, paste0("sex_comparison_", gene, ".svg")),
    make_sex_plot(gene),
    width = 5.5,
    height = 4.5,
    units = "cm"
  )
})


# Compare protein concentrations between female and male participants in the AICOVI cohort

full_report <- read.csv(
  file.path(processed_dir, "Conc_AICOVI.csv")
) %>%
  dplyr::select(-1) %>%
  pivot_longer(
    cols = -c(1:7),
    names_to = "Genes",
    values_to = "nM"
  ) %>%
  filter(!is.na(nM))

subject_map <- tibble(
  Subject_old = c("P01", "P03", "P04", "P05", "P06", "P09", "P12", "P13", "P19", "P22", "P26"),
  Subject_new = sprintf("P%02d", seq_along(Subject_old))
)

full_report <- full_report %>%
  left_join(subject_map, by = c("Subject" = "Subject_old")) %>%
  mutate(Subject = coalesce(Subject_new, Subject)) %>%
  dplyr::select(-Subject_new)

stats_all <- full_report %>%
  group_split(Genes) %>%
  map_dfr(fit_gene_kruskal) %>%
  mutate(padj = p.adjust(P.unadj, method = "BH"))

fc_table <- full_report %>%
  group_by(Genes, Sex) %>%
  summarise(mean_nM = mean(nM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Sex, values_from = mean_nM) %>%
  mutate(fold_change = female / male)

AICOVI_cohort_dir <- file.path(figure_dir, "AICOVI_cohort_sex_comparisons")
dir.create(AICOVI_cohort_dir, recursive = TRUE, showWarnings = FALSE)

walk(selected_genes, function(gene) {
  ggsave(
    file.path(AICOVI_cohort_dir, paste0("sex_comparison_", gene, ".svg")),
    make_sex_plot(gene),
    width = 5.5,
    height = 4.5,
    units = "cm"
  )
})


# Summarise longitudinal protein concentrations by participant

full_report <- read.csv(
  file.path(processed_dir, "Conc_TIME.csv")
) %>%
  dplyr::select(-1) %>%
  filter(!grepl("P09_TrP12_27", File.Name, ignore.case = TRUE)) %>%
  pivot_longer(
    cols = -c(1:7),
    names_to = "Genes",
    values_to = "nM"
  ) %>%
  filter(!is.na(nM))

selected_genes <- c(
  "PZP", "SHBG", "FETUB", "AGT",
  "SERPINA6", "SERPINA7", "CP", "APOL1", "KNG1"
)

gene_data <- full_report %>%
  mutate(Genes = factor(Genes, levels = selected_genes)) %>%
  filter(Genes %in% selected_genes) %>%
  group_by(Genes, Subject, Sex) %>%
  summarise(
    Mean_Cexp = mean(nM, na.rm = TRUE),
    SD_Cexp = sd(nM, na.rm = TRUE),
    n_timepoints = sum(!is.na(nM)),
    .groups = "drop"
  )

gene_data <- gene_data %>%
  group_by(Sex, Genes) %>%
  mutate(Subject = as.numeric(as.character(Subject))) %>%
  arrange(Subject) %>%
  mutate(Subject = factor(Subject, levels = unique(Subject))) %>%
  ungroup()

avg_data <- gene_data %>%
  group_by(Genes, Sex) %>%
  summarise(
    avg = mean(Mean_Cexp, na.rm = TRUE),
    sd = sd(Mean_Cexp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cv = 100 * sd / avg,
    label = paste0(
      "mean: ",
      formatC(avg, format = "f", big.mark = ",", digits = 0),
      " ± ",
      formatC(sd, format = "f", big.mark = ",", digits = 0),
      " nM; CV: ",
      round(cv, 1),
      "%"
    )
  )

n_labels <- gene_data %>%
  group_by(Genes) %>%
  mutate(
    y_label = -0.2 * max(Mean_Cexp, na.rm = TRUE),
    n_label = paste0("n=", n_timepoints)
  ) %>%
  ungroup()

longitudinal_barplot <- ggplot(
  gene_data,
  aes(x = Subject, y = Mean_Cexp, fill = Sex)
) +
  geom_col(color = "black", linewidth = 0.2) +
  geom_errorbar(
    aes(
      ymin = Mean_Cexp - SD_Cexp,
      ymax = Mean_Cexp + SD_Cexp
    ),
    linewidth = 0.2,
    width = 0.2
  ) +
  geom_hline(
    data = avg_data,
    aes(yintercept = avg),
    linetype = "dashed",
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  geom_label(
    data = avg_data,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.01,
    vjust = 1.3,
    size = 2.25,
    label.size = 0.2,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = n_labels,
    aes(x = Subject, y = y_label, label = n_label),
    angle = 90,
    size = 1.8,
    color = "grey50",
    inherit.aes = FALSE
  ) +
  ggh4x::facet_grid2(
    Genes ~ Sex,
    scales = "free",
    space = "free_x"
  ) +
  ggh4x::force_panelsizes(cols = c(1.9, 1.1)) +
  scale_fill_manual(
    values = c(male = "#04afbb", female = "#e5b822")
  ) +
  scale_y_continuous(
    labels = scales::scientific,
    expand = expansion(mult = c(0.15, 0.05))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "participant",
    y = "concentration [nM]"
  ) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 9),
    strip.text = element_text(size = 7, face = "bold"),
    panel.spacing = unit(0.2, "lines"),
    plot.margin = margin(10, 10, 45, 10)
  )

ggsave(
  file.path(figure_dir, "longitudinal_summary.svg"),
  longitudinal_barplot,
  width = 21,
  height = 28,
  units = "cm"
)


#Compare longitudinal protein concentrations in the AICOVI cohort

full_report <- read.csv(
  file.path(processed_dir, "Conc_AICOVI.csv")
) %>%
  dplyr::select(-1) %>%
  pivot_longer(
    cols = -c(1:7),
    names_to = "Genes",
    values_to = "nM"
  ) %>%
  filter(!is.na(nM))

full_report <- full_report %>%
  left_join(subject_map, by = c("Subject" = "Subject_old")) %>%
  mutate(Subject = coalesce(Subject_new, Subject)) %>%
  dplyr::select(-Subject_new)

gene_data <- full_report %>%
  filter(Genes %in% selected_genes) %>%
  mutate(Genes = factor(Genes, levels = selected_genes)) %>%
  group_by(Genes, Subject, Sex) %>%
  summarise(
    Mean_Cexp = mean(nM, na.rm = TRUE),
    SD_Cexp = sd(nM, na.rm = TRUE),
    n_timepoints = sum(!is.na(nM)),
    .groups = "drop"
  )

avg_data <- gene_data %>%
  group_by(Genes, Sex) %>%
  summarise(
    avg = mean(Mean_Cexp, na.rm = TRUE),
    sd = sd(Mean_Cexp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cv = 100 * sd / avg,
    label = paste0(
      "mean: ",
      formatC(avg, format = "f", big.mark = ",", digits = 0),
      " ± ",
      formatC(sd, format = "f", big.mark = ",", digits = 0),
      " nM; CV: ",
      round(cv, 1),
      "%"
    )
  )

n_labels <- gene_data %>%
  group_by(Genes) %>%
  mutate(
    y_label = -0.2 * max(Mean_Cexp, na.rm = TRUE),
    n_label = paste0("n=", n_timepoints)
  ) %>%
  ungroup()

AICOVI_longitudinal_barplot <- ggplot(
  gene_data,
  aes(x = Subject, y = Mean_Cexp, fill = Sex)
) +
  geom_col(color = "black", linewidth = 0.2) +
  geom_errorbar(
    aes(ymin = Mean_Cexp - SD_Cexp, ymax = Mean_Cexp + SD_Cexp),
    linewidth = 0.2,
    width = 0.2
  ) +
  geom_hline(
    data = avg_data,
    aes(yintercept = avg),
    linetype = "dashed",
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  geom_label(
    data = avg_data,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.01,
    vjust = 1.3,
    size = 2.25,
    label.size = 0.2,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = n_labels,
    aes(x = Subject, y = y_label, label = n_label),
    angle = 90,
    size = 1.8,
    color = "grey50",
    inherit.aes = FALSE
  ) +
  ggh4x::facet_grid2(Genes ~ Sex, scales = "free", space = "free_x") +
  ggh4x::force_panelsizes(cols = c(1.9, 1.1)) +
  scale_fill_manual(values = c(male = "#04afbb", female = "#e5b822")) +
  scale_y_continuous(
    labels = scales::scientific,
    expand = expansion(mult = c(0.16, 0.16))
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "participant", y = "concentration [nM]") +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 9),
    strip.text = element_text(size = 7, face = "bold"),
    panel.spacing = unit(0.2, "lines")
  )

ggsave(
  file.path(figure_dir, "AICOVI_cohort_longitudinal_summary.svg"),
  AICOVI_longitudinal_barplot,
  width = 18,
  height = 27,
  units = "cm"
)


# Plot selected longitudinal donor profiles

full_report <- read.csv(
  file.path(processed_dir, "Conc_TIMES.csv")
) %>%
  dplyr::select(-1) %>%
  filter(!grepl("P09_TrP12_27", File.Name, ignore.case = TRUE)) %>%
  pivot_longer(
    cols = -c(1:7),
    names_to = "Genes",
    values_to = "nM"
  ) %>%
  filter(!is.na(nM))

month_order <- c(
  "Jun", "Jul", "Aug", "Sept", "Oct", "Nov",
  "Dec", "Jan", "Feb", "Mar", "Apr", "May"
)

month_colors <- c(
  Jun = "#a04e4e",
  Jul = "#d25554",
  Aug = "#f37f6c",
  Sept = "#f99e82",
  Oct = "#febf9d",
  Nov = "#f9e6b6",
  Dec = "#e9e5cb",
  Jan = "#d8dac0",
  Feb = "#adbc91",
  Mar = "#95a971",
  Apr = "#70865a",
  May = "#527350"
)

full_report$Month <- month_order[full_report$Time.Point]

selected_subjects <- c(34, 35, 54, 7, 29, 3, 30, 32, 47)
bar.colors <- c("#e54c35ff","#4dbbd5FF","#00a087FF","#3c5488FF","#f39b7fFF","#8491b4FF","#91d1c2FF","#da1f26FF","#7f6148FF")


profile_data <- full_report %>%
  filter(
    Genes %in% selected_genes,
    Subject %in% selected_subjects
  ) %>%
  dplyr::select(Subject, Month, Sex, Genes, nM) %>%
  mutate(
    Subject = factor(Subject, levels = selected_subjects),
    Month = factor(Month, levels = month_order),
    Genes = factor(Genes, levels = selected_genes),
    Sex = factor(Sex, levels = c("female", "male"))
  )

selected_profile_plot <- ggplot(
  profile_data,
  aes(
    x = Month,
    y = nM,
    fill = Genes
  )
) +
  geom_col(width = 0.7, linewidth = 0.2, color = "black") +
  facet_grid2(
    Genes ~ Subject,
    scales = "free_y"
  ) +
  scale_y_continuous(labels = scales::scientific) +
  scale_fill_manual(values =bar.colors) +
  labs(
    x = "month",
    y = "concentration [nM]"
  ) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, size = 6, hjust = 1),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 8),
    strip.text = element_text(size = 8, face = "bold"),
    panel.spacing = unit(0.15, "lines")
  )

ggsave(
  file.path(figure_dir, "selected_donor_profiles.svg"),
  selected_profile_plot,
  width = 22,
  height = 16,
  units = "cm"
)

all_profile_data <- full_report %>%
  filter(Genes %in% selected_genes) %>%
  dplyr::select(Subject, Month, Sex, Genes, nM) %>%
  mutate(
    Subject = factor(Subject, levels = unique(sort(Subject))),
    Month = factor(Month, levels = month_order),
    Genes = factor(Genes, levels = selected_genes),
    Sex = factor(Sex, levels = c("female", "male"))
  )

all_profile_plot <- ggplot(
  all_profile_data,
  aes(x = Month, y = nM, fill = Genes)
) +
  geom_col(width = 0.6, linewidth = 0.2) +
  facet_grid2(
    Genes ~ Subject,
    scales = "free_y"
  ) +
  scale_y_continuous(labels = scales::scientific) +
  labs(
    x = "month",
    y = "concentration [nM]"
  ) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, size = 4, hjust = 1),
    axis.text.y = element_text(size = 6.5),
    axis.title = element_text(size = 11),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  file.path(figure_dir, "all_participant_profiles.pdf"),
  all_profile_plot,
  width = 34,
  height = 5,
  units = "in"
)


# Assess sex-specific correlations with AGT

full_report <- read.csv(
  file.path(processed_dir, "Conc_TIMES.csv")
) %>%
  dplyr::select(-1) %>%
  filter(!grepl("P09_TrP12_27", File.Name, ignore.case = TRUE)) %>%
  pivot_longer(
    cols = -c(1:7),
    names_to = "Genes",
    values_to = "nM"
  ) %>%
  filter(!is.na(nM))

selected_correlations <- c(
  "PZP", "SHBG", "FETUB", "SERPINA6",
  "SERPINA7", "CP", "APOL1", "KNG1"
)

correlation_data <- full_report %>%
  filter(Genes %in% c("AGT", selected_correlations)) %>%
  dplyr::select(Subject, Time.Point, Sex, Genes, nM) %>%
  pivot_wider(names_from = Genes, values_from = nM) %>%
  drop_na(AGT) %>%
  pivot_longer(
    cols = all_of(selected_correlations),
    names_to = "Genes",
    values_to = "Protein_value"
  ) %>%
  drop_na(Protein_value) %>%
  mutate(
    Genes = factor(Genes, levels = selected_correlations),
    Sex = factor(Sex, levels = c("female", "male"))
  )

correlation_stats <- correlation_data %>%
  group_by(Sex, Genes) %>%
  summarise(
    spearman = cor(
      Protein_value,
      AGT,
      method = "spearman",
      use = "complete.obs"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    label = glue(
      "r<sub>s</sub> = {format(round(spearman, 2), nsmall = 2)}"
    )
  )

correlation_plot <- ggplot(
  correlation_data,
  aes(x = Protein_value, y = AGT)
) +
  geom_point(alpha = 0.5, size = 1,colour="gray9") +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    linewidth = 0.4,
    colour = "hotpink3",
    alpha = 0.6
  ) +
  ggh4x::facet_grid2(
    rows = vars(Sex),
    cols = vars(Genes),
    scales = "free_x"
  ) +
  geom_richtext(
    data = filter(correlation_stats, Sex == "female"),
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1,
    vjust = 1.2,
    inherit.aes = FALSE,
    size = 2
  ) +
  labs(
    x = "concentration [nM]",
    y = "AGT concentration [nM]"
  ) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 7),
    axis.text.x = element_text(angle = 90, size = 6),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 7, face = "bold"),
    panel.spacing = unit(0.2, "lines")
  )

calc_female_in_male_ellipse <- function(data, level = 0.95) {
  male <- data %>% filter(Sex == "male")
  female <- data %>% filter(Sex == "female")
  
  if (nrow(male) < 3 || nrow(female) == 0) {
    return(tibble(
      pct_female_inside = NA_real_,
      n_female = nrow(female)
    ))
  }
  
  values <- male[, c("Protein_value", "AGT")]
  center <- colMeans(values)
  covariance <- cov(values)
  
  distance <- mahalanobis(
    female[, c("Protein_value", "AGT")],
    center = center,
    cov = covariance
  )
  
  tibble(
    pct_female_inside =
      mean(distance <= qchisq(level, df = 2)) * 100,
    n_female = nrow(female)
  )
}

ellipse_data <- correlation_data %>%
  filter(Sex == "male") %>%
  mutate(Sex = "female")

ellipse_results <- correlation_data %>%
  group_by(Genes) %>%
  group_modify(~ calc_female_in_male_ellipse(.x)) %>%
  ungroup() %>%
  mutate(
    label = sprintf("%.1f%%", pct_female_inside),
    Sex = "female"
  )

correlation_ellipse_plot <- correlation_plot +
  stat_ellipse(
    data = ellipse_data,
    aes(x = Protein_value, y = AGT),
    geom = "polygon",
    type = "norm",
    level = 0.95,
    linewidth = 0.3,
    fill = "#04afbb",
    colour = "#04afbb",
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  geom_richtext(
    data = ellipse_results,
    aes(x = Inf, y = -Inf, label = label),
    hjust = 1.4,
    vjust = -0.4,
    inherit.aes = FALSE,
    size = 2,
    label.colour = "#04afbb"
  )

ggsave(
  file.path(figure_dir, "AGT_correlations_male_ellipse.svg"),
  correlation_ellipse_plot,
  width = 22,
  height = 11,
  units = "cm"
)

# Identify highly abundant proteins without a detectable sex difference

full_report <- read.csv(
  file.path(processed_dir, "Conc_TIMES.csv")
) %>%
  dplyr::select(-1) %>%
  filter(!grepl("P09_TrP12_27", File.Name, ignore.case = TRUE)) %>%
  pivot_longer(
    cols = -c(1:7),
    names_to = "Genes",
    values_to = "nM"
  ) %>%
  filter(
    !is.na(nM),
    !grepl("\\.", Genes)
  )

stats_all <- full_report %>%
  group_split(Genes) %>%
  map_dfr(fit_gene_kruskal) %>%
  mutate(padj = p.adjust(P.unadj, method = "BH"))

genes_no_sex_diff <- stats_all %>%
  filter(padj >= 0.05) %>%
  pull(Genes) %>%
  unique()

top_stable_abundant <- full_report %>%
  group_by(Genes) %>%
  summarise(
    mean_nM = mean(nM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Genes %in% genes_no_sex_diff) %>%
  arrange(desc(mean_nM)) %>%
  slice_head(n = 3) %>%
  pull(Genes)

fc_table <- full_report %>%
  group_by(Genes, Sex) %>%
  summarise(mean_nM = mean(nM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Sex, values_from = mean_nM) %>%
  mutate(fold_change = female / male)

calculate_Sex_effect_size <- function(data) {
  female <- data %>% filter(Sex == "female") %>% pull(nM)
  male <- data %>% filter(Sex == "male") %>% pull(nM)
  
  result <- effsize::cliff.delta(
    female,
    male,
    conf.level = 0.95
  )
  
  tibble(
    effect_size = unname(result$estimate),
    ci_low = result$conf.int[1],
    ci_high = result$conf.int[2]
  )
}

make_stable_plot <- function(gene_name) {
  data <- full_report %>% filter(Genes == gene_name)
  
  n_table <- data %>%
    group_by(Sex) %>%
    summarise(n = sum(!is.na(nM)), .groups = "drop")
  
  effect <- calculate_Sex_effect_size(data)
  
  padj <- stats_all %>%
    filter(Genes == gene_name) %>%
    pull(padj)
  
  fold_change <- fc_table %>%
    filter(Genes == gene_name) %>%
    pull(fold_change)
  
  y_max <- max(data$nM, na.rm = TRUE)
  
  ggstatsplot::ggbetweenstats(
    data = data,
    x = Sex,
    y = nM,
    type = "nonparametric",
    pairwise.display = "none",
    bf.message = FALSE,
    results.subtitle = FALSE,
    xlab = "",
    ylab = "concentration [nM]",
    ggtheme = theme_light(),
    violin.args = list(width = 0.6, alpha = 0.2, linewidth = 0.2),
    boxplot.args = list(width = 0.3, alpha = 0.04, linewidth = 0.15),
    point.args = list(
      position = position_jitterdodge(dodge.width = 0.5),
      alpha = 0.4,
      size = 0.8,
      stroke = 0.5
    ),
    centrality.label.args = list(size = 1.8, nudge_y = y_max*0.25, direction = "y", segment.linetype = 4,
                                 min.segment.length = 0,label.size = 1.8,                      
                                 label.padding = unit(0.1, "lines")),
    centrality.point.args = list(
      size  = 1,
      colour = "darkred"
    )
  ) +
    scale_fill_manual(values = c(female = "#e5b822", male = "#04afbb")) +
    scale_color_manual(values = c(female = "#e5b822", male = "#04afbb")) +
    scale_y_continuous(
      labels = scales::scientific,
      expand = expansion(mult = c(0.05, 0.09))
    ) +
    scale_x_discrete(
      labels = c(
        female = paste0("n = ", n_table$n[n_table$Sex == "female"]),
        male = paste0("n = ", n_table$n[n_table$Sex == "male"])
      )
    ) +
    annotate(
      "text",
      x = 1.5,
      y = y_max * 1.12,
      label = paste0(
        "p[FDR-adj] == ",
        formatC(padj, format = "e", digits = 2)
      ),
      parse = TRUE,
      size = 1.9
    ) +
    labs(
      title = paste0(
        "FC = ", round(fold_change, 2),
        ", δ = ", round(effect$effect_size, 2),
        ", 95% CI = [",
        round(effect$ci_low, 2), ", ",
        round(effect$ci_high, 2), "]"
      )
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 5),
      axis.title.y = element_text(size = 5),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 5, hjust = 0)
    )
}

stable_dir <- file.path(figure_dir, "stable_abundant_proteins")
dir.create(stable_dir, recursive = TRUE, showWarnings = FALSE)

walk(top_stable_abundant, function(gene) {
  ggsave(
    file.path(stable_dir, paste0("stable_abundant_", gene, ".svg")),
    make_stable_plot(gene),
    width = 5.5,
    height = 4.5,
    units = "cm"
  )
})


# Plot individual longitudinal trajectories for selected proteins

full_report <- read.csv(
  file.path(processed_dir, "Conc_TIMES.csv")
) %>%
  dplyr::select(-1) %>%
  filter(!grepl("P09_TrP12_27", File.Name, ignore.case = TRUE)) %>%
  pivot_longer(
    cols = -c(1:7),
    names_to = "Genes",
    values_to = "nM"
  ) %>%
  filter(!is.na(nM))

plot_data <- full_report %>%
  filter(Genes %in% selected_genes) %>%
  mutate(
    Month = factor(
      c("Jun", "Jul", "Aug", "Sept", "Oct", "Nov",
        "Dec", "Jan", "Feb", "Mar", "Apr", "May")[Time.Point],
      levels = month_order
    ),
    Genes = factor(Genes, levels = selected_genes),
    Sex = factor(Sex, levels = c("female", "male"))
  )

female_subjects <- plot_data %>%
  filter(Sex == "female") %>%
  distinct(Subject) %>%
  arrange(Subject) %>%
  pull(Subject)

male_subjects <- plot_data %>%
  filter(Sex == "male") %>%
  distinct(Subject) %>%
  arrange(Subject) %>%
  pull(Subject)

plot_data <- plot_data %>%
  mutate(
    Subject = factor(
      Subject,
      levels = c(female_subjects, male_subjects)
    )
  )

trajectory_plot <- ggplot(
  plot_data,
  aes(
    x = Subject,
    y = nM,
    colour = Month,
    group = interaction(Subject, Genes)
  )
) +
  geom_line(
    colour = "grey65",
    linewidth = 0.25,
    alpha = 0.6
  ) +
  geom_point(size = 1.5, alpha = 0.9) +
  scale_colour_manual(values = month_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::scientific) +
  ggh4x::facet_grid2(
    Genes ~ Sex,
    scales = "free",
    independent = "x",
    space = "free_x"
  ) +
  ggh4x::force_panelsizes(cols = c(1.9, 1.1)) +
  labs(
    x = "participant",
    y = "concentration [nM]",
    colour = "Sampling month"
  ) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5.5),
    axis.title = element_text(size = 8),
    strip.text = element_text(size = 7, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.width = unit(0.4, "cm"),
    panel.spacing = unit(0.2, "lines")
  ) +
  guides(colour = guide_legend(nrow = 1))

ggsave(
  file.path(figure_dir, "individual_longitudinal_profiles.svg"),
  trajectory_plot,
  width = 20,
  height = 22,
  units = "cm"
)
