
## ================== Fig1 TITAN2 Threshold Plot + Heatmap ==================
rm(list = ls())

library(ggplot2)
library(dplyr)
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
library(tidyr)
library(tibble)
library(pheatmap)
library(cowplot)
library(grid)
library(RColorBrewer)
library(gridExtra)
library(colorspace)

load("all_TITAN2_data.RData")

global_abs_z <- abs(spp_all$zscore)
z_breaks <- pretty(c(0, max(global_abs_z, na.rm = TRUE)), n = 4)
z_breaks <- z_breaks[z_breaks >= 0]

fsumz <- sumz %>%
  filter(tao %in% c("fsumz-", "fsumz+")) %>%
  select(tao, cp_value = cp) %>%
  mutate(Direction = ifelse(tao == "fsumz-", "Decreaser", "Increaser"))
tau_dec <- fsumz %>% filter(Direction == "Decreaser") %>% pull(cp_value)
tau_inc <- fsumz %>% filter(Direction == "Increaser") %>% pull(cp_value)

tax_data$Phylum <- gsub("^p__", "", tax_data$Phylum)

soften_colors <- function(cols, desat = 0.20, lighten = 0.08) {
  cols %>%
    desaturate(amount = desat) %>%
    lighten(amount = lighten)
}

my_phylum_colors <- c(
  "#1b9e77", "#d95f02","#a65628", "#7570b3",  "#66a61e",
  "#fdc086","#a6761d", "#666666", "#8dd3c7", "#bebada",
  "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
  "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f", "#4daf4a",
  "#984ea3", "#fcc5c0", "#fa9fb5", "#f768a1", "#1f78b4",
  "#999999", "#a6cee3", "#7fc97f","#d73027",  "#beaed4",
  "#e6ab02"
)

my_colors_soft <- soften_colors(my_phylum_colors, desat = 0.06, lighten = 0.04)

all_phyla <- sort(unique(tax_data$Phylum))
if(length(my_colors_soft) < length(all_phyla)) {
  my_colors_soft <- rep(my_colors_soft, length.out = length(all_phyla))
}
phylum_colors <- setNames(my_colors_soft, all_phyla)
ann_colors <- list(Phylum = phylum_colors)

equalize_pheatmap_matrix_width <- function(p1, p2) {
  layout1 <- p1$gtable$layout
  layout2 <- p2$gtable$layout
  idx1 <- which(layout1$name == "matrix")
  idx2 <- which(layout2$name == "matrix")
  
  w1 <- convertWidth(sum(p1$gtable$widths[layout1$l[idx1]:layout1$r[idx1]]),"cm", valueOnly = TRUE)
  w2 <- convertWidth(sum(p2$gtable$widths[layout2$l[idx2]:layout2$r[idx2]]),"cm", valueOnly = TRUE)
  
  target <- max(w1, w2)
  
  p1$gtable$widths[layout1$l[idx1]:layout1$r[idx1]] <-
    p1$gtable$widths[layout1$l[idx1]:layout1$r[idx1]] * (target / w1)
  p2$gtable$widths[layout2$l[idx2]:layout2$r[idx2]] <-
    p2$gtable$widths[layout2$l[idx2]:layout2$r[idx2]] * (target / w2)
  
  list(p1 = p1, p2 = p2)
}

clean_colnames <- function(x) {
  x <- gsub("^Atom", "", x)
  x <- gsub("^0", "0.", x)
  return(x)
}

dec_data <- spp_df %>%
  filter(purity >= 0.95, reliability >= 0.95, maxgrp == 1) %>%
  left_join(tax_data %>% select(OTU_ID, Phylum),
            by = c("OTU" = "OTU_ID")) %>%
  arrange(desc(zenv.cp))

if(nrow(dec_data) == 0) {
  stop("No decreaser species found")
}

otu_z_dec <- t(apply(otu_dec[dec_data$OTU, ], 1,
                     function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))
otu_z_dec[is.na(otu_z_dec)] <- 0
otu_z_dec <- otu_z_dec[rev(rownames(otu_z_dec)), ]
heat_order_dec <- rownames(otu_z_dec)

row_annot_dec <- data.frame(Phylum = dec_data$Phylum)
rownames(row_annot_dec) <- dec_data$OTU
row_annot_dec <- row_annot_dec[heat_order_dec, , drop = FALSE]
colnames(otu_z_dec) <- clean_colnames(colnames(otu_z_dec))

pheat_dec <- pheatmap(
  otu_z_dec,
  cluster_rows = FALSE, cluster_cols = FALSE,
  annotation_row = row_annot_dec,
  annotation_colors = ann_colors,
  annotation_legend = FALSE,
  legend = FALSE,
  show_rownames = FALSE, show_colnames = TRUE,
  color = colorRampPalette(c("#2483C7","white","#EC681C"))(20),
  fontsize = 7, fontsize_col = 7, fontsize_row = 6, angle_col = 45,
  border_color = "white"
)

p_dec <- ggplot(dec_data, aes(y = factor(OTU, levels = dec_data$OTU))) +
  geom_hline(yintercept = seq_along(heat_order_dec), color = "grey90", size = 0.3) +
  geom_errorbarh(aes(xmin = X5., xmax = X95.), height = 0.3, color = "#5e92cd", size = 0.8) +
  geom_point(aes(x = zenv.cp, size = abs(zscore)), shape = 21, stroke = 0.8, fill = "#5e92cd", color = "black") +
  geom_vline(xintercept = tau_dec, color = "#5e92cd", linetype = "dashed", size = 0.8) +
  scale_x_continuous(name = "atom% 18O", limits = c(0,80), breaks = seq(0,80,10)) +
  theme(axis.title.x = element_blank()) +
  scale_y_discrete(position = "right") +
  scale_size_continuous(guide = "none", range = c(1.2, 5)) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size=7, color = "black"),
        axis.text.x = element_text(size=7, color = "black"),
        panel.border = element_rect(color = "black", size = 0.5),
        plot.margin = ggplot2::margin(10, 15, 10, 0))

inc_data <- spp_df %>%
  filter(purity >= 0.95, reliability >= 0.95, maxgrp == 2) %>%
  left_join(tax_data %>% select(OTU_ID, Phylum),
            by = c("OTU" = "OTU_ID")) %>%
  arrange(zenv.cp)

if(nrow(inc_data) == 0) {
  stop("No increaser species found")
}

otu_z_inc <- t(apply(otu_inc[inc_data$OTU, ], 1,
                     function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))
otu_z_inc[is.na(otu_z_inc)] <- 0
otu_z_inc <- otu_z_inc[rev(rownames(otu_z_inc)), ]
heat_order_inc <- rownames(otu_z_inc)

row_annot_inc <- data.frame(Phylum = inc_data$Phylum)
rownames(row_annot_inc) <- inc_data$OTU
row_annot_inc <- row_annot_inc[heat_order_inc, , drop = FALSE]
colnames(otu_z_inc) <- clean_colnames(colnames(otu_z_inc))

pheat_inc <- pheatmap(
  otu_z_inc,
  cluster_rows = FALSE, cluster_cols = FALSE,
  annotation_row = row_annot_inc,
  annotation_colors = ann_colors,
  annotation_legend = FALSE,
  legend = FALSE,
  show_rownames = FALSE, show_colnames = TRUE,
  color = colorRampPalette(c("#2483C7","white","#EC681C"))(20),
  fontsize = 8, fontsize_col = 8, fontsize_row = 6, angle_col = 45,
  border_color = "white"
)

p_inc <- ggplot(inc_data, aes(y = factor(OTU, levels = inc_data$OTU))) +
  geom_hline(yintercept = seq_along(inc_data$OTU), color = "grey90", size = 0.3) +
  geom_errorbarh(aes(xmin = X5., xmax = X95.), height = 0.3, color = "#f0904d", size = 0.8) +
  geom_point(aes(x = zenv.cp, size = abs(zscore)), shape = 21, stroke = 0.8, fill = "#f0904d", color = "black") +
  geom_vline(xintercept = tau_inc, color = "#f0904d", linetype = "dashed", size = 0.8) +
  scale_x_continuous(name = "atom% 18O", limits = c(0,80), breaks = seq(0,80,10)) +
  theme(axis.title.x = element_blank()) +
  scale_y_discrete(limits = inc_data$OTU) +
  scale_size_continuous(guide = "none", range = c(1.2, 5)) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size=8, color = "black"),
        axis.text.x = element_text(size=9),
        panel.border = element_rect(color = "black", size = 1),
        plot.margin = ggplot2::margin(10, 0, 10, 15))

eq_pheats <- equalize_pheatmap_matrix_width(pheat_dec, pheat_inc)
pheat_dec_eq <- eq_pheats$p1
pheat_inc_eq <- eq_pheats$p2

combo_dec <- ggdraw() +
  draw_plot(plot_grid(ggdraw(pheat_dec_eq$gtable), ggdraw(ggplotGrob(p_dec)), 
                      ncol = 2, align = "hv", rel_widths = c(1, 1)),
            x = 0, y = 0.03, width = 1, height = 0.97) +
  draw_label("atom% 18O", x = 0.5, y = 0.01, hjust = 0.5, vjust = 0, size = 8) +
  theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

combo_inc <- ggdraw() +
  draw_plot(plot_grid(ggdraw(ggplotGrob(p_inc)), ggdraw(pheat_inc_eq$gtable),
                      ncol = 2, align = "hv", rel_widths = c(1, 1)),
            x = 0, y = 0.03, width = 1, height = 0.97) +
  draw_label("atom% 18O", x = 0.5, y = 0.01, hjust = 0.5, vjust = 0, size = 12) +
  theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

phy_used <- sort(unique(c(row_annot_dec$Phylum, row_annot_inc$Phylum)))
phylum_colors_used <- phylum_colors[phy_used]

legend_phylum <- cowplot::get_legend(
  ggplot(data.frame(Phylum = factor(phy_used, levels = phy_used)),
         aes(x = 1, y = Phylum, fill = Phylum)) +
    geom_tile(width = 0.8, height = 0.8, color = "grey30", linewidth = 0.2) +
    scale_fill_manual(values = phylum_colors_used) +
    guides(fill = guide_legend(title = "Phylum", ncol = 1, byrow = TRUE)) +
    theme_void() +
    theme(legend.title = element_text(size = 11, hjust = 0),
          legend.text = element_text(size = 9),
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.35, "cm"),
          legend.spacing.y = unit(0.15, "cm"))
)

heat_legend <- cowplot::get_legend(
  ggplot(data.frame(x = 1:3, y = 1:3, z = c(-2, 0, 2))) +
    geom_tile(aes(x, y, fill = z)) +
    scale_fill_gradient2(low = "#2483c7", mid = "white", high = "#EC681C",
                         midpoint = 0, name = "Relative abundance (Z-score)") +
    theme_minimal() +
    theme(legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          legend.key.height = unit(0.6, "cm"),
          legend.key.width = unit(0.4, "cm"))
)

legend_compact <- plot_grid(legend_phylum, heat_legend, ncol = 1, rel_heights = c(1.2, 0.8))

two_panels <- plot_grid(combo_dec, combo_inc, ncol = 2, rel_widths = c(0.48, 0.48), align = "hv") +
  theme(plot.margin = ggplot2::margin(0, -25, 0, -25))

final_fig <- plot_grid(two_panels, legend_compact, ncol = 2, rel_widths = c(0.85, 0.15), align = "h")

print(final_fig)

## =========== Fig3 ======================= ##
rm(list = ls())

## ===============================================================
## Threshold analysis from integrated RData
## Variable definitions:
##   CO2            : Respiration
##   concentration  : DNA content
##   EAFDNA         : EAF-18O of DNA
## ===============================================================

## ---------------------------
## 1. Package management
## ---------------------------
required_pkgs <- c("segmented", "MuMIn")

install_if_missing <- function(pkgs) {
  missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing_pkgs) > 0) {
    install.packages(missing_pkgs, dependencies = TRUE)
  }
}

install_if_missing(required_pkgs)

suppressPackageStartupMessages({
  library(segmented)
  library(MuMIn)
})

## ---------------------------
## 2. Parameters
## ---------------------------
infile <- "Fig3_data.RData"

n_boot <- 500
min_n  <- 8
seed   <- 2025
set.seed(seed)

out_data_dir <- "data/processed"
out_tab_dir  <- "results"

for (dir_i in c(out_data_dir, out_tab_dir)) {
  if (!dir.exists(dir_i)) dir.create(dir_i, recursive = TRUE)
}

## ---------------------------
## 3. Load integrated RData
## ---------------------------
if (!file.exists(infile)) {
  stop(
    "Cannot find input RData file: ", infile,
    "\nPlease run 01_prepare_all_data.R first."
  )
}

load(infile)

## Expected objects from 01_prepare_all_data.R:
## figure4_df, figure4_x_var, figure4_vars
## figure3_df, figure3_x_var, figure3_vars
required_objects <- c(
  "figure4_df", "figure4_x_var", "figure4_vars",
  "figure3_df", "figure3_x_var", "figure3_vars"
)

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop(
    "The following required objects are missing from the RData file: ",
    paste(missing_objects, collapse = ", ")
  )
}

## ---------------------------
## 4. Core threshold analysis function
## ---------------------------
analyze_one <- function(d, x, y, n_boot = 500, min_n = 8) {
  dd <- d[, c(x, y)]
  names(dd) <- c("x", "y")
  
  dd$x <- suppressWarnings(as.numeric(dd$x))
  dd$y <- suppressWarnings(as.numeric(dd$y))
  dd <- dd[complete.cases(dd), , drop = FALSE]
  
  ## Result template: only keep requested statistics
  empty_row <- function(model_label = "too few data",
                        tau = NA_real_,
                        boot_lo = NA_real_,
                        boot_hi = NA_real_,
                        davies_p = NA_real_,
                        delta_aicc = NA_real_) {
    data.frame(
      Variable = y,
      Tau = tau,
      Boot_CI_lower = boot_lo,
      Boot_CI_upper = boot_hi,
      Davies_p = davies_p,
      Delta_AICc = delta_aicc,
      Model = model_label,
      N = nrow(dd),
      stringsAsFactors = FALSE
    )
  }
  
  if (nrow(dd) < min_n || sd(dd$y, na.rm = TRUE) == 0) {
    return(list(row = empty_row("too few data"), model = NULL))
  }
  
  ## 1) Linear model
  m0 <- lm(y ~ x, data = dd)
  AICc_linear <- MuMIn::AICc(m0)
  
  ## 2) Try segmented model using several initial breakpoint values
  seg_fit <- NULL
  init_pool <- unique(as.numeric(quantile(dd$x, c(0.25, 0.50, 0.75), na.rm = TRUE)))
  
  for (psi0 in init_pool) {
    seg_try <- try(
      segmented::segmented(
        m0,
        seg.Z = ~x,
        psi = list(x = psi0),
        control = segmented::seg.control(n.boot = 0, it.max = 500)
      ),
      silent = TRUE
    )
    
    if (!inherits(seg_try, "try-error")) {
      seg_fit <- seg_try
      break
    }
  }
  
  ## 3) If segmented model fails, retain linear model
  if (is.null(seg_fit)) {
    return(list(
      row = empty_row(
        model_label = "linear (segmented failed)",
        delta_aicc = NA_real_
      ),
      model = m0
    ))
  }
  
  ## 4) AICc comparison and Davies test
  AICc_segment <- MuMIn::AICc(seg_fit)
  Delta_AICc <- AICc_segment - AICc_linear
  ## Delta_AICc < 0 means segmented model has lower AICc than linear model.
  
  dv_p <- try(segmented::davies.test(m0, seg.Z = ~x)$p.value, silent = TRUE)
  dv_p <- if (inherits(dv_p, "try-error")) NA_real_ else as.numeric(dv_p)
  
  ## Keep segmented only when segmented has lower AICc and Davies test supports slope change.
  use_segmented <- is.finite(AICc_segment) &&
    is.finite(AICc_linear) &&
    (AICc_segment < AICc_linear) &&
    is.finite(dv_p) &&
    dv_p < 0.05
  
  ## 5) If segmented model is not supported, retain linear model.
  ## Tau and bootstrap CI are NA because no breakpoint is retained.
  if (!use_segmented) {
    return(list(
      row = empty_row(
        model_label = "linear (AICc/Davies)",
        davies_p = dv_p,
        delta_aicc = Delta_AICc
      ),
      model = m0
    ))
  }
  
  ## 6) Extract breakpoint
  tau_est <- as.numeric(seg_fit$psi[1, 2])
  
  ## 7) Bootstrap breakpoint CI only
  taus <- numeric(0)
  n <- nrow(dd)
  
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    db <- dd[idx, , drop = FALSE]
    
    ok <- try({
      mb <- lm(y ~ x, data = db)
      sb <- segmented::segmented(
        mb,
        seg.Z = ~x,
        psi = list(x = median(db$x, na.rm = TRUE)),
        control = segmented::seg.control(n.boot = 0, it.max = 500)
      )
      as.numeric(sb$psi[1, 2])
    }, silent = TRUE)
    
    if (!inherits(ok, "try-error") && is.finite(ok)) {
      taus <- c(taus, ok)
    }
  }
  
  if (length(taus) >= max(30, 0.3 * n_boot)) {
    boot_ci <- quantile(taus, c(0.025, 0.975), na.rm = TRUE)
    b_lo <- as.numeric(boot_ci[1])
    b_hi <- as.numeric(boot_ci[2])
    
    if (is.finite(b_lo) && is.finite(b_hi) && b_lo > b_hi) {
      tmp <- b_lo
      b_lo <- b_hi
      b_hi <- tmp
    }
  } else {
    b_lo <- b_hi <- NA_real_
  }
  
  list(
    row = empty_row(
      model_label = "segmented",
      tau = tau_est,
      boot_lo = b_lo,
      boot_hi = b_hi,
      davies_p = dv_p,
      delta_aicc = Delta_AICc
    ),
    model = seg_fit
  )
}

run_threshold_analysis <- function(data, x_var, vars, dataset_name,
                                   n_boot = 500, min_n = 8) {
  vars <- intersect(vars, names(data))
  
  if (length(vars) == 0) {
    stop("No valid response variables found for ", dataset_name)
  }
  
  res_list <- lapply(vars, function(v) {
    analyze_one(data, x_var, v, n_boot = n_boot, min_n = min_n)
  })
  names(res_list) <- vars
  
  res_tbl <- do.call(rbind, lapply(res_list, `[[`, "row"))
  res_tbl$Dataset <- dataset_name
  res_tbl <- res_tbl[, c("Dataset", setdiff(names(res_tbl), "Dataset"))]
  
  list(
    dataset = dataset_name,
    x_var = x_var,
    vars = vars,
    results = res_tbl,
    models = lapply(res_list, `[[`, "model")
  )
}

## ---------------------------
## 5. Run Figure4 and Figure3 analyses
## ---------------------------
figure4_analysis <- run_threshold_analysis(
  data = figure4_df,
  x_var = figure4_x_var,
  vars = figure4_vars,
  dataset_name = "Figure4",
  n_boot = n_boot,
  min_n = min_n
)

figure3_analysis <- run_threshold_analysis(
  data = figure3_df,
  x_var = figure3_x_var,
  vars = figure3_vars,
  dataset_name = "Figure3_EAF_DNA",
  n_boot = n_boot,
  min_n = min_n
)

figure4_results <- figure4_analysis$results
figure3_results <- figure3_analysis$results
all_threshold_results <- rbind(figure4_results, figure3_results)

print(all_threshold_results)
