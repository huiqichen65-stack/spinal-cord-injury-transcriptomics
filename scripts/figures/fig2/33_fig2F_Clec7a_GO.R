source("config.R", local = TRUE)

options(stringsAsFactors = FALSE)

required_pkgs <- c("ggplot2", "dplyr", "svglite")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

ROOT <- PROJECT_ROOT
RES_DIR <- file.path(ROOT, "results")

STEP22_DIR <- file.path(RES_DIR, "22_Clec7a_GO_KEGG_enrichment")

candidate_go_files <- c(file.path(STEP22_DIR, "22_Clec7a_GO_BP_up_v03.csv"), file.path(STEP22_DIR, "22_Clec7a_GO_BP_up_v02.csv"), 
    file.path(STEP22_DIR, "22_Clec7a_GO_BP_up_v01.csv"))

GO_FILE <- candidate_go_files[file.exists(candidate_go_files)][1]

if (is.na(GO_FILE) || is.null(GO_FILE)) {
    stop("No Step 22 GO BP result found under:\n", STEP22_DIR)
}

OUT_DIR <- file.path(RES_DIR, "33_Figure2F_Clec7a_GO_barplot")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_BASE <- file.path(OUT_DIR, "33_Figure2F_Clec7a_GO_barplot_v06")

go <- read.csv(GO_FILE, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("ID", "Description", "Count", "pvalue", "p.adjust")

missing_cols <- setdiff(required_cols, colnames(go))

if (length(missing_cols) > 0) {
    stop("Required GO column(s) missing: ", paste(missing_cols, collapse = ", "))
}

go$Description <- as.character(go$Description)
go$Count <- as.numeric(go$Count)
go$pvalue <- as.numeric(go$pvalue)
go$p.adjust <- as.numeric(go$p.adjust)

go <- go %>% filter(is.finite(Count), Count > 0, is.finite(p.adjust), p.adjust > 0) %>% mutate(neglog10_padj = -log10(p.adjust))

manual_terms <- c(

  "endocytosis",
  "import into cell",
  "receptor-mediated endocytosis",
  "small GTPase-mediated signal transduction",

  "glycerolipid metabolic process",
  "neutral lipid catabolic process",
  "acylglycerol catabolic process",

  "response to insulin",
  "insulin receptor signaling pathway",
  "cellular response to insulin stimulus",

  "cell activation",
  "leukocyte activation",
  "innate immune response",
  "microglial cell activation",
  "lipopolysaccharide-mediated signaling pathway"
)

plot_df <- go %>% filter(Description %in% manual_terms) %>% distinct(Description, .keep_all = TRUE)

missing_manual <- setdiff(manual_terms, plot_df$Description)

if (length(missing_manual) > 0) {
    warning("The following manually curated GO terms were not found:\n", paste(missing_manual, collapse = "\n"))
}

if (nrow(plot_df) < 8) {
    stop("Too few curated GO terms were found in the selected GO table. ", "Please verify that v03 is being read.")
}

plot_df <- plot_df %>% arrange(Count, p.adjust)

plot_df$Description_display <- factor(plot_df$Description, levels = plot_df$Description)

write.csv(plot_df %>% select(ID, Description, Count, pvalue, p.adjust, neglog10_padj, everything()), 
    paste0(OUT_BASE, "_terms.csv"), row.names = FALSE)

print(plot_df %>% select(Description, Count, pvalue, p.adjust) %>% arrange(desc(Count)))

max_count <- max(plot_df$Count, na.rm = TRUE)

label_x <- max_count * 0.015

p <- ggplot(plot_df, aes(x = Count, y = Description_display, fill = neglog10_padj)) + geom_col(width = 0.8, 
    linewidth = 0) + geom_text(aes(x = label_x, label = Description), hjust = 0, family = "Arial", size = 3.55, 
    color = "black") + scale_fill_gradientn(colours = c("#FFF7A8", "#F6D56A", "#F0A94E", "#E66A3C", "#C92F3A"), 
    name = expression(-log[10](p.adjust))) + scale_x_continuous(limits = c(0, max_count * 1.05), expand = c(0, 
    0)) + labs(title = "Clec7a+ Microglia GO EnrichmentBarplot", x = "Number of Gene", y = "pathway") + 
    theme_classic(base_family = "Arial", base_size = 12) + theme(text = element_text(family = "Arial", 
    color = "black"), plot.title = element_text(size = 14, face = "plain", hjust = 0.5, margin = margin(b = 5)), 
    axis.title.x = element_text(size = 12.5, margin = margin(t = 5)), axis.title.y = element_text(size = 12.5, 
        margin = margin(r = 6)), axis.text.x = element_text(size = 11.8, color = "black"), axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", 
        fill = NA, linewidth = 0.55), panel.grid = element_blank(), legend.position = "right", legend.title = element_text(size = 12), 
    legend.text = element_text(size = 11.2), legend.key.height = grid::unit(1.15, "cm"), legend.key.width = grid::unit(0.3, 
        "cm"), plot.margin = margin(5, 7, 5, 5))

W_IN <- 7.8
H_IN <- 6.3

ggsave(filename = paste0(OUT_BASE, ".svg"), plot = p, width = W_IN, height = H_IN, units = "in", device = svglite::svglite, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".pdf"), plot = p, width = W_IN, height = H_IN, units = "in", device = cairo_pdf, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".png"), plot = p, width = W_IN, height = H_IN, units = "in", dpi = 600, 
    bg = "white")
