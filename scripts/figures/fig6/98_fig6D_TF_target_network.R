options(stringsAsFactors = FALSE)
options(width = 180)

indir <- "D:/zz/results/98_Figure6D_Cytoscape_input_preparation"

nodes_file <- file.path(indir, "01_98_Figure6D_Cytoscape_nodes_v01.csv")

edges_file <- file.path(indir, "02_98_Figure6D_Cytoscape_edges_v01.csv")

outdir <- file.path(indir, "GraphML_v05")

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

graphml_file <- file.path(outdir, "01_98_Figure6D_Cytoscape_network_v05.graphml")

qc_file <- file.path(outdir, "02_98_GraphML_QC_v05.csv")

report_file <- file.path(outdir, "03_98_GraphML_REPORT_v05.txt")

color_audit_file <- file.path(outdir, "04_98_Figure6D_final_node_colors_v05.csv")

tf_fill_final <- "#D9D9D9"

target_fill_low <- "#FFF5F0"
target_fill_mid <- "#FB6A4A"
target_fill_high <- "#CB181D"

if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required. Install it with:\n", "install.packages(\"igraph\")")
}

required_files <- c(nodes_file, edges_file)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
    stop("Missing Step98 input file(s):\n", paste(missing_files, collapse = "\n"))
}

map_target_log2fc_to_hex <- function(
  x,
  low = target_fill_low,
  mid = target_fill_mid,
  high = target_fill_high
) {

  x <- as.numeric(x)

  out <- rep(
    "#BDBDBD",
    length(x)
  )

  ok <- is.finite(x)

  if (!any(ok)) {
    return(out)
  }

  xr <- range(
    x[ok],
    na.rm = TRUE
  )

  pal <- grDevices::colorRampPalette(
    c(
      low,
      mid,
      high
    )
  )(101)

  if (
    !all(is.finite(xr)) ||
    diff(xr) == 0
  ) {

    idx <- rep(
      51L,
      sum(ok)
    )

  } else {

    scaled <- (
      x[ok] - xr[1]
    ) / diff(xr)

    idx <- pmax(
      1L,
      pmin(
        101L,
        as.integer(
          round(
            scaled * 100
          ) + 1L
        )
      )
    )
  }

  out[ok] <- pal[idx]

  out
}

nodes <- read.csv(nodes_file, stringsAsFactors = FALSE, check.names = FALSE)

edges <- read.csv(edges_file, stringsAsFactors = FALSE, check.names = FALSE)

nodes$cy_fill_hex_final <- tf_fill_final

is_target_node <- nodes$node_type == "Target"

nodes$cy_fill_hex_final[is_target_node] <- map_target_log2fc_to_hex(nodes$log2FC[is_target_node])

if (anyNA(nodes$cy_fill_hex_final)) {
    stop("NA values detected in cy_fill_hex_final.")
}

required_node_cols <- c("id", "node_type", "regulon_size", "log2FC", "display_name", "RSS_rank", "target_RSS", 
    "cy_shape", "cy_node_size", "cy_fill_hex", "cy_border_hex", "cy_border_width", "cy_label_font_size", 
    "cy_label_font_style", "cy_fill_hex_final")

required_edge_cols <- c("edge_id", "source", "target_id", "interaction", "directed", "regulon_weight", 
    "cy_edge_width", "cy_edge_color", "cy_edge_transparency", "cy_target_arrow_shape")

if (!all(required_node_cols %in% colnames(nodes))) {
    stop("Node table lacks required column(s): ", paste(setdiff(required_node_cols, colnames(nodes)), 
        collapse = ", "))
}

if (!all(required_edge_cols %in% colnames(edges))) {
    stop("Edge table lacks required column(s): ", paste(setdiff(required_edge_cols, colnames(edges)), 
        collapse = ", "))
}

if (anyDuplicated(nodes$id)) {
  stop("Duplicated node IDs found.")
}

if (anyDuplicated(edges$edge_id)) {
  stop("Duplicated edge IDs found.")
}

edge_nodes <- unique(c(edges$source, edges$target_id))

missing_nodes <- setdiff(edge_nodes, nodes$id)

if (length(missing_nodes) > 0) {
    stop("Edge endpoint(s) absent from node table: ", paste(missing_nodes, collapse = ", "))
}

expected_tf <- c("Mef2a", "Nr3c1", "Pura", "Ikzf1", "Bach1", "Fli1", "Zfp148", "Ep300", "Bclaf1", "Chd1")

observed_tf <- nodes$id[nodes$node_type == "TF"]

if (length(observed_tf) != length(expected_tf) || !setequal(observed_tf, expected_tf)) {
    stop("TF identity mismatch.\nObserved: ", paste(observed_tf, collapse = ", "))
}

vertices <- data.frame(name = nodes$id, nodes[, setdiff(colnames(nodes), "id"), drop = FALSE], check.names = FALSE, 
    stringsAsFactors = FALSE)

vertices$node_id <- vertices$name

if (anyDuplicated(vertices$name)) {
    stop("Duplicated vertex IDs after GraphML preparation: ", paste(unique(vertices$name[duplicated(vertices$name)]), 
        collapse = ", "))
}

char_cols_v <- names(vertices)[vapply(vertices, is.character, logical(1))]

for (cc in char_cols_v) {
  vertices[[cc]][is.na(vertices[[cc]])] <- ""
}

edge_df <- data.frame(from = edges$source, to = edges$target_id, stringsAsFactors = FALSE)

edge_attr_cols <- setdiff(colnames(edges), c("source", "target_id"))

for (cc in edge_attr_cols) {
  edge_df[[cc]] <- edges[[cc]]
}

char_cols_e <- names(edge_df)[vapply(edge_df, is.character, logical(1))]

for (cc in char_cols_e) {
  edge_df[[cc]][is.na(edge_df[[cc]])] <- ""
}

color_audit <- nodes[, c("id", "node_type", "regulon_size", "log2FC", "cy_fill_hex_final"), drop = FALSE]

write.csv(color_audit, color_audit_file, row.names = FALSE)

g <- igraph::graph_from_data_frame(d = edge_df, directed = TRUE, vertices = vertices)

if (colnames(vertices)[1] != "name" || anyDuplicated(vertices$name)) {
    stop("Internal GraphML vertex-key audit failed.")
}

if (igraph::vcount(g) != nrow(nodes)) {
    stop("Vertex count mismatch: graph=", igraph::vcount(g), ", nodes table=", nrow(nodes))
}

if (igraph::ecount(g) != nrow(edges)) {
    stop("Edge count mismatch: graph=", igraph::ecount(g), ", edges table=", nrow(edges))
}

igraph::write_graph(g, file = graphml_file, format = "graphml")

if (!file.exists(graphml_file)) {
  stop("GraphML file was not created.")
}

qc_check <- c("GraphML written", "Unique GraphML vertex key", "Node count", "Edge count", "TF count", 
    "Target count", "All edge endpoints resolved", "node_type embedded", "regulon_size embedded", "log2FC embedded", 
    "display_name embedded", "Final fill color embedded", "TF final fill fixed gray", "Target final fill derived from log2FC", 
    "regulon_weight embedded")

qc_status <- c("PASS", if (colnames(vertices)[1] == "name" && !anyDuplicated(vertices$name)) "PASS" else "FAIL", 
    if (igraph::vcount(g) == 54L) "PASS" else "WARN", if (igraph::ecount(g) == 50L) "PASS" else "WARN", 
    if (sum(nodes$node_type == "TF") == 10L) "PASS" else "FAIL", if (sum(nodes$node_type == "Target") == 
        44L) "PASS" else "FAIL", if (length(missing_nodes) == 0L) "PASS" else "FAIL", if ("node_type" %in% 
        igraph::vertex_attr_names(g)) "PASS" else "FAIL", if ("regulon_size" %in% igraph::vertex_attr_names(g)) "PASS" else "FAIL", 
    if ("log2FC" %in% igraph::vertex_attr_names(g)) "PASS" else "FAIL", if ("display_name" %in% igraph::vertex_attr_names(g)) "PASS" else "FAIL", 
    if ("cy_fill_hex_final" %in% igraph::vertex_attr_names(g)) "PASS" else "FAIL", if (all(nodes$cy_fill_hex_final[nodes$node_type == 
        "TF"] == tf_fill_final)) "PASS" else "FAIL", if (length(unique(nodes$cy_fill_hex_final[nodes$node_type == 
        "Target"])) > 1L) "PASS" else "FAIL", if ("regulon_weight" %in% igraph::edge_attr_names(g)) "PASS" else "FAIL")

qc_detail <- c(graphml_file, paste0(length(unique(vertices$name)), " unique gene IDs; first vertex column = ", 
    colnames(vertices)[1]), paste0(igraph::vcount(g), " nodes"), paste0(igraph::ecount(g), " edges"), 
    paste0(sum(nodes$node_type == "TF"), " TFs"), paste0(sum(nodes$node_type == "Target"), " targets"), 
    if (length(missing_nodes) == 0L) "All resolved" else paste(missing_nodes, collapse = ", "), "Embedded directly in GraphML", 
    "Embedded directly in GraphML", "Embedded directly in GraphML", "Embedded directly in GraphML", "Embedded directly in GraphML", 
    paste0("All TFs = ", tf_fill_final), paste0("Target color range from log2FC; unique colors = ", length(unique(nodes$cy_fill_hex_final[nodes$node_type == 
        "Target"]))), "Embedded directly in GraphML")

if (length(qc_check) != length(qc_status) || length(qc_check) != length(qc_detail)) {
    stop("Internal QC vector-length mismatch: check=", length(qc_check), ", status=", length(qc_status), 
        ", detail=", length(qc_detail))
}

qc <- data.frame(check = qc_check, status = qc_status, detail = qc_detail, stringsAsFactors = FALSE)

write.csv(qc, qc_file, row.names = FALSE)

if (any(qc$status == "FAIL")) {
  stop("GraphML QC contains FAIL item(s).")
}

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    report_file)
