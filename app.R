# ====================== app.R ======================
library(shiny)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(plotly)
library(umap)
library(pheatmap)
library(enrichR)
library(randomForest)
library(pROC)

# ==== Load static miRNA-gene associations (CSV fallback) ====
mirna_static_targets <- read.csv("data/mirna_targets.csv", stringsAsFactors = FALSE)

get_targets_from_csv <- function(mirnas, static_df = mirna_static_targets) {
  static_df %>% filter(mirna %in% mirnas)
}

# ==== Helper: Enrichment Plot ====
renderEnrichPlot <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  top_terms <- df %>% arrange(Adjusted.P.value) %>% head(10)
  
  plot_ly(
    data = top_terms,
    x = ~reorder(Term, -Adjusted.P.value),
    y = ~-log10(Adjusted.P.value),
    type = "bar",
    text = ~paste("Term:", Term,
                  "<br>Adj. P:", signif(Adjusted.P.value, 3),
                  "<br>Overlap:", Overlap),
    hoverinfo = "text",
    marker = list(color = "steelblue")
  ) %>%
    layout(
      title = "-log10(Adj. P) of Top Enriched Pathways",
      xaxis = list(title = "Pathway", tickangle = -45),
      yaxis = list(title = "-log10(Adjusted P-value)"),
      margin = list(b = 120)
    )
}

# ==== Database choices ====
db_choices <- c(
  "GO_Biological_Process_2021",
  "GO_Molecular_Function_2021",
  "GO_Cellular_Component_2021",
  "KEGG_2021_Human",
  "Reactome_2022"
)

# ==== Safe wrapper around enrichR call ====
safe_enrichr <- function(genes, db = "KEGG_2021_Human") {
  if (length(genes) < 2) return(data.frame())
  res <- tryCatch({
    result <- enrichR::enrichr(genes, databases = db)
    if (is.null(result[[db]]) || nrow(result[[db]]) == 0) return(data.frame())
    result[[db]]
  }, error = function(e) {
    message("❌ Enrichment error: ", e$message)
    return(data.frame())
  })
  return(res)
}

# ==== UI ====
ui <- fluidPage(
  titlePanel("miRNA Differential Expression (DESeq2)"),
  sidebarLayout(
    sidebarPanel(
      fileInput("countsFile", "Upload Count Matrix CSV", accept = ".csv"),
      fileInput("metaFile", "Upload Sample Metadata CSV", accept = ".csv"),
      actionButton("runAnalysis", "Run DE Analysis"),
      hr(),
      actionButton("volcanoBtn", "Generate Volcano Plot"),
      actionButton("pcaBtn", "Generate PCA Plot"),
      actionButton("umapBtn", "Generate UMAP Plot"),
      actionButton("heatmapBtn", "Generate Heatmap"),
      actionButton("barplotBtn", "Generate Top miRNAs Barplot"),
      hr(),
      actionButton("enrichAllBtn", "Enrichment: All DE miRNAs"),
      actionButton("enrichUpBtn", "Enrichment: Upregulated"),
      actionButton("enrichDownBtn", "Enrichment: Downregulated"),
      selectInput(
        inputId = "selectedDB",
        label = "Select Enrichment Database",
        choices = db_choices,
        selected = "KEGG_2021_Human"
      ),
      actionButton("plotEnrichAllBtn", "Plot: All DE Enrichment"),
      actionButton("plotEnrichUpBtn", "Plot: Upregulated Enrichment"),
      actionButton("plotEnrichDownBtn", "Plot: Downregulated Enrichment"),
      hr(),
      actionButton("runRF", "Run Random Forest Classification"),
      hr(),
      downloadButton("downloadTop", "Download Top DE miRNAs")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Volcano", plotlyOutput("volcanoPlot")),
        tabPanel("PCA", plotlyOutput("pcaPlot")),
        tabPanel("UMAP", plotlyOutput("umapPlot")),
        tabPanel("Heatmap", plotlyOutput("heatmapPlot")),
        tabPanel("Top miRNAs", tableOutput("topTable")),
        tabPanel("Top miRNA Barplot", plotlyOutput("barplotPlot")),
        tabPanel("Enrich: All DE",
                 tabsetPanel(
                   tabPanel("Enrichment Table",
                            downloadButton("downloadEnrichAll", "Download CSV"),
                            br(), br(), tableOutput("enrichAllTable")),
                   tabPanel("Enrichment Barplot", plotlyOutput("enrichAllPlot"))
                 )
        ),
        tabPanel("Enrich: Upregulated",
                 tabsetPanel(
                   tabPanel("Enrichment Table",
                            downloadButton("downloadEnrichUp", "Download CSV"),
                            br(), br(), tableOutput("enrichUpTable")),
                   tabPanel("Enrichment Barplot", plotlyOutput("enrichUpPlot"))
                 )
        ),
        tabPanel("Enrich: Downregulated",
                 tabsetPanel(
                   tabPanel("Enrichment Table",
                            downloadButton("downloadEnrichDown", "Download CSV"),
                            br(), br(), tableOutput("enrichDownTable")),
                   tabPanel("Enrichment Barplot", plotlyOutput("enrichDownPlot"))
                 )
        ),
        tabPanel("Random Forest",
                 tabsetPanel(
                   tabPanel("Predictions", downloadButton("downloadRFpred", "Download Predictions"), br(), br(), tableOutput("rfPredictions")),
                   tabPanel("Metrics", downloadButton("downloadRFmetrics", "Download Metrics"), br(), br(), tableOutput("rfMetrics")),
                   tabPanel("Variable Importance", downloadButton("downloadRFimportance", "Download Importance Table"), br(), br(), tableOutput("rfImportance"))
                 )
        )
      )
    )
  )
)
# ==== SERVER ====
server <- function(input, output, session) {
  # Reactive values
  resultsDF <- reactiveVal(NULL)
  enrichment <- reactiveVal(NULL)
  ddsData <- reactiveVal(NULL)
  vsdData <- reactiveVal(NULL)
  heatmapMatrix <- reactiveVal(NULL)
  
  enrichAll <- reactiveVal(NULL)
  enrichUp <- reactiveVal(NULL)
  enrichDown <- reactiveVal(NULL)
  genesAll <- reactiveVal(NULL)
  genesUp <- reactiveVal(NULL)
  genesDown <- reactiveVal(NULL)
  
  rf_preds <- reactiveVal(NULL)
  rf_metrics <- reactiveVal(NULL)
  rf_importance <- reactiveVal(NULL)
  
  # 🚀 Analysis, plots, enrichment, and classification
  #source("server_core.R", local = TRUE)  # Optional if you want to modularize
  
  # -- Enrichment ALL
  observeEvent(input$enrichAllBtn, {
    req(resultsDF(), input$selectedDB)
    showNotification("⚙️ Starting enrichment for RF-ranked miRNAs...", type = "message")
    mirnas <- resultsDF()$miRNA
    tg_df <- get_targets_from_csv(mirnas, mirna_static_targets)
    tg <- unique(tg_df$target)
    genesAll(tg)
    
    if (length(tg) < 1) {
      showNotification("❌ No targets found", type = "error")
      enrichAll(data.frame())
      return(NULL)
    }
    
    result <- safe_enrichr(tg, input$selectedDB)
    enrichAll(result)
    showNotification("✅ Enrichment complete (All RF)", type = "message")
  })
  
  output$enrichAllTable <- renderTable({ req(enrichAll()); enrichAll() })
  
  output$downloadEnrichAll <- downloadHandler(
    filename = function() { "enrichment_all_rf_ranked.csv" },
    content = function(file) {
      write.csv(enrichAll(), file, row.names = FALSE)
    }
  )
  
  observeEvent(input$plotEnrichAllBtn, {
    req(enrichAll())
    output$enrichAllPlot <- renderPlotly({
      renderEnrichPlot(enrichAll())
    })
  })
  
  # -- Enrichment UP
  observeEvent(input$enrichUpBtn, {
    req(resultsDF(), input$selectedDB)
    showNotification("⚙️ Enrichment: Upregulated miRNAs...", type = "message")
    mirnas <- resultsDF() %>% filter(!is.na(padj), padj < 0.1, log2FoldChange > 0) %>% pull(miRNA)
    tg_df <- get_targets_from_csv(mirnas, mirna_static_targets)
    tg <- unique(tg_df$target)
    genesUp(tg)
    
    if (length(tg) < 1) {
      showNotification("❌ No targets found (upregulated)", type = "error")
      enrichUp(data.frame())
      return(NULL)
    }
    
    result <- safe_enrichr(tg, input$selectedDB)
    enrichUp(result)
    showNotification("✅ Enrichment done (Upregulated)", type = "message")
  })
  
  output$enrichUpTable <- renderTable({ req(enrichUp()); enrichUp() })
  
  output$downloadEnrichUp <- downloadHandler(
    filename = function() { "enrichment_upregulated_rf.csv" },
    content = function(file) {
      write.csv(enrichUp(), file, row.names = FALSE)
    }
  )
  
  observeEvent(input$plotEnrichUpBtn, {
    req(enrichUp())
    output$enrichUpPlot <- renderPlotly({
      renderEnrichPlot(enrichUp())
    })
  })
  
  # -- Enrichment DOWN
  observeEvent(input$enrichDownBtn, {
    req(resultsDF(), input$selectedDB)
    showNotification("⚙️ Enrichment: Downregulated miRNAs...", type = "message")
    mirnas <- resultsDF() %>% filter(!is.na(padj), padj < 0.1, log2FoldChange < 0) %>% pull(miRNA)
    tg_df <- get_targets_from_csv(mirnas, mirna_static_targets)
    tg <- unique(tg_df$target)
    genesDown(tg)
    
    if (length(tg) < 1) {
      showNotification("❌ No targets found (downregulated)", type = "error")
      enrichDown(data.frame())
      return(NULL)
    }
    
    result <- safe_enrichr(tg, input$selectedDB)
    enrichDown(result)
    showNotification("✅ Enrichment done (Downregulated)", type = "message")
  })
  
  output$enrichDownTable <- renderTable({ req(enrichDown()); enrichDown() })
  
  output$downloadEnrichDown <- downloadHandler(
    filename = function() { "enrichment_downregulated_rf.csv" },
    content = function(file) {
      write.csv(enrichDown(), file, row.names = FALSE)
    }
  )
  
  observeEvent(input$plotEnrichDownBtn, {
    req(enrichDown())
    output$enrichDownPlot <- renderPlotly({
      renderEnrichPlot(enrichDown())
    })
  })
  observeEvent(input$runAnalysis, {
    req(input$countsFile, input$metaFile)
    showNotification("🧬 Performing DESeq2 Analysis...", type = "message")
    
    # === Load & Prepare Count Data ===
    count_data <- read_csv(input$countsFile$datapath)
    count_matrix <- as.matrix(count_data[,-1])
    rownames(count_matrix) <- count_data[[1]]
    
    # === Load & Prepare Metadata ===
    meta_data <- read_csv(input$metaFile$datapath)
    meta_data <- as.data.frame(meta_data)
    rownames(meta_data) <- meta_data[[1]]
    meta_data <- meta_data[,-1]
    meta_data$condition <- as.factor(meta_data[[ncol(meta_data)]])
    meta_data <- meta_data[colnames(count_matrix), , drop = FALSE]
    
    # === DESeq2 Setup ===
    dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                  colData = meta_data,
                                  design = ~ condition)
    dds <- DESeq(dds)
    ddsData(dds)
    
    vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
    vsdData(vsd)
    
    # === Extract DE Results ===
    res <- results(dds)
    res$miRNA <- rownames(res)
    res_df <- as.data.frame(res) %>% mutate(sig = padj < 0.1)
    
    showNotification("✅ DE completed. Running RF feature ranking...", type = "message")
    
    # === Top 20 DE by padj
    top20 <- res_df %>%
      arrange(padj) %>%
      filter(!is.na(padj)) %>%
      slice_head(n = 20)
    
    vst_mat <- assay(vsd)
    
    # === Clean Names for RF
    vst_mat_renamed <- vst_mat
    rownames(vst_mat_renamed) <- make.names(rownames(vst_mat))
    top20$miRNA_clean <- make.names(top20$miRNA)
    
    vst_subset <- vst_mat_renamed[rownames(vst_mat_renamed) %in% top20$miRNA_clean, , drop = FALSE]
    df_rf <- as.data.frame(t(vst_subset))
    df_rf$condition <- factor(vsdData()$condition)
    
    # === Random Forest
    rf <- randomForest(condition ~ ., data = df_rf, importance = TRUE)
    imp <- importance(rf, type = 2)
    
    imp_df <- data.frame(miRNA_clean = rownames(imp), Importance = imp[, 1]) %>%
      arrange(desc(Importance)) %>%
      slice_head(n = 11)
    
    name_map <- data.frame(
      miRNA = top20$miRNA,
      miRNA_clean = top20$miRNA_clean,
      stringsAsFactors = FALSE
    )
    
    imp_df <- merge(imp_df, name_map, by = "miRNA_clean", all.x = TRUE)
    
    final_hits <- res_df %>%
      filter(miRNA %in% imp_df$miRNA)
    
    # 🛡️ Sanity check & fallback
    if (!is.data.frame(final_hits) || nrow(final_hits) == 0) {
      showNotification("⚠️ DE results empty or failed to map RF features", type = "error")
      return(NULL)
    }
    
    resultsDF(final_hits)
    
    showNotification("✅ DE pipeline complete with RF-ranked features", type = "message")
  })
  
observeEvent(input$heatmapBtn, {
  req(resultsDF(), vsdData())
  showNotification("🔥 Generating heatmap...", type = "message")
  
  mat <- assay(vsdData())
  
  sig_miRNAs <- resultsDF() %>%
    filter(!is.na(padj)) %>%
    pull(miRNA)
  
  sig_miRNAs <- intersect(sig_miRNAs, rownames(mat))
  
  if (length(sig_miRNAs) < 2) {
    showNotification("❌ Not enough DE miRNAs to plot heatmap", type = "error")
    return(NULL)
  }
  
  mat <- mat[sig_miRNAs, , drop = FALSE]
  mat <- t(scale(t(mat)))  # Z-score
  mat[is.na(mat)] <- 0
  
  row_order <- hclust(dist(mat))$order
  col_order <- hclust(dist(t(mat)))$order
  mat <- mat[row_order, col_order, drop = FALSE]
  
  output$heatmapPlot <- renderPlotly({
    plot_ly(
      z = mat,
      x = colnames(mat),
      y = rownames(mat),
      type = "heatmap",
      colorscale = "Viridis",
      colorbar = list(title = "Z-score"),
      hovertemplate = paste(
        "miRNA: %{y}<br>",
        "Sample: %{x}<br>",
        "Z-score: %{z:.2f}<extra></extra>"
      )
    ) %>%
      layout(
        title = paste0("🧬 Heatmap of ", nrow(mat), " DE miRNAs"),
        xaxis = list(title = "Samples"),
        yaxis = list(title = "miRNAs"),
        margin = list(l = 100, b = 100)
      )
  })
  
  showNotification("✅ Heatmap rendered", type = "message")
})
output$topTable <- renderTable({
  req(resultsDF())
  df <- resultsDF()
  df <- df[order(df$padj), ]
  df <- head(df, 20)
  df$padj <- signif(df$padj, 3)
  df$log2FoldChange <- signif(df$log2FoldChange, 3)
  df
})
observeEvent(input$runRF, {
  req(resultsDF(), vsdData())
  showNotification("🌲 Running Random Forest Classification...", type = "message")
  
  vst_mat <- assay(vsdData())
  final_hits <- resultsDF()
  
  sig_miRNAs <- final_hits$miRNA
  vst_mat_renamed <- vst_mat
  rownames(vst_mat_renamed) <- make.names(rownames(vst_mat))
  
  sig_miRNAs_clean <- make.names(sig_miRNAs)
  
  vst_subset <- vst_mat_renamed[rownames(vst_mat_renamed) %in% sig_miRNAs_clean, , drop = FALSE]
  df_rf <- as.data.frame(t(vst_subset))
  df_rf$condition <- factor(vsdData()$condition)
  
  rf <- randomForest(condition ~ ., data = df_rf, importance = TRUE)
  preds <- predict(rf, type = "response")
  
  # === Predictions
  rf_preds(data.frame(Sample = rownames(df_rf), Predicted = preds))
  
  # === Metrics
  actual <- df_rf$condition
  confusion <- table(Predicted = preds, Actual = actual)
  accuracy <- sum(diag(confusion)) / sum(confusion)
  rf_metrics(data.frame(Metric = "Accuracy", Value = accuracy))
  
  # === Importance
  imp <- importance(rf, type = 2)
  rf_importance(data.frame(miRNA = rownames(imp), Importance = imp[, 1]) %>%
                  arrange(desc(Importance)))
  
  showNotification("✅ Random Forest complete", type = "message")
})
# === Random Forest Outputs ===

output$rfPredictions <- renderTable({
  req(rf_preds())
  rf_preds()
})

output$rfMetrics <- renderTable({
  req(rf_metrics())
  rf_metrics()
})

output$rfImportance <- renderTable({
  req(rf_importance())
  rf_importance()
})

output$downloadTop <- downloadHandler(
  filename = function() {
    paste0("top_de_miRNAs_", Sys.Date(), ".csv")
  },
  content = function(file) {
    req(resultsDF())
    write.csv(resultsDF(), file, row.names = FALSE)
  }
)
output$downloadRFpred <- downloadHandler(
  filename = function() {
    paste0("rf_predictions_", Sys.Date(), ".csv")
  },
  content = function(file) {
    req(rf_preds())
    write.csv(rf_preds(), file, row.names = FALSE)
  }
)
output$downloadRFmetrics <- downloadHandler(
  filename = function() {
    paste0("rf_metrics_", Sys.Date(), ".csv")
  },
  content = function(file) {
    req(rf_metrics())
    write.csv(rf_metrics(), file, row.names = FALSE)
  }
)
output$downloadRFimportance <- downloadHandler(
  filename = function() {
    paste0("rf_variable_importance_", Sys.Date(), ".csv")
  },
  content = function(file) {
    req(rf_importance())
    write.csv(rf_importance(), file, row.names = FALSE)
  }
)


observeEvent(input$pcaBtn, {
  req(vsdData())
  showNotification("📊 Generating PCA plot...", type = "message")
  
  vst <- assay(vsdData())
  vst <- t(vst)
  
  # ⚠️ Remove zero-variance columns (miRNAs with no variation)
  vst <- vst[, apply(vst, 2, function(x) sd(x, na.rm = TRUE) > 1e-8)]
  
  validate(need(ncol(vst) >= 2, "Too few variable features for PCA"))
  
  pca <- prcomp(vst, scale. = TRUE)
  pc_df <- as.data.frame(pca$x[, 1:2])
  pc_df$Sample <- rownames(pc_df)
  pc_df$Group <- vsdData()$condition
  
  output$pcaPlot <- renderPlotly({
    plot_ly(
      data = pc_df,
      x = ~PC1,
      y = ~PC2,
      color = ~Group,
      text = ~Sample,
      type = "scatter",
      mode = "markers",
      marker = list(size = 10, opacity = 0.7),
      hoverinfo = "text"
    ) %>%
      layout(
        title = "PCA: Principal Component Analysis",
        xaxis = list(title = "PC1"),
        yaxis = list(title = "PC2")
      )
  })
  
  showNotification("✅ PCA plot rendered", type = "message")
})



observeEvent(input$umapBtn, {
  req(vsdData())
  showNotification("🔮 Generating UMAP plot...", type = "message")
  
  vst <- assay(vsdData())
  vst <- t(vst)
  
  umap_result <- umap::umap(vst)
  umap_df <- as.data.frame(umap_result$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$Sample <- rownames(vst)
  umap_df$Group <- vsdData()$condition
  
  output$umapPlot <- renderPlotly({
    plot_ly(
      data = umap_df,
      x = ~UMAP1,
      y = ~UMAP2,
      color = ~Group,
      text = ~Sample,
      type = "scatter",
      mode = "markers",
      marker = list(size = 10, opacity = 0.7),
      hoverinfo = "text"
    ) %>%
      layout(
        title = "UMAP: Uniform Manifold Approximation and Projection",
        xaxis = list(title = "UMAP1"),
        yaxis = list(title = "UMAP2")
      )
  })
  
  showNotification("✅ UMAP plot rendered", type = "message")
})
observeEvent(input$volcanoBtn, {
  req(resultsDF())
  showNotification("🌋 Rendering Volcano plot...", type = "message")
  
  df <- resultsDF() %>%
    mutate(
      log10padj = -log10(padj),
      sig_label = case_when(
        padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
        padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
        TRUE ~ "Not Sig"
      )
    )
  
  output$volcanoPlot <- renderPlotly({
    plot_ly(
      data = df,
      x = ~log2FoldChange,
      y = ~log10padj,
      text = ~paste("miRNA:", miRNA,
                    "<br>log2FC:", signif(log2FoldChange, 3),
                    "<br>padj:", signif(padj, 3)),
      type = "scatter",
      mode = "markers",
      color = ~sig_label,
      colors = c("gray", "firebrick", "steelblue"),
      marker = list(size = 8, opacity = 0.8),
      hoverinfo = "text"
    ) %>%
      layout(
        title = "Volcano Plot: log2FC vs -log10(padj)",
        xaxis = list(title = "log2(Fold Change)"),
        yaxis = list(title = "-log10(padj)"),
        showlegend = TRUE
      )
  })
  
  showNotification("✅ Volcano plot ready", type = "message")
})
observeEvent(input$barplotBtn, {
  req(resultsDF())
  showNotification("📊 Rendering barplot of top miRNAs...", type = "message")
  
  df <- resultsDF() %>%
    arrange(padj) %>%
    filter(!is.na(padj)) %>%
    slice_head(n = 20)
  
  df <- df %>%
    mutate(miRNA = factor(miRNA, levels = rev(miRNA)))
  
  output$barplotPlot <- renderPlotly({
    plot_ly(
      data = df,
      x = ~log2FoldChange,
      y = ~miRNA,
      type = "bar",
      orientation = "h",
      text = ~paste("padj:", signif(padj, 3)),
      hoverinfo = "text",
      marker = list(color = "mediumseagreen")
    ) %>%
      layout(
        title = "Top 20 DE miRNAs (log2 Fold Change)",
        xaxis = list(title = "log2(Fold Change)"),
        yaxis = list(title = "", tickangle = 0),
        margin = list(l = 120)
      )
  })
  
  showNotification("✅ Barplot ready", type = "message")
})













}
# ==== Final App Run ====
shinyApp(ui, server)


