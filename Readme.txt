miRNA Differential Expression Explorer - README
------------------------------------------------

🧬 Overview:
This interactive Shiny app enables exploratory and analytical workflows for miRNA expression data. From differential expression to machine learning and power analysis, the app is designed to support bioinformatics workflows in a user-friendly environment.

🔹 Step-by-Step Guide:

1️⃣ Upload Data:
   - Use the **sidebar** to upload:
     • A count matrix CSV (genes x samples)
     • A metadata CSV (samples x condition/attributes)
   - Make sure sample names match between both files.

2️⃣ Run Differential Expression:
   - Click **"Run DE Analysis"**
   - The app will run DESeq2, apply variance stabilizing transformation, and identify top DE miRNAs.
   - View results in:
     • **Volcano Plot**
     • **PCA / UMAP**
     • **Heatmap**
     • **Top miRNA table and barplot**

3️⃣ Functional Enrichment:
   - Select a database (e.g., KEGG, GO) and use enrichment buttons for:
     • All DE miRNAs
     • Upregulated
     • Downregulated
   - Results include tables + enrichment barplots with -log10(adj. p-value).

4️⃣ Classification with Random Forest:
   - Click **"Run Random Forest Classification"**
   - Based on the top DE features, a model is trained/tested (70/30 split).
   - Outputs:
     • Sample predictions
     • Accuracy, Sensitivity, Specificity, AUC
     • Variable importance

5️⃣ Power Analysis:
   - Click **"Run Power Analysis"**
   - Power is computed for each miRNA using a t-test across a range of sample sizes (n = 5 to 100).
   - Outputs:
     • Power table (downloadable)
     • Boxplot of sample sizes required to reach 80% power

📝 Downloadable Results:
- Differential expression tables
- Enrichment results
- Random forest predictions, metrics, importance
- Power analysis summary

📦 Requirements:
- Gene/miRNA names must be consistent across input files
- Ensure that `"Primary Tumor"` is the condition of interest (used as the positive class)

💡 Pro Tips:
- Use filtering (e.g., padj < 0.1) to focus on meaningful DE results
- Run PCA/UMAP to visually assess sample separability
- Use power analysis to justify experimental sample sizes

Happy exploring! 🚀

— JCAP MiRNA-SEQ by John Caperella
