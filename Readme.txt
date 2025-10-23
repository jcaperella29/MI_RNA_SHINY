miRNA Differential Expression Explorer - README
------------------------------------------------

🧬 Overview:
This interactive Shiny app enables exploratory and analytical workflows for miRNA expression data across multiple species (human, mouse, zebrafish, and fruit fly). 
It supports everything from differential expression to enrichment, machine learning, and power analysis — all in one intuitive interface.  

Built for speed, the app now includes session caching and automatic species detection for seamless cross-species exploration.

🔹 Step-by-Step Guide:

1️⃣ Upload Data:
   - Use the **sidebar** to upload:
     • A count matrix CSV (miRNAs × samples)
     • A metadata CSV (samples × condition/attributes)
   - Make sure sample names match between both files.
   - Supported species:
     • Human (hsa)
     • Mouse (mmu)
     • Zebrafish (dre)
     • Fruit Fly (dme)

2️⃣ Run Differential Expression:
   - Click **"Run DE Analysis"**
   - The app will run DESeq2, apply variance stabilizing transformation, and identify top DE miRNAs.
   - View results in:
     • Volcano Plot
     • PCA / UMAP
     • Heatmap
     • Top miRNA table and barplot

3️⃣ Functional Enrichment — Multi-Species Aware:
   - The app automatically detects species (from miRNA prefixes like hsa-, mmu-, dre-, dme-).
   - The enrichment database menu updates automatically for the right organism.
   - Buttons available:
     • Enrichment: All DE miRNAs
     • Enrichment: Upregulated
     • Enrichment: Downregulated
   - Each produces:
     • Enrichment tables
     • Interactive barplots (-log10 adjusted p-values)

   🔬 Backends used:
     • multiMiR + miRNAtap for validated/predicted targets
     • Enrichr, clusterProfiler, or g:Profiler for pathway analysis

   ⚡ Caching:
     • Results are cached per session — repeated enrichment for the same species + miRNA list loads instantly.

4️⃣ Classification with Random Forest:
   - Click **"Run Random Forest Classification"**
   - Based on top DE features, the model trains/test (70/30 split)
   - Outputs:
     • Sample predictions
     • Accuracy, Sensitivity, Specificity, AUC
     • Variable importance

5️⃣ Power Analysis:
   - Click **"Run Power Analysis"**
   - Power is computed per miRNA (t-test) across n = 5–100 samples
   - Outputs:
     • Power table (downloadable)
     • Boxplot of sample size needed for 80% power

📝 Downloadable Results:
- DE tables
- Enrichment results
- Random Forest predictions, metrics, importance
- Power Analysis summaries

⚙️ Performance and Caching:
- Target caching speeds up multiMiR and miRNAtap queries.
- Automatic species detection chooses the right enrichment databases.
- Falls back to clusterProfiler when network-based tools fail.

📦 Requirements:
- Gene/miRNA names must match between files.
- "Primary Tumor" is the positive class for classification.
- Supports Human, Mouse, Zebrafish, and Fruit Fly.

💡 Pro Tips:
- Use padj < 0.1 to filter significant DE miRNAs.
- Run PCA/UMAP to confirm sample separation.
- Use Power Analysis before experiments to justify sample sizes.
- Caching makes iterative enrichment testing fast.

🚀 Summary:
JCAP miRNA-SEQ App v2 — Multi-species ready, cached, and optimized for reproducible bioinformatics workflows.

— John Caperella (JCAP Bioinformatics Hub)
