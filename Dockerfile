# ---- Dockerfile ----
FROM rocker/shiny:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libtiff-dev \
    libjpeg-dev \
    libgit2-dev \
    libglpk-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    && apt-get clean

# Install required R packages
RUN R -e "install.packages(c('shiny', 'DESeq2', 'ggplot2', 'dplyr', 'readr', 'plotly', 'umap', 'pheatmap', 'clusterProfiler', 'org.Hs.eg.db', 'enrichR', 'randomForest', 'pROC', 'pwr', 'DT', 'BiocManager'), repos='http://cran.r-project.org')" \
    && R -e "BiocManager::install(c('DESeq2', 'org.Hs.eg.db', 'clusterProfiler'))"

# Copy app files
COPY . /srv/shiny-server/

# Change working directory
WORKDIR /srv/shiny-server/

# Expose Shiny port
EXPOSE 3838

# Start the app
CMD ["/usr/bin/shiny-server"]
