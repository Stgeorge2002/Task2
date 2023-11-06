#!/bin/bash

# Exit script if any command fails
set -e

# Function to check and install system tools if not present
check_and_install() {
    package=$1
    install_command=$2

    if ! command -v "$package" &> /dev/null; then
        echo "Package $package could not be found. Installing..."
        eval "$install_command"
    else
        echo "Package $package is already installed."
    fi
}

# Check and install system dependencies
check_and_install "wget" "sudo apt-get install wget -y"
check_and_install "unzip" "sudo apt-get install unzip -y"
check_and_install "gzip" "sudo apt-get install gzip -y"
check_and_install "bzip2" "sudo apt-get install bzip2 -y"
check_and_install "fastqc" "sudo apt-get install fastqc -y"
check_and_install "hisat2" "sudo apt-get install hisat2 -y"
check_and_install "samtools" "sudo apt-get install samtools -y"
check_and_install "subread" "sudo apt-get install subread -y"
check_and_install "trim-galore" "sudo apt-get install trim-galore -y"
check_and_install "r-base" "sudo apt-get install r-base -y"
check_and_install "r-base-dev" "sudo apt-get install r-base-dev -y"
check_and_install "cd-hit" "sudo apt-get install cd-hit -y"
check_and_install "clustalw" "sudo apt-get install clustalw -y"
check_and_install "emboss" "sudo apt-get install emboss -y"
# Add check for sra-toolkit (fastq-dump)
check_and_install "fastq-dump" "sudo apt-get install sra-toolkit -y"

# Set up R and install necessary packages
R_PACKAGES=("BiocManager" "DESeq2" "sva" "WGCNA" "ggplot2" "pheatmap" "ggtree" "KEGGREST" "Biostrings" "PopGenome" "ape")

# Install R packages
for pkg in "${R_PACKAGES[@]}"; do
    sudo Rscript -e "if (!requireNamespace('$pkg', quietly = TRUE)) { if('$pkg' == 'BiocManager') { install.packages('$pkg') } else { BiocManager::install('$pkg') } }"
done

# Define directories for raw and processed data, databases and results.
RAW_DATA_DIR="$HOME/Analysis/RawData"
PROCESSED_DATA_DIR="$HOME/Analysis/ProcessedData"
GENOME_DIR="$HOME/Analysis/Genome"
RESULTS_DIR="$HOME/Analysis/Results"

# Create directories if they don't exist
mkdir -p "$RAW_DATA_DIR"
mkdir -p "$PROCESSED_DATA_DIR"
mkdir -p "$GENOME_DIR"
mkdir -p "$RESULTS_DIR"

# Download reference genomes and annotations
echo "Downloading reference genomes and annotation files..."
# Replace the FTP URLs with the actual URLs from NCBI or other sources
wget -P "$GENOME_DIR" "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz"
wget -P "$GENOME_DIR" "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz"

# Download genome annotation files
# MC
wget -P "$GENOME_DIR" "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gtf.gz"
# Ecoli
wget -P "$GENOME_DIR" "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.gtf.gz"

# Unzip the downloaded files
gunzip -k "$GENOME_DIR"/*.gz

# Download short-read data using SRA Toolkit
# Hard-coded SRR accession numbers
echo "Downloading short-read data..."
declare -a SRR_IDS=("SRR15278553" "SRR12746045") # Replace with your actual SRR IDs

for srr_id in "${SRR_IDS[@]}"; do
    fastq-dump --split-files "$srr_id" --outdir "$RAW_DATA_DIR" --gzip
done

# Quality Control and Read Trimming
echo "Performing QC and read trimming..."
for fq in "$RAW_DATA_DIR"/*_1.fastq.gz; do
    # Assuming paired-end data: _1 and _2 files
    fq2=${fq/_1.fastq.gz/_2.fastq.gz}
    fastqc "$fq" "$fq2" -o "$PROCESSED_DATA_DIR"
    trim_galore --paired "$fq" "$fq2" --output_dir "$PROCESSED_DATA_DIR"
done

# Align Reads to Reference Genomes and Prepare for Feature Counting
echo "Aligning reads to reference genomes..."
hisat2-build "$GENOME_DIR"/GCF_000008865.2_ASM886v2_genomic.fna "$GENOME_DIR"/mg_index
hisat2-build "$GENOME_DIR"/GCF_000027325.1_ASM2732v1_genomic.fna "$GENOME_DIR"/ec_index

for fq in "$PROCESSED_DATA_DIR"/*_val_1.fq.gz; do
    base_name=$(basename "$fq" _val_1.fq.gz)
    hisat2 -p 8 -x "$GENOME_DIR"/mg_index -1 "$fq" -2 "${fq/_val_1.fq.gz/_val_2.fq.gz}" -S "$PROCESSED_DATA_DIR"/"${base_name}_mg.sam"
    hisat2 -p 8 -x "$GENOME_DIR"/ec_index -1 "$fq" -2 "${fq/_val_1.fq.gz/_val_2.fq.gz}" -S "$PROCESSED_DATA_DIR"/"${base_name}_ec.sam"
done

# Convert SAM to BAM, Sort, and Index
echo "Processing SAM files..."
for sam_file in "$PROCESSED_DATA_DIR"/*.sam; do
    bam_file=${sam_file/.sam/.bam}
    samtools view -S -b "$sam_file" > "$bam_file"
    samtools sort -o "${bam_file/.bam/_sorted.bam}" "$bam_file"
    samtools index "${bam_file/.bam/_sorted.bam}"
done

# Feature Counting
echo "Performing feature counting..."
featureCounts -a "$GENOME_DIR"/GCF_000027325.1_ASM2732v1_genomic.gtf -o "$PROCESSED_DATA_DIR"/mg_counts.txt "$PROCESSED_DATA_DIR"/*_mg_sorted.bam
featureCounts -a "$GENOME_DIR"/GCF_000008865.2_ASM886v2_genomic.gtf -o "$PROCESSED_DATA_DIR"/ec_counts.txt "$PROCESSED_DATA_DIR"/*_ec_sorted.bam

# The script now assumes all tools are installed and databases are in place
# Further analysis like gene duplication, biosynthetic pathway analysis, etc. are done here

run_r_analysis() {
    Rscript -e "
    # Load necessary libraries
    library(DESeq2)
    library(sva)
    library(WGCNA)
    library(ggplot2)
    library(pheatmap)
    library(ggtree)
    library(KEGGREST)
    library(Biostrings)
    library(PopGenome)
    library(ape)

    # Define the data directories
    raw_data_dir <- '$RAW_DATA_DIR'
    processed_data_dir <- '$PROCESSED_DATA_DIR'
    results_dir <- '$RESULTS_DIR'

    # Create results directory if it doesn't exist
    if (!dir.exists(results_dir)) {
        dir.create(results_dir)
    }

    # Read in the count data
    mg_counts <- read.csv(file.path(processed_data_dir, 'mg_counts.txt'), sep='\t', header=TRUE)
    ec_counts <- read.csv(file.path(processed_data_dir, 'ec_counts.txt'), sep='\t', header=TRUE)

    # Preprocessing for DESeq2
    dds_mg <- DESeqDataSetFromMatrix(countData = mg_counts[,4:ncol(mg_counts)],
                                     colData = DataFrame(condition = factor(c('control', 'treatment'))),
                                     design = ~ condition)

    dds_ec <- DESeqDataSetFromMatrix(countData = ec_counts[,4:ncol(ec_counts)],
                                     colData = DataFrame(condition = factor(c('control', 'treatment'))),
                                     design = ~ condition)

    # Differential expression analysis
    dds_mg <- DESeq(dds_mg)
    dds_ec <- DESeq(dds_ec)

    res_mg <- results(dds_mg)
    res_ec <- results(dds_ec)

    # Write results to file
    write.csv(as.data.frame(res_mg), file=file.path(results_dir, 'DE_results_mg.csv'))
    write.csv(as.data.frame(res_ec), file=file.path(results_dir, 'DE_results_ec.csv'))

    # Plotting results with ggplot2
    ggplot(as.data.frame(res_mg), aes(x=log2FoldChange, y=-log10(pvalue))) +
        geom_point() +
        theme_minimal() +
        ggtitle('M. genitalium Differential Expression')

    # Batch effect correction with sva
    mod <- model.matrix(~ condition, colData(dds_mg))
    mod0 <- model.matrix(~ 1, colData(dds_mg))
    svseq <- svaseq(counts(dds_mg), mod, mod0)

    # WGCNA analysis for co-expression network
    datExpr <- as.data.frame(t(assay(dds_mg)))
    names(datExpr) <- colnames(dds_mg)
    powers <- c(c(1:10), seq(from = 12, to=20, by=2))
    sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
    net <- blockwiseModules(datExpr, power = sft$powerEstimate, maxBlockSize = 5000)

    # Plotting a dendrogram of the genes
    plotDendroAndColors(net$labels, net$colors, 'Module colors')

    # KEGG Pathway Analysis
    kegg_genes <- substr(names(res_mg), start=1, stop=15)
    kegg_res <- keggLink('genes', kegg_genes)
    pathway <- keggGet(kegg_res)

    # Save KEGG results to a file
    saveRDS(pathway, file = file.path(results_dir, 'kegg_pathway_analysis.rds'))

    # Analysis with PopGenome
    gen <- readData(raw_data_dir, format='VCF')
    gen_stats <- neutrality.test(gen, test='Tajima.D')

    # Phylogenetic analysis
    dna <- read.dna(file.path(processed_data_dir, 'alignment.fasta'), format='fasta')
    phy_tree <- nj(dist.dna(dna))
    plot(phy_tree, main='Phylogenetic Tree')

    # Save the plotted phylogenetic tree as a PDF
    pdf(file.path(results_dir, 'Phylogenetic_Tree.pdf'))
    plot(phy_tree, main='Phylogenetic Tree')
    dev.off()

    # Ensure all results are saved and available for later
    save.image(file.path(results_dir, 'R_analysis_workspace.RData'))

    # Print session information for reproducibility
    sessionInfo()
    "
}
run_r_analysis

echo "Analysis pipeline completed successfully."
