#!/bin/bash
# 16S analysis script using QIIME2

# Create and activate conda environment for QIIME2
wget https://data.qiime2.org/distro/core/qiime2-2023.2-py38-linux-conda.yml
conda env create -n qiime2-2023.2 --file qiime2-2023.2-py38-linux-conda.yml
conda activate qiime2-2023.2

# Install greengenes2 plugin
pip install q2-greengenes2

# Download required sequences
bash ena-file-PRJNA801673_16s.sh

# Create manifest.csv using python script
python 16s_manifest.py

# Import sequences into QIIME2 format
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.csv \
  --output-path demux.qza \
  --input-format PairedEndFastqManifestPhred33

# Summarize and visualize sequence quality
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

# Execute DADA2 for denoising and chimera removal
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 19 \
  --p-trim-left-r 22 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

# Generate additional visualizations
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file metadata.tsv

# Download necessary references
wget http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.fna.qza
wget http://ftp.microbio.me/greengenes_release/2022.10/2022.10.taxonomy.asv.nwk.qza

# Map feature table and assign taxonomy
qiime greengenes2 non-v4-16s \
    --i-table table.qza \
    --i-sequences rep-seqs.qza \
    --i-backbone 2022.10.backbone.full-length.fna.qza \
    --o-mapped-table 16s.gg2.biom.qza \
    --o-representatives 16s.gg2.fna.qza

qiime greengenes2 taxonomy-from-table \
     --i-reference-taxonomy 2022.10.taxonomy.asv.nwk.qza \
     --i-table 16s.gg2.biom.qza \
     --o-classification 16s.gg2.taxonomy.qza

# Final visualizations
qiime metadata tabulate \
  --m-input-file 16s.gg2.taxonomy.qza \
  --o-visualization 16s.gg2.taxonomy.qzv

qiime taxa barplot --i-table 16s.gg2.biom.qza \
                   --i-taxonomy 16s.gg2.taxonomy.qza \
                   --m-metadata-file metadata.tsv \
                   --o-visualization taxa_barplot.qzv

qiime feature-table summarize \
  --i-table 16s.gg2.biom.qza \
  --o-visualization 16s.gg2.biom.qzv \
  --m-sample-metadata-file metadata.tsv

# Results available for visualization at: https://view.qiime2.org/
