#!/bin/bash
# WGS analysis script using Bowtie2 and Woltka
module load metawrap/1.3

# Download sequences
# Define base directory path
BASE_DIR="/home/proyectos/imdeaalim/alba"
cd "$BASE_DIR"
chmod +x ena-PRJNA798058.sh
./ena-PRJNA798058.sh

### 2. Quality check of sequences and adapter removal
mkdir "$BASE_DIR/READ_QC"

# Execute read quality control (QC)
for F in "$BASE_DIR"/*_1.fastq.gz; do
    R=${F%_*}_2.fastq.gz

    # Decompress files
    gunzip -c "$F" > "${F%.gz}" && gunzip -c "$R" > "${R%.gz}"

    FD="${F%.gz}"
    RD="${R%.gz}"

    BASE=${FD##*/}
    SAMPLE=${BASE%_*}

    # Check if the output directory already exists
    if [ ! -d "$BASE_DIR/READ_QC/$SAMPLE" ]; then
        metawrap read_qc -1 "$FD" -2 "$RD" -t 6 -o "$BASE_DIR/READ_QC/$SAMPLE"
    else
        echo "Directory $BASE_DIR/READ_QC/$SAMPLE already exists. Skipping processing."
    fi

    # Clean up decompressed files
    rm "$FD" "$RD"
done > "$BASE_DIR/read_qcERROR_01-10-2024-11-54.txt" 2>&1

# Move QC-passed files to a new directory (CLEAN_READS) and rename them
mkdir "$BASE_DIR/CLEAN_READS"

for i in $BASE_DIR/READ_QC/*; do
    BASE=${i##*/}
    SAMPLE=${BASE%_*}
    mv "${i}/final_pure_reads_1.fastq" "$BASE_DIR/CLEAN_READS/${SAMPLE}_1.fastq"
    mv "${i}/final_pure_reads_2.fastq" "$BASE_DIR/CLEAN_READS/${SAMPLE}_2.fastq"
done



# Load Bowtie2 module and execute Perl script to align samples
module load bowtie2/2.5.3
cd "$BASE_DIR/CLEAN_READS"
perl bowtie2.pl *_1.fastq

# Create directory for results and move SAM files
mkdir -p /home/proyectos/imdeaalim/alba/CLEAN_READS/alineamientos/SamFiles/results
mv *.sam /home/proyectos/imdeaalim/alba/CLEAN_READS/alineamientos/SamFiles/

# Execute Woltka for classification
module load woltka/0.1.6
woltka classify -i /home/proyectos/imdeaalim/alba/CLEAN_READS/alineamientos/SamFiles \
                -f sam \
                -o /home/proyectos/imdeaalim/alba/CLEAN_READS/alineamientos/SamFiles/results/ogu.biom

# Import Woltka-generated biom into Qiime2
module load qiime/2023.2
conda activate qiime2-2023.2
qiime tools import \
     --input-path /home/proyectos/imdeaalim/alba/CLEAN_READS/alineamientos/SamFiles/results/ogu.biom \
     --output-path ogu.biom.qza \
     --type FeatureTable[Frequency]

# Download taxonomy file
wget http://ftp.microbio.me/greengenes_release/2022.10/2022.10.taxonomy.asv.nwk.qza

# Filter and assign taxonomy
qiime greengenes2 filter-features \
     --i-feature-table ogu.biom.qza \
     --i-reference 2022.10.taxonomy.asv.nwk.qza \
     --o-filtered-feature-table woltka_gg2.biom.qza

qiime greengenes2 taxonomy-from-table \
     --i-reference-taxonomy 2022.10.taxonomy.asv.nwk.qza \
     --i-table woltka_gg2.biom.qza \
     --o-classification woltka_gg2.taxonomy.qza

# Generate final visualization
qiime metadata tabulate \
  --m-input-file woltka_gg2.taxonomy.qza \
  --o-visualization woltka_gg2.taxonomy.qzv

# Results available for visualization at: https://view.qiime2.org/