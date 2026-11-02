#!/bin/bash -l
#
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J HiSca
#SBATCH -o HiCScaffolding_%j.out
#SBATCH -e HiCScaffolding_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=300G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

# ===================================================================
# PREREQUISITES
# 1. Conda environments: 'tools_env', 'hicup_env', 'haphic_env'
# 2. Modules: 'juicer', 'bwa', '3d-dna'
# 3. External tools: Juicer scripts/JARs, HapHiC, HiCUP configuration file
# ===================================================================

# -------------------------------------------------------------------
# 1) ENABLE FAIL-FAST BEHAVIOR
# -------------------------------------------------------------------
set -e
set -o pipefail

# -------------------------------------------------------------------
# 2) UNIVERSAL VARIABLES DEFINITION
# -------------------------------------------------------------------

# -- Conda Initialization --
conda_sh="" 
mamba_sh=""

# -- Hardware Allocation --
threads=32
mem=300

# -- Pipeline Logic & Modes --
filt_only=""    # "1" = Filter reads only | "2" = Run complete pipeline
mode=""         # "1" = Generate .hic map for existing assembly | "2" = Scaffold with HapHiC first
filter=""       # "1" = Use Juicer for HiC filtering | "2" = Use HiCUP for HiC filtering

# -- Input Files & Directories --
genome=""         # Absolute path to the input genome FASTA
path=""           # Working directory
suffix=""         # Suffix/prefix for output files
index_dir=""      # Directory to store the Bowtie2/BWA genome index
path_hic_reads="" # Directory containing raw Hi-C reads (.fastq.gz)
sample_name=""    # Common prefix of the raw Hi-C read files (e.g., "Sample_1")

# -- Filtered Hi-C Data (If applicable) --
filt_hic=""       # Path to already filtered HiC BAM (if running HapHiC directly)

# -- Assembly Parameters --
exp_NumChrs="1"   # Expected chromosomes. "1" = Unknown/Unguided. >1 = Guided count.

# -- External Tools & Scripts --
JUICER=""         # Path to Juicer scripts directory
JUICER_TOOLS=""   # Path to juicer_tools.jar file
HAPHIC=""         # Path to HapHiC executable

# -------------------------------------------------------------------
# 3) LOGGING FUNCTION
# -------------------------------------------------------------------
log_file="progress_pipeline_HiC.txt"

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$log_file"
}

#####################################################################
#################### MAIN PIPELINE STARTS HERE ######################
#####################################################################

cd "$path" || { echo "Error: Working directory not found"; exit 1; }
mkdir -p "$index_dir"

log "Script execution started."
log "Mode: $mode | Filter: $filter"

source "$conda_sh"
source "$mamba_sh"

genome_name=$(basename "$genome")
genomeIdx_Name="${genome_name}.fai"

# ===================================================================
# BRANCH A: EXISTING ASSEMBLY + JUICER FILTERING (Mode 1, Filter 1)
# ===================================================================
if [[ "$mode" -eq 1 && "$filter" -eq 1 ]]; then
    log "--- Entering Branch A: Existing Assembly + Juicer ---"
    
    conda activate tools_env
    module load juicer bwa 3d-dna

    mkdir -p juicer_setup/{fastq,references,restriction_sites}
    cp "${path_hic_reads}"/"${sample_name}"*.gz juicer_setup/fastq/

    # Standardize FASTQ naming for Juicer
    cd juicer_setup/fastq/
    for f in ./${sample_name}*.gz; do
      base=$(basename "$f")
      read_num=$(echo "$base" | sed -nE 's/.*([12])\.(f(ast)?q)\.gz$/\1/p')
      if [ -n "$read_num" ]; then
        mv "$f" "${sample_name}_R${read_num}.fastq.gz"
      fi
    done
    cd ../../

    # Generate supporting files
    samtools faidx "$genome"
    cut -f1,2 "${genome}.fai" > juicer_setup/chrom_sizes
    echo -e "chr1\t1" > juicer_setup/restriction_sites/dummy_sites.txt

    # Index Genome
    log "Indexing genome with BWA..."
    cd juicer_setup/references
    cp "$genome" .
    bwa index "$genome_name"
    cd ../../

    # Run Juicer
    log "Running Juicer pipeline..."
    juicer.sh -d juicer_setup -p juicer_setup/chrom_sizes \
      -s none -z juicer_setup/references/"$genome_name" \
      -q long -Q 48:00:00 -l long -L 96:00:00 -t $threads

    conda deactivate 
    rm -rf juicer_setup/fastq/*
    log "Branch A completed successfully."
    exit 0
fi

# ===================================================================
# BRANCH B: HICUP DIGESTION & FILTERING 
# ===================================================================
log "--- Starting HiCUP Digestion & Filtering ---"
conda activate hicup_env

log "Building Bowtie2 index..."
bowtie2-build "$genome" "$index_dir/${suffix}" || { log "Bowtie2 index failed"; exit 1; }

log "Running HiCUP Digester (^GATC,DpnII)..."
hicup_digester -genome "$suffix" -re1 ^GATC,DpnII "$genome" --outdir "$path" || exit 1

mkdir -p filtered_reads

digest=$(ls "$path"/Digest_*.txt | head -n 1)
if [ -z "$digest" ]; then log "Error: No digest file found"; exit 1; fi

# Dynamically update the config file
r1=$(realpath $(ls "$path_hic_reads"/${sample_name}*1*.gz))
r2=$(realpath $(ls "$path_hic_reads"/${sample_name}*2*.gz))

sed -e "s|^Threads:.*|Threads: $threads|" \
    -e "s|^Index:.*|Index: $index_dir/$suffix|" \
    -e "s|^Digest:.*|Digest: $digest|" \
    -e "\|^/hic/reads1/path|c\\$r1" \
    -e "\|^/hic/reads2/path|c\\$r2" \
    pipeline_HiC.conf > pipeline_HiC_updated.conf
mv pipeline_HiC_updated.conf pipeline_HiC.conf

log "Running HiCUP..."
hicup --config pipeline_HiC.conf
conda deactivate

if [ "$filt_only" -eq 1 ]; then
  log "Filtering complete. 'filt_only' set to 1. Exiting early."
  exit 0 
fi

# ===================================================================
# BRANCH C: HAPHIC SCAFFOLDING (Mode 2)
# ===================================================================
if [ "$mode" -eq 2 ]; then
    log "--- Entering Branch C: HapHiC Scaffolding ---"
    conda activate haphic_env
    
    log "Scaffolding with HapHiC (Expected Chrs: $exp_NumChrs)..."
    "$HAPHIC" pipeline "$genome" "$filt_hic" "$exp_NumChrs"
    
    conda deactivate
    log "HapHiC Scaffolding completed."
fi

# ===================================================================
# BRANCH D: .hic FILE GENERATION FOR OMICSBOX
# ===================================================================

# Sub-Branch D1: Pre-existing assembly filtered with HiCUP
if [[ "$mode" -eq 1 && "$filter" -eq 2 ]]; then
  log "--- Sub-Branch D1: Generating .hic from HiCUP alignments ---"
  conda activate tools_env

  samtools faidx "$genome" 
  cut -f1,2 "${genome}.fai" > "${suffix}.chrom.sizes"
  awk '{print $1 "\t1\t" $2 "\t1\tW\t"$1 "\t1\t" $2 "\t+"}' "${suffix}.chrom.sizes" > "${suffix}.agp"

  log "Sorting BAM file..."
  samtools sort -n -o filtered_reads/"${sample_name}-pairsR1_2.hicup.namesort.bam" \
      filtered_reads/${sample_name}*.hicup.bam -@ $threads

  log "Generating Juicer pre-alignments..."
  "$JUICER" pre filtered_reads/"${sample_name}-pairsR1_2.hicup.namesort.bam" \
    "${suffix}.agp" "$genomeIdx_Name" \
    | sort -k2,2d -k6,6d -T ./ --parallel="$threads" -S"${mem}G" \
    | awk 'NF' > alignments_sorted.txt

  log "Creating .hic file..."
  java -Djava.awt.headless=true -jar -Xmx"${mem}G" "$JUICER_TOOLS" \
    pre alignments_sorted.txt "${suffix}_ContMap.hic" "${suffix}.chrom.sizes"

  conda deactivate
  log "Sub-Branch D1 completed successfully."
  exit 0
fi

# Sub-Branch D2: Generating map FOR the newly created HapHiC scaffolds
if [ "$mode" -eq 2 ]; then
  log "--- Sub-Branch D2: Post-HapHiC .hic Generation ---"
  genome="04.build/scaffolds.fa"

  conda activate tools_env
  seqkit seq -m 1000000 "$genome" > filtered_scaffolds_1Mb.fasta
  genome="filtered_scaffolds_1Mb.fasta"
  genomeIdx_Name="$(basename $genome).fai"
  suffix="${suffix}_hicScaff"
  conda deactivate

  # Repeat HiCUP logic for the new scaffolds
  if [ "$filter" -eq 2 ]; then
      log "Re-running HiCUP mapping on new scaffolds..."
      mv filtered_reads filtered_reads_FirstRun
      mkdir filtered_reads
      
      conda activate hicup_env
      bowtie2-build "$genome" "$index_dir/${suffix}"
      rm -f Digest_*.txt
      hicup_digester -genome "$suffix" -re1 ^GATC,DpnII "$genome" --outdir "$path"
      
      digest=$(ls "$path"/Digest_*.txt | head -n 1)
      sed -e "s|^Threads:.*|Threads: $threads|" \
          -e "s|^Index:.*|Index: $index_dir/$suffix|" \
          -e "s|^Digest:.*|Digest: $digest|" \
          pipeline_HiC.conf > pipeline_HiC_updated.conf
      mv pipeline_HiC_updated.conf pipeline_HiC.conf
      
      hicup --config pipeline_HiC.conf
      conda deactivate

      # Convert to .hic
      conda activate tools_env
      samtools faidx "$genome"
      cut -f1,2 "$genomeIdx_Name" > "${suffix}.chrom.sizes"
      awk '{print $1 "\t1\t" $2 "\t1\tW\t"$1 "\t1\t" $2 "\t+"}' "${suffix}.chrom.sizes" > "${suffix}.agp"

      samtools sort -n -o filtered_reads/"${sample_name}-pairsR1_2.hicup.namesort.bam" \
          filtered_reads/${sample_name}*.hicup.bam -@ $threads

      "$JUICER" pre filtered_reads/"${sample_name}-pairsR1_2.hicup.namesort.bam" \
        "${suffix}.agp" "$genomeIdx_Name" \
        | sort -k2,2d -k6,6d -T ./ --parallel="$threads" -S"${mem}G" \
        | awk 'NF' > alignments_sorted.txt

      java -Djava.awt.headless=true -jar -Xmx"${mem}G" "$JUICER_TOOLS" \
        pre alignments_sorted.txt "${suffix}_ContMap.hic" "${suffix}.chrom.sizes"
      
      conda deactivate
      log "Pipeline completed."
      exit 0
  fi
fi

log "Pipeline finished."