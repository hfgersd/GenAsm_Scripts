#!/bin/bash
#
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J kmer.prof
#SBATCH -o kmer_profiling_%j.out
#SBATCH -e kmer_profiling_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=800G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

# ===================================================================
# PREREQUISITES
# 1. Modules: 'jellyfish', 'smudgeplot'
# 2. Conda environment containing 'genomescope2' and 'smudgeplot'.
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
conda_env=""          # e.g., genomescope2 (should also contain smudgeplot dependencies)

# -- Hardware Allocation --
threads=32            # Number of threads for Jellyfish (must match --cpus-per-task)

# -- Input / Output Parameters --
reads=""              # Absolute path to the input raw reads (e.g., .fastq.gz)
out_dir=""            # Base output directory for all analyses
genome_prefix=""      # Name prefix for your sample (used in output filenames)

# -- Analysis Parameters --
jf_hash_size="5G"     # Initial hash size for Jellyfish (e.g., 5G, 10G)

# K-mer values to be evaluated (e.g., 21, 31, 41)
kmer_values=(
    ""
    ""
)

# Ploidy levels to evaluate in GenomeScope2 (e.g., 2, 4, 6)
ploidy_levels=(
    ""
    ""
)

# -------------------------------------------------------------------
# 3) LOGGING FUNCTION
# -------------------------------------------------------------------
log_file="progress_kmer_profiling.txt"

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$log_file"
}

#####################################################################
#################### MAIN PIPELINE STARTS HERE ######################
#####################################################################

log "K-mer Profiling Pipeline (Jellyfish + GenomeScope2 + Smudgeplot) Started"

# -------------------------------------------------------------------
# DIRECTORY SETUP
# -------------------------------------------------------------------
mkdir -p "$out_dir"
cd "$out_dir" || { echo "Directory not found"; exit 1; }
log "Working directory set to $out_dir"

mkdir -p jellyfish
mkdir -p genomescope
mkdir -p smudgeplot

# -------------------------------------------------------------------
# ENVIRONMENT SETUP
# -------------------------------------------------------------------
source "$conda_sh"
conda activate "$conda_env"
log "Conda environment $conda_env activated."

module load jellyfish smudgeplot
log "Modules 'jellyfish' and 'smudgeplot' loaded."

# ===================================================================
# ITERATE OVER SPECIFIED K-MER SIZES
# ===================================================================
for k in "${kmer_values[@]}"; do
    
    log "=========================================================="
    log "==== Processing k-mer size: $k ===="
    log "=========================================================="

    jf_db="jellyfish/${genome_prefix}_${k}mer.jf"
    histo_file="jellyfish/${genome_prefix}_${k}mer.histo"

    # ---------------------------------------------------------------
    # STEP 1: JELLYFISH COUNT & HISTOGRAM
    # ---------------------------------------------------------------
    log "--- STEP 1: Jellyfish Counting ---"
    
    # Run Jellyfish count
    jellyfish count -t $threads -C -m "$k" -s "$jf_hash_size" -o "$jf_db" "$reads"
    log "Jellyfish count completed."
    
    # Generate histogram
    log "Generating Jellyfish histogram..."
    jellyfish histo -o "$histo_file" "$jf_db"
    log "Jellyfish histogram generated."
    
    # ---------------------------------------------------------------
    # STEP 2: GENOMESCOPE2
    # ---------------------------------------------------------------
    log "--- STEP 2: GenomeScope2 ---"
    
    for p in "${ploidy_levels[@]}"; do
        log "Running GenomeScope2 for ploidy: $p..."
        gs_out="genomescope/${genome_prefix}_${k}mer_p${p}"
        
        genomescope2 -i "$histo_file" -o "$gs_out" -n "${genome_prefix}_${k}mer_p${p}" -k "$k" -p "$p"
        log "GenomeScope2 analysis completed for ploidy $p."
    done

    # ---------------------------------------------------------------
    # STEP 3: SMUDGEPLOT
    # ---------------------------------------------------------------
    log "--- STEP 3: Smudgeplot ---"

    # 3.1: Calculate Cutoffs
    log "Calculating Smudgeplot coverage cutoffs (L and U)..."
    L=$(smudgeplot.py cutoff "$histo_file" L)
    U=$(smudgeplot.py cutoff "$histo_file" U)
    log "Cutoffs defined: L=$L | U=$U"

    # 3.2: Dump k-mers inside the coverage boundaries
    jf_dump="smudgeplot/${genome_prefix}_${k}mer_L${L}_U${U}.fa"
    log "Dumping k-mers within coverage range [L, U] to $jf_dump..."
    jellyfish dump -c -L "$L" -U "$U" "$jf_db" > "$jf_dump"

    # 3.3: Convert to hetkmers
    sp_prefix="smudgeplot/${genome_prefix}_${k}mer_L${L}_U${U}"
    log "Converting to k-mer pairs (hetkmers)..."
    smudgeplot.py hetkmers -o "$sp_prefix" < "$jf_dump"

    # 3.4: Generate the Smudgeplot plot
    log "Plotting Smudgeplot..."
    smudgeplot.py plot \
      -o "${sp_prefix}_plot" \
      -t "${genome_prefix}_${k}mer" \
      -k "$k" \
      "${sp_prefix}_coverages.tsv"

    log "==== Finished all analyses for k=$k ===="
done

# -------------------------------------------------------------------
# CLEANUP
# -------------------------------------------------------------------
conda deactivate
log "Conda environment deactivated."
log "K-mer Profiling Pipeline Completed Successfully."