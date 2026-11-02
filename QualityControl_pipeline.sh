#!/bin/bash
#
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J Qua.Run
#SBATCH -o Quality_%j.out
#SBATCH -e Quality_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=250G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=

# ===================================================================
# PREREQUISITES
# 1. Conda environment containing 'quast' and 'busco'.
# 2. Local BUSCO lineage database downloaded.
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
conda_env="" # e.g., quality_env

# -- Processing Parameters --
threads=32
out_dir=""   # Output directory for results

# -- Run Mode Logic --
# 1 = Only run QUAST | 2 = Only run BUSCO | 3 = Run BOTH
mode=2

# -- QUAST Parameters --
ctg_thr=""   # Contig thresholds (comma separated, e.g. "0,1000,10000"). Leave blank for default.

# -- BUSCO Parameters --
buscoDB=""   # Path to the BUSCO lineage database directory

# -- Assemblies Array --
# Add the absolute paths to the FASTA assemblies you want to evaluate
assemblies=(
    ""
    ""
    ""
)

# -- Suffixes Array --
# Provide a matching label for each assembly path above
suffixes=(
    ""
    ""
    ""
)

# Number of assemblies to process (Should match the arrays above)
number_assemblies=3

# -------------------------------------------------------------------
# 3) LOGGING FUNCTION
# -------------------------------------------------------------------
log_file="progress_quality.txt"

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$log_file"
}

#####################################################################
#################### MAIN PIPELINE STARTS HERE ######################
#####################################################################

log "Quality assessment pipeline started"

source "$conda_sh"
source "$mamba_sh"
conda activate "$conda_env"

cd "$out_dir" || { echo "Directory not found"; exit 1; }
log "Working directory set to $out_dir"

# ===================================================================
# STEP 1: QUAST ANALYSIS
# ===================================================================
if [ "$mode" -eq 1 ] || [ "$mode" -eq 3 ]; then
    log "--- Starting QUAST Analysis ---"

    # Join suffixes and assemblies arrays into properly formatted strings
    label_string=$(IFS=, ; echo "${suffixes[*]}")
    assembly_string="${assemblies[*]}"

    quast.py \
        --labels "$label_string" \
        $assembly_string \
        --threads "$threads" \
        -o "$out_dir/quast_all" \
        --contig-thresholds "$ctg_thr"

    log "QUAST analysis completed."
fi

# ===================================================================
# STEP 2: BUSCO ANALYSIS
# ===================================================================
# Note: Original script had a condition limiting BUSCO if mode=3 and number_assemblies=1.
# I have generalized this so BUSCO runs as long as mode allows it.
if [ "$mode" -eq 2 ] || [ "$mode" -eq 3 ]; then
    log "--- Starting BUSCO Analysis ---"

    for ((i=0; i<number_assemblies; i++)); do
        assembly="${assemblies[$i]}"
        suffix="${suffixes[$i]}"

        log "Running BUSCO on assembly: $suffix"
        busco -i "$assembly" -c $threads -l "$buscoDB" -o "busco_${suffix}_tmp" -m genome -f
        
        # Rename standard output dir to custom suffix
        mv "busco_${suffix}_tmp" "busco_${suffix}"
        log "BUSCO completed for: $suffix"
    done
fi

log "Pipeline finished successfully."