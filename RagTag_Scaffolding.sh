#!/bin/bash 
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J ragtag.Scaf
#SBATCH -o ragtag_%j.out
#SBATCH -e ragtag_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

# ===================================================================
# PREREQUISITES
# 1. Modules: 'ragtag', 'samtools', 'bedtools'
# 2. Conda environment containing R/GenomeScope2.
# 3. Custom GS-viewer R script.
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
conda_env=""    # e.g., genomescope2

# -- Hardware Allocation --
threads=64

# -- Input Parameters --
assembly=""     # Path to the assembly FASTA to be scaffolded
GSviewer=""     # Path to the GS-viewer.r script

# -- References Array --
# Absolute paths to the reference genomes
references=(
    ""
    ""
    ""
)

# -- Suffixes Array --
# Matching short names for the references used above
suffixes=(
    ""
    ""
    ""
)

# Total number of references provided above
num_references=3  

# -------------------------------------------------------------------
# 3) LOGGING FUNCTION
# -------------------------------------------------------------------
log_file="progress_ragtag.txt"

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$log_file"
}

#####################################################################
#################### MAIN PIPELINE STARTS HERE ######################
#####################################################################

log "RagTag Pipeline Started"

# Validate array lengths
if [ "${#references[@]}" -ne "$num_references" ] || [ "${#suffixes[@]}" -ne "$num_references" ]; then
  log "Error: Number of references or suffixes array does not match num_references."
  exit 1
fi

genome_name=$(basename "$assembly")
joined_suffixes=$(IFS=', '; echo "${suffixes[*]}")
log "Scaffolding $genome_name against: $joined_suffixes"

# ===================================================================
# STEP 1: RAGTAG SCAFFOLDING (PARALLEL)
# ===================================================================
log "--- Starting RagTag Scaffolding ---"
module load ragtag

# Distribute CPU cores evenly among the parallel reference jobs
threads_job=$(( threads / num_references ))

for ((i=0; i<num_references; i++)); do
    ref="${references[$i]}"
    suffix="${suffixes[$i]}"
    
    (
        log "Launching RagTag for reference: $suffix"
        ragtag.py scaffold "$ref" "$assembly" \
            -o "scaffold_${suffix}" -t $threads_job -a 0.01 \
            --remove-small
        log "Finished RagTag for reference: $suffix"
    ) &
done

log "Waiting for all parallel RagTag jobs to complete..."
wait
log "All RagTag jobs finished."

# ===================================================================
# STEP 2: GS-VIEWER ANALYSIS (PARALLEL)
# ===================================================================
log "--- Starting GS-Viewer Analysis ---"

source "$conda_sh"
conda activate "$conda_env"
module load samtools bedtools

for ((i=0; i<num_references; i++)); do
    ref="${references[$i]}"
    suffix="${suffixes[$i]}"

    (
        log "Launching GS-Viewer for: $suffix"
        Rscript "$GSviewer" --path="scaffold_${suffix}/" \
            --ref="$ref" --out.suffix="scaffold_${suffix}"
        log "Finished GS-Viewer for: $suffix"
    ) &
done

log "Waiting for all parallel GS-Viewer jobs to complete..."
wait
log "All GS-Viewer jobs completed."

conda deactivate
log "Pipeline finished successfully."