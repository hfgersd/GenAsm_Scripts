#!/bin/bash
#
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J RepProfile
#SBATCH -o RepProfile_%j.out
#SBATCH -e RepProfile_%j.err
#SBATCH --time=14-00:00:00      # Combined time for all steps
#SBATCH --mem=400G              # Max memory required (by EDTA)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80       # Max CPUs required (by TRF)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=juan.rosascabrera@kaust.edu.sa

# -------------------------------------------------------------------
# 1) ENABLE FAIL-FAST BEHAVIOR
# -------------------------------------------------------------------
set -e
set -o pipefail

# -------------------------------------------------------------------
# 2) UNIVERSAL VARIABLES DEFINITION
# -------------------------------------------------------------------

# -- Input / Output Parameters --
genome="/path/to/your/genome.fasta"  # Absolute path to the input genome
path="/path/to/working/directory"    # Working directory
suffix="SpeciesName_v1"              # Base prefix for naming output files

# -- Workflow Toggles --
# 1 = Run EDTA annotation step | 2 = Skip annotation (run TE only)
edta_annotation=1 

# -- Thread & Processing Allocation --
# Note: Total SLURM cpus-per-task must be >= the highest value here.
EDTA_threads=70
rm_threads=50
trf_subfiles=80  # Number of splits/threads for parallel TRF

# -- Databases and Custom Scripts --
rDNA_DB="/ibex/project/c2016/rosascjd/rDNA_DB/rDNA_datepalm.fna"
rDNAfinder="/ibex/scratch/projects/c2016/celiim/tools/rDNA_finder.R"
Tel_finder="/ibex/scratch/projects/c2016/celiim/tools/Telomere_finder.R"

# -- Conda Initialization Paths --
conda_sh="/ibex/user/rosascjd/miniforge/etc/profile.d/conda.sh"
mamba_sh="/ibex/user/rosascjd/miniforge/etc/profile.d/mamba.sh"

# -------------------------------------------------------------------
# 3) LOGGING FUNCTION
# -------------------------------------------------------------------
log_file="progress_RepProfile.txt"

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$log_file"
}

#####################################################################
#################### MAIN PIPELINE STARTS HERE ######################
#####################################################################

cd "$path" || { echo "Error: Working directory not found"; exit 1; }
log "PIPELINE STARTED"
log "Working directory set to $path"

# ===================================================================
# STEP 1: EDTA (Extensive de-novo TE Annotator)
# ===================================================================
log "--- STEP 1: Starting EDTA ---"
module load edta

if [ "$edta_annotation" -eq 1 ]; then
  log "Running EDTA WITH annotation steps..."
  EDTA.pl --genome "$genome" --species others --anno 1 --t $EDTA_threads --force 1 --overwrite 0
elif [ "$edta_annotation" -eq 2 ]; then
  log "Running EDTA WITHOUT annotation steps..."
  EDTA.pl --genome "$genome" --species others --anno 0 --t $EDTA_threads --force 1 --overwrite 0 --sensitive 0
fi
log "EDTA execution finished."

# ===================================================================
# STEP 2: HELITRON REMOVAL & rDNA IDENTIFICATION
# ===================================================================
log "--- STEP 2: Filtering Helitrons & Identifying rDNA ---"
source "$conda_sh"
source "$mamba_sh"
conda activate tools_env

log "Creating FASTA index..."
samtools faidx ./*mod.EDTA.TElib.fa

mkdir -p helitron_removing

log "Filtering Helitron IDs..."
grep "Hel" ./*mod.EDTA.TElib.fa.fai | cut -f1 > helitron_removing/TElib_Helitron_only.txt
grep -v "Hel" ./*mod.EDTA.TElib.fa.fai | cut -f1 > helitron_removing/TElib_No_Helitron.txt

log "Extracting sequences..."
seqtk subseq ./*mod.EDTA.TElib.fa helitron_removing/TElib_Helitron_only.txt \
  > helitron_removing/"${suffix}".mod.EDTA.TElib_Helitron_only.fa
seqtk subseq ./*mod.EDTA.TElib.fa helitron_removing/TElib_No_Helitron.txt \
  > helitron_removing/"${suffix}".mod.EDTA.TElib_No_Helitron.fa

log "Running rDNA identification script..."
Rscript "$rDNAfinder" --ref.rDNA="$rDNA_DB" --genome="$genome" --out=helitron_removing/"${suffix}"_rDNA.fasta

log "Concatenating final TE library..."
cat helitron_removing/"${suffix}".mod.EDTA.TElib_No_Helitron.fa \
    helitron_removing/"${suffix}"_rDNA.fasta \
  > helitron_removing/"${suffix}".mod.EDTA.TElib_No_Helitron_plus_rDNA.fa

# Save the final library path for downstream steps
final_lib_file=$(ls helitron_removing/"${suffix}".mod.EDTA.TElib_No_Helitron_plus_rDNA.fa)

# ===================================================================
# STEP 3: 100KB FASTA FILTERING
# ===================================================================
log "--- STEP 3: Filtering Contigs/Scaffolds < 100kb ---"
samtools faidx "$genome"

# Check if any sequence is shorter than 100,000 bp
if awk '$2 <= 100000 {exit 1}' "$genome".fai; then
    log "No sequences below 100kb detected. Copying original genome."
    cp "$genome" "${suffix}.fasta"
    touch NoFilt.txt
else
    log "Sequences below 100kb detected. Filtering with seqkit..."
    suffix="${suffix}.100kb" # Update suffix to reflect filtered state dynamically
    seqkit seq -m 100000 "$genome" > "${suffix}.fasta"
fi
conda deactivate

# ===================================================================
# STEP 4: REPEATMASKER
# ===================================================================
log "--- STEP 4: Starting RepeatMasker ---"
module load repeatmasker

log "Running RepeatMasker with $rm_threads cores..."
RepeatMasker -qq -no_is --pa $rm_threads \
  -lib "$final_lib_file" \
  "${suffix}.fasta" \
  -dir "$path" \
  -gff -xsmall 

log "RepeatMasker finished."

# ===================================================================
# STEP 5: TANDEM REPEAT FINDER (TRF)
# ===================================================================
log "--- STEP 5: Starting TRF (Parallelized) ---"
module load trf seqkit

mkdir -p subseqs
seq_count=$(grep -c "^>" "${suffix}.fasta")
log "Found $seq_count sequences in ${suffix}.fasta"

# Determine how many subfiles to create
if [ "$seq_count" -lt "$trf_subfiles" ]; then
    nfiles=$seq_count
else
    nfiles=$trf_subfiles
fi

log "Splitting FASTA into $nfiles subfiles for TRF..."
seqkit split --by-part "$nfiles" "${suffix}.fasta" --out-dir subseqs

# Run TRF in parallel (backgrounding each process)
for file in subseqs/*.fasta; do
    (
      fname=$(basename "$file")
      trf "subseqs/$fname" 2 3 5 80 10 30 2000 -l 6 -dat -h > /dev/null 2>&1
    ) &
done

log "Waiting for all TRF parallel jobs to complete..."
wait 
log "All parallel TRF jobs finished. Merging .dat files..."

# Move generated .dat files into the subseqs folder to keep it clean
mv *2.3.5.80.10.30.2000.dat subseqs/ 2>/dev/null || true

# Merge TRF outputs
merged_output="${suffix}.fasta.2.3.5.80.10.30.2000.dat"
dat_files=(subseqs/*.fasta.2.3.5.80.10.30.2000.dat)

if [ ${#dat_files[@]} -eq 0 ]; then
    log "ERROR: No TRF output files found!"
    exit 1
fi

head -n 5 "${dat_files[0]}" > "$merged_output"
for f in "${dat_files[@]}"; do
    tail -n +6 "$f" >> "$merged_output"
done
log "Merged TRF output saved to $merged_output."

# ===================================================================
# STEP 6: TELOMERE FINDER
# ===================================================================
log "--- STEP 6: Starting Telomere Finder ---"
source "$conda_sh"
conda activate genomescope2
module load blast bedtools trf samtools ragtag repeatmasker

log "Running Telomere_finder.R..."
Rscript "$Tel_finder" \
  --fasta="${suffix}.fasta" \
  --chr.filt=1Mb \
  --sort.size=TRUE \
  --tel.seq=aaaccct \
  --lib="$final_lib_file"

conda deactivate

# ===================================================================
# STEP 7: CLEANUP
# ===================================================================
log "--- STEP 7: Cleaning up temporary files ---"

rm -rf subseqs

# If the FASTA was simply copied (no 100kb filter), delete the copy to save space
if [[ -f "NoFilt.txt" ]]; then
  log "NoFilt.txt found. Removing unfiltered FASTA copy (${suffix}.fasta) to save space."
  rm -f "${suffix}.fasta"
  rm -f NoFilt.txt
fi

log "PIPELINE COMPLETED SUCCESSFULLY"