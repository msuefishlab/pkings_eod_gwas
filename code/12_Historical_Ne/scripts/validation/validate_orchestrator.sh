#!/bin/bash
#SBATCH --account=general
#SBATCH --partition=general-long
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00

# smc++ recovery of simulated bottlenecks (msprime 1.4.4): 7 scenarios x 25 chromosomes,
# vcf2smc (5 DIs) and smc++ estimate as in the real run, then analyze_results.py.
# Submit from the repository root after `source config.env`:
#   sbatch -o output_data/slurm_logs/12_Historical_Ne/validation_%j.log \
#     code/12_Historical_Ne/scripts/validation/validate_orchestrator.sh
# MSPPY = a Python with msprime, numpy and matplotlib (default: python3 on PATH).
set -euo pipefail

REPO=${SLURM_SUBMIT_DIR:-$(git rev-parse --show-toplevel)}
cd "$REPO" && source config.env
WORK=${SMCPP_WORK:-${scratch_dir}/smcpp}
MSPPY=${MSPPY:-python3}
SMCSIF=$REPO/images/smcpp.sif
SCRIPTS=$REPO/code/12_Historical_Ne/scripts/validation   # the two worker files below are (re)written here
VALDIR=$WORK/validation
RESDIR=$WORK/results/validation
NJOBS_SIM=8
NJOBS_SMC=8
NJOBS_EST=7

echo "======================================"
echo "smc++ validation orchestrator"
echo "Start: $(date)"
echo "Host:  $(hostname)"
echo "======================================"

# HTSlib for bgzip + tabix
module load HTSlib/1.18-GCC-12.3.0 2>/dev/null || module load HTSlib 2>/dev/null \
    || { echo "ERROR: HTSlib module not found"; exit 1; }
which bgzip && echo "bgzip OK: $(bgzip --version 2>&1 | head -1)"
which tabix && echo "tabix OK"

# matplotlib
$MSPPY -c "import matplotlib" 2>/dev/null \
    || $MSPPY -m pip install --quiet matplotlib 2>&1

# directories
for SC in null geneflow Tb500 Tb1000 Tb2000 Tb5000 Tb10000; do
    mkdir -p $VALDIR/sim_$SC
    mkdir -p $VALDIR/estimates/$SC
done
mkdir -p $RESDIR

# ── write worker scripts to known paths ─────────────────────────────────────

# vcf2smc worker
cat > $SCRIPTS/w_vcf2smc.sh << 'WORKER_EOF'
#!/bin/bash
# Usage: bash w_vcf2smc.sh VCF_GZ OUTFILE CHR DI_SAMPLE SMCSIF SAMPLE_LIST
VCF=$1; OUTF=$2; CHR=$3; DI=$4; SMCSIF=$5; SAMPLE_LIST=$6
if [ -f "$OUTF" ]; then
    echo "skip $CHR $DI"
    exit 0
fi
/usr/bin/singularity exec "$SMCSIF" smc++ vcf2smc \
    --cores 1 \
    -d "$DI" "$DI" \
    "$VCF" "$OUTF" "$CHR" \
    "POP:$SAMPLE_LIST" 2>&1 | tail -3
RC=$?
if [ $RC -ne 0 ]; then
    echo "ERROR ($RC) in vcf2smc for $CHR $DI" >&2
    exit 0
fi
echo "done $CHR $DI"
WORKER_EOF
chmod +x $SCRIPTS/w_vcf2smc.sh

# estimate worker
cat > $SCRIPTS/w_estimate.sh << 'WORKER_EOF'
#!/bin/bash
# Usage: bash w_estimate.sh SC VALDIR SMCSIF
SC=$1; VALDIR=$2; SMCSIF=$3
SIMDIR="$VALDIR/sim_$SC"
ESTDIR="$VALDIR/estimates/$SC"
mkdir -p "$ESTDIR"
N_SMC=$(ls "$SIMDIR"/*.smc.gz 2>/dev/null | wc -l || echo 0)
if [ "$N_SMC" -eq 0 ]; then
    echo "ERROR: no smc.gz files for $SC" >&2
    exit 0
fi
echo "  $SC: estimating on $N_SMC smc.gz files at $(date)..."
/usr/bin/singularity run "$SMCSIF" estimate \
    --polarization-error 0.5 \
    --cores 4 \
    -o "$ESTDIR/" \
    1e-8 \
    "$SIMDIR"/*.smc.gz 2>&1 | tail -10
echo "  Done: $SC  $(date)"
if [ -f "$ESTDIR/model.final.json" ]; then
    echo "  model.final.json: OK"
else
    echo "  WARNING: model.final.json not found!"
fi
WORKER_EOF
chmod +x $SCRIPTS/w_estimate.sh

# ── STEP 1: simulate ────────────────────────────────────────────────────────
echo ""
echo "=== STEP 1: msprime simulations ==="
echo "Start: $(date)"

SIM_TASKS=$VALDIR/sim_tasks.txt
> $SIM_TASKS

CHR_LENS=(
    chr1:77107170 chr2:45065832 chr3:44126694 chr4:49774254 chr5:42734709
    chr6:41647949 chr7:36680006 chr8:40861117 chr9:37833278 chr10:35833763
    chr11:30920909 chr12:33386478 chr13:33314801 chr14:30390776 chr15:30147657
    chr16:33751348 chr17:27084498 chr18:28897476 chr19:31461873 chr20:24061261
    chr21:21860095 chr22:22622874 chr23:25458968 chr24:28052655 chr25:35261244
)

SEED_BASE=42
for SC in null geneflow Tb500 Tb1000 Tb2000 Tb5000 Tb10000; do
    case $SC in
        null)      SCEN="null";       TB="NA" ;;
        geneflow)  SCEN="geneflow";   TB="NA" ;;
        Tb*)       SCEN="bottleneck"; TB=${SC#Tb} ;;
    esac
    for CLPAIR in "${CHR_LENS[@]}"; do
        CHR=${CLPAIR%%:*}
        LEN=${CLPAIR##*:}
        SEED=$(( SEED_BASE + $(echo "$SC$CHR" | cksum | awk '{print $1 % 9999}') ))
        echo "$SC $SCEN $TB $CHR $LEN $SEED" >> $SIM_TASKS
    done
done

echo "Total simulation tasks: $(grep -c . $SIM_TASKS)"

export WORK MSPPY SCRIPTS
grep . $SIM_TASKS | xargs -L 1 -P $NJOBS_SIM bash -c '
    SC=$1; SCEN=$2; TB=$3; CHR=$4; LEN=$5; SEED=$6
    OUTDIR="$WORK/validation/sim_$SC"
    if [ "$TB" = "NA" ]; then
        "$MSPPY" "$SCRIPTS/simulate_chrom.py" \
            --scenario "$SCEN" --chrom "$CHR" --length "$LEN" \
            --outdir "$OUTDIR" --seed "$SEED" 2>&1
    else
        "$MSPPY" "$SCRIPTS/simulate_chrom.py" \
            --scenario "$SCEN" --T_b "$TB" --chrom "$CHR" --length "$LEN" \
            --outdir "$OUTDIR" --seed "$SEED" 2>&1
    fi
' _ || echo "WARNING: some simulations returned non-zero (continuing)"

echo "Simulations complete: $(date)"
echo "VCF counts per scenario:"
for SC in null geneflow Tb500 Tb1000 Tb2000 Tb5000 Tb10000; do
    N=$(ls $VALDIR/sim_$SC/*.vcf.gz 2>/dev/null | wc -l || echo 0)
    echo "  $SC: $N VCFs"
done

# ── STEP 2: vcf2smc ─────────────────────────────────────────────────────────
echo ""
echo "=== STEP 2: vcf2smc ==="
echo "Start: $(date)"

SMC_TASKS=$VALDIR/vcf2smc_tasks.txt
> $SMC_TASKS

SAMPLE_LIST=$(python3 -c "print(','.join('tsk_'+str(i) for i in range(20)))")

for SC in null geneflow Tb500 Tb1000 Tb2000 Tb5000 Tb10000; do
    SIMDIR=$VALDIR/sim_$SC
    for CLPAIR in "${CHR_LENS[@]}"; do
        CHR=${CLPAIR%%:*}
        VCF_GZ=$SIMDIR/${CHR}.vcf.gz
        if [ ! -f "$VCF_GZ" ]; then
            echo "  WARN: missing $VCF_GZ" >&2
            continue
        fi
        for DI in 0 1 2 3 4; do
            OUTF=$SIMDIR/${CHR}_di${DI}.smc.gz
            echo "$VCF_GZ $OUTF $CHR tsk_${DI} $SMCSIF $SAMPLE_LIST" >> $SMC_TASKS
        done
    done
done

echo "Total vcf2smc tasks: $(grep -c . $SMC_TASKS || echo 0)"

grep . $SMC_TASKS | xargs -L 1 -P $NJOBS_SMC bash "$SCRIPTS/w_vcf2smc.sh" \
    || echo "WARNING: some vcf2smc tasks returned non-zero (continuing)"

echo "vcf2smc complete: $(date)"
echo "smc.gz counts per scenario:"
for SC in null geneflow Tb500 Tb1000 Tb2000 Tb5000 Tb10000; do
    N=$(ls $VALDIR/sim_$SC/*.smc.gz 2>/dev/null | wc -l || echo 0)
    echo "  $SC: $N smc.gz files"
done

# ── STEP 3: smc++ estimate ───────────────────────────────────────────────────
echo ""
echo "=== STEP 3: smc++ estimate ==="
echo "Start: $(date)"

for SC in null geneflow Tb500 Tb1000 Tb2000 Tb5000 Tb10000; do
    bash "$SCRIPTS/w_estimate.sh" "$SC" "$VALDIR" "$SMCSIF" &
done
wait
echo "All estimate jobs finished"

echo "smc++ estimate complete: $(date)"
echo "Estimate outputs:"
for SC in null geneflow Tb500 Tb1000 Tb2000 Tb5000 Tb10000; do
    if [ -f $VALDIR/estimates/$SC/model.final.json ]; then
        echo "  $SC: model.final.json OK"
    else
        echo "  $SC: MISSING"
    fi
done

# ── STEP 4: analyze ──────────────────────────────────────────────────────────
echo ""
echo "=== STEP 4: Analyze results ==="
echo "Start: $(date)"

$MSPPY $SCRIPTS/analyze_results.py \
    --workdir $VALDIR \
    --outdir  $RESDIR

echo ""
echo "======================================"
echo "Validation DONE: $(date)"
echo "Results: $RESDIR"
ls -lh $RESDIR/ 2>/dev/null || echo "Results dir empty"
echo "======================================"
exit 0
