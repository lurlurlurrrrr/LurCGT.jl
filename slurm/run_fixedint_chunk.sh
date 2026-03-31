#!/usr/bin/env bash
#SBATCH --job-name=fixedint-chunk
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=medium
#SBATCH --partition=B2

set -euo pipefail

: "${PROJECT_ROOT:?PROJECT_ROOT is required}"
: "${PAIR_FILE:?PAIR_FILE is required}"
: "${NUMTYPE:?NUMTYPE is required}"
: "${SYMMETRY:?SYMMETRY is required}"

task_index="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
pair_line="$(sed -n "$((task_index + 1))p" "$PAIR_FILE")"

if [[ -z "$pair_line" ]]; then
    echo "No pair found for task index $task_index in $PAIR_FILE" >&2
    exit 1
fi

read -r chunk1 chunk2 <<< "$pair_line"

echo "Running $NUMTYPE $SYMMETRY chunk $chunk1,$chunk2"
cd "$PROJECT_ROOT"
"${JULIA_BIN:-julia}" --project=. test/fixedint_chunk_driver.jl "$NUMTYPE" "$SYMMETRY" "$D1MIN" "$D1MAX" "$D2MIN" "$D2MAX" "$M1" "$M2" "$chunk1" "$chunk2"
