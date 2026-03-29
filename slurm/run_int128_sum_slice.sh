#!/usr/bin/env bash
#SBATCH --job-name=int128-slice
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=medium
#SBATCH --partition=B1

set -euo pipefail

: "${PROJECT_ROOT:?PROJECT_ROOT is required}"
: "${PAIR_FILE:?PAIR_FILE is required}"
: "${SYMMETRY:?SYMMETRY is required}"

task_index="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
pair_line="$(sed -n "$((task_index + 1))p" "$PAIR_FILE")"

if [[ -z "$pair_line" ]]; then
    echo "No pair found for task index $task_index in $PAIR_FILE" >&2
    exit 1
fi

read -r n1 n2 <<< "$pair_line"

echo "Running $SYMMETRY slice with n1=$n1 n2=$n2"
cd "$PROJECT_ROOT"
"${JULIA_BIN:-julia}" --project=. test/int128_sum_slice_driver.jl "$SYMMETRY" "$n1" "$n2"
