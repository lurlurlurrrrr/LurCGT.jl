#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "Usage: $0 <SYMMETRY> <N1MIN> <N1MAX> <N2MIN> <N2MAX>" >&2
    echo "Example: $0 SU3 4 8 1 5" >&2
    exit 1
fi

symmetry="$1"
n1min="$2"
n1max="$3"
n2min="$4"
n2max="$5"

for value in "$n1min" "$n1max" "$n2min" "$n2max"; do
    [[ "$value" =~ ^[0-9]+$ ]] || {
        echo "All range arguments must be nonnegative integers." >&2
        exit 1
    }
done

if (( n1min > n1max )); then
    echo "N1MIN must be <= N1MAX." >&2
    exit 1
fi

if (( n2min > n2max )); then
    echo "N2MIN must be <= N2MAX." >&2
    exit 1
fi

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
generated_dir="$project_root/slurm/generated"
logs_dir="$project_root/slurm/logs"
mkdir -p "$generated_dir" "$logs_dir"

timestamp="$(date +%Y%m%d_%H%M%S)"
pair_file="$generated_dir/${symmetry}_${n1min}_${n1max}_${n2min}_${n2max}_${timestamp}.txt"
pair_count=0

for (( n1 = n1min; n1 <= n1max; n1++ )); do
    upper_n2=$(( n1 < n2max ? n1 : n2max ))
    for (( n2 = n2min; n2 <= upper_n2; n2++ )); do
        printf "%d %d\n" "$n1" "$n2" >> "$pair_file"
        pair_count=$((pair_count + 1))
    done
done

if (( pair_count == 0 )); then
    echo "No (n1, n2) pairs satisfy n1 >= n2 inside the requested ranges." >&2
    rm -f "$pair_file"
    exit 1
fi

array_spec="0-$((pair_count - 1))"
if [[ -n "${MAX_PARALLEL:-}" ]]; then
    array_spec="${array_spec}%${MAX_PARALLEL}"
fi

echo "Pair list written to $pair_file"
echo "Submitting $pair_count array tasks for $symmetry"

job_id="$(
    sbatch --parsable \
        --array="$array_spec" \
        --export=ALL,PROJECT_ROOT="$project_root",PAIR_FILE="$pair_file",SYMMETRY="$symmetry",JULIA_BIN="${JULIA_BIN:-julia}" \
        --output="$logs_dir/int128_sum_slice_%A_%a.out" \
        --error="$logs_dir/int128_sum_slice_%A_%a.err" \
        "$project_root/slurm/run_int128_sum_slice.sh"
)"

echo "Submitted job array $job_id"
