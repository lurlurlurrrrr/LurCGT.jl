#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 8 ]]; then
    echo "Usage: $0 <NUMTYPE> <SYMMETRY> <D1MIN> <D1MAX> <D2MIN> <D2MAX> <M1> <M2>" >&2
    echo "Example: $0 Int128 SU3 1 20 1 20 4 4" >&2
    exit 1
fi

numtype="$1"
symmetry="$2"
d1min="$3"
d1max="$4"
d2min="$5"
d2max="$6"
m1="$7"
m2="$8"

for value in "$d1min" "$d1max" "$d2min" "$d2max" "$m1" "$m2"; do
    [[ "$value" =~ ^[0-9]+$ ]] || {
        echo "Range and chunk arguments must be nonnegative integers." >&2
        exit 1
    }
done

if (( d1min > d1max || d2min > d2max )); then
    echo "Dimension minimum must be <= maximum." >&2
    exit 1
fi

if (( m1 == 0 || m2 == 0 )); then
    echo "Chunk counts must be positive." >&2
    exit 1
fi

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
generated_dir="$project_root/slurm/generated"
logs_dir="$project_root/slurm/logs"
mkdir -p "$generated_dir" "$logs_dir"

chunk_start() {
    local min="$1" max="$2" count="$3" index="$4"
    local len=$((max - min + 1))
    echo $((min + ((index - 1) * len) / count))
}

chunk_end() {
    local min="$1" max="$2" count="$3" index="$4"
    local len=$((max - min + 1))
    echo $((min + (index * len) / count - 1))
}

timestamp="$(date +%Y%m%d_%H%M%S)"
pair_file="$generated_dir/${numtype}_${symmetry}_${d1min}_${d1max}_${d2min}_${d2max}_${m1}_${m2}_${timestamp}.txt"
: > "$pair_file"

pair_count=0
for ((i = 1; i <= m1; i++)); do
    start1="$(chunk_start "$d1min" "$d1max" "$m1" "$i")"
    end1="$(chunk_end "$d1min" "$d1max" "$m1" "$i")"
    (( start1 <= end1 )) || continue
    for ((j = 1; j <= m2; j++)); do
        start2="$(chunk_start "$d2min" "$d2max" "$m2" "$j")"
        end2="$(chunk_end "$d2min" "$d2max" "$m2" "$j")"
        (( start2 <= end2 )) || continue
        (( start1 <= end2 )) || continue
        printf "%d %d\n" "$i" "$j" >> "$pair_file"
        pair_count=$((pair_count + 1))
    done
done

if (( pair_count == 0 )); then
    echo "No canonical chunk jobs to submit." >&2
    rm -f "$pair_file"
    exit 1
fi

array_spec="0-$((pair_count - 1))"
if [[ -n "${MAX_PARALLEL:-}" ]]; then
    array_spec="${array_spec}%${MAX_PARALLEL}"
fi

echo "Chunk list written to $pair_file"
echo "Submitting $pair_count array tasks for $numtype $symmetry"

job_id="$(
    sbatch --parsable \
        --array="$array_spec" \
        --export=ALL,PROJECT_ROOT="$project_root",PAIR_FILE="$pair_file",NUMTYPE="$numtype",SYMMETRY="$symmetry",D1MIN="$d1min",D1MAX="$d1max",D2MIN="$d2min",D2MAX="$d2max",M1="$m1",M2="$m2",JULIA_BIN="${JULIA_BIN:-julia}" \
        --output="$logs_dir/fixedint_chunk_%A_%a.out" \
        --error="$logs_dir/fixedint_chunk_%A_%a.err" \
        "$project_root/slurm/run_fixedint_chunk.sh"
)"

echo "Submitted job array $job_id"
