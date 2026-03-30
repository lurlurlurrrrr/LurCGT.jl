# Int128 Chunked CG Test Design

## Goal

Add a new Int128 test workflow to `LurCGT` that:
- builds a per-symmetry saved catalog of irreps sorted by dimension
- excludes irreps whose `Int128` and `BigInt` realizations differ
- runs CG generation tests in small remote-friendly chunks over dimension ranges
- saves chunk results for later aggregation and plotting
- keeps the implementation style close to the existing `src/inttest.jl` helpers

## Scope

In scope:
- Extend `src/inttest.jl` with reusable catalog, chunk-runner, and plot-preparation helpers.
- Add thin CLI drivers in `test/` for catalog updates, chunk execution, and plotting.
- Save one catalog file per symmetry under the repo so the data can be updated and reused.
- Keep chunk execution canonical so `generate_every_CGT` only receives `q1 <= q2`.
- Support remote execution by dividing requested dimension windows into `m1 x m2` chunk jobs.

Out of scope for this pass:
- Replacing the existing frontier and sum-slice helpers.
- Changing `generate_every_CGT` ordering requirements.
- Introducing a new external database format for these test results.
- Rich visualization beyond a pass/fail dimension plot.

## Existing Integration Points

The current Int128-related helpers already live in `src/inttest.jl`, including:
- irrep comparison between `BigInt` and `Int128`
- `run_int128_cgt_sum_slice`
- helper logic for enumerating qlabels and testing CG generation

The current remote-friendly execution style already exists in:
- `test/int128_sum_slice_driver.jl`
- `slurm/submit_int128_sum_slices.sh`
- `slurm/run_int128_sum_slice.sh`

The new workflow should follow the same pattern:
- keep most logic in `src/inttest.jl`
- keep `test/` entrypoints thin
- keep the server-facing unit of work small and resumable

## File Layout

### Code

Modify:
- `src/inttest.jl`
- `src/LurCGT.jl`
- `Project.toml`
- `test/int128_sum_slice_smoketest.jl`

Create:
- `test/int128_irep_catalog_driver.jl`
- `test/int128_chunk_driver.jl`
- `test/int128_plot_driver.jl`
- `test/int128_chunk_smoketest.jl`
- `slurm/submit_int128_chunks.sh`
- `slurm/run_int128_chunk.sh`

### Saved Data

Create per-symmetry directories under:
- `test/int128_data/<SYMMETRY>/`

Expected contents:
- `test/int128_data/<SYMMETRY>/catalog.jls`
- `test/int128_data/<SYMMETRY>/chunks/`
- `test/int128_data/<SYMMETRY>/plots/`

This keeps the generated artifacts versionable, easy to inspect, and easy to reuse on the remote server.

## Catalog Data Model

Each symmetry gets one saved catalog file:
- `test/int128_data/<SYMMETRY>/catalog.jls`

The file should store both accepted irreps and scan state:

```julia
(
    symmetry = "SU3",
    version = 1,
    accepted = [
        (qlabel=(0, 0), dim=1),
        (qlabel=(1, 0), dim=3),
        (qlabel=(0, 1), dim=3),
        (qlabel=(2, 0), dim=6),
    ],
    scanned = Dict(
        (0, 0) => (status=:accepted, dim=1, reason="match"),
        (1, 0) => (status=:accepted, dim=3, reason="match"),
        (0, 1) => (status=:accepted, dim=3, reason="match"),
        (3, 0) => (status=:rejected, dim=10, reason="Int128 and BigInt irreps differ"),
    ),
)
```

Why both parts are needed:
- `accepted` is the list used by CG testing and plotting.
- `scanned` allows updates to resume without redoing work.
- `scanned` also prevents rejected irreps from blocking exploration of larger irreps.

The `accepted` list must stay sorted by:
1. irrep dimension
2. `qlabel`

## Catalog Update Rules

The update path should remain close to the current `inttest.jl` exploration style:
- start from known scanned qlabels
- grow candidates by adding fundamental qlabels
- inspect candidates in ascending dimension order

Recommended API:

```julia
update_int128_irrep_catalog(::Type{S}; maxdim::Int, maxcount=nothing, save=true)
```

Behavior for each candidate qlabel:
1. Build the `BigInt` irrep.
2. Build the `Int128` irrep.
3. If the two irreps match except for `size_byte`, record `:accepted`.
4. If they differ, or `Int128` throws an expected overflow-like failure, record `:rejected`.
5. Mark the qlabel as scanned in either case.

The main stopping rule should be `maxdim`, because chunk jobs are dimension-range based. `maxcount` can remain optional for debugging or bounded experiments.

## Canonical Pair Rule

`generate_every_CGT` only accepts ordered inputs with `q1 <= q2`, so the new workflow must enforce canonical pairs by dimension and qlabel:
- test only when `dim(q1) <= dim(q2)`
- if `dim(q1) == dim(q2)`, test only when `q1 <= q2`

This rule must be applied consistently in:
- in-memory pair enumeration
- chunk execution
- server-side chunk submission
- saved results

## Chunked Execution Model

The server-facing input should be:

```text
<SYMMETRY> <d1min> <d1max> <d2min> <d2max> <m1> <m2>
```

Interpretation:
- `[d1min, d1max]` is the first dimension window
- `[d2min, d2max]` is the second dimension window
- split the first window into `m1` contiguous chunks
- split the second window into `m2` contiguous chunks
- assign one chunk-pair to one process

Example:

```text
SU3 1 20 1 20 4 4
```

Possible chunk ranges:
- `1:5`
- `6:10`
- `11:15`
- `16:20`

A chunk job such as `(6:10, 11:15)` should:
1. load the saved catalog for `SU3`
2. collect accepted irreps with dimensions in `6:10` and `11:15`
3. enumerate candidate pairs
4. keep only canonical pairs satisfying the ordering rule above
5. call `generate_every_CGT(S, Int128, Int128, (q1, q2), nothing; assertlev=1, save=false)`
6. save one chunk result file

The mirrored chunk `(11:15, 6:10)` should not be submitted, because it can only produce non-canonical work.

## Chunk Result Format

Each chunk job should save a self-contained file under:
- `test/int128_data/<SYMMETRY>/chunks/<range1>__<range2>.jls`

Recommended payload:

```julia
(
    symmetry = "SU3",
    dim_range1 = (6, 10),
    dim_range2 = (11, 15),
    total_pairs = 8,
    passed_pairs = 7,
    failed_pairs = 1,
    results = [
        (q1=(2,0), dim1=6, q2=(2,1), dim2=15, status=:passed, reason="passed"),
        (q1=(0,2), dim1=6, q2=(2,1), dim2=15, status=:cgt_failed, reason="OverflowError"),
    ],
)
```

Only canonical pairs should appear in `results`.

Failure classes for this pass only need to distinguish:
- `:passed`
- `:cgt_failed`

Rejected irreps should never appear here because they are already filtered out during catalog construction.

## Drivers And Remote Execution

The CLI layer should stay thin:

- `test/int128_irep_catalog_driver.jl`
  - updates or builds `catalog.jls` for one symmetry
  - takes a symmetry name and a maximum dimension

- `test/int128_chunk_driver.jl`
  - runs one chunk job for one symmetry and one chunk-pair
  - takes the server-facing input shape described above plus the selected chunk indices

- `test/int128_plot_driver.jl`
  - loads saved chunk files for one symmetry and writes plot output

The SLURM layer should mirror the existing sum-slice pattern:
- one submit script computes the canonical chunk-pair list
- one run script executes exactly one chunk job from that list
- job arrays stay small, resumable, and easy to rerun

## Plotting

Plotting is a separate pass over saved chunk files.

For one symmetry:
1. load all saved chunk results
2. aggregate by `(dim(q1), dim(q2))`
3. mark each dimension cell:
   - green if every tested pair in that cell passed
   - red if any tested pair in that cell failed

This keeps the plot simple and safety-oriented while preserving full pair detail in the serialized chunk results.

The plot output should be written under:
- `test/int128_data/<SYMMETRY>/plots/`

The first implementation should use `Plots.jl` to write a simple static image under `plots/`. The saved result format must stay independent of the plotting library so the backend can change later without rewriting chunk data.

## Testing Strategy

Add focused tests that cover the reusable logic rather than only the drivers.

Required coverage:
1. Catalog update behavior on a tiny symmetry such as `SU{2}`.
2. Comparison logic that accepts identical irreps except for `size_byte`.
3. Rejection logic when `Int128` and `BigInt` irreps differ.
4. Chunk-splitting logic for dimension windows and `m1`, `m2`.
5. Canonical pair filtering:
   - keep pairs with `dim(q1) < dim(q2)`
   - keep equal-dimension pairs only when `q1 <= q2`
   - reject mirrored duplicates
6. Chunk result summaries for a tiny runnable example.
7. Thin-driver smoke tests for argument parsing and output shape where practical.

The implementation should prefer small deterministic smoke tests and avoid large remote-only workflows in the regular local test suite.

## Risks And Mitigations

Risk: large dimension windows may still produce very uneven chunk workloads.

Mitigation: start with simple contiguous dimension chunks because they match the requested API and are easy to reason about. If load balance becomes a problem later, add a separate manifest/chunking strategy without changing the catalog format.

Risk: storing only accepted irreps might make updates non-resumable.

Mitigation: persist both `accepted` and `scanned` state in the catalog file.

Risk: duplicate work may reappear if chunk submission and in-process filtering disagree.

Mitigation: apply the same canonical ordering rule in both the job generator and the chunk runner.

Risk: plotting may become coupled to a specific graphics package too early.

Mitigation: keep serialized chunk results backend-agnostic and make plotting a separate final pass.

## Implementation Summary

This should be an extension of the current Int128 tooling, not a parallel system:
- keep the core logic in `src/inttest.jl`
- expose the new helpers from `src/LurCGT.jl`
- add thin drivers under `test/`
- save per-symmetry catalog and chunk files under `test/int128_data/`
- use chunked dimension windows for remote execution
- aggregate saved chunk results into a dimension-based red/green plot
