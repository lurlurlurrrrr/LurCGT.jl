# Refactor plan: SQLite finalize and server-mode headnode merge

## Goal

Change LurCGT server-mode SQLite behavior so Slurm jobs use node-local
databases during computation but merge generated rows back into the global
database on the headnode or shared filesystem.

Target server-mode data flow:

1. When CGT data for a concrete non-Abelian symmetry `S` is first requested,
   copy only the headnode global database for that `S` into the node-local
   global directory.
2. During computation, read from the node-local global database and write new
   rows to a process-local node-local database.
3. At finalization, merge the process-local node-local database directly into
   the headnode global database.
4. After a successful merge, delete node-local local/global database files and
   SQLite sidecar files.

The headnode local database is no longer part of the server-mode workflow. It
should not be copied to the compute node and should not receive writes.

Target local-mode data flow:

1. During computation, read from the active global database and write new rows
   to the active process-local database.
2. At finalization, merge the process-local database into the active global
   database.
3. After a successful merge, delete the process-local database and SQLite
   sidecar files.

## Current behavior

`src/sqlite_io.jl` already separates source and node directories:

- `LURCGT_GLOBALDB_DIR`: source global database directory.
- `LURCGT_LOCALDB_DIR`: source local database directory.
- `LURCGT_GLOBALDB_DIR_NODE`: node-local global database directory.
- `LURCGT_LOCALDB_DIR_NODE`: node-local process-local database directory.

In `LURCGT_RUN_MODE=server`, `sqlite_global_dir()` resolves to
`LURCGT_GLOBALDB_DIR_NODE`, and `sqlite_local_dir()` resolves to
`LURCGT_LOCALDB_DIR_NODE`.

This is correct for computation, but it means `merge_all_to_global` currently
merges into the node-local global copy, not into the headnode global database.
That does not persist new rows after node cleanup unless the whole node-global
database is copied back, which would be unsafe when multiple Slurm jobs finish
around the same time.

## Desired semantics

Use these terms in the implementation and documentation:

- source global: `sqlite_global_source_dir()`, normally on headnode or NFS.
- active global: `sqlite_global_dir()`, node-local in server mode.
- active local: `sqlite_local_dir()`, process-local and node-local in server
  mode.
- merge target global: source global in server mode, active global in local
  mode.

Server-mode global copying must be lazy and per symmetry. Starting Julia,
loading LurCGT, or opening a local DB must not copy any global database. A
server-mode copy is allowed only when a load path for a concrete
`S<:NonabelianSymm` needs that symmetry's active global database. The copy
source is:

```julia
joinpath(sqlite_global_source_dir(env), "$(totxt(S)).db")
```

and the copy target is:

```julia
joinpath(sqlite_global_node_dir(env), "$(totxt(S)).db")
```

No other symmetry's global database should be copied as a side effect.

Missing global databases are not errors. If the merge target global DB for `S`
does not exist, create it with the normal schema and PRAGMAs before merging. In
server mode this means creating `joinpath(sqlite_global_source_dir(env),
"$(totxt(S)).db")`; in local mode this means creating
`joinpath(sqlite_global_dir(env), "$(totxt(S)).db")`.

In local mode, runtime behavior should remain unchanged, but finalization should
also merge and cleanup:

```julia
load:  global source/active DB, then local source/active DB
save:  local source/active DB
merge: local source/active DB -> global source/active DB
cleanup: remove local process DB
```

In server mode, behavior should become:

```julia
load:  node-global DB, then node-local process DB
save:  node-local process DB
merge: node-local process DB -> source global DB
cleanup: remove node-local process DB and node-global copy
```

Do not use the source local DB in server mode. If `LURCGT_LOCALDB_DIR` is set
while `LURCGT_RUN_MODE=server`, it should not affect runtime reads, writes, or
merge source selection.

## Non-goals

- Do not introduce a new database type or wrapper object.
- Do not change table schemas or serialized object formats.
- Do not copy the whole node-global database back to the source global
  location.
- Do not make every save write through to the source global database.
- Do not make the headnode local database meaningful in server mode.
- Do not add a background synchronization process.
- Do not delete global databases in local mode or server mode.
- Do not copy all global databases at package initialization or job startup.
- Do not copy global databases for symmetries that are not requested in the
  current Julia process.

## Phase 1: separate active DB paths from merge target paths

Add small path helpers in `src/sqlite_io.jl`:

```julia
sqlite_merge_global_dir(env=ENV) =
    sqlite_run_mode(env) == "server" ? sqlite_global_source_dir(env) : sqlite_global_dir(env)

sqlite_merge_lock_dir(env=ENV) = joinpath(sqlite_merge_global_dir(env), "locks")
```

Add a path helper for the merge target database:

```julia
sqlite_merge_global_db_path(S; env=ENV) =
    joinpath(sqlite_merge_global_dir(env), "$(totxt(S)).db")
```

Keep `sqlite_global_dir()` unchanged. It should continue to mean the active
global database used for reads during computation.

Update comments and docstrings so "global" is not ambiguous. Use "active
global" for the node-local copy and "source global" for the persistent
headnode/NFS database.

## Phase 2: merge into the source global DB in server mode

Change merge operations so the destination DB is not opened through
`get_sqlite_db(S, :global)` in server mode. That function intentionally opens
the active global DB, which is the node-local copy.

Add an internal helper:

```julia
function get_sqlite_merge_global_db(::Type{S}) where {S<:NonabelianSymm}
    path = sqlite_merge_global_db_path(S)
    # Open/create the source global DB, initialize schema and pragmas.
end
```

This helper should use the existing `SQLITE_DBS` registry with a distinct key,
for example:

```julia
"merge_global:$(path)"
```

This avoids aliasing the active node-global connection with the source-global
merge connection.

`merge_table_to_global` should use:

- source DB: `get_sqlite_db(S, :local)`
- destination DB: `get_sqlite_merge_global_db(S)`
- lock directory: `sqlite_merge_lock_dir()`

Continue to use `INSERT OR IGNORE`. Existing global rows should win. This is
the safest behavior when two jobs compute the same object.

If the merge target DB file does not exist, `get_sqlite_merge_global_db` should
create its parent directory, open a new SQLite DB, and call `_init_sqlite_db`.
The same behavior applies in local mode, where the merge target is the active
global DB. A missing global DB is an empty cache, not an exceptional state.

## Phase 3: make server-mode global copy lazy and per symmetry

Keep `_prepare_server_sqlite_db(S, :global, path)` on the global DB open path,
not on package initialization or any broad setup function. The helper should
copy at most the SQLite file set for `S`:

```julia
source_path = sqlite_db_source_path(S, :global; env, process_id)
target_path = sqlite_db_path(S, :global; env, process_id)
```

For `S == SU{2}`, this means only `SU2.db`, `SU2.db-wal`, and `SU2.db-shm` are
eligible for copying. It must not scan `sqlite_global_source_dir()` for other
`*.db` files.

Opening the active local DB for `S` must not force the active global DB for
`S` to open or copy. A job that only writes local rows and never loads existing
global CGT data should not pay the global copy cost until a global lookup is
actually attempted.

If the source global DB for `S` does not exist, create an empty active
node-global DB for `S` only when the global load path asks for it.

## Phase 4: remove source-local copy behavior in server mode

Currently `_prepare_server_sqlite_db(S, :local, path)` can copy from
`sqlite_local_source_dir()` into `sqlite_local_node_dir()`. Remove this copy
path for local databases in server mode.

Desired behavior:

- `location == :global`: copy source global to active node-global on first open
  if the source exists.
- `location == :local`: create a fresh node-local process DB if it does not
  exist.

This makes `LURCGT_LOCALDB_DIR` irrelevant in server mode and avoids accidental
reuse of stale local rows from the headnode.

Keep per-process local filenames based on `process_local_id()`.

## Phase 5: add explicit finalization

Add a user-facing helper:

```julia
finalize_sqlite!(S; tables=collect(SQLITE_TABLES), cleanup=true, verbose=1)
```

Semantics:

1. Merge active local rows into the merge target global DB.
   - In local mode, this is the active global DB.
   - In server mode, this is the source global DB on the headnode or shared
     filesystem.
2. Clear successfully merged rows from the active local DB, preserving the
   current `merge_all_to_global(...; clear_local_after=true)` behavior.
3. Close all SQLite DBs.
4. If `cleanup=true`, delete temporary runtime DB files.
   - In local mode, delete the current process-local DB for `S`.
   - In server mode, delete the current process-local node DB for `S` and the
     node-local global copy for `S`.

The cleanup step should remove base `.db`, `-wal`, and `-shm` files. It should
not delete any global merge target files.

If the merge throws, do not cleanup. Leaving local files behind is better than
losing generated rows.

Keep `finalize_server_sqlite!` only as an optional backward-compatible alias if
needed. New documentation and examples should use `finalize_sqlite!`.

Add internal cleanup helpers instead of broad directory deletion:

```julia
delete_current_local_sqlite_db(S; env=ENV, verbose=0)
delete_active_global_sqlite_copy(S; env=ENV, verbose=0)
```

`delete_current_local_sqlite_db` should delete only the DB for
`process_local_id()`. The existing `delete_closed_local_sqlite_dbs()` scans a
whole local root and is too broad for finalization.

`delete_active_global_sqlite_copy` should do nothing unless
`sqlite_run_mode(env) == "server"`. In server mode it should delete only the
active node-global copy for `S`, never `sqlite_global_source_dir(env)`.

## Phase 6: checkpoint and file-copy safety

Before copying a source global database to a node-local global database, make
the copy robust with WAL sidecars:

- If the source global DB is opened by LurCGT in this process, run a passive or
  truncate checkpoint before copying.
- Copy the `.db`, `-wal`, and `-shm` files as a set, as the current code does.
- Prefer not to copy while holding the merge lock for a long time; copying a
  large database under the global lock can serialize all job startup. If
  consistency problems appear in practice, add a short source-global copy lock
  separate from the merge lock.

For the first implementation, the existing file-set copy is acceptable if the
source global DB is mostly append-only and merges are protected by the merge
lock. Document the assumption.

## Phase 7: tests

Add focused tests to the existing SQLite test area rather than creating large
fixtures.

Required regression tests:

- In server mode, `merge_all_to_global(SU{2})` inserts rows into
  `LURCGT_GLOBALDB_DIR/SU2.db`, not `LURCGT_GLOBALDB_DIR_NODE/SU2.db`.
- In server mode, `merge_all_to_global(SU{2})` creates
  `LURCGT_GLOBALDB_DIR/SU2.db` if the source global DB does not exist.
- In server mode, opening a local DB creates a fresh node-local local DB even
  if a same-named source local DB exists.
- In server mode, opening a global DB still copies source global rows to the
  node-global DB for fast reads.
- In server mode, opening or writing a local DB does not copy any global DB.
- In server mode, requesting global data for `SU{2}` copies only `SU2.db` and
  does not copy `SU3.db` or any other symmetry DB present in
  `LURCGT_GLOBALDB_DIR`.
- In local mode, merge behavior is unchanged.
- In local mode, finalization creates the active global DB for `S` if it does
  not exist, then merges the current local rows into it.
- In local mode, `finalize_sqlite!(S; cleanup=true)` merges rows into the
  active global DB and deletes only the current process-local DB for `S`.
- Concurrent merge behavior keeps `INSERT OR IGNORE` semantics. A minimal
  single-process duplicate-row test is enough unless a reliable multiprocess
  test already exists.
- In server mode, `finalize_sqlite!(S; cleanup=true)` deletes node-local files
  after a successful merge and leaves source global files intact.

Do not add broad CGT-generation tests for this refactor. The database behavior
can be tested by inserting small `UInt8` payload rows into existing tables.

## Validation

Run:

```bash
julia --project=LurCGT LurCGT/test/runtests.jl
```

For a quick manual Slurm-like check, use temporary directories:

```julia
withenv(
    "LURCGT_RUN_MODE" => "server",
    "LURCGT_GLOBALDB_DIR" => source_global,
    "LURCGT_GLOBALDB_DIR_NODE" => node_global,
    "LURCGT_LOCALDB_DIR_NODE" => node_local,
) do
    # create/load/save/finalize
end
```

Inspect both files after merge:

- `source_global/$(totxt(S)).db` contains newly generated rows.
- `node_global/$(totxt(S)).db` is only a compute-time read cache.

For local mode, check that `finalize_sqlite!` leaves the global DB in place and
removes only the current process-local DB for the finalized symmetry.

## Risks

The main risk is SQLite consistency when copying a source global DB while
another process is merging. `INSERT OR IGNORE` protects merge conflicts, but a
node-local copy could miss rows committed moments later. This is acceptable:
the local job may recompute an object and then skip it at merge time if another
job inserted it first.

The second risk is cleanup after failed merge. Cleanup must happen only after a
successful merge and connection close.

The third risk is ambiguity in names. Keep `sqlite_global_dir()` as the active
runtime global DB and add explicit `sqlite_merge_global_*` names for the
persistent merge target.

Expected error cases are limited to real filesystem or SQLite failures:

- the global or local directory cannot be created because of permissions;
- the filesystem is full or the quota is exceeded;
- an existing DB file is corrupt or is not a SQLite database;
- SQLite cannot acquire a lock because another process holds a stale or long
  running lock;
- WAL sidecar files are inconsistent because an external process copied or
  deleted only part of a SQLite file set;
- cleanup tries to remove a file that is still open by another process.

Do not catch and silently ignore these errors during merge. The safe behavior is
to leave the local DB files in place so the user can retry the merge.
