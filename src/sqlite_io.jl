# SQLite-based I/O for LurCGT
#
# Architecture:
# - ONE global database (shared read, WAL mode)
# - ONE local database PER PROCESS (no write contention)
# - Per-type LRU caches (sized to match typical object sizes)
# - Tables created at DB initialization
# - Performance PRAGMAs applied at connection time

using SQLite
using DBInterface
using Serialization
using LRUCache: LRU
using FileWatching: mkpidlock

# ============================================================================
# Configuration
# ============================================================================

const SQLITE_DBS = Dict{String, SQLite.DB}()
const SQLITE_DB_LOCK = ReentrantLock()

# Per-type LRU caches (sized to match typical object sizes)
const IREP_CACHE = LRU{String, Any}(maxsize=100*1024*1024, by=x -> x.size_byte)
const IREP_CACHE_LOCK = ReentrantLock()

const CGT_CACHE = LRU{String, Any}(maxsize=200*1024*1024, by=x -> x.size_byte)
const CGT_CACHE_LOCK = ReentrantLock()

const FSYMBOL_CACHE = LRU{String, Any}(maxsize=10*1024*1024, by=x -> x.size_byte)
const FSYMBOL_CACHE_LOCK = ReentrantLock()

const RSYMBOL_CACHE = LRU{String, Any}(maxsize=10*1024*1024, by=x -> x.size_byte)
const RSYMBOL_CACHE_LOCK = ReentrantLock()

const XSYMBOL_CACHE = LRU{String, Any}(maxsize=50*1024*1024, by=x -> x.size_byte)
const XSYMBOL_CACHE_LOCK = ReentrantLock()

const OMLIST_CACHE = LRU{String, Any}(maxsize=10*1024*1024, by=x -> x.size_byte)
const OMLIST_CACHE_LOCK = ReentrantLock()

const VALIDOUT_CACHE = LRU{String, Any}(maxsize=10*1024*1024, by=x -> x.size_byte)
const VALIDOUT_CACHE_LOCK = ReentrantLock()

const CGTPERM_CACHE = LRU{String, Any}(maxsize=20*1024*1024, by=x -> x.size_byte)
const CGTPERM_CACHE_LOCK = ReentrantLock()

const CONJPERM_CACHE = LRU{String, Any}(maxsize=20*1024*1024, by=x -> x.size_byte)
const CONJPERM_CACHE_LOCK = ReentrantLock()

const CGTSVD_CACHE = LRU{String, Any}(maxsize=20*1024*1024, by=x -> x.size_byte)
const CGTSVD_CACHE_LOCK = ReentrantLock()

const CG3FLIP_CACHE = LRU{String, Any}(maxsize=20*1024*1024, by=x -> x.size_byte)
const CG3FLIP_CACHE_LOCK = ReentrantLock()

function _nonempty_env_value(env::AbstractDict, key::AbstractString)
    haskey(env, key) || return nothing
    value = strip(String(env[key]))
    return isempty(value) ? nothing : value
end

_resolve_storage_dir(path::AbstractString) = normpath(abspath(expanduser(path)))

"""
    sqlite_run_mode([env])

Return the active SQLite storage mode, either `"local"` or `"server"`.

The mode is read from `LURCGT_RUN_MODE`; unset means `"local"`.  In `"local"`
mode, LurCGT reads and writes directly under `LURCGT_LOCALDB_DIR` and
`LURCGT_GLOBALDB_DIR`.  In `"server"` mode, active database files are placed
under the node-local directories and copied from source directories on first
open.
"""
function sqlite_run_mode(env::AbstractDict=ENV)
    mode = lowercase(something(_nonempty_env_value(env, "LURCGT_RUN_MODE"), "local"))
    mode in ("local", "server") ||
        throw(ArgumentError("LURCGT_RUN_MODE must be \"local\" or \"server\", got \"$(mode)\""))
    return mode
end

sqlite_local_source_dir(env::AbstractDict=ENV) =
    _resolve_storage_dir(something(_nonempty_env_value(env, "LURCGT_LOCALDB_DIR"),
                                   joinpath(homedir(), ".LurCGT_sqlite", "local")))

sqlite_global_source_dir(env::AbstractDict=ENV) =
    _resolve_storage_dir(something(_nonempty_env_value(env, "LURCGT_GLOBALDB_DIR"),
                                   joinpath(homedir(), ".LurCGT_sqlite", "global")))

sqlite_local_node_dir(env::AbstractDict=ENV) =
    _resolve_storage_dir(something(_nonempty_env_value(env, "LURCGT_LOCALDB_DIR_NODE"),
                                   sqlite_local_source_dir(env)))

sqlite_global_node_dir(env::AbstractDict=ENV) =
    _resolve_storage_dir(something(_nonempty_env_value(env, "LURCGT_GLOBALDB_DIR_NODE"),
                                   sqlite_global_source_dir(env)))

"""
    sqlite_local_dir([env])

Return the directory used for process-local SQLite databases in the active run
mode.

Local databases are stored as one database per symmetry and process id, for
example `joinpath(sqlite_local_dir(), "SU2", "\$(process_local_id()).db")`.
Set `LURCGT_LOCALDB_DIR` to choose the source/local directory, and
`LURCGT_LOCALDB_DIR_NODE` to choose the active node-local directory in
`LURCGT_RUN_MODE=server`.
"""
sqlite_local_dir(env::AbstractDict=ENV) =
    sqlite_run_mode(env) == "server" ? sqlite_local_node_dir(env) : sqlite_local_source_dir(env)

"""
    sqlite_global_dir([env])

Return the directory used for shared global SQLite databases in the active run
mode.

Global databases are stored as one database per symmetry, for example
`joinpath(sqlite_global_dir(), "SU2.db")`.  Set `LURCGT_GLOBALDB_DIR` to choose
the source/global directory, and `LURCGT_GLOBALDB_DIR_NODE` to choose the active
node-local directory in `LURCGT_RUN_MODE=server`.
"""
sqlite_global_dir(env::AbstractDict=ENV) =
    sqlite_run_mode(env) == "server" ? sqlite_global_node_dir(env) : sqlite_global_source_dir(env)

"""
    sqlite_lock_dir([env])

Return the directory used for cross-process SQLite merge locks.

Merge operations create lock files here so concurrent jobs can safely call
`merge_table_to_global` or `merge_all_to_global` for the same symmetry.
"""
sqlite_lock_dir(env::AbstractDict=ENV) = joinpath(sqlite_global_dir(env), "locks")

function sqlite_db_source_path(::Type{S}, location::Symbol=:local;
                               env::AbstractDict=ENV,
                               process_id::AbstractString=process_local_id()) where {S<:NonabelianSymm}
    @assert location in (:local, :global)
    symm_name = totxt(S)
    if location == :global
        return joinpath(sqlite_global_source_dir(env), "$(symm_name).db")
    end
    return joinpath(sqlite_local_source_dir(env), symm_name, "$(process_id).db")
end

function sqlite_db_path(::Type{S}, location::Symbol=:local;
                        env::AbstractDict=ENV,
                        process_id::AbstractString=process_local_id()) where {S<:NonabelianSymm}
    @assert location in (:local, :global)
    symm_name = totxt(S)
    if location == :global
        return joinpath(sqlite_global_dir(env), "$(symm_name).db")
    end
    return joinpath(sqlite_local_dir(env), symm_name, "$(process_id).db")
end

function _copy_sqlite_file_set(source_path::AbstractString, target_path::AbstractString)
    mkpath(dirname(target_path))
    for suffix in ("", "-wal", "-shm")
        source_file = string(source_path, suffix)
        target_file = string(target_path, suffix)
        isfile(source_file) || continue
        cp(source_file, target_file; force=true)
    end
    return nothing
end

function _prepare_server_sqlite_db(::Type{S}, location::Symbol, path::AbstractString;
                                   env::AbstractDict=ENV,
                                   process_id::AbstractString=process_local_id()) where {S<:NonabelianSymm}
    sqlite_run_mode(env) == "server" || return nothing
    isfile(path) && return nothing

    source_path = sqlite_db_source_path(S, location; env, process_id)
    path == source_path && return nothing

    if isfile(source_path)
        _copy_sqlite_file_set(source_path, path)
    end
    return nothing
end

function build_process_local_id(env::AbstractDict=ENV;
                                hostname::AbstractString=gethostname(),
                                pid::Integer=getpid())
    host_part = replace(String(hostname), r"[^A-Za-z0-9._-]" => "_")
    parts = String[host_part]

    if (job_id = _nonempty_env_value(env, "SLURM_JOB_ID")) !== nothing
        push!(parts, "job$(job_id)")
    end

    if (task_id = _nonempty_env_value(env, "SLURM_ARRAY_TASK_ID")) !== nothing
        push!(parts, "task$(task_id)")
    elseif (proc_id = _nonempty_env_value(env, "SLURM_PROCID")) !== nothing
        push!(parts, "proc$(proc_id)")
    end

    push!(parts, "pid$(pid)")
    return join(parts, "_")
end

const PROCESS_LOCAL_ID = Ref{String}("")

process_local_id() = PROCESS_LOCAL_ID[]
process_local_id(env::AbstractDict; hostname::AbstractString=gethostname(), pid::Integer=getpid()) =
    build_process_local_id(env; hostname, pid)

# All table names (created at DB init)
const SQLITE_TABLES = ("irreps", "cgt", "fsymbol", "rsymbol", "xsymbol",
                        "omlist", "validout", "cgtperm", "conjperm", "cgtsvd", "cg3flip")

# ============================================================================
# Database Management
# ============================================================================

"""Initialize a new SQLite DB with WAL mode, performance PRAGMAs, and all tables."""
function _init_sqlite_db(db::SQLite.DB)
    DBInterface.execute(db, "PRAGMA journal_mode=WAL")
    DBInterface.execute(db, "PRAGMA synchronous=NORMAL")
    DBInterface.execute(db, "PRAGMA cache_size=-64000")   # 64 MB
    DBInterface.execute(db, "PRAGMA mmap_size=268435456")  # 256 MB
    DBInterface.execute(db, "PRAGMA temp_store=MEMORY")
    for table in SQLITE_TABLES
        DBInterface.execute(db, "CREATE TABLE IF NOT EXISTS $(table) (key TEXT PRIMARY KEY, data BLOB NOT NULL)")
    end
end

"""
Get or create SQLite database for a symmetry type.

Locations:
- :global - Shared database, read by all processes, written only during merge
- :local  - Per-process database (one per process ID), no write contention
"""
function get_sqlite_db(::Type{S}, location::Symbol=:local) where {S<:NonabelianSymm}
    @assert location in (:local, :global)
    path = sqlite_db_path(S, location)

    Base.lock(SQLITE_DB_LOCK) do
        _prepare_server_sqlite_db(S, location, path)
        db_key = "$(location):$(path)"
        get!(SQLITE_DBS, db_key) do
            mkpath(dirname(path))
            db = SQLite.DB(path)
            _init_sqlite_db(db)
            db
        end
    end
end

# ============================================================================
# Serialization Helpers
# ============================================================================

function serialize_object(obj)
    buf = IOBuffer()
    serialize(buf, obj)
    return take!(buf)
end

function deserialize_object(data::Vector{UInt8})
    return deserialize(IOBuffer(data))
end

# Backward compatibility aliases
serialize_irep(rep) = serialize_object(rep)
deserialize_irep(data) = deserialize_object(data)

# ============================================================================
# Generic Save/Load (core DB operations)
# ============================================================================

"""
Save object to SQLite local database.
Tables are pre-created at DB init, so no CREATE TABLE here.
"""
function save_object_sqlite(::Type{S}, table_name::String, key::String, obj;
                           verbose=0) where {S<:NonabelianSymm}
    db = get_sqlite_db(S, :local)
    data = serialize_object(obj)
    DBInterface.execute(db,
        "INSERT OR REPLACE INTO $(table_name) (key, data) VALUES (?, ?)",
        [key, data])
    verbose > 0 && println("Saved $(table_name) to local DB ($(process_local_id())): $(key)")
    return nothing
end

"""
Load object from SQLite. Tries global DB first, then local.
"""
function load_object_sqlite(::Type{S}, table_name::String, key::String;
                           verbose=0) where {S<:NonabelianSymm}
    for location in (:global, :local)
        try
            db = get_sqlite_db(S, location)
            result_rows = DBInterface.execute(db,
                "SELECT data FROM $(table_name) WHERE key = ?", [key])
            for row in result_rows
                data = row.data
                if data !== nothing && !isempty(data)
                    result = deserialize_object(data)
                    verbose > 0 && println("✓ Found $(table_name) in $(location): $(key)")
                    return result
                end
            end
        catch e
            isa(e, SQLite.SQLiteException) ? continue : rethrow(e)
        end
    end
    verbose > 0 && println("✗ Not found $(table_name): $(key)")
    return nothing
end

# ============================================================================
# Generic Cached Save/Load (cache + DB)
# ============================================================================

"""Save object to DB and update the given cache."""
function _cached_save(::Type{S}, table_name::String, db_key::String, cache_key::String, obj,
                      cache::LRU, cache_lock::ReentrantLock;
                      verbose=0) where {S<:NonabelianSymm}
    save_object_sqlite(S, table_name, db_key, obj; verbose)
    Base.lock(cache_lock) do
        cache[cache_key] = obj
    end
    return nothing
end

"""Load object from the given cache, falling back to DB. Updates cache on DB hit."""
function _cached_load(::Type{S}, table_name::String, db_key::String, cache_key::String,
                      cache::LRU, cache_lock::ReentrantLock;
                      verbose=0) where {S<:NonabelianSymm}
    # 1. Check cache
    cached = Base.lock(cache_lock) do
        get(cache, cache_key, nothing)
    end
    cached !== nothing && return cached

    # 2. Load from DB
    result = load_object_sqlite(S, table_name, db_key; verbose)

    # 3. Update cache on hit
    if result !== nothing
        Base.lock(cache_lock) do
            cache[cache_key] = result
        end
    end
    return result
end

# ============================================================================
# Key Generation Functions
# ============================================================================

function irep_key(::Type{S}, ::Type{RT}, qlabel::NTuple{NZ, Int}) where {S<:Symmetry, RT, NZ}
    return "irep_$(totxt(S))_$(RT)_$(join(qlabel, "_"))"
end

function cgt_key(::Type{S}, ::Type{CT},
    qlabels::NTuple{N, NTuple{NZ,Int}}) where {S<:Symmetry, CT, NZ, N}
    @assert N == 2 || N == 3; @assert qlabels[1] <= qlabels[2]
    in1_str = join(qlabels[1], "_")
    in2_str = join(qlabels[2], "_")
    out_str = N == 2 ? "1j" : join(qlabels[3], "_")
    return "cgt_$(totxt(S))_$(CT)_$(in1_str)_$(in2_str)_$(out_str)"
end

function fsymbol_key(::Type{S}, ::Type{CT}, in1::NTuple{NZ,Int}, in2::NTuple{NZ,Int},
                     in3::NTuple{NZ,Int}, out::NTuple{NZ,Int}) where {S<:Symmetry, CT, NZ}
    return "fsymbol_$(totxt(S))_$(CT)_$(join(in1,"_"))_$(join(in2,"_"))_$(join(in3,"_"))_$(join(out,"_"))"
end

function rsymbol_key(::Type{S}, ::Type{CT}, in::NTuple{NZ,Int}, out::NTuple{NZ,Int}) where {S<:Symmetry, CT, NZ}
    return "rsymbol_$(totxt(S))_$(CT)_$(join(in,"_"))_$(join(out,"_"))"
end

function xsymbol_key(::Type{S},
    up1::NTuple{U1, NTuple{NZ, Int}},
    dn1::NTuple{D1, NTuple{NZ, Int}},
    up2::NTuple{U2, NTuple{NZ, Int}},
    dn2::NTuple{D2, NTuple{NZ, Int}},
    legs1::NTuple{M, Int},
    legs2::NTuple{M, Int}) where {S<:Symmetry, NZ, U1, D1, U2, D2, M}
    up1_str = join([join(q, "_") for q in up1], "_")
    dn1_str = join([join(q, "_") for q in dn1], "_")
    up2_str = join([join(q, "_") for q in up2], "_")
    dn2_str = join([join(q, "_") for q in dn2], "_")
    return "xsymbol_$(totxt(S))_$(up1_str):$(dn1_str);$(up2_str):$(dn2_str);$(join(legs1,"_")),$(join(legs2,"_"))"
end

function omlist_key(::Type{S},
    qs::NTuple{N, NTuple{NZ, Int}},
    outsp::NTuple{NZ, Int}) where {S<:Symmetry, N, NZ}
    qs_strs = [join(q, "_") for q in qs]
    return "omlist_$(totxt(S))_$(join(qs_strs, "_"))_$(join(outsp, "_"))"
end

function validout_key(::Type{S},
    qs::NTuple{N, NTuple{NZ, Int}}) where {S<:Symmetry, N, NZ}
    sorted_qs = sort(collect(qs))
    qs_strs = [join(q, "_") for q in sorted_qs]
    return "validout_$(totxt(S))_$(join(qs_strs, "_"))"
end

function cgtperm_key(::Type{S}, upsp::Tuple, dnsp::Tuple, perm::NTuple{N, Int}) where {S<:Symmetry, N}
    up_strs = [join(q, "_") for q in upsp]
    dn_strs = [join(q, "_") for q in dnsp]
    return "cgtperm_$(totxt(S))_$(join(up_strs, "_")),$(join(dn_strs, "_")),$(join(perm, "_"))"
end

function conjperm_key(::Type{S}, upsp::Tuple) where {S<:Symmetry}
    up_strs = [join(q, "_") for q in upsp]
    return "conjperm_$(totxt(S))_$(join(up_strs, "_"))"
end

function cgtsvd_key(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    leftlegs::NTuple{L, Int}) where {S<:Symmetry, U, D, NZ, L}
    up_strs = [join(q, "_") for q in upsp]
    dn_strs = [join(q, "_") for q in dnsp]
    return "cgtsvd_v2_$(totxt(S))_$(join(up_strs, "_")),$(join(dn_strs, "_")),$(join(leftlegs, "_"))"
end

function cg3flip_key(::Type{S},
    ::Type{CT},
    incom_spaces::NTuple{2, NTuple{NZ, Int}},
    out_space::NTuple{NZ, Int}) where {S<:Symmetry, CT<:Number, NZ}
    in1, in2 = incom_spaces
    return "cg3flip_$(totxt(S))_$(CT)_$(join(in1,"_"))_$(join(in2,"_"))_$(join(out_space,"_"))"
end

# ============================================================================
# Type-Specific Save/Load Functions
# ============================================================================

# --- Irep ---

function save_irep_sqlite(rep::Irep{S, NL, NZ, RT}; verbose=0) where {S, NL, NZ, RT}
    key = irep_key(S, RT, rep.qlabel)
    cache_key = "$(S)_$(RT)_$(join(rep.qlabel, "_"))"
    _cached_save(S, "irreps", key, cache_key, rep, IREP_CACHE, IREP_CACHE_LOCK; verbose)
end

function load_irep_sqlite(::Type{S}, ::Type{RT}, qlabel::NTuple{NZ, Int};
                          verbose=0) where {S<:NonabelianSymm, RT<:Number, NZ}
    key = irep_key(S, RT, qlabel)
    cache_key = "$(S)_$(RT)_$(join(qlabel, "_"))"
    return _cached_load(S, "irreps", key, cache_key, IREP_CACHE, IREP_CACHE_LOCK; verbose)
end

# --- CGT ---

function save_cgt_sqlite(::Type{S},
    cgt::CGT{S, CT, NZ, N};
    verbose=0) where {S<:NonabelianSymm, CT<:Number, NZ, N}
    key = cgt_key(S, CT, cgt.qlabels)
    cache_key = "$(S)_$(CT)_$(join([join(q, "_") for q in cgt.qlabels], "__"))"
    _cached_save(S, "cgt", key, cache_key, cgt, CGT_CACHE, CGT_CACHE_LOCK; verbose)
end

function load_cgt_sqlite(::Type{S},
    ::Type{CT},
    qlabels::NTuple{N, NTuple{NZ,Int}};
    verbose=0) where {S<:NonabelianSymm, CT, NZ, N}
    key = cgt_key(S, CT, qlabels)
    cache_key = "$(S)_$(CT)_$(join([join(q, "_") for q in qlabels], "__"))"
    return _cached_load(S, "cgt", key, cache_key, CGT_CACHE, CGT_CACHE_LOCK; verbose)
end

# --- Fsymbol ---

function save_Fsymbol_sqlite(fsym::Fsymbol{S, CT, NZ}; verbose=0) where {S<:NonabelianSymm, CT, NZ}
    save_Fsymbol_sqlite(S, fsym; verbose)
end

function save_Fsymbol_sqlite(::Type{S}, fsym; verbose=0) where {S<:NonabelianSymm}
    CT = eltype(fsym)
    key = fsymbol_key(S, CT, fsym.in1, fsym.in2, fsym.in3, fsym.out)
    cache_key = "$(S)_$(CT)_$(join(fsym.in1,"_"))_$(join(fsym.in2,"_"))_$(join(fsym.in3,"_"))_$(join(fsym.out,"_"))"
    _cached_save(S, "fsymbol", key, cache_key, fsym, FSYMBOL_CACHE, FSYMBOL_CACHE_LOCK; verbose)
end

function load_Fsymbol_sqlite(::Type{S}, ::Type{CT}, in1::NTuple{NZ,Int}, in2::NTuple{NZ,Int},
                            in3::NTuple{NZ,Int}, out::NTuple{NZ,Int};
                            verbose=0) where {S<:NonabelianSymm, CT, NZ}
    key = fsymbol_key(S, CT, in1, in2, in3, out)
    cache_key = "$(S)_$(CT)_$(join(in1,"_"))_$(join(in2,"_"))_$(join(in3,"_"))_$(join(out,"_"))"
    return _cached_load(S, "fsymbol", key, cache_key, FSYMBOL_CACHE, FSYMBOL_CACHE_LOCK; verbose)
end

# --- Rsymbol ---

function save_Rsymbol_sqlite(rsym::Rsymbol{S, CT, NZ}; verbose=0) where {S<:NonabelianSymm, CT, NZ}
    save_Rsymbol_sqlite(S, rsym; verbose)
end

function save_Rsymbol_sqlite(::Type{S}, rsym; verbose=0) where {S<:NonabelianSymm}
    CT = eltype(rsym)
    key = rsymbol_key(S, CT, rsym.in, rsym.out)
    cache_key = "$(S)_$(CT)_$(join(rsym.in,"_"))_$(join(rsym.out,"_"))"
    _cached_save(S, "rsymbol", key, cache_key, rsym, RSYMBOL_CACHE, RSYMBOL_CACHE_LOCK; verbose)
end

function load_Rsymbol_sqlite(::Type{S}, ::Type{CT}, in::NTuple{NZ,Int}, out::NTuple{NZ,Int};
                            verbose=0) where {S<:NonabelianSymm, CT, NZ}
    key = rsymbol_key(S, CT, in, out)
    cache_key = "$(S)_$(CT)_$(join(in,"_"))_$(join(out,"_"))"
    return _cached_load(S, "rsymbol", key, cache_key, RSYMBOL_CACHE, RSYMBOL_CACHE_LOCK; verbose)
end

# --- Xsymbol ---

function save_Xsymbol_sqlite(xsym::Xsymbol{S, U1, D1, U2, D2, NZ, M}; verbose=0) where {S<:NonabelianSymm, U1, D1, U2, D2, NZ, M}
    save_Xsymbol_sqlite(S, xsym; verbose)
end

function save_Xsymbol_sqlite(::Type{S}, xsym; verbose=0) where {S<:NonabelianSymm}
    key = xsymbol_key(S, xsym.up1sp, xsym.dn1sp, xsym.up2sp, xsym.dn2sp, xsym.legs1, xsym.legs2)
    _cached_save(S, "xsymbol", key, key, xsym, XSYMBOL_CACHE, XSYMBOL_CACHE_LOCK; verbose)
end

function load_Xsymbol_sqlite(::Type{S},
    up1::NTuple{U1, NTuple{NZ, Int}},
    dn1::NTuple{D1, NTuple{NZ, Int}},
    up2::NTuple{U2, NTuple{NZ, Int}},
    dn2::NTuple{D2, NTuple{NZ, Int}},
    legs1::NTuple{M, Int},
    legs2::NTuple{M, Int};
    verbose=0) where {S<:NonabelianSymm, NZ, U1, D1, U2, D2, M}
    key = xsymbol_key(S, up1, dn1, up2, dn2, legs1, legs2)
    return _cached_load(S, "xsymbol", key, key, XSYMBOL_CACHE, XSYMBOL_CACHE_LOCK; verbose)
end

# --- OMList ---

function save_omlist_sqlite(omlist::OMList{S, N, NZ, M1, M2}; verbose=0) where {S<:NonabelianSymm, N, NZ, M1, M2}
    save_omlist_sqlite(S, omlist; verbose)
end

function save_omlist_sqlite(::Type{S}, omlist; verbose=0) where {S<:NonabelianSymm}
    key = omlist_key(S, omlist.incom_spaces, omlist.out_space)
    _cached_save(S, "omlist", key, key, omlist, OMLIST_CACHE, OMLIST_CACHE_LOCK; verbose)
end

function load_omlist_sqlite(::Type{S},
    qs::NTuple{N, NTuple{NZ, Int}},
    outsp::NTuple{NZ, Int};
    verbose=0) where {S<:NonabelianSymm, N, NZ}
    key = omlist_key(S, qs, outsp)
    return _cached_load(S, "omlist", key, key, OMLIST_CACHE, OMLIST_CACHE_LOCK; verbose)
end

# --- ValidOuts ---

function save_validout_sqlite(validout::ValidOuts{S, N, NZ}; verbose=0) where {S<:NonabelianSymm, N, NZ}
    save_validout_sqlite(S, validout; verbose)
end

function save_validout_sqlite(::Type{S}, validout; verbose=0) where {S<:NonabelianSymm}
    key = validout_key(S, validout.incom_spaces)
    _cached_save(S, "validout", key, key, validout, VALIDOUT_CACHE, VALIDOUT_CACHE_LOCK; verbose)
end

function load_validout_sqlite(::Type{S},
    qs::NTuple{N, NTuple{NZ, Int}};
    verbose=0) where {S<:NonabelianSymm, N, NZ}
    key = validout_key(S, qs)
    return _cached_load(S, "validout", key, key, VALIDOUT_CACHE, VALIDOUT_CACHE_LOCK; verbose)
end

# --- CGTperm ---

function save_CGTperm_sqlite(cgtperm::CGTperm{S, U, D, N, NZ}; verbose=0) where {S<:NonabelianSymm, U, D, N, NZ}
    save_CGTperm_sqlite(S, cgtperm; verbose)
end

function save_CGTperm_sqlite(::Type{S}, cgtperm; verbose=0) where {S<:NonabelianSymm}
    key = cgtperm_key(S, cgtperm.upsp, cgtperm.dnsp, cgtperm.perm)
    _cached_save(S, "cgtperm", key, key, cgtperm, CGTPERM_CACHE, CGTPERM_CACHE_LOCK; verbose)
end

function load_CGTperm_sqlite(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    perm::NTuple{N, Int};
    verbose=0) where {S<:NonabelianSymm, U, D, N, NZ}
    key = cgtperm_key(S, upsp, dnsp, perm)
    return _cached_load(S, "cgtperm", key, key, CGTPERM_CACHE, CGTPERM_CACHE_LOCK; verbose)
end

# --- Conjperm ---

function save_Conjperm_sqlite(conjperm::Conjperm{S, U, NZ}; verbose=0) where {S<:NonabelianSymm, U, NZ}
    save_Conjperm_sqlite(S, conjperm; verbose)
end

function save_Conjperm_sqlite(::Type{S}, conjperm; verbose=0) where {S<:NonabelianSymm}
    key = conjperm_key(S, conjperm.upsp)
    _cached_save(S, "conjperm", key, key, conjperm, CONJPERM_CACHE, CONJPERM_CACHE_LOCK; verbose)
end

function load_Conjperm_sqlite(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}};
    verbose=0) where {S<:NonabelianSymm, U, NZ}
    key = conjperm_key(S, upsp)
    return _cached_load(S, "conjperm", key, key, CONJPERM_CACHE, CONJPERM_CACHE_LOCK; verbose)
end

# --- CGTSVD ---

function save_CGTSVD_sqlite(cgtsvd::CGTSVD{S, U, D, L, NZ}; verbose=0) where {S<:NonabelianSymm, U, D, L, NZ}
    save_CGTSVD_sqlite(S, cgtsvd; verbose)
end

function save_CGTSVD_sqlite(::Type{S}, cgtsvd; verbose=0) where {S<:NonabelianSymm}
    key = cgtsvd_key(S, cgtsvd.upsp, cgtsvd.dnsp, cgtsvd.leftlegs)
    _cached_save(S, "cgtsvd", key, key, cgtsvd, CGTSVD_CACHE, CGTSVD_CACHE_LOCK; verbose)
end

function load_CGTSVD_sqlite(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    leftlegs::NTuple{L, Int};
    verbose=0) where {S<:NonabelianSymm, U, D, L, NZ}
    key = cgtsvd_key(S, upsp, dnsp, leftlegs)
    return _cached_load(S, "cgtsvd", key, key, CGTSVD_CACHE, CGTSVD_CACHE_LOCK; verbose)
end

# --- CG3flip ---

function save_cg3flip_sqlite(::Type{S},
    cg3flip::CG3Flip{S, CT, NZ};
    verbose=0) where {S<:NonabelianSymm, CT<:Number, NZ}
    key = cg3flip_key(S, CT, cg3flip.incom_spaces, cg3flip.out_space)
    _cached_save(S, "cg3flip", key, key, cg3flip, CG3FLIP_CACHE, CG3FLIP_CACHE_LOCK; verbose)
end

function load_cg3flip_sqlite(::Type{S},
    ::Type{CT},
    incom_spaces::NTuple{2, NTuple{NZ, Int}},
    out_space::NTuple{NZ, Int};
    verbose=0) where {S<:NonabelianSymm, CT<:Number, NZ}
    key = cg3flip_key(S, CT, incom_spaces, out_space)
    return _cached_load(S, "cg3flip", key, key, CG3FLIP_CACHE, CG3FLIP_CACHE_LOCK; verbose)
end

# ============================================================================
# Merge Operations
# ============================================================================

"""
    merge_table_to_global(S, table_name[, filter_prefix]; clear_local_after=true, verbose=1)

Merge entries for one SQLite table from the current process-local database into
the shared global database for symmetry `S`.

`table_name` must be one of the internal cache tables, for example `"irreps"`,
`"cgt"`, `"xsymbol"`, `"cgtperm"`, or `"cgtsvd"`.  If `filter_prefix` is
provided, only rows whose key starts with that prefix are merged.

The merge is protected by a file lock under `sqlite_lock_dir()` so multiple
processes can merge safely.  Rows are inserted with `INSERT OR IGNORE`; existing
global rows are counted as skipped.  When `clear_local_after=true`, successfully
merged rows are removed from the local table.

Returns a named tuple `(merged=..., skipped=...)`.
"""
function merge_table_to_global(::Type{S}, table_name::String, filter_prefix::String="";
                               clear_local_after::Bool=true,
                               verbose=1) where {S<:NonabelianSymm}
    local_db = get_sqlite_db(S, :local)
    global_db = get_sqlite_db(S, :global)

    lockdir = sqlite_lock_dir()
    mkpath(lockdir)
    lockpath = joinpath(lockdir, "$(totxt(S))_merge.lock")

    mkpidlock(lockpath) do
        if verbose > 0
            println("Merging $(totxt(S)) $(table_name) from local → global ($(process_local_id()))")
        end

        # Read from local
        local_rows = try
            if isempty(filter_prefix)
                DBInterface.execute(local_db, "SELECT key, data FROM $(table_name)")
            else
                DBInterface.execute(local_db,
                    "SELECT key, data FROM $(table_name) WHERE key LIKE ?",
                    ["$(filter_prefix)%"])
            end
        catch
            verbose > 0 && println("Table $(table_name) not found in local DB")
            return (merged=0, skipped=0)
        end

        # Collect rows first (can't nest SQLite queries)
        rows = [(row.key, row.data) for row in local_rows]

        if isempty(rows)
            verbose > 0 && println("No rows to merge in $(table_name)")
            return (merged=0, skipped=0)
        end

        # Bulk insert with INSERT OR IGNORE in a transaction
        # Count rows before/after to determine how many were actually inserted
        count_before = first(DBInterface.execute(global_db,
            "SELECT COUNT(*) AS n FROM $(table_name)")).n
        DBInterface.execute(global_db, "BEGIN TRANSACTION")
        try
            for (key, data) in rows
                DBInterface.execute(global_db,
                    "INSERT OR IGNORE INTO $(table_name) (key, data) VALUES (?, ?)",
                    [key, data])
            end
            DBInterface.execute(global_db, "COMMIT")
        catch e
            DBInterface.execute(global_db, "ROLLBACK")
            rethrow(e)
        end
        count_after = first(DBInterface.execute(global_db,
            "SELECT COUNT(*) AS n FROM $(table_name)")).n
        merged_count = count_after - count_before

        skipped_count = length(rows) - merged_count
        if verbose > 0
            println("Merged: $(merged_count) | Skipped: $(skipped_count)")
        end

        # Clear local after successful merge
        if clear_local_after && merged_count > 0
            if isempty(filter_prefix)
                DBInterface.execute(local_db, "DELETE FROM $(table_name)")
            else
                DBInterface.execute(local_db,
                    "DELETE FROM $(table_name) WHERE key LIKE ?",
                    ["$(filter_prefix)%"])
            end
            verbose > 0 && println("Cleared $(table_name) from local DB")
        end

        return (merged=merged_count, skipped=skipped_count)
    end
end

"""
    merge_all_to_global(S; tables=collect(SQLITE_TABLES), clear_local_after=true, verbose=1)

Merge all selected SQLite cache tables from the current process-local database
into the shared global database for symmetry `S`.

This is the usual user-facing merge operation after a local run has generated
irreps, CGTs, F/R/X symbols, CGT permutations, or CGTSVD data.  Pass `tables` to
merge only selected tables.  Each table is merged by `merge_table_to_global`, so
cross-process locking and `clear_local_after` behavior are the same.

Returns a named tuple `(merged=..., skipped=...)` with totals across tables.
"""
function merge_all_to_global(::Type{S};
                             tables=collect(SQLITE_TABLES),
                             clear_local_after::Bool=true,
                             verbose=1) where {S<:NonabelianSymm}
    total_merged = 0
    total_skipped = 0

    for table in tables
        result = merge_table_to_global(S, table; clear_local_after, verbose)
        total_merged += result.merged
        total_skipped += result.skipped
    end

    if verbose > 0
        println("TOTAL: Merged $(total_merged) | Skipped $(total_skipped)")
    end
    return (merged=total_merged, skipped=total_skipped)
end

# ============================================================================
# Utility Functions
# ============================================================================

function list_ireps_sqlite(::Type{S}, ::Type{RT}; location::Symbol=:local) where {S, RT}
    db = get_sqlite_db(S, location)
    keys = String[]
    try
        result = DBInterface.execute(db,
            "SELECT key FROM irreps WHERE key LIKE ?",
            ["irep_$(totxt(S))_$(RT)_%"])
        for row in result
            push!(keys, row.key)
        end
    catch e
        isa(e, SQLite.SQLiteException) && return String[]
        rethrow(e)
    end
    return keys
end

"""
    sqlite_stats(S; location=:local)

Print a short summary of the SQLite cache database for symmetry `S`.

`location` may be `:local` for the current process-local database or `:global`
for the shared database.  The report lists non-empty table counts and the
database file size.  Calling this function opens the requested database if it is
not already open.
"""
function sqlite_stats(::Type{S}; location::Symbol=:local) where {S<:NonabelianSymm}
    db = get_sqlite_db(S, location)
    loc_str = location == :global ? "Global" : "Local ($(process_local_id()))"

    println("SQLite Statistics for $(totxt(S)) ($loc_str):")

    for table in SQLITE_TABLES
        try
            count_result = DBInterface.execute(db, "SELECT COUNT(*) as count FROM $(table)")
            entry_count = first(count_result).count
            entry_count > 0 && println("  $(table): $(entry_count) entries")
        catch
        end
    end

    path = db.file
    if isfile(path)
        size_mb = round(filesize(path) / 1024^2, digits=2)
        println("  File size: $(size_mb) MB")
    end
end

const _ALL_CACHES = (
    (IREP_CACHE, IREP_CACHE_LOCK),
    (CGT_CACHE, CGT_CACHE_LOCK),
    (FSYMBOL_CACHE, FSYMBOL_CACHE_LOCK),
    (RSYMBOL_CACHE, RSYMBOL_CACHE_LOCK),
    (XSYMBOL_CACHE, XSYMBOL_CACHE_LOCK),
    (OMLIST_CACHE, OMLIST_CACHE_LOCK),
    (VALIDOUT_CACHE, VALIDOUT_CACHE_LOCK),
    (CGTPERM_CACHE, CGTPERM_CACHE_LOCK),
    (CG3FLIP_CACHE, CG3FLIP_CACHE_LOCK),
)

"""
    clear_sqlite_cache!()

Clear all in-memory LRU caches used by the SQLite-backed object store.

This does not delete any SQLite database files and does not close database
connections.  It is useful when measuring memory use, forcing subsequent loads
to come from SQLite, or after manually changing database contents.
"""
function clear_sqlite_cache!()
    for (cache, cache_lock) in _ALL_CACHES
        Base.lock(cache_lock) do
            empty!(cache)
        end
    end
end

# Backward compatibility
clear_irep_sqlite_cache!() = clear_sqlite_cache!()

_sqlite_file_set(path::AbstractString) = (String(path), string(path, "-wal"), string(path, "-shm"))

function _sqlite_db_base_path(path::AbstractString)
    path_str = String(path)
    endswith(path_str, ".db") && return path_str
    endswith(path_str, ".db-wal") && return path_str[1:end-4]
    endswith(path_str, ".db-shm") && return path_str[1:end-4]
    return nothing
end

function _opened_local_sqlite_paths()
    opened_paths = Set{String}()
    Base.lock(SQLITE_DB_LOCK) do
        for (key, db) in SQLITE_DBS
            startswith(key, "local:") || continue
            push!(opened_paths, normpath(abspath(db.file)))
        end
    end
    return opened_paths
end

function _local_sqlite_db_paths(root::AbstractString)
    isdir(root) || return String[]

    paths = Set{String}()
    for (dir, _, files) in walkdir(root)
        for file in files
            base_path = _sqlite_db_base_path(joinpath(dir, file))
            base_path === nothing && continue
            push!(paths, normpath(abspath(base_path)))
        end
    end
    return sort!(collect(paths))
end

"""
    delete_closed_local_sqlite_dbs(; env=ENV, verbose=0)

Delete local SQLite database files that are not currently open in this Julia
process.

The function scans the configured local database root `sqlite_local_dir(env)`.
For each closed database, it removes the base `.db` file and any `-wal`/`-shm`
sidecar files.  Databases currently tracked by `get_sqlite_db(..., :local)` in
this process are skipped.

Returns the base `.db` paths that were removed.  This only considers local
databases; global databases are never deleted.
"""
function delete_closed_local_sqlite_dbs(; env::AbstractDict=ENV, verbose=0)
    root = sqlite_local_dir(env)
    opened_paths = _opened_local_sqlite_paths()
    deleted_paths = String[]

    for path in _local_sqlite_db_paths(root)
        path in opened_paths && continue

        deleted_any = false
        for file in _sqlite_file_set(path)
            isfile(file) || continue
            rm(file; force=true)
            deleted_any = true
        end

        if deleted_any
            push!(deleted_paths, path)
            verbose > 0 && println("Deleted local SQLite DB: $(path)")
        end
    end

    return deleted_paths
end

"""
    close_all_sqlite_dbs()

Close all SQLite database connections currently open in this Julia process.

This clears LurCGT's open-connection registry but does not clear the in-memory
object caches and does not delete any database files.  Use it before changing
SQLite-related environment variables in the same process, or before calling
`delete_closed_local_sqlite_dbs` when you want the current local databases to be
eligible for deletion.
"""
function close_all_sqlite_dbs()
    Base.lock(SQLITE_DB_LOCK) do
        for (_, db) in SQLITE_DBS
            try SQLite.close(db) catch end
        end
        empty!(SQLITE_DBS)
    end
end

function __init__()
    PROCESS_LOCAL_ID[] = build_process_local_id()
end
