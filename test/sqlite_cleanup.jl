using LurCGT
using Test

@testset "SQLite closed local DB cleanup" begin
    mktempdir() do tmp
        local_dir = joinpath(tmp, "local-db")
        global_dir = joinpath(tmp, "global-db")
        closed_path = joinpath(local_dir, LurCGT.totxt(SU{3}), "closed.db")

        withenv("LURCGT_LOCALDB_DIR" => local_dir,
                "LURCGT_GLOBALDB_DIR" => global_dir) do
            try
                LurCGT.close_all_sqlite_dbs()
                open_db = LurCGT.get_sqlite_db(SU{2}, :local)

                mkpath(dirname(closed_path))
                for file in (closed_path, string(closed_path, "-wal"), string(closed_path, "-shm"))
                    write(file, UInt8[])
                end

                deleted = LurCGT.delete_closed_local_sqlite_dbs()

                @test normpath(abspath(closed_path)) in deleted
                @test isfile(open_db.file)
                @test !isfile(closed_path)
                @test !isfile(string(closed_path, "-wal"))
                @test !isfile(string(closed_path, "-shm"))
            finally
                LurCGT.close_all_sqlite_dbs()
            end
        end
    end
end
