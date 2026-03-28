using Test
using LurCGT

@testset "SQLite local id builder" begin
    @test LurCGT.build_process_local_id(Dict{String, String}(); hostname="nodeA", pid=4321) ==
          "nodeA_pid4321"
    @test LurCGT.process_local_id(Dict{String, String}(); hostname="nodeA", pid=4321) ==
          "nodeA_pid4321"

    @test LurCGT.build_process_local_id(
        Dict("SLURM_JOB_ID" => "1001", "SLURM_ARRAY_TASK_ID" => "7");
        hostname="nodeB",
        pid=987,
    ) == "nodeB_job1001_task7_pid987"
    @test LurCGT.process_local_id(
        Dict("SLURM_JOB_ID" => "1001", "SLURM_ARRAY_TASK_ID" => "7");
        hostname="nodeB",
        pid=987,
    ) == "nodeB_job1001_task7_pid987"

    @test LurCGT.build_process_local_id(
        Dict("SLURM_JOB_ID" => "2002", "SLURM_PROCID" => "3");
        hostname="nodeC",
        pid=654,
    ) == "nodeC_job2002_proc3_pid654"
    @test LurCGT.process_local_id(
        Dict("SLURM_JOB_ID" => "2002", "SLURM_PROCID" => "3");
        hostname="nodeC",
        pid=654,
    ) == "nodeC_job2002_proc3_pid654"
end
