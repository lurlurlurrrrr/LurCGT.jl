# This function should be tested when there is no HDF5 file existing.
# Test HDF5 I/O functionality
function test_hdf5_io(::Type{S}) where {S<:NonabelianSymm}
    NZ = nzops(S)
    
    # Test Irep save/load
    println("Testing Irep save/load...")
    test_qlabel = ntuple(i -> i == 1 ? 2 : (i == 2 && NZ >= 2 ? 1 : 0), NZ)
    # Use getNsave to generate and save
    irep_orig = getNsave_irep(S, BigInt, test_qlabel)
    
    # Load directly to verify save worked
    irep_loaded = getNsave_irep(S, BigInt, test_qlabel)
    @assert !isnothing(irep_loaded)
    @assert irep_loaded.qlabel == irep_orig.qlabel
    @assert irep_loaded.dimension == irep_orig.dimension
    @assert length(irep_loaded.Sz) == length(irep_orig.Sz)
    println("  Irep save/load: PASSED")
    
    # Test CG3 save/load
    println("Testing CG3 save/load...")
    q1 = ntuple(i -> i == 2 && NZ >= 2 ? 1 : (i == 1 && NZ == 1 ? 1 : 0), NZ)
    q2 = ntuple(i -> i == 1 ? 1 : 0, NZ)
    # Use getNsave_cg3 from clebsch_io which handles generation
    LurCGT.generate_every_CGT(S, BigInt, BigInt, (q1, q2), nothing; verbose=0)
    
    # Get valid outputs using getNsave
    vo = getNsave_validout(S, (q1, q2))
    @assert !isnothing(vo) && length(vo.out_spaces) > 0
    possible_out = vo.out_spaces[1]
    
    # Use getNsave_cg3 to verify data accessibility
    cg3s_loaded = getNsave_cg3(S, BigInt, (q1, q2), [possible_out])
    @assert !isnothing(cg3s_loaded)
    @assert !isempty(cg3s_loaded[possible_out].blocks)
    @assert length(cg3s_loaded[possible_out].nfactor) > 0
    println("  CG3 save/load: PASSED")
    
    # Test F-symbol save/load
    println("Testing F-symbol save/load...")
    in1, in2, in3 = q1, q2, q1
    # Find a valid output for F-symbol
    vo_f = getNsave_validout(S, sort((in1, in2, in3)))
    @assert !isnothing(vo_f) && length(vo_f.out_spaces) > 0
    out = vo_f.out_spaces[1]
    
    # Use getNsave to generate and save
    fsym = getNsave_Fsymbol(S, BigInt, in1, in2, in3, out; verbose=0)
    @assert !isnothing(fsym)
    @assert fsym.in1 == in1
    @assert fsym.in2 == in2
    @assert fsym.in3 == in3
    @assert fsym.out == out
    
    # Call again to verify load works
    fsym_loaded = getNsave_Fsymbol(S, BigInt, in1, in2, in3, out; verbose=0)
    @assert !isnothing(fsym_loaded)
    @assert size(fsym_loaded.fsym_mat) == size(fsym.fsym_mat)
    println("  F-symbol save/load: PASSED")
    
    # Test R-symbol save/load
    println("Testing R-symbol save/load...")
    in_r = q1
    out_r = q1 .* 2
    # Use getNsave to generate and save
    rsym = getNsave_Rsymbol(S, BigInt, in_r, out_r; verbose=0)
    @assert !isnothing(rsym)
    @assert rsym.in == in_r
    @assert rsym.out == out_r
    
    # Call again to verify load works
    rsym_loaded = getNsave_Rsymbol(S, BigInt, in_r, out_r; verbose=0)
    @assert !isnothing(rsym_loaded)
    @assert size(rsym_loaded.rsym_mat) == size(rsym.rsym_mat)
    println("  R-symbol save/load: PASSED")
    
    # Test OMList save/load
    println("Testing OMList save/load...")
    incom_test = (q1, q2)
    out_test = possible_out
    # Use getNsave to generate and save
    omlist = getNsave_omlist(S, incom_test, out_test)
    @assert !isnothing(omlist)
    @assert omlist.incom_spaces == incom_test
    @assert omlist.out_space == out_test
    
    # Call again to verify load works
    omlist_loaded = getNsave_omlist(S, incom_test, out_test)
    @assert !isnothing(omlist_loaded)
    @assert omlist_loaded.totalOM == omlist.totalOM
    println("  OMList save/load: PASSED")
    
    # Test ValidOuts save/load
    println("Testing ValidOuts save/load...")
    # Use getNsave to generate and save
    vo = getNsave_validout(S, incom_test)
    @assert !isnothing(vo)
    @assert vo.incom_spaces == incom_test
    
    # Call again to verify load works
    vo_loaded = getNsave_validout(S, incom_test)
    @assert !isnothing(vo_loaded)
    @assert length(vo_loaded.out_spaces) == length(vo.out_spaces)
    @assert vo_loaded.out_spaces == vo.out_spaces
    println("  ValidOuts save/load: PASSED")
    
    println("All HDF5 I/O tests passed!")
end

# Test thread-safe HDF5 access
function test_hdf5_threadsafety(::Type{S}) where {S<:NonabelianSymm}
    NZ = nzops(S)
    println("Testing thread-safe concurrent access...")
    
    # Generate some test data first using getNsave
    test_qlabels = [ntuple(i -> i <= j ? 1 : 0, NZ) for j in 0:min(3, NZ)]
    for q in test_qlabels
        getNsave_irep(S, BigInt, q)
    end
    
    # Concurrent reads using getNsave (which will load from cache)
    results = Vector{Any}(undef, length(test_qlabels))
    println(Threads.nthreads(), " threads will be used for concurrent reads.")
    Threads.@threads for i in 1:length(test_qlabels)
        results[i] = getNsave_irep(S, BigInt, test_qlabels[i])
    end
    
    # Verify all reads succeeded
    for (i, irep) in enumerate(results)
        @assert !isnothing(irep)
        @assert irep.qlabel == test_qlabels[i]
    end
    
    println("  Thread-safe concurrent reads: PASSED")
    println("Thread-safety test completed!")
end

