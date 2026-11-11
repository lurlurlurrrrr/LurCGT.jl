"""
    validate_fundamental_ireps(S, RT=BigInt; verbose=0) -> Int

Request and save all fundamentals using recursive recipe dependencies, collecting all
incidental nonfundamental copies outside the cache. After every fundamental is
saved, compare each candidate with its canonical construction. Return the count
of validated copies; throw on the first unequal field set. Run with a fresh
database so canonical references cannot come from a previous backend.
"""
function validate_fundamental_ireps(::Type{S}, ::Type{RT}=BigInt; verbose=0) where {S,RT<:Integer}
    initialize_symmetry(S,RT;verbose)
    candidates = Irep{S,nlops(S),nzops(S),RT}[]
    for k in eachindex(fundamental_recipes(S))
        rep = construct_fundamental(S,RT,k; candidates,verbose)
        save_irep_sqlite(rep;verbose)
    end
    for candidate in candidates
        canonical = getNsave_irep(S,RT,candidate.qlabel;verbose)
        same_irep_fields(candidate,canonical) || error("Canonical Irep mismatch for $S $(candidate.qlabel)")
    end
    return length(candidates)
end

"""
    validate_fundamental_products(S, RT=BigInt; rules=irrep_validation_rules(S)) -> Int

Execute the supplied finite cross-product rules after fundamental construction.
Compare each independently extracted target with its canonical cached irrep;
return the number of passing target copies across all rules. Large exceptional rules are opt-in work.
"""
function validate_fundamental_products(::Type{S}, ::Type{RT}=BigInt;
        rules=irrep_validation_rules(S), verbose=0) where {S,RT<:Integer}
    checked = 0
    for rule in rules
        target = _recipe_label(S,rule.target)
        left,right = map(q -> _recipe_label(S,q),rule.sources)
        x = getNsave_irep(S,RT,left;verbose); y = getNsave_irep(S,RT,right;verbose)
        Sl,Sz,_,wi = tensor_prod_reps(x,y;verbose)
        copies = Irep{S,nlops(S),nzops(S),RT}[]
        get_irep_frommw(S,Sl,Sz,x.innerprod,y.innerprod,qlab2mwz(S,target),wi;copies,verbose)
        canonical = getNsave_irep(S,RT,target;verbose)
        for result in copies
            same_irep_fields(result,canonical) || error("Cross-product mismatch for $S $target from $left × $right")
            checked += 1
        end
    end
    return checked
end
