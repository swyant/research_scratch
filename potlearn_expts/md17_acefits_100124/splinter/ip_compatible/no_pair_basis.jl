using ACEpotentials

function no_pair_basis(;kwargs...)
    kwargs = ACE1x._clean_args(kwargs)
    rpiB   = ACE1x.mb_ace_basis(kwargs)
    return rpiB
end

function ACE1x.algebraic_smoothness_prior(basis::ACE1.RPI.RPIBasis; p = 4, wL = 1.0) 
   d = Float64[]
   append!(d, ACE1.scaling(basis, p, wL))
   return Diagonal(1 .+ d)
end



