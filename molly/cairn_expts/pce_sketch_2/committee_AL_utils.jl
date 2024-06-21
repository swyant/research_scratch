struct CommitteePotential 
  members::Vector{AbstractPotential} #Should be AbstractPotential
  leader::Integer
end 

function get_all_energies(sys::AbstractSystem, cmte_pot::CommitteePotential)
  all_potengs = Vector{Float64}()

  for pot in cmte_pot.members
      poteng = potential_energy(sys,pot)
      push!(res,poteng) 
  end

  all_potengs
end

function get_all_forces(sys::AbstractSystem, cmte_pot::CommitteePotential)
  return
end

