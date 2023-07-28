using KrylovKit
using ITensors.HDF5
using Plots
include("DMRG_iT.jl")
include("Observables_iT.jl")

function run_SW_DMRG(sites, params, H, ishermitian)::Tuple{Float64, MPS}

    energy::Float64, psi::MPS = DMRG(H, sites, params, ishermitian)

    return energy, psi

end
