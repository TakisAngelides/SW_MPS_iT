using ITensors
using Statistics
include("MPO_iT.jl")

mutable struct my_observer <: AbstractObserver
    energy_tol::Float64
    last_energy::Float64
    sites::Vector{Index{Int64}}
    params::Dict
end

function ITensors.measure!(o::my_observer; kwargs...)

    sweep = kwargs[:sweep]
    psi = kwargs[:psi]
    energy = kwargs[:energy]
    sweep_is_done = kwargs[:sweep_is_done]
    sites = o.sites
    params = o.params
    l_0 = params["l_0"]
    N = params["N"]
    sweep_observables_file_path = params["sweep_observables_file_path"]

    if sweep_is_done
        open(sweep_observables_file_path, "a+") do f
            charge_configuration_list = get_SW_charge_configuration(psi, sites)
            electric_field_configuration_list = get_SW_electric_field_configuration(charge_configuration_list, l_0)
            left_edge = floor(Int, N*0.48)
            right_edge = floor(Int, N*0.52)
            middle_efl = electric_field_configuration_list[left_edge:right_edge]
            avg_E_field = real(mean(middle_efl))
            total_charge = sum(charge_configuration_list)
            write(f, "$(sweep),$(total_charge),$(energy),$(avg_E_field)\n")
        end
    end

end

function DMRG(H, sites, params)

    """
    ns = number of sweeps
    """ 

    D = params["D"]
    ns = params["ns"]
    acc = params["acc"]
    
    noise_vector = zeros(ns)
    psi_0 = randomMPS(sites, D)
    sweeps = Sweeps(ns)
    observer = my_observer(acc, 1000.0, sites, params)
    
    energy, psi = dmrg(H, psi_0, sweeps, ishermitian = true, noise = noise_vector, observer = observer)

    return energy, psi

end