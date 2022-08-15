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
    silent = params["silent"]
    l_0 = params["l_0"]
    N = params["N"]
    sweep_observables_file_path = params["sweep_observables_file_path"]

    if sweep_is_done && !silent
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
    previous_mps_file_path = params["previous_mps_file_path"]
    initial_noise = params["initial_noise"]

    if isfile(previous_mps_file_path)
        f = h5open(previous_mps_file_path,"r")
        psi_0 = read(f, "MPS", MPS)
        close(f)
    else
        psi_0 = randomMPS(sites, D)
    end

    noise_vector = LinRange(initial_noise, 0.0, ns) # Noise to be added to the MPS during DMRG
    sweeps = Sweeps(ns) # This is the maximum number of sweeps to be done if accuracy (acc) is not reached
    observer = my_observer(acc, 1000.0, sites, params) # 1000.0 is for initial energy of the algorithm and should be set well above the estimated g.s. energy
    
    energy, psi = dmrg(H, psi_0, sweeps, ishermitian = true, noise = noise_vector, observer = observer, maxdim = D, write_when_maxdim_exceeds = 80)

    return energy, psi

end