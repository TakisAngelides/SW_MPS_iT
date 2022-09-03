using ITensors
using Statistics
include("MPO_iT.jl")

mutable struct my_observer <: AbstractObserver
    energy_tol::Float64
    last_energy::Float64
    sites::Vector{Index{Int64}}
    params::Dict
end

function ITensors.checkdone!(o::my_observer;kwargs...)::Bool
    
    sw::Int64 = kwargs[:sweep]
    energy::Float64 = kwargs[:energy]
    
    if abs(energy-o.last_energy)/abs(energy) < o.energy_tol
        
        println("Energy tolerance reached, stopping DMRG after sweep $sw")
        
        return true
    
    end
    
    # Otherwise, update last_energy and keep going
    
    o.last_energy = energy
    
    return false
end

function ITensors.measure!(o::my_observer; kwargs...)

    sweep::Int64 = kwargs[:sweep]
    psi::MPS = kwargs[:psi]
    energy::Float64 = kwargs[:energy]
    sweep_is_done::Bool = kwargs[:sweep_is_done]
    sites = o.sites
    params = o.params
    silent::Bool = params["silent"]
    l_0::Float64 = params["l_0"]
    N::Int64 = params["N"]
    sweep_observables_file_path = params["sweep_observables_file_path"]

    if sweep_is_done && !silent
        
        open(sweep_observables_file_path, "a+") do f
            
            z_config_list::Vector{ComplexF64} = get_z_configuration(psi, sites)
            
            charge_configuration_list::Vector{ComplexF64} = get_SW_charge_configuration(z_config_list)
            
            electric_field_configuration_list::Vector{ComplexF64} = get_SW_electric_field_configuration(charge_configuration_list, l_0)
            
            left_edge::Int64 = floor(Int, N*0.48)
            
            right_edge::Int64 = floor(Int, N*0.52)
            
            middle_efl::Vector{ComplexF64} = electric_field_configuration_list[left_edge:right_edge]
            
            avg_E_field::Float64 = real(mean(middle_efl))
            
            total_charge::ComplexF64 = sum(charge_configuration_list)
            
            write(f, "$(sweep),$(total_charge),$(energy),$(avg_E_field)\n")
        end
    end

end

function DMRG(H, sites::Vector{Index{Int64}}, params::Dict, ishermitian::Bool)::Tuple{Float64, MPS}

    """
    ns = number of sweeps
    """ 

    D::Int64 = params["D"]
    ns::Int64 = params["ns"]
    acc::Float64 = params["acc"]
    initial_noise::Float64 = params["initial_noise"]
    initial_ansatz::MPS = params["previous_psi"]
    first_excited::MPS = params["first_excited"]

    noise_vector = LinRange(initial_noise, 0.0, ns) # Noise to be added to the MPS during DMRG
    sweeps = Sweeps(ns, maxdim = D) # This is the maximum number of sweeps to be done if accuracy (acc) is not reached
    observer = my_observer(acc, 1000.0, sites, params) # 1000.0 is for initial energy of the algorithm and should be set well above the estimated g.s. energy
    
    if !first_excited
        energy::Float64, psi::MPS = dmrg(H, initial_ansatz, sweeps, ishermitian = ishermitian, noise = noise_vector, observer = observer, maxdim = D)
    else
        Ms::Vector{MPS} = params["Ms"]
        w::Float64 = params["w"]
        energy, psi = dmrg(H, Ms, initial_ansatz, sweeps, weight = w, ishermitian = ishermitian, noise = noise_vector, observer = observer, maxdim = D)
    end

    return energy, psi

end