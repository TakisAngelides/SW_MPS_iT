include("MPO_iT.jl")
using ITensors

function get_SW_charge_configuration(z_configuration_list)

    N_spin = length(z_configuration_list)
    N = Int(N_spin/2)

    charge_configuration_list = []
    for k=1:N
        charge_on_site_k = 0.5*(z_configuration_list[2*k-1] + z_configuration_list[2*k])
        append!(charge_configuration_list, charge_on_site_k)
    end

    return charge_configuration_list

end

function get_z_configuration(psi, sites)

    N_spin = length(psi)

    z_configuration_list = []
    for i=1:N_spin
        z_OpSum = get_SW_local_z_OpSum(i)
        z_mpo = get_MPO_from_OpSum(z_OpSum, sites)
        append!(z_configuration_list, inner(psi', z_mpo, psi))
    end

    return z_configuration_list

end

function get_SW_electric_field_configuration(charge_configuration_list, l_0)

    N = length(charge_configuration_list)

    electric_field_configuration_list = []

    for i in 1:N-1
        L_n = l_0 + sum(charge_configuration_list[1:i])
        append!(electric_field_configuration_list, L_n)
    end

    return electric_field_configuration_list

end

function get_SW_entanglement_entropy(psi)

    N_spin = length(psi)

    half_of_spin_chain = Int(N_spin/2)

    sites = siteinds(psi)

    orthogonalize!(psi, half_of_spin_chain)

    _, S = svd(psi[half_of_spin_chain], (linkind(psi, half_of_spin_chain-1), sites[half_of_spin_chain]))

    SvN = 0.0 # von Neumann entropy

    for n in 1:dim(S, 1)
        p = S[n,n]^2
        SvN -= p*log(p)
    end

    return SvN

end

function get_SW_chiral_condensate(psi, params)

    cc_opsum = get_SW_chiral_condensate_OpSum(params)
    sites = siteinds(psi)
    cc_mpo = get_MPO_from_OpSum(cc_opsum, sites)
    
    return inner(psi', cc_mpo, psi)

end