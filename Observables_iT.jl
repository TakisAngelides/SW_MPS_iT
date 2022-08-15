include("MPO_iT.jl")

function get_SW_charge_configuration(psi, sites)

    N_spin = length(psi)
    N = Int(N_spin/2)

    charge_configuration_list = []
    for i=1:N
        Q_local_OpSum = get_SW_local_charge_OpSum(i)
        Q_local_MPO = get_MPO_from_OpSum(Q_local_OpSum, sites)
        append!(charge_configuration_list, inner(psi', Q_local_MPO, psi))
    end

    return charge_configuration_list

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
