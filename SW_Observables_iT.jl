using ITensors.HDF5
using ITensors
using Statistics
include("Observables_iT.jl")

N = parse(Int, ARGS[1]) # Number of physical lattice sites
x = parse(Float64, ARGS[2]) # 1/(ag)^2
mg = parse(Float64, ARGS[3]) # m/g
D = parse(Int64, ARGS[4]) # Bond dimension of MPS to be computed
l_0 = parse(Float64, ARGS[5]) # l_0 = theta/(2*pi)
lambda = parse(Float64, ARGS[6]) # Lagrange multiplier to enforce total charge squared to be 0
acc = parse(Float64, ARGS[7]) # Tolerance for stopping condition of the variational algorithm
ns = parse(Int64, ARGS[8]) # Maximum number of sweeps of the variational algorithm
r = parse(Float64, ARGS[9]) # Wilson parameter

mps_file_path = "/lustre/fs23/group/nic/tangelides/SW_MPS_States_iT/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r).h5"

path = "/lustre/fs23/group/nic/tangelides/"

f_h5 = h5open(mps_file_path, "r")
psi = read(f_h5, "MPS", MPS)
gs_energy = read(f_h5, "gs_energy")
close(f_h5)

text_file_name = "/SW_Observables_iT/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r).txt"
path_to_text_file = path*text_file_name

open(path_to_text_file, "w") do f
    
    sites = siteinds(psi)
    
    z_configuration_list = get_z_configuration(psi, sites)
    
    charge_configuration_list = get_SW_charge_configuration(z_configuration_list)
    
    total_charge = sum(charge_configuration_list)

    electric_field_configuration_list = get_SW_electric_field_configuration(charge_configuration_list, l_0)
    
    left_edge = floor(Int, N*0.48)
    
    right_edge = floor(Int, N*0.52)
    
    middle_efl = electric_field_configuration_list[left_edge:right_edge]
    
    num_links = length(middle_efl)
    
    avg_E_field = real(mean(middle_efl))

    ee = get_SW_entanglement_entropy(psi)

    cc = get_SW_chiral_condensate(psi, params)

    write(f, "Energy, Average Electric Field, Number of Links, Total Charge, Chiral Condensate, Entanglement Entropy (First Line)\n")

    write(f, "z configuration list (Second Line)\n")

    write(f, "$(gs_energy),$(avg_E_field),$(num_links),$(total_charge),$(cc),$(ee)\n")

    l_tmp = length(z_configuration_list)

    for (index, element) in enumarate(z_configuration_list)
        if index != l_tmp
            write(f, "$(element),")
        else
            write(f, "$(element)")
        end
    end

end

