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
choice = parse(Int64, ARGS[10]) # Choice for observables

mps_file_path = "/lustre/fs23/group/nic/tangelides/SW_MPS_iT/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r).h5"

path = "/lustre/fs23/group/nic/tangelides/SW_Observables_iT"

f_h5 = h5open(mps_file_path, "r")
psi = read(f_h5, "MPS", MPS)
close(f_h5)

if choice == 1

    text_file_name = "/Electric Field/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r).txt"
    path_to_text_file = path*text_file_name
    
    open(path_to_text_file, "w") do f
        sites = siteinds(psi)
        charge_configuration_list = get_SW_charge_configuration(psi, sites)
        electric_field_configuration_list = get_SW_electric_field_configuration(charge_configuration_list, l_0)
        left_edge = floor(Int, N*0.48)
        right_edge = floor(Int, N*0.52)
        middle_efl = electric_field_configuration_list[left_edge:right_edge]
        num_links = length(middle_efl)
        avg_E_field = real(mean(middle_efl))
        total_charge = sum(charge_configuration_list)
        write(f, "$(avg_E_field),$(num_links),$(total_charge)\n")
    end
end
