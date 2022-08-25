include("run_DMRG_iT.jl")
include("Observables_iT.jl")

# ---------------------------------------------------------------------------------

# Testing the first excited state energy using exact diagonalization for Schwinger Wilson

# N = 4
# x = 1.0
# D = 80
# mg = 0.125
# l_0 = 0.125
# lambda = 100.0
# acc = 1e-12
# ns = 1000
# r = 1.0
# D_p = 0
# mg_p = 0.0
# silent = true
# initial_noise = 0.0

# sweep_observables_file_path = "N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r).txt"

# mps_file_path = "N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r).h5"

# previous_mps_file_path = "N_$(N)_x_$(x)_D_$(D_p)_l0_$(l_0)_mg_$(mg_p)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r).h5"

# if isfile(mps_file_path)
#     f = h5open(mps_file_path, "r")
#     previous_psi = read(f, "MPS", MPS)
#     sites = siteinds(previous_psi)
#     close(f)
# elseif isfile(previous_mps_file_path)
#     f = h5open(previous_mps_file_path, "r")
#     previous_psi = read(f, "MPS", MPS)
#     sites = siteinds(previous_psi)
#     close(f)
# else
#     sites = siteinds("S=1/2", 2*N)
#     previous_psi = randomMPS(sites, D)
# end

# # params = Dict("initial_noise" => initial_noise, "silent" => silent, "N" => N, "D" => D, "x" => x, "ns" => ns, "lambda" => lambda, "l_0" => l_0, "mg" => mg, "r" => r, "acc" => acc, "sweep_observables_file_path" => sweep_observables_file_path, "previous_mps_file_path" => previous_mps_file_path, "previous_psi" => previous_psi)

# # Compute the MPS

# # H = get_MPO_from_OpSum(get_SW_OpSum(params), sites)

# params = Dict("initial_noise" => initial_noise, "silent" => silent, "N" => 10, "J" => -1, "g_z" => -0.1, "g_x" => 1.5, "ns" => 5, "D" => 34, "acc" => acc, "sweep_observables_file_path" => sweep_observables_file_path, "previous_mps_file_path" => previous_mps_file_path, "previous_psi" => previous_psi)
# J = -1.0
# N = 10
# g_z = -0.1
# g_x = 1.5
# ns = 1000
# D = 80
# opsum_ising = get_Ising_OpSum(N, J, g_z, g_x)
# sites = siteinds("S=1/2", N)
# H = get_MPO_from_OpSum(opsum_ising, sites)

# energy_0, psi_0 = run_SW_DMRG(sites, params, H, true)

# # z_config = get_z_configuration(psi_0, sites)
# # charge_config = get_SW_charge_configuration(z_config)
# # electric_field_config = get_SW_electric_field_configuration(charge_config, l_0)
# # left_edge = floor(Int, N*0.48)
# # right_edge = floor(Int, N*0.52)
# # middle_efl = electric_field_config[left_edge:right_edge]
# # num_links = length(middle_efl)
# # avg_E_field = real(mean(middle_efl))
# # println(avg_E_field)

# initial_noise = 1e-5
# previous_psi = randomMPS(sites, D)
# # params = Dict("initial_noise" => initial_noise, "silent" => silent, "N" => N, "D" => D, "x" => x, "ns" => ns, "lambda" => lambda, "l_0" => l_0, "mg" => mg, "r" => r, "acc" => acc, "sweep_observables_file_path" => sweep_observables_file_path, "previous_mps_file_path" => previous_mps_file_path, "previous_psi" => previous_psi)

# params = Dict("initial_noise" => initial_noise, "silent" => silent, "N" => 10, "J" => -1, "g_z" => -0.1, "g_x" => 1.5, "ns" => 5, "D" => 34, "acc" => acc, "sweep_observables_file_path" => sweep_observables_file_path, "previous_mps_file_path" => previous_mps_file_path, "previous_psi" => previous_psi)

# # H = get_MPO_from_OpSum(get_SW_OpSum(params), sites)
# P = outer(psi_0', psi_0)
# Heff = H + energy_0.*P

# energy_1, psi_1 = run_SW_DMRG(sites, params, Heff, true)

# println(energy_0)
# println(energy_1)

# ---------------------------------------------------------------------------------

# Testing the first excited state energy using exact diagonalization for Ising model

# N = 20
# J = 0.0001
# g_z = 0.1
# g_x = 0.2
# ns = 500
# D = 20

# sites = siteinds("S=1/2", N)
# opsum_ising = get_Ising_OpSum(N, J, g_z, g_x)
# H = get_MPO_from_OpSum(opsum_ising, sites)
# initial_ansatz_0 = randomMPS(sites, D)
# sweeps = Sweeps(ns, maxdim = D)

# energy_0, psi_0 = dmrg(H, initial_ansatz_0, sweeps, ishermitian = true, maxdim = D)

# # Ms = [psi_0]
# # w = energy_0
# # initial_ansatz_1 = randomMPS(sites, D)

# # energy_1, psi_1 = dmrg(H, Ms, initial_ansatz_1, sweeps, weight = w, ishermitian = true, maxdim = D)

# Heff = H + energy_0.*outer(psi_0', psi_0)
# initial_ansatz_1 = randomMPS(sites, D)

# println("Overlap of ansatz with gs: ", inner(psi_0, initial_ansatz_1))

# initial_noise = 1e-2

# noise_vector = LinRange(initial_noise, 0.0, ns)

# energy_1, psi_1 = dmrg(Heff, initial_ansatz_1, sweeps, ishermitian = true, maxdim = D, noise = noise_vector)

# println("Overlap of 1st excited state with gs: ", inner(psi_0, psi_1))

# println(energy_0)
# println(energy_1)

# ---------------------------------------------------------------------------------

# Testing the first excited state of the Schwinger model

N = 6
x = 1.0
D = 80
mg = -0.2
l_0 = 0.125
lambda = 100.0
r = 1.0
ns = 500

params = Dict("N" => N, "l_0" => l_0, "N" => N, "x" => x, "mg" => mg, "r" => r, "lambda" => lambda)

sites = siteinds("S=1/2", 2*N)
H = get_MPO_from_OpSum(get_SW_OpSum(params), sites)
initial_ansatz_0 = randomMPS(sites, D)
sweeps = Sweeps(ns, maxdim = D)

energy_0, psi_0 = dmrg(H, initial_ansatz_0, sweeps, ishermitian = true, maxdim = D)

println("Norm of psi_0: ", norm(psi_0))

Heff = H + energy_0.*outer(psi_0', psi_0)
initial_ansatz_1 = randomMPS(sites, D)

println("Overlap of ansatz with gs: ", inner(psi_0, initial_ansatz_1))

initial_noise = 1e-2

noise_vector = LinRange(initial_noise, 0.0, ns)

# energy_1, psi_1 = dmrg(Heff, initial_ansatz_1, sweeps, ishermitian = true, maxdim = D, noise = noise_vector)

Ms = [psi_0]
w = 1.0
initial_ansatz_1 = randomMPS(sites, D)

energy_1, psi_1 = dmrg(Heff, Ms, initial_ansatz_1, sweeps, weight = w, ishermitian = true, maxdim = D)

println("Overlap of 1st excited state with gs: ", inner(psi_0, psi_1))

println(energy_0)
println(energy_1)

# ---------------------------------------------------------------------------------
