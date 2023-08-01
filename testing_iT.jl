# include("run_DMRG_iT.jl")
# include("Observables_iT.jl")
using ITensors
using Plots
# using DataStructures
# using StatsBase
include("MPO_iT.jl")

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

# # Testing the first excited state of the Schwinger model with Wilson fermions

# N = 20
# x = 1.0
# D = 10
# mg = -0.1
# l_0 = 0.01
# lambda = 10.0
# r = 1.0
# ns = 50

# params = Dict("N" => N, "l_0" => l_0, "N" => N, "x" => x, "mg" => mg, "r" => r, "lambda" => lambda)

# sites = siteinds("S=1/2", 2*N)
# opsum = get_SW_OpSum(params)
# H = get_MPO_from_OpSum(opsum, sites)
# initial_ansatz_0 = randomMPS(sites, D)
# sweeps = Sweeps(ns, maxdim = D)

# energy, psi = dmrg(H, initial_ansatz_0, sweeps, ishermitian = true, maxdim = D)

# z_configuration_list = get_z_configuration(psi, sites)
    
# charge_configuration_list = get_SW_charge_configuration(z_configuration_list)

# total_charge = sum(charge_configuration_list)

# electric_field_configuration_list = get_SW_electric_field_configuration(charge_configuration_list, l_0)

# left_edge = floor(Int, N*0.48)

# right_edge = floor(Int, N*0.52)

# middle_efl = electric_field_configuration_list[left_edge:right_edge]

# num_links = length(middle_efl)

# avg_E_field = real(mean(middle_efl))

# ee = get_SW_entanglement_entropy(psi)

# cc = get_SW_chiral_condensate(psi)

# println("Norm of psi_0: ", norm(psi_0))

# Ms = [psi_0]
# w = abs(energy_0)
# initial_ansatz_1 = randomMPS(sites, D)

# println(typeof(Ms))
# println(typeof(w))

# println("Overlap of ansatz with gs: ", inner(psi_0, initial_ansatz_1))

# energy_1, psi_1 = dmrg(H, Ms, initial_ansatz_1, sweeps, weight = w, ishermitian = true, maxdim = D)

# println("Overlap of 1st excited state with gs: ", inner(psi_0, psi_1))

# println(energy_0)
# println(energy_1)

# ---------------------------------------------------------------------------------

# Testing the staggered fermions Schwinger model

# N = 10
# x = 1.0
# D = 20
# mg = -0.1
# l_0 = 0.01
# lambda = 100.0
# r = 1.0
# ns = 1000

# params = Dict("N" => N, "l_0" => l_0, "N" => N, "x" => x, "mg" => mg, "r" => r, "lambda" => lambda)

# sites = siteinds("S=1/2", N)
# opsum = get_Staggered_OpSum(params)
# H = get_MPO_from_OpSum(opsum, sites)
# initial_ansatz_0 = randomMPS(sites, D)
# sweeps = Sweeps(ns, maxdim = D)

# energy_0, psi_0 = dmrg(H, initial_ansatz_0, sweeps, ishermitian = true, maxdim = D)

# println(energy_0)

# ---------------------------------------------------------------------------------

# Finding the mass shift at N = 4, x = 0.01: 
# theta = pi/8, r = 1 and theta = pi/8 + pi, r = -1

# N = 20
# x = 1.0
# D = 60
# l_0 = (pi/8)/(2*pi)
# lambda = 10.0
# r = 1.0
# ns = 50
# mg_list = LinRange(-2.0, -1.0, 10)
# avg_e_field_list = []
# file = open("file.txt", "w")
# for mg in mg_list
#     params = Dict("N" => N, "l_0" => l_0, "N" => N, "x" => x, "mg" => mg, "r" => r, "lambda" => lambda)
#     sites = siteinds("S=1/2", 2*N)
#     opsum = get_SW_OpSum(params)
#     H = get_MPO_from_OpSum(opsum, sites)
#     initial_ansatz_0 = randomMPS(sites, D)
#     sweeps = Sweeps(ns, maxdim = D)
#     energy, psi = dmrg(H, initial_ansatz_0, sweeps, ishermitian = true, maxdim = D)
#     z_configuration_list = get_z_configuration(psi, sites)
#     charge_configuration_list = get_SW_charge_configuration(z_configuration_list)
#     total_charge = sum(charge_configuration_list)
#     electric_field_configuration_list = get_SW_electric_field_configuration(charge_configuration_list, l_0)
#     display(real(electric_field_configuration_list))
#     println("The mass is: $(mg)")
#     avg_E_field = real(electric_field_configuration_list[9])
#     write(file, "$(mg),$(avg_E_field)\n")
# end
# close(file)

# N = 20
# x = 1.0
# D = 80
# l_0 = (pi/8+pi)/(2*pi)
# lambda = 100.0
# r = -1.0
# ns = 200
# mg_list = LinRange(-5.0, 5.0, 50)
# avg_e_field_list = []
# file = open("file_1.txt", "w")
# for mg in mg_list
#     params = Dict("N" => N, "l_0" => l_0, "N" => N, "x" => x, "mg" => mg, "r" => r, "lambda" => lambda)
#     sites = siteinds("S=1/2", 2*N)
#     opsum = get_SW_OpSum(params)
#     H = get_MPO_from_OpSum(opsum, sites)
#     initial_ansatz_0 = randomMPS(sites, D)
#     sweeps = Sweeps(ns, maxdim = D)
#     energy, psi = dmrg(H, initial_ansatz_0, sweeps, ishermitian = true, maxdim = D)
#     z_configuration_list = get_z_configuration(psi, sites)
#     charge_configuration_list = get_SW_charge_configuration(z_configuration_list)
#     total_charge = sum(charge_configuration_list)
#     electric_field_configuration_list = get_SW_electric_field_configuration(charge_configuration_list, l_0)
#     display(real(electric_field_configuration_list))
#     println("The mass is: $(mg)")
#     avg_E_field = real(electric_field_configuration_list[9])
#     write(file, "$(mg),$(avg_E_field)\n")
# end
# close(file)

# ---------------------------------------------------------------------------------

# # Testing the particle number operator

# N = 10
# x = 1.0
# D = 10
# l_0 = 0.125
# lambda = 100.0
# r = 1.0
# ns = 200
# # mg_list = LinRange(-5.0, 5.0, 1)
# mg_list = [-0.1]
# avg_e_field_list = []

# for mg in mg_list
#     params = Dict("N" => N, "l_0" => l_0, "N" => N, "x" => x, "mg" => mg, "r" => r, "lambda" => lambda)
#     sites = siteinds("S=1/2", 2*N)
#     opsum = get_SW_OpSum(params)
#     H = get_MPO_from_OpSum(opsum, sites)
#     initial_ansatz_0 = randomMPS(sites, D)
#     sweeps = Sweeps(ns, maxdim = D)
#     energy, psi = dmrg(H, initial_ansatz_0, sweeps, ishermitian = true, maxdim = D)
#     pn = get_particle_number(psi)
#     println(pn)
# end

# ---------------------------------------------------------------------------------

# s = siteinds("Qubit", 4)
# os = [("H", 1), ("X", 2), ("CNOT", 2, 4)]
# gates = ops(os, s)
# print(print(gates[1]))

# N = 2
# δτ = 0.1
# s = siteinds("S=1/2", N; conserve_qns=true)
# function ITensors.op(::OpName"expτSS", ::SiteType"S=1/2", s1::Index, s2::Index; τ)
#     h = -4 * op("Sz", s1) * op("Sz", s2)
#     return exp(τ * h)
# end
# gates = ops([("expτSS", (n, n + 1), (τ=-δτ / 2,)) for n in 1:(N - 1)], s)
# print(gates)

# ---------------------------------------------------------------------------------

# function ITensors.op(::OpName"expτSS", ::SiteType"S=1/2", s1::Index, s2::Index; τ)
#     h = -4 * op("Sz", s1) * op("Sz", s2)
#     return exp(τ * h)
#   end
  
#   function get_rho(N; cutoff = 1E-8, δτ = 0.1, beta_max = 1.0)

#     # Make an array of 'site' indices
#     s = siteinds("S=1/2", N)

#     # Make gates (1,2),(2,3),(3,4),...
#     gates = ops([("expτSS", (n, m), (τ=-δτ / 2,)) for n in 1:N for m in n+1:N], s)
#     # Include gates in reverse order to complete Trotter formula
#     append!(gates, reverse(gates))

#     # Initial state is infinite-temperature mixed state
#     rho = MPO(s, "Id") ./ √2

#     # Make H for measuring the energy        
#     function get_Hamiltonian(N)
#         H = OpSum()
#         for n=1:N-1
#             H += -4,"Sz",n,"Sz",n+1
#         end
#         return H
#     end
#     H = get_Hamiltonian(N)
#     H = MPO(H, s)

#     # Do the time evolution by applying the gates
#     # for Nsteps steps
#     for β in 0:δτ:beta_max
#         energy = inner(rho, H)
#         println("β ", β)
#         println("energy ", energy)
#         rho = apply(gates, rho; cutoff)
#         rho = rho / tr(rho)
#     end

#     return rho
# end

# rho = get_rho(3)
# sample_list = []
# number_of_samples = 3
# for _ in 1:number_of_samples
#     sample_element = ITensors.sample(rho)
#     shifted_sample = join(string.(sample_element.-1))
#     println(shifted_sample)
#     append!(sample_list, shifted_sample)
# end

# sample_dictionary = counter(sample_list)
# print(sample_dictionary)
# print(Dict{String, Float64}(sample_dictionary))
# map!(x -> x./number_of_samples, values(sample_dictionary))
# b = bar(sample_dictionary)
# savefig(b, "barplot.png")

# N = 20
# r = [bitstring(UInt32(element))[(32-N+1):32] for element in 0:2^N-1]
# print(r)


# ---------------------------------------------------------------------------------

# gamma = 0.5
# N = 3

# i,j,a,b = Index(2, "i"), Index(2, "j"), Index(2, "a"), Index(2, "b") 
# i1, jN = Index(1, "i1"), Index(1, "iN")

# tm = zeros(Float64, 2, 2, 2, 2) 
# tm[1,1,:,:] = [1.0 0.0; 0.0 1.0]
# tm[2,2,:,:] = [1.0 0.0; 0.0 1.0]
# tm[2,1,:,:] = [1.0 0.0; 0.0 sqrt(1-gamma)]

# tmns = zeros(Float64, 2, 2, 2, 2) 
# tmns[1,1,:,:] = [1.0 0.0; 0.0 1.0]
# tmns[2,2,:,:] = [1.0 0.0; 0.0 1.0]
# tmns[2,1,:,:] = [1.0 0.0; 0.0 sqrt(1-gamma)]

# t1 = zeros(Float64, 2, 2, 2) 
# t1[2,:,:] = [1.0 0.0; 0.0 1.0]

# tf = zeros(Float64, 2, 2, 2) 
# tf[2,:,:] = [1.0 0.0; 0.0 1.0]


# mpo = MPO(N)
# j = 2
# link_indices = [Index(2, "Link,l=$idx") for idx in 1:(N-1)]

# for i in 1:N 

#     site_index = Index(2, "S=1/2,Site,n=$(i)")

#     if i == 1

#         mpo[i] = ITensor(t1, link_indices[i], site_index, site_index')
    
#     elseif i == j
        
#         mpo[i] = ITensor(tm, dag(link_indices[i-1]), link_indices[i], site_index, site_index')

#     elseif i == N

#         mpo[i] = ITensor(tf, dag(link_indices[i-1]), site_index, site_index')
    
#     else

#         mpo[i] = ITensor(tmns, dag(link_indices[i-1]), link_indices[i], site_index, site_index')
        
#     end

# end

# # @show mpo

# tmp = mpo[1]
# tmp = contract(tmp, mpo[2])
# tmp = contract(tmp, mpo[3])

# display(NDTensors.dense(tmp[1,1,:,:,1,1]))

# # include("MPO_iT.jl")
# # s = siteinds("S=1/2", 4)
# # # mps = randomMPS(Float64, s)
# # mpo = MPO(get_SW_local_z_OpSum(2), s)

# # @show mpo

# ---------------------------------------------------------------------------------

# println(get_two_point_correlator(2, 5, 5, 0.5))
# s = siteinds("S=1/2", 10)
# display(get_MPO_from_OpSum(get_two_point_correlator(2, 5, 5, 0.5), s))

# ---------------------------------------------------------------------------------

# using Printf
# using HDF5

# function gather_data()

#     gap_list = []
    
#     N_list = [8]
#     l_0_list = [0.125] # LinRange(0, 2*pi, 5)
#     x = 1.0
#     lambda = 100.0
#     ns = 200
#     val = 1/sqrt(pi)
#     mg_list = LinRange(-0.5, 0.0, 5)
#     D = 100

#     f = h5open("mass_shift_data.hdf5", "w")

#     for N in N_list

#         sites = siteinds("S=1/2", N)

#         for (idx_l_0, l_0) in enumerate(l_0_list)

#             # efd_list = []

#             for mg in mg_list

#                 charge_configuration = []
#                 electric_field_configuration = []

#                 mu = 2*mg*sqrt(x)

#                 H = get_Schwinger_staggered_Hamiltonian_OpSum(N, x, mu, l_0, lambda)
#                 mpo = get_MPO_from_OpSum(H, sites)

#                 initial_ansatz_0 = randomMPS(sites)
#                 sweeps = Sweeps(ns, maxdim = D) # , maxdim = D

#                 energy_0, psi_0 = dmrg(mpo, initial_ansatz_0, ishermitian = true, sweeps, outputlevel = 1) # , maxdim = D

#                 for i in 1:N-1
#                     link_mpo = get_MPO_from_OpSum(get_electric_field_link_operator(i, l_0), sites)
#                     electric_field = inner(psi_0', link_mpo, psi_0)
#                     push!(electric_field_configuration, electric_field)
#                 end
#                 for i in 1:N
#                     charge_mpo = get_MPO_from_OpSum(get_staggered_site_charge_operator(i), sites)
#                     charge = inner(psi_0', charge_mpo, psi_0)
#                     push!(charge_configuration, charge)
#                 end
#                 # scatter(1:N-1, electric_field_configuration, label = "EF")
#                 # display(scatter!(1:N, charge_configuration, label = "Q"))
#                 # push!(efd_list, real((electric_field_configuration[Int(floor(N/2))-1] + electric_field_configuration[Int(floor(N/2))])/2))

#                 Ms = [psi_0]
#                 w = abs(energy_0)
#                 initial_ansatz_1 = randomMPS(sites)

#                 energy_1, psi_1 = dmrg(mpo, Ms, initial_ansatz_1, sweeps, weight = w, ishermitian = true, outputlevel = 1) # , maxdim = D

#                 push!(gap_list, energy_1 - energy_0 - val)

#                 f["$(N)_$(l_0)_$(mg)_gap"] = energy_1 - energy_0 - val
#                 f["$(N)_$(l_0)_$(mg)_efd"] = real((electric_field_configuration[Int(floor(N/2))-1] + electric_field_configuration[Int(floor(N/2))])/2)
#                 # f["$(N)_$(l_0)_$(mg)_q"] = charge_configuration
#                 # f["$(N)_$(l_0)_$(mg)_ef"] = electric_field_configuration

#                 println(N, " ", l_0, " ", mg, " total charge = ", sum(charge_configuration))

#             end

#             # scatter(mg_list, efd_list, label = "")
#             # plot!(mg_list, efd_list, label = "EFD")

#             plot(mg_list, gap_list, label = "Gap")
#             hline!([0], label = "")
#             title!(@sprintf("l_0/pi = %.5f", l_0 / pi))
#             display(scatter!(mg_list, gap_list, label = ""))

#         end

#         # hline!([0], label = "")

#     end

#     close(f)

# end

# gather_data()

# ---------------------------------------------------------------------------------

# function get_efd_vs_l_0()

#     N = 12

#     vol = 10

#     x = (N/vol)^2

#     D = 40

#     mass_shift = 1.5 - 1.185401

#     mg = 1 - mass_shift

#     l_0_list = LinRange(0.8, 0.9, 5)

#     lambda = 100.0

#     acc = 1e-9

#     ns = 50

#     r = 1.0

#     D_p = 0

#     mg_p = 0.0

#     output_level = 0

#     initial_noise = 0.0

#     sites = siteinds("S=1/2", 2*N)

#     sweeps = Sweeps(ns, maxdim = D) # , maxdim = D

#     efd_list = []

#     params = Dict()
#     params["N"] = N 
#     params["x"] = x
#     params["mg"] = mg
#     params["r"] = r
#     params["lambda"] = lambda

#     for l_0 in l_0_list

#         println(l_0)

#         params["l_0"] = l_0

#         H = get_SW_OpSum(params)
#         mpo = get_MPO_from_OpSum(H, sites)

#         initial_ansatz_0 = randomMPS(sites)

#         energy_0, psi_0 = dmrg(mpo, initial_ansatz_0, ishermitian = true, sweeps, outputlevel = output_level) # , maxdim = D

#         z_config = get_z_configuration(psi_0, sites)
#         charge_config = get_SW_charge_configuration(z_config)
#         electric_field_config = get_SW_electric_field_configuration(charge_config, l_0)
#         push!(efd_list, real(electric_field_config[3]))

#     end

#     scatter(l_0_list, efd_list)

# end

# get_efd_vs_l_0()


# ---------------------------------------------------------------------------------

# Check the Wilson fermions Schwinger model Hamiltonian against my python computed

# N = 4
# x = 1.0
# D = 20
# l_0 = 1.0
# lambda = 1.0
# r = 1.0
# ns = 100
# mg = 1.0

# params = Dict("N" => N, "l_0" => l_0, "N" => N, "x" => x, "mg" => mg, "r" => r, "lambda" => lambda)

# sites = siteinds("S=1/2", 2*N)

# opsum = get_Schwinger_Wilson_OpSum(params)
# H = get_MPO_from_OpSum(opsum, sites)
# initial_ansatz_0 = randomMPS(sites, D)
# sweeps = Sweeps(ns, maxdim = D)
# energy, psi = dmrg(H, initial_ansatz_0, sweeps, ishermitian = true, maxdim = D)

# println(energy)

# ---------------------------------------------------------------------------------
