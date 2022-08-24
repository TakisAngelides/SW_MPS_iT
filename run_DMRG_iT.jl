using KrylovKit
using ITensors.HDF5
using Plots
include("DMRG_iT.jl")
include("Observables_iT.jl")

function run_SW_DMRG(sites, params, H)

    energy, psi = DMRG(H, sites, params)

    return energy, psi

end

function run_Ising_DMRG()

    # Exact gs energy for N => 10, J => -1, g_z => -0.1, g_x => 1.5, ns => 5, D => 34 is: -16.706153497716304

    params = Dict("N" => 10, "J" => -1, "g_z" => -0.1, "g_x" => 1.5, "ns" => 5, "D" => 34)

    opsum_ising = get_Ising_OpSum(N, J, g_z, g_x)
    sites = siteinds("S=1/2", N)
    H_1 = get_MPO_from_OpSum(opsum_ising, sites)
    energy_1, _ = DMRG(H_1, sites, params)
    println(energy_1)

end

function main()

    # N = 4
    # D = 10
    # x = 1.0
    # ns = 100
    # lambda = 100.0
    # l_0 = 0.125
    # mg = 0.125
    # r = 1.0
    # acc = 1E-10

    # file_path = "sweep_observables_N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r).txt"
    # params = Dict("N" => N, "D" => D, "x" => x, "ns" => ns, "lambda" => lambda, "l_0" => l_0, "mg" => mg, "r" => r, "acc" => acc, "sweep_observables_file_path" => file_path)
    # sites = siteinds("S=1/2", 2*N)

    # energy, psi = run_SW_DMRG(sites, params)

    # sweep_list = []
    # total_charge = []
    # energy_list = []
    # avg_E_field = []

    # open(file_path, "r") do f
        
    #     for line in eachline(f)

    #         row = split(line, ",")
    #         append!(sweep_list, real(parse(Float64, row[1])))
    #         append!(total_charge, real(parse(Float64, split(row[2], " ")[1])))
    #         append!(energy_list, real(parse(Float64, split(row[3], " ")[1])))
    #         append!(avg_E_field, real(parse(Float64, split(row[4], " ")[1])))

    #     end
    
    # end

    # scatter(sweep_list, total_charge, label = "Total Charge")
    # plot!(sweep_list, energy_list, label = "Energy")
    # scatter!(sweep_list, avg_E_field, label = "Average Electric Field")

    # f = h5open("my_psi.h5", "w")
    # write(f, "my_psi", psi)
    # close(f)

    # f = h5open("my_psi.h5","r")
    # my_psi = read(f, "my_psi", MPS)
    # close(f)

end

# main()
