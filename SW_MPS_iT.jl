include("run_DMRG_iT.jl")

N = parse(Int, ARGS[1]) # Number of physical lattice sites
x = parse(Float64, ARGS[2]) # 1/(ag)^2
mg = parse(Float64, ARGS[3]) # m/g
D = parse(Int64, ARGS[4]) # Bond dimension of MPS to be computed
l_0 = parse(Float64, ARGS[5]) # l_0 = theta/(2*pi)
lambda = parse(Float64, ARGS[6]) # Lagrange multiplier to enforce total charge squared to be 0
acc = parse(Float64, ARGS[7]) # Tolerance for stopping condition of the variational algorithm
ms = parse(Int64, ARGS[8]) # Maximum number of sweeps of the variational algorithm
D_p = parse(Int64, ARGS[9]) # Bond dimension of ansatz
mg_p = parse(Float64, ARGS[10]) # m/g of ansatz
r = parse(Float64, ARGS[11]) # Wilson parameter

file_path = "sweep_observables_N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r).txt"
params = Dict("N" => N, "D" => D, "x" => x, "ns" => ns, "lambda" => lambda, "l_0" => l_0, "mg" => mg, "r" => r, "acc" => acc, "sweep_observables_file_path" => file_path)
sites = siteinds("S=1/2", 2*N)

energy, psi = run_SW_DMRG(sites, params)

f = h5open("my_psi.h5", "w")
write(f, "my_psi", psi)
close(f)
