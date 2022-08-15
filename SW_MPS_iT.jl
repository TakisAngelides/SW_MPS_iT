include("run_DMRG_iT.jl")

N = parse(Int, ARGS[1]) # Number of physical lattice sites
x = parse(Float64, ARGS[2]) # 1/(ag)^2
mg = parse(Float64, ARGS[3]) # m/g
D = parse(Int64, ARGS[4]) # Bond dimension of MPS to be computed
l_0 = parse(Float64, ARGS[5]) # l_0 = theta/(2*pi)
lambda = parse(Float64, ARGS[6]) # Lagrange multiplier to enforce total charge squared to be 0
acc = parse(Float64, ARGS[7]) # Tolerance for stopping condition of the variational algorithm
ns = parse(Int64, ARGS[8]) # Maximum number of sweeps of the variational algorithm
D_p = parse(Int64, ARGS[9]) # Bond dimension of ansatz
mg_p = parse(Float64, ARGS[10]) # m/g of ansatz
r = parse(Float64, ARGS[11]) # Wilson parameter
silent = parse(Bool, ARGS[12]) # Whether to compute total charge, energy, average electric field density per sweep
initial_noise = parse(Float64, ARGS[13]) # Inital noise that decays with each sweep for the DMRG

sweep_observables_file_path = "/lustre/fs23/group/nic/tangelides/SW_Sweep_Observables_iT/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r).txt"

mps_file_path = "/lustre/fs23/group/nic/tangelides/SW_MPS_iT/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r).h5"

previous_mps_file_path = "/lustre/fs23/group/nic/tangelides/SW_MPS_iT/N_$(N)_x_$(x)_D_$(D_p)_l0_$(l_0)_mg_$(mg_p)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r).h5"

params = Dict("initial_noise" => initial_noise, "silent" => silent, "N" => N, "D" => D, "x" => x, "ns" => ns, "lambda" => lambda, "l_0" => l_0, "mg" => mg, "r" => r, "acc" => acc, "sweep_observables_file_path" => sweep_observables_file_path, "previous_mps_file_path" => previous_mps_file_path)

sites = siteinds("S=1/2", 2*N)

# Compute the MPS

_, psi = run_SW_DMRG(sites, params)

# Save the MPS as h5 file

f = h5open(mps_file_path, "w")
write(f, "MPS", psi)
close(f)
