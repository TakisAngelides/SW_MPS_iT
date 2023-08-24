include("run_DMRG_iT.jl")

w_1_s_2 = parse(Int64, ARGS[1])
N = parse(Int, ARGS[2]) # Number of physical lattice sites
x = parse(Float64, ARGS[3]) # 1/(ag)^2
D = parse(Int64, ARGS[4]) # Bond dimension of MPS to be computed
mg = parse(Float64, ARGS[5]) # m/g
l_0 = parse(Float64, ARGS[6]) # l_0 = theta/(2*pi)
lambda = parse(Float64, ARGS[7]) # Lagrange multiplier to enforce total charge squared to be 0
acc = parse(Float64, ARGS[8]) # Tolerance for stopping condition of the variational algorithm
ns = parse(Int64, ARGS[9]) # Maximum number of sweeps of the variational algorithm
D_p = parse(Int64, ARGS[10]) # Bond dimension of ansatz
mg_p = parse(Float64, ARGS[11]) # m/g of ansatz
silent = parse(Bool, ARGS[12]) # Whether to compute total charge, energy, average electric field density per sweep
initial_noise = parse(Float64, ARGS[13]) # Inital noise that decays with each sweep for the DMRG
first_excited = parse(Bool, ARGS[14]) # If this is true we will compute the first excited state only and not the ground state
r = 1.0

path_to_tangelides = "/lustre/fs24/group/cqta/tangelides/"

if !first_excited

    sweep_observables_file_path = "$(path_to_tangelides)SW_Sweep_Observables_iT/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_w1s2_$(w_1_s_2)_fe_$(first_excited).txt"

    mps_file_path = "$(path_to_tangelides)SW_MPS_States_iT/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_w1s2_$(w_1_s_2)_fe_$(first_excited).hdf5"

    previous_mps_file_path = "$(path_to_tangelides)SW_MPS_States_iT/N_$(N)_x_$(x)_D_$(D_p)_l0_$(l_0)_mg_$(mg_p)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_w1s2_$(w_1_s_2)_fe_$(first_excited).hdf5"

    if isfile(previous_mps_file_path)
        f = h5open(previous_mps_file_path, "r")
        previous_psi = read(f, "MPS", MPS)
        sites = siteinds(previous_psi)
        close(f)
    else
        if w_1_s_2 == 1
            sites = siteinds("S=1/2", 2*N; conserve_qns = true)
        else
            sites = siteinds("S=1/2", N; conserve_qns = true)
        end
        state = [isodd(n) ? "0" : "1" for n = 1:length(sites)]
        previous_psi = randomMPS(sites, state; linkdims = D)
    end

    if w_1_s_2 == 1
        params = Dict("first_excited" => first_excited, "initial_noise" => initial_noise, "silent" => silent, "N" => N, "D" => D, "x" => x, "ns" => ns, "lambda" => lambda, "l_0" => l_0, "mg" => mg, "r" => r, "acc" => acc, "sweep_observables_file_path" => sweep_observables_file_path, "previous_mps_file_path" => previous_mps_file_path, "previous_psi" => previous_psi)
        H = get_MPO_from_OpSum(get_Schwinger_Wilson_OpSum(params), sites)
    else 
        mu = 2*mg*sqrt(x)
        params = Dict("first_excited" => first_excited, "initial_noise" => initial_noise, "silent" => silent, "N" => N, "D" => D, "x" => x, "ns" => ns, "lambda" => lambda, "l_0" => l_0, "mu" => mu, "acc" => acc, "sweep_observables_file_path" => sweep_observables_file_path, "previous_mps_file_path" => previous_mps_file_path, "previous_psi" => previous_psi)
        H = get_MPO_from_OpSum(get_Schwinger_staggered_Hamiltonian_OpSum(params), sites)
    end

    energy, psi = run_SW_DMRG(sites, params, H, true)

    # Save the MPS as h5 file including its sites object

    if isfile(mps_file_path)
        rm(mps_file_path)
    end

    f = h5open(mps_file_path, "w")
    write(f, "MPS", psi)
    write(f, "gs_energy", energy)
    close(f)

else

    sweep_observables_file_path = "$(path_to_tangelides)SW_Sweep_Observables_First_Excited_iT/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_w1s2_$(w_1_s_2)_fe_$(first_excited).txt"

    previous_mps_file_path = "$(path_to_tangelides)SW_MPS_First_Excited_States_iT/N_$(N)_x_$(x)_D_$(D_p)_l0_$(l_0)_mg_$(mg_p)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_w1s2_$(w_1_s_2)_fe_$(first_excited).hdf5"

    mps_file_path = "$(path_to_tangelides)SW_MPS_First_Excited_States_iT/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_w1s2_$(w_1_s_2)_fe_$(first_excited).hdf5"

    psi_0_file_path = "$(path_to_tangelides)SW_MPS_States_iT/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_w1s2_$(w_1_s_2)_fe_false.hdf5"
    
    f = h5open(psi_0_file_path, "r")
    psi_0 = read(f, "MPS", MPS)
    energy_0 = read(f, "gs_energy")
    sites = siteinds(psi_0)
    close(f)

    if isfile(previous_mps_file_path)
        f = h5open(previous_mps_file_path, "r")
        previous_psi = read(f, "first_excited_MPS", MPS)
        close(f)
    else
        state = [isodd(n) ? "1" : "0" for n = 1:length(sites)]
        previous_psi = randomMPS(sites, state; linkdims = D)
    end

    if w_1_s_2 == 1
        params = Dict("Ms" => [psi_0], "w" => abs(energy_0), "first_excited" => first_excited, "initial_noise" => initial_noise, "silent" => silent, "N" => N, "D" => D, "x" => x, "ns" => ns, "lambda" => lambda, "l_0" => l_0, "mg" => mg, "r" => r, "acc" => acc, "sweep_observables_file_path" => sweep_observables_file_path, "previous_mps_file_path" => previous_mps_file_path, "previous_psi" => previous_psi)
        H = get_MPO_from_OpSum(get_Schwinger_Wilson_OpSum(params), sites)
    else 
        mu = 2*mg*sqrt(x)
        params = Dict("Ms" => [psi_0], "w" => abs(energy_0), "first_excited" => first_excited, "initial_noise" => initial_noise, "silent" => silent, "N" => N, "D" => D, "x" => x, "ns" => ns, "lambda" => lambda, "l_0" => l_0, "mu" => mu, "acc" => acc, "sweep_observables_file_path" => sweep_observables_file_path, "previous_mps_file_path" => previous_mps_file_path, "previous_psi" => previous_psi)
        H = get_MPO_from_OpSum(get_Schwinger_staggered_Hamiltonian_OpSum(params), sites)
    end

    energy, psi = run_SW_DMRG(sites, params, H, true)

    # Save the MPS as h5 file including its sites object

    if isfile(mps_file_path)
        rm(mps_file_path)
    end

    f = h5open(mps_file_path, "w")
    write(f, "first_excited_MPS", psi)
    write(f, "first_excited_energy", energy)
    close(f)

end