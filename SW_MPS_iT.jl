include("run_DMRG_iT.jl")

N = parse(Int, ARGS[1]) # Number of physical lattice sites
x = parse(Float64, ARGS[2]) # 1/(ag)^2
D = parse(Int64, ARGS[3]) # Bond dimension of MPS to be computed
mg = parse(Float64, ARGS[4]) # m/g
l_0 = parse(Float64, ARGS[5]) # l_0 = theta/(2*pi)
lambda = parse(Float64, ARGS[6]) # Lagrange multiplier to enforce total charge squared to be 0
acc = parse(Float64, ARGS[7]) # Tolerance for stopping condition of the variational algorithm
ns = parse(Int64, ARGS[8]) # Maximum number of sweeps of the variational algorithm
D_p = parse(Int64, ARGS[9]) # Bond dimension of ansatz
mg_p = parse(Float64, ARGS[10]) # m/g of ansatz
r = parse(Float64, ARGS[11]) # Wilson parameter
silent = parse(Bool, ARGS[12]) # Whether to compute total charge, energy, average electric field density per sweep
initial_noise = parse(Float64, ARGS[13]) # Inital noise that decays with each sweep for the DMRG
w_1_s_2 = parse(Int64, ARGS[14]) # If this is equal to 1 do Wilson, else do staggered
first_excited = parse(Bool, ARGS[15]) # If this is true we will compute the first excited state only and not the ground state

if !first_excited

    sweep_observables_file_path = "/lustre/fs23/group/nic/tangelides/SW_Sweep_Observables_iT/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r)_w1s2_$(w_1_s_2)_fe_$(first_excited).txt"

    mps_file_path = "/lustre/fs23/group/nic/tangelides/SW_MPS_States_iT/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r)_w1s2_$(w_1_s_2)_fe_$(first_excited).h5"

    previous_mps_file_path = "/lustre/fs23/group/nic/tangelides/SW_MPS_States_iT/N_$(N)_x_$(x)_D_$(D_p)_l0_$(l_0)_mg_$(mg_p)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r)_w1s2_$(w_1_s_2)_fe_$(first_excited).h5"

    if isfile(previous_mps_file_path)
        f = h5open(previous_mps_file_path, "r")
        previous_psi = read(f, "MPS", MPS)
        sites = siteinds(previous_psi)
        close(f)
    else
        if w_1_s_2 == 1
            sites = siteinds("S=1/2", 2*N)
        else
            sites = siteinds("S=1/2", N)
        end
        previous_psi = randomMPS(sites, D)
    end

    params = Dict("initial_noise" => initial_noise, "silent" => silent, "N" => N, "D" => D, "x" => x, "ns" => ns, "lambda" => lambda, "l_0" => l_0, "mg" => mg, "r" => r, "acc" => acc, "sweep_observables_file_path" => sweep_observables_file_path, "previous_mps_file_path" => previous_mps_file_path, "previous_psi" => previous_psi)

    # Compute the MPS

    if w_1_s_2 == 1
        H = get_MPO_from_OpSum(get_SW_OpSum(params), sites)
    else
        H = get_MPO_from_OpSum(get_SW_Staggered_OpSum(params), sites)
    end

    energy, psi = run_SW_DMRG(sites, params, H)

    # Save the MPS as h5 file including its sites object

    if isfile(mps_file_path)
        rm(mps_file_path)
    end

    f = h5open(mps_file_path, "w")
    write(f, "MPS", psi)
    write(f, "gs_energy", energy)
    close(f)

else

    sweep_observables_file_path = "/lustre/fs23/group/nic/tangelides/SW_Sweep_Observables_First_Excited_iT/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r)_w1s2_$(w_1_s_2)_fe_$(first_excited).txt"

    previous_mps_file_path = "/lustre/fs23/group/nic/tangelides/SW_MPS_First_Excited_States_iT/N_$(N)_x_$(x)_D_$(D_p)_l0_$(l_0)_mg_$(mg_p)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r)_w1s2_$(w_1_s_2)_fe_$(first_excited).h5"

    mps_file_path = "/lustre/fs23/group/nic/tangelides/SW_MPS_First_Excited_States_iT/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r)_w1s2_$(w_1_s_2)_fe_$(first_excited).h5"

    psi_0_file_path = "/lustre/fs23/group/nic/tangelides/SW_MPS_States_iT/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ns_$(ns)_acc_$(acc)_lam_$(lambda)_r_$(r)_w1s2_$(w_1_s_2)_fe_false.h5"
    f = h5open(psi_0_file_path, "r")
    psi_0 = read(f, "MPS", MPS)
    energy_0 = read(f, "gs_energy")
    sites = siteinds(psi_0)
    close(f)

    if isfile(previous_mps_file_path)
        f = h5open(previous_mps_file_path, "r")
        previous_psi = read(f, "MPS", MPS)
        close(f)
    else
        previous_psi = randomMPS(sites, D)
    end

    params = Dict("initial_noise" => initial_noise, "silent" => silent, "N" => N, "D" => D, "x" => x, "ns" => ns, "lambda" => lambda, "l_0" => l_0, "mg" => mg, "r" => r, "acc" => acc, "sweep_observables_file_path" => sweep_observables_file_path, "previous_mps_file_path" => previous_mps_file_path, "previous_psi" => previous_psi)

    # Compute the MPS

    if w_1_s_2 == 1
        P = outer(psi_0', psi_0)
        H = get_MPO_from_OpSum(get_SW_OpSum(params), sites)
        Heff = H + energy_0.*P
    else
        P = outer(psi_0', psi_0)
        H = get_MPO_from_OpSum(get_Staggered_OpSum(params), sites)
        Heff = H + energy_0.*P
    end

    energy, psi = run_SW_DMRG(sites, params, Heff)

    # Save the MPS as h5 file including its sites object

    if isfile(mps_file_path)
        rm(mps_file_path)
    end

    f = h5open(mps_file_path, "w")
    write(f, "first_excited_MPS", psi)
    write(f, "first_excited_energy", energy)
    close(f)

end