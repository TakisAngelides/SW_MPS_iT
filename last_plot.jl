using ITensors
using Plots
using LinearAlgebra
using Statistics
using CurveFit

function get_Staggered_Hamiltonian_MPO_paper(N, l_0, x, mg, lambda, sites)

    ampo = OpSum()

    for n=1:N-1
        ampo += x,"S+",n,"S-",n+1
        ampo += x,"S-",n,"S+",n+1
    end

    for n=0:N-1
        for k=n+1:N-1    
            ampo += (0.5*(N-k-1+lambda)),"Z",n+1,"Z",k+1
        end
    end

    for n=0:N-2
        ampo += (N/4 - 0.5*ceil(n/2) + l_0*(N-n-1)),"Z",n+1
    end

    for n=0:N-1
        ampo += (mg*sqrt(x)*(-1)^n),"Z",n+1
    end

    ampo += ((N-1)*l_0^2 + 0.5*l_0*N + (1/8)*N^2 + lambda*N/4),"Id",1

    return MPO(ampo, sites)

end

function get_efd(N, l_0, x, mg, lambda, nsweeps, Dmax, tol)

    sites = siteinds("S=1/2", N; conserve_qns = true)
    state = [isodd(n) ? "Dn" : "Up" for n=1:N]
    mps1 = randomMPS(sites, state)
    sweeps = Sweeps(nsweeps, maxdim = Dmax)

    H = get_Staggered_Hamiltonian_MPO_paper(N, l_0, x, mg, lambda, sites)

    sweeps = Sweeps(nsweeps, maxdim = Dmax)
    observer = DMRGObserver(;energy_tol = tol) # This allows to track the energy during DMRG and if the change in energy is smaller than tol it stops
    energy1, gs1 = dmrg(H, mps1, sweeps; ishermitian = true, outputlevel = 1, observer = observer)

    # gs = gs1
    state2 = [isodd(n) ? "Up" : "Dn" for n=1:N]
    mps2 = randomMPS(sites, state2)
    Ms = [mps1]
    w = abs(real(energy1))
    energy2, gs2 = dmrg(H, Ms, mps2, sweeps, weight = w, ishermitian = true, outputlevel = 1; observer = observer)
    if energy1 < energy2
        gs = gs1
    else
        println("SWITCHED")
        gs = gs2
    end
    
    z_config = expect(gs, "Z")
    efd = l_0 + 0.25 + 0.5*sum(z_config[1:div(N, 2)-1]) + 0.25*z_config[div(N, 2)]

    return efd

end

function get_pn(N, l_0, x, mg, lambda, nsweeps, Dmax, tol)

    sites = siteinds("S=1/2", N; conserve_qns = true)
    state = [isodd(n) ? "Dn" : "Up" for n=1:N]
    mps1 = randomMPS(sites, state)
    sweeps = Sweeps(nsweeps, maxdim = Dmax)

    H = get_Staggered_Hamiltonian_MPO_paper(N, l_0, x, mg, lambda, sites)

    sweeps = Sweeps(nsweeps, maxdim = Dmax)
    observer = DMRGObserver(;energy_tol = tol) # This allows to track the energy during DMRG and if the change in energy is smaller than tol it stops
    energy1, gs1 = dmrg(H, mps1, sweeps; ishermitian = true, outputlevel = 1, observer = observer)

    # gs = gs1
    state2 = [isodd(n) ? "Up" : "Dn" for n=1:N]
    mps2 = randomMPS(sites, state2)
    Ms = [mps1]
    w = abs(real(energy1))
    energy2, gs2 = dmrg(H, Ms, mps2, sweeps, weight = w, ishermitian = true, outputlevel = 1; observer = observer)
    if energy1 < energy2
        gs = gs1
    else
        println("SWITCHED")
        gs = gs2
    end
    
    stag = 0.5 .* [(-1)^n for n in 0:N-1]
    z_config = expect(gs, "Z")
    pn = 0.5*N + sum(stag .* z_config)

    return pn

end

function get_mass_shift(N, x, l_0, lambda, nsweeps, Dmax, tol, mg_star_tol, mg_left_initial, mg_right_initial)

    mg_values = []
    efd = []

    mg_left = mg_left_initial
    mg_right = mg_right_initial

    left = get_efd(N, l_0, x, mg_left, lambda, nsweeps, Dmax, tol)
    left_initial = left
    right = get_efd(N, l_0, x, mg_right, lambda, nsweeps, Dmax, tol)
    
    push!(mg_values, mg_left)
    push!(efd, left)
    push!(mg_values, mg_right)
    push!(efd, right)

    compare = abs(mg_right - mg_left)
    while compare > mg_star_tol

        println("The distance between the mgs is $(compare)")

        mg_middle = (mg_right + mg_left)/2
        middle = get_efd(N, l_0, x, mg_middle, lambda, nsweeps, Dmax, tol)

        push!(mg_values, mg_middle)
        push!(efd, middle)

        if left_initial < 0
            if middle < 0
                mg_left = mg_middle
                left = middle
            else
                mg_right = mg_middle
                right = middle
            end
        else
            if middle > 0
                mg_left = mg_middle
                left = middle
            else
                mg_right = mg_middle
                right = middle
            end
        end

        compare = abs(mg_right - mg_left)

    end

    return -(mg_left + mg_right)/2, mg_values, efd

end

function get_l0_star_pn(N, x, mg_r, lambda, nsweeps, Dmax, tol, l0_star_tol, l0_right_initial, l0_left_initial, mg_left_initial, mg_right_initial)

    l_0_values = []
    pn = []

    l0_left = l0_left_initial
    l0_right = l0_right_initial

    mass_shift_left, _, _ = get_mass_shift(N, x, l0_left, lambda, nsweeps, Dmax, tol, mg_star_tol, mg_left_initial, mg_right_initial)
    mass_shift_right, _, _ = get_mass_shift(N, x, l0_right, lambda, nsweeps, Dmax, tol, mg_star_tol, mg_left_initial, mg_right_initial)

    mg_lat_left = mg_r - mass_shift_left
    mg_lat_right = mg_r - mass_shift_right

    left = get_pn(N, l0_left, x, mg_lat_left, lambda, nsweeps, Dmax, tol)
    right = get_pn(N, l0_right, x, mg_lat_right, lambda, nsweeps, Dmax, tol)

    push!(l_0_values, l0_left)
    push!(pn, left)
    push!(l_0_values, l0_right)
    push!(pn, right)

    compare = abs(l0_right - l0_left)
    while compare > l0_star_tol

        println("The distance between l0s is $(compare)")

        l0_middle = (l0_right + l0_left)/2
        mass_shift_middle, _, _ = get_mass_shift(N, x, l0_middle, lambda, nsweeps, Dmax, tol, mg_star_tol, mg_left_initial, mg_right_initial)
        mg_lat = mg_r - mass_shift_middle
        middle = get_pn(N, l0_middle, x, mg_lat, lambda, nsweeps, Dmax, tol)

        push!(l_0_values, l0_middle)
        push!(pn, middle)

        if middle < 1
            l0_left = l0_middle
            left = middle
        else
            l0_right = l0_middle
            right = middle
        end

        compare = abs(l0_right - l0_left)

    end

    return (l0_left + l0_right)/2, l_0_values, pn

end

function get_l0_star_efd(N, x, mg_r, lambda, nsweeps, Dmax, tol, l0_star_tol, l0_right_initial, l0_left_initial, mg_left_initial, mg_right_initial)

    l_0_values = []
    efd = []

    l0_left = l0_left_initial
    l0_right = l0_right_initial

    mass_shift_left, _, _ = get_mass_shift(N, x, l0_left, lambda, nsweeps, Dmax, tol, mg_star_tol, mg_left_initial, mg_right_initial)
    mass_shift_right, _, _ = get_mass_shift(N, x, l0_right, lambda, nsweeps, Dmax, tol, mg_star_tol, mg_left_initial, mg_right_initial)

    mg_lat_left = mg_r - mass_shift_left
    mg_lat_right = mg_r - mass_shift_right

    left = get_efd(N, l0_left, x, mg_lat_left, lambda, nsweeps, Dmax, tol)
    left_initial = left
    right = get_efd(N, l0_right, x, mg_lat_right, lambda, nsweeps, Dmax, tol)

    push!(l_0_values, l0_left)
    push!(efd, left)
    push!(l_0_values, l0_right)
    push!(efd, right)

    compare = abs(l0_right - l0_left)
    while compare > l0_star_tol

        println("The distance between l0s is $(compare)")

        l0_middle = (l0_right + l0_left)/2
        mass_shift_middle, _, _ = get_mass_shift(N, x, l0_middle, lambda, nsweeps, Dmax, tol, mg_star_tol, mg_left_initial, mg_right_initial)
        mg_lat = mg_r - mass_shift_middle
        middle = get_efd(N, l0_middle, x, mg_lat, lambda, nsweeps, Dmax, tol)

        push!(l_0_values, l0_middle)
        push!(efd, middle)

        if left_initial < 0
            if middle < 0
                l0_left = l0_middle
                left = middle
            else
                l0_right = l0_middle
                right = middle
            end
        else
            if middle > 0
                l0_left = l0_middle
                left = middle
            else
                l0_right = l0_middle
                right = middle
            end
        end

        compare = abs(l0_right - l0_left)

    end

    return (l0_left + l0_right)/2, l_0_values, efd

end

function main()

    mg_lat = 10.0
    l_0 = 1.375

    lambda = 0
    
    nsweeps = 50
    Dmax = 100
    tol = 1e-9

    mg_star_tol = 1e-4
    mg_left_initial = -0.4
    mg_right_initial = 0.0

    l0_star_tol = 1e-4
    l0_left_initial = 0.4
    l0_right_initial = 1.1

    mg_r = 10.0

    N = 8
    vol = 30
    x = (N/vol)^2

    # Short test to see the phase transition with EFD for different values of N and volume
    p = plot()
    vol_list = [10, 20]
    for vol in vol_list

        x = 10
        N = round(Int64, sqrt(x)*vol)
        efd_list = []
        l_0_list = LinRange(0, 1.5, 5)
        for l_0 in l_0_list
            mass_shift, _, _ = get_mass_shift(N, x, l_0, lambda, nsweeps, Dmax, tol, mg_star_tol, mg_left_initial, mg_right_initial)
            mg_lat = mg_r - mass_shift
            push!(efd_list, get_efd(N, l_0, x, mg_lat, lambda, nsweeps, Dmax, tol))
        end
        plot!(p, l_0_list, efd_list, label = "N = $(N), x = $(x), vol = $(vol)")
    
    end
    savefig(p, "efd_final_plot.png")

    # Short test to see the phase transition with EFD for different values of N and x
    # p = plot()
    # N_list = [8, 10, 12]
    # vol = 30
    # for N in N_list

    #     x = (N/vol)^2
    #     efd_list = []
    #     l_0_list = LinRange(0, 2.5, 5)
    #     for l_0 in l_0_list
    #         mass_shift, _, _ = get_mass_shift(N, x, l_0, lambda, nsweeps, Dmax, tol, mg_star_tol, mg_left_initial, mg_right_initial)
    #         mg_lat = mg_r - mass_shift
    #         push!(efd_list, get_efd(N, l_0, x, mg_lat, lambda, nsweeps, Dmax, tol))
    #     end
    #     plot!(p, l_0_list, efd_list, label = "N = $(N), x = $(x), vol = $(vol)")
    
    # end
    # savefig(p, "efd_final_plot.png")

    # Short test to see that the mass shift is working
    # mass_shift, mg_values, pn = get_mass_shift(N, x, l_0, lambda, nsweeps, Dmax, tol, mg_star_tol, mg_left_initial, mg_right_initial)
    # perm = sortperm(mg_values)
    # mg_values = mg_values[perm]
    # pn = pn[perm]
    # p = plot()
    # plot!(p, mg_values, pn)
    # scatter!(p, mg_values, pn)
    # vline!([-mass_shift])
    # savefig(p, "mass_shift_final_plot.png")
    # println(mass_shift)

    # Short test to see the phase transition with PN
    # let
    # pn_list = []
    # l_0_list = LinRange(0.1, 1.1, 4)
    # for l_0 in l_0_list
    #     mass_shift, _, _ = get_mass_shift(N, x, l_0, lambda, nsweeps, Dmax, tol, mg_star_tol, mg_left_initial, mg_right_initial)
    #     mg_lat = mg_r - mass_shift
    #     push!(pn_list, get_pn(N, l_0, x, mg_lat, lambda, nsweeps, Dmax, tol))
    # end
    # p = plot()
    # plot!(p, l_0_list, pn_list)
    # savefig(p, "pt_final_plot.png")
    # end

    # Short test to see the phase transition with PN
    # p = plot()
    # efd_list = []
    # l_0_list = LinRange(0, 2.5, 5)
    # for l_0 in l_0_list
    #     mass_shift, _, _ = get_mass_shift(N, x, l_0, lambda, nsweeps, Dmax, tol, mg_star_tol, mg_left_initial, mg_right_initial)
    #     mg_lat = mg_r - mass_shift
    #     push!(efd_list, get_efd(N, l_0, x, mg_lat, lambda, nsweeps, Dmax, tol))
    # end
    # plot!(p, l_0_list, efd_list, label = "N = $(N), x = $(x), vol = $(vol)")
    # savefig(p, "efd_final_plot_N_$(N)_x_$(x)_vol_$(vol).png")

    # Short test to see that the getting the l_0_star is working with PN
    # l_0_star_pn, l_0_values_pn, pn = get_l0_star_pn(N, x, mg_r, lambda, nsweeps, Dmax, tol, l0_star_tol, l0_right_initial, l0_left_initial, mg_left_initial, mg_right_initial)
    # perm = sortperm(l_0_values_pn)
    # l_0_values_pn = l_0_values_pn[perm]
    # pn = pn[perm]
    # p = plot()
    # plot!(p, l_0_values_pn, pn)
    # scatter!(p, l_0_values_pn, pn)
    # vline!([l_0_star_pn])
    # savefig(p, "l_0_star_final_plot_pn.png")
    # println(l_0_star_pn)

    # Short test to see that the getting the l_0_star is working with EFD
    # l_0_star_efd, l_0_values_efd, efd = get_l0_star_efd(N, x, mg_r, lambda, nsweeps, Dmax, tol, l0_star_tol, l0_right_initial, l0_left_initial, mg_left_initial, mg_right_initial)
    # perm = sortperm(l_0_values_efd)
    # l_0_values_efd = l_0_values_efd[perm]
    # efd = efd[perm]
    # p = plot()
    # plot!(p, l_0_values_efd, efd)
    # scatter!(p, l_0_values_efd, efd)
    # vline!([l_0_star_efd])
    # savefig(p, "l_0_star_final_plot_efd.png")
    # println(l_0_star_efd)

    # Extrapolate to zero lattice spacing at fixed volume
    # N = [8, 10, 20]
    # ag = vol ./ N
    # l_0_star = [0.675732421875, 0.6639160156250001, 0.6476074218750001]
    # res = linear_fit(ag, l_0_star)
    # println(res)

end

main()
