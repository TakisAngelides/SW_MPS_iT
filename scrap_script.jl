
using ITensors
using Plots
using DataStructures
using StatsBase
using Dates

println("Program started running and the time is: ", now())

N = 2

function ITensors.op(::OpName"expτSS", ::SiteType"S=1/2", s1::Index, s2::Index; tau)
    h = -4 * op("Sz", s1) * op("Sz", s2)
    return exp(tau * h)
  end
  
function get_rho(N, dt, beta; cutoff = 1E-8)

    # Make an array of 'site' indices
    s = siteinds("S=1/2", N)

    t = time()
    # Make gates (1,2),(2,3),(3,4),...
    gates = ops([("expτSS", (n, m), (tau = -dt / 2,)) for n in 1:N for m in n+1:N], s)
    dt = time() - t
    println("Writing gates as ITensor objects took: ", dt, " seconds")

    # Include gates in reverse order to complete Trotter formula
    append!(gates, reverse(gates))

    # Initial state is infinite-temperature mixed state
    rho = MPO(s, "Id") ./ √2

    # Make H for measuring the energy        
    function get_Hamiltonian(N)
        H = OpSum()
        for n=1:N-1
            H += -4,"Sz",n,"Sz",n+1
        end
        return H
    end
    H = get_Hamiltonian(N)
    t = time()
    H = MPO(H, s)
    dt = time() - t
    println("Getting the MPO of the Hamiltonian took: ", dt, " seconds")

    # Do the time evolution by applying the gates
    # for Nsteps steps
    for beta in 0:dt:beta
        t = time()
        energy = inner(rho, H)
        println("------------------------------------------------------------------------")
        println("beta = ", beta)
        println("energy = ", energy)
        rho = apply(gates, rho; cutoff)
        rho = rho / tr(rho)
        dt = time() - t
        println("This evolution iteration took ", dt, " seconds")
        println("The maximum bond dimension of the density matrix MPO is: ", maxlinkdim(rho))    
    end
    println("------------------------------------------------------------------------")
    return rho
end

beta = 0.02
dt = 0.01
rho = get_rho(N, dt, beta)

# bitstrings = [bitstring(UInt64(element)) for element in 0:2^N-1]

# sample_dictionary = Dict{String, Float64}(element => 0 for element in bitstrings)

number_of_samples = 100

open("samples_N_$(N)_beta_$(beta)_dt_$(dt).txt", "w") do f

    for _ in 1:number_of_samples
        
        sample_element = ITensors.sample(rho)
        
        shifted_sample = join(string.(sample_element.-1))

        write(f, "$(shifted_sample)\n")

        # println("Sample: ", shifted_sample)
        
        # sample_dictionary[shifted_sample] += 1
    end

end

# map!(x -> x./number_of_samples, values(sample_dictionary))

# b = bar(sample_dictionary);
# savefig(b, "barplot.png")

