function get_SW_charge_configuration(z_configuration_list::Vector{ComplexF64})::Vector{ComplexF64}

    N_spin::Int64 = length(z_configuration_list)
    N::Int64 = Int(N_spin/2)

    charge_configuration_list::Vector{ComplexF64} = []
    
    for k=1:N
    
        charge_on_site_k::ComplexF64 = 0.5*(z_configuration_list[2*k-1] + z_configuration_list[2*k])
    
        append!(charge_configuration_list, charge_on_site_k)
    end

    return charge_configuration_list

end

function get_z_configuration(psi::MPS, sites::Vector{Index{Int64}})::Vector{ComplexF64}

    N_spin::Int64 = length(psi)

    z_configuration_list::Vector{ComplexF64} = []
    
    for i=1:N_spin
    
        z_OpSum::Sum{Scaled{ComplexF64, Prod{Op}}} = get_SW_local_z_OpSum(i)
    
        z_mpo::MPO = get_MPO_from_OpSum(z_OpSum, sites)
    
        append!(z_configuration_list, inner(psi', z_mpo, psi))
    end

    return z_configuration_list

end

function get_SW_electric_field_configuration(charge_configuration_list::Vector{ComplexF64}, l_0::Float64)::Vector{ComplexF64}

    N::Int64 = length(charge_configuration_list)

    electric_field_configuration_list::Vector{ComplexF64} = []

    for i in 1:N-1
        
        L_n::ComplexF64 = l_0 + sum(charge_configuration_list[1:i])
        
        append!(electric_field_configuration_list, L_n)
    end

    return electric_field_configuration_list

end

function get_SW_local_charge_OpSum(site_idx::Int64)::Sum{Scaled{ComplexF64, Prod{Op}}}

    ampo::Sum{Scaled{ComplexF64, Prod{Op}}} = OpSum()

    ampo += "Sz",2*site_idx-1
    ampo += "Sz",2*site_idx

    return ampo

end

function get_SW_local_z_OpSum(site_idx::Int64)::Sum{Scaled{ComplexF64, Prod{Op}}}

    ampo::Sum{Scaled{ComplexF64, Prod{Op}}} = OpSum()

    ampo += 2,"Sz",site_idx

    return ampo

end

function get_SW_OpSum(params)::Sum{Scaled{ComplexF64, Prod{Op}}}

    """
    N = number of physical sites
    """

    N::Int64 = params["N"]
    l_0::Float64 = params["l_0"]
    x::Float64 = params["x"]
    mg::Float64 = params["mg"]
    r::Float64 = params["r"]
    lambda::Float64 = params["lambda"]

    A = 1im*x*(r - 1)
    B = 1im*x*(r + 1)
    C = 2*1im*(sqrt(x)*mg + x*r)
    F = (N-1)*l_0^2 + N*(N-1)/4 + lambda*N/2

    ampo::Sum{Scaled{ComplexF64, Prod{Op}}} = OpSum()

    for n=1:N-1
        ampo += A,"S-",2*n,"S+",2*n+1
        ampo += -A,"S+",2*n,"S-",2*n+1
        ampo += 4*B,"S+",2*n-1,"Sz",2*n,"Sz",2*n+1,"S-",2*n+2
        ampo += -4*B,"S-",2*n-1,"Sz",2*n,"Sz",2*n+1,"S+",2*n+2
    end

    for n=1:N
        ampo += C,"S-",2*n-1,"S+",2*n
        ampo += -C,"S+",2*n-1,"S-",2*n
    end

    for k=1:2*N-2
        D = l_0*(N-ceil(k/2))
        ampo += 2*D,"Sz",k
    end

    for k=1:2*N
        for p=k+1:2*N
            E = 0.5*(N-ceil(p/2)+lambda)
            ampo += 4*E,"Sz",k,"Sz",p
        end
    end

    ampo += F,"Id",2

    return ampo

end

function get_SW_OpSum_alternative(params)::Sum{Scaled{ComplexF64, Prod{Op}}}

    """
    N = number of physical sites
    """

    N::Int64 = params["N"]
    l_0::Float64 = params["l_0"]
    x::Float64 = params["x"]
    mg::Float64 = params["mg"]
    r::Float64 = params["r"]
    lambda::Float64 = params["lambda"]

    A = 1im*x*(r - 1)
    B = 1im*x*(r + 1)
    C = 2*1im*(sqrt(x)*mg + x*r)
    F = (N-1)*l_0^2 + N*(N-1)/4 + lambda*N/2

    ampo::Sum{Scaled{ComplexF64, Prod{Op}}} = OpSum()

    for n=1:N-1
        ampo += A,"S-",2*n,"S+",2*n+1
        ampo += A,"S+",2*n,"S-",2*n+1
        ampo += 4*B,"S+",2*n-1,"Sz",2*n,"Sz",2*n+1,"S-",2*n+2
        ampo += -4*B,"S-",2*n-1,"Sz",2*n,"Sz",2*n+1,"S+",2*n+2
    end

    for n=1:N
        ampo += C,"S-",2*n-1,"S+",2*n
        ampo += -C,"S+",2*n-1,"S-",2*n
    end

    for k=1:N-1
        J = l_0*(N-k)
        ampo += 2*J,"Sz",2*k
        ampo += 2*J,"Sz",2*k-1
    end

    for k=1:N-1
        for p=k+1:N
            E = 0.5*(N-p+lambda)
            ampo += 4*E,"Sz",2*k-1,"Sz",2*p-1
            ampo += 4*E,"Sz",2*k-1,"Sz",2*p
            ampo += 4*E,"Sz",2*k,"Sz",2*p-1
            ampo += 4*E,"Sz",2*k,"Sz",2*p
        end
    end

    for k=1:N
        G = 0.5*(N - k + lambda)
        ampo += 4*G,"Sz",2*k-1,"Sz",2*k
    end

    ampo += F,"Id",1

    return ampo

end

function get_Ising_OpSum(N::Int64, J::Float64, g_z::Float64, g_x::Float64)::Sum{Scaled{ComplexF64, Prod{Op}}}
        
    ampo = OpSum()

    for n=1:N-1
        ampo += 4*J,"Sz",n,"Sz",n+1
    end

    for n=1:N
        ampo += 2*g_x,"Sx",n
        ampo += 2*g_z,"Sz",n
    end

    return ampo

end

function get_MPO_from_OpSum(OpSum, sites::Vector{Index{Int64}})::MPO

    return MPO(OpSum, sites)

end

function get_two_point_correlator(n, m, N, x)

    """
    Inputs:

    n = left index of correlator
    m = right index of correlator
    N = number of lattice sites (=> we have 2*N spin sites)
    x = 1/(ag)^2

    Outputs:

    The OpSum corresponding to the operator of the two point correlation function psi_bar_n*psi_m/g

    """

    @assert(n <= m, "The condition n <= m was not satisfied, the left index of the correlator needs to be less than the right index")
    @assert(m <= N, "The condition m <= N was not satisfied, the right index of the correlator needs to be within the boundaries of the lattice")

    os = OpSum()

    C = ((-1)^(n+m))*sqrt(x)

    term_1 = []
    term_2 = []

    push!(term_1, C*2^((2*m-1)-(2*n-1)+1)) # "Sz" is Z/2 so we multiply a factor here to make it Z
    push!(term_2, C*2^((2*m-2)-(2*n)+1)) 
        
    for i in 1:2*N
        if i < (2*n-1)
            push!(term_1, "I", i)
        elseif i == 2*n-1
            push!(term_1, "S+", i)
        elseif i < 2*m
            push!(term_1, "Sz", i)
        elseif i == 2*m
            push!(term_1, "S-", i)
        else
            push!(term_1, "I", i)
        end
    end

    for i in 1:2*N
        if i < (2*n)
            push!(term_2, "I", i)
        elseif i == 2*n
            push!(term_2, "S+", i)
        elseif i < (2*m-1)
            push!(term_2, "Sz", i)
        elseif i == 2*m-1
            push!(term_2, "S-", i)
        else
            push!(term_2, "I", i)
        end
    end

    os += tuple(term_1...)
    os += tuple(term_2...)

    return os

end

function get_Staggered_OpSum(params)::Sum{Scaled{ComplexF64, Prod{Op}}}

    N::Int64 = params["N"]
    l_0::Float64 = params["l_0"]
    x::Float64 = params["x"]
    mg::Float64 = params["mg"]
    lambda::Float64 = params["lambda"]

    ampo::Sum{Scaled{ComplexF64, Prod{Op}}} = OpSum()

    for n=1:N-1
        ampo += x,"S+",n,"S-",n+1
        ampo += x,"S-",n,"S+",n+1
        ampo += 2*sqrt(x)*mg*(-1)^n,"Sz",n
    end

    for n=1:N-1
        for k=1:n
            ampo += 2*l_0,"Sz",k
            ampo += l_0*(-1)^k,"Id",1
            for p=1:n
                ampo += "Sz",k,"Sz",p
                ampo += 0.5*(-1)^p,"Sz",k
                ampo += 0.5*(-1)^k,"Sz",p
                ampo += 0.25*(-1)^(p+k),"Id",1
            end
        end
    end

    for k=1:N-1
        for p=1:N-1
            ampo += lambda,"Sz",k,"Sz",p
            ampo += 0.5*lambda*(-1)^p,"Sz",k
            ampo += 0.5*lambda*(-1)^k,"Sz",p
            ampo += 0.25*lambda*(-1)^(p+k),"Id",1
        end
    end

    ampo += N*sqrt(x)*mg,"Id",1
    ampo += (l_0^2)*(N-1),"Id",1

    return ampo

end

function get_Staggered_local_charge_Opsum(site_idx::Int64)::Sum{Scaled{ComplexF64, Prod{Op}}}

    ampo::Sum{Scaled{ComplexF64, Prod{Op}}} = OpSum()

    ampo += "Sz",site_idx
    ampo += 0.5*(-1)^site_idx,"Id",site_idx

    return ampo

end

function get_SW_chiral_condensate_OpSum(N::Int64)::Sum{Scaled{ComplexF64, Prod{Op}}}

    ampo::Sum{Scaled{ComplexF64, Prod{Op}}} = OpSum()

    for n=1:N
        ampo += 1im,"S-",2*n-1,"S+",2*n
        ampo += -1im,"S+",2*n-1,"S-",2*n
    end

    return ampo

end

function get_particle_number_OpSum(N::Int64)::Sum{Scaled{ComplexF64, Prod{Op}}}

    ampo::Sum{Scaled{ComplexF64, Prod{Op}}} = OpSum()

    for n=1:N
        ampo += 1im,"S-",2*n-1,"S+",2*n
        ampo += -1im,"S+",2*n-1,"S-",2*n
    end

    ampo += N,"Id",1

    return ampo

end

function get_Schwinger_staggered_Hamiltonian_OpSum(params)::Sum{Scaled{ComplexF64, Prod{Op}}}
    
    """
    
    Inputs:
    
    N = number of spin lattice sites
    
    x = 1 / (a g)^2
    
    mu = 2 m / (g^2 a) = 2 m sqrt(x) / g
    
    l_0 = theta / (2 pi) topological theta term - background electric field
    
    lambd = penalty term for total charge
    
    Outputs: 
    
    H = the staggered Hamiltonian for the Schwinger model with type qiskit.opflow.primitive_ops.pauli_sum_op.PauliSumOp
    
    """

    N = params["N"]
    x = params["x"]
    mu = params["mu"]
    l_0 = params["l_0"]
    lambd = params["lambda"]
    
    # Build terms of H and at the end we will add them all together
    H::Sum{Scaled{ComplexF64, Prod{Op}}} = OpSum()
    
    # Kinetic term
    for n in 1:N-1
        H += 2*x,"Sy",n,"Sy",n+1 # this is 0.5 * 2 * 2 * x where the 2s account for Sy = 0.5*Y
        H += 2*x,"Sx",n,"Sx",n+1
    end
    
    # Mass term
    for n in 1:N
        H += mu * (-1)^(n-1), "Sz", n 
    end
     
    # Electric term
    for n in 1:N-1 
        for k in 1:n
            H += l_0 * (-1)^(k-1), "Id", 1   
            H += 2*l_0, "Sz", k
        end
    end
            
    # Electric term
    for n in 1:N-1
        for k in 1:n
            for k_dash in 1:n
                H += 0.25 * (-1)^(k + k_dash - 2), "Id", 1   
                H += 0.5 * (-1)^(k - 1), "Sz", k_dash 
                H += 0.5 * (-1)^(k_dash - 1), "Sz", k
            end
        end
    end
    
    # Electric term
    for n in 1:N-1
        for k in 1:n
            for k_dash in k+1:n
                H += 2, "Sz", k, "Sz", k_dash
            end
        end
    end
                
    # Penalty term
    for k in 1:N
        for k_dash in 1:N
            H += lambd * 0.25 * (-1)^(k + k_dash - 2), "Id", 1
            H += lambd * 0.5 * (-1)^(k - 1), "Sz", k_dash
            H += lambd * 0.5 * (-1)^(k_dash - 1), "Sz", k
        end
    end
    
    # Penalty term
    for k in 1:N
        for k_dash in k+1:N
            H += 2 * lambd, "Sz", k, "Sz", k_dash
        end
    end
        
    # Constant terms    
    H += lambd*N/4 + (l_0^2)*(N-1) + (1/8)*(N-2)*(N-1) + (1/4)*(N-1), "Id", 1

    return H

end

function get_staggered_site_charge_operator(site_idx)
        
    op = OpSum()
    op += "Sz", site_idx
    op += 0.5*(-1)^(site_idx-1), "Id", 1
    
    return op

end

function get_electric_field_link_operator(link, l_0)

    op = OpSum()
    for i in 1:link
        op += get_staggered_site_charge_operator(i)
    end
    op += l_0, "Id", 1

end

