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

function get_z_configuration(psi::MPS, sites)::Vector{ComplexF64}

    N_spin::Int64 = length(psi)

    z_configuration_list::Vector{ComplexF64} = []
    
    for i=1:N_spin
    
        z_OpSum::Sum{Scaled{ComplexF64, Prod{Op}}} = get_SW_local_z_OpSum(i)
    
        z_mpo::MPO = get_MPO_from_OpSum(z_OpSum, sites)
    
        append!(z_configuration_list, inner(psi', z_mpo, psi))
    end

    return z_configuration_list

end

function get_SW_electric_field_configuration(charge_configuration_list::Vector, l_0::Float64)::Vector{ComplexF64}

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

function get_Schwinger_Wilson_OpSum(params)::Sum{Scaled{ComplexF64, Prod{Op}}}

    """
    N = number of physical sites
    """

    N::Int64 = params["N"]
    l_0::Float64 = params["l_0"]
    x::Float64 = params["x"]
    mg::Float64 = params["mg"]
    r::Float64 = params["r"]
    lambda::Float64 = params["lambda"]

    H::Sum{Scaled{ComplexF64, Prod{Op}}} = OpSum()

    # These have been commented and we follow the S+, S- convention to make it obvious that Sz is conserved and then we can use conserved quantum number MPS

    # for n=1:N-1
    #     H += 0.5*x*(r-1),"X",2*n,"X",2*n+1
    #     H += 0.5*x*(r-1),"Y",2*n,"Y",2*n+1
    # end

    # for n=1:N-1
    #     H += 0.5*x*(r+1),"X",2*n,"X",2*n+1
    #     H += 0.5*x*(r+1),"Y",2*n,"Y",2*n+1
    # end

    # for n=1:N
    #     H += (mg*sqrt(x) + x*r),"X",2*n-1,"X",2*n
    #     H += (mg*sqrt(x) + x*r),"Y",2*n-1,"Y",2*n
    # end

    for n=1:N-1
        H += x*(r-1),"S-",2*n,"S+",2*n+1
        H += x*(r-1),"S+",2*n,"S-",2*n+1
    end

    for n=1:N-1
        H += x*(r+1),"S-",2*n,"S+",2*n+1
        H += x*(r+1),"S+",2*n,"S-",2*n+1
    end

    for n=1:N
        H += 2*(mg*sqrt(x) + x*r),"S-",2*n-1,"S+",2*n
        H += 2*(mg*sqrt(x) + x*r),"S+",2*n-1,"S-",2*n
    end

    for n=1:2*N-2
        H += l_0*(N-ceil(n/2)),"Z",n
    end

    for n=1:2*N
        for k in n+1:2*N
            H += 0.5*(N-ceil(k/2)+lambda),"Z",n,"Z",k
        end
    end

    H += ((l_0^2)*(N-1) + 0.25*N*(N-1) + lambda*N/2),"Id",1

    return H

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

function get_MPO_from_OpSum(OpSum, sites)::MPO

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
        ampo += "S-",2*n-1,"S+",2*n
        ampo += "S+",2*n-1,"S-",2*n
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
    # for n in 1:N-1
    #     H += 2*x,"Sy",n,"Sy",n+1 # this is 0.5 * 2 * 2 * x where the 2s account for Sy = 0.5*Y
    #     H += 2*x,"Sx",n,"Sx",n+1
    # end
    for n in 1:N-1
        H += x,"S+",n,"S-",n+1 # this is 0.5 * 2 * 2 * x where the 2s account for Sy = 0.5*Y
        H += x,"S-",n,"S+",n+1
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

function get_Schwinger_staggered_Hamiltonian_MPO(params, sites)::MPO

    N = params["N"]
    x = params["x"]
    mu = params["mu"]
    l_0 = params["l_0"]

    ## construct the MPO QN subspace ##
    H = MPO(N)
    Link = []
    ## 1 generate all the link space
    for n = 1:N-1
        ln = Index(QN() => 1, QN("Sz", -2) => 1, QN("Sz", 2) => 1, QN("Sz", 0) => 1, QN("Sz", 0) => 1; tags = join(["Link,l=", string(n)]))
        push!(Link, ln)
    end
    ## 1 generate the local tensor subspace of MPO, where lr is the right link and ll is the left link
    lr = Link[1]
    lp = prime(sites[1])
    lpdag = dag(sites[1])
    H[1] = ITensor(lr, lp, lpdag)
    h02 = 0.5*mu + 0.25 + (0.25+l_0)*(N-1)
    ## (1) H[1][1] : h_{2}*\sigma^{z} + const*I
    H[1][lr=>1,lp=>1,lpdag=>1] = h02*1 + (0.5*l_0+0.125*N)*N+(N-1)*l_0*l_0
    H[1][lr=>1,lp=>2,lpdag=>2] = -1*h02 + (0.5*l_0+0.125*N)*N+(N-1)*l_0*l_0
    ## (2) H[1][2] : \sigma^{+}
    H[1][lr=>2,lp=>1,lpdag=>2] = 1.
    ## (3) H[1][3] : \sigma^{-}
    H[1][lr=>3,lp=>2,lpdag=>1] = 1.
    ## (4) H[1][4] : \sigma^{z}
    H[1][lr=>4,lp=>1,lpdag=>1] = 1.
    H[1][lr=>4,lp=>2,lpdag=>2] = -1.
    ## (5) H[1][5] = I
    H[1][lr=>5,lp=>1,lpdag=>1] = 1.
    H[1][lr=>5,lp=>2,lpdag=>2] = 1.

    for n = 2:N-1
        ll = dag(Link[n-1])
        lr = Link[n]
        lp = prime(sites[n])
        lpdag = dag(sites[n])
        H[n] = ITensor(ll, lr, lp, lpdag)
        hn1 = 0.5*(N-n)
        hn2 = -0.5*mu*(-1)^n + 0.125*(1-(-1)^n) + (0.25+l_0)*(N-n)
        ## (1) H[n][1,1] = I
        H[n][ll=>1,lr=>1,lp=>1,lpdag=>1] = 1.
        H[n][ll=>1,lr=>1,lp=>2,lpdag=>2] = 1.
        ## (2) H[n][2,1] = x * \sigma^{-}
        H[n][ll=>2,lr=>1,lp=>2,lpdag=>1] = x
        ## (3) H[n][3,1] = x * \sigma^{+}
        H[n][ll=>3,lr=>1,lp=>1,lpdag=>2] = x
        ## (4) H[n][4,1] = h_{1} * \sigma^{z}
        H[n][ll=>4,lr=>1,lp=>1,lpdag=>1] = hn1
        H[n][ll=>4,lr=>1,lp=>2,lpdag=>2] = -hn1
        ## (5) H[n][5,1] = h_{2} * \sigma^{z}
        H[n][ll=>5,lr=>1,lp=>1,lpdag=>1] = hn2
        H[n][ll=>5,lr=>1,lp=>2,lpdag=>2] = -hn2
        ## (6) H[n][5,2] = \sigma^{+}
        H[n][ll=>5,lr=>2,lp=>1,lpdag=>2] = 1.
        ## (7) H[n][5,3] = \sigma^{-}
        H[n][ll=>5,lr=>3,lp=>2,lpdag=>1] = 1.
        ## (8) H[n][5,4] = \sigma^{z}
        H[n][ll=>5,lr=>4,lp=>1,lpdag=>1] = 1.
        H[n][ll=>5,lr=>4,lp=>2,lpdag=>2] = -1.
        ## (9) H[n][5,5] = Id
        H[n][ll=>5,lr=>5,lp=>1,lpdag=>1] = 1.
        H[n][ll=>5,lr=>5,lp=>2,lpdag=>2] = 1.
        ## (10) H[n][4,4] = Id
        H[n][ll=>4,lr=>4,lp=>1,lpdag=>1] = 1.
        H[n][ll=>4,lr=>4,lp=>2,lpdag=>2] = 1.
    end

    ll = dag(Link[N-1])
    lp = prime(sites[N])
    lpdag = dag(sites[N])
    H[N] = ITensor(ll, lp, lpdag)
    hN2 = -0.5*mu
    ## (1) H[N][1] = I
    H[N][ll=>1,lp=>1,lpdag=>1] = 1.
    H[N][ll=>1,lp=>2,lpdag=>2] = 1.
    ## (2) H[N][2] = x * \sigma^{-}
    H[N][ll=>2,lp=>2,lpdag=>1] = x
    ## (3) H[N][3] = x * \sigma^{+}
    H[N][ll=>3,lp=>1,lpdag=>2] = x
    ## (4) H[N][5] = h_{2} * \sigma^{z}
    H[N][ll=>5,lp=>1,lpdag=>1] = hN2
    H[N][ll=>5,lp=>2,lpdag=>2] = -hN2

    return H

end

function get_staggered_particle_number_MPO(N, sites)

    res = OpSum()

    for n in 1:N
        res += 0.5*(-1)^(n-1), "Z", n 
    end

    res += 0.5*N,"Id",1

    return MPO(res, sites)

end

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
