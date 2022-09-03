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