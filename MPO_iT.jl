function get_SW_local_charge_OpSum(site_idx)

    ampo = OpSum()

    ampo += "Sz",2*site_idx-1
    ampo += "Sz",2*site_idx
    ampo += -"Id",1

    return ampo

end

function get_SW_OpSum(params)

    """
    N = number of physical sites
    """

    N = params["N"]
    l_0 = params["l_0"]
    x = params["x"]
    mg = params["mg"]
    r = params["r"]
    lambda = params["lambda"]

    A = 1im*x*(r - 1)
    B = 1im*x*(r + 1)
    C = 2*1im*(sqrt(x)*mg + x*r)
    F = (N-1)*l_0^2 + N*(N-1)/4 + lambda*N/2

    ampo = OpSum()

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

function get_SW_OpSum_alternative(params)

    """
    N = number of physical sites
    """

    N = params["N"]
    l_0 = params["l_0"]
    x = params["x"]
    mg = params["mg"]
    r = params["r"]
    lambda = params["lambda"]

    A = 1im*x*(r - 1)
    B = 1im*x*(r + 1)
    C = 2*1im*(sqrt(x)*mg + x*r)
    F = (N-1)*l_0^2 + N*(N-1)/4 + lambda*N/2

    ampo = OpSum()

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

function get_Ising_OpSum(N, J, g_z, g_x)
        
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

function get_MPO_from_OpSum(OpSum, sites)

    return MPO(OpSum, sites)

end
