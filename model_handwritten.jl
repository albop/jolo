import Base

function g(s,x,e,p)

    a = s[:,1]
    k = s[:,2]
    i = x[:,1]

    delta = p[1]
    rho = p[2]

    K = k*(1-delta) + i
    A = a*rho + e
    
    S = hcat(A, K)

    return S
                
end

function f(s,x,S,X,e,p)

    a = s[:,1]
    k = s[:,2]
    i = x[:,1]

    A = S[:,1]
    K = S[:,2]
    I = X[:,1]

    delta = p[1]
    gamma = p[3]
    theta = p[4]
    beta = p[5]

    c = exp(a).*k.^theta - i
    C = exp(A).*K.^theta - I

    res = beta*(C./c).^(-gamma).*( (1-delta) + exp(A).*theta.*K.^(theta-1) ) - 1

    return res

end


function steady_state(p)
    (delta, rho, gamma, theta,beta) = p
    a_ss = 0.0
    k_ss = ((1/beta-(1-delta))/theta)^(1/(theta-1))
    i_ss = delta*k_ss

    s_ss = [a_ss, k_ss]
    x_ss = [i_ss]

    return (s_ss, x_ss)
    
end


params = [0.1, 0.9, 4, 0.3,  0.96]
(s_ss, x_ss) = steady_state(params)

type Model
        controls
        states
        auxiliaries
        parameters
        shocks
        g:: Function #(Any, Any, Any , Array{Float64,1})
        f:: Function
        a:: Function
        params:: Array{Float64,1}
        s_ss
        x_ss
end

model = Model( 
        ["i"],
        ["a", "k"],
        ["delta", "rho", "gamma", "theta", "beta"],
        g,
        f,
        params,
        s_ss,
        x_ss
)

nodes = [0.0]
weights = [1.0]


