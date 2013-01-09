load("rbc_model.jl")
load("decision_rules.jl")
load("newton.jl")

smin = [1-0.001, 8]
smax = [1+0.001, 10]
orders = [3, 10]

nodes = zeros( (1,1) ) 
weights = [1.0]

dr = DecisionRule(smin, smax, orders)

grid = dr.grid
N = size(grid,1)

s_ss = model["s_ss"]'
x_ss = model["x_ss"]'
X_s = [[ 1.25268335 -0.02185392 ;  0.20023482 -0.00591657]];

init = repmat(x_ss, N, 1) + (grid - repmat(s_ss, N, 1) ) * X_s'

type ApproximationSpace
    smin:: Array{Float64, 1}
    smax:: Array{Float64, 1}
    orders:: Array{Int64, 1}
    nodes:: Array{Float64, 2}
    weights:: Array{Float64, 1}
end

type Model
        states
        controls
        auxiliaries
        parameters
        shocks
        g:: Function #(Any, Any, Any , Array{Float64,1})
        f:: Function
        a:: Function
        params:: Array{Float64,1}
        s_ss
        x_ss
        a_ss
end

rbc_model = Model(
        model["states"],
        model["controls"],
        model["auxiliaries"],
        model["parameters"],
        model["shocks"],
        model["transition"],
        model["arbitrage"],
        model["auxiliary"],
        model["params"],
        model["s_ss"],
        model["x_ss"],
        model["a_ss"]
)


approx = ApproximationSpace(smin,smax,orders,nodes,weights)




function step_residuals(f, g, aux, s::Array{Float64,2}, x::Array{Float64,2}, p::Array{Float64,1}, dr::DecisionRule, nodes::Array{Float64,2}, weights::Array{Float64,1})

   
    n = (size(s, 1))
    res = zeros(size(x))
    a = aux(s,x,p)
    for i = 1:size(weights,1)
        e = nodes[i,:]
        e = repmat(e, n, 1)
        S = g(s,x,a,e,p)
        X = evaluate(dr,S)
        A = aux(S,X,p)
        ff =  f(s,x,a,S,X,A,p)
        res += ff*weights[i]
    end
    return res
end



function solve_model(model::Model, approx::ApproximationSpace, x0::Array{Float64,2})

    f = model.f
    g = model.g
    a = model.a
    
    p = model.params

 
    verbose = false
    smin = approx.smin
    smax = approx.smax
    orders = approx.orders
    nodes = approx.nodes
    weights = approx.weights
    
    maxit = 500
    err = 100
    tol = 1e-6
    it = 0

    dr = DecisionRule(smin,smax,orders)

    print("Starting iterations\n")
    while (err>tol)&&(it<maxit)
        
            dr.values = x0
    
            fobj(t::Array{Float64,2}) = step_residuals(f, g, a, dr.grid, t, p, dr, nodes, weights)

            res = fobj(x0)
            (x,nit) = serial_solver(fobj, x0, 5)
        
            err = max(max(abs(x-x0)))
            it +=1

            x0 = x
        if verbose    
            print((it, err, nit))
            print("\n")
        end
    end

end


tic()
solve_model(rbc_model, approx, init)
toc()

tic()
solve_model(rbc_model, approx, init)
toc()
