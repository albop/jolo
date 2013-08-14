require("newton.jl")

type Model
    symbols:: Dict{String,Array{String,1}}
    functions:: Dict{String,Function}
    calibration:: Dict{String,Array{Float64,1}}
end


type ApproximationSpace
    smin:: Array{Float64, 1}
    smax:: Array{Float64, 1}
    orders:: Array{Int64, 1}
    nodes:: Array{Float64, 2}
    weights:: Array{Float64, 1}
end

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


function solve_model(model::Model, approx::ApproximationSpace, dr::DecisionRule)

    grid = mlinspace(approx.smin, approx.smax, approx.orders)
    init = evaluate(dr, grid)
    dr = solve_model(model, approx, init)
    return dr

end

function solve_model(model::Model, approx::ApproximationSpace, x0::Array{Float64,2}, verbose=false)

    f = model.functions["arbitrage"]
    g = model.functions["transition"]
    a = model.functions["auxiliary"]
    
    p = model.calibration["parameters"]

 
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
    grid = dr.grid

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
    println("Finished in ", it, " iterations.")
    return dr
end
