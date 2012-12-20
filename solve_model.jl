load("model_test.jl")
load("decision_rules.jl")
load("newton.jl")

smin = [-0.001, 1]
smax = [0.001, 10]
orders = [2, 10]

dr = DecisionRule(smin, smax, orders)

grid = dr.grid

s_ss = model.s_ss
x_ss = model.x_ss

values = hcat( grid[:,2]*0.1 ) # stupid starting value for investment
#reshape(values, (size( values,1),1) )

function step_residuals(model::Model, s::Array{Float64,2}, x::Array{Float64,2}, dr::DecisionRule, nodes::Array{Float64,1}, weights::Array{Float64,1})
    p = model.params
    f = model.f
    g = model.g
    res = zeros(size(x))
    for i = 1:size(weights,1)
        e = nodes[i]
        S = g(s,x,e,p)
        X = evaluate(dr,S)
        res += f(s,x,S,X,e,p)*weights[i]
    end
    return res
end


#maxit = 0
#dr.values = x0
#    
# fobj(t::Array{Float64,2}) = step_residuals(model, dr.grid, t, dr, nodes, weights)
#    (x,nit) = serial_solver(fobj, x0, 5)

type ApproximationSpace
    smin:: Array{Float64, 1}
    smax:: Array{Float64, 1}
    orders:: Array{Int64, 1}
    nodes:: Array{Float64, 1}
    weights:: Array{Float64, 1}
end

approx = ApproximationSpace(smin,smax,orders,nodes,weights)

function solve_model(model::Model, approx::ApproximationSpace, x0::Array{Float64,2})

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
    
            fobj(t::Array{Float64,2}) = step_residuals(model, dr.grid, t, dr, nodes, weights)
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
solve_model(model, approx, values)
toc()

tic()
solve_model(model, approx, values)
toc()
