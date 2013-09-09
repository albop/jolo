#require("ndgrid.jl")
#require("mlininterp.jl")
#require("interpolation.jl")


type DecisionRule

    smin::Array{Float64,1}
    smax::Array{Float64,1}
    orders::Array{Int64,1}
    grid::Array{Float64,2}
    values::Array{Float64,2}

end


function evaluate(dr::DecisionRule, s::Array{Float64,2}) # return type : Array{Float64,2}
    if length(dr.orders) == 2
        n_x = size(dr.values,2)
        N = size(s,1)
        output = multilinear_interpolation_2d(dr.smin, dr.smax, dr.orders, dr.values, s)
        return output
    else
        return multilinear_interpolation(dr.smin, dr.smax, dr.orders, dr.values, s)
    end
end

function mlinspace( smin, smax, orders)
    d = size(smin,1)
    if d == 1
        grid = linspace(smin[1],smax[1],orders[1])
    elseif d == 2
        (x,y) = ndgrid(
                linspace(smin[1],smax[1],orders[1]),
                linspace(smin[2],smax[2],orders[2])
               )
        grid = hcat( x[:] , y[:] )
    elseif d == 3
        grid = ndgrid(
                linspace(smin[3],smax[3],orders[3])
               )
    elseif d == 4
        grid = ndgrid(
                linspace(smin[1],smax[1],orders[1]),
                linspace(smin[2],smax[2],orders[2]),
                linspace(smin[3],smax[3],orders[3]),
                linspace(smin[4],smax[4],orders[4])
               )
    end
    return grid
end


function DecisionRule(smin::Array{Float64,1}, smax::Array{Float64,1}, orders::Array{Int64,1})
    grid = mlinspace(smin, smax, orders)
    # we fill the d.r. with dummy values
    n = size(grid,1)
    values = zeros(n,1)
    dr =  DecisionRule(smin, smax, orders, grid, values)
    return dr
end


function DecisionRule(smin::Array{Float64,1}, smax::Array{Float64,1}, orders::Array{Int64,1}, values::Array{Float64,2})
    grid = mlinspace(smin, smax, orders)
    dr =  DecisionRule(smin, smax, orders, grid, values)
    return dr
end
