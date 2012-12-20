load("ndgrid.jl")
load("mlininterp.jl")

type DecisionRule

    smin::Array{Float64,1}
    smax::Array{Float64,1}
    orders::Array{Int64,1}
    grid::Array{Float64,2}
    values::Array{Float64,2}

end

function evaluate(dr::DecisionRule, s::Array{Float64,2}) # return type : Array{Float64,2}
        return multilinear_interpolation(dr.smin, dr.smax, dr.orders, dr.values, s)
end


function DecisionRule(smin::Array{Float64,1}, smax::Array{Float64,1}, orders::Array{Int64,1})
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
                linspace(smin[1],smax[1],orders[1]),
                linspace(smin[2],smax[2],orders[2]),
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

    # we fill the d.r. with dummy values
    n = size(grid,1)
    values = zeros(n,1)
    dr =  DecisionRule(smin, smax, orders, grid, values)

    return dr

end

