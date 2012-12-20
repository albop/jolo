#import Base.isless

#isless(b::Array, a::Number) = a*ones(size(b)) .< b
#isless(a::Array{Int64,1},b::Int64) = a.< ones(size(b));
#isless(Array{Int64,1},Int64)


load("mlininterp")


NN = 1000000

orders = [100,100]

N = prod(orders)

grid = randn( N, 2)
fine_grid = randn( NN, 2 )

fun(x) = [ x[:,1].*x[:,2]  x[:,1] ]

values = fun(grid)

smin = [-1,-1]
smax = [1,1]

times = []
for i = 1:10
    tic()
    t = multilinear_interpolation(smin,smax,orders,values, fine_grid)
    toc()
    #times = [times, elapsed]
end


#print(times)
