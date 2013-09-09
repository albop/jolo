#isless(a::Number, b::Array) = a*ones(size(b)) .< b 

import Base.isless

#    using Base
#    isless(b::Array, a::Number) = a*ones(size(b)) .< b
#    isless(a::Number, b::Array) = a*ones(size(b)) .< b
#function isless(a::Array{Int64,1},b::Int64)
#    lhs = a
#    rhs = b*ones(Int64,size(a))
#    print("\nlhs\n")
#    print(lhs)
#    print("\nrhs\n")
##    print(rhs)
#    res = lhs .< rhs
#    return res
#end
    isless(a::Array{Int64,1},b::Int64) = a .< b*ones(Int64,size(a))
#    isless(Array{Int64,1},Int64)


function multilinear_interpolation(smin,smax,orders,x,y)

	d = length(orders)
	N = size(y,1)
        nx = size(x,2)
	qq = zeros(Int,N,d)
	mm = zeros(N,d)

    
    ss = zeros(size(y))
    for i = 1:d
        ss[:,i] = (y[:,i]-smin[i])/(smax[i]-smin[i])
    end

	for i=1:d
		s = ss[:,i]
		n = orders[i]
		delta = 1/(n-1)
		r = s/delta
		q = ifloor(r) # why doesn't floor return an integer ?
		q = (q<0) .* 0 + (q>= n-2) .* (n-2) + (0<=q).*(q<n-2).*q
                q = min(q, orders[i]-1)
                q = max(q, 0)
		m = r-q
		mm[:,i] = m
		qq[:,i] = q
	end
#	qq = qq+1


	(b,g) = strange_construction( x, qq, orders)

    z = b + recursive_evaluation(d,[], mm, g)

    return z

end

function recursive_evaluation(d,ind,mm,g)
    if length(ind) == d
        e = g[:,:,ind...]
    else
        j = length(ind)
        ind1 = [ind, 1]
        ind2 = [ind, 2]
        nx = size(g,2)
        tmm = repmat( mm[:,j+1], 1, nx )
        a =  recursive_evaluation(d,ind1,mm,g)
        b = recursive_evaluation(d,ind2,mm,g)
        e = (1-tmm) .* a + tmm .* b
    end

    return e
end

function strange_construction(a, q, dims)

	N = size(q,1)
    nx = size(a,2) 

    d = length(dims)
	k = length(dims)

#	q = q-1
	cdims = cumprod(dims)

    lin_q = q[:,1]
	for i = 2:k
		lin_q = lin_q + q[:,i]*cdims[i-1]
	end

    cart_prod = binary_cart_prod(d)

	lin_cp = cart_prod[:,1]
	for i = 2:k
		lin_cp = lin_cp + cart_prod[:,i]*cdims[i-1]
	end

	lin_q = lin_q + 1
	b = a[lin_q,:]

	g = zeros( N, nx, size(cart_prod,1) )
	for i = 1:size(cart_prod,1)
		g[:,:,i] = a[lin_q+lin_cp[i],:] - b
    end

    ndims = [[ N, nx ], 2*ones(Int64,d)]
    g = reshape(g,ndims...)

	return b,g

end


function binary_cart_prod(d)
    if d == 1
        return [0, 1]
    else
        p = binary_cart_prod(d-1)
        n = size(p,1)
        a = [p zeros(Int,n) ]
        b = [p ones(Int,n) ]
        return [ a, b ]
    end
end

