function multilinear_interpolation_2d(smin, smax, orders, x, s)

    n_x = size(x,2)
    N = size(s,1)

    output = zeros(N, n_x)
    d = 2

    order_0 = orders[1]
    order_1 = orders[2]
    
    M_1 = order_0

    for j in 1:n_x

    for i in 1:N

            # (s_1, ..., s_d) : evaluation point
            s_0 = s[ i, 1 ]
            s_1 = s[ i, 2 ]

            # (sn_1, ..., sn_d) : normalized evaluation point (in [0,1] inside the grid)
            sn_0 = (s_0-smin[1])/(smax[1]-smin[1])
            sn_1 = (s_1-smin[2])/(smax[2]-smin[2])

            # q_k : index of the interval "containing" s_k
            q_0 = max( min( ifloor(sn_0 *(order_0-1)), (order_0-2) ), 0 )
            q_1 = max( min( ifloor(sn_1 *(order_1-1)), (order_1-2) ), 0 )

            # lam_k : barycentric coordinate in interval k
            lam_0 = sn_0*(order_0-1) - q_0
            lam_1 = sn_1*(order_1-1) - q_1

            # v_ij: values on vertices of hypercube "containing" the point
            v_00 = x[M_1*(q_1) + (q_0) + 1, j]
            v_10 = x[M_1*(q_1) + (q_0+1) + 1, j]
            v_01 = x[M_1*(q_1+1) + (q_0) + 1, j]
            v_11 = x[M_1*(q_1+1) + (q_0+1) + 1, j]

            # interpolated/extrapolated value
            output[i,j] = (1-lam_0)*((1-lam_1)*(v_00) + (lam_1)*(v_01)) + (lam_0)*((1-lam_1)*(v_10) + (lam_1)*(v_11))

    end
    end 
    return output
end


