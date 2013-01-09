

function transition(states_p, controls_p, auxiliary_p, shocks, p)

    n = size(states_p,1)

    val = zeros( n, 2 )
    val[:,1] = p[7].*states_p[:,1] + p[8].*(-p[7] + 1) + shocks[:,1]
    a =  controls_p[:,1]
    b = states_p[:,2].*(-p[5] + 1)

    val[:,2] = a + b

    return val

end
    
function arbitrage(states, controls, auxiliary, states_f, controls_f, auxiliary_f, p)

    n = size(states,1)

    val = zeros( n, 2 )
    val[:,1] = -p[1].*(auxiliary[:,1]./auxiliary_f[:,1]).^p[2].*(-p[5] + auxiliary_f[:,2] + 1) + 1
    val[:,2] = -p[4].*auxiliary[:,1].^p[2].*controls[:,2].^p[3] + auxiliary[:,3]


    return val

end
    
function auxiliary(states, controls, p)

    n = size(states,1)

    val = zeros( n, 3 )
    val[:,1] = -controls[:,1] + states[:,2].^p[6].*controls[:,2].^(-p[6] + 1).*states[:,1]
    val[:,2] = p[6].*states[:,1].*(controls[:,2]./states[:,2]).^(-p[6] + 1)
    val[:,3] = states[:,1].*(states[:,2]./controls[:,2]).^p[6].*(-p[6] + 1)


    return val

end
    

model = {
    "states" => ["z","k"],
    "controls" => ["i","n"],
    "auxiliaries" => ["c","rk","w"],
    "parameters" => ["beta","sigma","eta","chi","delta","alpha","rho","zbar"],
    "shocks" => ["e_z"],
    "transition" => transition,
    "arbitrage" => arbitrage,
    "auxiliary" => auxiliary,
    "params" => [0.99, 1, 1, 8.04277481517292, 0.025, 0.33, 0.8, 1],
    "s_ss" => [1, 9.35497829014598],
    "x_ss" => [0.233874457253650, 0.33],
    "a_ss" => [0.233874457253650, 0.33]
}
