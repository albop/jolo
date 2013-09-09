

function transition(states_p, controls_p, auxiliary_p, shocks, p)

    n = size(states_p,1)

    val = zeros( n, 2 )
    val[:,1] = p[7].*states_p[:,1] + p[8].*(-p[7] + 1) + shocks[:,1]
    val[:,2] = controls_p[:,1] + states_p[:,2].*(-p[5] + 1)


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
    

functions = {
	"transition" => transition,
	"arbitrage" => arbitrage,
	"auxiliary" => auxiliary,
}

steady_state = {
	"states" => [1.0, 9.35497829015],
	"controls" => [0.233874457254, 0.33],
	"auxiliary" => [0.761183686556, 0.035101010101, 2.0202695647],
}

symbols = {
	"states" => ["z","k"],
	"controls" => ["i","n"],
	"auxiliary" => ["c","rk","w"],
	"shocks" => ["e_z"],
	"parameters" => ["beta","sigma","eta","chi","delta","alpha","rho","zbar"],
}

calibration = {
    "steady_state" => steady_state,
    "parameters" => [0.99, 1.0, 1.0, 8.04277481517, 0.025, 0.33, 0.8, 1.0],
    "covariances" => [[ 0.0015]]
}

model = {
    "symbols" => symbols,
    "functions" => functions,
    "calibration" => calibration
}
