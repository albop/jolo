

function yaml_import( filename::String )
    
    command = `dolo-julia $filename.yaml`
    run(command)
    
    txt = string( pwd(), "/$(filename)_model.jl")
    include( txt )

symbols = Dict{String, Array{String,1}}()  
symbols["states"] = model["symbols"]["states"]
symbols["controls"] = model["symbols"]["controls"]
symbols["auxiliary"] = model["symbols"]["parameters"]
symbols["parameters"] = model["symbols"]["parameters"]
symbols["shocks"] = model["symbols"]["shocks"]

functions = Dict{String, Function}()
functions["arbitrage"] = model["functions"]["arbitrage"]
functions["transition"] = model["functions"]["transition"]
functions["auxiliary"] = model["functions"]["auxiliary"]

calib = Dict{String, Array{Float64,1}}()
calib["states"] = model["calibration"]["steady_state"]["states"]
calib["controls"] = model["calibration"]["steady_state"]["controls"]
calib["auxiliaries"] = model["calibration"]["steady_state"]["auxiliary"]
calib["parameters"] = model["calibration"]["parameters"]

mymodel = Model( symbols, functions, calib )

return mymodel

end
