
function yaml_import( filename::String )

    command = `dolo-julia $filename.yaml`
    run(command)
    
    path = pwd()
    txt = string( pwd(), "/$(filename)_model.jl")

    model = evalfile(txt)

    mymodel = Model( model["symbols"], model["functions"], model["calibration"])

    return mymodel

end
