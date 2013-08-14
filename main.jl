require("decision_rules.jl")
require("yaml_import")
require("solve_model")


## import the model

model = yaml_import("rbc")


## construct approximation space

# boundaries for the state-space
smin = [1-0.001, 5]
smax = [1+0.001, 10]
orders = [10, 50]



# nodes and weight to compute expectations
nodes = zeros( (1,1) )  # all expected shocks are zero
weights = [1.0]

approx = ApproximationSpace(smin,smax,orders,nodes,weights)


## initial decision rule

# construct cartesian grid
grid = mlinspace(smin,smax,orders)
s_ss = model.calibration["states"]
x_ss = model.calibration["controls"]
X_s = [[ 1.25268335 -0.02185392 ;  0.20023482 -0.00591657]];
N = prod(approx.orders)
init = repmat(x_ss', N, 1) + (grid - repmat(s_ss', N, 1) ) * X_s'

# create decision rule object
dr = DecisionRule(smin,smax,orders,init)

## solve the model


tic()
dr = solve_model(model, approx, init)
toc()

tic()
dr = solve_model(model, approx, init)
toc()

tic()
# use an initial d.r.
dr = solve_model(model, approx, dr)
toc()

