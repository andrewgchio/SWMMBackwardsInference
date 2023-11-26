module SWMMBackwardsInference

# Enable keyword arguments in function definitions/structs
import Base.@kwdef

# Import optimization software
import JuMP
import Ipopt

# Import graph structures and visualization
import Graphs: DiGraph
import MetaGraphsNext: MetaGraph, label_for, code_for, vertices, edges, 
                        inneighbors, outneighbors, induced_subgraph,
                        nv, ne
import GraphMakie: graphplot
import CairoMakie
import GeometryBasics: Point

# Config Files
import ConfParser: ConfParse, parse_conf!, retrieve

# Date Utils
import Dates: Hour

# DataFrame / CSV Utils
import DataFrames: DataFrame
import CSV: read

# Stormwater Network
export construct_graph, SWNode, SWEdge
export get_local_branch, get_upstream_nodes, nv, ne, label_for, code_for
export gen_SWnodes, gen_SWedges, gen_SWvnames, gen_SWenames
export print_graph, draw_graph

# Simulation and Inference
export all_source_inference, single_source_inference
export differentiable_swmm

# Visualization support
export plot_edgeflows, plot_inflow

# Utils
export dict2df

# Load pyswmm and hymo (SWMM parser) in path
import PyCall
const hymo = PyCall.PyNULL()
const plt  = PyCall.PyNULL()

# Globals
MAXTRIALS   = 8
Ω           = 0.5
FUDGE       = 0.0001
A_Smin      = 12.566
MAXVELOCITY = 50

LINEAR_SOLVER = "ma57"
PRINT_LEVEL   = 5
TOL           = 0.001
MAX_CPU_TIME  = 600
Qmin          = 0.0
Qmax          = 1.0

function __init__()
    # Load python libraries: hymo
    hymopath = "./hymo"
    pushfirst!(PyCall.PyVector(PyCall.pyimport("sys")."path"), hymopath)
    copy!(hymo, PyCall.pyimport("hymo"))

    # Load matplotlib.pyplot
    copy!(plt, PyCall.pyimport("matplotlib.pyplot"))

    # Load Configuration files
    conf = ConfParse("./config/config.ini")
    parse_conf!(conf)
    global MAXTRIALS   = parse(Int64, retrieve(conf, "routing", "MAXTRIALS"))
    global Ω           = parse(Float64, retrieve(conf, "routing", "OMEGA"))
    global FUDGE       = parse(Float64, retrieve(conf, "routing", "FUDGE"))
    global A_Smin      = parse(Float64, retrieve(conf, "routing", "A_Smin"))
    global MAXVELOCITY = parse(Float64, retrieve(conf, "routing", "MAXVELOCITY"))

    global LINEAR_SOLVER = retrieve(conf, "optimization", "LINEAR_SOLVER")
    global PRINT_LEVEL   = parse(Int64, retrieve(conf, "optimization", "PRINT_LEVEL"))
    global TOL           = parse(Float64, retrieve(conf, "optimization", "TOL"))
    global MAX_CPU_TIME  = parse(Float64, retrieve(conf, "optimization", "MAX_CPU_TIME"))
    global Qmin          = parse(Float64, retrieve(conf, "optimization", "Qmin"))
    global Qmax          = parse(Float64, retrieve(conf, "optimization", "Qmax"))

end

# General
include("Globals.jl")
include("Utils.jl")

# Network
include("StormwaterNetwork.jl")
include("XSect.jl")
include("XSectdat.jl")

# Inference
include("DynwaveRouting.jl")
include("Inference.jl")

end
