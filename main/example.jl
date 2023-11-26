import Distributions: Normal, pdf

using SWMMBackwardsInference

function main()

    # Example EPA SWMM input file
    FSWMM_INP = "./data/swmminp/tutorial/tutorial.inp"

    # Initialize the stormwater graph
    G = construct_graph(FSWMM_INP)

    # Set time step and range
    Δt = 30
    T = 0:Δt:2*60*60

    # Initialize a stormwater anomaly Qdwf
    src = "J1"
    st, et, q = 60, 100, 0.2
    Qdwf = zeros(length(T))
    μ,σ = st + 0.5*(et-st), (et-st)/3
    for (i,t) in enumerate(T); Qdwf[i] = pdf(Normal(μ,σ), t) end
    Qdwf /= (maximum(Qdwf) / q)

    # Run forwards simulation
    Q = differentiable_swmm(G, T, Δt, src, Qdwf...)
    println(Q)

    # # Take "J2" as a sensor
    S = ["J2-J4"]
    Qdf = dict2df(Q, T, S)
    Qobs = Qdf[:, ["time", S...]]
    println(Qobs)

    # Backwards inference can be done on the single source - but this is
    # repeated in when iterating across all sources
    # Qinf, objval, status, pstatus = single_source_inference(G, src, S, Qobs, T, Δt)

    # # And try backwards inference on all sources
    Qinfs = all_source_inference(G, S, Qobs, T, Δt)
    println("Results from all source inference")
    println(Qinfs)
    println()

end

main()