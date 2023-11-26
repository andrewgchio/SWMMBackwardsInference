# Inference.jl
#
# Backwards inference models

function all_source_inference(G, S, Qobs, T, Δt)
    Q_dwfs = Dict()

    for src in gen_SWvnames(G)

        # ignore outfalls
        if G[src].nodetype == "OUTFALL"; continue end
                                    
        Q_dwf, objval, status, pstatus = 
            single_source_inference(G, src, S, Qobs, T, Δt)
                
        println("Value = $(objval) with q = (OMITTED) [$(status) - $(pstatus)]")
        println()
                
        # Save values of Q for plotting later
        Q_dwfs[src] = (objval, JuMP.value.(Q_dwf).data)
                
    end

    return Q_dwfs
end

function single_source_inference(G, src, S, Qobs, T, Δt)

    println("Starting inference with src=$(src)")

    # Initialize the model
    model = JuMP.Model(Ipopt.Optimizer)

    # Set a few optimizer attributes
    JuMP.set_optimizer_attribute(model, "linear_solver", LINEAR_SOLVER) 
    JuMP.set_optimizer_attribute(model, "print_level", PRINT_LEVEL)
    JuMP.set_optimizer_attribute(model, "tol", TOL)

    # JuMP.set_optimizer_attribute(model, "max_iter", MAX_ITER)
    JuMP.set_optimizer_attribute(model, "max_cpu_time", MAX_CPU_TIME) 
    
    # Declare the Qdwf variable to optimize, ranging from 0.0 to 1.0
    JuMP.@variable(model, Qmin <= Q_dwf[items=T] <= Qmax)
        
    function wrapper(Q_dwf::TR...) where {TR <: Real}
        Q = differentiable_swmm(G, T,Δt, src, Q_dwf...)

        # Multiply by Q factor to "normalize"
        return sum( (Qfactor*Q[i,s]-(Qfactor*Qobs[i,s])) ^ 2 
                for i in 1:length(T) for s in S)
    end
        
    JuMP.register(model, :wrapper, length(T), wrapper; autodiff = true)
    JuMP.@NLobjective(model, Min, wrapper(Q_dwf...))
    
    # println(model)
    JuMP.optimize!(model)
        
    return Q_dwf, JuMP.objective_value(model), 
            JuMP.termination_status(model), JuMP.primal_status(model)
end
    
function my_diff_simulate_df_wrapper(; src, q, st, et)
    Q_dwf = zeros(length(T))
#     Q_dwf[time2index(st,T):time2index(et,T)] .= q
    for (i,t) in enumerate(T)
        Q_dwf[i] = my_swmm_dwf(t+Δt, q=q, st=st, et=et)
    end
    
    Qdiff_simu = my_diff_simulate(src,Q_dwf...)
    return dict2df(Qdiff_simu, T, Ename)
end

################################################################################
# Visualization
################################################################################

function plot_inflow(T, Q; label=nothing, ϵ=0.025, ax=plt, 
        title="", xlabel="", ylabel="", ylimits=[0.0,0.5])

    ax.plot(T, Q, label=label)

    ax.set_title(title) 
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # PyPlot.grid(grid) 
    ax.set_ylim(ylimits...)

    ax.legend(bbox_to_anchor=[1.05,1],loc=2,borderaxespad=0)

end
