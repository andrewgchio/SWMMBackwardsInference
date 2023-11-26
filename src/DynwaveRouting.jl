# DynwaveRouting.jl
#
# Contains routing functions

################################################################################
# Differentiable SWMM Implementation
################################################################################

function differentiable_swmm(G, T, Δt, src, Q_dwf::TR...) where {TR <: Real}
    
    ############################################################################
    # Initialize relevant dataframes for simulation
    ############################################################################

    l(x) = label_for(G,x)
    c(x) = code_for(G,x)
    
    # Inflow rate Q (cfs)
    Q = Dict() # DataFrame(time=T)
    
    # Hydraulic head H (ft)
    H = Dict() # DataFrame(time=T)
    
    # Initialize for node
    for v in gen_SWnodes(G)
        # IPOPT EDIT: Changed to start at v.invertElev, so that max(h1,z1) does
        # not need to occur later on.
        for (ti,t) in enumerate(T); H[ti,v.name] = v.invertElev end
        # H[1,v.name]  = v.initDepth + v.invertElev
    end
    
    # Initialize for edge
    for e in gen_SWedges(G)
        for (ti,t) in enumerate(T); Q[ti,e.name] = Real(0.0) end
        # Q[1,e.name]  = e.initFlow 
    end
    
    ############################################################################
    # Run simulation 
    ############################################################################
    
    for e in gen_SWedges(G); e.a1 = 0.0 end
    for v in gen_SWnodes(G); v.Qin, v.Qout = 0.0, 0.0 end
    
    for (ti,t) in enumerate(T)
        
        # skip last row
        if ti >= length(T); continue end
        
        # Only start printing debug statements right before flow starts
        
        ######################################################################## 
        # Get simulation indices
        ######################################################################## 
        
        # The next time step that we are aiming to compute
        # IPOPT EDIT TRY: Is using findfirst an issue?
        tj = ti+1 # findfirst(x -> x == t+Δt, T)
        
        ######################################################################## 
        # Set inflow for src according to dwf
        ######################################################################## 
        
        for v in gen_SWnodes(G);  v.newLatFlow = 0.0 end
        G[src].newLatFlow = Q_dwf[tj] # my_swmm_dwf_std(t+Δt, q=q, st=st, et=et)
        
        ######################################################################## 
        # Execute dynamic wave analysis
        ######################################################################## 
        
        ######################################################################## 
        # DW Step 1: 
        # Initially let Qlast and Hlast be the flow in each link and the head at 
        # each node, respectively, computed at time t. At time 0 these values 
        # are provided by the user-supplied initial conditions.
        ######################################################################## 
        
        for e in gen_SWedges(G)
            e.Qlast, e.Qnew = Q[ti,e.name], 0.0
            e.a2 = e.a1 # conduit area from solution at last time step
        end

        for v in gen_SWnodes(G)
            v.Hlast, v.Hnew = H[ti,v.name], 0.0
            v.oldNetInflow = v.Qin - v.Qout
        end
        
        # keep iterating until convergence
        # IPOPT EDIT: Continue iterating until MaxTrials is reached
        for Steps in 0:MAXTRIALS-1

            #################################################################### 
            # DW Step 2: 
            # Solve Eqn. 3-14 for each link producing a new flow estimate Qnew 
            # for time t + Δt, basing the values of A, Ā, Ū, and R̄ on Hlast.
            #################################################################### 
            
            # init nodes with inflow / outflow
            for v in gen_SWnodes(G)
                v.Qin, v.Qout = 0.0, 0.0
                # IPOPT EDIT: Instead of modifying Qin and Qout depending on the
                # value of newLatFlow at each node, we assume that the only new
                # lateral flow comes at the node with the dry weather flow.
                # Since this new flow > 0 , we only take the first branch.
                # if v.newLatFlow >= 0.0;  v.Qin  += v.newLatFlow
                # else                     v.Qout -= v.newLatFlow
                # end
            end
            G[src].Qin = G[src].newLatFlow

            # find link flows
            for e in gen_SWedges(G)
                
                # get most current heads at upstream/downstream ends of conduit
                n1, n2 = G[e.u], G[e.v]
                z1, z2 = n1.invertElev + e.offset1, n2.invertElev + e.offset2
                h1, h2 = n1.Hlast, n2.Hlast # n1.newDepth + n1.invertElev, n2.newDepth + n2.invertElev
                h1, h2 = max(h1, z1), max(h2, z2)
                
                # get unadjusted upstream/downstream flow depths in conduit
                #   (flow depth) = head in conduit - elev. of conduit invert
                # IPOPT ISSUE: Ignore max; This causes a domain error with a complex number
                # y1, y2 = h1-z1, h2-z2
                y1, y2 = max(h1-z1, FUDGE), max(h2-z2, FUDGE)

                # get area from solution at previous time step
                # IPOPT EDIT: Ignore max
                aOld = e.a2
                # aOld = max(e.a2, FUDGE)

                # use Courant-modified length instead of conduit's actual length
                length = e.length # e.modLength

                # find surface area contributions to upstream/downstream nodes
                # h1, h2, y1, y2 = findSurfArea(e, qLast, length, h1, h2, y1, y2)

                # compute area at each end of conduit & hyd. radius at upstream end
                # wSlot = _get_slot_width(e, y1) # Assume 0 since not using slot
                a1    = _get_area(e, y1)
                r1    = _get_hyd_rad(e, y1)
                a2    = _get_area(e, y2)
                # a1 = y1 * e.base # assume rectangular pipe
                # r1 = (e.base * y1) / (e.base + 2 * y1)
                # a2 = y2 * e.base

                # compute area & hyd. radius at midpoint
                Ȳ = 0.5 * (y1 + y2)
                Ā = _get_area(e, Ȳ) # Ȳ * e.base
                R̄ = _get_hyd_rad(e, Ȳ) # (e.base * Ȳ) / (e.base + 2 * Ȳ)

                # if edge is dry, then just record that
                # IPOPT EDIT: Ignore whether the node is dry or not. This
                # marginally increases the average difference
                # if Ā <= FUDGE
                    
                #     e.a1 = 0.5 * (a1 + a2)
                #     e.Qnew = 0.0
                    
                #     # TODO: head will change even if pipe is dry? (Ā ≤ FUDGE)
                    
                #     continue
                # end

                # compute velocity from last flow estimate
                # Ā = max(Ā, FUDGE)
                Ū = e.Qlast / Ā
                # Ū = min(MAXVELOCITY, Ū)
                # if abs(Ū) > MAXVELOCITY; Ū = MAXVELOCITY * sign(e.Qlast) end

                # compute Froude No. - link_getFroude(e, Ū, Ȳ)
                froude = abs(Ū) / sqrt(GRAVITY * (Ā / _get_width(e, Ȳ))) 
                # IPOPT EDIT: Igmore abs
                # froude = Ū / sqrt(GRAVITY * Ȳ) # link_getFroude(e, Ū, Ȳ)

                # find inertial damping factor (σ)
                # if     froude <= 0.5; σ = 1.0
                # elseif froude >= 1.0; σ = 0.0
                # else                  σ = 2.0 * (1.0 - froude)
                # end
                # IPOPT EDIT: Approximation; MathConstants.e is not a Real
                σ = 1.0/(1 + 2.718281828459045^(10 * (froude - 0.75))) 
                
                # get upstream-weighted area & hyd. rad. based on damping factor
                # ρ = 1.0
                # if e.Qlast > 0.0 && h1 >= h2; ρ = σ end
                # IPOPT EDIT: Use this approximation 
                ρ = σ
                Awtd = a1 + (Ā - a1) * ρ
                Rwtd = r1 + (R̄ - r1) * ρ

                # compute terms of momentum eqn.
                # 1. friction slope term
                ΔQ_fric = Δt * GRAVITY * 
                            (e.η / 1.486)^2 / (Rwtd ^ (1.33333)) * abs(Ū)

                # 2. energy slope term
                ΔQ_pres = -Δt * GRAVITY * Awtd * (h2 - h1) / length

                # 3 & 4. inertial terms, weighted with σ
                # dq3, dq4 = 0.0, 0.0
                # IPOPT EDIT: Ignore the if statement and always set dq3 and dq4
                # to be the values found below
                # if σ > 0.0
                #     dq3 = 2.0 * Ū * (Ā - aOld) * σ
                #     dq4 = Δt * Ū * Ū * (a2 - a1) / length * σ
                # end
                dq3 = 2.0 * Ū * (Ā - aOld) * σ
                dq4 = Δt * Ū * Ū * (a2 - a1) / length * σ
                ΔQ_iner = dq3 + dq4

                # Eqn 3-14
                e.Qnew = (Q[ti,e.name] + ΔQ_iner + ΔQ_pres) / (1.0 + ΔQ_fric)
                
                ################################################################ 
                # DW Step 3: 
                # Combine Qnew and Qlast together using a relaxation factor θ 
                # to produce a weighted value of Qnew
                ################################################################ 
                
                β = 1.486 / e.η * 
                    sqrt(abs(z1-z2) / sqrt(e.length^2 - abs(z1-z2)^2))
                Qnorm = β * a1 * r1 ^ (2.0/3.0)
                e.Qnew = min(Qnorm, e.Qnew)
                
                # IPOPT UNDO: remove reliance on steps
                if Steps > 0
                    e.Qnew = (1.0-Ω) * e.Qlast + Ω * e.Qnew
                    # IPOPT EDIT: Assume Qnew and Qlast will not change signs
                    # if e.Qnew*e.Qlast < 0.0; e.Qnew = 0.001*sign(e.Qnew) end
                end
                
                # if e.Qnew > FUDGE && n1.Hlast - n1.invertElev <= FUDGE; e.Qnew = FUDGE end
                # if e.Qnew < -FUDGE && n2.Hlast - n1.invertElev <= FUDGE; e.Qnew = -FUDGE end
                
                # TODO: The new depth should be set here and used in the
                # computation for the new head
                
                # Save new values
                e.a1 = Ā
                
            end
            
            # update node flows
            # for e in E
            #     n1,n2 = G[e.u], G[e.v]
            #     if e.Qnew >= 0.0
            #         n1.Qout += e.Qnew
            #         n2.Qin  += e.Qnew
            #     else
            #         n1.Qin  -= e.Qnew
            #         n2.Qout -= e.Qnew
            #     end
            # end
            
            # IPOPT EDIT: Assume that all Qnew values are positive for now, and 
            # sum them as shown. 
            for v in gen_SWnodes(G)
                # Update Qin
                for u in inneighbors(G,c(v.name))
                    v.Qin += G[l(u),v.name].Qnew
                end
            
                # Update Qout
                for w in outneighbors(G, c(v.name))
                    v.Qout += G[v.name,l(w)].Qnew
                end
            end
            
            ####################################################################
            # DW Step 4: 
            # Compute a value for Hnew at each node from Eqn 3-15 using the 
            # flows Qnew for Q^(t+Δt) and heads Hlast to evaluate A_S^(t+Δt).
            ####################################################################

            # Compute surfArea from contributing nodes
            surfArea = Dict(vname => A_Smin for vname in gen_SWvnames(G))
            for e in gen_SWedges(G)
                # IPOPT EDIT: Ignore dry vs not dry. This makes the difference 
                # jump up a bit.
                # if e.Qnew > FUDGE # Not a dry node
                # widths assume rectangular pipes
                width1   = _get_width(e, 0) 
                width2   = _get_width(e, 0) 
                widthMid = _get_width(e, 0) 
                surfArea[e.u] += (width1 + widthMid) * e.length / 4.0
                surfArea[e.v] += (widthMid + width2) * e.length / 4.0
                # end
            end
            
            for v in gen_SWnodes(G)
                
                # Skip outfalls
                # if v.nodetype == "OUTFALL"; continue end
                
                # net flow at node = Qin - Qout
                ΣQold, ΣQnew = 0.0, 0.0
                ΣQold = v.oldNetInflow
                ΣQnew = v.Qin - v.Qout
                
                A_SN   = 0 # no area from storage nodes
                ΣA_SL = surfArea[v.name]
                
                # IPOPT EDIT: Changed surfArea to start at A_Smin as a default
                A_S = A_SN + ΣA_SL
                # A_S = max(A_SN + ΣA_SL, A_Smin)
                Δh = (0.5 * Δt) * (ΣQold + ΣQnew) / A_S
                v.Hnew = H[ti,v.name] + Δh
                
                ################################################################ 
                # DW Step 5: 
                # As with flows, apply a relaxation factor to combine Hlast and 
                # Hnew
                ################################################################ 
                
                if Steps > 0; v.Hnew = (1.0-Ω) * v.Hlast + Ω * v.Hnew end
                
            end
            
            ####################################################################
            # DW Step 6: 
            # If Hnew is close enough to Hlast for each node then the process 
            # stops with Qnew and Hnew as the solution for time t+Δt. 
            # Otherwise, Hlast, Qlast are set equal to Hnew and Qnew, 
            # respectively, and the process returns to step 2.
            ####################################################################
            
            for v in gen_SWnodes(G); v.Hlast = v.Hnew end
            for e in gen_SWedges(G); e.Qlast = e.Qnew end
            
        end
        
        ######################################################################## 
        # Save into dataframes after convergence
        ######################################################################## 
        
        # Save hydraulic head
        # TODO: Set head at outfall to a constant. SWMM sets it to the 
        # elevation of the critical or normal flow depth.
        for v in gen_SWnodes(G)
            if v.nodetype == "OUTFALL"
                H[tj, v.name] = H[ti, v.name] # Let the outfall head not change
            else
                H[tj,v.name] = v.Hnew  
            end
        end
        
        for e in gen_SWedges(G);  Q[tj,e.name] = e.Qnew end
    end
    
    # For Ipopt to print out status
    flush(stdout)
    
    return Q
    
end

function _get_area(e, y)
    return y * e.xsect.wMax # Assuming closed rectangular pipe
    # FUDGE = 0.0001
    # return max(xsect_getAofY(e.xsect, y), FUDGE)
end

function _get_hyd_rad(e, y)
    a = y * e.xsect.wMax
    return a / (e.xsect.wMax + 2.0*a/e.xsect.wMax)
    # return (e.base * Ȳ) / (e.base + 2 * Ȳ)
    # FUDGE = 0.0001
    # return max(xsect_getRofY(e.xsect, y), FUDGE)
end

function _get_width(e, y)
    return e.xsect.wMax # Assuming closed rectangular pipe
    # return xsect_getWofY(e.xsect, y)
end

################################################################################
# Visualization
################################################################################

function plot_edgeflows(G, T, Q; ϵ=0.025, ax=plt, 
                    title="", xlabel="", ylabel="", ylimits=[0.0,0.5], 
                    legend=true)
    
    for e in gen_SWenames(G)
        if e in names(Q) && any(Q[:,e] .> ϵ) # Used to avoid plotting 0 lines
            ax.plot(Q[:,"t"], Q[:,e], label=e) 
        end
    end

    ax.set_title(title) 
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # PyPlot.grid(grid) 
    ax.set_ylim(ylimits...)

    if legend
        ax.legend(bbox_to_anchor=[1.05,1],loc=2,borderaxespad=0)
    end

end
