# StormwaterNetwork.jl
#
# Contains related functions for working with Stormwater Nodes and Edges

################################################################################
# Stormwater Network Models
################################################################################

# Construct a graph from a SWMM input file
function construct_graph(finp::String)
    
    # create graph object
    G = MetaGraph(DiGraph(), 
                  label_type=String, 
                  vertex_data_type=SWNode, 
                  edge_data_type=SWEdge)

    # Hydrological Model Reader (to read coords)
    f = hymo.SWMMInpFile(finp)
    
    # Read coordinates
    coords_map = Dict(
        convert(String, r["Node"]) => 
            (convert(Float64, r["X_Coord"]), convert(Float64, r["Y_Coord"])) 
        for (_,r) in f.coordinates.reset_index().iterrows())
        
    shapes_map = Dict(
        convert(String, r["Link"]) => (convert(String, r["Shape"]), 
            [
                convert(Float64, r["Geom1"]),
                convert(Float64, r["Geom2"]),
                convert(Float64, r["Geom3"]),
                convert(Float64, r["Geom4"])
            ])
            for (_,r) in f.xsections.reset_index().iterrows()
    )
    xsect_map = Dict()
    for (e,xsect_params) in shapes_map
        xsect = TXsect()
        xsect_type = xsect_stringTypeToCode(xsect_params[1])
        xsect_setParams(xsect, xsect_type, xsect_params[2], 1.0)
        xsect_map[e] = xsect
    end
    
    # Add all junctions as nodes
    for (_,node) in f.junctions.reset_index().iterrows()
        name = convert(String, node["Name"])
        x,y = coords_map[name]
        G[name] = SWNode(name=name, nodetype="JUNCTION",
                         x=x, y=y,
                         invertElev=node["Invert_Elev"],
                         fullDepth=node["Max_Depth"],
                         initDepth=node["Init_Depth"],
                         surDepth=node["Surcharge_Depth"],
                         pondedArea=node["Ponded_Area"])
    end
    
    # Also add all outlets as nodes
    for (_,node) in f.outfalls.reset_index().iterrows()
        name = convert(String, node["Name"])
        x,y = coords_map[name]
        G[name] = SWNode(name=name, nodetype="OUTFALL",
                         x=x, y=y,
                         invertElev=node["Invert_Elev"])
    end
    
    # Add conduits as edges
    for (_,link) in f.conduits.reset_index().iterrows()
        u = convert(String, link["Inlet_Node"])
        v = convert(String, link["Outlet_Node"])

        if !haskey(G, u) || !haskey(G, v); continue end

        # Conduit properties
        name = convert(String, link["Name"])
        G[u,v] = SWEdge(u=u, v=v, name=u*"-"*v,
                        length=link["Length"],
                        η=link["Manning_N"],
                        offset1=link["Inlet_Offset"],
                        offset2=link["Outlet_Offset"],
                        initFlow=link["Init_Flow"],
                        qMax=link["Max_Flow"],
                        xsect=xsect_map[name]
        )
    end

    return G
    
end
    
# Used to declare as SWNode(name="A", x=1, y=2)
@kwdef mutable struct SWNode
    name::String
    nodetype::String # "JUNCTION" or "OUTFALL"
    x::Real
    y::Real
    
    invertElev::Real = 0.0
    fullDepth::Real = 0.0
    initDepth::Real = 0.0
    surDepth::Real = 0.0
    pondedArea::Real = 0.0
    
    # simulation
    converged::Bool = false
    dYdT::Real = 0.0
    oldSurfArea::Real = 0.0
    newSurfArea::Real = 0.0
    inflow::Real = 0.0
    outflow::Real = 0.0
    losses::Real = 0.0
    sumdqdh::Real = 0.0
    
    oldDepth::Real = 0.0
    newDepth::Real = 0.0
    crownElev::Real = 0.0
    overflow::Real = 0.0
    newVolume::Real = 0.0
    fullVolume::Real = 0.0
    
    # version based on manual
    Hlast::Real = 0.0
    Hnew::Real = 0.0
    newLatFlow::Real = 0.0
    Qin::Real = 0.0
    Qout::Real = 0.0
    oldNetInflow::Real = 0.0
end

@kwdef mutable struct SWEdge
    u::String
    v::String
    name::String
    
    length::Real
    η::Real
    
    offset1::Real
    offset2::Real
    
    initFlow::Real
    qMax::Real
    
    # Dynwave
    base::Real = 1.0
    
    # simulation
    bypassed::Bool = false
    surfArea1::Real = 0.0
    surfArea2::Real = 0.0
    a1::Real = 0.0
    a2::Real = 0.0
    
    newDepth::Real = 0.0
    newVolume::Real = 0.0
    newFlow::Real = 0.0
    oldFlow::Real = 0.0
    q1::Real = 0.0
    q2::Real = 0.0
    modLength::Real = 0.0
    dqdh::Real = 0.0
    froude::Real = 0.0
    
    # version based on manual
    Qlast::Real = 0.0
    Qnew::Real = 0.0
    
    flowClass::String = ""

    xsect
end

################################################################################
# Graph Queries
################################################################################

# Return the subgraph that is local to the provided node
function get_local_branch(v_start, G; max_hops_down=5, max_hops_up=5)
    # Utility functions
    l(v) = label_for(G,v)
    c(v) = code_for(G,v)
        
    # Find the node that is max_hops downstream of v_start
    for i in 1:max_hops_down
        out = outneighbors(G, c(v_start))
        if !isempty(out)
            v_start = l(outneighbors(G, c(v_start))[1])
        end
    end
        
    # Compute upstream branch within max_hops
    to_visit = [(v_start,0)]
    visited = Set{String}()
    downstream = Vector{Integer}()
        
    while !isempty(to_visit)
        vname, hops = popfirst!(to_visit)
        if vname in visited; continue end

        push!(visited, vname)
        push!(downstream, c(vname))

        if hops <= max_hops_up
            for w in inneighbors(G, c(vname))
                push!(to_visit, (l(w), hops+1))
            end
        end
    end
            
    sg,vmap = induced_subgraph(G, downstream)
            
    return sg
end

function get_upstream_nodes(v_start, G; τ_d=typemax(Int))
    c(x) = code_for(G, x)
    l(x) = label_for(G, x)
    v = G[v_start]
    Vup = Set([v_start])
    to_visit = [(0.0,v)]
    
    while !isempty(to_visit)
        d,v = pop!(to_visit)
        
        for cu in inneighbors(G, c(v.name))
            u = G[l(cu)]
            new_d = d + G[u.name,v.name].length
            if new_d <= τ_d && !(u.name in Vup)
                push!(to_visit, (new_d, u))
                push!(Vup, u.name)
            end
        end
    end
    
    return Vup
end

################################################################################
# Iterator Utils
################################################################################

# Return an iterable of SWNodes
function gen_SWnodes(G)
    return (G[label_for(G,v)] for v in vertices(G))
end

# Return an iterable of SWEdges
function gen_SWedges(G)
    return (G[label_for(G,e.src), label_for(G, e.dst)] for e in edges(G))
end

# Return an iterable of node names
function gen_SWvnames(G; sort=false)
    if sort; return sort(v.name for v in gen_SWnodes(G))
    else return (v.name for v in gen_SWnodes(G))
    end
end

# Return an iterable of edge names
function gen_SWenames(G; sort=false)
    if sort; return sort(e.name for e in gen_SWedges(G))
    else return (e.name for e in gen_SWedges(G))
    end
end

################################################################################
# Graph I/O and Visualization Utilities 
################################################################################

# Print the nodes/edges of the graph
function print_graph(G::MetaGraph)
    for v in gen_SWnodes(G); println(v) end
    for e in gen_SWedges(G); println(e) end
end

# Plot the graph
function draw_graph(G::MetaGraph; with_labels=true)
    # Layout specified by coordinates in metagraph nodes
    function StaticLayout()
        function layout(G)
            xs = [v.x for v in gen_SWnodes(G)]
            ys = [v.y for v in gen_SWnodes(G)]
            return Point.(zip(xs, ys))
        end
        return layout
    end
    
    if with_labels
        nlabels = collect(gen_SWvnames(G))
    else
        nlabels = nothing
    end
    
    fig, ax, p = graphplot(G, layout=StaticLayout(), nlabels=nlabels)
    display(fig)
    
end

function print_graph_branch(G::MetaGraph, v::String; max_hops=5)
    Gsub = get_local_branch(v, G; max_hops_down=max_hops)
    print_graph(Gsub)
end

function draw_graph_branch(G::MetaGraph, v::String; 
            max_hops=5, with_labels=true)
    Gsub = get_local_branch(v, G; max_hops_down=max_hops)
    draw_graph(Gsub; with_labels=with_labels)
end
