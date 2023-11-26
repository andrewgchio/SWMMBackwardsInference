# SWMMBackwardsInference

[![Build Status](https://github.com/andrewgchio/SWMMBackwardsInference.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/andrewgchio/SWMMBackwardsInference.jl/actions/workflows/CI.yml?query=branch%3Amain)

SWMMBackwardsInference.jl is a Julia/JuMP package that implements a physics-informed backwards inference model to determine potential origins of *dry weather flows (DWFs)* in stormwater networks using sensor observations. Our code applies several approximations to the dynamic wave analysis method provided in [EPA SWMM](https://www.epa.gov/water-research/storm-water-management-model-swmm) to satisfy the computational requirements of backwards inference. 

Two key components are provided by SWMMBackwardsInference.jl: 

- *DynwaveRouting*, which re-implements the dynamic wave analysis method for the propagation of flow in [EPA SWMM](https://www.epa.gov/water-research/storm-water-management-model-swmm), and adds physics approximations to obtain a differentiable model.

- *BackwardsInference*, which leverages auto-differentiation with the flow routing model above to produce a set of potential origins with their respective flow curves and likelihood of occurrence, given a set of sensor observations. 

## Getting Started

This section details the installation of SWMMBackwardsInference.jl, and a quick tutorial with a simple example. 

### Dependencies

The dependencies for SWMMBackwardsInference.jl are the following: 

**Julia** packages:

- [JuMP](https://github.com/jump-dev/JuMP.jl) (version >= v1.15.0)
- [Ipopt](https://github.com/jump-dev/Ipopt.jl) (version >= 1.4.2)
- [Graphs](https://github.com/JuliaGraphs/Graphs.jl) (version >= 1.8.0)
- [MetaGraphsNext](https://github.com/JuliaGraphs/MetaGraphsNext.jl) (version >= 0.6.0)
- [GraphMakie](https://github.com/MakieOrg/GraphMakie.jl) (version >= 0.5.6)
- [CairoMakie](https://github.com/JuliaPlots/CairoMakie.jl) (version >= 0.10.9)
- [GeometryBasics](https://github.com/JuliaGeometry/GeometryBasics.jl) (version >= 0.4.9)
- [PyCall](https://github.com/JuliaPy/PyCall.jl) (version >= 1.96.1)
- [DataFrame](https://github.com/JuliaData/DataFrames.jl/tree/main) (version >= 1.61)
- [CSV](https://github.com/JuliaData/CSV.jl) (version >= 0.10.11)
- [COIN-HSL](https://licences.stfc.ac.uk/product/coin-hsl) (version >= 2023.5.26+0, optional, but recommended)

    The backwards inference model provided in SWMMBackwardsInference.jl works best when using the COIN-HSL package, which includes faster, non-linear solvers (as compared to the default solvers in Ipopt). We recommend installing and using these solvers if possible, but note that the default solver will be used automatically if none of these solvers are found.

**Python** packages:

- [hymo](https://github.com/lucashtnguyen/hymo) (version >= 0.1.1b)
- [pyswmm](https://github.com/pyswmm/pyswmm) (version >= 0.6.2, optional)

### Installation

The latest release of SWMMBackwardsInference.jl can by cloning the project. 

```
git clone https://github.com/andrewgchio/SWMMBackwardsInference
```

> ***Note:*** It is highly recommended to install and use COIN-HSL (listed above) to enable faster backwards inference. Instructions for installation depend on the machine and architecture. More details can be found [here](https://github.com/coin-or-tools/ThirdParty-HSL).

### Quickstart Tutorial

#### Initial Setup

To start working with SWMMBackwardsInference.jl, an EPA SWMM input file is required. This file can be created using the [EPA SWMM](https://www.epa.gov/water-research/storm-water-management-model-swmm) GUI. 

For this tutorial, we use the sample network in `examples/SampleNetwork.inp`. This small network consists of 4 junctions connected by 4 conduits. We assume that a __dry weather flow__ anomaly is injected at node `J1`, which then propagates to nodes `J2` and `J4` before reaching the outfall `Out1`. The flow anomaly has a peak flow of `Q=0.2` cfs at time `t=80`. The network and anomaly injection profile can be visualized as follows:

<network and anomaly>

To use SWMMBackwardsInference.jl, the package must first be brought into the current namespace. We first construct a graph of the stormwater network using the EPA SWMM input file, as well as the time range/time step over which simulation and inference is run, as shown below:

```
using SWMMBackwardsInference

# initialize sample network
G = construct_graph("./examples/tutorial/tutorial.inp") 

dt = 30           # 30 second time step
T  = 0:dt:2*60*60 # 2 hours, in seconds
```

#### Running a Forwards Simulation

Running the forwards simulation requires a source junction node `src` and a time-series inflow `Qdwf`. These represent the location, and amount of inflow introduced as part of the anomaly, respectively. 

Note that `Qdwf` is only constrained to fall within `Qmin` and `Qmax`, which are set in the configuration file `config.txt`. 

```
src = "J1"

# Initialize Qdwf 
Qdwf = ... 
```

For demonstration purposes, we model `Qdwf` to follow the normal distribution curve shown above; this can be done as follows: 

```
import Distributions: Normal, pdf

Qdwf = zeros(length(T))
st,et,q = 60, 100, 0.2
μ,σ = st + 0.5 * (et-st), (et-st)/3
for (i,t) in enumerate(T)
    Qdwf[i] = pdf(Normal(μ,σ), t)
end

# Normalize Qdwf height to around q
Qdwf /= (maximum(Qdwf) / q)
```

Then, the differentiable forwards model can be run as follows:

```
Q = differentiable_swmm(G, T, dt, src, Qdwf...) 
```

Note that here, `Q` is a dictionary of edge flows at varying timestep *indices*. 

#### Running the Backwards Inference

The backwards inference model can be run either with a specific source in mind, or with all nodes as possible. 

First, a set of sensors, and their corresponding observations must be initialized to observe parts of the network. Since flows are primarily defined on edges (*not nodes*), initialize sensors `S` to be a list of edges in the network, e.g., 

```
S = ["J2-J4"]
```

The corresponding observations made by sensors should be in a `DataFrame` with `length(T)` rows and columns with names of edges. `SWMMBackwardsInference.jl` provides a utility function to turn the dictionary returned by `differentiable_swmm`into a `DataFrame` as follows: 

```
Qobs = dict2df(Q, T, S)[:, ["time", S...]]
```

To run it with a single source, use: 

```
Qinf, objval, status, pstatus = single_source_inference(G, src, S, Qobs, T, dt)
```

The function `single_source_inference` returns 4 values: 
- the inferred array of inflows `Qinf`
- the mean squared error between observations and simulated ground truth `objval`
- the termination status of the inference model `status`
- the primal termination status of the inference model `pstatus`

Or, to iterate through multiple different sources, use: 

```
Qinfs = all_source_inference(G, S, Qobs, T, dt)
```

Note that this returns a dictionary of `src => (objval, Qinf)`. 

The full code can be seen at `main/example.jl`. 

## License

```
The MIT License (MIT)

Copyright (c) 2023 Andrew Chio

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
