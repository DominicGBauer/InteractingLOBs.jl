__precompile__()


module InteractingLOBs

using LinearAlgebra
using Statistics
using Distributions
using Random
using SharedArrays
using Distributed
using SpecialFunctions
using Logging
using IOLogging
using ProgressMeter
using TSSM

include("source_function.jl")
include("coupling_function.jl")
include("rl_push_function.jl")
include("randomness_function.jl")
include("data_passer.jl")
include("reaction_diffusion_path.jl")
include("parse_params.jl")
include("reaction_diffusion_spde.jl")
include("objective_surface.jl")
include("StylizedFacts.jl") # for now, this actually means that InteractingLOBs needs to declare that it depends on all the same things
                            # StylizedFacts depends on (i.e. anything that StylizedFacts.jl calls with "using"). Should one day
                            # be made into its own package
include("plot_price_impact.jl")
__version__ = "Interacting LOB"

foo(x) = 2*x

export  SLOB,
        SourceTerm,
        CouplingTerm,
        RLPushTerm,
        RandomnessTerm,
        DataPasser,
        parse_commandline,
        ObjectiveSurface,
        InteractOrderBooks,
        StylizedFacts,
        to_simulation_time,
        to_real_time,
        clear_double_dict,
        plot_price_impact,
        foo,
        SibuyaKernelModified

end # module
