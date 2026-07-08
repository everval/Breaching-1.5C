cd(@__DIR__)
using Pkg
Pkg.activate(pwd())
Pkg.instantiate()

include("AdditionalFunctions.jl")

total_time = @elapsed begin

# Regression model
println("Running Analysis-ERF.jl..."); @time include("Analysis-ERF.jl")
println("Running Projections-ERF.jl..."); @time include("Projections-ERF.jl")

# Time series model
println("Running Analysis-ONI.jl..."); @time include("Analysis-ONI.jl")
println("Running Analysis-TSModel.jl..."); @time include("Analysis-TSModel.jl")
println("Running Coverage-TSModel.jl..."); @time include("Coverage-TSModel.jl")
println("Running Plot-Temperature-ONI.jl..."); @time include("Plot-Temperature-ONI.jl")
println("Running ProbabilityPaths-PreParis.jl..."); @time include("ProbabilityPaths-PreParis.jl")
println("Running ProbabilityPaths.jl..."); @time include("ProbabilityPaths.jl")
println("Running ProbabilityPathsPlots.jl..."); @time include("ProbabilityPathsPlots.jl")

end

println("Total runtime: ", total_time, " seconds")