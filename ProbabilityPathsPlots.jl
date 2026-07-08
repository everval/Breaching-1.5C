cd(@__DIR__)
using Pkg
Pkg.activate(pwd())
using CSV, DataFrames, Dates, Random, StatsPlots, HypothesisTests, JLD2, MarSwitching
include("AdditionalFunctions.jl")


hadcrut = CSV.read("results/ProbabilityPaths-HadCRUT.csv", DataFrame)
gistemp = CSV.read("results/ProbabilityPaths-GISTEMP.csv", DataFrame)
berkeley = CSV.read("results/ProbabilityPaths-Berkeley.csv", DataFrame)

xls = (Date(2023, 10, 1), Date(2065, 03, 1))

# Figure 8: Probability of breaching the thresholds
plot(hadcrut."Date (month)", hadcrut."1.5°C Threshold", label="HadCRUT(1.5°C)", xlabel="Date (monthly)", ylabel="Probability of breaching", title="Probability of breaching the thresholds", legend=:inside, linewidth=3, linestyle= :dashdot, color=:orange)
plot!(gistemp."Date (month)", gistemp."1.5°C Threshold", label="GISTEMP(1.5°C)", linewidth=3, linestyle= :dash, color=:darkorange2)
plot!(berkeley."Date (month)", berkeley."1.5°C Threshold", label="Berkeley   (1.5°C)", linewidth=3, linestyle= :dot, color=:darkorange)

plot!(hadcrut."Date (month)", hadcrut."2°C Threshold", label="HadCRUT  (2°C)", linewidth=3, linestyle= :dashdot, color=:red2)
plot!(gistemp."Date (month)", gistemp."2°C Threshold", label="GISTEMP  (2°C)", linewidth=3, linestyle= :dash, color=:crimson)
plot!(berkeley."Date (month)", berkeley."2°C Threshold", label="Berkeley     (2°C)", linewidth=3, linestyle= :dot, color=:firebrick)

plot!(xlims=xls, xticks=(hadcrut."Date (month)"[223:90:end-100], Dates.format.(hadcrut."Date (month)"[223:90:end-100], "mm/yyyy")), ylims=(0, 1))

savefig("figures/Figure8.pdf")



hadcrut_pre_pa = CSV.read("results/ProbabilityPaths-HadCRUT-PrePA.csv", DataFrame)
gistemp_pre_pa = CSV.read("results/ProbabilityPaths-GISTEMP-PrePA.csv", DataFrame)
berkeley_pre_pa = CSV.read("results/ProbabilityPaths-Berkeley-PrePA.csv", DataFrame)

xls = (Date(2025, 10, 1), Date(2050, 07, 1))

# Figure 9: Probability of breaching the thresholds before the Paris Agreement
plot(hadcrut."Date (month)", hadcrut."1.5°C Threshold", label="HadCRUT", xlabel="Date (monthly)", ylabel="Probability of breaching", title="Probability of breaching the 1.5°C threshold", linewidth=3, linestyle= :dashdot, color=:orange)
plot!(gistemp."Date (month)", gistemp."1.5°C Threshold", label="GISTEMP", linewidth=3, linestyle= :dash, color=:darkorange2)
plot!(berkeley."Date (month)", berkeley."1.5°C Threshold", label="Berkeley", linewidth=3, linestyle= :dot, color=:darkorange)

plot!(hadcrut_pre_pa."Date (month)", hadcrut_pre_pa."1.5°C Threshold", label="HadCRUT - PA-Start", xlabel="Date (monthly)", ylabel="Probability of breaching", title="Probability of breaching the 1.5°C threshold", linewidth=3, linestyle= :dashdot, color=:pink2, xlims=xls)
plot!(gistemp_pre_pa."Date (month)", gistemp_pre_pa."1.5°C Threshold", label="GISTEMP - PA-Start", linewidth=3, linestyle= :dash, color=:deeppink3)
plot!(berkeley_pre_pa."Date (month)", berkeley_pre_pa."1.5°C Threshold", label="Berkeley    - PA-Start", linewidth=3, linestyle= :dot, color=:hotpink)

plot!(xlims=xls, xticks=(hadcrut."Date (month)"[223:60:end-100], Dates.format.(hadcrut."Date (month)"[223:60:end-100], "mm/yyyy")), ylims=(0, 1), legend = :bottomright)

savefig("figures/Figure9.pdf")
