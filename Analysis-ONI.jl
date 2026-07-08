cd(@__DIR__)
using Pkg
Pkg.activate(pwd())
using CSV, MarSwitching, JLD2, DataFrames, Dates, Random
include("AdditionalFunctions.jl")

Random.seed!(1234)

all_data = CSV.read("data/Compiled_Global_Temperature_Data.csv", DataFrame)

nino = select(all_data, :Date, :ONI_Anomaly)
rename!(nino, :ONI_Anomaly => :N34 )
dropmissing!(nino)

nino_model3 = MSModel(nino[!, :N34], 3);

nino_model5 = MSModel(nino[!, :N34], 5);

nino_model7 = MSModel(nino[!, :N34], 7);

nino_model9 = MSModel(nino[!, :N34], 9);;

nino_model11 = MSModel(nino[!, :N34], 11);

nino_model13 = MSModel(nino[!, :N34], 13);

# Table 9: AIC for different number of regimes
msm_aic(model) = 2 * length(model.raw_params) - 2 * model.Likelihood

regime_models = [nino_model3, nino_model5, nino_model7, nino_model9, nino_model11, nino_model13];
n_regimes = [3, 5, 7, 9, 11, 13];

aic_table = DataFrame(Regimes = n_regimes, AIC = msm_aic.(regime_models))
CSV.write("results/Table9.csv", aic_table)

# Tables 10 and 11: Estimation results of the best model (7 states)
summary_msm(nino_model7)

# Save the best model (7 states) to a JLD2 file for later use
@save "results/nino_model7.jld2" nino_model7
