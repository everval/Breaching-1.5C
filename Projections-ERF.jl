cd(@__DIR__)
using Pkg
Pkg.activate(pwd())
using StatsPlots, Dates, CSV, DataFrames, HypothesisTests
include("AdditionalFunctions.jl")


rawtemp = CSV.read("data/HadCRUT 5.0 Global Annual Summary.csv", DataFrame)
temp = rawtemp[!, [:Time, :Anomaly]]
rename!(temp, :Time => :Year, :Anomaly => :Temp)

oldbase = mean(temp[(temp.Year.>=1850).&(temp.Year.<1900), :Temp])
temp[!, :Temp] = temp[!, :Temp] .- oldbase;

rawerf = CSV.read("data/ERF_best_aggregates_1750-2024.csv", DataFrame)
rawerf_05 = CSV.read("data/ERF_p05_aggregates_1750-2024.csv", DataFrame)
rawerf_95 = CSV.read("data/ERF_p95_aggregates_1750-2024.csv", DataFrame)
erf = rawerf[!, [:timebound_lower, :total]]
rename!(erf, :timebound_lower => :Year, :total => :ERF)
erf.ERF_p05 = rawerf_05[!, :total];
erf.ERF_p95 = rawerf_95[!, :total];

calstart = 1980
projstart = 2015
predstart = 2025
predend = 2100;

temperf_all = leftjoin(temp, erf, on = :Year)
dropmissing!(temperf_all)

temperf = temperf_all[temperf_all.Year .>= calstart, :]

reg = robust_est(temperf.Temp, [ones(length(temperf.ERF)) temperf.ERF]; verbose=true);

erfproj = DataFrame();
erfpref = DataFrame();

ssp119_best = CSV.read("data/SSPs/table_A3.4a_ssp119_ERF_1750-2500_best_estimate.csv", DataFrame);
ssp119_p05 = CSV.read("data/SSPs/table_A3.4a_ssp119_ERF_1750-2500_5th_percentile.csv", DataFrame);
ssp119_p95 = CSV.read("data/SSPs/table_A3.4a_ssp119_ERF_1750-2500_95th_percentile.csv", DataFrame);

erfproj.Year = ssp119_best[!, :year];
erfproj.SSP119_Best = ssp119_best[!, :total];
erfproj.SSP119_p05 = ssp119_p05[!, :total];
erfproj.SSP119_p95 = ssp119_p95[!, :total];

ssp126_best = CSV.read("data/SSPs/table_A3.4b_ssp126_ERF_1750-2500_best_estimate.csv", DataFrame);
ssp126_p05 = CSV.read("data/SSPs/table_A3.4b_ssp126_ERF_1750-2500_5th_percentile.csv", DataFrame);
ssp126_p95 = CSV.read("data/SSPs/table_A3.4b_ssp126_ERF_1750-2500_95th_percentile.csv", DataFrame);

erfproj.SSP126_Best = ssp126_best[!, :total];
erfproj.SSP126_p05 = ssp126_p05[!, :total];
erfproj.SSP126_p95 = ssp126_p95[!, :total];

ssp245_best = CSV.read("data/SSPs/table_A3.4c_ssp245_ERF_1750-2500_best_estimate.csv", DataFrame);
ssp245_p05 = CSV.read("data/SSPs/table_A3.4c_ssp245_ERF_1750-2500_5th_percentile.csv", DataFrame);
ssp245_p95 = CSV.read("data/SSPs/table_A3.4c_ssp245_ERF_1750-2500_95th_percentile.csv", DataFrame);

erfproj.SSP245_Best = ssp245_best[!, :total];
erfproj.SSP245_p05 = ssp245_p05[!, :total];
erfproj.SSP245_p95 = ssp245_p95[!, :total];

ssp370_best = CSV.read("data/SSPs/table_A3.4d_ssp370_ERF_1750-2500_best_estimate.csv", DataFrame);
ssp370_p05 = CSV.read("data/SSPs/table_A3.4d_ssp370_ERF_1750-2500_5th_percentile.csv", DataFrame);
ssp370_p95 = CSV.read("data/SSPs/table_A3.4d_ssp370_ERF_1750-2500_95th_percentile.csv", DataFrame);

erfproj.SSP370_Best = ssp370_best[!, :total];
erfproj.SSP370_p05 = ssp370_p05[!, :total];
erfproj.SSP370_p95 = ssp370_p95[!, :total];

ssp585_best = CSV.read("data/SSPs/table_A3.4e_ssp585_ERF_1750-2500_best_estimate.csv", DataFrame);
ssp585_p05 = CSV.read("data/SSPs/table_A3.4e_ssp585_ERF_1750-2500_5th_percentile.csv", DataFrame);
ssp585_p95 = CSV.read("data/SSPs/table_A3.4e_ssp585_ERF_1750-2500_95th_percentile.csv", DataFrame);

erfproj.SSP585_Best = ssp585_best[!, :total];
erfproj.SSP585_p05 = ssp585_p05[!, :total];
erfproj.SSP585_p95 = ssp585_p95[!, :total];

erfproj = erfproj[erfproj.Year .>= projstart .&& erfproj.Year .<= predend, :]
erfpred = erfproj[erfproj.Year .>= predstart .&& erfproj.Year .<= predend, :]

## Figure 3(a): Observed ERF and projections under SSPs
pj = plot(temperf.Year, temperf.ERF, label="Observed ERF", color=:black, linewidth=2, title="Effective Radiative Forcing Observed and Projections", xlabel="Year", ylabel="ERF (W/m2)", legend=:topleft)
plot!(pj, erfpred.Year, erfpred.SSP119_Best, label="SSP1-1.9", xlabel="Year", ylabel="ERF (W/m2)", legend=:topleft, linewidth=2, linestyle=:dash, color=:darkorange)
plot!(pj, erfpred.Year, erfpred.SSP126_Best, label="SSP1-2.6", linewidth=2, linestyle=:dash, color=:green)
plot!(pj, erfpred.Year, erfpred.SSP245_Best, label="SSP2-4.5", linewidth=2, linestyle=:dash, color=:red)
plot!(pj, erfpred.Year, erfpred.SSP370_Best, label="SSP3-7.0", linewidth=2, linestyle=:dash, color=:purple)
plot!(pj, erfpred.Year, erfpred.SSP585_Best, label="SSP5-8.5", linewidth=2, xlimit=(1980, 2100), linestyle=:dash, color=:brown)

savefig("figures/Figure3(a).pdf")


dist_ss119_p05 = -erfpred.SSP119_p05 .+ erfpred.SSP119_Best
dist_ss126_p05 = -erfpred.SSP126_p05 .+ erfpred.SSP126_Best
dist_ss245_p05 = -erfpred.SSP245_p05 .+ erfpred.SSP245_Best
dist_ss370_p05 = -erfpred.SSP370_p05 .+ erfpred.SSP370_Best
dist_ss585_p05 = -erfpred.SSP585_p05 .+ erfpred.SSP585_Best

dist_ss119_p95 = erfpred.SSP119_p95 .- erfpred.SSP119_Best
dist_ss126_p95 = erfpred.SSP126_p95 .- erfpred.SSP126_Best
dist_ss245_p95 = erfpred.SSP245_p95 .- erfpred.SSP245_Best
dist_ss370_p95 = erfpred.SSP370_p95 .- erfpred.SSP370_Best
dist_ss585_p95 = erfpred.SSP585_p95 .- erfpred.SSP585_Best

σ_ssp119 = (dist_ss119_p05 .+ dist_ss119_p95) ./ (2 * 1.6448536269)
σ_ssp126 = (dist_ss126_p05 .+ dist_ss126_p95) ./ (2 * 1.6448536269)
σ_ssp245 = (dist_ss245_p05 .+ dist_ss245_p95) ./ (2 * 1.6448536269)
σ_ssp370 = (dist_ss370_p05 .+ dist_ss370_p95) ./ (2 * 1.6448536269)
σ_ssp585 = (dist_ss585_p05 .+ dist_ss585_p95) ./ (2 * 1.6448536269)

function predict_temperature_anomalies(scenario, regression, σ=0.0)
    T = length(scenario)

    if σ != 0.0
        erf_errors = rand.(Normal.(0, σ))
    else
        erf_errors = zeros(T)
    end

    # Sample the coefficients from a normal distribution with the estimated coefficients and standard deviation
    sigma = regression.betavar
    sigma_pd = Symmetric(sigma + 1e-8 * I)  # small nugget for robustness
    mvn_dist = MvNormal(regression.β, sigma_pd)
    betas = rand(mvn_dist)

    updated_scenario = scenario .+ erf_errors

    X = [ones(T) updated_scenario]
    Yfit = X * betas
    Yerr = Yfit .+ rand(Normal(0, sqrt(regression.σ²)), T)
    return (Yerr = Yerr, Yfit = Yfit)
end

npaths = 10^5
fullcalendar = calstart:predend
ncal = length(fullcalendar)
n_test = ncal - size(temperf, 1)
dates_forecast = collect(predstart:predend)

# Figure 3(b): Forecasting temperature anomalies under SSP1-1.9 with prediction intervals
matforecasts_ssp119 = zeros(n_test, npaths)

for ii = 1:npaths
    matforecasts_ssp119[:, ii] = predict_temperature_anomalies(erfpred.SSP119_Best, reg, σ_ssp119).Yerr
end

quantilesforecasts_ssp119 = zeros(n_test, 9)

for ii = 1:n_test
    quantilesforecasts_ssp119[ii, :] = quantile(matforecasts_ssp119[ii, :], [0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995 ])
end

p_pi = plot(temperf.Year, temperf.Temp, label="Observed", xlabel="Date", ylabel="Temperature Anomaly (°C)", title="Forecasted Temperature with PI (SSP1-1.9)")

plot!(dates_forecast, quantilesforecasts_ssp119[1:end, 5], label="Median", color=:red, linewidth=2, linestyle=:dot)
plot!(dates_forecast, quantilesforecasts_ssp119[1:end, 5], fillrange=(quantilesforecasts_ssp119[1:end, 4],quantilesforecasts_ssp119[1:end, 6]), fillalpha=0.3, label="50% PI", color=:purple, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts_ssp119[1:end, 5], fillrange=(quantilesforecasts_ssp119[1:end, 3],quantilesforecasts_ssp119[1:end, 7]), fillalpha=0.2, label="90% PI", color=:green, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts_ssp119[1:end, 5], fillrange=(quantilesforecasts_ssp119[1:end, 1],quantilesforecasts_ssp119[1:end, 9]), fillalpha=0.1, label="99% PI", color=:red, linezalpha=0, linealpha=0)
vline!([predstart], label="Predictions", linestyle=:dashdot, color=:black, legend=:topleft, legendfontsize=10)

savefig("figures/Figure3(b).pdf")


dummies1 = zeros(length(fullcalendar), npaths)

for ii = 1:npaths
    fullpath = [temperf.Temp; matforecasts_ssp119[1:end, ii]]
    for jj = 10:(ncal-10)
        dummies1[jj, ii] = mean(fullpath[(jj-9):(jj+10)])
    end
end

pa151 = dropdims(mean(dummies1[:, :] .>= 1.5, dims=2), dims=2);
pa201 = dropdims(mean(dummies1[:, :] .>= 2, dims=2), dims=2);

plot(fullcalendar, pa151, label="Probability of breaching 1.5°C", xlabel="Year", ylabel="Probability", title="Probability of Breaching 1.5°C and 2°C", legend=:bottomright, linewidth=4, color=:darkorange, linestyle=:dash)
plot!(fullcalendar, pa201, label="Probability of breaching 2°C", color=:red3, xlims=(2016, 2087), ylims=(0, 1), xticks=2016:10:2100, yticks=0:0.1:1, linewidth=4, linestyle=:dashdot)

probabilities_df = DataFrame(Year=fullcalendar, Probability_15_SSP1=pa151, Probability_20_SSP1=pa201)

matforecasts_ssp126 = zeros(n_test, npaths)

for ii = 1:npaths
    matforecasts_ssp126[:, ii] = predict_temperature_anomalies(erfpred.SSP126_Best, reg, σ_ssp126).Yerr
end

quantilesforecasts_ssp126 = zeros(n_test, 9)

for ii = 1:n_test
    quantilesforecasts_ssp126[ii, :] = quantile(matforecasts_ssp126[ii, :], [0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995 ])
end

# Figure 3(c): Forecasting temperature anomalies under SSP1-2.6 with prediction intervals
p_pi2 = plot(temperf.Year, temperf.Temp, label="Observed", xlabel="Date", ylabel="Temperature Anomaly (°C)", title="Forecasted Temperature with PI (SSP1-2.6)")
plot!(dates_forecast, quantilesforecasts_ssp126[1:end, 5], label="Median", color=:red, linewidth=2, linestyle=:dot)
plot!(dates_forecast, quantilesforecasts_ssp126[1:end, 5], fillrange=(quantilesforecasts_ssp126[1:end, 4],quantilesforecasts_ssp126[1:end, 6]), fillalpha=0.3, label="50% PI", color=:purple, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts_ssp126[1:end, 5], fillrange=(quantilesforecasts_ssp126[1:end, 3],quantilesforecasts_ssp126[1:end, 7]), fillalpha=0.2, label="90% PI", color=:green, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts_ssp126[1:end, 5], fillrange=(quantilesforecasts_ssp126[1:end, 1],quantilesforecasts_ssp126[1:end, 9]), fillalpha=0.1, label="99% PI", color=:red, linezalpha=0, linealpha=0)
vline!([predstart], label="Predictions", linestyle=:dashdot, color=:black, legend=:topleft, legendfontsize=10)

savefig("figures/Figure3(c).pdf")


dummies2 = zeros(length(fullcalendar), npaths)

for ii = 1:npaths
    fullpath = [temperf.Temp; matforecasts_ssp126[1:end, ii]]
    for jj = 10:(ncal-10)
        dummies2[jj, ii] = mean(fullpath[(jj-9):(jj+10)])
    end
end

pa152 = dropdims(mean(dummies2[:, :] .>= 1.5, dims=2), dims=2);
pa202 = dropdims(mean(dummies2[:, :] .>= 2, dims=2), dims=2);

plot(fullcalendar, pa152, label="Probability of breaching 1.5°C", xlabel="Year", ylabel="Probability", title="Probability of Breaching 1.5°C and 2°C", legend=:bottomright, linewidth=4, color=:darkorange, linestyle=:dash)
plot!(fullcalendar, pa202, label="Probability of breaching 2°C", color=:red3, xlims=(2016, 2087), ylims=(0, 1), xticks=2016:10:2100, yticks=0:0.1:1, linewidth=4, linestyle=:dashdot)

probabilities_df.Probability_15_SSP2 = pa152
probabilities_df.Probability_20_SSP2 = pa202

matforecasts_ssp245 = zeros(n_test, npaths)

for ii = 1:npaths
    matforecasts_ssp245[:, ii] = predict_temperature_anomalies(erfpred.SSP245_Best, reg, σ_ssp245).Yerr
end

quantilesforecasts_ssp245 = zeros(n_test, 9)

for ii = 1:n_test
    quantilesforecasts_ssp245[ii, :] = quantile(matforecasts_ssp245[ii, :], [0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995 ])
end

# Figure 3(d): Forecasting temperature anomalies under SSP2-4.5 with prediction intervals
p_pi3 = plot(temperf.Year, temperf.Temp, label="Observed", xlabel="Date", ylabel="Temperature Anomaly (°C)", title="Forecasted Temperature with PI (SSP2-4.5)")
plot!(dates_forecast, quantilesforecasts_ssp245[1:end, 5], label="Median", color=:red, linewidth=2, linestyle=:dot)
plot!(dates_forecast, quantilesforecasts_ssp245[1:end, 5], fillrange=(quantilesforecasts_ssp245[1:end, 4],quantilesforecasts_ssp245[1:end, 6]), fillalpha=0.3, label="50% PI", color=:purple, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts_ssp245[1:end, 5], fillrange=(quantilesforecasts_ssp245[1:end, 3],quantilesforecasts_ssp245[1:end, 7]), fillalpha=0.2, label="90% PI", color=:green, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts_ssp245[1:end, 5], fillrange=(quantilesforecasts_ssp245[1:end, 1],quantilesforecasts_ssp245[1:end, 9]), fillalpha=0.1, label="99% PI", color=:red, linezalpha=0, linealpha=0)
vline!([predstart], label="Predictions", linestyle=:dashdot, color=:black, legend=:topleft, legendfontsize=10)

savefig("figures/Figure3(d).pdf")

dummies3 = zeros(length(fullcalendar), npaths)

for ii = 1:npaths
    fullpath = [temperf.Temp; matforecasts_ssp245[1:end, ii]]
    for jj = 10:(ncal-10)
        dummies3[jj, ii] = mean(fullpath[(jj-9):(jj+10)])
    end
end

pa153 = dropdims(mean(dummies3[:, :] .>= 1.5, dims=2), dims=2);
pa203 = dropdims(mean(dummies3[:, :] .>= 2, dims=2), dims=2);

plot(fullcalendar, pa153, label="Probability of breaching 1.5°C", xlabel="Year", ylabel="Probability", title="Probability of Breaching 1.5°C and 2°C", legend=:bottomright, linewidth=4, color=:darkorange, linestyle=:dash)
plot!(fullcalendar, pa203, label="Probability of breaching 2°C", color=:red3, xlims=(2016, 2087), ylims=(0, 1), xticks=2016:10:2100, yticks=0:0.1:1, linewidth=4, linestyle=:dashdot)

probabilities_df.Probability_15_SSP3 = pa153
probabilities_df.Probability_20_SSP3 = pa203

matforecasts_ssp370 = zeros(n_test, npaths)

for ii = 1:npaths
    matforecasts_ssp370[:, ii] = predict_temperature_anomalies(erfpred.SSP370_Best, reg, σ_ssp370).Yerr
end

quantilesforecasts_ssp370 = zeros(n_test, 9)

for ii = 1:n_test
    quantilesforecasts_ssp370[ii, :] = quantile(matforecasts_ssp370[ii, :], [0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995 ])
end

# Figure 3(e): Forecasting temperature anomalies under SSP3-7.0 with prediction intervals
p_pi4 = plot(temperf.Year, temperf.Temp, label="Observed", xlabel="Date", ylabel="Temperature Anomaly (°C)", title="Forecasted Temperature with PI (SSP3-7.0)")

plot!(dates_forecast, quantilesforecasts_ssp370[1:end, 5], label="Median", color=:red, linewidth=2, linestyle=:dot)
plot!(dates_forecast, quantilesforecasts_ssp370[1:end, 5], fillrange=(quantilesforecasts_ssp370[1:end, 4],quantilesforecasts_ssp370[1:end, 6]), fillalpha=0.3, label="50% PI", color=:purple, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts_ssp370[1:end, 5], fillrange=(quantilesforecasts_ssp370[1:end, 3],quantilesforecasts_ssp370[1:end, 7]), fillalpha=0.2, label="90% PI", color=:green, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts_ssp370[1:end, 5], fillrange=(quantilesforecasts_ssp370[1:end, 1],quantilesforecasts_ssp370[1:end, 9]), fillalpha=0.1, label="99% PI", color=:red, linezalpha=0, linealpha=0)
vline!([predstart], label="Predictions", linestyle=:dashdot, color=:black, legend=:topleft, legendfontsize=10)

display(p_pi4)

savefig("figures/Figure3(e).pdf")


dummies4 = zeros(length(fullcalendar), npaths)

for ii = 1:npaths
    fullpath = [temperf.Temp; matforecasts_ssp370[1:end, ii]]
    for jj = 10:(ncal-10)
        dummies4[jj, ii] = mean(fullpath[(jj-9):(jj+10)])
    end
end

pa154 = dropdims(mean(dummies4[:, :] .>= 1.5, dims=2), dims=2);
pa204 = dropdims(mean(dummies4[:, :] .>= 2, dims=2), dims=2);

plot(fullcalendar, pa154, label="Probability of breaching 1.5°C", xlabel="Year", ylabel="Probability", title="Probability of Breaching 1.5°C and 2°C", legend=:bottomright, linewidth=4, color=:darkorange, linestyle=:dash)
plot!(fullcalendar, pa204, label="Probability of breaching 2°C", color=:red3, xlims=(2016, 2087), ylims=(0, 1), xticks=2016:10:2100, yticks=0:0.1:1, linewidth=4, linestyle=:dashdot)

probabilities_df.Probability_15_SSP4 = pa154
probabilities_df.Probability_20_SSP4 = pa204

matforecasts_ssp585 = zeros(n_test, npaths)

for ii = 1:npaths
    matforecasts_ssp585[:, ii] = predict_temperature_anomalies(erfpred.SSP585_Best, reg, σ_ssp585).Yerr
end

quantilesforecasts_ssp585 = zeros(n_test, 9)

for ii = 1:n_test
    quantilesforecasts_ssp585[ii, :] = quantile(matforecasts_ssp585[ii, :], [0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995 ])
end

# Figure 3(f): Forecasting temperature anomalies under SSP5-8.5 with prediction intervals
p_pi5 = plot(temperf.Year, temperf.Temp, label="Observed", xlabel="Date", ylabel="Temperature Anomaly (°C)", title="Forecasted Temperature with PI (SSP5-8.5)")

plot!(dates_forecast, quantilesforecasts_ssp585[1:end, 5], label="Median", color=:red, linewidth=2, linestyle=:dot)
plot!(dates_forecast, quantilesforecasts_ssp585[1:end, 5], fillrange=(quantilesforecasts_ssp585[1:end, 4],quantilesforecasts_ssp585[1:end, 6]), fillalpha=0.3, label="50% PI", color=:purple, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts_ssp585[1:end, 5], fillrange=(quantilesforecasts_ssp585[1:end, 3],quantilesforecasts_ssp585[1:end, 7]), fillalpha=0.2, label="90% PI", color=:green, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts_ssp585[1:end, 5], fillrange=(quantilesforecasts_ssp585[1:end, 1],quantilesforecasts_ssp585[1:end, 9]), fillalpha=0.1, label="99% PI", color=:red, linezalpha=0, linealpha=0)
vline!([predstart], label="Predictions", linestyle=:dashdot, color=:black, legend=:topleft, legendfontsize=10)

savefig("figures/Figure3(f).pdf")

dummies5 = zeros(length(fullcalendar), npaths)

for ii = 1:npaths
    fullpath = [temperf.Temp; matforecasts_ssp585[1:end, ii]]
    for jj = 10:(ncal-10)
        dummies5[jj, ii] = mean(fullpath[(jj-9):(jj+10)])
    end
end

pa155 = dropdims(mean(dummies5[:, :] .>= 1.5, dims=2), dims=2);
pa205 = dropdims(mean(dummies5[:, :] .>= 2, dims=2), dims=2);

plot(fullcalendar, pa155, label="Probability of breaching 1.5°C", xlabel="Year", ylabel="Probability", title="Probability of Breaching 1.5°C and 2°C", legend=:bottomright, linewidth=4, color=:darkorange, linestyle=:dash)
plot!(fullcalendar, pa205, label="Probability of breaching 2°C", color=:red3, xlims=(2016, 2087), ylims=(0, 1), xticks=2016:10:2100, yticks=0:0.1:1, linewidth=4, linestyle=:dashdot)

probabilities_df.Probability_15_SSP5 = pa155
probabilities_df.Probability_20_SSP5 = pa205

CSV.write("results/ProbabilityPaths-ERF.csv", probabilities_df)

# Figure 4: Probability of breaching 1.5°C under different SSPs
probabilities_df = CSV.read("results/ProbabilityPaths-ERF.csv", DataFrame)
ppsp = plot(probabilities_df.Year, probabilities_df.Probability_15_SSP1, label="SSP1-1.9", xlabel="Year", ylabel="Probability", title="Probability of Breaching 1.5°C", legend=:right, linewidth=2, color=:darkorange, linestyle=:dash)
plot!(ppsp, probabilities_df.Year, probabilities_df.Probability_15_SSP2, label="SSP1-2.6", color=:green, linewidth=2, linestyle=:dash)
plot!(ppsp, probabilities_df.Year, probabilities_df.Probability_15_SSP3, label="SSP2-4.5", color=:red, linewidth=2, linestyle=:dash)
plot!(ppsp, probabilities_df.Year, probabilities_df.Probability_15_SSP4, label="SSP3-7.0", color=:purple, linewidth=2, linestyle=:dash)
plot!(ppsp, probabilities_df.Year, probabilities_df.Probability_15_SSP5, label="SSP5-8.5", color=:brown, linewidth=2, linestyle=:dash)
plot!(ppsp, xlims=(2016, 2087), ylims=(0, 1), xticks=2015:10:2100, yticks=0:0.1:1)

savefig("figures/Figure4.pdf")


# Table 3: Probability of breaching 1.5°C under different SSPs (1980 is used when no probability is found)
table_probabilities = DataFrame("Probability Level" => String[], "SSP1" => Int64[], "SSP2" => Int64[], "SSP3" => Int64[], "SSP4" => Int64[], "SSP5" => Int64[])

push!(table_probabilities, ["Above 1%", probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP1 .> 0.01), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP2 .> 0.01), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP3 .> 0.01), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP4 .> 0.01), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP5 .> 0.01), 1)]])
push!(table_probabilities, ["Above 10%", probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP1 .> 0.1), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP2 .> 0.1), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP3 .> 0.1), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP4 .> 0.1), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP5 .> 0.1), 1)]])
push!(table_probabilities, ["Above 25%", probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP1 .> 0.25), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP2 .> 0.25), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP3 .> 0.25), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP4 .> 0.25), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP5 .> 0.25), 1)]])
push!(table_probabilities, ["Above 50%", probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP1 .> 0.5), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP2 .> 0.5), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP3 .> 0.5), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP4 .> 0.5), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP5 .> 0.5), 1)]])
push!(table_probabilities, ["Above 95%", probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP1 .> 0.95), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP2 .> 0.95), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP3 .> 0.95), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP4 .> 0.95), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP5 .> 0.95), 1)]])
push!(table_probabilities, ["Above 99%", probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP1 .> 0.99), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP2 .> 0.99), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP3 .> 0.99), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP4 .> 0.99), 1)], probabilities_df.Year[something(findfirst(probabilities_df.Probability_15_SSP5 .> 0.99), 1)]])

CSV.write("results/Table3.csv", table_probabilities)
