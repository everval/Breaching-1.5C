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
erf.ERF_p05 = rawerf_05[!, :total]
erf.ERF_p95 = rawerf_95[!, :total]

temperf_all = leftjoin(temp, erf, on = :Year)
dropmissing!(temperf_all)

temperf = temperf_all[temperf_all.Year .>= 1980, :]

# Table 1: Regression of temperature anomalies on ERF
reg = robust_est(temperf.Temp, [ones(length(temperf.ERF)) temperf.ERF]; verbose=true);

table1 = DataFrame(
    Term = ["Intercept", "ERF"],
    Estimate = reg.β,
    StdError = reg.stdev,
    tStat = reg.t_stat,
    CI_lower = reg.confint[:, 1],
    CI_upper = reg.confint[:, 2]
)

CSV.write("results/Table1.csv", table1)

# Table 8: Dickey-Fuller tests for stationarity of temperature anomalies, ERF, and residuals
df_temp = dickey_fuller_tests(temperf.Temp)
df_df = DataFrame("Lag" => df_temp.lag, "Specification" => string.(df_temp.deterministic), "Temperature" => df_temp.pvalue)

df_erf = dickey_fuller_tests(temperf.ERF)
df_df.ERF = df_erf.pvalue;

res_adf = dickey_fuller_tests(reg.res, det_specs=[:none, :constant], maxlags=3)
df_res = zeros(8)
for ii in 5:8
    df_res[ii] = eg_pvalue(res_adf.statistic[ii], 2, length(temperf.Temp); model=:constant)
end
df_df.Residuals = df_res

df_df

CSV.write("results/Table8.csv", df_df)

# Figure 1: Temperature anomalies, ERF, and fitted model
plot(temperf.Year, temperf.Temp, label="Temperature Anomalies", xlabel="", ylabel="", legend=:topleft, linewidth=2, color=:blue)
plot!(temperf.Year, reg.Yfit, linestyle=:dashdot, label="Fitted Model", linewidth=2, color=:red, xlabel="Year", ylabel="Temperature Anomaly (°C)", title="Temperature Anomalies, ERF, and Fitted Model")
plot!(twinx(), temperf.Year, temperf.ERF, label="Effective Radiative Forcing (Best)", ylabel="W/m2", linestyle=:dash, color=:black, linewidth=2, legend=:bottomright)
savefig("figures/Figure1.pdf")

# Figure 10: Residual diagnostics
res_raw = reg.res;
res_std = (res_raw .- mean(res_raw))./std(res_raw);

p1 = plot(temperf.Year, res_std, xlabel="Year", ylabel="Residuals", title="Standardised residuals", legend=false, linewidth=2)

p2 = scatter(reg.Yfit, res_std, xlabel="Fitted", ylabel="Standardised residuals ", title="Residuals vs fitted", legend=false)

p3 = histogram(res_std, xlabel="Residuals", ylabel="Frequency", title="Histogram", legend=false)

p4 = qqplot(res_std, Normal(0, 1), xlabel="Theoretical Quantiles", ylabel="Sample Quantiles", title="QQ-plot", linewidth=2)

mlag = 10
lags = 0:mlag
acs = autocor(res_std, lags)
n = length(res_std)
ci = 1.96 ./ sqrt.(n .- lags)

p5 = plot(lags, acs, seriestype=:stem, marker=:circle, markersize=3, xlabel="Lag", ylabel="ACF", title="Autocorrelation function", label="", linewidth=2)
plot!(lags, ci, linestyle=:dash, color=:red, label="95% CI", linewidth=2)
plot!(lags, -ci, linestyle=:dash, color=:red, label="", xticks=0:mlag, xminorticks=false, linewidth=2)

pacfs = pacf(res_std, lags)
p6 = plot(lags, pacfs, seriestype=:stem, marker=:circle, markersize=3, xlabel="Lag", ylabel="PACF", title="Partial autocorrelation function", label="", linewidth=2)
plot!(lags, ci, linestyle=:dash, color=:red, label="95% CI", linewidth=2)
plot!(lags, -ci, linestyle=:dash, color=:red, label="", xticks=0:mlag, xminorticks=false, linewidth=2)

plot(p1, p2, p3, p4, p5, p6, layout=(3,2), size=(650,800))
savefig("figures/Figure10.pdf")

## Table 12: Diagnostic tests for residuals (R-squared is from the robust regression output)
normality_test = ShapiroWilkTest(res_std)

dw_test = DurbinWatsonTest(reg.X, res_std)

lb_lags = ceil(Int, log(length(res_std)))
lb_test = LjungBoxTest(res_std, lb_lags, 3)

bp_test = BreuschPaganTest(reg.X, res_std .^ 2)

results_diagnostic_tests = DataFrame("Test" => String[], "Statistic" => Float64[], "p-value" => Float64[])
push!(results_diagnostic_tests, ["Shapiro-Wilk Normality Test", normality_test.W, pvalue(normality_test)])
push!(results_diagnostic_tests, ["Durbin-Watson Test", dw_test.DW, pvalue(dw_test)])
push!(results_diagnostic_tests, ["Ljung-Box Test ($(lb_lags) lags)", lb_test.Q, pvalue(lb_test)])
push!(results_diagnostic_tests, ["Breusch-Pagan Test", bp_test.lm, pvalue(bp_test)])
push!(results_diagnostic_tests, ["R-squared", reg.rsquared, NaN])

CSV.write("results/Table12.csv", results_diagnostic_tests)

# Figure 2: Forecasting temperature anomalies with prediction intervals
PA_start = 2015

temperf_sub = temperf[temperf.Year .< PA_start, :]
temperf_test = temperf[temperf.Year .>= PA_start, :]
n_test = size(temperf_test, 1)

erf_p95_p50 = temperf_test.ERF_p95 .- temperf_test.ERF
erf_p50_p05 = temperf_test.ERF .- temperf_test.ERF_p05
erf_p95_p05 = temperf_test.ERF_p95 .- temperf_test.ERF_p05

[ [mean(erf_p95_p50), mean(erf_p50_p05), mean(erf_p95_p05)],
  [std(erf_p95_p50), std(erf_p50_p05), std(erf_p95_p05)] ,
  [median(erf_p95_p50), median(erf_p50_p05), median(erf_p95_p05)]
  ]

[erf_p50_p05 erf_p95_p50 erf_p95_p05]

erf_std_p95 = erf_p95_p50 ./ 1.6448536269
erf_std_p05 = erf_p50_p05 ./ 1.6448536269

erf_measurements_std = (erf_std_p95 .+ erf_std_p05) ./ 2

reg_sub = robust_est(temperf_sub.Temp, [ones(length(temperf_sub.ERF)) temperf_sub.ERF]; verbose=false);

temperf_test = temperf[temperf.Year .>= PA_start, :]
n_test = size(temperf_test, 1)

nsim = 10^5

temp_pred = zeros(n_test, nsim)

for ii = 1:nsim
    #erf_sim = rand.(Uniform.(temperf_test.ERF_p05, temperf_test.ERF_p95), 1)
    erf_sim = rand.(Normal.(temperf_test.ERF, erf_measurements_std), 1)
    sigma = reg_sub.betavar
    sigma_pd = Symmetric(sigma + 1e-8 * I)
    betas_sim = rand(MvNormal(reg_sub.β, sigma_pd))
    temp_pred[:, ii] = betas_sim[1] .+ betas_sim[2] .* vcat(erf_sim...) + rand.(Normal(0, sqrt(reg_sub.σ²)), n_test)
end

quantilesforecasts = zeros(n_test, 9)

for ii = 1:n_test
    quantilesforecasts[ii, :] = quantile(temp_pred[ii, :], [0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995 ])
end

dates_forecast = collect(PA_start:temperf.Year[end])
p_pi = plot(temperf.Year, temperf.Temp, label="Observed", xlabel="Date", ylabel="Temperature Anomaly (°C)", title="Forecasted Temperature with Prediction Intervals", legend=:topleft)
plot!(dates_forecast, quantilesforecasts[:, 5], label="Median Forecast", color=:red, linewidth=2, linestyle=:dot)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 4],quantilesforecasts[:, 6]), fillalpha=0.3, label="50% Prediction Interval", color=:purple, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 3],quantilesforecasts[:, 7]), fillalpha=0.2, label="90% Prediction Interval", color=:green, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 1],quantilesforecasts[:, 9]), fillalpha=0.1, label="99% Prediction Interval", color=:red, linezalpha=0, linealpha=0)

savefig("figures/Figure2.pdf")

# Table 2: Coverage probabilities of prediction intervals
coverage_probs = zeros(5)
levels = [0.99, 0.95, 0.90, 0.50]

results = DataFrame(Level = String[], Coverage = Float64[])
actuals = temperf_test.Temp

res99 = sum((actuals .>= quantilesforecasts[:, 1]) .& (actuals .<= quantilesforecasts[:, 9])) / n_test
push!(results, ("99%", res99))
res95 = sum((actuals .>= quantilesforecasts[:, 2]) .& (actuals .<= quantilesforecasts[:, 8])) / n_test
push!(results, ("95%", res95))
res90 = sum((actuals .>= quantilesforecasts[:, 3]) .& (actuals .<= quantilesforecasts[:, 7])) / n_test
push!(results, ("90%", res90))
res50 = sum((actuals .>= quantilesforecasts[:, 4]) .& (actuals .<= quantilesforecasts[:, 6])) / n_test
push!(results, ("50%", res50))

CSV.write("results/Table2.csv", results)
