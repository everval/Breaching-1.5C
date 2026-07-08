cd(@__DIR__)
using Pkg
Pkg.activate(pwd())
using CSV, DataFrames, Dates, Random, StatsPlots, HypothesisTests
include("AdditionalFunctions.jl")

all_data = CSV.read("data/Compiled_Global_Temperature_Data.csv", DataFrame)


############# HadCRUT
dataselection = "HadCRUT"  # Change here to select the dataset

column1 = dataselection*"_RawTemperature"
column2 = dataselection*"_Temp"

data = select(all_data, :Date, column1, column2, :ONI_Anomaly)
rename!(data, Symbol(column1) => :RawTemperature, Symbol(column2) => :Temp, :ONI_Anomaly => :N34 )

tempnino = dropmissing(data)

lmodel_exo = robust_est(tempnino.Temp, [ones(size(tempnino, 1)) 1:size(tempnino, 1) tempnino.N34]; verbose=true);

mband = 0.6;
d_est = exact_whittle_est(lmodel_exo.res; m=mband)
d_std = sqrt(exact_whittle_est_variance(length(lmodel_exo.res); m = mband))
d_tstat = d_est/d_std

# Table 4: Estimation results of the long memory parameter d
param_names = ["Intercept", "Trend", "El Niño"]

results_coefficients = DataFrame(
    "Parameter" => param_names,
    "Coefficient" => lmodel_exo.β,
    "StdError" => lmodel_exo.stdev,
    "tStat" => lmodel_exo.t_stat
)
push!(results_coefficients, ["d", d_est, d_std, d_tstat])

CSV.write("results/Table4-"*dataselection*".csv", results_coefficients)


# Figure 11: Residual diagnostics HadCRUT
res_detrended = fracdiff(lmodel_exo.res, d_est);
res_std = (res_detrended .- mean(res_detrended))./std(res_detrended);

p1 = plot(tempnino.Date, res_std, xlabel="Date", ylabel="Residuals", title="Standardised residuals", legend=false,
xticks=(tempnino.Date[1:120:end], Dates.format.(tempnino.Date[1:120:end], "Y-m")) )

p2 = scatter(lmodel_exo.Yfit, res_std, xlabel="Fitted", ylabel="Standardised residuals ", title="Residuals vs fitted", legend=false)

p3 = histogram(res_std, xlabel="Residuals", ylabel="Frequency", title="Histogram", legend=false)

p4 = qqplot(res_std, Normal(0, 1), xlabel="Theoretical Quantiles", ylabel="Sample Quantiles", title="QQ-plot")

mlag = 30
lags = 0:mlag
acs = autocor(res_std, lags)

n = length(res_std)
ci = 1.96 ./ sqrt.(n .- lags)

p5 = plot(lags, acs, seriestype=:stem, marker=:circle, markersize=3, xlabel="Lag", ylabel="ACF", title="Autocorrelation function", label="")
plot!(lags, ci, linestyle=:dash, color=:red, label="95% CI")
plot!(lags, -ci, linestyle=:dash, color=:red, label="")

pacfs = pacf(res_std, lags)
p6 = plot(lags, pacfs, seriestype=:stem, marker=:circle, markersize=3, xlabel="Lag", ylabel="PACF", title="Partial autocorrelation function", label="")
plot!(lags, ci, linestyle=:dash, color=:red, label="95% CI")
plot!(lags, -ci, linestyle=:dash, color=:red, label="")

plot(p1, p2, p3, p4, p5, p6, layout=(3,2), size=(650,800))

savefig("figures/Figure11.pdf")

# Table 13: Residual diagnostics tests

normality_test = ShapiroWilkTest(res_std);

dw_test = DurbinWatsonTest(lmodel_exo.X, res_std);

lb_lags = ceil(Int, log(length(res_std)));
lb_test = LjungBoxTest(res_std, lb_lags, 3);

bp_test = BreuschPaganTest(lmodel_exo.X, res_std .^ 2);

results_diagnostic_tests = DataFrame("Test" => String[], "Statistic" => Float64[], "p-value" => Float64[])
push!(results_diagnostic_tests, ["Shapiro-Wilk Normality Test", normality_test.W, pvalue(normality_test)])
push!(results_diagnostic_tests, ["Durbin-Watson Test", dw_test.DW, pvalue(dw_test)])
push!(results_diagnostic_tests, ["Ljung-Box Test ("*string(lb_lags)*" lags)", lb_test.Q, pvalue(lb_test)])
push!(results_diagnostic_tests, ["Breusch-Pagan Test", bp_test.lm, pvalue(bp_test)])
push!(results_diagnostic_tests, ["R-squared", lmodel_exo.rsquared, NaN])

CSV.write("results/Table13-"*dataselection*".csv", results_diagnostic_tests)

####### Berkeley
dataselection = "Berkeley"  # Change here to select the dataset

column1 = dataselection*"_RawTemperature"
column2 = dataselection*"_Temp"

data = select(all_data, :Date, column1, column2, :ONI_Anomaly)
rename!(data, Symbol(column1) => :RawTemperature, Symbol(column2) => :Temp, :ONI_Anomaly => :N34 )

tempnino = dropmissing(data)

lmodel_exo = robust_est(tempnino.Temp, [ones(size(tempnino, 1)) 1:size(tempnino, 1) tempnino.N34]; verbose=true);

mband = 0.6;
d_est = exact_whittle_est(lmodel_exo.res; m=mband)
d_std = sqrt(exact_whittle_est_variance(length(lmodel_exo.res); m = mband))
d_tstat = d_est/d_std

# Table 4: Estimation results of the long memory parameter d
param_names = ["Intercept", "Trend", "El Niño"]

results_coefficients = DataFrame(
    "Parameter" => param_names,
    "Coefficient" => lmodel_exo.β,
    "StdError" => lmodel_exo.stdev,
    "tStat" => lmodel_exo.t_stat
)
push!(results_coefficients, ["d", d_est, d_std, d_tstat])

CSV.write("results/Table4-"*dataselection*".csv", results_coefficients)


# Figure 12: Residual diagnostics Berkeley

res_detrended = fracdiff(lmodel_exo.res, d_est);
res_std = (res_detrended .- mean(res_detrended))./std(res_detrended);

p1 = plot(tempnino.Date, res_std, xlabel="Date", ylabel="Residuals", title="Standardised residuals", legend=false,
xticks=(tempnino.Date[1:120:end], Dates.format.(tempnino.Date[1:120:end], "Y-m")) )

p2 = scatter(lmodel_exo.Yfit, res_std, xlabel="Fitted", ylabel="Standardised residuals ", title="Residuals vs fitted", legend=false)

p3 = histogram(res_std, xlabel="Residuals", ylabel="Frequency", title="Histogram", legend=false)

p4 = qqplot(res_std, Normal(0, 1), xlabel="Theoretical Quantiles", ylabel="Sample Quantiles", title="QQ-plot")

mlag = 30
lags = 0:mlag
acs = autocor(res_std, lags)

n = length(res_std)
ci = 1.96 ./ sqrt.(n .- lags)

p5 = plot(lags, acs, seriestype=:stem, marker=:circle, markersize=3, xlabel="Lag", ylabel="ACF", title="Autocorrelation function", label="")
plot!(lags, ci, linestyle=:dash, color=:red, label="95% CI")
plot!(lags, -ci, linestyle=:dash, color=:red, label="")

pacfs = pacf(res_std, lags)
p6 = plot(lags, pacfs, seriestype=:stem, marker=:circle, markersize=3, xlabel="Lag", ylabel="PACF", title="Partial autocorrelation function", label="")
plot!(lags, ci, linestyle=:dash, color=:red, label="95% CI")
plot!(lags, -ci, linestyle=:dash, color=:red, label="")

plot(p1, p2, p3, p4, p5, p6, layout=(3,2), size=(650,800))

savefig("figures/Figure12.pdf")

# Table 13: Residual diagnostics tests

normality_test = ShapiroWilkTest(res_std);

dw_test = DurbinWatsonTest(lmodel_exo.X, res_std);

lb_lags = ceil(Int, log(length(res_std)));
lb_test = LjungBoxTest(res_std, lb_lags, 3);

bp_test = BreuschPaganTest(lmodel_exo.X, res_std .^ 2);

results_diagnostic_tests = DataFrame("Test" => String[], "Statistic" => Float64[], "p-value" => Float64[])
push!(results_diagnostic_tests, ["Shapiro-Wilk Normality Test", normality_test.W, pvalue(normality_test)])
push!(results_diagnostic_tests, ["Durbin-Watson Test", dw_test.DW, pvalue(dw_test)])
push!(results_diagnostic_tests, ["Ljung-Box Test ("*string(lb_lags)*" lags)", lb_test.Q, pvalue(lb_test)])
push!(results_diagnostic_tests, ["Breusch-Pagan Test", bp_test.lm, pvalue(bp_test)])
push!(results_diagnostic_tests, ["R-squared", lmodel_exo.rsquared, NaN])

CSV.write("results/Table13-"*dataselection*".csv", results_diagnostic_tests)


####### GISTEMP
dataselection = "GISTEMP"  # Change here to select the dataset

column1 = dataselection*"_RawTemperature"
column2 = dataselection*"_Temp"

data = select(all_data, :Date, column1, column2, :ONI_Anomaly)
rename!(data, Symbol(column1) => :RawTemperature, Symbol(column2) => :Temp, :ONI_Anomaly => :N34 )

tempnino = dropmissing(data)

lmodel_exo = robust_est(tempnino.Temp, [ones(size(tempnino, 1)) 1:size(tempnino, 1) tempnino.N34]; verbose=true);

mband = 0.6;
d_est = exact_whittle_est(lmodel_exo.res; m=mband)
d_std = sqrt(exact_whittle_est_variance(length(lmodel_exo.res); m = mband))
d_tstat = d_est/d_std

# Table 4: Estimation results of the long memory parameter d
param_names = ["Intercept", "Trend", "El Niño"]

results_coefficients = DataFrame(
    "Parameter" => param_names,
    "Coefficient" => lmodel_exo.β,
    "StdError" => lmodel_exo.stdev,
    "tStat" => lmodel_exo.t_stat
)
push!(results_coefficients, ["d", d_est, d_std, d_tstat])

CSV.write("results/Table4-"*dataselection*".csv", results_coefficients)


# Figure 13: Residual diagnostics GISTEMP

res_detrended = fracdiff(lmodel_exo.res, d_est);
res_std = (res_detrended .- mean(res_detrended))./std(res_detrended);

p1 = plot(tempnino.Date, res_std, xlabel="Date", ylabel="Residuals", title="Standardised residuals", legend=false,
xticks=(tempnino.Date[1:120:end], Dates.format.(tempnino.Date[1:120:end], "Y-m")) )

p2 = scatter(lmodel_exo.Yfit, res_std, xlabel="Fitted", ylabel="Standardised residuals ", title="Residuals vs fitted", legend=false)

p3 = histogram(res_std, xlabel="Residuals", ylabel="Frequency", title="Histogram", legend=false)

p4 = qqplot(res_std, Normal(0, 1), xlabel="Theoretical Quantiles", ylabel="Sample Quantiles", title="QQ-plot")

mlag = 30
lags = 0:mlag
acs = autocor(res_std, lags)

n = length(res_std)
ci = 1.96 ./ sqrt.(n .- lags)

p5 = plot(lags, acs, seriestype=:stem, marker=:circle, markersize=3, xlabel="Lag", ylabel="ACF", title="Autocorrelation function", label="")
plot!(lags, ci, linestyle=:dash, color=:red, label="95% CI")
plot!(lags, -ci, linestyle=:dash, color=:red, label="")

pacfs = pacf(res_std, lags)
p6 = plot(lags, pacfs, seriestype=:stem, marker=:circle, markersize=3, xlabel="Lag", ylabel="PACF", title="Partial autocorrelation function", label="")
plot!(lags, ci, linestyle=:dash, color=:red, label="95% CI")
plot!(lags, -ci, linestyle=:dash, color=:red, label="")

plot(p1, p2, p3, p4, p5, p6, layout=(3,2), size=(650,800))

savefig("figures/Figure13.pdf")

# Table 13: Residual diagnostics tests

normality_test = ShapiroWilkTest(res_std);

dw_test = DurbinWatsonTest(lmodel_exo.X, res_std);

lb_lags = ceil(Int, log(length(res_std)));
lb_test = LjungBoxTest(res_std, lb_lags, 3);

bp_test = BreuschPaganTest(lmodel_exo.X, res_std .^ 2);

results_diagnostic_tests = DataFrame("Test" => String[], "Statistic" => Float64[], "p-value" => Float64[])
push!(results_diagnostic_tests, ["Shapiro-Wilk Normality Test", normality_test.W, pvalue(normality_test)])
push!(results_diagnostic_tests, ["Durbin-Watson Test", dw_test.DW, pvalue(dw_test)])
push!(results_diagnostic_tests, ["Ljung-Box Test ("*string(lb_lags)*" lags)", lb_test.Q, pvalue(lb_test)])
push!(results_diagnostic_tests, ["Breusch-Pagan Test", bp_test.lm, pvalue(bp_test)])
push!(results_diagnostic_tests, ["R-squared", lmodel_exo.rsquared, NaN])

CSV.write("results/Table13-"*dataselection*".csv", results_diagnostic_tests)