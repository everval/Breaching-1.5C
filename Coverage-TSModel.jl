cd(@__DIR__)
using Pkg
Pkg.activate(pwd())
using CSV, DataFrames, Dates, Random, StatsPlots, HypothesisTests, JLD2, MarSwitching
include("AdditionalFunctions.jl")

all_data = CSV.read("data/Compiled_Global_Temperature_Data.csv", DataFrame)

results = DataFrame(Dataset = String[], StartDate = String[], Level = String[], Coverage = Float64[])

###################### Coverage start date: 2000-01-01
PA_start = Date("2000-01-01")

######### HadCRUT
dataselection = "HadCRUT"  # Change here to select the dataset

column1 = dataselection*"_RawTemperature"
column2 = dataselection*"_Temp"

data = select(all_data, :Date, column1, column2, :ONI_Anomaly)
rename!(data, Symbol(column1) => :RawTemperature, Symbol(column2) => :Temp, :ONI_Anomaly => :SST )

tempnino = dropmissing(data)

nino_model7 = MSModel(tempnino[!, :SST], 7);

tempnino_pre_pa = filter(row -> row.Date < PA_start, tempnino)
last(tempnino_pre_pa, 5)

h = size(tempnino.Date, 1)-size(tempnino_pre_pa.Date, 1)+1

T = size(tempnino_pre_pa.Temp, 1)
lmodel_exo = robust_est(tempnino_pre_pa.Temp, [ones(T, 1) collect(1:T) tempnino_pre_pa.SST]; verbose=false);

mband = 0.6;
d_est = exact_whittle_est(lmodel_exo.res; m=mband)
d_std = sqrt(exact_whittle_est_variance(length(lmodel_exo.res); m = mband))
d_tstat = d_est/d_std


# label: forecast-simulations
nsim = 10^5
matforecasts = zeros(h, nsim)

for ii = 1:nsim

  simul_nino = generate_msm(nino_model7, h)[1]
  matforecasts[:, ii] = trend_forecast_simulation(lmodel_exo, h, simul_nino, d_est, d_std).Yforecasterr

end

quantilesforecasts = zeros(h, 9)

for ii = 1:h
    quantilesforecasts[ii, :] = quantile(matforecasts[ii, :], [0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995 ])
end

dates_forecast = collect(PA_start:Month(1):PA_start + Month(h-1))
plot(tempnino.Date, tempnino.Temp, label="Observed", xlabel="Date", ylabel="Temperature Anomaly (°C)", title="Forecasted Temperature with PI ("*dataselection*")", legend=:topleft,
xticks=(tempnino.Date[1:60:end], Dates.format.(tempnino.Date[1:60:end], "Y-m")) )
plot!(dates_forecast, quantilesforecasts[:, 5], label="Median Forecast", color=:red, linewidth=2, linestyle=:dot)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 4],quantilesforecasts[:, 6]), fillalpha=0.3, label="50% PI", color=:purple, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 3],quantilesforecasts[:, 7]), fillalpha=0.2, label="90% PI", color=:green, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 1],quantilesforecasts[:, 9]), fillalpha=0.1, label="99% PI", color=:red, linezalpha=0, linealpha=0)
plot!(xlims=(Date(1982,1,1),Date(2025,10,1)))

savefig("figures/Figure14(a).pdf")

actuals = tempnino.Temp[(end - h + 1):end]

res99 = sum((actuals .>= quantilesforecasts[:, 1]) .& (actuals .<= quantilesforecasts[:, 9])) / h
push!(results, (dataselection, string(PA_start), "99%", res99))
res95 = sum((actuals .>= quantilesforecasts[:, 2]) .& (actuals .<= quantilesforecasts[:, 8])) / h
push!(results, (dataselection, string(PA_start), "95%", res95))
res90 = sum((actuals .>= quantilesforecasts[:, 3]) .& (actuals .<= quantilesforecasts[:, 7])) / h
push!(results, (dataselection, string(PA_start), "90%", res90))
res50 = sum((actuals .>= quantilesforecasts[:, 4]) .& (actuals .<= quantilesforecasts[:, 6])) / h
push!(results, (dataselection, string(PA_start), "50%", res50))

######### GISTEMP
dataselection = "GISTEMP"  # Change here to select the dataset

column1 = dataselection*"_RawTemperature"
column2 = dataselection*"_Temp"

data = select(all_data, :Date, column1, column2, :ONI_Anomaly)
rename!(data, Symbol(column1) => :RawTemperature, Symbol(column2) => :Temp, :ONI_Anomaly => :SST )

tempnino = dropmissing(data)

tempnino_pre_pa = filter(row -> row.Date < PA_start, tempnino)
last(tempnino_pre_pa, 5)

h = size(tempnino.Date, 1)-size(tempnino_pre_pa.Date, 1)+1

T = size(tempnino_pre_pa.Temp, 1)
lmodel_exo = robust_est(tempnino_pre_pa.Temp, [ones(T, 1) collect(1:T) tempnino_pre_pa.SST]; verbose=false);

mband = 0.6;
d_est = exact_whittle_est(lmodel_exo.res; m=mband)
d_std = sqrt(exact_whittle_est_variance(length(lmodel_exo.res); m = mband))
d_tstat = d_est/d_std


# label: forecast-simulations
nsim = 10^5
matforecasts = zeros(h, nsim)

for ii = 1:nsim

  simul_nino = generate_msm(nino_model7, h)[1]
  matforecasts[:, ii] = trend_forecast_simulation(lmodel_exo, h, simul_nino, d_est, d_std).Yforecasterr

end

quantilesforecasts = zeros(h, 9)

for ii = 1:h
    quantilesforecasts[ii, :] = quantile(matforecasts[ii, :], [0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995 ])
end

dates_forecast = collect(PA_start:Month(1):PA_start + Month(h-1))
plot(tempnino.Date, tempnino.Temp, label="Observed", xlabel="Date", ylabel="Temperature Anomaly (°C)", title="Forecasted Temperature with PI ("*dataselection*")", legend=:topleft,
xticks=(tempnino.Date[1:60:end], Dates.format.(tempnino.Date[1:60:end], "Y-m")) )
plot!(dates_forecast, quantilesforecasts[:, 5], label="Median Forecast", color=:red, linewidth=2, linestyle=:dot)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 4],quantilesforecasts[:, 6]), fillalpha=0.3, label="50% PI", color=:purple, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 3],quantilesforecasts[:, 7]), fillalpha=0.2, label="90% PI", color=:green, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 1],quantilesforecasts[:, 9]), fillalpha=0.1, label="99% PI", color=:red, linezalpha=0, linealpha=0)
plot!(xlims=(Date(1982,1,1),Date(2025,10,1)))

savefig("figures/Figure14(b).pdf")

actuals = tempnino.Temp[(end - h + 1):end]

res99 = sum((actuals .>= quantilesforecasts[:, 1]) .& (actuals .<= quantilesforecasts[:, 9])) / h
push!(results, (dataselection, string(PA_start), "99%", res99))
res95 = sum((actuals .>= quantilesforecasts[:, 2]) .& (actuals .<= quantilesforecasts[:, 8])) / h
push!(results, (dataselection, string(PA_start), "95%", res95))
res90 = sum((actuals .>= quantilesforecasts[:, 3]) .& (actuals .<= quantilesforecasts[:, 7])) / h
push!(results, (dataselection, string(PA_start), "90%", res90))
res50 = sum((actuals .>= quantilesforecasts[:, 4]) .& (actuals .<= quantilesforecasts[:, 6])) / h
push!(results, (dataselection, string(PA_start), "50%", res50))


######### Berkeley
dataselection = "Berkeley"  # Change here to select the dataset

column1 = dataselection*"_RawTemperature"
column2 = dataselection*"_Temp"

data = select(all_data, :Date, column1, column2, :ONI_Anomaly)
rename!(data, Symbol(column1) => :RawTemperature, Symbol(column2) => :Temp, :ONI_Anomaly => :SST )

tempnino = dropmissing(data)

tempnino_pre_pa = filter(row -> row.Date < PA_start, tempnino)
last(tempnino_pre_pa, 5)

h = size(tempnino.Date, 1)-size(tempnino_pre_pa.Date, 1)+1

T = size(tempnino_pre_pa.Temp, 1)
lmodel_exo = robust_est(tempnino_pre_pa.Temp, [ones(T, 1) collect(1:T) tempnino_pre_pa.SST]; verbose=false);

mband = 0.6;
d_est = exact_whittle_est(lmodel_exo.res; m=mband)
d_std = sqrt(exact_whittle_est_variance(length(lmodel_exo.res); m = mband))
d_tstat = d_est/d_std


# label: forecast-simulations
nsim = 10^5
matforecasts = zeros(h, nsim)

for ii = 1:nsim

  simul_nino = generate_msm(nino_model7, h)[1]
  matforecasts[:, ii] = trend_forecast_simulation(lmodel_exo, h, simul_nino, d_est, d_std).Yforecasterr

end

quantilesforecasts = zeros(h, 9)

for ii = 1:h
    quantilesforecasts[ii, :] = quantile(matforecasts[ii, :], [0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995 ])
end

dates_forecast = collect(PA_start:Month(1):PA_start + Month(h-1))
plot(tempnino.Date, tempnino.Temp, label="Observed", xlabel="Date", ylabel="Temperature Anomaly (°C)", title="Forecasted Temperature with PI ("*dataselection*")", legend=:topleft,
xticks=(tempnino.Date[1:60:end], Dates.format.(tempnino.Date[1:60:end], "Y-m")) )
plot!(dates_forecast, quantilesforecasts[:, 5], label="Median Forecast", color=:red, linewidth=2, linestyle=:dot)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 4],quantilesforecasts[:, 6]), fillalpha=0.3, label="50% PI", color=:purple, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 3],quantilesforecasts[:, 7]), fillalpha=0.2, label="90% PI", color=:green, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 1],quantilesforecasts[:, 9]), fillalpha=0.1, label="99% PI", color=:red, linezalpha=0, linealpha=0)
plot!(xlims=(Date(1982,1,1),Date(2025,10,1)))

savefig("figures/Figure14(c).pdf")

actuals = tempnino.Temp[(end - h + 1):end]

res99 = sum((actuals .>= quantilesforecasts[:, 1]) .& (actuals .<= quantilesforecasts[:, 9])) / h
push!(results, (dataselection, string(PA_start), "99%", res99))
res95 = sum((actuals .>= quantilesforecasts[:, 2]) .& (actuals .<= quantilesforecasts[:, 8])) / h
push!(results, (dataselection, string(PA_start), "95%", res95))
res90 = sum((actuals .>= quantilesforecasts[:, 3]) .& (actuals .<= quantilesforecasts[:, 7])) / h
push!(results, (dataselection, string(PA_start), "90%", res90))
res50 = sum((actuals .>= quantilesforecasts[:, 4]) .& (actuals .<= quantilesforecasts[:, 6])) / h
push!(results, (dataselection, string(PA_start), "50%", res50))


###################### Coverage start date: 2016-11-01
PA_start = Date("2016-11-01")

######### HadCRUT
dataselection = "HadCRUT"  # Change here to select the dataset

column1 = dataselection*"_RawTemperature"
column2 = dataselection*"_Temp"

data = select(all_data, :Date, column1, column2, :ONI_Anomaly)
rename!(data, Symbol(column1) => :RawTemperature, Symbol(column2) => :Temp, :ONI_Anomaly => :SST )

tempnino = dropmissing(data)

nino_model7 = MSModel(tempnino[!, :SST], 7);

tempnino_pre_pa = filter(row -> row.Date < PA_start, tempnino)
last(tempnino_pre_pa, 5)

h = size(tempnino.Date, 1)-size(tempnino_pre_pa.Date, 1)+1

T = size(tempnino_pre_pa.Temp, 1)
lmodel_exo = robust_est(tempnino_pre_pa.Temp, [ones(T, 1) collect(1:T) tempnino_pre_pa.SST]; verbose=false);

mband = 0.6;
d_est = exact_whittle_est(lmodel_exo.res; m=mband)
d_std = sqrt(exact_whittle_est_variance(length(lmodel_exo.res); m = mband))
d_tstat = d_est/d_std


# label: forecast-simulations
nsim = 10^5
matforecasts = zeros(h, nsim)

for ii = 1:nsim

  simul_nino = generate_msm(nino_model7, h)[1]
  matforecasts[:, ii] = trend_forecast_simulation(lmodel_exo, h, simul_nino, d_est, d_std).Yforecasterr

end

quantilesforecasts = zeros(h, 9)

for ii = 1:h
    quantilesforecasts[ii, :] = quantile(matforecasts[ii, :], [0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995 ])
end

dates_forecast = collect(PA_start:Month(1):PA_start + Month(h-1))
plot(tempnino.Date, tempnino.Temp, label="Observed", xlabel="Date", ylabel="Temperature Anomaly (°C)", title="Forecasted Temperature with PI ("*dataselection*")", legend=:topleft,
xticks=(tempnino.Date[1:60:end], Dates.format.(tempnino.Date[1:60:end], "Y-m")) )
plot!(dates_forecast, quantilesforecasts[:, 5], label="Median Forecast", color=:red, linewidth=2, linestyle=:dot)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 4],quantilesforecasts[:, 6]), fillalpha=0.3, label="50% PI", color=:purple, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 3],quantilesforecasts[:, 7]), fillalpha=0.2, label="90% PI", color=:green, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 1],quantilesforecasts[:, 9]), fillalpha=0.1, label="99% PI", color=:red, linezalpha=0, linealpha=0)
plot!(xlims=(Date(1982,1,1),Date(2025,10,1)))

savefig("figures/Figure15(a).pdf")

actuals = tempnino.Temp[(end - h + 1):end]

res99 = sum((actuals .>= quantilesforecasts[:, 1]) .& (actuals .<= quantilesforecasts[:, 9])) / h
push!(results, (dataselection, string(PA_start), "99%", res99))
res95 = sum((actuals .>= quantilesforecasts[:, 2]) .& (actuals .<= quantilesforecasts[:, 8])) / h
push!(results, (dataselection, string(PA_start), "95%", res95))
res90 = sum((actuals .>= quantilesforecasts[:, 3]) .& (actuals .<= quantilesforecasts[:, 7])) / h
push!(results, (dataselection, string(PA_start), "90%", res90))
res50 = sum((actuals .>= quantilesforecasts[:, 4]) .& (actuals .<= quantilesforecasts[:, 6])) / h
push!(results, (dataselection, string(PA_start), "50%", res50))

######### GISTEMP
dataselection = "GISTEMP"  # Change here to select the dataset

column1 = dataselection*"_RawTemperature"
column2 = dataselection*"_Temp"

data = select(all_data, :Date, column1, column2, :ONI_Anomaly)
rename!(data, Symbol(column1) => :RawTemperature, Symbol(column2) => :Temp, :ONI_Anomaly => :SST )

tempnino = dropmissing(data)

tempnino_pre_pa = filter(row -> row.Date < PA_start, tempnino)
last(tempnino_pre_pa, 5)

h = size(tempnino.Date, 1)-size(tempnino_pre_pa.Date, 1)+1

T = size(tempnino_pre_pa.Temp, 1)
lmodel_exo = robust_est(tempnino_pre_pa.Temp, [ones(T, 1) collect(1:T) tempnino_pre_pa.SST]; verbose=false);

mband = 0.6;
d_est = exact_whittle_est(lmodel_exo.res; m=mband)
d_std = sqrt(exact_whittle_est_variance(length(lmodel_exo.res); m = mband))
d_tstat = d_est/d_std


# label: forecast-simulations
nsim = 10^5
matforecasts = zeros(h, nsim)

for ii = 1:nsim

  simul_nino = generate_msm(nino_model7, h)[1]
  matforecasts[:, ii] = trend_forecast_simulation(lmodel_exo, h, simul_nino, d_est, d_std).Yforecasterr

end

quantilesforecasts = zeros(h, 9)

for ii = 1:h
    quantilesforecasts[ii, :] = quantile(matforecasts[ii, :], [0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995 ])
end

dates_forecast = collect(PA_start:Month(1):PA_start + Month(h-1))
plot(tempnino.Date, tempnino.Temp, label="Observed", xlabel="Date", ylabel="Temperature Anomaly (°C)", title="Forecasted Temperature with PI ("*dataselection*")", legend=:topleft,
xticks=(tempnino.Date[1:60:end], Dates.format.(tempnino.Date[1:60:end], "Y-m")) )
plot!(dates_forecast, quantilesforecasts[:, 5], label="Median Forecast", color=:red, linewidth=2, linestyle=:dot)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 4],quantilesforecasts[:, 6]), fillalpha=0.3, label="50% PI", color=:purple, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 3],quantilesforecasts[:, 7]), fillalpha=0.2, label="90% PI", color=:green, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 1],quantilesforecasts[:, 9]), fillalpha=0.1, label="99% PI", color=:red, linezalpha=0, linealpha=0)
plot!(xlims=(Date(1982,1,1),Date(2025,10,1)))

savefig("figures/Figure15(b).pdf")

actuals = tempnino.Temp[(end - h + 1):end]

res99 = sum((actuals .>= quantilesforecasts[:, 1]) .& (actuals .<= quantilesforecasts[:, 9])) / h
push!(results, (dataselection, string(PA_start), "99%", res99))
res95 = sum((actuals .>= quantilesforecasts[:, 2]) .& (actuals .<= quantilesforecasts[:, 8])) / h
push!(results, (dataselection, string(PA_start), "95%", res95))
res90 = sum((actuals .>= quantilesforecasts[:, 3]) .& (actuals .<= quantilesforecasts[:, 7])) / h
push!(results, (dataselection, string(PA_start), "90%", res90))
res50 = sum((actuals .>= quantilesforecasts[:, 4]) .& (actuals .<= quantilesforecasts[:, 6])) / h
push!(results, (dataselection, string(PA_start), "50%", res50))


######### Berkeley
dataselection = "Berkeley"  # Change here to select the dataset

column1 = dataselection*"_RawTemperature"
column2 = dataselection*"_Temp"

data = select(all_data, :Date, column1, column2, :ONI_Anomaly)
rename!(data, Symbol(column1) => :RawTemperature, Symbol(column2) => :Temp, :ONI_Anomaly => :SST )

tempnino = dropmissing(data)

tempnino_pre_pa = filter(row -> row.Date < PA_start, tempnino)
last(tempnino_pre_pa, 5)

h = size(tempnino.Date, 1)-size(tempnino_pre_pa.Date, 1)+1

T = size(tempnino_pre_pa.Temp, 1)
lmodel_exo = robust_est(tempnino_pre_pa.Temp, [ones(T, 1) collect(1:T) tempnino_pre_pa.SST]; verbose=false);

mband = 0.6;
d_est = exact_whittle_est(lmodel_exo.res; m=mband)
d_std = sqrt(exact_whittle_est_variance(length(lmodel_exo.res); m = mband))
d_tstat = d_est/d_std


# label: forecast-simulations
nsim = 10^5
matforecasts = zeros(h, nsim)

for ii = 1:nsim

  simul_nino = generate_msm(nino_model7, h)[1]
  matforecasts[:, ii] = trend_forecast_simulation(lmodel_exo, h, simul_nino, d_est, d_std).Yforecasterr

end

quantilesforecasts = zeros(h, 9)

for ii = 1:h
    quantilesforecasts[ii, :] = quantile(matforecasts[ii, :], [0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995 ])
end

dates_forecast = collect(PA_start:Month(1):PA_start + Month(h-1))
plot(tempnino.Date, tempnino.Temp, label="Observed", xlabel="Date", ylabel="Temperature Anomaly (°C)", title="Forecasted Temperature with PI ("*dataselection*")", legend=:topleft,
xticks=(tempnino.Date[1:60:end], Dates.format.(tempnino.Date[1:60:end], "Y-m")) )
plot!(dates_forecast, quantilesforecasts[:, 5], label="Median Forecast", color=:red, linewidth=2, linestyle=:dot)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 4],quantilesforecasts[:, 6]), fillalpha=0.3, label="50% PI", color=:purple, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 3],quantilesforecasts[:, 7]), fillalpha=0.2, label="90% PI", color=:green, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 1],quantilesforecasts[:, 9]), fillalpha=0.1, label="99% PI", color=:red, linezalpha=0, linealpha=0)
plot!(xlims=(Date(1982,1,1),Date(2025,10,1)))

savefig("figures/Figure15(c).pdf")

actuals = tempnino.Temp[(end - h + 1):end]

res99 = sum((actuals .>= quantilesforecasts[:, 1]) .& (actuals .<= quantilesforecasts[:, 9])) / h
push!(results, (dataselection, string(PA_start), "99%", res99))
res95 = sum((actuals .>= quantilesforecasts[:, 2]) .& (actuals .<= quantilesforecasts[:, 8])) / h
push!(results, (dataselection, string(PA_start), "95%", res95))
res90 = sum((actuals .>= quantilesforecasts[:, 3]) .& (actuals .<= quantilesforecasts[:, 7])) / h
push!(results, (dataselection, string(PA_start), "90%", res90))
res50 = sum((actuals .>= quantilesforecasts[:, 4]) .& (actuals .<= quantilesforecasts[:, 6])) / h
push!(results, (dataselection, string(PA_start), "50%", res50))


###################### Coverage start date: 2010-01-01
PA_start = Date("2010-01-01")

######### HadCRUT
dataselection = "HadCRUT"  # Change here to select the dataset

column1 = dataselection*"_RawTemperature"
column2 = dataselection*"_Temp"

data = select(all_data, :Date, column1, column2, :ONI_Anomaly)
rename!(data, Symbol(column1) => :RawTemperature, Symbol(column2) => :Temp, :ONI_Anomaly => :SST )

tempnino = dropmissing(data)

nino_model7 = MSModel(tempnino[!, :SST], 7);

tempnino_pre_pa = filter(row -> row.Date < PA_start, tempnino)
last(tempnino_pre_pa, 5)

h = size(tempnino.Date, 1)-size(tempnino_pre_pa.Date, 1)+1

T = size(tempnino_pre_pa.Temp, 1)
lmodel_exo = robust_est(tempnino_pre_pa.Temp, [ones(T, 1) collect(1:T) tempnino_pre_pa.SST]; verbose=false);

mband = 0.6;
d_est = exact_whittle_est(lmodel_exo.res; m=mband)
d_std = sqrt(exact_whittle_est_variance(length(lmodel_exo.res); m = mband))
d_tstat = d_est/d_std


# label: forecast-simulations
nsim = 10^5
matforecasts = zeros(h, nsim)

for ii = 1:nsim

  simul_nino = generate_msm(nino_model7, h)[1]
  matforecasts[:, ii] = trend_forecast_simulation(lmodel_exo, h, simul_nino, d_est, d_std).Yforecasterr

end

quantilesforecasts = zeros(h, 9)

for ii = 1:h
    quantilesforecasts[ii, :] = quantile(matforecasts[ii, :], [0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995 ])
end

dates_forecast = collect(PA_start:Month(1):PA_start + Month(h-1))
plot(tempnino.Date, tempnino.Temp, label="Observed", xlabel="Date", ylabel="Temperature Anomaly (°C)", title="Forecasted Temperature with PI ("*dataselection*")", legend=:topleft,
xticks=(tempnino.Date[1:60:end], Dates.format.(tempnino.Date[1:60:end], "Y-m")) )
plot!(dates_forecast, quantilesforecasts[:, 5], label="Median Forecast", color=:red, linewidth=2, linestyle=:dot)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 4],quantilesforecasts[:, 6]), fillalpha=0.3, label="50% PI", color=:purple, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 3],quantilesforecasts[:, 7]), fillalpha=0.2, label="90% PI", color=:green, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 1],quantilesforecasts[:, 9]), fillalpha=0.1, label="99% PI", color=:red, linezalpha=0, linealpha=0)
plot!(xlims=(Date(1982,1,1),Date(2025,10,1)))

savefig("figures/Figure6(a).pdf")

actuals = tempnino.Temp[(end - h + 1):end]

res99 = sum((actuals .>= quantilesforecasts[:, 1]) .& (actuals .<= quantilesforecasts[:, 9])) / h
push!(results, (dataselection, string(PA_start), "99%", res99))
res95 = sum((actuals .>= quantilesforecasts[:, 2]) .& (actuals .<= quantilesforecasts[:, 8])) / h
push!(results, (dataselection, string(PA_start), "95%", res95))
res90 = sum((actuals .>= quantilesforecasts[:, 3]) .& (actuals .<= quantilesforecasts[:, 7])) / h
push!(results, (dataselection, string(PA_start), "90%", res90))
res50 = sum((actuals .>= quantilesforecasts[:, 4]) .& (actuals .<= quantilesforecasts[:, 6])) / h
push!(results, (dataselection, string(PA_start), "50%", res50))


######### GISTEMP
dataselection = "GISTEMP"  # Change here to select the dataset

column1 = dataselection*"_RawTemperature"
column2 = dataselection*"_Temp"

data = select(all_data, :Date, column1, column2, :ONI_Anomaly)
rename!(data, Symbol(column1) => :RawTemperature, Symbol(column2) => :Temp, :ONI_Anomaly => :SST )

tempnino = dropmissing(data)

tempnino_pre_pa = filter(row -> row.Date < PA_start, tempnino)
last(tempnino_pre_pa, 5)

h = size(tempnino.Date, 1)-size(tempnino_pre_pa.Date, 1)+1

T = size(tempnino_pre_pa.Temp, 1)
lmodel_exo = robust_est(tempnino_pre_pa.Temp, [ones(T, 1) collect(1:T) tempnino_pre_pa.SST]; verbose=false);

mband = 0.6;
d_est = exact_whittle_est(lmodel_exo.res; m=mband)
d_std = sqrt(exact_whittle_est_variance(length(lmodel_exo.res); m = mband))
d_tstat = d_est/d_std


# label: forecast-simulations
nsim = 10^5
matforecasts = zeros(h, nsim)

for ii = 1:nsim

  simul_nino = generate_msm(nino_model7, h)[1]
  matforecasts[:, ii] = trend_forecast_simulation(lmodel_exo, h, simul_nino, d_est, d_std).Yforecasterr

end

quantilesforecasts = zeros(h, 9)

for ii = 1:h
    quantilesforecasts[ii, :] = quantile(matforecasts[ii, :], [0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995 ])
end

dates_forecast = collect(PA_start:Month(1):PA_start + Month(h-1))
plot(tempnino.Date, tempnino.Temp, label="Observed", xlabel="Date", ylabel="Temperature Anomaly (°C)", title="Forecasted Temperature with PI ("*dataselection*")", legend=:topleft,
xticks=(tempnino.Date[1:60:end], Dates.format.(tempnino.Date[1:60:end], "Y-m")) )
plot!(dates_forecast, quantilesforecasts[:, 5], label="Median Forecast", color=:red, linewidth=2, linestyle=:dot)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 4],quantilesforecasts[:, 6]), fillalpha=0.3, label="50% PI", color=:purple, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 3],quantilesforecasts[:, 7]), fillalpha=0.2, label="90% PI", color=:green, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 1],quantilesforecasts[:, 9]), fillalpha=0.1, label="99% PI", color=:red, linezalpha=0, linealpha=0)
plot!(xlims=(Date(1982,1,1),Date(2025,10,1)))

savefig("figures/Figure6(b).pdf")


actuals = tempnino.Temp[(end - h + 1):end]

res99 = sum((actuals .>= quantilesforecasts[:, 1]) .& (actuals .<= quantilesforecasts[:, 9])) / h
push!(results, (dataselection, string(PA_start), "99%", res99))
res95 = sum((actuals .>= quantilesforecasts[:, 2]) .& (actuals .<= quantilesforecasts[:, 8])) / h
push!(results, (dataselection, string(PA_start), "95%", res95))
res90 = sum((actuals .>= quantilesforecasts[:, 3]) .& (actuals .<= quantilesforecasts[:, 7])) / h
push!(results, (dataselection, string(PA_start), "90%", res90))
res50 = sum((actuals .>= quantilesforecasts[:, 4]) .& (actuals .<= quantilesforecasts[:, 6])) / h
push!(results, (dataselection, string(PA_start), "50%", res50))

######### Berkeley
dataselection = "Berkeley"  # Change here to select the dataset

column1 = dataselection*"_RawTemperature"
column2 = dataselection*"_Temp"

data = select(all_data, :Date, column1, column2, :ONI_Anomaly)
rename!(data, Symbol(column1) => :RawTemperature, Symbol(column2) => :Temp, :ONI_Anomaly => :SST )

tempnino = dropmissing(data)

tempnino_pre_pa = filter(row -> row.Date < PA_start, tempnino)
last(tempnino_pre_pa, 5)

h = size(tempnino.Date, 1)-size(tempnino_pre_pa.Date, 1)+1

T = size(tempnino_pre_pa.Temp, 1)
lmodel_exo = robust_est(tempnino_pre_pa.Temp, [ones(T, 1) collect(1:T) tempnino_pre_pa.SST]; verbose=false);

mband = 0.6;
d_est = exact_whittle_est(lmodel_exo.res; m=mband)
d_std = sqrt(exact_whittle_est_variance(length(lmodel_exo.res); m = mband))
d_tstat = d_est/d_std


# label: forecast-simulations
nsim = 10^5
matforecasts = zeros(h, nsim)

for ii = 1:nsim

  simul_nino = generate_msm(nino_model7, h)[1]
  matforecasts[:, ii] = trend_forecast_simulation(lmodel_exo, h, simul_nino, d_est, d_std).Yforecasterr

end

quantilesforecasts = zeros(h, 9)

for ii = 1:h
    quantilesforecasts[ii, :] = quantile(matforecasts[ii, :], [0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995 ])
end

dates_forecast = collect(PA_start:Month(1):PA_start + Month(h-1))
plot(tempnino.Date, tempnino.Temp, label="Observed", xlabel="Date", ylabel="Temperature Anomaly (°C)", title="Forecasted Temperature with PI ("*dataselection*")", legend=:topleft,
xticks=(tempnino.Date[1:60:end], Dates.format.(tempnino.Date[1:60:end], "Y-m")) )
plot!(dates_forecast, quantilesforecasts[:, 5], label="Median Forecast", color=:red, linewidth=2, linestyle=:dot)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 4],quantilesforecasts[:, 6]), fillalpha=0.3, label="50% PI", color=:purple, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 3],quantilesforecasts[:, 7]), fillalpha=0.2, label="90% PI", color=:green, linezalpha=0, linealpha=0)
plot!(dates_forecast, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 1],quantilesforecasts[:, 9]), fillalpha=0.1, label="99% PI", color=:red, linezalpha=0, linealpha=0)
plot!(xlims=(Date(1982,1,1),Date(2025,10,1)))

savefig("figures/Figure6(c).pdf")


actuals = tempnino.Temp[(end - h + 1):end]

res99 = sum((actuals .>= quantilesforecasts[:, 1]) .& (actuals .<= quantilesforecasts[:, 9])) / h
push!(results, (dataselection, string(PA_start), "99%", res99))
res95 = sum((actuals .>= quantilesforecasts[:, 2]) .& (actuals .<= quantilesforecasts[:, 8])) / h
push!(results, (dataselection, string(PA_start), "95%", res95))
res90 = sum((actuals .>= quantilesforecasts[:, 3]) .& (actuals .<= quantilesforecasts[:, 7])) / h
push!(results, (dataselection, string(PA_start), "90%", res90))
res50 = sum((actuals .>= quantilesforecasts[:, 4]) .& (actuals .<= quantilesforecasts[:, 6])) / h
push!(results, (dataselection, string(PA_start), "50%", res50))


## Table 5: Coverage results for the three datasets across all coverage start dates
CSV.write("results/Table5.csv", results)