cd(@__DIR__)
using Pkg
Pkg.activate(pwd())
using CSV, DataFrames, Dates, Random, StatsPlots, HypothesisTests, JLD2, MarSwitching
include("AdditionalFunctions.jl")

all_data = CSV.read("data/Compiled_Global_Temperature_Data.csv", DataFrame)
h = 800
nsim = 10^5

####### HadCRUT
dataselection = "HadCRUT"  # Change here to select the dataset

column1 = dataselection*"_RawTemperature"
column2 = dataselection*"_Temp"

data = select(all_data, :Date, column1, column2, :ONI_Anomaly)
rename!(data, Symbol(column1) => :RawTemperature, Symbol(column2) => :Temp, :ONI_Anomaly => :SST )

tempnino = dropmissing(data)

nino_model7 = MSModel(tempnino[!, :SST], 7);

T = size(tempnino.Temp, 1)

lmodel_exo = robust_est(tempnino.Temp, [ones(T, 1) collect(1:T) tempnino.SST]; verbose=false);

mband = 0.6;
d_est = exact_whittle_est(lmodel_exo.res; m=mband)
d_std = sqrt(exact_whittle_est_variance(length(lmodel_exo.res); m = mband))
d_tstat = d_est/d_std

T = length(tempnino.Date);
date_for = collect((tempnino.Date[1]+Dates.Month(T)):Month(1):(tempnino.Date[1]+Dates.Month(T - 1)+Dates.Month(h)));
fulldate = [tempnino.Date; date_for]


matforecasts = zeros(h, nsim)

for ii = 1:nsim

  simul_nino = generate_msm(nino_model7, h)[1]
  matforecasts[:, ii] = trend_forecast_simulation(lmodel_exo, h, simul_nino, d_est, d_std).Yforecasterr

end


ninodate = [data.Date; date_for]

quantilesforecasts = zeros(h, 9)

for ii = 1:h
    quantilesforecasts[ii, :] = quantile(matforecasts[ii, :], [0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995 ])
end

plot(tempnino.Date, tempnino.Temp, label="Observed", xlabel="Date", ylabel="Temperature Anomaly (°C)", title="Forecasted Temperature with PI ("*dataselection*")", legend=:topleft,
xticks=(tempnino.Date[1:60:end], Dates.format.(tempnino.Date[1:60:end], "Y-m")) )
plot!(date_for, quantilesforecasts[:, 5], label="Median Forecast", color=:red, linewidth=2, linestyle=:dot)
plot!(date_for, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 4],quantilesforecasts[:, 6]), fillalpha=0.3, label="50% PI", color=:purple, linezalpha=0, linealpha=0)
plot!(date_for, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 3],quantilesforecasts[:, 7]), fillalpha=0.2, label="90% PI", color=:green, linezalpha=0, linealpha=0)
plot!(date_for, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 1],quantilesforecasts[:, 9]), fillalpha=0.1, label="99% PI", color=:red, linezalpha=0, linealpha=0)
plot!(xlims=(Date(1982,1,1),Date(2082,10,1)) , xticks=(fulldate[37:120:end], Dates.format.(fulldate[37:120:end], "Y")) )

savefig("figures/Figure7(a).pdf")

PA_start = Date("2016-11-01") # PA enters into force
inicio = findfirst(isequal(PA_start),tempnino.Date)-119
fin = T + h - 120 # 10 years

datejoin = collect(PA_start-Month(119):Month(1):date_for[h-120])
meantemp = tempnino.Temp[inicio:T]
tsize = fin-inicio+1
dummies = zeros(tsize, nsim)

for jj = 1:nsim
    completo = [meantemp; matforecasts[:, jj]]
    for ii = 120:tsize
        dummies[ii, jj] = mean(completo[ii-119:ii+120])
    end
end

pa15 = mean(dummies[:, :] .>= 1.5, dims=2);
pa20 = mean(dummies[:, :] .>= 2, dims=2);


results_probability_paths = DataFrame("Date (month)" => datejoin, "1.5°C Threshold" => vec(pa15), "2°C Threshold" => vec(pa20))
CSV.write("results/ProbabilityPaths-"*dataselection*".csv", results_probability_paths)

results_probabilities = DataFrame("Dataset" => String[],"Probability level and period" => String[], "1.5°C Threshold" => Date[], "2°C Threshold" => Date[])

push!(results_probabilities, [dataselection, "Above 1%, 20-years avg.", datejoin[something(findfirst(pa15 .> 0.01), 1)], datejoin[something(findfirst(pa20 .> 0.01), 1)]])
push!(results_probabilities, [dataselection, "Above 50%, 20-years avg.", datejoin[something(findfirst(pa15 .>= 0.5), 1)], datejoin[something(findfirst(pa20 .>= 0.5), 1)]])
push!(results_probabilities, [dataselection, "Above 99%, 20-years avg.", datejoin[something(findfirst(pa15 .>= 0.99), 1)], datejoin[something(findfirst(pa20 .>= 0.99), 1)]])

PA_start = Date("2016-11-01") # PA enters into force
inicio30 = findfirst(isequal(PA_start),tempnino.Date)-179
fin30 = T + h - 180 # 15 years

datejoin30 = collect(PA_start-Month(179):Month(1):date_for[h-180])
meantemp30 = tempnino.Temp[inicio30:T]
tsize30 = fin30-inicio30+1
dummies30 = zeros(tsize30, nsim)

for jj = 1:nsim
    completo30 = [meantemp30; matforecasts[:, jj]]
    for ii = 180:tsize30
        dummies30[ii, jj] = mean(completo30[ii-179:ii+180])
    end
end

pa15_30 = mean(dummies30[:, :] .>= 1.5, dims=2);
pa20_30 = mean(dummies30[:, :] .>= 2, dims=2);

push!(results_probabilities, [dataselection, "Above 1%, 30-years avg.", datejoin30[something(findfirst(pa15_30 .> 0.01), 1)], datejoin30[something(findfirst(pa20_30 .> 0.01), 1)]])
push!(results_probabilities, [dataselection, "Above 50%, 30-years avg.", datejoin30[something(findfirst(pa15_30 .>= 0.5), 1)], datejoin30[something(findfirst(pa20_30 .>= 0.5), 1)]])
push!(results_probabilities, [dataselection, "Above 99%, 30-years avg.", datejoin30[something(findfirst(pa15_30 .>= 0.99), 1)], datejoin30[something(findfirst(pa20_30 .>= 0.99), 1)]])




####### GISTEMP
dataselection = "GISTEMP"  # Change here to select the dataset

column1 = dataselection*"_RawTemperature"
column2 = dataselection*"_Temp"

data = select(all_data, :Date, column1, column2, :ONI_Anomaly)
rename!(data, Symbol(column1) => :RawTemperature, Symbol(column2) => :Temp, :ONI_Anomaly => :SST )

tempnino = dropmissing(data)

T = size(tempnino.Temp, 1)

lmodel_exo = robust_est(tempnino.Temp, [ones(T, 1) collect(1:T) tempnino.SST]; verbose=false);

mband = 0.6;
d_est = exact_whittle_est(lmodel_exo.res; m=mband)
d_std = sqrt(exact_whittle_est_variance(length(lmodel_exo.res); m = mband))
d_tstat = d_est/d_std

T = length(tempnino.Date);
date_for = collect((tempnino.Date[1]+Dates.Month(T)):Month(1):(tempnino.Date[1]+Dates.Month(T - 1)+Dates.Month(h)));
fulldate = [tempnino.Date; date_for]

matforecasts = zeros(h, nsim)

for ii = 1:nsim

  simul_nino = generate_msm(nino_model7, h)[1]
  matforecasts[:, ii] = trend_forecast_simulation(lmodel_exo, h, simul_nino, d_est, d_std).Yforecasterr

end


ninodate = [data.Date; date_for]

quantilesforecasts = zeros(h, 9)

for ii = 1:h
    quantilesforecasts[ii, :] = quantile(matforecasts[ii, :], [0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995 ])
end

plot(tempnino.Date, tempnino.Temp, label="Observed", xlabel="Date", ylabel="Temperature Anomaly (°C)", title="Forecasted Temperature with PI ("*dataselection*")", legend=:topleft,
xticks=(tempnino.Date[1:60:end], Dates.format.(tempnino.Date[1:60:end], "Y-m")) )
plot!(date_for, quantilesforecasts[:, 5], label="Median Forecast", color=:red, linewidth=2, linestyle=:dot)
plot!(date_for, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 4],quantilesforecasts[:, 6]), fillalpha=0.3, label="50% PI", color=:purple, linezalpha=0, linealpha=0)
plot!(date_for, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 3],quantilesforecasts[:, 7]), fillalpha=0.2, label="90% PI", color=:green, linezalpha=0, linealpha=0)
plot!(date_for, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 1],quantilesforecasts[:, 9]), fillalpha=0.1, label="99% PI", color=:red, linezalpha=0, linealpha=0)
plot!(xlims=(Date(1982,1,1),Date(2082,10,1)) , xticks=(fulldate[37:120:end], Dates.format.(fulldate[37:120:end], "Y")) )

savefig("figures/Figure7(b).pdf")

PA_start = Date("2016-11-01") # PA enters into force
inicio = findfirst(isequal(PA_start),tempnino.Date)-119
fin = T + h - 120 # 10 years

datejoin = collect(PA_start-Month(119):Month(1):date_for[h-120])
meantemp = tempnino.Temp[inicio:T]
tsize = fin-inicio+1
dummies = zeros(tsize, nsim)

for jj = 1:nsim
    completo = [meantemp; matforecasts[:, jj]]
    for ii = 120:tsize
        dummies[ii, jj] = mean(completo[ii-119:ii+120])
    end
end

pa15 = mean(dummies[:, :] .>= 1.5, dims=2);
pa20 = mean(dummies[:, :] .>= 2, dims=2);


results_probability_paths = DataFrame("Date (month)" => datejoin, "1.5°C Threshold" => vec(pa15), "2°C Threshold" => vec(pa20))
CSV.write("results/ProbabilityPaths-"*dataselection*".csv", results_probability_paths)


push!(results_probabilities, [dataselection, "Above 1%, 20-years avg.", datejoin[something(findfirst(pa15 .> 0.01), 1)], datejoin[something(findfirst(pa20 .> 0.01), 1)]])
push!(results_probabilities, [dataselection, "Above 50%, 20-years avg.", datejoin[something(findfirst(pa15 .>= 0.5), 1)], datejoin[something(findfirst(pa20 .>= 0.5), 1)]])
push!(results_probabilities, [dataselection, "Above 99%, 20-years avg.", datejoin[something(findfirst(pa15 .>= 0.99), 1)], datejoin[something(findfirst(pa20 .>= 0.99), 1)]])

PA_start = Date("2016-11-01") # PA enters into force
inicio30 = findfirst(isequal(PA_start),tempnino.Date)-179
fin30 = T + h - 180 # 15 years

datejoin30 = collect(PA_start-Month(179):Month(1):date_for[h-180])
meantemp30 = tempnino.Temp[inicio30:T]
tsize30 = fin30-inicio30+1
dummies30 = zeros(tsize30, nsim)

for jj = 1:nsim
    completo30 = [meantemp30; matforecasts[:, jj]]
    for ii = 180:tsize30
        dummies30[ii, jj] = mean(completo30[ii-179:ii+180])
    end
end

pa15_30 = mean(dummies30[:, :] .>= 1.5, dims=2);
pa20_30 = mean(dummies30[:, :] .>= 2, dims=2);

push!(results_probabilities, [dataselection, "Above 1%, 30-years avg.", datejoin30[something(findfirst(pa15_30 .> 0.01), 1)], datejoin30[something(findfirst(pa20_30 .> 0.01), 1)]])
push!(results_probabilities, [dataselection, "Above 50%, 30-years avg.", datejoin30[something(findfirst(pa15_30 .>= 0.5), 1)], datejoin30[something(findfirst(pa20_30 .>= 0.5), 1)]])
push!(results_probabilities, [dataselection, "Above 99%, 30-years avg.", datejoin30[something(findfirst(pa15_30 .>= 0.99), 1)], datejoin30[something(findfirst(pa20_30 .>= 0.99), 1)]])



####### Berkeley
dataselection = "Berkeley"  # Change here to select the dataset

column1 = dataselection*"_RawTemperature"
column2 = dataselection*"_Temp"

data = select(all_data, :Date, column1, column2, :ONI_Anomaly)
rename!(data, Symbol(column1) => :RawTemperature, Symbol(column2) => :Temp, :ONI_Anomaly => :SST )

tempnino = dropmissing(data)

T = size(tempnino.Temp, 1)

lmodel_exo = robust_est(tempnino.Temp, [ones(T, 1) collect(1:T) tempnino.SST]; verbose=false);

mband = 0.6;
d_est = exact_whittle_est(lmodel_exo.res; m=mband)
d_std = sqrt(exact_whittle_est_variance(length(lmodel_exo.res); m = mband))
d_tstat = d_est/d_std

T = length(tempnino.Date);
date_for = collect((tempnino.Date[1]+Dates.Month(T)):Month(1):(tempnino.Date[1]+Dates.Month(T - 1)+Dates.Month(h)));
fulldate = [tempnino.Date; date_for]

matforecasts = zeros(h, nsim)

for ii = 1:nsim

  simul_nino = generate_msm(nino_model7, h)[1]
  matforecasts[:, ii] = trend_forecast_simulation(lmodel_exo, h, simul_nino, d_est, d_std).Yforecasterr

end


ninodate = [data.Date; date_for]

quantilesforecasts = zeros(h, 9)

for ii = 1:h
    quantilesforecasts[ii, :] = quantile(matforecasts[ii, :], [0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995 ])
end

plot(tempnino.Date, tempnino.Temp, label="Observed", xlabel="Date", ylabel="Temperature Anomaly (°C)", title="Forecasted Temperature with PI ("*dataselection*")", legend=:topleft,
xticks=(tempnino.Date[1:60:end], Dates.format.(tempnino.Date[1:60:end], "Y-m")) )
plot!(date_for, quantilesforecasts[:, 5], label="Median Forecast", color=:red, linewidth=2, linestyle=:dot)
plot!(date_for, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 4],quantilesforecasts[:, 6]), fillalpha=0.3, label="50% PI", color=:purple, linezalpha=0, linealpha=0)
plot!(date_for, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 3],quantilesforecasts[:, 7]), fillalpha=0.2, label="90% PI", color=:green, linezalpha=0, linealpha=0)
plot!(date_for, quantilesforecasts[:, 5], fillrange=(quantilesforecasts[:, 1],quantilesforecasts[:, 9]), fillalpha=0.1, label="99% PI", color=:red, linezalpha=0, linealpha=0)
plot!(xlims=(Date(1982,1,1),Date(2082,10,1)) , xticks=(fulldate[37:120:end], Dates.format.(fulldate[37:120:end], "Y")) )

savefig("figures/Figure7(c).pdf")

PA_start = Date("2016-11-01") # PA enters into force
inicio = findfirst(isequal(PA_start),tempnino.Date)-119
fin = T + h - 120 # 10 years

datejoin = collect(PA_start-Month(119):Month(1):date_for[h-120])
meantemp = tempnino.Temp[inicio:T]
tsize = fin-inicio+1
dummies = zeros(tsize, nsim)

for jj = 1:nsim
    completo = [meantemp; matforecasts[:, jj]]
    for ii = 120:tsize
        dummies[ii, jj] = mean(completo[ii-119:ii+120])
    end
end

pa15 = mean(dummies[:, :] .>= 1.5, dims=2);
pa20 = mean(dummies[:, :] .>= 2, dims=2);


results_probability_paths = DataFrame("Date (month)" => datejoin, "1.5°C Threshold" => vec(pa15), "2°C Threshold" => vec(pa20))
CSV.write("results/ProbabilityPaths-"*dataselection*".csv", results_probability_paths)


push!(results_probabilities, [dataselection, "Above 1%, 20-years avg.", datejoin[something(findfirst(pa15 .> 0.01), 1)], datejoin[something(findfirst(pa20 .> 0.01), 1)]])
push!(results_probabilities, [dataselection, "Above 50%, 20-years avg.", datejoin[something(findfirst(pa15 .>= 0.5), 1)], datejoin[something(findfirst(pa20 .>= 0.5), 1)]])
push!(results_probabilities, [dataselection, "Above 99%, 20-years avg.", datejoin[something(findfirst(pa15 .>= 0.99), 1)], datejoin[something(findfirst(pa20 .>= 0.99), 1)]])

PA_start = Date("2016-11-01") # PA enters into force
inicio30 = findfirst(isequal(PA_start),tempnino.Date)-179
fin30 = T + h - 180 # 15 years

datejoin30 = collect(PA_start-Month(179):Month(1):date_for[h-180])
meantemp30 = tempnino.Temp[inicio30:T]
tsize30 = fin30-inicio30+1
dummies30 = zeros(tsize30, nsim)

for jj = 1:nsim
    completo30 = [meantemp30; matforecasts[:, jj]]
    for ii = 180:tsize30
        dummies30[ii, jj] = mean(completo30[ii-179:ii+180])
    end
end

pa15_30 = mean(dummies30[:, :] .>= 1.5, dims=2);
pa20_30 = mean(dummies30[:, :] .>= 2, dims=2);

push!(results_probabilities, [dataselection, "Above 1%, 30-years avg.", datejoin30[something(findfirst(pa15_30 .> 0.01), 1)], datejoin30[something(findfirst(pa20_30 .> 0.01), 1)]])
push!(results_probabilities, [dataselection, "Above 50%, 30-years avg.", datejoin30[something(findfirst(pa15_30 .>= 0.5), 1)], datejoin30[something(findfirst(pa20_30 .>= 0.5), 1)]])
push!(results_probabilities, [dataselection, "Above 99%, 30-years avg.", datejoin30[something(findfirst(pa15_30 .>= 0.99), 1)], datejoin30[something(findfirst(pa20_30 .>= 0.99), 1)]])


### Table 6: Probabilities of breaching 1.5°C threshold at the last observation
CSV.write("results/Table6.csv", results_probabilities)