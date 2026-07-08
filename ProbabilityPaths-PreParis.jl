cd(@__DIR__)
using Pkg
Pkg.activate(pwd())
using CSV, DataFrames, Dates, Random, StatsPlots, HypothesisTests, JLD2, MarSwitching
include("AdditionalFunctions.jl")

all_data = CSV.read("data/Compiled_Global_Temperature_Data.csv", DataFrame)
h = 1000
nsim = 10^5

######## HadCRUT
dataselection = "HadCRUT" # Change here to select the dataset

column1 = dataselection*"_RawTemperature"
column2 = dataselection*"_Temp"

data = select(all_data, :Date, column1, column2, :ONI_Anomaly)
rename!(data, Symbol(column1) => :RawTemperature, Symbol(column2) => :Temp, :ONI_Anomaly => :SST )

startPA = Date("2016-11-01") # Paris Agreement enters into force
tempnino = dropmissing(data)
tempnino = tempnino[tempnino.Date .< startPA, :]

nino_model7 = MSModel(tempnino[!, :SST], 7);

T = size(tempnino.Temp, 1)

lmodel_exo = robust_est(tempnino.Temp, [ones(T, 1) collect(1:T) tempnino.SST]; verbose=false);

mband = 0.6;
d_est = exact_whittle_est(lmodel_exo.res; m=mband)
d_std = sqrt(exact_whittle_est_variance(length(lmodel_exo.res); m = mband))
d_tstat = d_est/d_std

T = length(tempnino.Date);
date_for = collect((tempnino.Date[1]+Dates.Month(T)):Month(1):(tempnino.Date[1]+Dates.Month(T - 1)+Dates.Month(h)));

matforecasts = zeros(h, nsim)
for ii = 1:nsim

  simul_nino = generate_msm(nino_model7, h)[1]
  matforecasts[:, ii] = trend_forecast_simulation(lmodel_exo, h, simul_nino, d_est, d_std).Yforecasterr

end

PA_start = Date("2006-11-01") # PA enters into force
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

results_probability_paths = DataFrame("Date (month)" => datejoin, "1.5°C Threshold" => vec(pa15))
CSV.write("results/ProbabilityPaths-"*dataselection*"-prePA.csv", results_probability_paths)

results_probabilities = DataFrame("Dataset" => String[], "Probability level and period" => String[], "1.5°C Threshold" => Date[])
push!(results_probabilities, [dataselection,  "Above 1%, 20-years avg.", datejoin[something(findfirst(pa15 .> 0.01), 1)]])
push!(results_probabilities, [dataselection, "Above 50%, 20-years avg.", datejoin[something(findfirst(pa15 .>= 0.5), 1)]])
push!(results_probabilities, [dataselection, "Above 99%, 20-years avg.", datejoin[something(findfirst(pa15 .>= 0.99), 1)]])

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

push!(results_probabilities, [dataselection, "Above 1%, 30-years avg.", datejoin30[something(findfirst(pa15_30 .> 0.01), 1)]])
push!(results_probabilities, [dataselection, "Above 50%, 30-years avg.", datejoin30[something(findfirst(pa15_30 .>= 0.5),1)]])
push!(results_probabilities, [dataselection, "Above 99%, 30-years avg.", datejoin30[something(findfirst(pa15_30 .>= 0.99),1)]])


######## GISTEMP
dataselection = "GISTEMP" # Change here to select the dataset

column1 = dataselection*"_RawTemperature"
column2 = dataselection*"_Temp"

data = select(all_data, :Date, column1, column2, :ONI_Anomaly)
rename!(data, Symbol(column1) => :RawTemperature, Symbol(column2) => :Temp, :ONI_Anomaly => :SST )

startPA = Date("2016-11-01") # Paris Agreement enters into force
tempnino = dropmissing(data)
tempnino = tempnino[tempnino.Date .< startPA, :]

T = size(tempnino.Temp, 1)

lmodel_exo = robust_est(tempnino.Temp, [ones(T, 1) collect(1:T) tempnino.SST]; verbose=false);

mband = 0.6;
d_est = exact_whittle_est(lmodel_exo.res; m=mband)
d_std = sqrt(exact_whittle_est_variance(length(lmodel_exo.res); m = mband))
d_tstat = d_est/d_std

T = length(tempnino.Date);
date_for = collect((tempnino.Date[1]+Dates.Month(T)):Month(1):(tempnino.Date[1]+Dates.Month(T - 1)+Dates.Month(h)));
matforecasts = zeros(h, nsim)

for ii = 1:nsim

  simul_nino = generate_msm(nino_model7, h)[1]
  matforecasts[:, ii] = trend_forecast_simulation(lmodel_exo, h, simul_nino, d_est, d_std).Yforecasterr

end

PA_start = Date("2006-11-01") # PA enters into force
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

results_probability_paths = DataFrame("Date (month)" => datejoin, "1.5°C Threshold" => vec(pa15))
CSV.write("results/ProbabilityPaths-"*dataselection*"-prePA.csv", results_probability_paths)

push!(results_probabilities, [dataselection,  "Above 1%, 20-years avg.", datejoin[something(findfirst(pa15 .> 0.01), 1)]])
push!(results_probabilities, [dataselection, "Above 50%, 20-years avg.", datejoin[something(findfirst(pa15 .>= 0.5), 1)]])
push!(results_probabilities, [dataselection, "Above 99%, 20-years avg.", datejoin[something(findfirst(pa15 .>= 0.99), 1)]])

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

push!(results_probabilities, [dataselection, "Above 1%, 30-years avg.", datejoin30[something(findfirst(pa15_30 .> 0.01), 1)]])
push!(results_probabilities, [dataselection, "Above 50%, 30-years avg.", datejoin30[something(findfirst(pa15_30 .>= 0.5),1)]])
push!(results_probabilities, [dataselection, "Above 99%, 30-years avg.", datejoin30[something(findfirst(pa15_30 .>= 0.99),1)]])


######## HadCRUT
dataselection = "Berkeley" # Change here to select the dataset

column1 = dataselection*"_RawTemperature"
column2 = dataselection*"_Temp"

data = select(all_data, :Date, column1, column2, :ONI_Anomaly)
rename!(data, Symbol(column1) => :RawTemperature, Symbol(column2) => :Temp, :ONI_Anomaly => :SST )

startPA = Date("2016-11-01") # Paris Agreement enters into force
tempnino = dropmissing(data)
tempnino = tempnino[tempnino.Date .< startPA, :]

T = size(tempnino.Temp, 1)

lmodel_exo = robust_est(tempnino.Temp, [ones(T, 1) collect(1:T) tempnino.SST]; verbose=false);

mband = 0.6;
d_est = exact_whittle_est(lmodel_exo.res; m=mband)
d_std = sqrt(exact_whittle_est_variance(length(lmodel_exo.res); m = mband))
d_tstat = d_est/d_std

T = length(tempnino.Date);
date_for = collect((tempnino.Date[1]+Dates.Month(T)):Month(1):(tempnino.Date[1]+Dates.Month(T - 1)+Dates.Month(h)));
matforecasts = zeros(h, nsim)

for ii = 1:nsim

  simul_nino = generate_msm(nino_model7, h)[1]
  matforecasts[:, ii] = trend_forecast_simulation(lmodel_exo, h, simul_nino, d_est, d_std).Yforecasterr

end

PA_start = Date("2006-11-01") # PA enters into force
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

results_probability_paths = DataFrame("Date (month)" => datejoin, "1.5°C Threshold" => vec(pa15))
CSV.write("results/ProbabilityPaths-"*dataselection*"-prePA.csv", results_probability_paths)

push!(results_probabilities, [dataselection,  "Above 1%, 20-years avg.", datejoin[something(findfirst(pa15 .> 0.01), 1)]])
push!(results_probabilities, [dataselection, "Above 50%, 20-years avg.", datejoin[something(findfirst(pa15 .>= 0.5), 1)]])
push!(results_probabilities, [dataselection, "Above 99%, 20-years avg.", datejoin[something(findfirst(pa15 .>= 0.99), 1)]])

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

push!(results_probabilities, [dataselection, "Above 1%, 30-years avg.", datejoin30[something(findfirst(pa15_30 .> 0.01), 1)]])
push!(results_probabilities, [dataselection, "Above 50%, 30-years avg.", datejoin30[something(findfirst(pa15_30 .>= 0.5),1)]])
push!(results_probabilities, [dataselection, "Above 99%, 30-years avg.", datejoin30[something(findfirst(pa15_30 .>= 0.99),1)]])


## Table 7: Probabilities of breaching 1.5°C threshold at the start of the Paris Agreement (2016-11-01) and at the start of the Paris Agreement (2006-11-01) for 20-year and 30-year averages, respectively.
CSV.write("results/Table7.csv", results_probabilities)