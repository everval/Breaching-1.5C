cd(@__DIR__)
using Pkg
Pkg.activate(pwd())
using CSV, DataFrames, Dates, Random, StatsPlots, HypothesisTests
include("AdditionalFunctions.jl")

all_data = CSV.read("data/Compiled_Global_Temperature_Data.csv", DataFrame)

# Options: HadCRUT, Berkeley, GISTEMP
dataselection = "HadCRUT"  # Change here to select the dataset

column1 = "HadCRUT_Temp"
column2 = "Berkeley_Temp"
column3 = "GISTEMP_Temp"
column4 = "ONI_Anomaly"

data = select(all_data, :Date, column1, column2, column3, column4)

rename!(data, :HadCRUT_Temp => :HadCRUT, :Berkeley_Temp => :Berkeley, :GISTEMP_Temp => :GISTEMP, :ONI_Anomaly => :ONI )
dropmissing!(data)

T = length(data.Date)
indicator = zeros(Int, T)

for i in 1:T
    if data.ONI[i] >= 0.5
        indicator[i] = 1
    elseif data.ONI[i] <= -0.5
        indicator[i] = -1
    else
        indicator[i] = 0
    end
end

data[!, :Indicator] = indicator

p = plot(data.Date, data.HadCRUT, label="HadCRUT", xlabel="Date (monthly)", ylabel="°C", title="Temperature Anomalies", legend=:topleft, linewidth=1.5, linestyle=:dash)
plot!(data.Date, data.Berkeley, label="Berkeley", linewidth=1.5, linestyle=:dot)
plot!(data.Date, data.GISTEMP, label="GISTEMP", linewidth=1.5, linestyle=:dashdot)

function shade_indicator_spans!(p, dates, indicator)
    i = 1
    while i <= length(indicator)
        current_val = indicator[i]
        if current_val in (-1, 1)
            start_idx = i
            while i <= length(indicator) && indicator[i] == current_val
                i = i + 1
            end
            stop_idx = i - 1
            if (stop_idx - start_idx) >= 4
                vspan!(p, [dates[start_idx], dates[stop_idx]], color=current_val == 1 ? :red : :blue, alpha=0.1, label="")
            end
        else
            i = i + 1
        end
    end
end

shade_indicator_spans!(p, data.Date, indicator)

plot!(p, xlims=(Date(1982, 1, 1), Date(2025, 1, 1)), xticks=(data.Date[61:90:end], Dates.format.(data.Date[61:90:end], "mm/yyyy")))

savefig("figures/Figure5.pdf")

