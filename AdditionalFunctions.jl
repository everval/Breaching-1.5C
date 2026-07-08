using LongMemory
using Statistics
using StatsBase
using Distributions
using GLM
using CovarianceMatrices
using LinearAlgebra
using Random
using Plots
using Plots.PlotMeasures

## Set default plot size and font
Random.seed!(123456)
theme(:ggplot2)
default(
        size = (600, 400),
        fontfamily = "Computer Modern",
        tickfontsize = 10,
        legendfontsize = 10,
        titlefontsize = 12,
        xlabelfontsize = 10,
        ylabelfontsize = 10,
        titlefontfamily = "Computer Modern",
        legendfontfamily = "Computer Modern",
        tickfontfamily = "Computer Modern",
        framestyle = :box,
        dpi = 600,
        margin = 3mm
) # Set default plot size and font

## Estimation functions

function robust_est(Y, X; verbose=false)
    fit = lm(X, Y)
    betavar = vcov(Bartlett{Andrews}(), fit)
    params = coef(fit)
    res = residuals(fit)
    rss = sum(res .^ 2)
    σ² = dispersion(fit)#rss / (length(Y) - size(X, 2))
    std_fit = sqrt.(diag(betavar))
    conf_int = [params .- 1.96*std_fit params .+ 1.96*std_fit]
    t_st = params ./ std_fit

    if verbose
        println("-----Trend Estimation Results-----")
        println("Estimated coefficients: ", params)
        println("Standard errors: ", std_fit)
        println("t-statistics: ", t_st)
        println("Residual sum of squares: ", rss)
        println("R-squared: ", r2(fit))
        println("σ²: ", σ²)
        println("AIC: ", aic(fit))
        println("BIC: ", bic(fit))
        println("R-squared: ", r2(fit))
    end

    return (β = params, σ² = σ², confint = conf_int, stdev = std_fit, t_stat = t_st, u=residuals(fit), Yfit = GLM.predict(fit), betavardiag = std_fit.^2, X = X, res=res, rsquared = r2(fit), aic=aic(fit), bic=bic(fit), betavar=betavar)
end

## Forecasting and simulation

### Forecast uncertainty
function trend_forecast_simulation(params, h, exog, dest, dva)
    T = size(params.res, 1)
    X = [ones(h, 1) collect((T+1):(T+h)) exog]
    tol = 1e-6

    defrac_raw = fracdiff(params.res, dest)
    defrac = defrac_raw .- mean(defrac_raw)
    stderr = std(defrac)
    nyerr = rand(Normal(0, stderr), h)
    fulerr = [defrac; nyerr] 
    
    drel_raw = rand(Normal(dest, dva)) 
    drel = clamp(drel_raw, tol, 0.5-tol)
    fracerr = fracdiff(fulerr, -drel)
    forecastfracerr = fracerr[(T+1):(T+h)]

    Σ = params.betavar
    Σ_pd = Symmetric(Σ + 1e-8 * I)  # small nugget for robustness
    mvndist = MvNormal(params.β, Σ_pd) 
    betas = rand(mvndist)

    Yforecastmean = X * betas
    Yforecasterr = Yforecastmean + forecastfracerr

    return (Yforecastmean=Yforecastmean, Yforecasterr=Yforecasterr)

end


## Dickey-Fuller test

function dickey_fuller_tests(series; det_specs = [:none, :constant], maxlags=3)
    results = DataFrame(
        lag = Int[],
        deterministic = Symbol[],
        statistic = Float64[],
        pvalue = Float64[]
    )

    # Loop through lags and deterministic specifications
    for det in det_specs
        for lag in 0:maxlags
            adf = ADFTest(series, det, lag)
            push!(results, (lag=lag, deterministic=det, statistic=adf.stat, pvalue=pvalue(adf)))
        end
    end

    return results
end


## Engle-Granger cointegration test

using HypothesisTests
using Distributions

const EG_CV_LEVELS = [0.01, 0.05, 0.10]

const EG_CV_COEFF = Dict(
    # ── K=1 (unit root / ADF, included for completeness) ────────────
    (1, :none,     0.01) => (-2.5658,  -1.960, -10.04),
    (1, :none,     0.05) => (-1.9393,  -0.398,   0.00),
    (1, :none,     0.10) => (-1.6156,  -0.181,   0.00),

    (1, :constant, 0.01) => (-3.4335,  -5.999, -29.25),
    (1, :constant, 0.05) => (-2.8621,  -2.738,  -8.36),
    (1, :constant, 0.10) => (-2.5671,  -1.438,  -4.48),

    # ── K=2 ─────────────────────────────────────────────────────────
    (2, :constant, 0.01) => (-3.9001, -10.534, -30.03),
    (2, :constant, 0.05) => (-3.3377,  -5.967,  -8.98),
    (2, :constant, 0.10) => (-3.0462,  -4.069,  -5.73),

    # ── K=3 ─────────────────────────────────────────────────────────
    (3, :constant, 0.01) => (-4.2981, -13.790, -46.37),
    (3, :constant, 0.05) => (-3.7429,  -8.352, -13.41),
    (3, :constant, 0.10) => (-3.4518,  -6.241,  -2.79),

    # ── K=4 ─────────────────────────────────────────────────────────
    (4, :constant, 0.01) => (-4.6493, -17.188, -59.20),
    (4, :constant, 0.05) => (-4.1000, -10.745, -21.57),
    (4, :constant, 0.10) => (-3.8110,  -8.317,  -5.19),
)


"""
    eg_critical_value(n, T; model=:constant, level=0.05)

MacKinnon (2010) response surface critical value for the EG test.
"""
function eg_critical_value(n::Int, T::Int;
                            model::Symbol = :constant,
                            level::Float64 = 0.05)
    key = (n, model, level)
    haskey(EG_CV_COEFF, key) || error(
        "No EG critical value for n=$n, model=$model, level=$level")
    β∞, β1, β2 = EG_CV_COEFF[key]
    return β∞ + β1/T + β2/T^2
end

"""
    eg_pvalue(τ, n, T; model=:constant)

Compute a p-value for an Engle-Granger test statistic `τ` by interpolating
across MacKinnon (2010) critical values at multiple quantile levels.

Strategy:
  1. Compute CVs at 7 quantile levels spanning [0.01, 0.90].
  2. The CV grid gives us (cv, p) pairs — i.e. points on the inverse CDF.
  3. Fit a monotone cubic interpolant and evaluate at τ.
  4. Extrapolate in the tails via a shifted/scaled normal approximation.
"""
function eg_pvalue(τ::Float64, n::Int, T::Int; model::Symbol = :constant)
    # Collect quantile levels that are tabulated for this (n, model)
    levels  = filter(l -> haskey(EG_CV_COEFF, (n, model, l)), EG_CV_LEVELS)
    cvs     = [eg_critical_value(n, T; model=model, level=l) for l in levels]

    # cvs are in DECREASING order (more negative = lower quantile)
    # we need increasing cv → increasing p
    # cvs[1] corresponds to level 0.01 (most negative), cvs[end] to 0.90
    # so cvs is already increasing (less negative toward the right tail)

    # --- tail extrapolation ---
    if τ <= cvs[1]
        # Deep left tail: use a stretched normal below the 1% point
        # Calibrated so that p(cv_0.01) ≈ 0.01
        z = (τ - cvs[1]) / abs(cvs[1])   # normalized distance into tail
        return max(0.001, levels[1] * exp(2.5 * z))
    end

    if τ >= cvs[end]
        # Right tail: not rejecting at all, p approaches 1
        z = (τ - cvs[end]) / abs(cvs[end])
        return min(0.999, levels[end] + (1.0 - levels[end]) * (1 - exp(-3.0 * z)))
    end

    # --- interior: linear interpolation in log-probability space ---
    # Find the bracket [cvs[i], cvs[i+1]] containing τ
    i = findlast(c -> c <= τ, cvs)
    i = clamp(i, 1, length(cvs) - 1)

    cv_lo, cv_hi = cvs[i],   cvs[i+1]
    p_lo,  p_hi  = levels[i], levels[i+1]

    # Interpolate linearly in (cv, log(p)) space for better tail behaviour
    log_p = log(p_lo) + (log(p_hi) - log(p_lo)) * (τ - cv_lo) / (cv_hi - cv_lo)
    return exp(log_p)
end;