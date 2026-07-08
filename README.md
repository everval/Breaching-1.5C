[![DOI](https://zenodo.org/badge/1292616026.svg)](https://doi.org/10.5281/zenodo.21267065)

# README

## Title and Context
This folder contains the replicable code and data for the manuscript:

**Breaching 1.5°C: Give me the odds** by *J. Eduardo Vera-Valdés, Olivia Kvist*  

The scripts implement the statistical/econometric models and generate the figures and tables presented in the manuscript. They are organized to facilitate reproducibility of the results.

---

## Contents
The repository includes the following scripts:

- `Analysis-ERF.jl` (~11 seconds): Implements the regression model for the Effective Radiative Forcing (ERF) data and generates related figures and tables. Inputs are the yearly temperature anomalies datasets from HadCRUT5 [`HadCRUT 5.0 Global Annual Summary.csv`] and ERF data from Indicators of Global Climate Change 2024 [`ERF_best_aggregates_1750-2024.csv`, `ERF_p05_aggregates_1750-2024.csv`, `ERF_p95_aggregates_1750-2024.csv`]. Outputs include:
    - `Figure1.pdf` - Temperature anomalies, ERF, and fitted model.
    - `Figure2.pdf` - Forecasting temperature anomalies with prediction intervals.
    - `Figure10.pdf` - Diagnostic plots for the ERF model.
    - `Table1.csv` - Regression results for the ERF model.
    - `Table2.csv` - Coverage probabilities of prediction intervals.
    - `Table8.csv` - Stationarity tests for the ERF data.
    - `Table12.csv` - Diagnostic tests for the ERF model residuals.

- `Projections-ERF.jl` (~10 seconds): Implements the climate model and estimates temperature probability paths. Inputs are the yearly temperature anomalies datasets from HadCRUT5 [`HadCRUT 5.0 Global Annual Summary.csv`] and ERF data from Indicators of Global Climate Change 2024 [`ERF_best_aggregates_1750-2024.csv`, `ERF_p05_aggregates_1750-2024.csv`, `ERF_p95_aggregates_1750-2024.csv`]. For the projections, the code uses the Shared Socioeconomic Pathways (SSPs) scenarios [`data/SSPs` folder]. Outputs include:
    - `Figure3.pdf` - ERF projections and temperature anomaly prediction intervals for each SSP. 6 subplots are generated: projected ERF, plus projected temperature anomalies for each of the five SSP scenarios.
    - `Figure4.pdf` - Probability of breaching 1.5°C and 2°C for each SSP.
    - `Table3.csv` - Years when the probability of breaching 1.5°C exceeds 1%, 10%, 25%, 50%, 95%, and 99% for each SSP.
    - `ProbabilityPaths-ERF.csv` - Yearly probability of breaching 1.5°C and 2°C for each SSP, underlying Figure4.pdf.

- `Analysis-ONI.jl` (~11 seconds): Implements the Markov-switching model for the Oceanic Niño Index (ONI) data and generates related figures and tables. Inputs are the monthly ONI data from NOAA, loaded through the `Compiled_Global_Temperature_Data.csv` file. Outputs include:
    - `Table9.csv` - AIC values for different numbers of regimes in the Markov-switching model.
    - Estimation results of the best model (7 states) including transition probabilities are printed to the console and saved in `results/nino_model7.jld2` for later use.

- `Analysis-TSModel.jl` (~3 seconds): Implements the time series model for the temperature anomalies data and generates related figures and tables. Inputs are the yearly temperature anomalies datasets from HadCRUT5, GISTEMP, and Berkeley, loaded through the `Compiled_Global_Temperature_Data.csv` file. Outputs include:
    - `Figure11.pdf` - Residual diagnostics for HadCRUT dataset.
    - `Figure12.pdf` - Residual diagnostics for Berkeley dataset.
    - `Figure13.pdf` - Residual diagnostics for GISTEMP dataset.
    - `Table4-HadCRUT.csv`, `Table4-GISTEMP.csv`, `Table4-Berkeley.csv` - Regression results for each dataset.
    - `Table13-HadCRUT.csv`, `Table13-GISTEMP.csv`, `Table13-Berkeley.csv` - Residual diagnostics tests for each dataset.

- `Plot-Temperature-ONI.jl` (<1 second): Generates the plot of the temperature anomalies and ONI data used in the manuscript. Inputs are the yearly temperature anomalies datasets from HadCRUT5, GISTEMP, and Berkeley, and the monthly ONI data from NOAA, loaded through the `Compiled_Global_Temperature_Data.csv` file. Outputs include:
    - `Figure5.pdf` - Plot of temperature anomalies and ONI data.

- `Coverage-TSModel.jl` (~5.9 minutes): Implements the time series model for the temperature anomalies data and generates coverage figures. Inputs are the yearly temperature anomalies datasets from HadCRUT5, GISTEMP, and Berkeley, loaded through the `Compiled_Global_Temperature_Data.csv` file. Outputs include:
    - `Figure14(a).pdf`, `Figure14(b).pdf`, `Figure14(c).pdf`: Coverage plots for each dataset starting from 2000-01-01.
    - `Figure15(a).pdf`, `Figure15(b).pdf`, `Figure15(c).pdf`: Coverage plots for each dataset starting from 2016-11-01.
    - `Figure6(a).pdf`, `Figure6(b).pdf`, `Figure6(c).pdf`: Coverage plots for each dataset starting from 2010-01-01.
    - `Table5.csv`: Coverage probabilities for each dataset for all three starting dates.

- `ProbabilityPaths-PreParis.jl` (~25.2 minutes): Implements the climate model and estimates temperature probability paths for the period before the Paris Agreement. Inputs are the yearly temperature anomalies datasets from HadCRUT5, GISTEMP, and Berkeley, and the monthly ONI data from NOAA, loaded through the `Compiled_Global_Temperature_Data.csv` file. Outputs include:
    - `Table7.csv` - Years when the probability of breaching 1.5°C exceeds 1%, 50%, and 99% for each dataset before the Paris Agreement.
    - Probability paths for each dataset are generated and saved in the `results` folder. They are used to create the figures in the manuscript in the `ProbabilityPathsPlots.jl` script.

- `ProbabilityPaths.jl` (~17.4 minutes): Implements the climate model and estimates temperature probability paths for the full sample. Inputs are the yearly temperature anomalies datasets from HadCRUT5, GISTEMP, and Berkeley, and the monthly ONI data from NOAA, loaded through the `Compiled_Global_Temperature_Data.csv` file. Outputs include:
    - `Table6.csv` - Years when the probability of breaching 1.5°C and 2°C exceeds 1%, 50%, and 99% for each dataset after the Paris Agreement.
    - `Figure7(a).pdf`, `Figure7(b).pdf`, `Figure7(c).pdf`: Prediction interval plots for each dataset, extended forecast for the full sample.
    - Probability paths for each dataset are generated and saved in the `results` folder. They are used to create the figures in the manuscript in the `ProbabilityPathsPlots.jl` script.

- `ProbabilityPathsPlots.jl` (<1 second): Generates the plots of the temperature probability paths for the full sample and for the period before the Paris Agreement. Inputs are the probability paths generated by the `ProbabilityPaths.jl` and `ProbabilityPaths-PreParis.jl` scripts. Outputs include:
    - `Figure8.pdf` - Probability of breaching 1.5°C and 2°C for each dataset.
    - `Figure9.pdf` - Probability of breaching 1.5°C for each dataset before and after the Paris Agreement.

- `AdditionalFunctions.jl`: Contains additional functions used in the other scripts.

---

## Requirements
- **Language**: Julia (v1.12.6 was used for the reproduction of the results)
- **Key packages**: Listed in `Project.toml` and `Manifest.toml`

To reproduce the environment, use the provided `Project.toml` and `Manifest.toml` files.

---

## Reproducing the Results

1. Clone/unzip the repository to your local machine.
2. Open Julia and activate the environment using the command:
   ```julia
   using Pkg
   Pkg.activate("path/to/repository")
   Pkg.instantiate()
   ```
3. Open each julia script (`.jl` files) in a Julia IDE or Jupyter notebook, ensure that the required packages are installed, and execute the scripts in the order they are listed in the Contents section to generate the figures and tables. Alternatively, run all scripts in sequence with the `main.jl` script (~49 minutes total), which executes them in the correct order.
4. The output figures will be saved in the `figures` folder, and the output tables will be saved in the `results` folder.

---

### Author and System Information

The reproducibility package was created by [*J. Eduardo Vera-Valdés*](https://everval.github.io) on *2026-07-08*. Contact: [eduardo@math.aau.dk](mailto:eduardo@math.aau.dk).

The codes were run in a MacBook Pro M4 Pro with 32GB RAM and 8-core CPU. The codes were written for ease of understanding, not for computational efficiency. Hence, no parallelization or GPU acceleration was implemented. The codes were run in a single thread, and the execution time for each script are listed in their description. The total execution time for all scripts (run via `main.jl`) is approximately 49 minutes.

---

### Data Sources

The data used in this repository are publicly available and can be downloaded from the sources listed in the manuscript. Nevertheless, the data files are included in the `data` folder for convenience.

Additionally, the main data files for the time series model are obtained from the Quarto manuscript [Global Temperature Anomalies in Practice](https://everval.github.io/Global-Temperature-Anomalies/), where the data are updated, compiled and cleaned. 

---

### Updates and Dashboard

Updates to the results and figures using the latest data can be found in the [dashboard](https://everval.github.io/Dashboard-Breaching-1.5C/). The dashboard is updated every 2 to 3 months with the latest data and results. 


--- 

### Citation

If you use this code or results in your research, please cite the manuscript:

```bibtex
@article{Vera-Valdes2026,
  title={Breaching 1.5°C: Give me the odds},
  author={Vera-Valdés, J. Eduardo and Kvist, Olivia},
  journal={International Journal of Forecasting},
  year={2026},
  note={Forthcoming},
}
```

---

## License

This repository is licensed under the MIT License.
You are free to use, modify, and distribute the code provided that proper credit is given to the original authors.
