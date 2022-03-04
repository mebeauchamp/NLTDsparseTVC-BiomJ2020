# NLTDsparseTVC-BiomJ2020
R code for simulations in article "Wang Y, Beauchamp ME, Abrahamowicz M. Non-linear and time-dependent effects of sparsely measured continuous time-varying covariates in time-to-event analysis. *Biom J* 2020;62(2):492-515. doi: 10.1002/bimj.201900042".
## Description
Original code developed by Yishu Wang, with some earlier parts developped by Willy Wynant.
Modifications, corrections, and extensions of the code by Marie-Eve Beauchamp.

For questions or comments about the code please contact Marie-Eve Beauchamp (marie-eve.beauchamp at rimuhc.ca).
 
The code has been written using R with the following version information:<br/>
- R version 3.3.1 (2016-06-21)<br/> 
- Platform x86_64-w64-mingw32/x64 (64-bit)<br/> 
- Using R packages:<br/> 
  - PermAlgo version 1.1<br/>
  - survival version 2.39-5
## Content
#### Main programs for simulations 
Code to reproduce the simulations of sections 3.2 and 3.3 of the manuscript, and section S.4 of the Supporting Information:

- `Scenario_A1.R`<br/>
- `Scenario_A2.R`<br/>
- `Scenario_A3.R`<br/>
- `Scenario_A4.R`<br/>
- `Scenario_A4_cubic.R`, for scenario A4 with models (1) and (2) estimated with cubic B-Splines (instead of quadratic splines)<br/>
- `Scenario_A2_withME.R`, for scenario A2 with measurement errors (ME) added and including the SIMEX-based correction

Each program above runs the simulations for the corresponding scenario (e.g., A1).
The programs call the core functions in `Functions.R` to implement the proposed models (1) and (2) of the manuscript.
The programs read the data in `SBPdata.csv`, which were used to generate the simulated data. 
Each program creates a workspace .RData, with the corresponding name (e.g. `Scenario_A1.RData`), including the simulation results.
#### `Functions.R`
Includes the functions to implement models (1) and (2) of the manuscript, to calculate the NL, TD and TEL effect estimates, as well as the internal functions required. This program is called by the other programs. 
#### `SBPdata.csv` 
Contains a dataset derived from data from the Framingham Heart Study and that is used in the simulations. For the simulated data, the independent variables were fixed, and taken from SBPdata, and the outcomes were generated.
#### `Figures.R` 
The program to reproduce Figures 1-5 of the manuscript and Figures S1-S5 of the Supporting Information,
from the simulation results saved by each simulation program.
