# A Residual Prediction Test for the Well-Specification of Linear Instrumental Variable Models
This repository contains the code for reproducing the plots and the tables in the paper 
<i>Cyrill Scheidegger, Malte Londschien and Peter Bühlmann (2025). A residual prediction test for the well-specification of linear instrumental variable models, arXiv:2506.12771</i>.


- The files `SimulationSetup.R` and `AdditionalImplementations.R` contain code that is sourced in the simulations.

- The files `SimulationsH0.R`, `SimulationsHA.R` and `SimulationsPower.R` contain code to run the simulations.

- The file `Plots_for_Main.R` produces the figures in Sections 5.2 and 5.3 of the paper based on the simulation results.

- The file `RealDataAnalysis.R` conducts two real-data analyses and produces the figures in Section 5.4 of the paper and Section E of the supplement.

- The files `PlotsH0_for_Supplement.R`, `PlotsHA_for_Supplement.R` and `PlotsPower_for_Supplement.R` produce the plots in Sections D.3, D.4 and D.5 of the supplement based on the simulation results.

- The files `SimulationSetupCluster.R`, `SimulationsCluster.R` and `PlotsCluster.R` contain code to run the simulations and produce the plots in Section F of the supplement.
