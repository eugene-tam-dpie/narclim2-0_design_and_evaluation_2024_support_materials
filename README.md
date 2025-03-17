# Supporting information for "Design, evaluation and future projections of the NARCliM2.0 CMIP6-CORDEX Australasia regional climate ensemble" #

## Executive Summary ##

This repository contains supporting information for the manuscript titled "Design, evaluation and future projections of the NARCliM2.0 CMIP6-CORDEX Australasia regional climate ensemble", by Giovanni Di Virgilio, Jason Evans, Fei Ji, Eugene Tam, Jatin Kala, Julia Andrys, Christopher Thomas, Dipayan Choudhury, Carlos Rocha, Stephen White, Yue Li, Moutassem El Rafei, Rishav Goyal and Matthew L. Riley, submitted to Geoscientific Model Development on 2024-05-05.

The manuscript is available with the DOI [10.5194/gmd-18-671-2025]( https://doi.org/10.5194/gmd-18-671-2025)

We ask everyone who uses or adapts our scripts to acknowledge the NSW Department of Climate Change, Energy, Environment and Water 2024 (NSW DCCEEW), not to be mistaken with our Federal counterpart.

Additional information can be found on the [NSW Climate Data Portal](https://climatedata-beta.environment.nsw.gov.au/)

## Repository Structure ##

The included supported materials include

* **WRF code and configuration**
  * **Snapshots of the version of WRF/WPS source code** (within the `WRF/repo_snaphots` subdirectory) used for GCMs with
    * standard and no leap calendars
    * 360 day calendars
  * **Namelists (configuration files)** used to generate boundary/initial conditions from metgrid outputs as well as downscale the GCMs using the boundary/initial condition files generated within the `WRF/namelists` subdirectory.
* **GCM pre-processing scripts and configuration** including
  * **NCL preprocessing scripts** that directly processed the GCM outputs and generate "intermediate" files within the `gcm_preprocessing/ncl_prepare_intermediate_files` subdirectory
  * **Namelists (configuration files)** used to convert the intermediate files into metgrid files within the `gcm_preprocessing/make_metgrids`
* Boundary/Initial file checker which ensures the boundary/initial conditions generated does have periods that contain significant outliers which can be found within the `miscellaneous` subdirectory. This isn't strictly necessary to run WRF but provides additional confidence in the downloaded/generated data before performing the more time consuming steps.

## Acknowledgements ##

We acknowledge

* Climate Change Fund
* Other funding bodies

## References ##

1. Di Virgilio Giovanni, F. Ji, E. Tam, N. Nishant, J. P. Evans, C. Thomas, M. Riley, K. Beyer, M. Grose, S Narsey, F Delage. (2022). Selecting CMIP6 GCMs for CORDEX dynamical downscaling: model performance, independence, and climate change signals, Earth's Future, DOI: [10.1029/2021EF002625](https://dx.doi.org/10.1029/2021EF002625)

2. Di Virgilio Giovanni, F. Ji, E. Tam, J. P. Evans, J. Kala, J. Andrys C. Thomas, D. Choudhury, C. Rocha, Y. Li, M. Riley (Under review) Evaluation of CORDEX ERA5-forced ‘NARCliM2.0’ regional climate models over Australia using the Weather Research and Forecasting (WRF) model version 4.1.2, Geoscientific Model Development, DOI: [10.5194/gmd-2024-41](https://dx.doi.org/10.5194/gmd-2024-41)

3. Di Virgilio Giovanni, J. P. Evans, F. Ji, E. Tam, J. Kala, J. Andrys, C. Thomas, D. Choudhury, C. Rocha, S. White, Y. Li, M. El Rafei, R. Goyal and M. Riley (Under review) Design, evaluation and future projections of the NARCliM2.0 CORDEX-CMIP6 Australasia regional climate ensemble, Geoscientific Model Development, DOI: [10.5194/gmd-2024-87](https://doi.org/10.5194/gmd-2024-87)
