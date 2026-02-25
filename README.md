# DSC-Wolbachia
This repository contains codes in R files for "Heterogeneous suppressive effect of Wolbachia incompatible insect technique coupled with sterile insect technique across time and historical Ae. aegypti abundance - using distributional synthetic controls" only.

The data preprocessing and analysis were run using **R version 4.4.1**. A list of package dependencies and their versions is provided below.

| Package name | Version |
| --- | --- |
|dplyr|1.2.0|
|ggplot2|3.5.2|
|parallel|4.5.2|
|CVXR|1.0.15|
|osqp|0.6.3.3|
|pbmcapply|1.5.1|
|ECOSolverR|0.5.5|
|frechet|0.3.0|
|msm|1.8.2|
|patchwork|1.3.2|
|RColorBrewer|1.1.3|

Code File description:
- "**Preprocessing_surveillance_data.R**". This file is used for data aggregation and preprocessing.
- "**DSC_aegy_sequential.R**". The main file for analysis. Serial computation for sectors is proceeded in this file.
- "**DSC_aegy_sector.R**". The source file for analysis. This file shows how the analysis is proceeded for each sector step by step.
- "**DSC_aegy_aux_func.R**". The source file containing all auxiliary functions.
- (Other files will be updated later)

Data File description: (**Data is not provided in this repository**)
- "**YYYY_QQ_weekly_trap_albo.xlsx**" which "**YYYY**" is replaced by the year of data recording and "**QQ**" is replaced by the quarter of data recording. These files contains the number of wild-type ***Ae. aegypti*** and ***Ae. albopictus*** caught in Gravitraps in an epidemiological week.
- "**Release_matrix_2016_2022ew26.xlsx**". This file contains weeks of releasing <i>Wolbachia</i>-infected male mosquitoes in intervened areas.
- "**4-levels_Sectors.csv**". This file contains the type of sectors (intervened or not; adjacent to the other type of sectors or not) every week throughout the study period.
