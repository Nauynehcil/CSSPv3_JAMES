This repository provides the newly developed modules within CSSPv3, i.e., irrigation, reservoir and glacier modules. And the processed outputs and scripts associated with the study: “Development and evaluation of the Conjunctive Surface–Subsurface Process Model version 3 (CSSPv3) at global scale”.

Reference:
[1] X. Yuan*., C. Li, P. Ji, XZ. Xi, SW. Li, Y. Jiao, H. Yi, XZ. Liang, QC. Zeng, 2026:  Development and evaluation of Conjunctive Surface-Subsurface Process Model version 3 (CSSPv3) at global scale. Journal of Advances in Modeling Earth Systems. Revised.
[2] Li, C.,  X. Yuan*,  Y. Jiao,  P. Ji,  Z. Huang,  2024: High-resolution Land Surface Modeling of the Irrigation Effects on Evapotranspiration over the Yellow River Basin. Journal of Hydrology, 633, 130986. 
[3] Xi, X., X. Yuan*, & Yi. Hao, 2026: High-resolution land surface modeling of climate and CO2 effects on ecosystem carbon-water coupling across the Qinghai-Tibet Plateau. Agricultural and Forest Meteorology, 385, 111195. https://doi.org/10.1016/j.agrformet.2026.111195

The main folders and files are described as follows:
1)	Matlab_toolbox_2_CSSPv3
Contains the required MATLAB toolboxes for running the provided scripts. 
Note that external toolboxes m_map and tight_subplot are not included and should be installed separately. 
Detailed information is provided in the Matlab Requirements section below.

2)	Scripts_and_Data
Contains the scripts and processed datasets used for model evaluation and figure generation in the manuscript.
Subfolders named 01_IA to 06_GPP_ER correspond to different evaluation sections of the paper. Each folder includes:
	plotting/analysis scripts
	processed data files (MAT format) used to generate the figures
Subfolder named Other contains additional scripts for key analysis procedures and intermediate processing steps.
The MATLAB data file mask.mat contains auxiliary information required for plotting, including: spatial masks, grid latitude and longitude, and so on.
 
1. Matlab Requirements
1.1 Software Version
The scripts were developed and tested in Matlab-R2024a. It is recommended to use this version to ensure all built-in functions and graphics features work as intended.
1.2 External Dependencies
The following toolboxes are authored by other researchers. Users are required to download them via the links provided below and add them to the MATLAB search path before execution:
1)	m_map: Used for geographic mapping and projections. 
2)	tight_subplot: Used for precise control over subplot margins and spacing. 
1.3 Custom Toolboxes 
The custom functions used in this study are provided in the 'Matlab_toolbox_2_CSSPv3' directory. Before running the code, please ensure that this folder and all its subfolders are added to your MATLAB path, or append the path settings at the end of the corresponding scripts:
1)	KGE.m: Evaluates model performance using the Kling–Gupta Efficiency (KGE), which balances correlation, bias, and variability.
2)	CC.m: Computes the Pearson correlation coefficient to measure the linear relationship between simulated and observed variables.
3)	RMSE.m: Calculates the Root Mean Square Error (RMSE) to quantify the magnitude of prediction errors.
4)	BIAS.m: Estimates the mean bias between simulations and observations, indicating systematic overestimation or underestimation.
5)	Colorbar: Provides color palette files (in .mat format) designed to mimic NCL (NCAR Command Language) color tables for high-quality spatial visualization.
6)	weightmean.m: Performs weighted mean calculations, such as area-weighted spatial averages.
7)	checkLeapYear.m: Identifies leap years to ensure temporal consistency in long-term time series analyses.

2. Scripts and Data
All plotting scripts are named using the format: Figure*_* where each script corresponds to a specific figure in the main manuscript and is used to reproduce the associated results and visualizations.
The data required for plotting are provided in MATLAB Data format (.mat) and can be loaded using the MATLAB command: load *.mat.

Note:
Due to GitHub file size limits, some .mat files are split into multiple volumes.
To Extract: Download all .z01, .z02 ... and .zip files into the same folder, then extract the .zip file.
