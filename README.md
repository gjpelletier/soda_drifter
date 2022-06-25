soda_drifter.m

Matlab scripts to track the trajectories and biogeochemical exposure history of passive drifters for any location in the global oceans, and any time between 1980-2017, using velocity vectors from SODA3.12.2 and biogeochemical tracers from the World Ocean Atlas, TCO2_NNGv2, AT_NNGv2, and OceanSODA-ETHZ

The soda_drifter.m script calculates the trajectories of drifters using either forward or reverse particle tracking using Euler's method.

In addition to the soda_drifter.m script, the user also needs to have a copy of the following netcdf files that are read by soda_drifter.m:

% SODAv3122_2D_K100_1980_2017_v20220606.nc

% SODAv3122_2D_K200_1980_2017_v20220606.nc

% SODAv3122_2D_KT100_1980_2017_v20220606.nc

% SODAv3122_2D_KT200_1980_2017_v20220606.nc

% woa_2D_monthly_no3_o2_aou_v20220610.nc

% woa_nng_co2sys_2D_monthly_climatology_2005-2017_v20220608.nc

% woa_mobo_nng_co2sys_2D_monthly_climatology_2005-2017_v20220604.nc

% OceanSODA-ETHZ_GRaCER_v2021a_1982-2020.nc

% ncei_ocads_integrated_aragonite_v20220512.nc

These netcdf files are available to download from this link:

https://drive.google.com/drive/folders/1IEhqOa8RSszXh-aHRGV7t2JSzICBKce7?usp=sharing

The repository also includes a soda_drifter_utilities.zip file with several functions that are needed by soda_drifter.m

INSTALLATION INSTRUCTIONS

To install soda_drifter, do the following steps:

1) Copy the soda_drifter.m file into a working directory on your computer. 

2) Unzip the utilities scripts into that same working folder where you put the soda_drifter.m file. 

3) Create a subfolder named 'nc' in your working folder and unzip the netcdf files into that 'nc' folder. 

4) Edit line 161 in soda_drifter.m so it represents your working directory. For example, if your working directory is 'c:\soda_drifter', then line 161 should be edited so it is 'cd 'C:\soda_drifter';'

5) The user also needs to install the following three free matlab toolboxes from the followng links (m_map, gif, and export_fig) and edit their paths on line's 157-159 to reflect where they are stored on the user's computer:

addpath(genpath('c:/matlab/m_map/')) % functions for making maps in Mercator projection https://www.eoas.ubc.ca/~rich/map.html 

addpath(genpath('c:/matlab/gif/')) % make animated gif https://www.mathworks.com/matlabcentral/fileexchange/63239-gif 

addpath(genpath('c:/matlab/export_fig/')) % needed by gif https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig

Output files from soda_drifter will be stored in subfolders named png, gif, and csv in your working directory

If you have any questions, please contact Greg Pelletier at gjpelletier@gmail.com

