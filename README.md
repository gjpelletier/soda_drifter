# soda_drifter
Matlab scripts to track the trajectories of passive drifters for any location and any time between 1980-2017 in the North Pacific Ocean using velocity vectors from SODA3.12.2 

The soda_drifter.m script calculates the trajectories of drifters using either forward or reverse particle tracking using Euler's method. 

In addition to the soda_drifter.m script, the user also needs to have a copy of the following two netcdf files that are read by soda_drifter.m:

1) map1_velocity_2D_soda3_12_2_1980_2017_v20220421.nc

This netcdf file contains the vertically integrated velocity vectors. The user may select either 0-100m, 0-200m, 0-500m, or 0-1000m depth-integrated velocity vectors. These data were derived from 3D velocity vectors from the SODA3.12.2 1980-2017 monthly netcdf files that are available for download from https://www.soda.umd.edu/ . The scripts calculate vertically-integrated velocity vectors for 0-100m, 0-200m, 0-500m, and 0-1000m, and save them to a netcdf file that is read by soda_drifer.m

2) map1_2D_PlankTOM12_A_1980-2018_v20220419.nc

This netcdf file contains the optional map overlay data of biogeochemical tracers derived from output from the PlankTOM12 model. The user may select from a variety of biogeochemical tracer variables. These data were derived from the PlankTOM12 1980-2018 monthly netcdf files that are available for download from http://greenocean-data.uea.ac.uk/RECCAP/2021/ . The scripts calcualte various vertically-integrated biogeochemical tracers for 0-100m that may be used for the animation map overlays

These two netcdf files are available to download from this link:

https://drive.google.com/drive/folders/1IEhqOa8RSszXh-aHRGV7t2JSzICBKce7?usp=sharing

The repository also includes a soda_drifter_colormaps.zip file with two colormap functions that are needed to create the animations

To install soda_drifter, copy the soda_drifter.m file to a working directory on your computer. Then unzip the colormap files to that same working folder. Then create a subfolder named 'nc' in your working folder and unzip the two netcdf files in that 'nc' folder. Then edit line 87 in soda_drifter.m so it represents your working directory. For example, if your working directory is 'c:\soda_drifter', then line 87 should be edited so it is 'cd 'C:\soda_drifter';'
