# soda_drifter
Matlab scripts to track the trajectories of passive drifters for any location and any time between 1980-2017 in the North Pacific Ocean using velocity vectors from SODA3.12.2 

The soda_drifter.m script calculates the trajectories of drifters using either forward or reverse particle tracking using Euler's method. 

In addition to the soda_drifter.m script, the user also needs to have a copy of the following two netcdf files that are read by soda_drifter.m:

1) map1_velocity_2D_soda3_12_2_1980_2017_v20220421.nc

This netcdf file contains the vertically integrated velocity vectors. The user may select either 0-100m, 0-200m, 0-500m, or 0-1000m depth-integrated velocity vectors. These data were derived from 3D velocity vectors from the SODA3.12.2 1980-2017 monthly netcdf files that are available for download from https://www.soda.umd.edu/ . The scripts calculate vertically-integrated velocity vectors for 0-100m, 0-200m, 0-500m, and 0-1000m, and save them to a netcdf file that is read by soda_drifer.m

2) map1_2D_PlankTOM12_A_1980-2018_v20220429.nc

This netcdf file contains the optional map overlay data of biogeochemical tracers derived from output from the PlankTOM12 model. The user may select from a variety of biogeochemical tracer variables. These data were derived from the PlankTOM12 1980-2018 monthly netcdf files that are available for download from http://greenocean-data.uea.ac.uk/RECCAP/2021/ . The scripts calcualte various vertically-integrated biogeochemical tracers for 0-100m that may be used for the animation map overlays

These two netcdf files are available to download from this link:

https://drive.google.com/drive/folders/1IEhqOa8RSszXh-aHRGV7t2JSzICBKce7?usp=sharing

The repository also includes a soda_drifter_utilities.zip file with several functions that are needed by soda_drifter.m

To install soda_drifter, copy the soda_drifter.m file to a working directory on your computer. Then unzip the utilities files to that same working folder. Then create a subfolder named 'nc' in your working folder and unzip the two netcdf files in that 'nc' folder. Then edit line 87 in soda_drifter.m so it represents your working directory. For example, if your working directory is 'c:\soda_drifter', then line 87 should be edited so it is 'cd 'C:\soda_drifter';'

The user also needs to install the following three free matlab toolboxes from the followng links (m_map, gif, and export_fig) and edit their paths on line's 82-84 to reflect where they are stored on the user's computer:

addpath(genpath('c:/matlab/m_map/'))  		% functions for making maps in Mercator projection https://www.eoas.ubc.ca/~rich/map.html

addpath(genpath('c:/matlab/gif/'))  		% make animated gif https://www.mathworks.com/matlabcentral/fileexchange/63239-gif

addpath(genpath('c:/matlab/export_fig/'))  	% needed by gif https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig

The soda_drifter.m package includes a function named drifter.m that was developed as a stand-alone function to track the trajectories of passive drifter particles in the ocean. The drifter.m function is used by soda_drifter.m, and it may also be used outside of soda_drifter.m. It calculates the trajectories of drifters using either forward or reverse particle tracking using Euler's method.

Inputs to and outputs from the drifter.m function are as follows:

INPUTS

map_lon = vector of longitudes for the map extent (degE), must be equally spaced increments of degrees

map_lat = vector of latitudes for the map extent (degN), must be equally spaced increments of degrees

time = vector of continuous times for time series (datenum)

map_u = 3-d array (size of map_lon x map_lat x time) of east-west velocity vectors within the map extent (m/s positive east)

map_v = 3-d array (size of map_lon x map_lat x time) of north-south velocity vectors within the map extent (m/s positive north)

drifter_lon0 = vector of initial or final position longitudes of drifters (degE) 

drifter_lat0 = vector of initial or final position latitudes of drifters (degN)

direction = either 'forward' or 'reverse' for particle tracking:

direction = 'forward' calculates future trajectories starting at drifter_lon0 and drifter_lat0 intial drifter locations

direction = 'reverse' calculates past trajectories leading to the final drifter locations at drifter_lon0 and drifter_lat0

OUTPUTS

drifter_lon = 2-d array (size of numel(time) x numel(drifter_lon0)) of time series of trajectoy longitudes for each drifter (degE)

drifter_lat = 2-d array (size of numel(time) x numel(drifter_lon0)) of time series of trajectory latitudes for each drifter (degN)

The soda_drifter.m package is capable of extracting and interpolating the velocity data required by drifter.m from the SODA3.12.2 ocean model reanalysis data for any location in the North Pacific Ocean for any time period between 1980-2017. The required input arrays of velocity u and v vectors for drifter.m may also be created by the user from the output of any other available ocean model. 


