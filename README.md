# soda_drifter
Matlab scripts to track the trajectories of passive drifters for any locations in the North Pacific Ocean using velocity vectors from SODA3.12.2 

The soda_drifter.m script calculates the trajectories of drifters using either forward or reverse particle tracking. In addition to the soda_drifter.m script, the user also needs to have a copy of the following two netcdf files that are read by soda_drifter.m:

map1_velocity_2D_soda3_12_2_1980_2017_v20220421.nc
map1_2D_PlankTOM12_A_1980-2018_v20220419.nc

These two netcdf files are available to download from this link:



The u and v velocity vectors used by soda_driter.m are read from a netcdf file that was made by running the following scripts:

soda_v20.m 
soda_write_netcdf_v20.m 

These scripts extract 3D velocity vectors from the SODA3.12.2 1980-2017 monthly netcdf files that are available for download from https://www.soda.umd.edu/ . The scripts calculate vertically-integrated velocity vectors for 0-100m, 0-200m, 0-500m, and 0-1000m, and save them to a netcdf file that is read by soda_drifer.m

The animations produced by soda_drifter.m also include an optional map overlay of user selected biogeochemical tracers that are read by soda_drifter.m from a netcdf file that was produced by running the following scripts:

RECCAP_2021_3D_v42.m
RECCAP_2021_2D_v44_write_netcdf.m

These scripts extract various biogeochemical tracers from the PlankTOM12 1980-2018 monthly netcdf files that are available for download from http://greenocean-data.uea.ac.uk/RECCAP/2021/ . The scripts calcualte various vertically-integrated biogeochemical tracers for 0-100m that may be used for the animation map overlays

