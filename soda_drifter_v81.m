
% - - - - -
% soda_drifter.m
% Trajectories of passive drifters using monthly velocity u and v from SODA3.12.2 
% using Euler's method with the option to use either forward or reverse tracking,
% and optional overlay of biogeochemical tracers
% - - - - -
% Copyright (c) 2022 Greg Pelletier

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% - - - - -
%
% version history

% v81 flip jet2 in southern hemisphere for animations
% v80 made extra plots optional
% v79 added woa climatology of no3, o2, and aou for optional overlay
% v78 added choice of either woa_nng or woa_mobo_nng climatology for overlay biogeochem 
% v77 updated drifter_v04.m
% v76 use netcdf of entire OceanSODA model domain
% v75 new soda netcdf
% v74 use new 2D soda netcdf file of u and v at various depths and depth-averages covering the entire soda model domain
% v73 removed PlankTOM12 and replaced with new gridded co2sys climatology woa_co2sys_2D_monthly_climatology_2005-2017_v20220604.nc
% v72 updated PlankTOM12 nc files to v60 using corrected po4 divided by 122 pers comm Rebecca Wright 5/31/2022
% v71 added OceanSODA-ETHZ biogeochem tracers as options for map overlay data
% v70 added overlay option to use ncei ocads aragonite downloaded from https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.nodc%3A0139360/html
% v69 use map3 v20220511 (159E to 305.5E, -79.5N to 0.5N)
% v68 added option for different ocean basin map extents (map3 v20220510 = South Pacific Ocean)
% v67 updated drifter.m for less verbose display of progress to increase runtime speed
% v66 save nc file of drifter trajectories and time series of exposure values extracted from map overlay data
% v65 use new monthly_anomaly function to calculate anomalies for overlay data, added stations from Bednarsek et al 2021
% v64 added pteropod stations from Bednarsek et al 2021 paper on calcification regression, updated options to animate overlay without drifters
% v63 added option to save cropped interpolated u and v and initial drifter data to a mat file
% v62 added soda.whichdrifters='nearshore', same as Feely but with added nodes along coast of Alaska and Mexico, updated drifter_v03.m
% v61 cleaned up deprecated code
% v60 added time interpolation, and separated the integration and animations into separate functions
% v59 started adding time interpolation of u, v, and overlay data
% v58 added flexible inputs for soda.nsubsteps and divded each month into 28*8 substeps for approx 3-hour time step of spatial integration
% v57 added initial drifters 5degE of the Feely transect
% v56 new map1 extent for overlay variables from PlankTOM12, same as soda map1 extent: 115E-265E, 0N-66N
% v55 added station line for 2007 West Coast cruise from Feely et al 2016 Fig 1
% v54 added anomaly option for map overlay data
% v53 added option for single initial node
% v52 added another map extent PacificNE3 to use with drifters at 153W 
% v51 moved selection of map extent latitude and longitude ranges to user input section
% v50 more efficient plotting of drifter markers on animation frames for runtime speed
% v49 added more options for cmocean colormaps
% v48 more efficient integration code for runtime speed
% v47 added option to change marker edge color to jet or any other color
% v46 save csv files of monthly soda.d01.lon and soda.d01.lat
% v45 added initial drifters along north-south trnsect lon 135W from 20N-55N and larger PacificNE2 map extent
% v44 optional markerface colors (jet, black, or white), and Osborne et al 1989 map extent option (to use with Kuroshio transport transect)
% 		and added transport transects for Kuroshio, Oyashio, and Kuroshio Extension
% v43 initialized rand with matlab default seed=0 Mersenne Twister, extended PacificNE map extent to 170W
% v42 coded for multiple drifter particles for all drifter experiments
% v41 debug animation for multiple particles per node
% v40 start adding code for multiple particles at each node
% v39 fixed stuck drifter code for reverse particle tracking, v39 is the last version with 1 particle per initial drifter node
% v37 bug fix for stuck drifter
% v36 for stuck drifters find lon lat at last ggcount with non zero u and v, and then back up 1 day to and interp new u and v at current j
% v35 for stuck drifters find lon lat at last ggcount with non zero u and v and interp new u and v at current j
% v32 use interp2 to interpolate u and v at drifter positions
% v32 read all input data from netcdf files instead of using previously run soda_v20.m
% v31 changed nc input file for overlay variables
% v29 move stuck drifters to half way back to i-1 and i-2 (25% of the way back was leaving too many still stuck)
% v27-28 added more optional overlay layers, and reads the overlay data from a netcdf file instead of requiring pre-running RECCAP script
% v26 move stuck drifters to 25% of the way back to i-1 and i-2
% v25 move stuck drifters to half way back to i-1 and i-2
% v24 added option for movmean moving average u and v
% v22-23 adjusted figure size for gif
% v21 added optional weekly saving of frames in the gif
% v20 added optional weekly saving of lon and lat positions saved from g loop
% v19 added matlab code to write an animated gif to the gif folder
% v18 added option to map drifters in the region of Japanese pteropod stations, also fixed a bug in the title year-month of animation frames
% v17 added options for quiver width, scale, color
% v16 added option to pick start and end months
% v15 added option for an east-west drifter line along latitude 33N from 140E to 118W
% v14 added option to select either PacificN with drifters along 153W from 58N to 20N, or PacificNE with drifters along Nina's sample line
% v12 added option to calculate reverse drifter tracking for exposure history track prior to specified final position
% v11 added option to overlay depth to aragonite saturation on the animation frames (requires running RECCAP_2021_3D_v40.m)
% v10 fixed code to get drifters unstuck from shore
% v08 moved k loop inside of g loop
% v07 bug fix for stuck drifters to reset them to position at the previous timestep
% v06 added substeps to divide the timestep into increments for improved accuracy
% v05 fixed bug with drifters getting stuck on shore, added east-west line of drifters to the north-south line in NE Pacific
% v04 added the option to animate any number of years between a chosen start and end year
% v03 added option to select any number of drifters, and made a line of 21 drifters from 43N to 53N at 153W
% v02 changed time period to one selected year (pickyear)
% v01 One drifter starting at 153W 53N in map1 for all soda.time

% ----------
% The following files are required to be in a subfolder named 'nc' in the working directory:
% 
% SODAv3122_2D_K100_1980_2017_v20220606.nc
% SODAv3122_2D_K200_1980_2017_v20220606.nc
% SODAv3122_2D_KT100_1980_2017_v20220606.nc
% SODAv3122_2D_KT200_1980_2017_v20220606.nc
% woa_2D_monthly_no3_o2_aou_v20220610.nc
% woa_nng_co2sys_2D_monthly_climatology_2005-2017_v20220608.nc
% woa_mobo_nng_co2sys_2D_monthly_climatology_2005-2017_v20220604.nc
% OceanSODA-ETHZ_GRaCER_v2021a_1982-2020.nc
% ncei_ocads_integrated_aragonite_v20220512.nc
%
% In addition, the following matlab scripts are required to be in the working directory,
% or in paths that are added with matlab addpath statements:
%
% cmocean2.m
% drifter_v04.m
% interp2_drifter_overlay.m
% jet2.m
% monthly_anomaly.m
% nanmax2.m
% nanmean2.m
% nanmedian2.m
% nanmin2.m
% nanstd2.m
% oceansoda2planktom.m
% soda_drifter_animation_v05.m
% tight_subplot.m
%
% The user also needs to install the following free matlab toolboxes from the followng links 
% (m_map, gif, and export_fig) and edit their addpath statements below to reflect where they are stored on the user's computer:
%
% addpath(genpath('c:/matlab/m_map/')) % functions for making maps in Mercator projection https://www.eoas.ubc.ca/~rich/map.html 
% addpath(genpath('c:/matlab/gif/')) % make animated gif https://www.mathworks.com/matlabcentral/fileexchange/63239-gif 
% addpath(genpath('c:/matlab/export_fig/')) % needed by gif https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig
%
% Outputs from soda_drifter.m will be put in the png, gif, and csv folders
%
 
% ----------
% ----------
% ----------
% START OF USER INPUTS
tic
clear all
clc
disp(['Reading input data ...'])

% Edit the following addpath statements as needed to reflect where you have stored the required 
% matlab scripts listed above for functions used in this main script
% addpath(genpath('c:/matlab/greg/'))  		% cmocean2 modified from original cmocean at https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps?s_tid=srchtitle
% addpath(genpath('c:/matlab/tight_subplot/'))  
addpath(genpath('c:/data/soda3.12.2/soda_drifter_utilities/'))  
addpath(genpath('c:/matlab/m_map/'))  		% functions for making maps in Mercator projection https://www.eoas.ubc.ca/~rich/map.html
addpath(genpath('c:/matlab/gif/'))  		% make animated gif https://www.mathworks.com/matlabcentral/fileexchange/63239-gif
addpath(genpath('c:/matlab/export_fig/'))  	% needed by gif https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig

% Edit the following default directories as needed to reflect your working directory where you store this main script
cd 'C:\data\soda3.12.2';

if not(isfolder('png'))
    mkdir('png')
end
if not(isfolder('gif'))
    mkdir('gif')
end
if not(isfolder('csv'))
    mkdir('csv')
end
if not(isfolder('nc'))
    mkdir('nc')
end

% default random seed
rng('default');	

soda.ver = 'v81';					% version number for output file names
soda.movmean_months = 1;			% number of months for moving average u and v (1= no smoothing, 3= 3-month moving average, etc)
soda.pickstartyear = 2016;			% the year when the drifter experiment starts			
soda.pickstartmonth = 12;			% the month to start the drifter experiment in the start year (e.g. 1=Jan, 2=Feb, etc)
soda.pickendyear = 2017;  			% the year when the drifter experiment ends. for one year simulation input the end year with same value as the start year 
soda.pickendmonth = 12;				% the month to end the drifter experiment in the end year (e.g. 12=Dec, 11=Nov, etc)
soda.gifdelay = 0.25;				% delay between frames for animated gif (suggest 0.25 for 12-month, 0.125 for 5 year)
soda.gifpause = 4;					% scale factor for gifdelay to pause first and last frame (suggest 4 for 12-month, 8 for 5 year)

soda.whichmap = 'Humboldt';		% Map extent to use for animated gif frames (add more options as needed in the code for making the animation) 
									% 'Humboldt' = Humboldt current region of South Pacific along west coast of South America
									% 'PacificN' = entire North Pacific from 0-65N and 128E-110W
									% 'PacificNE' = region of pteropod stations in the Northeast Pacific from 30N-61N and 165W-120W
									% 'PacificNE2' = region of 135W transect line in the Northeast Pacific from 10N-61N and 175E-110W
									% 'PacificNE3' = region of 153W transect line in the Northeast Pacific
									% 'PacificNE4' = region of Feely et al 2016 nearshore stations in the Northeast Pacific
									% 'PacificNE5' = region of Bednarsek et al 2021 pteropod stations in the Northeast Pacific
									% 'PacificNW' = region of Japanese pteropod stations from 125E-180E and 10N-60N
									% 'Osborne' = region of Osborne et al (1989) paper on trajectories for 12 months from start off SW coast of Japan
									% NOTE: any of the choices above may be used for the map extent, 
									% but if any other map extent is needed instead
									% it can be customized in the next if/elseif whichmap block of code below

soda.whichdrifters = 'shelf';		% Initial position of drifters (add more choices as needed in next if/elseif block of code):  
									% 'shelf' = slope of continental shelf along the Humboldt Current region of South Pacific
									% 'blob' = 160-140W and 43N in the blob of low TA anomaly at the start of the marine heat wave on 2014-01
									% 'Bednarsek2021' = N Pacific pteropod stations in Table 1 of Bednarsek et al 2021 https://www.frontiersin.org/articles/10.3389/fmars.2021.671497/full
									% 'Bednarsek2021b' = same as Bednarsek2021 except exclude Bering Sea
									% 'Feely2016' = offshore pteropod stations 1-13 along the US West Coast cruise 2007, Fig 1 in Feely,R.A.,etal.,ChemicalandbiologicalimpactsofoceanacidiﬁcationalongthewestcoastofNorthAmerica, Estuarine, Coastal and Shelf Science (2016), http://dx.doi.org/10.1016/j.ecss.2016.08.043
									% 'nearshore' = same as 'Feely' with added nodes along Alaska and Mexico coasts
									% 'offshore' = Parallel to 'Feely' stations, but moved 5 deg east
									% 'pteropods' = L-shaped line of pteropod stations in the NE Pacific,
									% '153W' = north south at 153W from 58N-20N,
									% '135W' = north south at 135W from 55N-20N,
									% '130W' = north south at 130W from 52N-10N,
									% '33N' = east-west at 33N from 140E to 118W (latitude of the Kuroshio Extension)
									% '50N' = east-west at 50N from 155E to 128W 
									% 'Japan' = L-shaped line surrounding Japanese pteropod stations 25N from 137.5E to 160E and 160E from 25N to 50N
									% 'Japan2' = Japanese pteropod stations 1-15
									% 'Kuroshio' = Kuroshio Current transport transect line near the SE coast of Japan in the Kuroshio Current
									% 'Oyashio' = Oyashio Current transport transect line off of East Kamchatka
									% 'KE' = Kuroshio Extension transport transect at lon 160E from 28N-409N
									% 'Spot' = Specify a single spot location below for intial drifters
									% NOTE: any of the choices above may be used for the location of initial drifter positions, 
									% but if any other drifter initial positions are needed instead
									% they can be customized in the next if/elseif whichdrifters block of code below

soda.drifter_particles = 25;		% number of drifter particles at each node to be initially spaced randomly between width and height of each node
soda.drifter_width = 0.5;			% east-west width of each drifter node (degrees)
soda.drifter_height = 0.5;			% north-south height of each drifter node (degrees)
soda.drifter_spot_lon = 360-173;	% longitude of initial drifters (degE) (only used if soda.whichdrifters='Spot')
soda.drifter_spot_lat = 60;			% latitude of initial drifters (degN) (only used if soda.whichdrifters='Spot')

soda.markeredgecolor = 'jet';		% color for the outline of the drifter markers on the animation frames 
									% 	('jet'=jet colormap, otherwise any matlab color code such as 'k'=black, 'm'=magenta, 'c'=cyan, 'y'=yellow, 'w'=white, 'b'= blue, 'r'=red, etc) 
soda.markerfacecolor = 'k';			% color for the face of the drifter markers on the animation frames 
									% 	('jet'=jet colormap, otherwise any matlab color code such as 'k'=black, 'm'=magenta, 'c'=cyan, 'y'=yellow, 'w'=white, 'b'= blue, 'r'=red, etc) 
soda.markersize = 4;				% size of the marker for the drifter trajectory animation (points, e.g. matlab default size is 6 points)

soda.whichdic = 'nng';				% 'nng' or 'mobo' (suggested default is 'nng') - Select which overlay netcdf file of climatology to use, either 'mobo' for MOBO-DIC based, or 'nng' for TCO2-NNGv2 based 

soda.overlay = 'omara_K100';		% Overlay data layer to show as background on the maps of drifter trajectories:
									% 'none' = do not overlay a data layer on the animation map frames
									% 'omara_satdep' = depth to aragonite saturation horizon
									% 'omcal_satdep' = depth to calcite saturation horizon
									% 'omara_KT' = Surface 0m omega aragonite
									% 'omara_K100' = Depth 100m omega aragonite
									% 'omara_KT100' = Depth-integrated 0-100m omega aragonite
									% 'omara_K200' = Depth 200m omega aragonite
									% 'omara_KT200' = Depth-integrated 0-200m omega aragonite
									% 'omcal_KT' = Surface 0m omega calcite
									% 'omcal_K100' = Depth 100m omega calcite
									% 'omcal_KT100' = Depth-integrated 0-100m omega calcite
									% 'omcal_K200' = Depth 200m omega calcite
									% 'omcal_KT200' = Depth-integrated 0-200m omega calcite
									% 'phtot_KT' = Surface 0m pH (total scale)
									% 'phtot_K100' = Depth 100m pH (total scale)
									% 'phtot_KT100' = Depth-integrated 0-100m pH (total scale)
									% 'phtot_K200' = Depth 200m pH (total scale)
									% 'phtot_KT200' = Depth-integrated 0-200m pH (total scale)
									% 'talk_KT100' = Depth-integrated 0-100m total alkalinity
									% 'dic_KT100' = Depth-integrated 0-100m DIC
									% 'pco2_KT100' = Depth-integrated 0-100m pCO2
									% 'revelle_KT100' = Depth-integrated 0-100m Revelle Factor
									% 'temp_KT100' = Depth-integrated 0-100m temperature
									% 'sal_KT100' = Depth-integrated 0-100m salinity
									% 'po4_KT100' = Depth-integrated 0-100m total inorganic phosphorus
									% 'si_KT100' = Depth-integrated 0-100m total inorganic silica
									% 'no3_KT100' = Depth-integrated 0-100m nitrate
									% 'o2_KT100' = Depth-integrated 0-100m oxygen
									% 'aou_KT100' = Depth-integrated 0-100m Apparent Oxygen Utilization
									% 'omara' - OceanSODA surface aragonite saturation state omega
									% 'omcal' - OceanSODA surface calcite saturation state omega
									% 'temp' - OceanSODA surface temperature
									% 'sal' - OceanSODA surface salinity
									% 'dic' - OceanSODA surface dissolved inorganic carbon
									% 'talk' - OceanSODA surface total alkalinity
									% 'phtot' - OceanSODA surface pH (total scale)
									% 'spco2' - OceanSODA surface partial pressure of carbon dioxide
									% 'revelle' - OceanSODA surface Revelle Factor
									% 'hco3' - OceanSODA surface HCO3
									% 'co3' - OceanSODA surface CO3
									% 'co2' - OceanSODA surface CO2
									% 'fgco2' - OceanSODA Sea-air flux of CO2 calculated using the SeaFlux dataset 
									% 'ncei_omara_satdep' = Jiang et al (2015) depth to aragonite saturation horizon (gridded data from NCEI OCADS https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.nodc%3A0139360/html)
									% 'ncei_omara_KT100' = Jiang et al (2015) 0-100m depth-integrated aragonite saturation (gridded data from NCEI OCADS https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.nodc%3A0139360/html)
									% NOTE: overlay data must be lon x lat x time and have same time origin of 1980 and be a continuous sequence of monthly data

soda.velocity = 'K100';			% Specify which depth or depth-integrated average to use for u and v velocity vectors
									% 'KT100' = Drifters at depth-integrated average velocity from 0-100m
									% 'KT200' = Drifters at depth-integrated average velocity from 0-200m
									% 'K100' = Drifters at depth of 100m
									% 'K200' = Drifters at depth of 200m

soda.quiverwidth = .1;			 	% line width of quiver arrows (suggested default 0.1)
soda.quiverscale = 5;			 	% scale factor for quiver length (suggested default 5)
soda.quivercolor = 'k';			 	% color for quiver arrows 

soda.nsubdayspermonth = 30;			% intermediate increment to divide each month into approximately days
soda.nsubhoursperday = 8;			% intermediate increment to divide each day into sub-daily increments, (e.g. 4=6 hours, 8=3 hours, 12=2 hours, 24=1 hour, etc.)
soda.skipdays = 7;					% number of days to skip between animation frames (for example soda.skipdays=7 shows weekly animation frames)

soda.forward_tracking = false;		% drifter tracking forward in time (true), or reverse particle tracking for exposure history (false)
soda.makegif = true;				% true or false to make animated gif
soda.overlay_quiver = true;			% overlay the current vector quiver arrows on the animation frames (used if soda.makegif=true)
soda.overlay_drifter = true;		% show the drifters animation frames (used if soda.makegif=true)
soda.integrate = true;				% true or false to calculate drifter trajectories 
soda.save_drifter_track = false;	% true or false to save drifter trajectory lon, lat, and interpolated overlay data as netcdf file
soda.plot_drifter_track = true;		% true or false to plot map, time-series, and boxplot of drifter tracks and interpolated overlay values along drifter tracks 
soda.overlay_anomaly = false;		% use the monthly anomaly compared with the long-term monthly averages for the selected map overlay (only appropriate for OceanSODA)
soda.savemat = false;				% true or false to save the interpolated cropped u and v and intial drifter location data to a mat file for external calculation of drifter tracks
soda.extraplots = false;			% true or false to output extra plots with pH and aragonite thresholds and cropped regions

% END OF USER INPUTS
% ----------
% ----------
% ----------

% calculated values using the inputs above
soda.picknyears = soda.pickendyear-soda.pickstartyear+1;	% intermediate calc to find the number of months for duration of the experiment (picknmonths)
soda.picknmonths = (soda.picknyears-2) * 12 + (12-soda.pickstartmonth+1) + soda.pickendmonth;
soda.nsubsteps = soda.nsubdayspermonth * soda.nsubhoursperday;		% number of increments to divide each time step (nsubsteps=30*24=720 is approx 1 hour substep increment

% Edit the following as needed to change the map extent of the animation frames 
% select the map extent for the animation frames (add more options as needed)
if strcmp(soda.whichmap,'PacificN') 
	soda.maplat = [0 65];					% south and north latitudes of map extent (degN)
	soda.maplon = [(115-360) -100];			% west and east longitudes of map extent (-degW)
	soda.croplat = [];
	soda.croplon = [];
elseif strcmp(soda.whichmap,'PacificNE') 
	soda.maplat = [30 61];					% south and north latitudes of map extent (degN)
	soda.maplon = [-175 -120];				% west and east longitudes of map extent (-degW)
	soda.croplat = [];
	soda.croplon = [];
elseif strcmp(soda.whichmap,'PacificNE2') 
	soda.maplat = [10 61];					% south and north latitudes of map extent (degN)
	soda.maplon = [-185 -110];				% west and east longitudes of map extent (-degW)
	soda.croplat = [];
	soda.croplon = [];
elseif strcmp(soda.whichmap,'PacificNE3') 
	soda.maplat = [10 61];					% south and north latitudes of map extent (degN)
	soda.maplon = [-195 -120];				% west and east longitudes of map extent (-degW)
	soda.croplat = [];
	soda.croplon = [];
elseif strcmp(soda.whichmap,'PacificNE4') 
	soda.maplat = [0 61];					% south and north latitudes of map extent (degN)
	soda.maplon = [-190 -94.5];				% west and east longitudes of map extent (-degW)
	soda.croplat = [];
	soda.croplon = [];
elseif strcmp(soda.whichmap,'PacificNE5') 
	soda.maplat = [20 65];					% south and north latitudes of map extent (degN)
	% soda.maplon = [-190 -94.5];				% west and east longitudes of map extent (-degW)
	soda.maplon = [-190-10 -94.5-20];				% west and east longitudes of map extent (-degW)
	soda.croplat = [];
	soda.croplon = [];
elseif strcmp(soda.whichmap,'PacificNW') 
	soda.maplat = [10 60];					% south and north latitudes of map extent (degN)
	soda.maplon = [(115-360) -160];		% west and east longitudes of map extent (-degW)
	soda.croplat = [];
	soda.croplon = [];
elseif strcmp(soda.whichmap,'Osborne') 
	soda.maplat = [10 60];					% south and north latitudes of map extent (degN)
	soda.maplon = [(120-360) -155];		% west and east longitudes of map extent (-degW)
	soda.croplat = [];
	soda.croplon = [];
elseif strcmp(soda.whichmap,'Humboldt') 
	soda.maplat = [-61 0];					% south and north latitudes of map extent (degN)
	soda.maplon = [-160 -55];		% west and east longitudes of map extent (-degW)
	soda.croplat = [-45 -5];					% south and north latitudes of map extent (degN)
	soda.croplon = [-100 -65];		% west and east longitudes of map extent (-degW)
	soda.croporient = 'portrait';	% 'portrait' or 'landscape' for inset crop map (default 'landscape' if not specified)
end

% Edit the following as needed to change the initial/final position of the drifter node centers 
% or add more expermients with drifter nodes in different locations. Each drifter node will have 
% the number of soda.drifter_particles specified above randomly spaced within the specifed width and height above

% v60
soda.drifter_lon0 = [];
soda.drifter_lat0 = [];
if strcmp(soda.whichdrifters,'Spot')	% User-specified spot location
	soda.ndrift = 1;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	soda.drifter_lon0 = soda.drifter_spot_lon;					% degE
	soda.drifter_lat0 = soda.drifter_spot_lat;					% degN
	soda.drifter.sta = [];
elseif strcmp(soda.whichdrifters,'shelf') 	% slope of continental shelf along the Humboldt Current region of South Pacific
	soda.ndrift = 24;					% number of drifters
	% initialize vectors for drifter lat and lon
	soda.drifter.lon0 = [280.4	281.7	282.7	284.4	287.3	288.5	288.6	288.4	288.4	287.8	287.5	287.2	286.9	285.8	285.3	284.9	284.4	284.2	283.9	283.4	283	283.6	284.2	288.3];		% degE
	soda.drifter.lat0 = [-10	-12	-14	-16	-18	-20	-22	-24	-26	-28	-30	-32	-34	-36	-38	-40	-42	-44	-46	-48	-50	-52	-54	-56];			% degN
	soda.drifter.sta = {['10',char(176),'S'],	['12',char(176),'S'],	['14',char(176),'S'],	['16',char(176),'S'],	['18',char(176),'S'],	['20',char(176),'S'],	['22',char(176),'S'],	['24',char(176),'S'],	['26',char(176),'S'],	['28',char(176),'S'],	['30',char(176),'S'],	['32',char(176),'S'],	['34',char(176),'S'],	['36',char(176),'S'],	['38',char(176),'S'],	['40',char(176),'S'],	['42',char(176),'S'],	['44',char(176),'S'],	['46',char(176),'S'],	['48',char(176),'S'],	['50',char(176),'S'],	['52',char(176),'S'],	['54',char(176),'S'],	['56',char(176),'S']};
elseif strcmp(soda.whichdrifters,'blob') 	% 160-140W and 43N in the blob of low TA anomaly at the start of the marine heat wave on 2014-01
	soda.ndrift = 31;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	soda.drifter.lon0 = [(360-160):1:(360-130)];		% degE
	soda.drifter.lat0(1:soda.ndrift) = 45;			% degN
	soda.drifter.sta = [];
elseif strcmp(soda.whichdrifters,'Bednarsek2021') 	% N Pacific pteropod stations in Table 1 of Bednarsek et al 2021 https://www.frontiersin.org/articles/10.3389/fmars.2021.671497/full
	soda.ndrift = 20;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	soda.drifter.lon0 = [207.999	208.401	207.402	208	208	208.402	208.403	207.496	206.889	208.933	211.966	216.566	219.703	222.817	196.399	194.913	194.854	194.084	194.798	186.356];		% degE
	soda.drifter.lat0 = [41.301	46.3	48.001	49.23	51	52.599	54.396	55.466	56.144	54.043	54.33	56.218	55.581	56.33	57.52	55.056	54.426	54.36	54.711	60.573];			% degN
	soda.drifter.sta = ["151",	"161",	"164"	"167"	"170"	"174"	"178"	"182"	"187"	"191"	"194"	"198"	"201"	"204"	"24"	"45"	"50"	"53"	"56"	"45H2"];
elseif strcmp(soda.whichdrifters,'Bednarsek2021b') 	% same as Bednarsek et al 2021 except exclude Bering Sea except for 45H2
	soda.ndrift = 15;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	soda.drifter.lon0 = [207.999	208.401	207.402	208	208	208.402	208.403	207.496	206.889	208.933	211.966	216.566	219.703	222.817	186.356];		% degE
	soda.drifter.lat0 = [41.301	46.3	48.001	49.23	51	52.599	54.396	55.466	56.144	54.043	54.33	56.218	55.581	56.33	60.573];			% degN
	soda.drifter.sta = ["151",	"161",	"164"	"167"	"170"	"174"	"178"	"182"	"187"	"191"	"194"	"198"	"201"	"204"	"45H2"];
elseif strcmp(soda.whichdrifters,'Feely2016') 	% offshore stations at 2007 West Coast cruise stations in Fig 1 by Feely,R.A.,etal.,ChemicalandbiologicalimpactsofoceanacidiﬁcationalongthewestcoastofNorthAmerica, Estuarine, Coastal and Shelf Science (2016), http://dx.doi.org/10.1016/j.ecss.2016.08.043
	soda.ndrift = 13;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	soda.drifter.lon0 = [227.135	230.517	232.456	233.069	232.771	233.566	234.892	236.135	237.412	238.953	240.478	243.412	245.218];		% degE
	soda.drifter.lat0 = [50.669	48.277	46.034	44.606	41.383	39.29	37.662	35.884	34.173	32.578	30.319	27.196	25.053];			% degN
	soda.drifter.sta = {['50.7',char(176),'N']	['48.3',char(176),'N']	['46.0',char(176),'N']	['44.6',char(176),'N']	['41.4',char(176),'N']	['39.3',char(176),'N']	['37.7',char(176),'N']	['35.9',char(176),'N']	['34.2',char(176),'N']	['32.6',char(176),'N']	['30.3',char(176),'N']	['27.2',char(176),'N']	['25.1',char(176),'N']};
elseif strcmp(soda.whichdrifters,'nearshore') 	% same as Feely with added nodes along coast of Alaska and Mexico
	soda.ndrift = 24;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	soda.drifter.lon0 = [196.55	202.78	208.49	212.66	217.42	220.33	222.14	224.15	227.135	230.517	232.456	233.069	232.771	233.566	234.892	236.135	237.412	238.953	240.478	243.412	245.218	248.64	253.56	257.55];		% degE
	soda.drifter.lat0 = [52.91	54.07	55.03	57.24	58.62	57.43	55.72	53.67	50.669	48.277	46.034	44.606	41.383	39.29	37.662	35.884	34.173	32.578	30.319	27.196	25.053	22.29	19.33	17.09];			% degN
	soda.drifter.sta = {['52.9',char(176),'N']	['54.1',char(176),'N']	['55.0',char(176),'N']	['57.2',char(176),'N']	['58.6',char(176),'N']	['57.4',char(176),'N']	['55.7',char(176),'N']	['53.7',char(176),'N']	['50.7',char(176),'N']	['48.3',char(176),'N']	['46.0',char(176),'N']	['44.6',char(176),'N']	['41.4',char(176),'N']	['39.3',char(176),'N']	['37.7',char(176),'N']	['35.9',char(176),'N']	['34.2',char(176),'N']	['32.6',char(176),'N']	['30.3',char(176),'N']	['27.2',char(176),'N']	['25.1',char(176),'N']	['22.3',char(176),'N']	['19.3',char(176),'N']	['17.1',char(176),'N']};
elseif strcmp(soda.whichdrifters,'offshore') 	% offshore stations at 2007 West Coast cruise stations in Fig 1 by Feely,R.A.,etal.,ChemicalandbiologicalimpactsofoceanacidiﬁcationalongthewestcoastofNorthAmerica, Estuarine, Coastal and Shelf Science (2016), http://dx.doi.org/10.1016/j.ecss.2016.08.043
	soda.ndrift = 14;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	soda.drifter.lon0 = [220.6	222.135	225.517	227.456	228.069	227.771	228.566	229.892	231.135	232.412	233.953	235.478	238.412	240.218];		% degE
	soda.drifter.lat0 = [53.7	50.669	48.277	46.034	44.606	41.383	39.29	37.662	35.884	34.173	32.578	30.319	27.196	25.053];			% degN
	soda.drifter.sta = {['53.7',char(176),'N']	['50.7',char(176),'N']	['48.3',char(176),'N']	['46.0',char(176),'N']	['44.6',char(176),'N']	['41.4',char(176),'N']	['39.3',char(176),'N']	['37.7',char(176),'N']	['35.9',char(176),'N']	['34.2',char(176),'N']	['32.6',char(176),'N']	['30.3',char(176),'N']	['27.2',char(176),'N']	['25.1',char(176),'N']};
elseif strcmp(soda.whichdrifters,'pteropods') 
	soda.ndrift = 51;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	% initial position of drifters
	soda.drifter.lat0(1:36) = 50;		% degN
	soda.drifter.lon0(1) = 360-135;	% degE
	for i = 2:36
		soda.drifter.lon0(i) = soda.drifter.lon0(i-1) - 0.5;
	end
	soda.drifter.lon0(37:51) = 360-153;
	soda.drifter.lat0(37) = 50;
	for i = 38:51
		soda.drifter.lat0(i) = soda.drifter.lat0(i-1) - 0.5;
	end
	soda.drifter.sta = [];
elseif strcmp(soda.whichdrifters,'153W') 
	soda.ndrift = 75;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	% initial position of drifters
	soda.drifter.lon0(1:soda.ndrift) = 360-153;	% degE
	soda.drifter.lat0(1) = 57;						% degN
	for i = 2:soda.ndrift
		soda.drifter.lat0(i) = soda.drifter.lat0(i-1) - 0.5;
	end
	soda.drifter.sta = [];
elseif strcmp(soda.whichdrifters,'135W') 
	soda.ndrift = 71;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	% initial position of drifters
	soda.drifter.lon0(1:soda.ndrift) = 360-135;		% degE
	soda.drifter.lat0(1) = 55;							% degN
	for i = 2:soda.ndrift
		soda.drifter.lat0(i) = soda.drifter.lat0(i-1) - 0.5;
	end
	soda.drifter.sta = [];
elseif strcmp(soda.whichdrifters,'130W') 
	soda.ndrift = 85;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	% initial position of drifters
	soda.drifter.lon0(1:soda.ndrift) = 360-130;		% degE
	soda.drifter.lat0(1) = 52;							% degN
	for i = 2:soda.ndrift
		soda.drifter.lat0(i) = soda.drifter.lat0(i-1) - 0.5;
	end
	soda.drifter.sta = [];
elseif strcmp(soda.whichdrifters,'33N') 
	soda.ndrift = 205;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	% initial position of drifters
	soda.drifter.lat0(1:soda.ndrift) = 33;			% degN
	soda.drifter.lon0(1) = 360-118;				% degE
	for i = 2:soda.ndrift
		soda.drifter.lon0(i) = soda.drifter.lon0(i-1) - 0.5;
	end
	soda.drifter.sta = [];
elseif strcmp(soda.whichdrifters,'50N') 
	soda.ndrift = 155;			% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	% initial position of drifters
	soda.drifter.lat0(1:soda.ndrift) = 50;		% degN
	soda.drifter.lon0(1) = 360-128;			% degE
	for i = 2:soda.ndrift
		soda.drifter.lon0(i) = soda.drifter.lon0(i-1) - 0.5;
	end
	soda.drifter.sta = [];
elseif strcmp(soda.whichdrifters,'Japan') 
	soda.ndrift = 96;					% number of drifters
	% initialize vectors for initial drifter nodes
	% v60
	% initial position of drifters
	soda.drifter.lon0(1:51) = 160;		% degE
	soda.drifter.lat0(1) = 50;			% degN
	for i = 2:51
		soda.drifter.lat0(i) = soda.drifter.lat0(i-1) - 0.5;
	end
	soda.drifter.lat0(52:96) = 25;		% degN
	soda.drifter.lon0(52) = 159.5;		% degE
	for i = 53:96
		soda.drifter.lon0(i) = soda.drifter.lon0(i-1) - 0.5;
	end
	soda.drifter.sta = [];
elseif strcmp(soda.whichdrifters,'Japan2') 	% Japanese pteropod stations 1-15
	soda.ndrift = 15;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	soda.drifter.lon0 = [137.5	137.5	137.5	145	144	141	145	150	157.5	159	160	167.5	155	147.5	142];		% degE
	soda.drifter.lat0 = [32.5	30	25	25	32	35	35	35	35	42.5	47	50	44.5	41	42];			% degN
	soda.drifter.sta = ["1"	"2"	"3"	"4"	"5"	"6"	"7"	"8"	"9"	"10"	"K2"	"12"	"13"	"14"	"15"];
elseif strcmp(soda.whichdrifters,'Kuroshio')	% Kuroshio Current transport transect off the SW coast of Japan
	soda.ndrift = 7;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	% initial position of drifters
	soda.drifter.lon0 = [141	141.375	141.75	142.125	142.5 142.875 143.25];		% degE
	soda.drifter.lat0 = [36	35.5	35	34.5	34 33.5 33];						% degN
	soda.drifter.sta = [];
elseif strcmp(soda.whichdrifters,'Oyashio') 	% Oyashio current transport transect at East Kamchatka
	soda.ndrift = 6;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	% initial position of drifters
	soda.drifter.lon0 = [160	160.625	161.25	161.875	162.5	163.125];		% degE
	soda.drifter.lat0 = [54	53.6875	53.375	53.0625	52.75	52.4375];						% degN
	soda.drifter.sta = [];
elseif strcmp(soda.whichdrifters,'KE')		% Kuroshio Extension transport transect 
	soda.ndrift = 43;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	% initial position of drifters
	soda.drifter.lon0(1:soda.ndrift) = 160;			% degE
	soda.drifter.lat0 = [49:-0.5:28];	% degN
	soda.drifter.sta = [];
end

% map sub-region, and file names for netcdf input files of soda lon, lat, time, u, v, and overlay data 
soda.mapx.lon_idx = [];
soda.mapx.lat_idx = [];
if strcmp(soda.whichmap, 'Humboldt')
	% map3 - South Pacific index values for sub-region of interest
	soda.mapx.lon_idx = ((160*2-1):1:(306*2-1))';	% 159.5E to 305.5E at 0.5 deg increments for SODA3.12.2 u and v arrays
	soda.mapx.lat_idx = (1:1:151)';					% -74.5N to 0.5N at 0.5 deg increments for SODA3.12.2 u and v arrays
	ethz.mapx.lon_idx = (160:1:306)';				% 159.5E to 305.5 degE at 1 deg increments for OceanSODA-ETHZ overlay arrays
	ethz.mapx.lat_idx = (11:1:91)';					% -79.5 to 0.5 degN at 1 deg increments for OceanSODA-ETHZ overlay arrays	
else
	% map1 - North Pacific index values for sub-region of interest
	soda.mapx.lon_idx = (115*2:1:265*2)';			% 114.5 to 264.5 degE at 0.5 deg increments for SODA3.12.2 u and v arrays
	soda.mapx.lat_idx = 75*2 + (0:1:66*2)';			% -0.5 to 65.5 degN at 0.5 deg increments for SODA3.12.2 u and v arrays
	ethz.mapx.lon_idx = (115:1:265)';				% 114.5 to 264.5 degE at 1 deg increments for OceanSODA-ETHZ overlay arrays
	ethz.mapx.lat_idx = 90 + (0:1:66)';				% -0.5 to 65.5 degN at 1 deg increments for OceanSODA-ETHZ overlay arrays
end
if strcmp(soda.whichdic, 'mobo')
	soda.overlay_fn = [pwd, '\nc\woa_mobo_nng_co2sys_2D_monthly_climatology_2005-2017_v20220604.nc'];		% edit as needed to correspond with the path to this nc file
elseif strcmp(soda.whichdic, 'nng')
	soda.overlay_fn = [pwd, '\nc\woa_nng_co2sys_2D_monthly_climatology_2005-2017_v20220608.nc'];		% edit as needed to correspond with the path to this nc file
else
	error('Select either mobo or nng for which climatology to use for the overlay data');
end
soda.overlay_fn2 = [pwd, '\nc\ncei_ocads_integrated_aragonite_v20220512.nc'];		% edit as needed to correspond with the path to this nc file
soda.overlay_fn3 = [pwd, '\nc\OceanSODA-ETHZ_GRaCER_v2021a_1982-2020.nc'];		% edit as needed to correspond with the path to this nc file
soda.overlay_fn4 = [pwd, '\nc\woa_2D_monthly_no3_o2_aou_v20220610.nc'];		% edit as needed to correspond with the path to this nc file

% - - -
% read velocity data from SODA netcdf file

if strcmp(soda.velocity,'KT100')
	soda.fn = [pwd '\nc\SODAv3122_2D_KT100_1980_2017_v20220606.nc'];	% edit as needed to correspond with the path to this nc file
	soda.plot_title = 'Drifters at depth of 0-100m';			
	strdep1 = ' at 0-100m';
elseif strcmp(soda.velocity,'KT200')
	soda.fn = [pwd '\nc\SODAv3122_2D_KT200_1980_2017_v20220606.nc'];	% edit as needed to correspond with the path to this nc file
	soda.plot_title = 'Drifters at depth of 0-200m';			
	strdep1 = ' at 0-200m';
elseif strcmp(soda.velocity,'K100')
	soda.fn = [pwd '\nc\SODAv3122_2D_K100_1980_2017_v20220606.nc'];	% edit as needed to correspond with the path to this nc file
	soda.plot_title = 'Drifters at depth of 100m';			
	strdep1 = ' at 100m';
elseif strcmp(soda.velocity,'K200')
	soda.fn = [pwd '\nc\SODAv3122_2D_K200_1980_2017_v20220606.nc'];	% edit as needed to correspond with the path to this nc file
	soda.plot_title = 'Drifters at depth of 200m';			
	strdep1 = ' at 200m';
else
	error('Enter either KT100, KT200, K100, or K200 for velocity');
end

% read all lon and lat from entire soda model
soda.lon = ncread(soda.fn,'lon');			% 0.5:0.5:359.5 degE
soda.lat = ncread(soda.fn,'lat');			% -74.5:0.5:89.5 degN

% extract the lon and lat for the sub-region of interest
soda.mapx.lon = soda.lon(soda.mapx.lon_idx);	% extract sub-region of interest		
soda.mapx.lat = soda.lat(soda.mapx.lat_idx);	% extract sub-region of interest	
soda.mapx.lonW = soda.mapx.lon - 360;			% extract sub-region of interest

soda.time = ncread(soda.fn,'time');
soda.datenum = soda.time + datenum(datetime(1980,1,1,0,0,0));
soda.datetime = datetime(soda.datenum,'ConvertFrom','datenum');

% read all u and v from the the entire soda model
soda.u_name = ['u_',soda.velocity];
soda.v_name = ['v_',soda.velocity];
soda.u = [];
soda.v = [];
soda.u = ncread(soda.fn,soda.u_name);
soda.v = ncread(soda.fn,soda.v_name);

% extract the u and v for the sub-region of interest
soda.mapx.u = [];
soda.mapx.v = [];
soda.mapx.u = soda.u(soda.mapx.lon_idx,soda.mapx.lat_idx,:);
soda.mapx.v = soda.v(soda.mapx.lon_idx,soda.mapx.lat_idx,:);

% - - -

% start and end year of the soda netcdf file data
soda.startyear = 1980;	% soda.startyear 1980 was used for the netcdf file above for vertically integrated u and v and must be 1980 to correctly align the map overlay the outputs from PlankTOM12 specified below
soda.endyear = 2017;	% 2017 was used for the netcdf above for vertically integrated u and v
soda.nyears = soda.endyear-soda.startyear+1;
soda.nmonths = soda.nyears * 12;
soda.yearnum = [soda.startyear:soda.nyears/soda.nmonths:(soda.endyear+1)]';   
soda.yearnum = soda.yearnum(1:soda.nmonths);									% decimal years for x-axis for plots

% year and month of each row for soda_drifter animation labels, calculated from the inputs above
soda.year(1:numel(soda.time),1) = NaN;
soda.month(1:numel(soda.time),1) = NaN;
imonth = 0;
iyear = soda.startyear;
for i = 1:numel(soda.time)
	imonth = imonth + 1;
	soda.year(i,1) = iyear;
	soda.month(i,1) = imonth;
	if imonth == 12
		imonth = 0;
		iyear = iyear +1;
	end
end

% calculate movmean u and v moving average
if soda.movmean_months > 1
	for i = 1:numel(soda.mapx.lon)
		for j = 1:numel(soda.mapx.lat)
			soda.mapx.u(i,j,:) = movmean(soda.mapx.u(i,j,:),soda.movmean_months);
			soda.mapx.v(i,j,:) = movmean(soda.mapx.v(i,j,:),soda.movmean_months);
		end
	end
end

% map overlay data
if ~strcmp(soda.overlay,'none') ...
	& ~strcmp(soda.overlay,'ncei_omara_satdep') ...
	& ~strcmp(soda.overlay,'ncei_omara_KT100') ...
	& ~strcmp(soda.overlay,'omara') ...
	& ~strcmp(soda.overlay,'omcal') ...
	& ~strcmp(soda.overlay,'temp') ...
	& ~strcmp(soda.overlay,'sal') ...
	& ~strcmp(soda.overlay,'dic') ...
	& ~strcmp(soda.overlay,'talk') ...
	& ~strcmp(soda.overlay,'phtot') ...
	& ~strcmp(soda.overlay,'spco2') ...
	& ~strcmp(soda.overlay,'revelle') ...
	& ~strcmp(soda.overlay,'hco3') ...
	& ~strcmp(soda.overlay,'co3') ...
	& ~strcmp(soda.overlay,'co2') ...
	& ~strcmp(soda.overlay,'fgco2')
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to -359.5:1:0.5 degE  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');		% -89.5:1:89.5 degN  
	if strcmp(soda.overlay,'no3_KT100') | strcmp(soda.overlay,'o2_KT100') | strcmp(soda.overlay,'aou_KT100')
		soda.overlay_data = ncread(soda.overlay_fn4,soda.overlay);  
		soda.overlay_nc_units = ncreadatt(soda.overlay_fn4,soda.overlay,'units');  
		soda.overlay_nc_long_name = ncreadatt(soda.overlay_fn4,soda.overlay,'long_name');  
	else
		soda.overlay_data = ncread(soda.overlay_fn,soda.overlay);  
		soda.overlay_nc_units = ncreadatt(soda.overlay_fn,soda.overlay,'units');  
		soda.overlay_nc_long_name = ncreadatt(soda.overlay_fn,soda.overlay,'long_name');  
	end
	% soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_data = repmat(soda.overlay_data,1,1,soda.nyears);
elseif strcmp(soda.overlay,'omara') ...
	| strcmp(soda.overlay,'omcal') ...
	| strcmp(soda.overlay,'temp') ...
	| strcmp(soda.overlay,'sal') ...
	| strcmp(soda.overlay,'dic') ...
	| strcmp(soda.overlay,'talk') ...
	| strcmp(soda.overlay,'phtot') ...
	| strcmp(soda.overlay,'spco2') ...
	| strcmp(soda.overlay,'revelle') ...
	| strcmp(soda.overlay,'hco3') ...
	| strcmp(soda.overlay,'co3') ...
	| strcmp(soda.overlay,'co2') ...
	| strcmp(soda.overlay,'fgco2')

	% - - -
	% OceanSODA-ETHZ

	ethz.lon_raw = ncread(soda.overlay_fn3,'lon');		% degE (-179.5:1:179.5)    
	ethz.lat_raw = ncread(soda.overlay_fn3,'lat');  	% degN (-89.5:1:89.5)
	ethz.time_raw = ncread(soda.overlay_fn3,'time');	% 'days since 1982-01-01'

	% read overlay data from OceanSODA-ETHZ nc file 
	if strcmp(soda.overlay,'omara')
		ethz.varname = 'omega_ar';				% aragonite saturation state
		soda.overlay_nc_units = 'omega_ar';  
		soda.overlay_nc_long_name = 'Aragonite saturation state omega';  
	elseif strcmp(soda.overlay,'omcal')
		ethz.varname = 'omega_ca';				% aragonite saturation state
		soda.overlay_nc_units = 'omega_ca';  
		soda.overlay_nc_long_name = 'Calcite saturation state omega';  
	elseif strcmp(soda.overlay,'temp')
		ethz.varname = 'temperature';			% degC
		soda.overlay_nc_units = 'degC';  
		soda.overlay_nc_long_name = 'Temperature';  
	elseif strcmp(soda.overlay,'dic')
		ethz.varname = 'dic';					% umol/kg
		soda.overlay_nc_units = 'umol/kg';  
		soda.overlay_nc_long_name = 'Dissolved Inorganic Carbon';  
	elseif strcmp(soda.overlay,'talk')
		ethz.varname = 'talk';					% umol/kg
		soda.overlay_nc_units = 'umol/kg';  
		soda.overlay_nc_long_name = 'Total Alkalinity';  
	elseif strcmp(soda.overlay,'phtot')
		ethz.varname = 'ph_total';	
		soda.overlay_nc_units = 'pH (total scale)';  
		soda.overlay_nc_long_name = 'pH (total scale)';  
	elseif strcmp(soda.overlay,'revelle')
		ethz.varname = 'revelle_factor';	
		soda.overlay_nc_units = 'Revelle Factor';  
		soda.overlay_nc_long_name = 'Revelle Factor';  
	elseif strcmp(soda.overlay,'spco2')
		ethz.varname = 'spco2';					% uatm
		soda.overlay_nc_units = 'uatm';  
		soda.overlay_nc_long_name = 'surface_partial_pressure_of_carbon_dioxide_in_sea_water';  
	elseif strcmp(soda.overlay,'hco3')
		ethz.varname = 'hco3';
		soda.overlay_nc_units = 'umol/kg';  
		soda.overlay_nc_long_name = 'HCO3';  
	elseif strcmp(soda.overlay,'co3')
		ethz.varname = 'co3';	
		soda.overlay_nc_units = 'umol/kg';  
		soda.overlay_nc_long_name = 'CO3';  
	elseif strcmp(soda.overlay,'co2')
		ethz.varname = 'co2';	
		soda.overlay_nc_units = 'umol/kg';  
		soda.overlay_nc_long_name = 'CO2';  
	elseif strcmp(soda.overlay,'sal')
		ethz.varname = 'salinity';	
		soda.overlay_nc_units = 'psu';  
		soda.overlay_nc_long_name = 'Salinity';  
	elseif strcmp(soda.overlay,'fgco2')
		ethz.varname = 'fgco2';	
		soda.overlay_nc_units = 'mol/m2/yr';  
		soda.overlay_nc_long_name = 'Sea-air flux calculated using the SeaFlux dataset.';  
	end
	ethz.overlay_raw = ncread(soda.overlay_fn3,ethz.varname);
	ethz.overlay_raw(ethz.overlay_raw==0)=NaN;			% set non-water cells to nan

	% new time, datenum, and datetime starting in 1980 to synch with PlankTOM12 and SODA start time
	ethz.datenum = datenum(1980,1:(24+468),1)';		% vector of months from 1980-2020
	ethz.datetime = datetime(ethz.datenum,'ConvertFrom','datenum');
	ethz.time = [];
	ethz.time(1:numel(ethz.datenum),1)=NaN;
	ethz.time = ethz.datenum - ethz.datenum(1,1);

	% new array for biogeochem overlay tracer, convert to ethz.lon 0.5:1:359.5 degE format
	% with time starting 1980,1,1 to synchronize with SODA u and v arrays
	ethz.overlay(1:numel(ethz.lon_raw),1:numel(ethz.lat_raw),1:numel(ethz.time)) = NaN;
	[ethz.lon, ethz.lat, ethz.overlay(:,:,25:numel(ethz.time))] = oceansoda2planktom(ethz.lon_raw,ethz.lat_raw,ethz.overlay_raw);

	% extract overlay data for mapx subset region
	soda.overlay_lon = ethz.lon(ethz.mapx.lon_idx) - 360; 		% extract and convert to -359.5:1:0.5 degE
	soda.overlay_lat = ethz.lat(ethz.mapx.lat_idx);				% extract degN
	soda.overlay_data = ethz.overlay(ethz.mapx.lon_idx, ethz.mapx.lat_idx, :);

end
if strcmp(soda.overlay,'omara')
	soda.overlay_name = '\Omega_A';
	soda.overlay_pivot = 1;
	soda.overlay_min = 0;
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_min = min(soda.overlay_pivot,soda.overlay_min);
	soda.overlay_max = max(soda.overlay_pivot,soda.overlay_max);
	soda.overlay_cm = '-delta';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'omcal')
	soda.overlay_name = '\Omega_C';
	soda.overlay_pivot = 1;
	soda.overlay_min = 0;
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_min = min(soda.overlay_pivot,soda.overlay_min);
	soda.overlay_max = max(soda.overlay_pivot,soda.overlay_max);
	soda.overlay_cm = '-delta';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'phtot')
	soda.overlay_name = 'pH';
	soda.overlay_pivot = 7.6;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_min = min(soda.overlay_pivot,soda.overlay_min);
	soda.overlay_max = max(soda.overlay_pivot,soda.overlay_max);
	soda.overlay_cm = '-delta';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'talk')
	soda.overlay_name = 'TA (\mumol kg^-^1)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'dic')
	soda.overlay_name = 'DIC (\mumol kg^-^1)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'temp') 
	soda.overlay_name = 'Temperature (\circC)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = 'thermal';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'sal')
	soda.overlay_name = 'Salinity (psu)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = 'haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'spco2')
	soda.overlay_name = 'pCO_2 (\muatm)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'revelle')
	soda.overlay_name = 'Revelle Factor';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'hco3')
	soda.overlay_name = 'HCO3 \mumol kg^-^1';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'co3')
	soda.overlay_name = 'CO3 \mumol kg^-^1';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'co2')
	soda.overlay_name = 'CO2 \mumol kg^-^1';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'fgco2')
	soda.overlay_name = 'CO2 flux mol m^-^2 yr^-^1';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'omara_satdep')
	soda.overlay_name = '\Omega_A saturation depth (m)';
	soda.overlay_pivot = 100;
	soda.overlay_min = 0;
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_min = min(soda.overlay_pivot,soda.overlay_min);
	soda.overlay_max = max(soda.overlay_pivot,soda.overlay_max);
	soda.overlay_cm = '-delta';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'omcal_satdep')
	soda.overlay_name = '\Omega_C saturation depth (m)';
	soda.overlay_pivot = 100;
	soda.overlay_min = 0;
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_min = min(soda.overlay_pivot,soda.overlay_min);
	soda.overlay_max = max(soda.overlay_pivot,soda.overlay_max);
	soda.overlay_cm = '-delta';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'omara_KT100') | strcmp(soda.overlay,'omara_K100') | strcmp(soda.overlay,'omara_KT') | strcmp(soda.overlay,'omara_KT200') | strcmp(soda.overlay,'omara_K200')  
	soda.overlay_name = '\Omega_A';
	soda.overlay_pivot = 1;
	soda.overlay_min = 0;
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_min = min(soda.overlay_pivot,soda.overlay_min);
	soda.overlay_max = max(soda.overlay_pivot,soda.overlay_max);
	soda.overlay_cm = '-delta';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'omcal_KT100') | strcmp(soda.overlay,'omcal_K100') | strcmp(soda.overlay,'omcal_KT') | strcmp(soda.overlay,'omcal_KT200') | strcmp(soda.overlay,'omcal_K200')  
	soda.overlay_name = '\Omega_C';
	soda.overlay_pivot = 1;
	soda.overlay_min = 0;
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_min = min(soda.overlay_pivot,soda.overlay_min);
	soda.overlay_max = max(soda.overlay_pivot,soda.overlay_max);
	soda.overlay_cm = '-delta';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'phtot_KT100') | strcmp(soda.overlay,'phtot_K100') | strcmp(soda.overlay,'phtot_KT100') | strcmp(soda.overlay,'phtot_KT200') | strcmp(soda.overlay,'phtot_K200') 
	soda.overlay_name = 'pH';
	soda.overlay_pivot = 7.6;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_min = min(soda.overlay_pivot,soda.overlay_min);
	soda.overlay_max = max(soda.overlay_pivot,soda.overlay_max);
	soda.overlay_cm = '-delta';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'talk_KT100')
	soda.overlay_name = 'TA (\mumol kg^-^1)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'dic_KT100')
	soda.overlay_name = 'DIC (\mumol kg^-^1)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'pco2_KT100')
	soda.overlay_name = 'pCO_2 (\muatm)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'revelle_KT100')
	soda.overlay_name = 'Revelle Factor';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'temp_KT100') 
	soda.overlay_name = 'Temperature (\circC)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = 'thermal';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'sal_KT100')
	soda.overlay_name = 'Salinity (psu)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = 'haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'po4_KT100')
	soda.overlay_name = 'PO4 (\mumol kg^-^1)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'si_KT100')
	soda.overlay_name = 'Si (\mumol kg^-^1)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'no3_KT100')
	soda.overlay_name = 'NO_3 (\mumol kg^-^1)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'o2_KT100')
	soda.overlay_name = 'O_2 (\mumol kg^-^1)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'aou_KT100')
	soda.overlay_name = 'AOU (\mumol kg^-^1)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
% elseif strcmp(soda.overlay,'epcalc100')
	% soda.overlay_name = 'epcalc100 (mol m^-^2 sec^-^1)';
	% soda.overlay_pivot = NaN;
	% soda.overlay_min = prctile(soda.overlay_data(:),10);
	% soda.overlay_max = prctile(soda.overlay_data(:),90);
	% soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	% soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
% elseif strcmp(soda.overlay,'intphyc')
	% soda.overlay_name = 'intphyc (mol m^-^2)';
	% soda.overlay_pivot = NaN;
	% soda.overlay_min = prctile(soda.overlay_data(:),10);
	% soda.overlay_max = prctile(soda.overlay_data(:),90);
	% soda.overlay_cm = 'speed';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	% soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
% elseif strcmp(soda.overlay,'intzooc')
	% soda.overlay_name = 'intzooc (mol m^-^2)';
	% soda.overlay_pivot = NaN;
	% soda.overlay_min = prctile(soda.overlay_data(:),10);
	% soda.overlay_max = prctile(soda.overlay_data(:),90);
	% soda.overlay_cm = 'speed';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	% soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
% elseif strcmp(soda.overlay,'intpp')
	% soda.overlay_name = 'intpp (mol m^-^2 sec^-^1)';
	% soda.overlay_pivot = NaN;
	% soda.overlay_min = prctile(soda.overlay_data(:),10);
	% soda.overlay_max = prctile(soda.overlay_data(:),90);
	% soda.overlay_cm = 'speed';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	% soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'ncei_omara_satdep')
	soda.overlay = ['omara_satdep'];
	soda.overlay_nc_units = ncreadatt(soda.overlay_fn2,soda.overlay,'units');  
	soda.overlay_nc_long_name = ncreadatt(soda.overlay_fn2,soda.overlay,'long_name');  
	soda.overlay_lon = ncread(soda.overlay_fn2,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn2,'lat');  
	soda.overlay_data = ncread(soda.overlay_fn2,soda.overlay);  
	soda.overlay_data = repmat(soda.overlay_data,1,1,size(soda.mapx.u,3));
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_name = '\Omega_A saturation depth (m)';
	soda.overlay_pivot = 100;
	soda.overlay_min = 0;
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_min = min(soda.overlay_pivot,soda.overlay_min);
	soda.overlay_max = max(soda.overlay_pivot,soda.overlay_max);
	soda.overlay_cm = '-delta';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'ncei_omara_KT100')
	soda.overlay = ['omara_KT100'];
	soda.overlay_nc_units = ncreadatt(soda.overlay_fn2,soda.overlay,'units');  
	soda.overlay_nc_long_name = ncreadatt(soda.overlay_fn2,soda.overlay,'long_name');  
	soda.overlay_lon = ncread(soda.overlay_fn2,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn2,'lat');  
	soda.overlay_data = ncread(soda.overlay_fn2,soda.overlay);  
	soda.overlay_data = repmat(soda.overlay_data,1,1,size(soda.mapx.u,3));
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_name = '\Omega_A';
	soda.overlay_pivot = 1;
	soda.overlay_min = 0;
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_min = min(soda.overlay_pivot,soda.overlay_min);
	soda.overlay_max = max(soda.overlay_pivot,soda.overlay_max);
	soda.overlay_cm = '-delta';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif ~strcmp(soda.overlay,'none')
	error('Program terminated because overlay were not selected correctly');
end
% use optional monthly anomaly of map overlay data
if soda.overlay_anomaly
	soda.overlay_data = monthly_anomaly(soda.overlay_data);		% convert the soda.overlay_data to the monthly anomaly
	soda.overlay_name = ['\Delta',soda.overlay_name];
	soda.overlay = [soda.overlay,'_anomaly'];
	soda.overlay_pivot = 0;
	soda.overlay_min = prctile(soda.overlay_data(:),1);
	soda.overlay_max = prctile(soda.overlay_data(:),99);
	soda.overlay_min = min(soda.overlay_pivot,soda.overlay_min);
	soda.overlay_max = max(soda.overlay_pivot,soda.overlay_max);
	soda.overlay_cm = 'balance';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
end
% % mask non-water cells as NaN
% if ~strcmp(soda.overlay,'none')
	% soda.overlay_mask = ncread(soda.overlay_fn,'mask');  
	% for i = 1:size(soda.overlay_data,3)
		% soda.datamask = squeeze(soda.overlay_data(:,:,i));
		% soda.datamask(soda.overlay_mask==0)=NaN;
		% soda.overlay_data(:,:,i) = soda.datamask;
	% end
% end

% assign gif file name for animation based on inputs above
if soda.forward_tracking
	soda.gifname = [pwd '\gif\Drifter_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_',soda.overlay,'_mov',num2str(soda.movmean_months),'mo_forward_',soda.ver,'.gif'];
else
	soda.gifname = [pwd '\gif\Drifter_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_',soda.overlay,'_mov',num2str(soda.movmean_months),'mo_reverse_',soda.ver,'.gif'];
end

% assign csv file name for output of drifter trajectory lon and lat based on inputs above
if soda.forward_tracking
	soda.csvname_lon_monthly = [pwd '\csv\Drifter_lon_monthly_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_',soda.overlay,'_mov',num2str(soda.movmean_months),'mo_forward_',soda.ver,'.csv'];
	soda.csvname_lat_monthly = [pwd '\csv\Drifter_lat_monthly_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_',soda.overlay,'_mov',num2str(soda.movmean_months),'mo_forward_',soda.ver,'.csv'];
	soda.csvname_lon_weekly = [pwd '\csv\Drifter_lon_weekly_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_',soda.overlay,'_mov',num2str(soda.movmean_months),'mo_forward_',soda.ver,'.csv'];
	soda.csvname_lat_weekly = [pwd '\csv\Drifter_lat_weekly_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_',soda.overlay,'_mov',num2str(soda.movmean_months),'mo_forward_',soda.ver,'.csv'];
else
	soda.csvname_lon_monthly = [pwd '\csv\Drifter_lon_monthly_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_',soda.overlay,'_mov',num2str(soda.movmean_months),'mo_reverse_',soda.ver,'.csv'];
	soda.csvname_lat_monthly = [pwd '\csv\Drifter_lat_monthly_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_',soda.overlay,'_mov',num2str(soda.movmean_months),'mo_reverse_',soda.ver,'.csv'];
	soda.csvname_lon_weekly = [pwd '\csv\Drifter_lon_monthly_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_',soda.overlay,'_mov',num2str(soda.movmean_months),'mo_reverse_',soda.ver,'.csv'];
	soda.csvname_lat_weekly = [pwd '\csv\Drifter_lat_weekly_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_',soda.overlay,'_mov',num2str(soda.movmean_months),'mo_reverse_',soda.ver,'.csv'];
end

% check to see if the soda simulation start year of 1980 matches the RECCAP PlankTOM12 simulation start year of 1980 for omara satdep
if ~strcmp(soda.overlay,'none') & ~(soda.startyear==1980)
	error('Program terminated because soda.startyear must be 1980 to overlay the map data correctly');
end

% ----------

% ----------
% ----------

% EXTRACT SELECTED MONTHS

% ----------
% ----------

% ----------
% v60

% % extract soda.overlay_time (different time stamp than soda.time but same months and indices, and also the overlay data runs through 2018 

soda.overlay_time = [];
% soda.overlay_time = ncread(soda.overlay_fn,'time');  
soda.overlay_datetime = datetime(1980,1:size(soda.overlay_data,3),15,0,0,0)';
soda.overlay_datenum = datenum(soda.overlay_datetime);
soda.overlay_time = soda.overlay_datenum - datenum(datetime(1980,1,1,0,0,0));	% days since 1980,1,1,0,0,0

% % extract data for the selected months 
soda.extract.time = [];
soda.extract.datenum = [];
soda.extract.datetime = [];
soda.extract.overlay_time = [];
soda.extract.overlay_datenum = [];
soda.extract.overlay_datetime = [];
soda.extract.overlay_data = [];
soda.extract.u = [];
soda.extract.v = [];
soda.extract.time(1:soda.picknmonths,1) = NaN;
soda.extract.overlay_time(1:soda.picknmonths,1) = NaN;
soda.extract.datenum(1:soda.picknmonths,1) = NaN;
soda.extract.overlay_datenum(1:soda.picknmonths,1) = NaN;
soda.extract.overlay_data(1:numel(soda.overlay_lon),1:numel(soda.overlay_lat),1:soda.picknmonths) = NaN;
soda.extract.u(1:numel(soda.mapx.lon),1:numel(soda.mapx.lat),1:soda.picknmonths) = NaN;
soda.extract.v(1:numel(soda.mapx.lon),1:numel(soda.mapx.lat),1:soda.picknmonths) = NaN;
iskip = (soda.pickstartyear - soda.startyear) * 12 + soda.pickstartmonth - 1;
for i = 1:soda.picknmonths
	soda.extract.time(i,1) = soda.time(i+iskip,1);
	soda.extract.overlay_time(i,1) = soda.overlay_time(i+iskip,1);
	soda.extract.overlay_data(:,:,i) = soda.overlay_data(:,:,i+iskip);
	soda.extract.u(:,:,i) = soda.mapx.u(:,:,i+iskip);
	soda.extract.v(:,:,i) = soda.mapx.v(:,:,i+iskip);
end
% % adjust first and last times for overlay data for interpolation range
soda.extract.overlay_time(1,1) = min(soda.extract.overlay_time(1,1),soda.extract.time(1,1));
soda.extract.overlay_time(numel(soda.extract.overlay_time),1) = ...
								max(soda.extract.overlay_time(numel(soda.extract.overlay_time),1),soda.extract.time(numel(soda.extract.time),1));
% % extract datenum and datetime
soda.extract.datenum = double(soda.extract.time + datenum(datetime(1980,1,1,0,0,0)));
soda.extract.datetime = datetime(soda.extract.datenum,'ConvertFrom','datenum');
soda.extract.overlay_datenum = soda.extract.overlay_time + datenum(datetime(1980,1,1,0,0,0));
soda.extract.overlay_datetime = datetime(soda.extract.overlay_datenum,'ConvertFrom','datenum');

% ----------

% ----------
% ----------

% INTERPOLATE TIME SERIES OF U, V, AND MAP OVERLAY

% ----------
% ----------

% ----------
% v60

% % interpolate times to time increments of soda.nsubsteps
soda.interp.time = [];
soda.interp.time(1:(soda.picknmonths-1)*soda.nsubsteps+1,1) = NaN;
soda.extract.time_diff = [];
soda.extract.time_diff = diff(soda.extract.time);
ijcount=0;
for i = 1:soda.picknmonths-1
	for j = 1:soda.nsubsteps
		ijcount = ijcount + 1;
		if ijcount == 1
			soda.interp.time(ijcount,1) = soda.extract.time(1,1);
			soda.interp.time(ijcount+1,1) = soda.interp.time(ijcount,1) + soda.extract.time_diff(i,1) / soda.nsubsteps;
			% soda.interp.overlay_time(ijcount,1) = soda.extract.overlay_time(1,1);
			% soda.interp.overlay_time(ijcount+1,1) = soda.interp.overlay_time(ijcount,1) + soda.extract.overlay_time_diff(i,1) / soda.nsubsteps;
		else
			soda.interp.time(ijcount+1,1) = soda.interp.time(ijcount,1) + soda.extract.time_diff(i,1) / soda.nsubsteps;
			% soda.interp.overlay_time(ijcount+1,1) = soda.interp.overlay_time(ijcount,1) + soda.extract.overlay_time_diff(i,1) / soda.nsubsteps;
		end
	end
end

soda.interp.time(1,1) = max(soda.interp.time(1,1),soda.extract.time(1,1));
soda.interp.time(numel(soda.interp.time),1) = ...
								min(soda.interp.time(numel(soda.interp.time),1),soda.extract.time(numel(soda.extract.time),1));

soda.interp.datenum = double(soda.interp.time + datenum(datetime(1980,1,1,0,0,0)));
soda.interp.datetime = datetime(soda.interp.datenum,'ConvertFrom','datenum');
% % interpolate overlay_data and u and v to time increments of soda.nsubsteps
soda.interp.overlay_data = [];
soda.interp.u = [];
soda.interp.v = [];
soda.interp.overlay_data(1:numel(soda.overlay_lon),1:numel(soda.overlay_lat),1:numel(soda.interp.time)) = NaN;
soda.interp.u(1:numel(soda.mapx.lon),1:numel(soda.mapx.lat),1:numel(soda.interp.time)) = NaN;
soda.interp.v(1:numel(soda.mapx.lon),1:numel(soda.mapx.lat),1:numel(soda.interp.time)) = NaN;
%
disp(['Interpolating velocity vectors ...'])
for i = 1:numel(soda.mapx.lon)
	for j = 1: numel(soda.mapx.lat)
		clear x u v xq uq vq
		x = soda.extract.datenum;
		u = squeeze(soda.extract.u(i,j,:));
		v = squeeze(soda.extract.v(i,j,:));
		xq = soda.interp.datenum;
		uq = interp1(x,u,xq);
		vq = interp1(x,v,xq);
		soda.interp.u(i,j,:) = uq;
		soda.interp.v(i,j,:) = vq;
	end
end	

%
disp(['Interpolating map overlay data ...'])
for i = 1:numel(soda.overlay_lon)
	for j = 1: numel(soda.overlay_lat)
		clear x y xq yq
		x = soda.extract.overlay_datenum;
		y = squeeze(soda.extract.overlay_data(i,j,:));
		xq = soda.interp.datenum;
		yq = interp1(x,y,xq);
		soda.interp.overlay_data(i,j,:) = yq;
	end
end	

%
% year, month, and day of interp time series
[soda.interp.year,soda.interp.month,soda.interp.day] = ymd(soda.interp.datetime);

% v60
% initial position of particles randomly placed within each initial drifter node
soda.drifter.lon0_kk = []; soda.drifter.lat0_kk = []; 
soda.drifter.lon0_kk_rsp = []; soda.drifter.lat0_kk_rsp = []; 
soda.drifter.lon0_kk(1:soda.ndrift,1:soda.drifter_particles) = NaN;
soda.drifter.lat0_kk(1:soda.ndrift,1:soda.drifter_particles) = NaN;
for i = 1:soda.ndrift
	lon_a = soda.drifter.lon0(i) - soda.drifter_width / 2;		% east side of drifter node
	lon_b = soda.drifter.lon0(i) + soda.drifter_width / 2;		% west side of drifter node
	lat_a = soda.drifter.lat0(i) - soda.drifter_height / 2;	% south side of drifter node
	lat_b = soda.drifter.lat0(i) + soda.drifter_height / 2;	% north side of drifter node
	for j = 1:soda.drifter_particles
		soda.drifter.lon0_kk(i,j) = lon_a + (lon_b - lon_a) * rand();
		soda.drifter.lat0_kk(i,j) = lat_a + (lat_b - lat_a) * rand();
	end
end
% reshape 2d arrays to 1d vectors of initial position
soda.drifter.lon0_kk_rsp = reshape(soda.drifter.lon0_kk,1,[]);
soda.drifter.lat0_kk_rsp = reshape(soda.drifter.lat0_kk,1,[]);

% error('Program stopped for debugging');

% ----------

% ----------
% ----------

% SAVE CROPPED INTERPOLATED DATA TO MAT FILE

% ----------
% ----------

% ----------

if soda.savemat

	% soda.matname = [pwd '\Drifter_crop_interp_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_mov',num2str(soda.movmean_months),'mo_',soda.ver,'.mat'];
	soda.matname = [pwd '\crop.mat'];
	soda.crop_lon_vec = ((soda.maplon(1)+360):.5:(soda.maplon(2)+360));		% west and east longitudes of cropped map extent in the animation frame (degE)
	soda.crop_lat_vec = (soda.maplat(1):.5:soda.maplat(2));					% south and north latitudes of cropped map extent (degN)
	[val,soda.crop_lon_idx] = min(abs(soda.mapx.lon-soda.crop_lon_vec));					   
	[val,soda.crop_lat_idx] = min(abs(soda.mapx.lat-soda.crop_lat_vec));	
	% create structure of cropped interpolated data to save as a mat file
	crop.map_u = soda.interp.u(soda.crop_lon_idx,soda.crop_lat_idx,:);		% 
	crop.map_v = soda.interp.v(soda.crop_lon_idx,soda.crop_lat_idx,:);		% 
	crop.map_lon = soda.crop_lon_vec;		% vector of longitudes for the full map extent (degE), must be equally spaced increments of degrees
	crop.map_lat = soda.crop_lat_vec;		% vector of latitudes for the full map extent (degN), must be equally spaced increments of degrees
	crop.datenum = soda.interp.datenum;		% vector of continuous times (datenum) for time series	
	crop.drifter_lon0 = soda.drifter.lon0_kk_rsp;		% 2-d array (size of time x numel(crop.drifter_lon0)) of initial drifter longitudes (degE)
	crop.drifter_lat0 = soda.drifter.lat0_kk_rsp;		% 2-d array (size of time x numel(crop.drifter_lon0)) of initial drifter latitudes (degN)
	save(soda.matname,'crop','-v7.3')

end

% ----------

% ----------
% ----------

% CALCULATE TRAJECTORIES

% ----------
% ----------

% ----------

if soda.integrate

	% calculate trajectories
	if soda.forward_tracking
		soda.direction = 'forward';
	else
		soda.direction = 'reverse';
	end
	soda.drifter.lon_rsp = [];
	soda.drifter.lat_rsp = [];

	% * * * * *

	[soda.drifter.lon_rsp, soda.drifter.lat_rsp] = drifter_v04(soda.mapx.lon, soda.mapx.lat, soda.interp.time, ...
					soda.interp.u, soda.interp.v, ...
					soda.drifter.lon0_kk_rsp, soda.drifter.lat0_kk_rsp, soda.direction);

	% * * * * *

	% reshape drifter.lon and drifter.lat to 3d (time x ndrift x drifter_particles)
	soda.drifter.lon = [];
	soda.drifter.lat = [];
	soda.drifter.lon = reshape(soda.drifter.lon_rsp,[],soda.ndrift,soda.drifter_particles);
	soda.drifter.lat = reshape(soda.drifter.lat_rsp,[],soda.ndrift,soda.drifter_particles);

end

% ----------

% ----------
% ----------

% INTERPOLATE OVERLAY DATA TO MATCH DRIFTER TRAJECTORIES

% ----------
% ----------

% ----------

if ~strcmp(soda.overlay,'none')

	disp(['Interpolating overlay data along drifter trajectories ...'])

	% input arguments for interp2_drifter_overlay.m
	soda.drifter.overlay = [];
	soda.drifter.overlay_rsp = [];
	drifter_lon = soda.drifter.lon_rsp;			% degE
	drifter_lat = soda.drifter.lat_rsp;			% degN
	overlay_lon = soda.overlay_lon+360; 		% convert from -degW to degE
	overlay_lat = soda.overlay_lat;				% degN
	overlay_data = soda.interp.overlay_data;

	% * * * * *

	[soda.drifter.overlay_rsp] = interp2_drifter_overlay(drifter_lon, drifter_lat, overlay_lon, overlay_lat, overlay_data);

	% * * * * *

	% reshape to 3d array (time x soda.ndrift x soda.drifter_particles)
	soda.drifter.overlay = reshape(soda.drifter.overlay_rsp,[],soda.ndrift,soda.drifter_particles);
	% summary stats by soda.ndrift across soda.drifter_particles
	soda.drifter.overlay_mean = squeeze(nanmean2(soda.drifter.overlay,3));
	soda.drifter.overlay_std = squeeze(nanstd2(soda.drifter.overlay,3));
	soda.drifter.overlay_min = squeeze(nanmin2(soda.drifter.overlay,3));
	soda.drifter.overlay_max = squeeze(nanmax2(soda.drifter.overlay,3));
	% summary stats by soda.ndrift across soda.drifter_particles
	soda.drifter.lon_mean = squeeze(nanmean2(soda.drifter.lon,3));
	soda.drifter.lat_mean = squeeze(nanmean2(soda.drifter.lat,3));
	
	% save drifter lon, lat, and interpolated overlay data as nc file
	if soda.save_drifter_track
	disp(['Saving netcdf and csv files of overlay data along drifter trajectories ...'])

		if soda.forward_tracking
			soda.drifter_nc_name = [pwd '\nc\Drifter_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_forward_',soda.ver,'.nc'];
			soda.drifter_csv_name_mean = [pwd '\csv\Drifter_mean_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_forward_',soda.ver,'.csv'];
			soda.drifter_csv_name_stdev = [pwd '\csv\Drifter_stdev_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_forward_',soda.ver,'.csv'];
			soda.drifter_csv_name_min = [pwd '\csv\Drifter_min_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_forward_',soda.ver,'.csv'];
			soda.drifter_csv_name_max = [pwd '\csv\Drifter_max_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_forward_',soda.ver,'.csv'];
			soda.drifter_csv_name_mean_lon = [pwd '\csv\Drifter_mean_lon_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_forward_',soda.ver,'.csv'];
			soda.drifter_csv_name_mean_lat = [pwd '\csv\Drifter_mean_lat_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_forward_',soda.ver,'.csv'];
		else
			soda.drifter_nc_name = [pwd '\nc\Drifter_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_reverse_',soda.ver,'.nc'];
			soda.drifter_csv_name_mean = [pwd '\csv\Drifter_mean_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_reverse_',soda.ver,'.csv'];
			soda.drifter_csv_name_stdev = [pwd '\csv\Drifter_stdev_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_reverse_',soda.ver,'.csv'];
			soda.drifter_csv_name_min = [pwd '\csv\Drifter_min_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_reverse_',soda.ver,'.csv'];
			soda.drifter_csv_name_max = [pwd '\csv\Drifter_max_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_reverse_',soda.ver,'.csv'];
			soda.drifter_csv_name_mean_lon = [pwd '\csv\Drifter_mean_lon_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_reverse_',soda.ver,'.csv'];
			soda.drifter_csv_name_mean_lat = [pwd '\csv\Drifter_mean_lat_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_reverse_',soda.ver,'.csv'];
		end
		%
		% write nc file
		fout = soda.drifter_nc_name;
		if exist(fout)
			delete(fout);
		end
		NX = numel(overlay_lon);
		NY = numel(overlay_lat);
		NT = numel(soda.interp.time);
		ND = soda.ndrift;
		NP = soda.drifter_particles;
		% create
		nccreate(fout,'time','Dimensions',{'time',NT},'Format','netcdf4','FillValue',NaN);   % create nc file with variables with needed dimensions
		nccreate(fout,'drifter_lon0','Dimensions',{'drifter',ND,'particle',NP},'Format','netcdf4','FillValue',NaN);   % create nc file with variables with needed dimensions
		nccreate(fout,'drifter_lat0','Dimensions',{'drifter',ND,'particle',NP},'Format','netcdf4','FillValue',NaN);   % create nc file with variables with needed dimensions
		nccreate(fout,'drifter_lon','Dimensions',{'time',NT,'drifter',ND,'particle',NP},'Format','netcdf4','FillValue',NaN);   % create nc file with variables with needed dimensions
		nccreate(fout,'drifter_lat','Dimensions',{'time',NT,'drifter',ND,'particle',NP},'Format','netcdf4','FillValue',NaN);   % create nc file with variables with needed dimensions
		nccreate(fout,'drifter_overlay','Dimensions',{'time',NT,'drifter',ND,'particle',NP},'Format','netcdf4','FillValue',NaN);   % create nc file with variables with needed dimensions
		nccreate(fout,'drifter_overlay_mean','Dimensions',{'time',NT,'drifter',ND},'Format','netcdf4','FillValue',NaN);   % create nc file with variables with needed dimensions
		nccreate(fout,'drifter_overlay_std','Dimensions',{'time',NT,'drifter',ND},'Format','netcdf4','FillValue',NaN);   % create nc file with variables with needed dimensions
		nccreate(fout,'drifter_overlay_min','Dimensions',{'time',NT,'drifter',ND},'Format','netcdf4','FillValue',NaN);   % create nc file with variables with needed dimensions
		nccreate(fout,'drifter_overlay_max','Dimensions',{'time',NT,'drifter',ND},'Format','netcdf4','FillValue',NaN);   % create nc file with variables with needed dimensions
		% write
		ncwrite(fout,'time',soda.interp.time);   % write the data into the variable in the nc file
		ncwrite(fout,'drifter_lon0',soda.drifter.lon0_kk);   % write the data into the variable in the nc file
		ncwrite(fout,'drifter_lat0',soda.drifter.lat0_kk);   % write the data into the variable in the nc file
		ncwrite(fout,'drifter_lon',soda.drifter.lon);   % write the data into the variable in the nc file
		ncwrite(fout,'drifter_lat',soda.drifter.lon);   % write the data into the variable in the nc file
		ncwrite(fout,'drifter_overlay',soda.drifter.overlay);   % write the data into the variable in the nc file
		ncwrite(fout,'drifter_overlay_mean',soda.drifter.overlay_mean);   % write the data into the variable in the nc file
		ncwrite(fout,'drifter_overlay_std',soda.drifter.overlay_std);   % write the data into the variable in the nc file
		ncwrite(fout,'drifter_overlay_min',soda.drifter.overlay_min);   % write the data into the variable in the nc file
		ncwrite(fout,'drifter_overlay_max',soda.drifter.overlay_max);   % write the data into the variable in the nc file
		% writeatt
		ncwriteatt(fout,'time','units','days')
		ncwriteatt(fout,'time','long_name','days since 1980-01-01 00:00:00')
		if soda.forward_tracking
			ncwriteatt(fout,'drifter_lon0','units','degrees_east')
			ncwriteatt(fout,'drifter_lon0','long_name','initial longitude (positve east)')
			ncwriteatt(fout,'drifter_lat0','units','degrees_north')
			ncwriteatt(fout,'drifter_lat0','long_name','initial latitude (positve north)')
		else
			ncwriteatt(fout,'drifter_lon0','units','degrees_east')
			ncwriteatt(fout,'drifter_lon0','long_name','final longitude (positve east)')
			ncwriteatt(fout,'drifter_lat0','units','degrees_north')
			ncwriteatt(fout,'drifter_lat0','long_name','final latitude (positve north)')
		end
		ncwriteatt(fout,'drifter_lon','units','degrees_east')
		ncwriteatt(fout,'drifter_lon','long_name','trajectory longitude (positve east)')
		ncwriteatt(fout,'drifter_lat','units','degrees_north')
		ncwriteatt(fout,'drifter_lat','long_name','trajectory latitude (positve north)')
		ncwriteatt(fout,'drifter_overlay','units',soda.overlay_nc_units)
		ncwriteatt(fout,'drifter_overlay','long_name',soda.overlay_nc_long_name)
		ncwriteatt(fout,'drifter_overlay_mean','units',soda.overlay_nc_units)
		ncwriteatt(fout,'drifter_overlay_mean','long_name',[soda.overlay_nc_long_name, ', mean of ',num2str(NP),' particle tracks'])
		ncwriteatt(fout,'drifter_overlay_std','units',soda.overlay_nc_units)
		ncwriteatt(fout,'drifter_overlay_std','long_name',[soda.overlay_nc_long_name, ', std dev of ',num2str(NP),' particle tracks'])
		ncwriteatt(fout,'drifter_overlay_min','units',soda.overlay_nc_units)
		ncwriteatt(fout,'drifter_overlay_min','long_name',[soda.overlay_nc_long_name, ', min of ',num2str(NP),' particle tracks'])
		ncwriteatt(fout,'drifter_overlay_max','units',soda.overlay_nc_units)
		ncwriteatt(fout,'drifter_overlay_max','long_name',[soda.overlay_nc_long_name, ', max of ',num2str(NP),' particle tracks'])
		% write global attributes
		ncwriteatt(fout,'/','title','Drifter trajectories and overlay biogeochemical exposure data')
		ncwriteatt(fout,'/','long_title','Drifter trajectories calculated from SODA3.12.2 u and v, and biogeochemical exposure data from PlankTOM12')
		ncwriteatt(fout,'/','institution','Greg Pelletier')
		ncwriteatt(fout,'/','source','Calculated using soda_drifter.m by Greg Pelletier, available at https://github.com/gjpelletier/soda_drifter')
		%
		% write csv file
		% delete old csv with this name if it exists
		if exist(soda.drifter_csv_name_mean, 'file')==2
		  delete(soda.drifter_csv_name_mean);
		end
		if exist(soda.drifter_csv_name_stdev, 'file')==2
		  delete(soda.drifter_csv_name_stdev);
		end
		if exist(soda.drifter_csv_name_min, 'file')==2
		  delete(soda.drifter_csv_name_min);
		end
		if exist(soda.drifter_csv_name_max, 'file')==2
		  delete(soda.drifter_csv_name_max);
		end
		% column names
		clear colname
		colname = strings();
		colname(1)='excel_datenum';
		kkcount = 1;
		for j = 1:soda.ndrift
			kkcount=kkcount+1;
			if numel(soda.drifter.sta)>0
				colname(kkcount) = strcat('sta_',soda.drifter.sta(j));
			else
				colname(kkcount) = ['sta_',num2str(j,'%03.f')];
			end
		end
		col1 = soda.interp.datenum - 693960;	% convert matlab datenum to excel datenum, difference between 1/1/1900 vs 1/1/0000 = 693960 days
		% mean overlay
		clear Mtx_mean T_mean
		Mtx_mean = [col1, soda.drifter.overlay_mean];  
		T_mean=array2table(Mtx_mean);
		T_mean.Properties.VariableNames(1:numel(colname))=colname;
		writetable(T_mean,soda.drifter_csv_name_mean);
		% stdev overlay
		clear Mtx_stdev T_stdev
		Mtx_stdev = [col1, soda.drifter.overlay_std];  
		T_stdev=array2table(Mtx_stdev);
		T_stdev.Properties.VariableNames(1:numel(colname))=colname;
		writetable(T_stdev,soda.drifter_csv_name_stdev);
		% min overlay
		clear Mtx_min T_min
		Mtx_min = [col1, soda.drifter.overlay_min];  
		T_min=array2table(Mtx_min);
		T_min.Properties.VariableNames(1:numel(colname))=colname;
		writetable(T_min,soda.drifter_csv_name_min);
		% max overlay
		clear Mtx_max T_max
		Mtx_max = [col1, soda.drifter.overlay_max];  
		T_max=array2table(Mtx_max);
		T_max.Properties.VariableNames(1:numel(colname))=colname;
		writetable(T_max,soda.drifter_csv_name_max);
		% mean lon
		clear Mtx_mean T_mean
		Mtx_mean = [col1, soda.drifter.lon_mean];  
		T_mean=array2table(Mtx_mean);
		T_mean.Properties.VariableNames(1:numel(colname))=colname;
		writetable(T_mean,soda.drifter_csv_name_mean_lon);
		% mean overlay
		clear Mtx_mean T_mean
		Mtx_mean = [col1, soda.drifter.lat_mean];  
		T_mean=array2table(Mtx_mean);
		T_mean.Properties.VariableNames(1:numel(colname))=colname;
		writetable(T_mean,soda.drifter_csv_name_mean_lat);

		% save csv of summary stats of overlay data on drifter trajectory track lines by station
		% column names
		clear colname col1 col2 col3 col4 col5 col6
		colname = strings();
		colname(1)='station';
		colname(2)='mean';
		colname(3)='stdev';
		colname(4)='min';
		colname(5)='max';
		colname(6)='median';
		col1 = soda.drifter.sta';
		col2 = nanmean2(soda.drifter.overlay_mean)';	% mean	
		col3 = nanstd2(soda.drifter.overlay_mean)';		% stdev
		col4 = nanmin2(soda.drifter.overlay_mean)';		% min
		col5 = nanmax2(soda.drifter.overlay_mean)';		% max
		col6 = nanmedian2(soda.drifter.overlay_mean)';	% median
		clear Mtx Tx Tname
		if soda.forward_tracking
			Tname = [pwd '\csv\Drifter_summary_stats_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_forward_',soda.ver,'.csv'];
		else	
			Tname = [pwd '\csv\Drifter_summary_stats_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_reverse_',soda.ver,'.csv'];
		end
		Mtx = [col1, col2, col3, col4 col5 col6];  
		Tx=array2table(Mtx);
		Tx.Properties.VariableNames(1:numel(colname))=colname;
		writetable(Tx,Tname);

	end

	if soda.plot_drifter_track
	disp(['Plotting exposure tracks of overlay data along drifter trajectories ...'])

		if strcmp(soda.overlay,'omara_K200') | strcmp(soda.overlay,'omcal_K200') | strcmp(soda.overlay,'phtot_K200')
			strdep = ' at 200m';
		elseif strcmp(soda.overlay,'omara_KT200') | strcmp(soda.overlay,'omcal_KT200') | strcmp(soda.overlay,'phtot_KT200')
			strdep = ' at 0-200m';
		elseif strcmp(soda.overlay,'omara_K100') | strcmp(soda.overlay,'omcal_K100') | strcmp(soda.overlay,'phtot_K100')
			strdep = ' at 100m';
		elseif strcmp(soda.overlay,'omara_KT100') | strcmp(soda.overlay,'omcal_KT100') | strcmp(soda.overlay,'phtot_KT100')
			strdep = ' at 0-100m';
		else
			strdep = '';
		end

		% - - -
		% a) map of tracklines and b) time series of mean exposure on trackline
		ck=colormap(jet2(soda.ndrift));
		% if strcmp(soda.whichmap,'Humboldt')
		if soda.maplat(2) <= 0
			ck = flipud(ck);
		end
		figure(1)
		hFig=gcf;
		clf(hFig);
		% [ (inner vert) (inner horiz) ]. [ (Bottom margin) (Top margin) ]. [ (L margin) (R margin) ] 
		% ha = tight_subplot(2,1,[.06 .035],[.06 .06],[.08 .04]);
		ha = tight_subplot(2,1,[.08 .035],[.08 .06],[.10 .04]);
		axes(ha(2))
		% subplot(2,1,1)
		hold on
		for i=1:soda.ndrift
			plot(soda.interp.datetime,movmean(soda.drifter.overlay_mean(:,i),soda.nsubsteps),'color',ck(i,:),'linewidth',1.5)
		end
		set(gca,'TickDir','out');	% the only other option is 'in'
		ylabel([soda.overlay_name,strdep]);
		if soda.forward_tracking
			title(['b._ ', soda.overlay_name,strdep,' on forward tracks',strdep1,' from ',num2str(soda.pickstartyear),'-',num2str(soda.pickstartmonth,'%02.f'),' to ',num2str(soda.pickendyear),'-',num2str(soda.pickendmonth,'%02.f')]);
		else
			title(['b._ ', soda.overlay_name,strdep,' on reverse tracks',strdep1,' from ',num2str(soda.pickstartyear),'-',num2str(soda.pickstartmonth,'%02.f'),' to ',num2str(soda.pickendyear),'-',num2str(soda.pickendmonth,'%02.f')]);
		end
		if numel(soda.drifter.sta)>0
			legend(soda.drifter.sta,'location','bestoutside')
		else
			% legend(colname(2:end),'location','bestoutside')
			clear legname
			legname = strings;
			for j = 1:soda.ndrift
				legname(j) = ['sta',num2str(j,'%02.f')];
			end
			legend(legname,'location','bestoutside')
		end
		ax = gca;
		ax.XAxis.TickValues = datetime(1980,1,1)+calquarters(0:(2018-1980)*4);
		% ax.XAxis.TickLabelFormat = 'QQQ-yyyy';
		ax.XAxis.TickLabelFormat = 'MMM-yyyy';
		xlim([min(soda.interp.datetime) max(soda.interp.datetime)])
		% yticks([0 .5 1 1.5 2 2.5 3 3.5 4 4.5 5]);
		% yticklabels([0:.5:5]);
		yticks('auto')
		yt=yticks;	% returns the current y-axis tick values as vector yt
		yticklabels(yt);
		grid on
		hold off
		%
		axes(ha(1))
		% subplot(2,1,2)
		hold on
		m_proj('mercator','lat',soda.maplat,'long',soda.maplon);
		if ~strcmp(soda.overlay,'none')
			clear X Y Z
			X = soda.overlay_lon;
			Y = soda.overlay_lat;
			Z = squeeze(nanmean2(soda.interp.overlay_data,3))';
			[C,h] = m_contourf(X,Y,Z,256);
			set(h,'LineColor','none')
			caxis([soda.overlay_min soda.overlay_max]);
			h=colorbar;
			h.Label.String=['Average ',soda.overlay_name,strdep];
			if strcmp(soda.overlay_cm,'balance') | strcmp(soda.overlay_cm,'-balance') | ...
					strcmp(soda.overlay_cm,'delta') | strcmp(soda.overlay_cm,'-delta') | ...
					strcmp(soda.overlay_cm,'curl') | strcmp(soda.overlay_cm,'-curl')
				if ~isnan(soda.overlay_pivot)
					cmocean2(soda.overlay_cm,'pivot',soda.overlay_pivot,'opacity',soda.overlay_opacity);
				else
					cmocean2(soda.overlay_cm,'opacity',soda.overlay_opacity);
				end
			else 
				cmocean2(soda.overlay_cm,'opacity',soda.overlay_opacity);
			end
		end
		m_gshhs_i('patch',[.9 .9 .9],'edgecolor',[.9 .9 .9]);		% intermediate res coast
		for i=1:soda.ndrift
			% initial or final position
			m_plot(soda.drifter.lon0(i)-360, soda.drifter.lat0(i),'o','linewidth',1,'markeredgecolor',ck(i,:),'markerfacecolor','k','MarkerSize',soda.markersize)
			% track line
			m_plot(soda.drifter.lon_mean(:,i)-360,soda.drifter.lat_mean(:,i),'color',ck(i,:),'linewidth',1.5)
		end
		if soda.forward_tracking
			title(['a._ Forward tracks of drifters',strdep1,' from_ ',num2str(soda.pickstartyear),'-',num2str(soda.pickstartmonth,'%02.f'),' to ',num2str(soda.pickendyear),'-',num2str(soda.pickendmonth,'%02.f')]);
		else
			title(['a._ Reverse tracks of drifters',strdep1,' from_ ',num2str(soda.pickstartyear),'-',num2str(soda.pickstartmonth,'%02.f'),' to ',num2str(soda.pickendyear),'-',num2str(soda.pickendmonth,'%02.f')]);
		end
		m_grid('box','off','linestyle','none','XaxisLocation','bottom','YaxisLocation','left','xtick',([-240 -230 -220 -210 -200 -190 -180 -170 -160 -150 -140 -130 -120 -110 -100 -90 -80 -70 -60]),'ytick',([-65 -60 -50 -40 -30 -20 -10 0 10 20 30 40 50 60 65]),'FontSize',8);	% N Pacific
		hold off
		%
		set(gcf, 'PaperPosition', [0 0 6 9])   
		if soda.forward_tracking
			print(gcf, [pwd '/png/Drifter_tracks_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_forward_',soda.ver,'.png'], '-dpng', '-r600' );   %save file as PNG w/ 300dpi
		else
			print(gcf, [pwd '/png/Drifter_tracks_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_reverse_',soda.ver,'.png'], '-dpng', '-r600' );   %save file as PNG w/ 300dpi
		end

		% - - -
		% box plot of overlay variable by station
		% ck=colormap(jet2(soda.ndrift));
		figure(1)
		hFig=gcf;
		clf(hFig);
		% T=bplot(soda.drifter.overlay_mean);
		hold on
		boxplot(soda.drifter.overlay_mean,'labels',soda.drifter.sta);
		h = findobj(gcf,'tag','Outliers');
		set(h,'MarkerSize',.5)		
		% overlay mean values as black diamonds
		plot(nanmean2(soda.drifter.overlay_mean), 'dk')
		ylim([nanmin2(soda.drifter.overlay_mean(:)) nanmax2(soda.drifter.overlay_mean(:))]);
		% boxplot(soda.drifter.overlay_mean,'labels',soda.drifter.sta,'symbol','w+');
		if strcmp(soda.overlay,'omara_KT100') | strcmp(soda.overlay,'omcal_KT100') 
			hline(1,'k:','linewidth',1,'color',[.5 .5 .5]);
		elseif strcmp(soda.overlay,'omara_satdep') | strcmp(soda.overlay,'omcal_satdep') 
			hline(100,'k:','linewidth',1,'color',[.5 .5 .5]);
		end
		if soda.forward_tracking
			title(['c._ Boxplot of ',soda.overlay_name,strdep,' on forward tracks',strdep1,' from ',num2str(soda.pickstartyear),'-',num2str(soda.pickstartmonth,'%02.f'),' to ',num2str(soda.pickendyear),'-',num2str(soda.pickendmonth,'%02.f')]);
		else
			title(['c._ Boxplot of ',soda.overlay_name,strdep,' on reverse tracks',strdep1,' from ',num2str(soda.pickstartyear),'-',num2str(soda.pickstartmonth,'%02.f'),' to ',num2str(soda.pickendyear),'-',num2str(soda.pickendmonth,'%02.f')]);
		end
		ylabel([soda.overlay_name,strdep]);
		% grid on
		set(gca, 'YGrid', 'on', 'XGrid', 'off');
		hold off
		set(gcf, 'PaperPosition', [0 0 7.5 5])   
		if soda.forward_tracking
			print(gcf, [pwd '/png/Boxplot_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_forward_',soda.ver,'.png'], '-dpng', '-r600' );   %save file as PNG w/ 300dpi
		else
			print(gcf, [pwd '/png/Boxplot_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_reverse_',soda.ver,'.png'], '-dpng', '-r600' );   %save file as PNG w/ 300dpi
		end

		% - - -
		% Inset map of drifter tracks zoomed into soda.croplat soda.croplon
		if soda.extraplots & numel(soda.croplat)==2 & numel(soda.croplon)==2
			% - - -
			% inset crop map of mean drifter tracks
			figure(1)
			hFig=gcf;
			clf(hFig);
			%
			hold on
			m_proj('mercator','lat',soda.croplat,'long',soda.croplon);
			if ~strcmp(soda.overlay,'none')
				clear X Y Z
				X = soda.overlay_lon;
				Y = soda.overlay_lat;
				Z = squeeze(nanmean2(soda.interp.overlay_data,3))';
				[C,h] = m_contourf(X,Y,Z,256);
				set(h,'LineColor','none')
				caxis([soda.overlay_min soda.overlay_max]);
				h=colorbar;
				h.Label.String=['Average ',soda.overlay_name,strdep];
				% v49
				if strcmp(soda.overlay_cm,'balance') | strcmp(soda.overlay_cm,'-balance') | ...
						strcmp(soda.overlay_cm,'delta') | strcmp(soda.overlay_cm,'-delta') | ...
						strcmp(soda.overlay_cm,'curl') | strcmp(soda.overlay_cm,'-curl')
					if ~isnan(soda.overlay_pivot)
						cmocean2(soda.overlay_cm,'pivot',soda.overlay_pivot,'opacity',soda.overlay_opacity);
					else
						cmocean2(soda.overlay_cm,'opacity',soda.overlay_opacity);
					end
				else 
					cmocean2(soda.overlay_cm,'opacity',soda.overlay_opacity);
				end
			end
			m_gshhs_i('patch',[.9 .9 .9],'edgecolor',[.9 .9 .9]);		% intermediate res coast
			for i=1:soda.ndrift
				% initial or final position
				m_plot(soda.drifter.lon0(i)-360, soda.drifter.lat0(i),'o','linewidth',1,'markeredgecolor',ck(i,:),'markerfacecolor','k','MarkerSize',soda.markersize)
				% track line
				m_plot(soda.drifter.lon_mean(:,i)-360,soda.drifter.lat_mean(:,i),'color',ck(i,:),'linewidth',1.5)
			end
			if soda.forward_tracking
				title(['d._ Forward tracks',strdep1,' from_ ',num2str(soda.pickstartyear),'-',num2str(soda.pickstartmonth,'%02.f'),' to ',num2str(soda.pickendyear),'-',num2str(soda.pickendmonth,'%02.f')]);
			else
				title(['d._ Reverse tracks',strdep1,' from_ ',num2str(soda.pickstartyear),'-',num2str(soda.pickstartmonth,'%02.f'),' to ',num2str(soda.pickendyear),'-',num2str(soda.pickendmonth,'%02.f')]);
			end
			m_grid('box','off','linestyle','none','XaxisLocation','bottom','YaxisLocation','left','xtick',([-240 -230 -220 -210 -200 -190 -180 -170 -160 -150 -140 -130 -120 -110 -100 -90 -80 -70 -60]),'ytick',([-65 -60 -50 -40 -30 -20 -10 0 10 20 30 40 50 60 65]),'FontSize',8);	% N Pacific
			hold off
			%
			if strcmp(soda.croporient,'portrait')
				set(gcf, 'PaperPosition', [0 0 6 8])   
			else
				set(gcf, 'PaperPosition', [0 0 8 6])   
			end
			if soda.forward_tracking
				print(gcf, [pwd '/png/Drifter_tracks_crop_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_forward_',soda.ver,'.png'], '-dpng', '-r600' );   %save file as PNG w/ 300dpi
			else
				print(gcf, [pwd '/png/Drifter_tracks_crop_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_reverse_',soda.ver,'.png'], '-dpng', '-r600' );   %save file as PNG w/ 300dpi
			end
		end 	% if print zoomed inset map
		
	end		% if soda.plot_drifter_track

end 	% if ~strcmp(soda.overlay,'none')

% ----------

% ----------
% ----------

% FILTER WEEKLY FOR ANIMATION FRAMES

% ----------
% ----------

% ----------

% if soda.integrate
if ~strcmp(soda.overlay,'none') & soda.makegif

	soda.drifter.weekly_lon_rsp = [];
	soda.drifter.weekly_lat_rsp = [];
	soda.drifter.weekly_lon = [];
	soda.drifter.weekly_lat = [];
	soda.drifter.weekly_time = [];
	soda.drifter.weekly_overlay_data = [];
	soda.drifter.weekly_u = [];
	soda.drifter.weekly_v = [];
	soda.timestep = nanmean2(diff(soda.interp.time));	% average time step in days 
	iskip = floor(soda.skipdays / soda.timestep);
	nweeks = floor(numel(soda.interp.time)/iskip);
	soda.drifter.weekly_lon(1:nweeks,1:soda.ndrift,1:soda.drifter_particles) = NaN;
	soda.drifter.weekly_lat(1:nweeks,1:soda.ndrift,1:soda.drifter_particles) = NaN;
	soda.drifter.weekly_time(1:nweeks,1) = NaN;
	soda.drifter.weekly_overlay_data(1:numel(soda.overlay_lon),1:numel(soda.overlay_lat),1:nweeks) = NaN;
	soda.drifter.weekly_u(1:numel(soda.mapx.lon),1:numel(soda.mapx.lat),1:nweeks) = NaN;
	soda.drifter.weekly_v(1:numel(soda.mapx.lon),1:numel(soda.mapx.lat),1:nweeks) = NaN;
	iweek = 0;
	icount = iskip - 1;
	for i = 1:numel(soda.interp.time)
		icount = icount + 1;

		%v61 include the last row of original data in the weekly filter
		% if icount == iskip
		if (icount == iskip) | (i == numel(soda.interp.time) & icount > 0)

			iweek = iweek + 1;
			icount = 0;

			if soda.integrate
				soda.drifter.weekly_lon_rsp(iweek,:) = soda.drifter.lon_rsp(i,:);
				soda.drifter.weekly_lat_rsp(iweek,:) = soda.drifter.lat_rsp(i,:);
				soda.drifter.weekly_lon(iweek,:,:) = soda.drifter.lon(i,:,:);
				soda.drifter.weekly_lat(iweek,:,:) = soda.drifter.lat(i,:,:);
			end

			soda.drifter.weekly_time(iweek,1) = soda.interp.time(i,1);
			soda.drifter.weekly_overlay_data(:,:,iweek) = soda.interp.overlay_data(:,:,i);
			soda.drifter.weekly_u(:,:,iweek) = soda.interp.u(:,:,i);
			soda.drifter.weekly_v(:,:,iweek) = soda.interp.v(:,:,i);
		end
	end

	% year, month, and day of weekly time series
	soda.drifter.weekly_datenum = double(soda.drifter.weekly_time + datenum(datetime(1980,1,1,0,0,0)));
	soda.drifter.weekly_datetime = datetime(soda.drifter.weekly_datenum,'ConvertFrom','datenum');
	[soda.drifter.weekly_year,soda.drifter.weekly_month,soda.drifter.weekly_day] = ymd(soda.drifter.weekly_datetime);

end		% if ~strcmp(soda.overlay,'none') & soda.makegif

% error('Program stopped for debugging');

% ----------

% ----------
% ----------

% MAKE ANIMATION

% ----------
% ----------

% ----------
% v60

% * * * * *

if soda.makegif

	if soda.integrate
		soda_drifter_animation_v05(soda.mapx.lon-360, soda.mapx.lat, soda.drifter.weekly_u, soda.drifter.weekly_v, ...
									soda.quiverscale, soda.quivercolor, soda.quiverwidth, soda.plot_title, soda.overlay_quiver, ...
									soda.drifter.weekly_lon-360, soda.drifter.weekly_lat, soda.overlay_drifter, ...
									soda.markerfacecolor, soda.markeredgecolor, soda.markersize, soda.integrate, ...
									soda.overlay, soda.overlay_name, ...
									soda.overlay_lon, soda.overlay_lat, soda.drifter.weekly_overlay_data, ...
									soda.overlay_min, soda.overlay_max, soda.overlay_cm, soda.overlay_pivot, soda.overlay_opacity, ...
									soda.maplon, soda.maplat, soda.whichmap, ...
									soda.drifter.weekly_time, soda.gifname, soda.gifdelay, soda.gifpause)
	else
		soda_drifter_animation_v05(soda.mapx.lon-360, soda.mapx.lat, [], [], ...
									soda.quiverscale, soda.quivercolor, soda.quiverwidth, soda.plot_title, soda.overlay_quiver, ...
									[], [], soda.overlay_drifter, ...
									soda.markerfacecolor, soda.markeredgecolor, soda.markersize, soda.integrate, ...
									soda.overlay, soda.overlay_name, ...
									soda.overlay_lon, soda.overlay_lat, soda.drifter.weekly_overlay_data, ...
									soda.overlay_min, soda.overlay_max, soda.overlay_cm, soda.overlay_pivot, soda.overlay_opacity, ...
									soda.maplon, soda.maplat, soda.whichmap, ...
									soda.drifter.weekly_time, soda.gifname, soda.gifdelay, soda.gifpause)
	end

end


toc



% ----------

% ----------
% ----------

% EXTRA PLOTS WITH THRESHOLDS IF OVERLAY IS omara

% ----------
% ----------

% ----------
% Overlay omega_ar threholds on drifter track exposure plot
if soda.extraplots

	if strcmp(soda.overlay,'omara_K200') | ...
		strcmp(soda.overlay,'omara_K100') | ...
		strcmp(soda.overlay,'omara_KT200') | ...
		strcmp(soda.overlay,'omara_KT100') | ...
		strcmp(soda.overlay,'omara_KT') | ...
		strcmp(soda.overlay,'phtot_K200') | ...
		strcmp(soda.overlay,'phtot_K100') | ...
		strcmp(soda.overlay,'phtot_KT200') | ...
		strcmp(soda.overlay,'phtot_KT100') | ...
		strcmp(soda.overlay,'phtot_KT')

		if strcmp(soda.overlay,'omara_K200') | strcmp(soda.overlay,'omcal_K200') | strcmp(soda.overlay,'phtot_K200')	
			strdep = ' at 200m';
		elseif strcmp(soda.overlay,'omara_KT200') | strcmp(soda.overlay,'omcal_KT200') | strcmp(soda.overlay,'phtot_KT200')
			strdep = ' at 0-200m';
		elseif strcmp(soda.overlay,'omara_K100') | strcmp(soda.overlay,'omcal_K100') | strcmp(soda.overlay,'phtot_K100')
			strdep = ' at 100m';
		elseif strcmp(soda.overlay,'omara_KT100') | strcmp(soda.overlay,'omcal_KT100') | strcmp(soda.overlay,'phtot_KT100')
			strdep = ' at 0-100m';
		else
			strdep = '';
		end

		% - - -
		% a) map of tracklines and b) time series of mean exposure on trackline
		ck=colormap(jet2(soda.ndrift));
		if strcmp(soda.whichmap,'Humboldt')
			ck = flipud(ck);
		end
		figure(1)
		hFig=gcf;
		clf(hFig);
		% [ (inner vert) (inner horiz) ]. [ (Bottom margin) (Top margin) ]. [ (L margin) (R margin) ] 
		% ha = tight_subplot(2,1,[.06 .035],[.06 .06],[.08 .04]);
		ha = tight_subplot(2,1,[.08 .035],[.08 .06],[.10 .04]);
		axes(ha(2))
		% subplot(2,1,1)
		hold on
		%
		if strcmp(soda.overlay,'omara_K200') | ...
			strcmp(soda.overlay,'omara_K100') | ...
			strcmp(soda.overlay,'omara_KT200') | ...
			strcmp(soda.overlay,'omara_KT100') | ...
			strcmp(soda.overlay,'omara_KT')
				plot(soda.interp.datetime,ones(size(soda.interp.datetime))*1.2,'k-','linewidth',1)
				plot(soda.interp.datetime,ones(size(soda.interp.datetime))*.95,'k-','linewidth',1.5)
		elseif	strcmp(soda.overlay,'phtot_K200') | ...
			strcmp(soda.overlay,'phtot_K100') | ...
			strcmp(soda.overlay,'phtot_KT200') | ...
			strcmp(soda.overlay,'phtot_KT100') | ...
			strcmp(soda.overlay,'phtot_KT')
				plot(soda.interp.datetime,ones(size(soda.interp.datetime))*7.75,'k-','linewidth',1)
				plot(soda.interp.datetime,ones(size(soda.interp.datetime))*7.4,'k-','linewidth',1.5)
		end
		%
		for i=1:soda.ndrift
			plot(soda.interp.datetime,movmean(soda.drifter.overlay_mean(:,i),soda.nsubsteps),'color',ck(i,:),'linewidth',1.5)
		end
		set(gca,'TickDir','out');	% the only other option is 'in'
		ylabel([soda.overlay_name,strdep]);
		if soda.forward_tracking
			title(['b._ ', soda.overlay_name,strdep,' on forward tracks',strdep1,' from ',num2str(soda.pickstartyear),'-',num2str(soda.pickstartmonth,'%02.f'),' to ',num2str(soda.pickendyear),'-',num2str(soda.pickendmonth,'%02.f')]);
		else
			title(['b._ ', soda.overlay_name,strdep,' on reverse tracks',strdep1,' from ',num2str(soda.pickstartyear),'-',num2str(soda.pickstartmonth,'%02.f'),' to ',num2str(soda.pickendyear),'-',num2str(soda.pickendmonth,'%02.f')]);
		end
		if numel(soda.drifter.sta)>0
			legend(soda.drifter.sta,'location','bestoutside')
		else
			% legend(colname(2:end),'location','bestoutside')
			clear legname
			legname = strings;
			for j = 1:soda.ndrift
				legname(j) = ['sta',num2str(j,'%02.f')];
			end
			legend(legname,'location','bestoutside')
		end
		ax = gca;
		ax.XAxis.TickValues = datetime(1980,1,1)+calquarters(0:(2018-1980)*4);
		% ax.XAxis.TickLabelFormat = 'QQQ-yyyy';
		ax.XAxis.TickLabelFormat = 'MMM-yyyy';
		xlim([min(soda.interp.datetime) max(soda.interp.datetime)])
		% yticks([0 .5 1 1.5 2 2.5 3 3.5 4 4.5 5]);
		% yticklabels([0:.5:5]);
		% ylim([0.5 2]);
		%
		% ylim([nanmin2(soda.drifter.overlay_mean(:)) nanmax2(soda.drifter.overlay_mean(:))]);
		if strcmp(soda.overlay,'omara_K200') | ...
			strcmp(soda.overlay,'omara_K100') | ...
			strcmp(soda.overlay,'omara_KT200') | ...
			strcmp(soda.overlay,'omara_KT100') | ...
			strcmp(soda.overlay,'omara_KT')
				ylim([0 3.5]);
		elseif	strcmp(soda.overlay,'phtot_K200') | ...
			strcmp(soda.overlay,'phtot_K100') | ...
			strcmp(soda.overlay,'phtot_KT200') | ...
			strcmp(soda.overlay,'phtot_KT100') | ...
			strcmp(soda.overlay,'phtot_KT')
				ylim([7.3 8.1]);
		end
		%
		yticks('auto')
		yt=yticks;	% returns the current y-axis tick values as vector yt
		yticklabels(yt);
		grid on
		hold off
		%
		axes(ha(1))
		% subplot(2,1,2)
		hold on
		m_proj('mercator','lat',soda.maplat,'long',soda.maplon);
		if ~strcmp(soda.overlay,'none')
			clear X Y Z
			X = soda.overlay_lon;
			Y = soda.overlay_lat;
			Z = squeeze(nanmean2(soda.interp.overlay_data,3))';
			[C,h] = m_contourf(X,Y,Z,256);
			set(h,'LineColor','none')
			caxis([soda.overlay_min soda.overlay_max]);
			h=colorbar;
			h.Label.String=['Average ',soda.overlay_name,strdep];
			if strcmp(soda.overlay_cm,'balance') | strcmp(soda.overlay_cm,'-balance') | ...
					strcmp(soda.overlay_cm,'delta') | strcmp(soda.overlay_cm,'-delta') | ...
					strcmp(soda.overlay_cm,'curl') | strcmp(soda.overlay_cm,'-curl')
				if ~isnan(soda.overlay_pivot)
					cmocean2(soda.overlay_cm,'pivot',soda.overlay_pivot,'opacity',soda.overlay_opacity);
				else
					cmocean2(soda.overlay_cm,'opacity',soda.overlay_opacity);
				end
			else 
				cmocean2(soda.overlay_cm,'opacity',soda.overlay_opacity);
			end
		end
		m_gshhs_i('patch',[.9 .9 .9],'edgecolor',[.9 .9 .9]);		% intermediate res coast
		for i=1:soda.ndrift
			% initial or final position
			m_plot(soda.drifter.lon0(i)-360, soda.drifter.lat0(i),'o','linewidth',1,'markeredgecolor',ck(i,:),'markerfacecolor','k','MarkerSize',soda.markersize)
			% track line
			m_plot(soda.drifter.lon_mean(:,i)-360,soda.drifter.lat_mean(:,i),'color',ck(i,:),'linewidth',1.5)
		end
		if soda.forward_tracking
			title(['a._ Forward tracks of drifters',strdep1,' from_ ',num2str(soda.pickstartyear),'-',num2str(soda.pickstartmonth,'%02.f'),' to ',num2str(soda.pickendyear),'-',num2str(soda.pickendmonth,'%02.f')]);
		else
			title(['a._ Reverse tracks of drifters',strdep1,' from_ ',num2str(soda.pickstartyear),'-',num2str(soda.pickstartmonth,'%02.f'),' to ',num2str(soda.pickendyear),'-',num2str(soda.pickendmonth,'%02.f')]);
		end
		m_grid('box','off','linestyle','none','XaxisLocation','bottom','YaxisLocation','left','xtick',([-240 -230 -220 -210 -200 -190 -180 -170 -160 -150 -140 -130 -120 -110 -100 -90 -80 -70 -60]),'ytick',([-65 -60 -50 -40 -30 -20 -10 0 10 20 30 40 50 60 65]),'FontSize',8);	% N Pacific
		hold off
		%
		set(gcf, 'PaperPosition', [0 0 6 9])   
		if soda.forward_tracking
			print(gcf, [pwd '/png/Drifter_tracks_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_forward_',soda.ver,'_threshold.png'], '-dpng', '-r600' );   %save file as PNG w/ 300dpi
		else
			print(gcf, [pwd '/png/Drifter_tracks_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_reverse_',soda.ver,'_threshold.png'], '-dpng', '-r600' );   %save file as PNG w/ 300dpi
		end

		% - - -
		% box plot of overlay variable by station
		% ck=colormap(jet2(soda.ndrift));
		figure(1)
		hFig=gcf;
		clf(hFig);
		% T=bplot(soda.drifter.overlay_mean);
		hold on
		boxplot(soda.drifter.overlay_mean,'labels',soda.drifter.sta);
		% ploplot(soda.interp.datetime,ones(size(soda.interp.datetime))*1.2,'k-','linewidth',1)
		% plot(soda.interp.datetime,ones(size(soda.interp.datetime))*.95,'k-','linewidth',1.5)
		%
		if strcmp(soda.overlay,'omara_K200') | ...
			strcmp(soda.overlay,'omara_K100') | ...
			strcmp(soda.overlay,'omara_KT200') | ...
			strcmp(soda.overlay,'omara_KT100') | ...
			strcmp(soda.overlay,'omara_KT')
				hline(1.2,'k-','linewidth',1)
				hline(.95,'k-','linewidth',1.5)
		elseif	strcmp(soda.overlay,'phtot_K200') | ...
			strcmp(soda.overlay,'phtot_K100') | ...
			strcmp(soda.overlay,'phtot_KT200') | ...
			strcmp(soda.overlay,'phtot_KT100') | ...
			strcmp(soda.overlay,'phtot_KT')
				hline(7.75,'k-','linewidth',1)
				hline(7.4,'k-','linewidth',1.5)
		end
		%
		h = findobj(gcf,'tag','Outliers');
		set(h,'MarkerSize',.5)		
		% overlay mean values as black diamonds
		plot(nanmean2(soda.drifter.overlay_mean), 'dk')
		%
		% ylim([nanmin2(soda.drifter.overlay_mean(:)) nanmax2(soda.drifter.overlay_mean(:))]);
		if strcmp(soda.overlay,'omara_K200') | ...
			strcmp(soda.overlay,'omara_K100') | ...
			strcmp(soda.overlay,'omara_KT200') | ...
			strcmp(soda.overlay,'omara_KT100') | ...
			strcmp(soda.overlay,'omara_KT')
				ylim([0 3.5]);
		elseif	strcmp(soda.overlay,'phtot_K200') | ...
			strcmp(soda.overlay,'phtot_K100') | ...
			strcmp(soda.overlay,'phtot_KT200') | ...
			strcmp(soda.overlay,'phtot_KT100') | ...
			strcmp(soda.overlay,'phtot_KT')
				ylim([7.3 8.1]);
		end
		%
		% boxplot(soda.drifter.overlay_mean,'labels',soda.drifter.sta,'symbol','w+');
		if strcmp(soda.overlay,'omara_KT100') | strcmp(soda.overlay,'omcal_KT100') 
			hline(1,'k:','linewidth',1,'color',[.5 .5 .5]);
		elseif strcmp(soda.overlay,'omara_satdep') | strcmp(soda.overlay,'omcal_satdep') 
			hline(100,'k:','linewidth',1,'color',[.5 .5 .5]);
		end
		if soda.forward_tracking
			title(['c._ Boxplot of ',soda.overlay_name,strdep,' on forward tracks',strdep1,' from ',num2str(soda.pickstartyear),'-',num2str(soda.pickstartmonth,'%02.f'),' to ',num2str(soda.pickendyear),'-',num2str(soda.pickendmonth,'%02.f')]);
		else
			title(['c._ Boxplot of ',soda.overlay_name,strdep,' on reverse tracks',strdep1,' from ',num2str(soda.pickstartyear),'-',num2str(soda.pickstartmonth,'%02.f'),' to ',num2str(soda.pickendyear),'-',num2str(soda.pickendmonth,'%02.f')]);
		end
		ylabel([soda.overlay_name,strdep]);
		% grid on
		set(gca, 'YGrid', 'on', 'XGrid', 'off');
		hold off
		set(gcf, 'PaperPosition', [0 0 7.5 5])   
		if soda.forward_tracking
			print(gcf, [pwd '/png/Boxplot_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_forward_',soda.ver,'_threshold.png'], '-dpng', '-r600' );   %save file as PNG w/ 300dpi
		else
			print(gcf, [pwd '/png/Boxplot_',soda.overlay,'_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),num2str(soda.pickstartmonth,'%02.f'),'_',num2str(soda.pickendyear),num2str(soda.pickendmonth,'%02.f'),'_reverse_',soda.ver,'_threshold.png'], '-dpng', '-r600' );   %save file as PNG w/ 300dpi
		end

	end 		% if overlay data is omara | phtot

end 	% if soda.extraplots





















