
% - - - - -
% Trajectories of passive drifters using monthly velocity u and v from SODA3.12.2
% with the option to use either forward or reverse tracking
% - - - - -
% Greg Pelletier (gjpelletier@gmail.com)
% - - - - -

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
% 
% input of u and v data requires the following netcdf file: 
% [pwd '\nc\map1_velocity_2D_soda3_12_2_1980_2017_v20220421.nc']
% that was obtained by running the following scripts:
%
% soda_v20.m
% soda_write_netcdf_v20.m
%
% ----------
%
% overlay data requires the following netcdf file: 
% [pwd, '\nc\map1_2D_PlankTOM12_A_1980-2018_v20220419.nc']
% that was obtained by running the following scripts:
% 
% RECCAP_2021_3D_v42.m
% RECCAP_2021_2D_v44_write_netcdf.m
% 

% ----------
% ----------
% ----------
% START OF USER INPUTS

% folders containing matlab scripts for functions used in this main script
addpath(genpath('c:/matlab/greg/'))  		% cmocean2 modified from original cmocean at https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps?s_tid=srchtitle
addpath(genpath('c:/matlab/m_map/'))  		% functions for making maps in Mercator projection https://www.eoas.ubc.ca/~rich/map.html
addpath(genpath('c:/matlab/gif/'))  		% make animated gif https://www.mathworks.com/matlabcentral/fileexchange/63239-gif
addpath(genpath('c:/matlab/export_fig/'))  	% needed by gif https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig

% default random seed and directories
cd 'C:\data\soda3.12.2';
if not(isfolder('png'))
    mkdir('png')
end
if not(isfolder('ani'))
    mkdir('ani')
end
if not(isfolder('gif'))
    mkdir('gif')
end
if not(isfolder('csv'))
    mkdir('csv')
end
rng('default');	

tic
clc

clear soda;

soda.ver = 'v52';					% version number for output file names
soda.forward_tracking = true;		% drifter tracking forward in time (true), or reverse drifter tracking for exposure history (false)
soda.movmean_months = 3;			% number of months for moving average u and v (1= no smoothing, 3= 3-month moving average, etc)
soda.pickstartyear = 2017;			% the year when the drifter experiment starts			
soda.pickstartmonth = 1;			% the month to start the drifter experiment in the start year (e.g. 1=Jan, 2=Feb, etc)
soda.pickendyear = 2017;			% the year when the drifter experiment ends. for one year simulation input the end year with same value as the start year 
soda.pickendmonth = 12;				% the month to end the drifter experiment in the end year (e.g. 12=Dec, 11=Nov, etc)
soda.gifweekly = true;				% write weekly animation frames from g loop to gif 
soda.gifdelay = 0.25;				% delay between frames for animated gif (suggest 0.25 for 12-month, 0.125 for 5 year)
soda.gifpause = 4;					% scale factor for gifdelay to pause first and last frame (suggest 4 for 12-month, 8 for 5 year)

soda.whichmap = 'PacificNW';		% Map extent to use for animated gif frames (add more options as needed in the code for making the animation) 
									% 'PacificN' = entire North Pacific from 0-65N and 128E-110W
									% 'PacificNE' = region of pteropod stations in the Northeast Pacific from 30N-61N and 165W-120W
									% 'PacificNE2' = region of 135W transect line in the Northeast Pacific from 10N-61N and 175E-110W
									% 'PacificNE3' = region of 153W transect line in the Northeast Pacific
									% 'PacificNW' = region of Japanese pteropod stations from 125E-180E and 10N-60N
									% 'Osborne' = region of Osborne et al (1989) paper on trajectories for 12 months from start off SW coast of Japan
									% NOTE: any of the choices above may be used for the map extent, 
									% but if any other map extent is needed instead
									% it can be customized in the next if/elseif whichmap block of code below

soda.whichdrifters = 'Japan';		% Initial position of drifters (add more choices as needed in next if/elseif block of code):  
									% 'pteropods' = L-shaped line of pteropod stations in the NE Pacific,
									% '153W' = north south at 153W from 58N-20N,
									% '135W' = north south at 135W from 55N-20N,
									% '33N' = east-west at 33N from 140E to 118W (latitude of the Kuroshio Extension)
									% '50N' = east-west at 50N from 155E to 128W 
									% 'Japan' = L-shaped line surrounding Japanese pteropod stations 25N from 137.5E to 160E and 160E from 25N to 50N
									% 'Kuroshio' = Kuroshio Current transport transect line near the SE coast of Japan in the Kuroshio Current
									% 'Oyashio' = Oyashio Current transport transect line off of East Kamchatka
									% 'KE' = Kuroshio Extension transport transect at lon 160E from 28N-409N
									% NOTE: any of the choices above may be used for the location of initial drifter positions, 
									% but if any other drifter initial positions are needed instead
									% they can be customized in the next if/elseif whichdrifters block of code below

soda.drifter_width = 0.5;			% east-west width of each drifter node (degrees)
soda.drifter_height = 0.5;			% north-south height of each drifter node (degrees)
soda.drifter_particles = 25;		% number of drifter particles at each node to be initially spaced randomly between width and height of each node

soda.markeredgecolor = 'jet';		% color for the outline of the drifter markers on the animation frames 
									% 	('jet'=jet colormap, otherwise any matlab color code such as 'k'=black, 'm'=magenta, 'c'=cyan, 'y'=yellow, 'w'=white, 'b'= blue, 'r'=red, etc) 
soda.markerfacecolor = 'k';			% color for the face of the drifter markers on the animation frames 
									% 	('jet'=jet colormap, otherwise any matlab color code such as 'k'=black, 'm'=magenta, 'c'=cyan, 'y'=yellow, 'w'=white, 'b'= blue, 'r'=red, etc) 
soda.markersize = 4;				% size of the marker for the drifter trajectory animation (points, e.g. matlab default size is 6 points)
soda.savepng = false;				% save png file for each animation frame

soda.overlay = 'omara_KT100';		% Overlay data layer to show as background on the animated gif map frames:
									% 'none' = do not overlay a data layer on the animation map frames
									% 'omara_satdep' = overlay depth to aragonite saturation horizon on the animation frames
									% 'omara_KT100' = overlay 0-100m depth-integrated omega aragonite
									% 'omcal_satdep' = overlay depth to calcite saturation horizon on the animation frames
									% 'omcal_KT100' = overlay 0-100m depth-integrated omega calcite
									% 'phtot_KT100' = overlay 0-100m depth-integrated pH (total scale)
									% 'talk_KT100' = overlay 0-100m depth-integrated total alkalinity
									% 'dic_KT100' = overlay 0-100m depth-integrated DIC
									% 'pco2_KT100' = overlay 0-100m depth-integrated pCO2
									% 'revelle_KT100' = overlay 0-100m depth-integrated Revelle Factor
									% 'temp_KT100' = overlay 0-100m depth-integrated temperature
									% 'sal_KT100' = overlay 0-100m depth-integrated salinity
									% 'po4_KT100' = overlay 0-100m depth-integrated total inorganic phosphorus
									% 'si_KT100' = overlay 0-100m depth-integrated total inorganic silica
									% 'o2_KT100' = overlay 0-100m depth-integrated dissolved oxygen
									% 'aou_KT100' = overlay 0-100m depth-integrated apparent oxygen utilization
									% 'epcalc100' = export flux of CaCO3 at 100m 
									% 'intphyc' = water column depth-integrated total phytoplankton carbon
									% 'intzooc' = water column depth-integrated total zooplankton carbon
									% 'intpp' = water column depth-integrated net primary production
									% NOTE: overlay data must be lon x lat x time and have same time origin of 1980 and be a continuous sequence of monthly data

soda.overlay_quiver = true;			% overlay the current vector quiver arrows on the animation frames
soda.quiverwidth = .1;			 	% line width of quiver arrows 
soda.quiverscale = 5;			 	% scale factor for quiver length
soda.quivercolor = 'k';			 	% color for quiver arrows 

soda.depthlabel = '0-100m';			% label for plots to label the depth or range of depths that is used for the velocity vectors

% calculated values using the inputs above
soda.picknyears = soda.pickendyear-soda.pickstartyear+1;	% intermediate calc to find the number of months for duration of the experiment (picknmonths)
soda.picknmonths = (soda.picknyears-2) * 12 + (12-soda.pickstartmonth+1) + soda.pickendmonth;

soda.nsubsteps = 30*24;				% needs to be a multiple of 4 and 30, number of increments to divide each time step (nsubsteps=720 is approx 1 hour substep increment (30/720*24)=1hr

% Edit the following as needed to change the map extent of the animation frames 
% select the map extent for the animation frames (add more options as needed)
if strcmp(soda.whichmap,'PacificN') 
	soda.maplat = [0 65];					% south and north latitudes of map extent (degN)
	soda.maplon = [(120-360) -110];		% west and east longitudes of map extent (-degW)
elseif strcmp(soda.whichmap,'PacificNE') 
	soda.maplat = [30 61];					% south and north latitudes of map extent (degN)
	soda.maplon = [-175 -120];				% west and east longitudes of map extent (-degW)
elseif strcmp(soda.whichmap,'PacificNE2') 
	soda.maplat = [10 61];					% south and north latitudes of map extent (degN)
	soda.maplon = [-185 -110];				% west and east longitudes of map extent (-degW)
elseif strcmp(soda.whichmap,'PacificNE3') 
	soda.maplat = [10 61];					% south and north latitudes of map extent (degN)
	soda.maplon = [-195 -120];				% west and east longitudes of map extent (-degW)
elseif strcmp(soda.whichmap,'PacificNW') 
	soda.maplat = [10 60];					% south and north latitudes of map extent (degN)
	soda.maplon = [(120-360) -170];		% west and east longitudes of map extent (-degW)
elseif strcmp(soda.whichmap,'Osborne') 
	soda.maplat = [10 60];					% south and north latitudes of map extent (degN)
	soda.maplon = [(120-360) -155];		% west and east longitudes of map extent (-degW)
end

% Edit the following as needed to change the initial/final position of the drifter node centers 
% or add more expermients with drifter nodes in different locations. Each drifter node will have 
% the number of soda.drifter_particles specified above randomly spaced within the specifed width and height above
if strcmp(soda.whichdrifters,'pteropods') 
	soda.ndrift = 51;					% number of drifters
	% initialize vectors for drifter lat and lon
	soda.d01.lon = []; soda.d01.lat = []; 
	soda.d01.lon(1:soda.picknmonths,1:soda.ndrift) = NaN;
	soda.d01.lat(1:soda.picknmonths,1:soda.ndrift) = NaN;
	% initialize index relative to map1 for drifter lat and lon
	soda.d01.map1_lon_idx = []; soda.d01.map1_lat_idx = []; 
	soda.d01.map1_lon_idx(1:soda.picknmonths,1:soda.ndrift) = NaN;
	soda.d01.map1_lat_idx(1:soda.picknmonths,1:soda.ndrift) = NaN;
	% v12_PacificNE using Nina's line of stations
	if soda.forward_tracking
		% initial position of drifters
		soda.d01.lat(1,1:36) = 50;
		soda.d01.lon(1,1) = 360-135;
		for i = 2:36
			soda.d01.lon(1,i) = soda.d01.lon(1,i-1) - 0.5;
		end
		soda.d01.lon(1,37:51) = 360-153;
		soda.d01.lat(1,37) = 50;
		for i = 38:51
			soda.d01.lat(1,i) = soda.d01.lat(1,i-1) - 0.5;
		end
	else		% reverse drifter tracking
		% final position of drifters
		soda.d01.lat(soda.picknmonths,1:36) = 50;
		soda.d01.lon(soda.picknmonths,1) = 360-135;
		for i = 2:36
			soda.d01.lon(soda.picknmonths,i) = soda.d01.lon(soda.picknmonths,i-1) - 0.5;
		end
		soda.d01.lon(soda.picknmonths,37:51) = 360-153;
		soda.d01.lat(soda.picknmonths,37) = 50;
		for i = 38:51
			soda.d01.lat(soda.picknmonths,i) = soda.d01.lat(soda.picknmonths,i-1) - 0.5;
		end
	end
elseif strcmp(soda.whichdrifters,'153W') 
	% soda.ndrift = 77;					% number of drifters
	soda.ndrift = 75;					% number of drifters
	% initialize vectors for drifter lat and lon
	soda.d01.lon = []; soda.d01.lat = []; 
	soda.d01.lon(1:soda.picknmonths,1:soda.ndrift) = NaN;
	soda.d01.lat(1:soda.picknmonths,1:soda.ndrift) = NaN;
	% initialize index relative to map1 for drifter lat and lon
	soda.d01.map1_lon_idx = []; soda.d01.map1_lat_idx = []; 
	soda.d01.map1_lon_idx(1:soda.picknmonths,1:soda.ndrift) = NaN;
	soda.d01.map1_lat_idx(1:soda.picknmonths,1:soda.ndrift) = NaN;
	%% PacificN_v14_153W_from_58N_to_20N
	if soda.forward_tracking
		% initial position of drifters
		soda.d01.lon(1,1:soda.ndrift) = 360-153;
		% soda.d01.lat(1,1) = 58;
		soda.d01.lat(1,1) = 57;
		for i = 2:soda.ndrift
			soda.d01.lat(1,i) = soda.d01.lat(1,i-1) - 0.5;
		end
	else		% reverse drifter tracking
		% final position of drifters
		soda.d01.lon(soda.picknmonths,1:soda.ndrift) = 360-153;
		% soda.d01.lat(soda.picknmonths,1) = 58;
		soda.d01.lat(soda.picknmonths,1) = 57;
		for i = 2:soda.ndrift
			soda.d01.lat(soda.picknmonths,i) = soda.d01.lat(soda.picknmonths,i-1) - 0.5;
		end
	end
elseif strcmp(soda.whichdrifters,'135W') 
	soda.ndrift = 71;					% number of drifters
	% initialize vectors for drifter lat and lon
	soda.d01.lon = []; soda.d01.lat = []; 
	soda.d01.lon(1:soda.picknmonths,1:soda.ndrift) = NaN;
	soda.d01.lat(1:soda.picknmonths,1:soda.ndrift) = NaN;
	% initialize index relative to map1 for drifter lat and lon
	soda.d01.map1_lon_idx = []; soda.d01.map1_lat_idx = []; 
	soda.d01.map1_lon_idx(1:soda.picknmonths,1:soda.ndrift) = NaN;
	soda.d01.map1_lat_idx(1:soda.picknmonths,1:soda.ndrift) = NaN;
	%% PacificN_v14_153W_from_58N_to_20N
	if soda.forward_tracking
		% initial position of drifters
		soda.d01.lon(1,1:soda.ndrift) = 360-135;
		soda.d01.lat(1,1) = 55;
		for i = 2:soda.ndrift
			soda.d01.lat(1,i) = soda.d01.lat(1,i-1) - 0.5;
		end
	else		% reverse drifter tracking
		% final position of drifters
		soda.d01.lon(soda.picknmonths,1:soda.ndrift) = 360-135;
		% soda.d01.lat(soda.picknmonths,1) = 58;
		soda.d01.lat(soda.picknmonths,1) = 55;
		for i = 2:soda.ndrift
			soda.d01.lat(soda.picknmonths,i) = soda.d01.lat(soda.picknmonths,i-1) - 0.5;
		end
	end
elseif strcmp(soda.whichdrifters,'33N') 
	soda.ndrift = 205;					% number of drifters
	% initialize vectors for drifter lat and lon
	soda.d01.lon = []; soda.d01.lat = []; 
	soda.d01.lon(1:soda.picknmonths,1:soda.ndrift) = NaN;
	soda.d01.lat(1:soda.picknmonths,1:soda.ndrift) = NaN;
	% initialize index relative to map1 for drifter lat and lon
	soda.d01.map1_lon_idx = []; soda.d01.map1_lat_idx = []; 
	soda.d01.map1_lon_idx(1:soda.picknmonths,1:soda.ndrift) = NaN;
	soda.d01.map1_lat_idx(1:soda.picknmonths,1:soda.ndrift) = NaN;
	%% PacificN_v14_153W_from_58N_to_20N
	if soda.forward_tracking
		% initial position of drifters
		soda.d01.lat(1,1:soda.ndrift) = 33;
		soda.d01.lon(1,1) = 360-118;
		for i = 2:soda.ndrift
			soda.d01.lon(1,i) = soda.d01.lon(1,i-1) - 0.5;
		end
	else		% reverse drifter tracking
		% final position of drifters
		soda.d01.lat(soda.picknmonths,1:soda.ndrift) = 33;
		soda.d01.lon(soda.picknmonths,1) = 360-118;
		for i = 2:soda.ndrift
			soda.d01.lon(soda.picknmonths,i) = soda.d01.lon(soda.picknmonths,i-1) - 0.5;
		end
	end
elseif strcmp(soda.whichdrifters,'50N') 
	soda.ndrift = 155;			% number of drifters
	% initialize vectors for drifter lat and lon
	soda.d01.lon = []; soda.d01.lat = []; 
	soda.d01.lon(1:soda.picknmonths,1:soda.ndrift) = NaN;
	soda.d01.lat(1:soda.picknmonths,1:soda.ndrift) = NaN;
	% initialize index relative to map1 for drifter lat and lon
	soda.d01.map1_lon_idx = []; soda.d01.map1_lat_idx = []; 
	soda.d01.map1_lon_idx(1:soda.picknmonths,1:soda.ndrift) = NaN;
	soda.d01.map1_lat_idx(1:soda.picknmonths,1:soda.ndrift) = NaN;
	%% PacificN_v14_153W_from_58N_to_20N
	if soda.forward_tracking
		% initial position of drifters
		soda.d01.lat(1,1:soda.ndrift) = 50;
		soda.d01.lon(1,1) = 360-128;
		for i = 2:soda.ndrift
			soda.d01.lon(1,i) = soda.d01.lon(1,i-1) - 0.5;
		end
	else		% reverse drifter tracking
		% final position of drifters
		soda.d01.lat(soda.picknmonths,1:soda.ndrift) = 50;
		soda.d01.lon(soda.picknmonths,1) = 360-128;
		for i = 2:soda.ndrift
			soda.d01.lon(soda.picknmonths,i) = soda.d01.lon(soda.picknmonths,i-1) - 0.5;
		end
	end
elseif strcmp(soda.whichdrifters,'Japan') 
	soda.ndrift = 96;					% number of drifters
	% initialize vectors for initial drifter nodes
	soda.d01.lon = []; soda.d01.lat = []; 
	soda.d01.lon(1:soda.picknmonths,1:soda.ndrift) = NaN;
	soda.d01.lat(1:soda.picknmonths,1:soda.ndrift) = NaN;
	% initialize index relative to map1 for drifter lat and lon
	soda.d01.map1_lon_idx = []; soda.d01.map1_lat_idx = []; 
	soda.d01.map1_lon_idx(1:soda.picknmonths,1:soda.ndrift) = NaN;
	soda.d01.map1_lat_idx(1:soda.picknmonths,1:soda.ndrift) = NaN;
	% L-shaped line surrounding Japanese pteropod stations 25N from 137.5E to 160E and 160E from 25N to 50N
	% initial position of each drifter node
	if soda.forward_tracking
		% initial position of drifters
		soda.d01.lon(1,1:51) = 160;
		soda.d01.lat(1,1) = 50;
		for i = 2:51
			soda.d01.lat(1,i) = soda.d01.lat(1,i-1) - 0.5;
		end
		soda.d01.lat(1,52:96) = 25;
		soda.d01.lon(1,52) = 159.5;
		for i = 53:96
			soda.d01.lon(1,i) = soda.d01.lon(1,i-1) - 0.5;
		end
	else		% reverse drifter tracking
		soda.d01.lon(soda.picknmonths,1:51) = 160;
		soda.d01.lat(soda.picknmonths,1) = 50;
		for i = 2:51
			soda.d01.lat(soda.picknmonths,i) = soda.d01.lat(soda.picknmonths,i-1) - 0.5;
		end
		soda.d01.lat(soda.picknmonths,52:96) = 25;
		soda.d01.lon(soda.picknmonths,52) = 159.5;
		for i = 53:96
			soda.d01.lon(soda.picknmonths,i) = soda.d01.lon(soda.picknmonths,i-1) - 0.5;
		end
	end
elseif strcmp(soda.whichdrifters,'Kuroshio')	% Kuroshio Current transport transect off the SW coast of Japan
	soda.ndrift = 7;					% number of drifters
	% soda.x_kuroshio.lat = [36	35.5	35	34.5	34 33.5 33];	% extend to the southwest 2 more nodes
	% soda.x_kuroshio.lon = [141	141.375	141.75	142.125	142.5 142.875 143.25];
	% initialize vectors for drifter lat and lon
	soda.d01.lon = []; soda.d01.lat = []; 
	soda.d01.lon(1:soda.picknmonths,1:soda.ndrift) = NaN;
	soda.d01.lat(1:soda.picknmonths,1:soda.ndrift) = NaN;
	% initialize index relative to map1 for drifter lat and lon
	soda.d01.map1_lon_idx = []; soda.d01.map1_lat_idx = []; 
	soda.d01.map1_lon_idx(1:soda.picknmonths,1:soda.ndrift) = NaN;
	soda.d01.map1_lat_idx(1:soda.picknmonths,1:soda.ndrift) = NaN;
	if soda.forward_tracking
		% initial position of drifters
		soda.d01.lon(1,:) = [141	141.375	141.75	142.125	142.5 142.875 143.25];
		soda.d01.lat(1,:) = [36	35.5	35	34.5	34 33.5 33];
	else		% reverse drifter tracking
		% final position of drifters
		soda.d01.lon(soda.picknmonths,:) = [141	141.375	141.75	142.125	142.5 142.875 143.25];
		soda.d01.lat(soda.picknmonths,:) = [36	35.5	35	34.5	34 33.5 33];
	end
elseif strcmp(soda.whichdrifters,'Oyashio') 	% Oyashio current transport transect at East Kamchatka
	soda.ndrift = 6;					% number of drifters
	% soda.x_kamchatka.lat = [54	53.6875	53.375	53.0625	52.75	52.4375];
	% soda.x_kamchatka.lon = [160	160.625	161.25	161.875	162.5	163.125];
	% initialize vectors for drifter lat and lon
	soda.d01.lon = []; soda.d01.lat = []; 
	soda.d01.lon(1:soda.picknmonths,1:soda.ndrift) = NaN;
	soda.d01.lat(1:soda.picknmonths,1:soda.ndrift) = NaN;
	% initialize index relative to map1 for drifter lat and lon
	soda.d01.map1_lon_idx = []; soda.d01.map1_lat_idx = []; 
	soda.d01.map1_lon_idx(1:soda.picknmonths,1:soda.ndrift) = NaN;
	soda.d01.map1_lat_idx(1:soda.picknmonths,1:soda.ndrift) = NaN;
	if soda.forward_tracking
		% initial position of drifters
		soda.d01.lon(1,:) = [160	160.625	161.25	161.875	162.5	163.125];
		soda.d01.lat(1,:) = [54	53.6875	53.375	53.0625	52.75	52.4375]
	else		% reverse drifter tracking
		% final position of drifters
		soda.d01.lon(soda.picknmonths,:) = [160	160.625	161.25	161.875	162.5	163.125];
		soda.d01.lat(soda.picknmonths,:) = [54	53.6875	53.375	53.0625	52.75	52.4375];
	end
elseif strcmp(soda.whichdrifters,'KE')		% Kuroshio Extension transport transect 
	soda.ndrift = 43;					% number of drifters
	% soda.x_KE.lat = [49:-0.5:28];
	% soda.x_KE.lon(1:numel(soda.x_KE.lat)) = 160;
	% initialize vectors for drifter lat and lon
	soda.d01.lon = []; soda.d01.lat = []; 
	soda.d01.lon(1:soda.picknmonths,1:soda.ndrift) = NaN;
	soda.d01.lat(1:soda.picknmonths,1:soda.ndrift) = NaN;
	% initialize index relative to map1 for drifter lat and lon
	soda.d01.map1_lon_idx = []; soda.d01.map1_lat_idx = []; 
	soda.d01.map1_lon_idx(1:soda.picknmonths,1:soda.ndrift) = NaN;
	soda.d01.map1_lat_idx(1:soda.picknmonths,1:soda.ndrift) = NaN;
	if soda.forward_tracking
		% initial position of drifters
		soda.d01.lon(1,:) = 160;
		soda.d01.lat(1,:) = [49:-0.5:28];
	else		% reverse drifter tracking
		% final position of drifters
		soda.d01.lon(soda.picknmonths,:) = 160;
		soda.d01.lat(soda.picknmonths,:) = [49:-0.5:28];
	end
end

% read netcdf file of soda lon, lat, time, u, v 
soda.map1_fn = [pwd '\nc\map1_velocity_2D_soda3_12_2_1980_2017_v20220421.nc'];	% edit as needed to correspond with the path to this nc file
soda.map1.lon = ncread(soda.map1_fn,'lon');
soda.map1.lat = ncread(soda.map1_fn,'lat');
soda.time = ncread(soda.map1_fn,'time');
soda.map1.u_d01 = ncread(soda.map1_fn,'u_KT100');	% KT100 is vertically integrated from 0-100m. KT200, KT500, and KT1000 are also available to select from this netcdf file
soda.map1.v_d01 = ncread(soda.map1_fn,'v_KT100');
soda.map1.lonW = soda.map1.lon - 360;

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

% meshgrid x and y arrays corresponding to the input u and v arrays, used by interp2 and m_map quiver plots
% x_mesh must use lon degW
[soda.map1.x_mesh,soda.map1.y_mesh]=meshgrid([115-360:.5:265-360],[0:.5:66]);		% this is the map1 extent that was extracted from SODA3.12.2 for the default 2D lon x lat x time arrays of vertically integrated u and v
																					
% calculate movmean u and v moving average
if soda.movmean_months > 1
	for i = 1:numel(soda.map1.lon)
		for j = 1:numel(soda.map1.lat)
			soda.map1.u_d01(i,j,:) = movmean(soda.map1.u_d01(i,j,:),soda.movmean_months);
			soda.map1.v_d01(i,j,:) = movmean(soda.map1.v_d01(i,j,:),soda.movmean_months);
		end
	end
end

% map overlay data for animation frames
if ~strcmp(soda.overlay,'none')
	soda.overlay_fn = [pwd, '\nc\map1_2D_PlankTOM12_A_1980-2018_v20220419.nc'];		% edit as needed to correspond with the path to this nc file
end
if strcmp(soda.overlay,'omara_satdep')
	soda.overlay_data = ncread(soda.overlay_fn,'omara_satdep_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = '\Omega_A saturation depth (m)';
	soda.overlay_pivot = 100;
	soda.overlay_min = 0;
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-delta';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'omcal_satdep')
	soda.overlay_data = ncread(soda.overlay_fn,'omcal_satdep_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = '\Omega_C saturation depth (m)';
	soda.overlay_pivot = 100;
	soda.overlay_min = 0;
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-delta';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'omara_KT100')
	soda.overlay_data = ncread(soda.overlay_fn,'omara_KT100_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = '\Omega_A';
	soda.overlay_pivot = 1;
	soda.overlay_min = 0;
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-delta';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'omcal_KT100')
	soda.overlay_data = ncread(soda.overlay_fn,'omcal_KT100_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = '\Omega_C';
	soda.overlay_pivot = 1;
	soda.overlay_min = 0;
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-delta';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'phtot_KT100')
	soda.overlay_data = ncread(soda.overlay_fn,'phtot_KT100_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = 'pH';
	soda.overlay_pivot = 7.6;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-delta';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'talk_KT100')
	soda.overlay_data = ncread(soda.overlay_fn,'talk_KT100_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = 'TA (\mumol kg^-^1)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'dic_KT100')
	soda.overlay_data = ncread(soda.overlay_fn,'dic_KT100_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = 'DIC (\mumol kg^-^1)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'pco2_KT100')
	soda.overlay_data = ncread(soda.overlay_fn,'pco2_KT100_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = 'pCO_2 (\muatm)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'revelle_KT100')
	soda.overlay_data = ncread(soda.overlay_fn,'revelle_KT100_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = 'Revelle Factor';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'temp_KT100')
	soda.overlay_data = ncread(soda.overlay_fn,'temp_KT100_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = 'Temperature (\circC)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'sal_KT100')
	soda.overlay_data = ncread(soda.overlay_fn,'sal_KT100_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = 'Salinity (psu)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'po4_KT100')
	soda.overlay_data = ncread(soda.overlay_fn,'po4_KT100_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = 'PO4 (\mumol kg^-^1)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'si_KT100')
	soda.overlay_data = ncread(soda.overlay_fn,'si_KT100_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = 'Si (\mumol kg^-^1)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'o2_KT100')
	soda.overlay_data = ncread(soda.overlay_fn,'o2_KT100_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = 'O_2 (\mumol kg^-^1)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'aou_KT100')
	soda.overlay_data = ncread(soda.overlay_fn,'aou_KT100_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = 'AOU (\mumol kg^-^1)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'epcalc100')
	soda.overlay_data = ncread(soda.overlay_fn,'epcalc100_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = 'epcalc100 (mol m^-^2 sec^-^1)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'intphyc')
	soda.overlay_data = ncread(soda.overlay_fn,'intphyc_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = 'intphyc (mol m^-^2)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'intzooc')
	soda.overlay_data = ncread(soda.overlay_fn,'intzooc_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = 'intzooc (mol m^-^2)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
elseif strcmp(soda.overlay,'intpp')
	soda.overlay_data = ncread(soda.overlay_fn,'intpp_a');  
	soda.overlay_data(soda.overlay_data==0)=NaN;
	soda.overlay_lon = ncread(soda.overlay_fn,'lon')-360;	% convert to degW  
	soda.overlay_lat = ncread(soda.overlay_fn,'lat');  
	soda.overlay_name = 'intpp (mol m^-^2 sec^-^1)';
	soda.overlay_pivot = NaN;
	soda.overlay_min = prctile(soda.overlay_data(:),10);
	soda.overlay_max = prctile(soda.overlay_data(:),90);
	soda.overlay_cm = '-haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
end
% mask non-water cells as NaN
if ~strcmp(soda.overlay,'none')
	soda.overlay_fn = [pwd, '\nc\map1_2D_PlankTOM12_A_1980-2018_v20220419.nc'];
	soda.overlay_mask = ncread(soda.overlay_fn,'mask');  
	for i = 1:size(soda.overlay_data,3)
		soda.datamask = squeeze(soda.overlay_data(:,:,i));
		soda.datamask(soda.overlay_mask==0)=NaN;
		soda.overlay_data(:,:,i) = soda.datamask;
	end
end

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
	disp(['Program terminated because soda.startyear must be 1980 to overlay the aragonite saturation depths correctly'])
	exit()
end

% END OF USER INPUTS
% ----------
% ----------
% ----------

% initialize arrays for particles within initial drifter nodes
soda.d01.lon_kk = []; soda.d01.lat_kk = []; 
soda.d01.lon_kk(1:soda.picknmonths,1:soda.ndrift,1:soda.drifter_particles) = NaN;
soda.d01.lat_kk(1:soda.picknmonths,1:soda.ndrift,1:soda.drifter_particles) = NaN;
soda.d01.map1_lon_idx_kk = []; soda.d01.map1_lat_idx_kk = []; 
soda.d01.map1_lon_idx_kk(1:soda.picknmonths,1:soda.ndrift,1:soda.drifter_particles) = NaN;
soda.d01.map1_lat_idx_kk(1:soda.picknmonths,1:soda.ndrift,1:soda.drifter_particles) = NaN;
% initial position of particles randomly placed within each initial drifter node
if soda.forward_tracking
	for i = 1:soda.ndrift
		lon_a = soda.d01.lon(1,i) - soda.drifter_width / 2;		% east side of drifter node
		lon_b = soda.d01.lon(1,i) + soda.drifter_width / 2;		% west side of drifter node
		lat_a = soda.d01.lat(1,i) - soda.drifter_height / 2;	% south side of drifter node
		lat_b = soda.d01.lat(1,i) + soda.drifter_height / 2;	% north side of drifter node
		for j = 1:soda.drifter_particles
			soda.d01.lon_kk(1,i,j) = lon_a + (lon_b - lon_a) * rand();
			soda.d01.lat_kk(1,i,j) = lat_a + (lat_b - lat_a) * rand();
		end
	end
else
	for i = 1:soda.ndrift
		lon_a = soda.d01.lon(soda.picknmonths,i) - soda.drifter_width / 2;		% east side of drifter node
		lon_b = soda.d01.lon(soda.picknmonths,i) + soda.drifter_width / 2;		% west side of drifter node
		lat_a = soda.d01.lat(soda.picknmonths,i) - soda.drifter_height / 2;	% south side of drifter node
		lat_b = soda.d01.lat(soda.picknmonths,i) + soda.drifter_height / 2;	% north side of drifter node
		for j = 1:soda.drifter_particles
			soda.d01.lon_kk(soda.picknmonths,i,j) = lon_a + (lon_b - lon_a) * rand();
			soda.d01.lat_kk(soda.picknmonths,i,j) = lat_a + (lat_b - lat_a) * rand();
		end
	end
end
% reshape 3d particle arrays to 2d
soda.d01.lon_kk_rsp = reshape(soda.d01.lon_kk,soda.picknmonths,[]);
soda.d01.lat_kk_rsp = reshape(soda.d01.lat_kk,soda.picknmonths,[]);

% optional weekly arrays for save weekly frames
if soda.gifweekly
	soda.d01.lon_weekly = []; soda.d01.lat_weekly = []; 
	soda.d01.lon_weekly(1:soda.picknmonths*4,1:soda.ndrift*soda.drifter_particles) = NaN;
	soda.d01.lat_weekly(1:soda.picknmonths*4,1:soda.ndrift*soda.drifter_particles) = NaN;
	if soda.forward_tracking
		% initial position of drifters
		soda.d01.lat_weekly(1,:) = soda.d01.lat_kk_rsp(1,:);
		soda.d01.lat_weekly(2,:) = soda.d01.lat_kk_rsp(1,:);
		soda.d01.lat_weekly(3,:) = soda.d01.lat_kk_rsp(1,:);
		soda.d01.lat_weekly(4,:) = soda.d01.lat_kk_rsp(1,:);
		soda.d01.lon_weekly(1,:) = soda.d01.lon_kk_rsp(1,:);
		soda.d01.lon_weekly(2,:) = soda.d01.lon_kk_rsp(1,:);
		soda.d01.lon_weekly(3,:) = soda.d01.lon_kk_rsp(1,:);
		soda.d01.lon_weekly(4,:) = soda.d01.lon_kk_rsp(1,:);
	else		% reverse drifter tracking
		% final position of drifters
		soda.d01.lat_weekly(soda.picknmonths*4,:) = soda.d01.lat_kk_rsp(soda.picknmonths,:);
		soda.d01.lat_weekly(soda.picknmonths*4-1,:) = soda.d01.lat_kk_rsp(soda.picknmonths,:);
		soda.d01.lat_weekly(soda.picknmonths*4-2,:) = soda.d01.lat_kk_rsp(soda.picknmonths,:);
		soda.d01.lat_weekly(soda.picknmonths*4-3,:) = soda.d01.lat_kk_rsp(soda.picknmonths,:);
		soda.d01.lon_weekly(soda.picknmonths*4,:) = soda.d01.lon_kk_rsp(soda.picknmonths,:);
		soda.d01.lon_weekly(soda.picknmonths*4-1,:) = soda.d01.lon_kk_rsp(soda.picknmonths,:);
		soda.d01.lon_weekly(soda.picknmonths*4-2,:) = soda.d01.lon_kk_rsp(soda.picknmonths,:);
		soda.d01.lon_weekly(soda.picknmonths*4-3,:) = soda.d01.lon_kk_rsp(soda.picknmonths,:);
	end
end

% - - - - -
% loop to calculate the drifter positions over time
% 

soda.d01.u_gg = [];
soda.d01.v_gg = [];
soda.d01.uv_gg = [];
soda.d01.lon_gg = [];
soda.d01.lat_gg = [];
soda.d01.u_gg(1:soda.picknmonths*soda.nsubsteps,1:soda.ndrift*soda.drifter_particles) = NaN;
soda.d01.v_gg(1:soda.picknmonths*soda.nsubsteps,1:soda.ndrift*soda.drifter_particles) = NaN;
soda.d01.uv_gg(1:soda.picknmonths*soda.nsubsteps,1:soda.ndrift*soda.drifter_particles) = NaN;
soda.d01.lon_gg(1:soda.picknmonths*soda.nsubsteps,1:soda.ndrift*soda.drifter_particles) = NaN;
soda.d01.lat_gg(1:soda.picknmonths*soda.nsubsteps,1:soda.ndrift*soda.drifter_particles) = NaN;

xq = [];
yq = [];
uq = [];
vq = [];
uvq = [];
xq(1:soda.ndrift*soda.drifter_particles,1) = NaN;
yq(1:soda.ndrift*soda.drifter_particles,1) = NaN;
uq(1:soda.ndrift*soda.drifter_particles,1) = NaN;
vq(1:soda.ndrift*soda.drifter_particles,1) = NaN;
uvq(1:soda.ndrift*soda.drifter_particles,1) = NaN;

soda.d01.u = [];
soda.d01.v = [];
soda.d01.uv = [];
soda.d01.dx = [];
soda.d01.dy = [];
soda.d01.dlon = [];
soda.d01.dlat = [];
soda.d01.u(1:soda.ndrift*soda.drifter_particles,1) = NaN;
soda.d01.v(1:soda.ndrift*soda.drifter_particles,1) = NaN;
soda.d01.uv(1:soda.ndrift*soda.drifter_particles,1) = NaN;
soda.d01.dx(1:soda.ndrift*soda.drifter_particles,1) = NaN;
soda.d01.dy(1:soda.ndrift*soda.drifter_particles,1) = NaN;
soda.d01.dlon(1:soda.ndrift*soda.drifter_particles,1) = NaN;
soda.d01.dlat(1:soda.ndrift*soda.drifter_particles,1) = NaN;

soda.d01.lon_g = [];
soda.d01.lat_g = [];
soda.d01.lon_g(1:soda.nsubsteps+1,1:soda.ndrift*soda.drifter_particles) = NaN;
soda.d01.lat_g(1:soda.nsubsteps+1,1:soda.ndrift*soda.drifter_particles) = NaN;

ggcount = 0;
if soda.forward_tracking		% drifter tracking forward in time from specified initial position

	% iweek=1;
	iweek=1+4;
	gcount=0;

	iskip = (soda.pickstartyear - soda.startyear) * 12 + soda.pickstartmonth - 1;

	for i = 1:soda.picknmonths-1
	disp(['Calculating drifter position for month ', num2str(i), ' ...'])

		j = i + iskip;	 % index for pickyear to select u and v in the selected years
		soda.timestep = (soda.time(j+1) - soda.time(j)) * 86400;		% time step in seconds

		for g = 1:soda.nsubsteps	% divide the timestep into nsubsteps increments

			gcount = gcount + 1;	% count substeps for this month i
			ggcount = ggcount + 1;	% cumulative across all months and substeps

			if g == 1
				soda.d01.lon_g(g,:) = soda.d01.lon_kk_rsp(i,:);
				soda.d01.lat_g(g,:) = soda.d01.lat_kk_rsp(i,:);
			end

			% use interp2 to interpolate uq and vq for all current drifter positions xq and yq
			x = soda.map1.x_mesh;
			y = soda.map1.y_mesh;
			u = squeeze(soda.map1.u_d01(:,:,j))';
			v = squeeze(soda.map1.v_d01(:,:,j))';
			xq = soda.d01.lon_g(g,:)'-360;	% x_mesh is degW, therefore xq needs to be degW
			yq = soda.d01.lat_g(g,:)';
			uq = interp2(x,y,u,xq,yq);
			vq = interp2(x,y,v,xq,yq);
			uvq = sqrt(uq.^2 + vq.^2); 	% resultant velocity of u and v
			soda.d01.u_gg(ggcount,:) = uq';
			soda.d01.v_gg(ggcount,:) = vq';
			soda.d01.uv_gg(ggcount,:) = uvq';
			soda.d01.lon_gg(ggcount,:) = soda.d01.lon_g(g,:)';
			soda.d01.lat_gg(ggcount,:) = soda.d01.lat_g(g,:)';

			% v48
			soda.d01.u = uq;
			soda.d01.v = vq;
			soda.d01.uv = uvq;
			% calc change dx and dy between current substep and next substep 
			% dx is change east-west position and dy is north-south position, both in Km
			soda.d01.dx = soda.timestep .* soda.d01.u ./ 1000 / soda.nsubsteps;			
			soda.d01.dy = soda.timestep .* soda.d01.v ./ 1000 / soda.nsubsteps;		
			% calc change dlon and dlat between current substep increment and next substep (dlon is change east-west position and dlat is north-south position, both in degrees)
			soda.d01.dlon = km2deg(soda.d01.dx);		
			soda.d01.dlat = km2deg(soda.d01.dy);

			for k = 1:soda.ndrift*soda.drifter_particles		% loop through drifter particles

				% % v48
				% soda.d01.u = uq(k,1);
				% soda.d01.v = vq(k,1);
				% soda.d01.uv = uvq(k,1);

				% if u and v are both 0 that means the drifter has drifted over the model boundary
				% therefore reset the position to the previous substep to find u and v
				if i>1 & soda.d01.uv == 0

					gglast = find(soda.d01.uv_gg(:,k),1,'last');	% gg index of last non-zero velocity for this drifter
					if gglast == size(soda.d01.uv_gg,1)
						gglast = max(1,ggcount - 1);
					end
					ggback = gglast - soda.nsubsteps / 30;			% back up to find lon lat 1 day before the last non-zero velocity
					if ggback < 1
						ggback = gglast;
					end
					soda.d01.lon_g(g,k) = soda.d01.lon_gg(ggback,k);
					soda.d01.lat_g(g,k) = soda.d01.lat_gg(ggback,k);
					x = soda.map1.x_mesh;
					y = soda.map1.y_mesh;
					u = squeeze(soda.map1.u_d01(:,:,j))';
					v = squeeze(soda.map1.v_d01(:,:,j))';
					xq2 = soda.d01.lon_g(g,k)-360;	% x_mesh is degW, therefore xq needs to be degW
					yq2 = soda.d01.lat_g(g,k);
					uq2 = interp2(x,y,u,xq2,yq2);
					vq2 = interp2(x,y,v,xq2,yq2);

					% % v48
					% soda.d01.u = uq2;
					% soda.d01.v = vq2;
					soda.d01.u(k,1) = uq2;
					soda.d01.v(k,1) = vq2;
					% calc change dx and dy between current substep and next substep 
					% dx is change east-west position and dy is north-south position, both in Km
					soda.d01.dx(k,1) = soda.timestep .* soda.d01.u(k,1) ./ 1000 ./ soda.nsubsteps;			
					soda.d01.dy(k,1) = soda.timestep .* soda.d01.v(k,1) ./ 1000 ./ soda.nsubsteps;		
					% calc change dlon and dlat between current substep increment and next substep (dlon is change east-west position and dlat is north-south position, both in degrees)
					soda.d01.dlon(k,1) = km2deg(soda.d01.dx(k,1));		
					soda.d01.dlat(k,1) = km2deg(soda.d01.dy(k,1));

				end

				% % v48
				% % calc change dx and dy between current substep and next substep 
				% % dx is change east-west position and dy is north-south position, both in Km
				% soda.d01.dx = soda.timestep * soda.d01.u / 1000 / soda.nsubsteps;			
				% soda.d01.dy = soda.timestep * soda.d01.v / 1000 / soda.nsubsteps;		
				% % calc change dlon and dlat between current substep increment and next substep (dlon is change east-west position and dlat is north-south position, both in degrees)
				% soda.d01.dlon = km2deg(soda.d01.dx);		
				% soda.d01.dlat = km2deg(soda.d01.dy);
				% % find new postion of drifter at the end of this substep
				% soda.d01.lon_g(g+1,k) = soda.d01.lon_g(g,k) + soda.d01.dlon;
				% soda.d01.lat_g(g+1,k) = soda.d01.lat_g(g,k) + soda.d01.dlat;
				soda.d01.lon_g(g+1,k) = soda.d01.lon_g(g,k) + soda.d01.dlon(k,1);
				soda.d01.lat_g(g+1,k) = soda.d01.lat_g(g,k) + soda.d01.dlat(k,1);

				if g == soda.nsubsteps

					% if u and v are both 0 that means the drifter has drifted over the model boundary
					% therefore reset the final position at the end of this timestep 
					% to what it was at the start of this timestep (same as the end of the previous timestep)
					if i>2 & soda.d01.uv == 0

						% soda.d01.stuck(k,1) = true;
						disp(['Drifter stuck onshore at i=', num2str(i), ', ','j=', num2str(j), ', ','k=', num2str(k), ', ','g=', num2str(g)])

						gglast = find(soda.d01.uv_gg(:,k),1,'last');	% gg index of last non-zero velocity for this drifter
						if gglast == size(soda.d01.uv_gg,1)
							gglast = max(1,ggcount - 1);
						end
						ggback = gglast - soda.nsubsteps / 30;			% back up to find lon lat 1 day before the last non-zero velocity
						if ggback < 1
							ggback = gglast;
						end
						soda.d01.lon_g(g+1,k) = soda.d01.lon_gg(ggback,k);
						soda.d01.lat_g(g+1,k) = soda.d01.lat_gg(ggback,k);

						x = soda.map1.x_mesh;
						y = soda.map1.y_mesh;
						u = squeeze(soda.map1.u_d01(:,:,j))';
						v = squeeze(soda.map1.v_d01(:,:,j))';
						xq2 = soda.d01.lon_g(g+1,k)-360;	% x_mesh is degW, therefore xq needs to be degW
						yq2 = soda.d01.lat_g(g+1,k);
						uq2 = interp2(x,y,u,xq2,yq2);
						vq2 = interp2(x,y,v,xq2,yq2);

						% % v48
						% soda.d01.u = uq2;
						% soda.d01.v = vq2;
						soda.d01.u(k,1) = uq2;
						soda.d01.v(k,1) = vq2;

					end

					soda.d01.lon_kk_rsp(i+1,k) = soda.d01.lon_g(g+1,k);
					soda.d01.lat_kk_rsp(i+1,k) = soda.d01.lat_g(g+1,k);

				end	

				if gcount >= soda.nsubsteps / 4
					if soda.gifweekly
						if iweek > 1
							soda.d01.lon_weekly(iweek,k) = soda.d01.lon_g(g+1,k);
							soda.d01.lat_weekly(iweek,k) = soda.d01.lat_g(g+1,k);
						end
					end

					if k == soda.ndrift * soda.drifter_particles

						iweek = iweek + 1;
						gcount = 0;
					end
				end
			end		% for k
		end			% for g
	end				% for i

else 	% reverse drifter tracking backward in time from specified final position

	% reverse the sign of u and v to use the opposite direction of velocity vectors for reverse drifter tracking
	soda.map1.u_d01 = -soda.map1.u_d01;
	soda.map1.v_d01 = -soda.map1.v_d01;

	iweek=soda.picknmonths*4-4;
	gcount=0;

	iskip = (soda.pickstartyear - soda.startyear) * 12 + soda.pickstartmonth - 1;

	for i = soda.picknmonths:-1:2

	disp(['Calculating drifter position for month ', num2str(i), ' ...'])

		j = i + iskip;	 % index for pickyear to select u and v in the selected years

		soda.timestep = (soda.time(j) - soda.time(j-1)) * 86400;		% time step in seconds

		for g = 1:soda.nsubsteps	% divide the timestep into nsubsteps increments

			gcount = gcount + 1;
			ggcount = ggcount + 1;	% cumulative across all months and substeps

			if g == 1
				soda.d01.lon_g(g,:) = soda.d01.lon_kk_rsp(i,:);
				soda.d01.lat_g(g,:) = soda.d01.lat_kk_rsp(i,:);
			end

			% use interp2 to interpolate uq and vq for all current drifter positions xq and yq
			x = soda.map1.x_mesh;
			y = soda.map1.y_mesh;
			u = squeeze(soda.map1.u_d01(:,:,j))';
			v = squeeze(soda.map1.v_d01(:,:,j))';
			xq = soda.d01.lon_g(g,:)'-360;	% x_mesh is degW, therefore xq needs to be degW
			yq = soda.d01.lat_g(g,:)';
			uq = interp2(x,y,u,xq,yq);
			vq = interp2(x,y,v,xq,yq);
			uvq = sqrt(uq.^2 + vq.^2); 	% resultant velocity of u and v
			soda.d01.u_gg(ggcount,:) = uq';
			soda.d01.v_gg(ggcount,:) = vq';
			soda.d01.uv_gg(ggcount,:) = uvq';
			soda.d01.lon_gg(ggcount,:) = soda.d01.lon_g(g,:)';
			soda.d01.lat_gg(ggcount,:) = soda.d01.lat_g(g,:)';

			% v48
			soda.d01.u = uq;
			soda.d01.v = vq;
			soda.d01.uv = uvq;
			% calc change dx and dy between current substep and next substep 
			% dx is change east-west position and dy is north-south position, both in Km
			soda.d01.dx = soda.timestep .* soda.d01.u ./ 1000 / soda.nsubsteps;			
			soda.d01.dy = soda.timestep .* soda.d01.v ./ 1000 / soda.nsubsteps;		
			% calc change dlon and dlat between current substep increment and next substep (dlon is change east-west position and dlat is north-south position, both in degrees)
			soda.d01.dlon = km2deg(soda.d01.dx);		
			soda.d01.dlat = km2deg(soda.d01.dy);

			for k = 1:soda.ndrift*soda.drifter_particles		% loop through drifter particles

				% % v48
				% soda.d01.u = uq(k,1);
				% soda.d01.v = vq(k,1);
				% soda.d01.uv = uvq(k,1);

				% if u and v are both 0 that means the drifter has drifted over the model boundary
				% therefore reset the position to the previous substep to find u and v
				if i<soda.picknmonths & soda.d01.uv==0

					gglast = find(soda.d01.uv_gg(:,k),1,'last');	% gg index of last non-zero velocity for this drifter
					if gglast == size(soda.d01.uv_gg,1)
						gglast = max(1,ggcount - 1);
					end
					ggback = gglast - soda.nsubsteps / 30;			% back up to find lon lat 1 day before the last non-zero velocity
					if ggback < 1
						ggback = gglast;
					end
					soda.d01.lon_g(g,k) = soda.d01.lon_gg(ggback,k);
					soda.d01.lat_g(g,k) = soda.d01.lat_gg(ggback,k);
					x = soda.map1.x_mesh;
					y = soda.map1.y_mesh;
					u = squeeze(soda.map1.u_d01(:,:,j))';
					v = squeeze(soda.map1.v_d01(:,:,j))';
					xq2 = soda.d01.lon_g(g,k)-360;	% x_mesh is degW, therefore xq needs to be degW
					yq2 = soda.d01.lat_g(g,k);
					uq2 = interp2(x,y,u,xq2,yq2);
					vq2 = interp2(x,y,v,xq2,yq2);

					% % v48
					% soda.d01.u = uq2;
					% soda.d01.v = vq2;
					soda.d01.u(k,1) = uq2;
					soda.d01.v(k,1) = vq2;
					% calc change dx and dy between current substep and next substep 
					% dx is change east-west position and dy is north-south position, both in Km
					soda.d01.dx(k,1) = soda.timestep .* soda.d01.u(k,1) ./ 1000 ./ soda.nsubsteps;			
					soda.d01.dy(k,1) = soda.timestep .* soda.d01.v(k,1) ./ 1000 ./ soda.nsubsteps;		
					% calc change dlon and dlat between current substep increment and next substep (dlon is change east-west position and dlat is north-south position, both in degrees)
					soda.d01.dlon(k,1) = km2deg(soda.d01.dx(k,1));		
					soda.d01.dlat(k,1) = km2deg(soda.d01.dy(k,1));

				end

				% % v48
				% % calc change dx and dy between current substep and next substep 
				% % dx is change east-west position and dy is north-south position, both in Km
				% soda.d01.dx = soda.timestep * soda.d01.u / 1000 / soda.nsubsteps;			
				% soda.d01.dy = soda.timestep * soda.d01.v / 1000 / soda.nsubsteps;		
				% % calc change dlon and dlat between current substep increment and next substep (dlon is change east-west position and dlat is north-south position, both in degrees)
				% soda.d01.dlon = km2deg(soda.d01.dx);		
				% soda.d01.dlat = km2deg(soda.d01.dy);
				% % find new postion of drifter at the end of this substep
				% soda.d01.lon_g(g+1,k) = soda.d01.lon_g(g,k) + soda.d01.dlon;
				% soda.d01.lat_g(g+1,k) = soda.d01.lat_g(g,k) + soda.d01.dlat;
				soda.d01.lon_g(g+1,k) = soda.d01.lon_g(g,k) + soda.d01.dlon(k,1);
				soda.d01.lat_g(g+1,k) = soda.d01.lat_g(g,k) + soda.d01.dlat(k,1);

				if g == soda.nsubsteps
					% if u and v are both 0 that means the drifter has drifted over the model boundary
					% therefore reset the final position at the end of this timestep 
					% to what it was at the start of this timestep (same as the end of the previous timestep)
					if i<soda.picknmonths-1 & soda.d01.uv==0

						disp(['Drifter stuck onshore at i=', num2str(i), ', ','j=', num2str(j), ', ','k=', num2str(k), ', ','g=', num2str(g)])
						gglast = find(soda.d01.uv_gg(:,k),1,'last');	% gg index of last non-zero velocity for this drifter
						if gglast == size(soda.d01.uv_gg,1)
							gglast = max(1,ggcount - 1);
						end
						ggback = gglast - soda.nsubsteps / 30;			% back up to find lon lat 1 day before the last non-zero velocity
						if ggback < 1
							ggback = gglast;
						end
						soda.d01.lon_g(g+1,k) = soda.d01.lon_gg(ggback,k);
						soda.d01.lat_g(g+1,k) = soda.d01.lat_gg(ggback,k);
						%
						x = soda.map1.x_mesh;
						y = soda.map1.y_mesh;
						u = squeeze(soda.map1.u_d01(:,:,j))';
						v = squeeze(soda.map1.v_d01(:,:,j))';
						xq2 = soda.d01.lon_g(g+1,k)-360;	% x_mesh is degW, therefore xq needs to be degW
						yq2 = soda.d01.lat_g(g+1,k);
						uq2 = interp2(x,y,u,xq2,yq2);
						vq2 = interp2(x,y,v,xq2,yq2);

						% % v48
						% soda.d01.u = uq2;
						% soda.d01.v = vq2;
						soda.d01.u(k,1) = uq2;
						soda.d01.v(k,1) = vq2;

					end

					soda.d01.lon_kk_rsp(i-1,k) = soda.d01.lon_g(g+1,k);
					soda.d01.lat_kk_rsp(i-1,k) = soda.d01.lat_g(g+1,k);

				end	% if g== soda.nsubsteps

				if gcount >= soda.nsubsteps / 4
					if soda.gifweekly
						if iweek < soda.picknmonths*4
							soda.d01.lon_weekly(iweek,k) = soda.d01.lon_g(g+1,k);
							soda.d01.lat_weekly(iweek,k) = soda.d01.lat_g(g+1,k);
						end
					end

					if k == soda.ndrift * soda.drifter_particles
						gcount = 0;
						iweek = iweek - 1;
					end
				end

			end		% k loop
		end			% g loop
	end				% i loop

	% put the u and v velocity vectors back into the correct sign after calculationg the reverse drifter tracking
	soda.map1.u_d01 = -soda.map1.u_d01;
	soda.map1.v_d01 = -soda.map1.v_d01;

end

% - - - - -
% save csv files of 

% delete old csv with this name if it exists
if exist(soda.csvname_lon_monthly, 'file')==2
  delete(soda.csvname_lon_monthly);
end
if exist(soda.csvname_lat_monthly, 'file')==2
  delete(soda.csvname_lat_monthly);
end
if exist(soda.csvname_lon_weekly, 'file')==2
  delete(soda.csvname_lon_weekly);
end
if exist(soda.csvname_lat_weekly, 'file')==2
  delete(soda.csvname_lat_weekly);
end
%
% monthly
% monthly column names
clear colname_monthly
colname_monthly = strings();
colname_monthly(1)='year';
colname_monthly(2)='month';
colname_monthly(3)='time_idx';	% index of soda.time days since 1/1/1980
kkcount = 3;
% for i = 1:soda.ndrift
	% for j = 1:soda.drifter_particles
for i = 1:soda.drifter_particles
	for j = 1:soda.ndrift
		kkcount=kkcount+1;
		% colname_monthly(kkcount) = ['d_',num2str(i),'_',num2str(j)];
		colname_monthly(kkcount) = ['d_',num2str(j),'_',num2str(i)];
	end
end
% monthly values for year and month columns
clear col1 col2 col3
iskip = (soda.pickstartyear - soda.startyear) * 12 + soda.pickstartmonth - 1;
for i = 1:soda.picknmonths
	col1(i,1) = soda.year(i+iskip);
	col2(i,1) = soda.month(i+iskip);
	col3(i,1) = i+iskip;	% time index of soda.time days since 1/1/1980 corresponding to this selected month
end
% monthly lon
clear Mtx_lon T_lon_monthly
Mtx_lon = [col1, col2, col3, soda.d01.lon_kk_rsp];  
T_lon_monthly=array2table(Mtx_lon);
T_lon_monthly.Properties.VariableNames(1:numel(colname_monthly))=colname_monthly;
writetable(T_lon_monthly,soda.csvname_lon_monthly);
% monthly lat
clear Mtx_lat T_lat_monthly
Mtx_lat = [col1, col2, col3, soda.d01.lat_kk_rsp];  
T_lat_monthly=array2table(Mtx_lat);
T_lat_monthly.Properties.VariableNames(1:numel(colname_monthly))=colname_monthly;
writetable(T_lat_monthly,soda.csvname_lat_monthly);
%
% weekly
% weekly column names
clear colname_weekly
colname_weekly = strings();
colname_weekly(1)='year';
colname_weekly(2)='month';
colname_weekly(3)='week';
kkcount = 3;
% for i = 1:soda.ndrift
	% for j = 1:soda.drifter_particles
for i = 1:soda.drifter_particles
	for j = 1:soda.ndrift
		kkcount=kkcount+1;
		% colname_weekly(kkcount) = ['d_',num2str(i),'_',num2str(j)];
		colname_weekly(kkcount) = ['d_',num2str(j),'_',num2str(i)];
	end
end
% weekly values for year and month columns
clear col1 col2 col3
iskip = (soda.pickstartyear - soda.startyear) * 12 + soda.pickstartmonth - 1;
icount = 0;
for i = 1:soda.picknmonths
	for j = 1:4
		icount = icount+1;
		col1(icount,1) = soda.year(i+iskip);
		col2(icount,1) = soda.month(i+iskip);
		col3(icount,1) = j;
	end
end
% weekly lon
clear Mtx_lon T_lon_weekly
Mtx_lon = [col1, col2, col3, soda.d01.lon_weekly];  
T_lon_weekly=array2table(Mtx_lon);
T_lon_weekly.Properties.VariableNames(1:numel(colname_weekly))=colname_weekly;
writetable(T_lon_weekly,soda.csvname_lon_weekly);
% weekly lat
clear Mtx_lat T_lat_weekly
Mtx_lat = [col1, col2, col3, soda.d01.lat_weekly];  
T_lat_weekly=array2table(Mtx_lat);
T_lat_weekly.Properties.VariableNames(1:numel(colname_weekly))=colname_weekly;
writetable(T_lat_weekly,soda.csvname_lat_weekly);

% - - - - -
% reshape drifter particle arrays back from 2d to 3d for plotting

soda.d01.lon_kk = reshape(soda.d01.lon_kk_rsp,soda.picknmonths,soda.ndrift,[]);
soda.d01.lat_kk = reshape(soda.d01.lat_kk_rsp,soda.picknmonths,soda.ndrift,[]);

soda.d01.lon_weekly_kk = [];
soda.d01.lat_weekly_kk = [];
soda.d01.lon_weekly_kk = reshape(soda.d01.lon_weekly,soda.picknmonths*4,soda.ndrift,[]);
soda.d01.lat_weekly_kk = reshape(soda.d01.lat_weekly,soda.picknmonths*4,soda.ndrift,[]);

% - - - - -
% - - - - -

% ANIMATION

% - - - - -
% - - - - -

% - - - - -
% make animation frames and animated gif 
%
if ~soda.gifweekly

	% ----------
	% monthly frames
	
	iskip = (soda.pickstartyear - soda.startyear) * 12 + soda.pickstartmonth - 1;
	soda.pickyear = soda.pickstartyear;
	clear colormap
	ck=colormap(jet2(soda.ndrift));
	for i = 1:soda.picknmonths	% loop thru months
	disp(['Making animation of drifters for month ', num2str(i), ' ...'])
		figure(1)
		hFig=gcf;
		clf(hFig);
		% 
		hold on
		%
		% v51
		m_proj('mercator','lat',soda.maplat,'long',soda.maplon);
		%
		% overlay omara satdep or other data
		if ~strcmp(soda.overlay,'none')
			clear X Y Z
			X = soda.overlay_lon;
			Y = soda.overlay_lat;
			Z = squeeze(soda.overlay_data(:,:,i+iskip))';
			[C,h] = m_contourf(X,Y,Z,256);
			set(h,'LineColor','none')
			caxis([soda.overlay_min soda.overlay_max]);
			h=colorbar;
			h.Label.String=soda.overlay_name;
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
		%
		clear x y u v scale
		x = soda.map1.x_mesh;
		y = soda.map1.y_mesh;
		u = squeeze(soda.map1.u_d01(:,:,i+iskip))';
		v = squeeze(soda.map1.v_d01(:,:,i+iskip))';
		m_gshhs_i('patch',[.9 .9 .9],'edgecolor',[.9 .9 .9]);		% intermediate res coast
		% overlay current vector quiver arrows
		if soda.overlay_quiver
			% m_quiver(x,y,u,v,scale,'color','k','linewidth',.5);
			m_quiver(x,y,u,v,soda.quiverscale,'color',soda.quivercolor,'linewidth',soda.quiverwidth);
		end
		% overlay drifters

		% % % v40 - - - - -
		% for k = 1:soda.ndrift
			% for kk = 1:soda.drifter_particles
				% m_plot(soda.d01.lon_kk(i,k,kk)-360, soda.d01.lat_kk(i,k,kk),'o','linewidth',1,'color',soda.markercolor,'markerfacecolor',ck(k,:),'MarkerSize',3)
			% end
		% end

		% % v47 - - - - -
		if strcmp(soda.markerfacecolor,'jet')
			for k = 1:soda.ndrift
				% v50
				if strcmp(soda.markeredgecolor,'jet')
					m_plot(squeeze(soda.d01.lon_kk(i,k,:))-360, squeeze(soda.d01.lat_kk(i,k,:)),'o','linewidth',1,'color',ck(k,:),'markerfacecolor',ck(k,:),'MarkerSize',soda.markersize)
				else
					m_plot(squeeze(soda.d01.lon_kk(i,k,:))-360, squeeze(soda.d01.lat_kk(i,k,:)),'o','linewidth',1,'color',soda.markeredgecolor,'markerfacecolor',ck(k,:),'MarkerSize',soda.markersize)
				end
			end
		else
			for k = 1:soda.ndrift
				% v50
				if strcmp(soda.markeredgecolor,'jet')
					m_plot(squeeze(soda.d01.lon_kk(i,k,:))-360, squeeze(soda.d01.lat_kk(i,k,:)),'o','linewidth',1,'color',ck(k,:),'markerfacecolor',soda.markerfacecolor,'MarkerSize',soda.markersize)
				else
					m_plot(squeeze(soda.d01.lon_kk(i,k,:))-360, squeeze(soda.d01.lat_kk(i,k,:)),'o','linewidth',1,'color',soda.markeredgecolor,'markerfacecolor',soda.markerfacecolor,'MarkerSize',soda.markersize)
				end
			end
		end
		title(['Drifters at depth of_ ',soda.depthlabel,', ',num2str(soda.year(i+iskip,1)),'-',num2str(soda.month(i+iskip,1),'%02.f')])
		m_grid('box','off','linestyle','none','XaxisLocation','bottom','YaxisLocation','left','xtick',([-240 -230 -220 -210 -200 -190 -180 -170 -160 -150 -140 -130 -120 -110 -100]),'ytick',([0 10 20 30 40 50 60 65]),'FontSize',8);	% N Pacific
		%
		set(gca,'TickDir','out');	% the only other option is 'in'
		hold off
		%
		if soda.savepng
			% print uses figure 'PaperPosition' for the size
			if strcmp(soda.whichmap,'PacificN') 
				set(gcf, 'PaperPosition', [0 0 10 6])    % N Pacific
			elseif strcmp(soda.whichmap,'PacificNE') 
				set(gcf, 'PaperPosition', [0 0 7 7])    % NE Pacific
				elseif strcmp(soda.whichmap,'PacificNE2') 
					set(gcf, 'PaperPosition', [0 0 8 6])    % NE Pacific region of Dick's 135W line
				elseif strcmp(soda.whichmap,'PacificNW') 
					set(gcf, 'PaperPosition', [0 0 8 6])    % NW Pacific region of Japanese pteropod stations
			end
			% print(gcf, [pwd '/ani/Drifter_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),'_',num2str(soda.pickendyear),'_',num2str(i,'%02.f'),'.png'], '-dpng', '-r300' );   %save file as PNG w/ 300dpi
			print(gcf, [pwd '/ani/Drifter_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),'_',num2str(soda.pickendyear),'_',soda.overlay,'_mov',num2str(soda.movmean_months),'mo_',soda.ver,'_',num2str(i,'%02.f'),'.png'], '-dpng', '-r300' );   %save file as PNG w/ 300dpi
		end
		% write animated gif
		% set(0,'DefaultFigureColor',[1 1 1])
		set(gcf,'color','w');
			% gif/export_fig uses figure 'Position' for the size
			if strcmp(soda.whichmap,'PacificN') 
				set(gcf,'Position',[100 120 960 540]);	% 16:9 ratio for wide figs
			else 
				set(gcf,'Position',[100 120 560 420]);	% 4:3 ratio for normal figs
			end
		if i == 1
			gif(soda.gifname,'DelayTime',soda.gifdelay*soda.gifpause,'resolution',300,'overwrite',true,'nodither')
		elseif i == soda.picknmonths
			gif('DelayTime',soda.gifdelay*soda.gifpause)
		else
			gif('DelayTime',soda.gifdelay)
		end
	end 	% loop through months
	clear x y u v scale

else 		% weekly gif

	% ----------
	% weekly frames

	iweekcum=0;
	iskip = (soda.pickstartyear - soda.startyear) * 12 + soda.pickstartmonth - 1;
	soda.pickyear = soda.pickstartyear;
	clear colormap
	ck=colormap(jet2(soda.ndrift));
	for i = 1:soda.picknmonths	% loop thru months
		for iweek = 1:4
			iweekcum = iweekcum + 1;
			disp(['Making animation of drifters for week ', num2str(iweekcum), ' ...'])
			figure(1)
			hFig=gcf;
			clf(hFig);
			% 
			hold on
			%
			% v51
			m_proj('mercator','lat',soda.maplat,'long',soda.maplon);
			%
			% overlay omara satdep or other data
			if ~strcmp(soda.overlay,'none')
				clear X Y Z
				X = soda.overlay_lon;
				Y = soda.overlay_lat;
				Z = squeeze(soda.overlay_data(:,:,i+iskip))';
				[C,h] = m_contourf(X,Y,Z,256);
				set(h,'LineColor','none')
				caxis([soda.overlay_min soda.overlay_max]);
				h=colorbar;
				h.Label.String=soda.overlay_name;
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
			%
			clear x y u v scale
			x = soda.map1.x_mesh;
			y = soda.map1.y_mesh;
			u = squeeze(soda.map1.u_d01(:,:,i+iskip))';
			v = squeeze(soda.map1.v_d01(:,:,i+iskip))';
			m_gshhs_i('patch',[.9 .9 .9],'edgecolor',[.9 .9 .9]);		% intermediate res coast
			% overlay current vector quiver arrows
			if soda.overlay_quiver
				% m_quiver(x,y,u,v,scale,'color','k','linewidth',.5);
				m_quiver(x,y,u,v,soda.quiverscale,'color',soda.quivercolor,'linewidth',soda.quiverwidth);
			end
			% overlay drifters
			% % v47 - - - - -
			if strcmp(soda.markerfacecolor,'jet')
				for k = 1:soda.ndrift
					% v50
					if strcmp(soda.markeredgecolor,'jet')
						m_plot(squeeze(soda.d01.lon_weekly_kk(iweekcum,k,:))-360, squeeze(soda.d01.lat_weekly_kk(iweekcum,k,:)),'o','linewidth',1,'color',ck(k,:),'markerfacecolor',ck(k,:),'MarkerSize',soda.markersize)
					else
						m_plot(squeeze(soda.d01.lon_weekly_kk(iweekcum,k,:))-360, squeeze(soda.d01.lat_weekly_kk(iweekcum,k,:)),'o','linewidth',1,'color',soda.markeredgecolor,'markerfacecolor',ck(k,:),'MarkerSize',soda.markersize)
					end
				end
			else
				for k = 1:soda.ndrift
					% v50
					if strcmp(soda.markeredgecolor,'jet')
						m_plot(squeeze(soda.d01.lon_weekly_kk(iweekcum,k,:))-360, squeeze(soda.d01.lat_weekly_kk(iweekcum,k,:)),'o','linewidth',1,'color',ck(k,:),'markerfacecolor',soda.markerfacecolor,'MarkerSize',soda.markersize)
					else
						m_plot(squeeze(soda.d01.lon_weekly_kk(iweekcum,k,:))-360, squeeze(soda.d01.lat_weekly_kk(iweekcum,k,:)),'o','linewidth',1,'color',soda.markeredgecolor,'markerfacecolor',soda.markerfacecolor,'MarkerSize',soda.markersize)
					end
				end
			end
			title(['Drifters at depth of_ ',soda.depthlabel,', ',num2str(soda.year(i+iskip,1)),'-',num2str(soda.month(i+iskip,1),'%02.f')])
			m_grid('box','off','linestyle','none','XaxisLocation','bottom','YaxisLocation','left','xtick',([-240 -230 -220 -210 -200 -190 -180 -170 -160 -150 -140 -130 -120 -110 -100]),'ytick',([0 10 20 30 40 50 60 65]),'FontSize',8);	% N Pacific
			%
			set(gca,'TickDir','out');	% the only other option is 'in'
			hold off
			%
			% - - - - -
			% save animation frame
			if soda.savepng
				% print uses figure 'PaperPosition' for the size
				if strcmp(soda.whichmap,'PacificN') 
					set(gcf, 'PaperPosition', [0 0 10 6])    % N Pacific
				elseif strcmp(soda.whichmap,'PacificNE') 
					set(gcf, 'PaperPosition', [0 0 7 7])    % NE Pacific region of Nina's stations
				elseif strcmp(soda.whichmap,'PacificNE2') 
					set(gcf, 'PaperPosition', [0 0 8 6])    % NE Pacific region of Dick's 135W line
				elseif strcmp(soda.whichmap,'PacificNW') 
					set(gcf, 'PaperPosition', [0 0 8 6])    % NW Pacific region of Japanese pteropod stations
				end
				if iweek == 4	% pring monthly png files using last week of each month
					print(gcf, [pwd '/ani/Drifter_',soda.whichmap,'_',soda.whichdrifters,'_',num2str(soda.pickstartyear),'_',num2str(soda.pickendyear),'_',soda.overlay,'_mov',num2str(soda.movmean_months),'mo_',soda.ver,'_',num2str(i,'%02.f'),'.png'], '-dpng', '-r300' );   %save file as PNG w/ 300dpi
				end
			end
			% write animated gif
			% set(0,'DefaultFigureColor',[1 1 1])
			set(gcf,'color','w');
			% gif/export_fig uses figure 'Position' for the size
			if strcmp(soda.whichmap,'PacificN') 
				set(gcf,'Position',[100 120 960 540]);	% 16:9 ratio for wide figs
			else 
				set(gcf,'Position',[100 120 560 420]);	% 4:3 ratio for normal figs
			end
			if i == 1
				gif(soda.gifname,'DelayTime',soda.gifdelay*soda.gifpause,'resolution',300,'overwrite',true,'nodither')
			elseif i == soda.picknmonths*4
				gif('DelayTime',soda.gifdelay*soda.gifpause)
			else
				gif('DelayTime',soda.gifdelay)
			end

		end 	% loop through weeks

	end 	% loop through months

	clear x y u v scale

end 	% if soda.gifweekly


toc








