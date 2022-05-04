
% - - - - -
% soda_drifter.m
% Trajectories of passive drifters using monthly velocity u and v from SODA3.12.2 
% using Euler's method with the option to use either forward or reverse tracking,
% and optional overlay of biogeochemical tracers from PlankTOM12
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
% [pwd, '\nc\map1_2D_PlankTOM12_A_1980-2018_v20220429.nc']
% that was obtained by running the following scripts:
% 
% RECCAP_2021_3D_v43.m
% RECCAP_2021_2D_v45_write_netcdf.m
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
disp(['Reading input data ...'])

clear soda;

soda.ver = 'v63';					% version number for output file names
soda.forward_tracking = true;		% drifter tracking forward in time (true), or reverse drifter tracking for exposure history (false)
soda.movmean_months = 1;			% number of months for moving average u and v (1= no smoothing, 3= 3-month moving average, etc)
soda.pickstartyear = 2010;			% the year when the drifter experiment starts			
soda.pickstartmonth = 4;			% the month to start the drifter experiment in the start year (e.g. 1=Jan, 2=Feb, etc)
soda.pickendyear = 2014;  			% the year when the drifter experiment ends. for one year simulation input the end year with same value as the start year 
soda.pickendmonth = 4;				% the month to end the drifter experiment in the end year (e.g. 12=Dec, 11=Nov, etc)
soda.gifweekly = true;				% write weekly animation frames from g loop to gif 
soda.gifdelay = 0.25;				% delay between frames for animated gif (suggest 0.25 for 12-month, 0.125 for 5 year)
soda.gifpause = 4;					% scale factor for gifdelay to pause first and last frame (suggest 4 for 12-month, 8 for 5 year)

soda.whichmap = 'PacificNE4';		% Map extent to use for animated gif frames (add more options as needed in the code for making the animation) 
									% 'PacificN' = entire North Pacific from 0-65N and 128E-110W
									% 'PacificNE' = region of pteropod stations in the Northeast Pacific from 30N-61N and 165W-120W
									% 'PacificNE2' = region of 135W transect line in the Northeast Pacific from 10N-61N and 175E-110W
									% 'PacificNE3' = region of 153W transect line in the Northeast Pacific
									% 'PacificNE4' = region of 130W transect line in the Northeast Pacific
									% 'PacificNW' = region of Japanese pteropod stations from 125E-180E and 10N-60N
									% 'Osborne' = region of Osborne et al (1989) paper on trajectories for 12 months from start off SW coast of Japan
									% NOTE: any of the choices above may be used for the map extent, 
									% but if any other map extent is needed instead
									% it can be customized in the next if/elseif whichmap block of code below

soda.whichdrifters = 'nearshore';	% Initial position of drifters (add more choices as needed in next if/elseif block of code):  
									% 'Feely' = offshore pteropod stations 1-13 along the US West Coast cruise 2007, Fig 1 in Feely,R.A.,etal.,ChemicalandbiologicalimpactsofoceanacidiﬁcationalongthewestcoastofNorthAmerica, Estuarine, Coastal and Shelf Science (2016), http://dx.doi.org/10.1016/j.ecss.2016.08.043
									% 'nearshore' = same as 'Feely' with added nodes along Alaska and Mexico coasts
									% 'offshore' = Parallel to 'Feely' stations, but moved 5 deg east
									% 'pteropods' = L-shaped line of pteropod stations in the NE Pacific,
									% '153W' = north south at 153W from 58N-20N,
									% '135W' = north south at 135W from 55N-20N,
									% '130W' = north south at 130W from 52N-10N,
									% '33N' = east-west at 33N from 140E to 118W (latitude of the Kuroshio Extension)
									% '50N' = east-west at 50N from 155E to 128W 
									% 'Japan' = L-shaped line surrounding Japanese pteropod stations 25N from 137.5E to 160E and 160E from 25N to 50N
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
soda.drifter_spot_lon = 360-153;	% longitude of initial drifters (degE) (only used if soda.whichdrifters='Spot')
soda.drifter_spot_lat = 53;			% latitude of initial drifters (degN) (only used if soda.whichdrifters='Spot')

soda.markeredgecolor = 'jet';		% color for the outline of the drifter markers on the animation frames 
									% 	('jet'=jet colormap, otherwise any matlab color code such as 'k'=black, 'm'=magenta, 'c'=cyan, 'y'=yellow, 'w'=white, 'b'= blue, 'r'=red, etc) 
soda.markerfacecolor = 'k';			% color for the face of the drifter markers on the animation frames 
									% 	('jet'=jet colormap, otherwise any matlab color code such as 'k'=black, 'm'=magenta, 'c'=cyan, 'y'=yellow, 'w'=white, 'b'= blue, 'r'=red, etc) 
soda.markersize = 4;				% size of the marker for the drifter trajectory animation (points, e.g. matlab default size is 6 points)
soda.savepng = false;				% save png file for each animation frame

soda.overlay_anomaly = false;		% use the monthly anomaly compared with the long-term monthly averages for the selected map overlay
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
soda.quiverwidth = .1;			 	% line width of quiver arrows (suggested default 0.1)
soda.quiverscale = 5;			 	% scale factor for quiver length (suggested default 5)
soda.quivercolor = 'k';			 	% color for quiver arrows 

soda.depthlabel = '0-100m';			% label for plots to label the depth or range of depths that is used for the velocity vectors

% calculated values using the inputs above
soda.picknyears = soda.pickendyear-soda.pickstartyear+1;	% intermediate calc to find the number of months for duration of the experiment (picknmonths)
soda.picknmonths = (soda.picknyears-2) * 12 + (12-soda.pickstartmonth+1) + soda.pickendmonth;

soda.nsubdayspermonth = 28;			% intermediate increment to divide each month into approximately days
soda.nsubhoursperday = 8;			% intermediate increment to divide each day into approximateiy hours
soda.skipdays = 7;					% number of days to skip between animation frames (for example soda.skipdays=7 shows weekly animation frames)
soda.nsubsteps = soda.nsubdayspermonth * soda.nsubhoursperday;		% number of increments to divide each time step (nsubsteps=30*24=720 is approx 1 hour substep increment

soda.savemat = false;				% true or false to save the cropped u and v and intial drifter location data to a mat file
soda.integrate = true;				% true or false to calculate drifter trajectories 
soda.makegif = true;				% true or false to make animated gif

% END OF USER INPUTS
% ----------
% ----------
% ----------

% Edit the following as needed to change the map extent of the animation frames 
% select the map extent for the animation frames (add more options as needed)
if strcmp(soda.whichmap,'PacificN') 
	soda.maplat = [0 65];					% south and north latitudes of map extent (degN)
	% soda.maplon = [(120-360) -110];		% west and east longitudes of map extent (-degW)
	soda.maplon = [(115-360) -100];			% west and east longitudes of map extent (-degW)
elseif strcmp(soda.whichmap,'PacificNE') 
	soda.maplat = [30 61];					% south and north latitudes of map extent (degN)
	soda.maplon = [-175 -120];				% west and east longitudes of map extent (-degW)
elseif strcmp(soda.whichmap,'PacificNE2') 
	soda.maplat = [10 61];					% south and north latitudes of map extent (degN)
	soda.maplon = [-185 -110];				% west and east longitudes of map extent (-degW)
elseif strcmp(soda.whichmap,'PacificNE3') 
	soda.maplat = [10 61];					% south and north latitudes of map extent (degN)
	soda.maplon = [-195 -120];				% west and east longitudes of map extent (-degW)
elseif strcmp(soda.whichmap,'PacificNE4') 
	soda.maplat = [0 61];					% south and north latitudes of map extent (degN)
	% soda.maplon = [-190 -100];				% west and east longitudes of map extent (-degW)
	soda.maplon = [-190 -94.5];				% west and east longitudes of map extent (-degW)
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

% v60
soda.drifter_lon0 = [];
soda.drifter_lat0 = [];
if strcmp(soda.whichdrifters,'Spot')	% User-specified spot location
	soda.ndrift = 1;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	soda.drifter_lon0 = soda.drifter_spot_lon;					% degE
	soda.drifter_lat0 = soda.drifter_spot_lat;					% degN
elseif strcmp(soda.whichdrifters,'Feely') 	% offshore stations at 2007 West Coast cruise stations in Fig 1 by Feely,R.A.,etal.,ChemicalandbiologicalimpactsofoceanacidiﬁcationalongthewestcoastofNorthAmerica, Estuarine, Coastal and Shelf Science (2016), http://dx.doi.org/10.1016/j.ecss.2016.08.043
	soda.ndrift = 13;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	soda.drifter.lon0 = [227.135	230.517	232.456	233.069	232.771	233.566	234.892	236.135	237.412	238.953	240.478	243.412	245.218];		% degE
	soda.drifter.lat0 = [50.669	48.277	46.034	44.606	41.383	39.29	37.662	35.884	34.173	32.578	30.319	27.196	25.053];			% degN
elseif strcmp(soda.whichdrifters,'nearshore') 	% same as Feely with added nodes along coast of Alaska and Mexico
	soda.ndrift = 24;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	soda.drifter.lon0 = [196.55	202.78	208.49	212.66	217.42	220.33	222.14	224.15	227.135	230.517	232.456	233.069	232.771	233.566	234.892	236.135	237.412	238.953	240.478	243.412	245.218	248.64	253.56	257.55];		% degE
	soda.drifter.lat0 = [52.91	54.07	55.03	57.24	58.62	57.43	55.72	53.67	50.669	48.277	46.034	44.606	41.383	39.29	37.662	35.884	34.173	32.578	30.319	27.196	25.053	22.29	19.33	17.09];			% degN
elseif strcmp(soda.whichdrifters,'offshore') 	% offshore stations at 2007 West Coast cruise stations in Fig 1 by Feely,R.A.,etal.,ChemicalandbiologicalimpactsofoceanacidiﬁcationalongthewestcoastofNorthAmerica, Estuarine, Coastal and Shelf Science (2016), http://dx.doi.org/10.1016/j.ecss.2016.08.043
	soda.ndrift = 14;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	soda.drifter.lon0 = [220.6	222.135	225.517	227.456	228.069	227.771	228.566	229.892	231.135	232.412	233.953	235.478	238.412	240.218];		% degE
	soda.drifter.lat0 = [53.7	50.669	48.277	46.034	44.606	41.383	39.29	37.662	35.884	34.173	32.578	30.319	27.196	25.053];			% degN
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
elseif strcmp(soda.whichdrifters,'Kuroshio')	% Kuroshio Current transport transect off the SW coast of Japan
	soda.ndrift = 7;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	% initial position of drifters
	soda.drifter.lon0 = [141	141.375	141.75	142.125	142.5 142.875 143.25];		% degE
	soda.drifter.lat0 = [36	35.5	35	34.5	34 33.5 33];						% degN
elseif strcmp(soda.whichdrifters,'Oyashio') 	% Oyashio current transport transect at East Kamchatka
	soda.ndrift = 6;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	% initial position of drifters
	soda.drifter.lon0 = [160	160.625	161.25	161.875	162.5	163.125];		% degE
	soda.drifter.lat0 = [54	53.6875	53.375	53.0625	52.75	52.4375];						% degN
elseif strcmp(soda.whichdrifters,'KE')		% Kuroshio Extension transport transect 
	soda.ndrift = 43;					% number of drifters
	% initialize vectors for drifter lat and lon
	% v60
	% initial position of drifters
	soda.drifter.lon0(1:soda.ndrift) = 160;			% degE
	soda.drifter.lat0 = [49:-0.5:28];	% degN
end

% read netcdf file of soda lon, lat, time, u, v 
soda.map1_fn = [pwd '\nc\map1_velocity_2D_soda3_12_2_1980_2017_v20220421.nc'];	% edit as needed to correspond with the path to this nc file
soda.map1.lon = ncread(soda.map1_fn,'lon');
soda.map1.lat = ncread(soda.map1_fn,'lat');
soda.time = ncread(soda.map1_fn,'time');
soda.datenum = soda.time + datenum(datetime(1980,1,1,0,0,0));
soda.datetime = datetime(soda.datenum,'ConvertFrom','datenum');
soda.map1.u = ncread(soda.map1_fn,'u_KT100');	% KT100 is vertically integrated from 0-100m. KT200, KT500, and KT1000 are also available to select from this netcdf file
soda.map1.v = ncread(soda.map1_fn,'v_KT100');
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

% calculate movmean u and v moving average
if soda.movmean_months > 1
	for i = 1:numel(soda.map1.lon)
		for j = 1:numel(soda.map1.lat)
			soda.map1.u(i,j,:) = movmean(soda.map1.u(i,j,:),soda.movmean_months);
			soda.map1.v(i,j,:) = movmean(soda.map1.v(i,j,:),soda.movmean_months);
		end
	end
end

% map overlay data for animation frames
if ~strcmp(soda.overlay,'none')
	soda.overlay_fn = [pwd, '\nc\map1_2D_PlankTOM12_A_1980-2018_v20220429.nc'];		% edit as needed to correspond with the path to this nc file
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
	soda.overlay_cm = 'thermal';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
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
	soda.overlay_cm = 'haline';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
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
	soda.overlay_cm = 'speed';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
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
	soda.overlay_cm = 'speed';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
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
	soda.overlay_cm = 'speed';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
end
% use optional monthly anomaly of map overlay data
if soda.overlay_anomaly
	soda.overlay_data_rsp = reshape(soda.overlay_data, size(soda.overlay_data,1), size(soda.overlay_data,2), 12, []); 
	soda.overlay_data_monthly_avg = squeeze(nanmean2(soda.overlay_data_rsp, 4));		% monthly avg across all years
	soda.overlay_data_monthly_climatology = repmat(soda.overlay_data_monthly_avg, 1, 1, size(soda.overlay_data_rsp,4));		% monthly climatology, repeating time series of 12 monthly means for 39 years
	soda.overlay_data = soda.overlay_data - soda.overlay_data_monthly_climatology;		% convert the soda.overlay_data to the monthly anomaly
	soda.overlay_name = ['\Delta',soda.overlay_name];
	soda.overlay = [soda.overlay,'_anomaly'];
	soda.overlay_pivot = 0;
	soda.overlay_min = prctile(soda.overlay_data(:),1);
	soda.overlay_max = prctile(soda.overlay_data(:),99);
	% soda.overlay_min = nanmin2(soda.overlay_data(:));
	% soda.overlay_max = nanmax2(soda.overlay_data(:));
	soda.overlay_cm = 'balance';		% which ccmocean colormap per https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
	soda.overlay_opacity = 0.8;		% lighten or darken colormap (default opacity=1)
end
% mask non-water cells as NaN
if ~strcmp(soda.overlay,'none')
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
soda.overlay_time = ncread(soda.overlay_fn,'time');  
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
soda.extract.u(1:numel(soda.map1.lon),1:numel(soda.map1.lat),1:soda.picknmonths) = NaN;
soda.extract.v(1:numel(soda.map1.lon),1:numel(soda.map1.lat),1:soda.picknmonths) = NaN;
iskip = (soda.pickstartyear - soda.startyear) * 12 + soda.pickstartmonth - 1;
for i = 1:soda.picknmonths
	soda.extract.time(i,1) = soda.time(i+iskip,1);
	soda.extract.overlay_time(i,1) = soda.overlay_time(i+iskip,1);
	soda.extract.overlay_data(:,:,i) = soda.overlay_data(:,:,i+iskip);
	soda.extract.u(:,:,i) = soda.map1.u(:,:,i+iskip);
	soda.extract.v(:,:,i) = soda.map1.v(:,:,i+iskip);
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
soda.interp.u(1:numel(soda.map1.lon),1:numel(soda.map1.lat),1:numel(soda.interp.time)) = NaN;
soda.interp.v(1:numel(soda.map1.lon),1:numel(soda.map1.lat),1:numel(soda.interp.time)) = NaN;
%
disp(['Interpolating velocity vectors ...'])
for i = 1:numel(soda.map1.lon)
	for j = 1: numel(soda.map1.lat)
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
	[val,soda.crop_lon_idx] = min(abs(soda.map1.lon-soda.crop_lon_vec));					   
	[val,soda.crop_lat_idx] = min(abs(soda.map1.lat-soda.crop_lat_vec));	
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

	[soda.drifter.lon_rsp, soda.drifter.lat_rsp] = drifter_v02(soda.map1.lon, soda.map1.lat, soda.interp.time, ...
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

% FILTER WEEKLY FOR ANIMATION FRAMES

% ----------
% ----------

% ----------

if soda.integrate

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
	soda.drifter.weekly_u(1:numel(soda.map1.lon),1:numel(soda.map1.lat),1:nweeks) = NaN;
	soda.drifter.weekly_v(1:numel(soda.map1.lon),1:numel(soda.map1.lat),1:nweeks) = NaN;
	iweek = 0;
	icount = iskip - 1;
	for i = 1:numel(soda.interp.time)
		icount = icount + 1;

		%v61 include the last row of original data in the weekly filter
		% if icount == iskip
		if (icount == iskip) | (i == numel(soda.interp.time) & icount > 0)

			iweek = iweek + 1;
			icount = 0;
			soda.drifter.weekly_lon_rsp(iweek,:) = soda.drifter.lon_rsp(i,:);
			soda.drifter.weekly_lat_rsp(iweek,:) = soda.drifter.lat_rsp(i,:);
			soda.drifter.weekly_lon(iweek,:,:) = soda.drifter.lon(i,:,:);
			soda.drifter.weekly_lat(iweek,:,:) = soda.drifter.lat(i,:,:);
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

end

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

if soda.integrate & soda.makegif

	soda_drifter_animation_v01(soda.map1.lon-360, soda.map1.lat, soda.drifter.weekly_u, soda.drifter.weekly_v, ...
								soda.quiverscale, soda.quivercolor, soda.quiverwidth, soda.depthlabel, soda.overlay_quiver, ...
								soda.drifter.weekly_lon-360, soda.drifter.weekly_lat, ...
								soda.markerfacecolor, soda.markeredgecolor, soda.markersize, ...
								soda.overlay, soda.overlay_name, ...
								soda.overlay_lon, soda.overlay_lat, soda.drifter.weekly_overlay_data, ...
								soda.overlay_min, soda.overlay_max, soda.overlay_cm, soda.overlay_pivot, soda.overlay_opacity, ...
								soda.maplon, soda.maplat, soda.whichmap, ...
								soda.drifter.weekly_time, soda.gifname, soda.gifdelay, soda.gifpause)

end


toc


