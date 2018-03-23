% Date: 2017-05-31 | Author: Y.Qin
% Update: 17-06-15 | added: CC, PBIAS, RMSE
% Update: 18-02-07 | added: output of 'E_avg'
% static_stefan.m:
%	use the simplified-Stefan solution to calculate the MTSFG
%   (adapted from 'post_valid.m')
clear
% ***
% Set workspace direction
root_dir = '..\';
matindir = [root_dir 'mat_input\'];
% Input Constant
IniYear  = 1961;
EndYear  = 2016;
yr_num   = EndYear-IniYear+1;
ObsIniY  = 1967;
ObsEndY  = 2015;
yr_obs   = ObsEndY-ObsIniY+1;
% Load the stations climatic forcing and land surface parameters
%	- stn_list : info about site row/col and station No.
%	- DDF  : Degree Days of Freezing of air temper.
load([matindir 'Site_data_input.mat'],'stn_list','DDF');
% Get the number of stations
[st_num, ~] = size(stn_list);
% Load the stations prior observation data
%	- ONLY mtsfg_obs : obs mtsfg at each site(1967-2015)
load([matindir 'Site_obs_mtsfg.mat'],'mtsfg_obs');
% ***
% Initialize of metrics matrix
%	- simulation of each year in each loop
Z_static = zeros(st_num, yr_num);
metrics  = zeros(st_num, 3);
E_avg    = zeros(st_num, 1);
% ***
% Loop of stations (SiteNo)
for st = 1:st_num
% GET station input (_st) DDF data in *.mat
	DDF_st  = DDF(st, :);
% GET station input (_st) MTSFG_obs data in *.mat
	Zobs_st = mtsfg_obs(st, :);
% LOOP OF simulating years (1961-2016)
	E_sum = 0; E_count = 0;
	for iobs = 1:yr_obs
		iddf = iobs +ObsIniY-IniYear;
% Check VALID DDF & MTSFG_obs (non-NaN)
		if DDF_st(iddf) > 0 && Zobs_st(iobs) > 0;
			E_count = E_count + 1;
			E_sum = E_sum + Zobs_st(iobs)/sqrt(DDF_st(iddf));
		end
	end
	E_avg(st) = E_sum / E_count;
	Z_static(st, :) = E_avg(st) * sqrt(DDF_st);
end
% Set zero-data as 'NaN'
Z_static(Z_static==0) = NaN;
% ***
% CALCULATE metrics(R2/BIAS/RMSE)
for st = 1:st_num
	line = 0;
	for tobs = 1:yr_obs
		tsim = tobs + ObsIniY-IniYear;
		if mtsfg_obs(st,tobs)>0 && ...
		isnan(Z_static(st,tsim))==0
			line = line+1;
			cal_obs(line) = mtsfg_obs(st,tobs);
			cal_sim(line) = Z_static(st,tsim);
		end
	end
	corr = corrcoef(cal_obs,cal_sim);
	bias = mean(cal_sim)/mean(cal_obs)-1;
	mse  = mean((cal_sim-cal_obs).^2);
	metrics(st,1) = corr(1,2);
	metrics(st,2) = bias;
	metrics(st,3) = mse^0.5;
	clear cal_obs;
	clear cal_sim;
end
% ***
% OUTPUT of xls:
xlswrite([matindir 'static_stefan.xlsx'], Z_static,1,'C2')
xlswrite([matindir 'static_stefan.xlsx'], metrics,1,'C17')
xlswrite([matindir 'static_stefan.xlsx'], E_avg,1,'M17')
