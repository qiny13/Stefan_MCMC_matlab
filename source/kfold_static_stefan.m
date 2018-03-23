% Date: 2018-02-08 | Author: Y.Qin
% kfold_static_stefan.m:
%	k-fold cross validation to calculate the predictive error
%	using the simplified-Stefan solution to calculate the MTSFG
%   (adapted from 'static_stefan.m')
clear
% ***
% Set workspace direction
root_dir = '..\';
matindir = [root_dir 'mat_input\'];
% Input Constant
%	- SET: number of 'k' in k-fold cross validation
fold_n   = 5;
%	- SET: start/end year of calculation
IniYear  = 1961;
EndYear  = 2016;
yr_num   = EndYear-IniYear+1;
%	- SET: start/end year of observation
ObsIniY  = 1967;
ObsEndY  = 2015;
yr_obs   = ObsEndY-ObsIniY+1;
%	- SET: year gap of obs(1967-2015) and sim(1961-2016)
yrgp     = ObsIniY - IniYear;
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
%	- record RMSE at each station for each k-fold
rmse     = zeros(st_num, fold_n);
% ***
% Loop of stations (SiteNo)
for st = 1:st_num
% GET station input (_st) DDF data in *.mat
	DDF_st  = DDF(st, 1+yrgp:end-1);
% GET station input (_st) MTSFG_obs data in *.mat
	Zobs_st = mtsfg_obs(st, :);
% GET k-fold cross validation index
	indices = crossvalind('Kfold', yr_obs, fold_n);
% Loop of train folds (k-fold validation)
for ifold = 1:fold_n
	valid = (indices == ifold); % test subsamples
	train = ~valid; % subsamples for training
% select the training data (*_train)
	DDF_train  = DDF_st(:,train);
	Zobs_train = Zobs_st(:,train);
	k_train    = size(Zobs_train,2);
% select the testing data (*_valid)
	DDF_valid  = DDF_st(:,valid);
	Zobs_valid = Zobs_st(:,valid);
	k_valid    = size(Zobs_valid,2);	
	Zsim_valid = zeros(1,k_valid);
% Loop of training
	E_sum = 0; E_count = 0;
	for itra = 1:k_train
		%iddf = itra + yrgp;
% Check VALID DDF & MTSFG_obs (non-NaN)
		if DDF_train(itra) > 0 && Zobs_train(itra) > 0
			E_count = E_count + 1;
			E_sum = E_sum + Zobs_train(itra)/sqrt(DDF_train(itra));
		end
	end
	E_avg = E_sum / E_count;
% Loop of testing samples
	for ival = 1:k_valid
		if DDF_valid(ival) > 0 && Zobs_valid(ival) > 0
			Zsim_valid(ival) = E_avg * sqrt(DDF_valid(ival));
		else
			Zobs_valid(ival) = NaN;
			Zsim_valid(ival) = NaN;
		end
	end
	ind_nonan = ~isnan(Zobs_valid);
	rmse(st,ifold) = ...
		sqrt(mean((Zsim_valid(ind_nonan)-Zobs_valid(ind_nonan)).^2));
end % end of k-fold (ifold)
end % end of station (st)
% ***
% OUTPUT of xls:
xlswrite([matindir 'static_stefan.xlsx'], rmse,2,'C2')
