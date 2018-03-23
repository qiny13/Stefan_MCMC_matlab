% Date: 2017-03-31 | Author: Y.Qin
% Update: 17-04-05
% Update: 17-04-17 | solo-result ->> results and metrics-pdf
% updata: 17-05-26 | solo chain ->> ensemble para dataset
% update: 18-02-11 | pa_num = 6->7, API(alpha,p) ->> API(alpha,p,Ds)
% post_valid.m:
%	Post-process load main_stefan_mcmc.mat and validate the
%   "optimized" parameters from MCMC results
clear
% ***
% Set workspace direction
root_dir = '..\';
matindir = [root_dir 'mat_input\'];
out_dir  = [root_dir 'results\'];
% Input Constant
IniYear  = 1961;
EndYear  = 2016;
yr_num   = EndYear-IniYear+1;
ObsIniY  = 1967;
ObsEndY  = 2015;
yr_obs   = ObsEndY-ObsIniY+1;
%	- SET: significance level of quantile(0.05->>95%)
qt       = 0.05;
% Load the stations climatic forcing and land surface parameters
%	- stn_list : info about site row/col and station No.
%   - AP_mnum  : AP number of months
%	- AP   : the antecedent precip index
%	- DDF  : Degree Days of Freezing of air temper.
% 	- nf   : n-factor value based on CMA site obs (Wang,YH. TgTa 201608)
% 	- rhod : soil bulk density (kg/m3) (HWSD, China Soil Map)
% 	- sand : sand content (%) (HWSD, China Soil Map)
% 	- ths  : saturated soil water constant (cm3/cm3) (Dai,YJ.,2013)
load([matindir 'Site_data_input.mat']);
% Load the stations prior observation data
%	- ONLY mtsfg_obs : obs mtsfg at each site(1967-2015)
load([matindir 'Site_obs_mtsfg.mat'],'mtsfg_obs');
% Load the ensemble para dataset (from post_para_plot.m)
%	- Ensemb : row = LoopNum/2*Chain_num (lhos) | col = (alp,bet,gam,p,q,w,Ds) x stations
%	- SiteNo : row = st_num
load([out_dir 'post_para.mat']);
[LoopNum, pa_col] = size(Ensemb);
[st_num, ~] = size(SiteNo);
pa_num = pa_col/st_num; % number of total parameters
% Initialize of metrics matrix
%	- record values of parameters in each loop
metrics = zeros(st_num, 9);
metloop = zeros(LoopNum, 3);
%	- simulation of each year in each loop
Zsim_st = zeros(LoopNum, yr_num);
%	- quantiles of simulation (low median high)
Zsim_lo = zeros(st_num, yr_num);
Zsim_me = zeros(st_num, yr_num);
Zsim_hi = zeros(st_num, yr_num);
% Loop of stations (SiteNo)
tic
for st = 1:st_num
% GET station row(strow) in *.mat
	[strow,~] = find(stn_list==SiteNo(st));
% GET station input (_st) forcing data or parameters in *.mat
	AP_st   = AP((strow-1)*AP_mnum+1 : strow*AP_mnum, :);
	DDF_st  = DDF(strow, :);
	nf_st   = nf(strow);	
	rhod_st = rhod(strow);
	sand_st = sand(strow);
	ths_st  = ths(strow);
% GET optimized (*_opt) parameter-sets from MCMC results
	alp_opt = Ensemb(:, pa_num*(st-1)+1);
	bet_opt = Ensemb(:, pa_num*(st-1)+2);
	gam_opt = Ensemb(:, pa_num*(st-1)+3);
	p_opt   = Ensemb(:, pa_num*(st-1)+4);
	q_opt   = Ensemb(:, pa_num*(st-1)+5);
	w_opt   = Ensemb(:, pa_num*(st-1)+6);
	Ds_opt  = Ensemb(:, pa_num*(st-1)+7);
% GET station input (_st) obs data in *.mat
	Zobs_st = mtsfg_obs(strow, :);
	Zobs_st(Zobs_st==0) = NaN;  % get rid of blank value
% LOOP OF mcmc chains (last half of samples)
	for ilp = 1 : LoopNum
% LOOP OF simulating years (1961-2016)
		for iy = 1:yr_num
			if DDF_st(iy) > 0
				Zsim_st(ilp,iy) = StefanFunc( ...
					alp_opt(ilp),bet_opt(ilp),gam_opt(ilp), ...
					p_opt(ilp),q_opt(ilp),w_opt(ilp), ...
					AP_st(:,iy),DDF_st(iy), ...
					Ds_opt(ilp),nf_st,rhod_st,sand_st,ths_st);
			else % DDF_st(iy) == 0 (nodata)
				Zsim_st(ilp,iy) = NaN;
			end
		end
% CALCULATE metrics(R2/BIAS/RMSE)
		line = 0;
		for tobs = 1:yr_obs
			tsim = tobs + ObsIniY-IniYear;
			if isnan(Zobs_st(tobs))==0 && ...
			isnan(Zsim_st(ilp,tsim))==0
				line = line+1;
				cal_obs(line) = Zobs_st(tobs);
				cal_sim(line) = Zsim_st(ilp,tsim);
			end
		end
		corr = corrcoef(cal_obs,cal_sim);
		bias = mean(cal_sim)/mean(cal_obs)-1;
		mse  = mean((cal_sim-cal_obs).^2);
		metloop(ilp,1) = corr(1,2);
		metloop(ilp,2) = bias;
		metloop(ilp,3) = mse^0.5;
		line = 0;
		clear cal_obs;
		clear cal_sim;
% Process display
		if mod(ilp, LoopNum/10) == 0
			disp(['Validated: ' num2str(ilp) ' of ' num2str(LoopNum) ...
			' Loops | st-' num2str(st) ' of ' num2str(st_num) ' Sites'])
		end
	end % end of i-loop (ilp)
	toc
% Quantile of metrics(R2/BIAS/RMSE)
	quant = [qt/2, 0.50, 1-qt/2];
	metrics(st,1:3) = quantile(metloop(:,1),quant);
	metrics(st,4:6) = quantile(metloop(:,2),quant);
	metrics(st,7:9) = quantile(metloop(:,3),quant);
% Quantile of mcmc-simul (95%-range median)
	for iy = 1:yr_num
		quansim = quantile(Zsim_st(:,iy),quant);
		Zsim_lo(st,iy) = quansim(1);
		Zsim_me(st,iy) = quansim(2);
		Zsim_hi(st,iy) = quansim(3);
	end
end
% ***
% OUTPUT of xls: save mid and std of each parameter
xlswrite([out_dir 'validmat.xlsx'],metrics,1,'B2')
% OUTPUT of xls: save mode, mid, and std of each parameter
xlswrite([out_dir 'post_valid.xlsx'],Zsim_lo,1, 'C2')
xlswrite([out_dir 'post_valid.xlsx'],Zsim_me,1,'C14')
xlswrite([out_dir 'post_valid.xlsx'],Zsim_hi,1,'C26')
% OUTPUT of *.mat: for post_valid_plot.m
save([out_dir 'post_valid.mat'],'Zsim_lo','Zsim_me','Zsim_hi',...
	'ObsIniY','ObsEndY','IniYear','EndYear','SiteNo');
