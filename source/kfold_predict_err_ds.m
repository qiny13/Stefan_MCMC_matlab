% Date: 2017-06-23 | Author: Y.Qin
% update: 17-09-04 | 3-term model ->> 1-term model
% update: 18-02-10 | API(alpha,p) ->> API(alpha,p, Ds)
% kfold_predict_err.m:
%	k-fold cross validation to calculate the predictive error
%   (adapted from: 'main_stefan_mcmc_solo.m')
% Function(s) Needed:
%   - StefanFunc.m
%   - LikelihoodFunc2.m
% Input Data Needed:
%   - .\mat_input\Site_data_input.mat
%    >> climatic forcing (DDF & AP) and land surface parameters
%   - .\mat_input\Site_obs_mtsfg.mat
%    >> prior observation data (Zobs & k)
clc
clear
% ***
% Set workspace direction
root_dir = '..\';
% Input Constant
%	- SET: number of 'k' in k-fold cross validation
fold_n     = 5;
%	- SET: number of maximum number of observation
k_st       = 2015 - 1967 +1;
%	- SET: year gap of obs(1967-2015) and sim(1961-2016)
yrgp       = 1967 - 1961;
%	- SET: number of loops in MCMC
LoopNum    = 10 *10^4;
% 	- SET: lhos: last half of samples
lhos       = LoopNum/2+1 : LoopNum;
%	- SET: standard deviation '/sigma' * Ds (mark-180210)
sig        = 0.05;
%	- SET: stno of calc stations
SiteNo     = [52943;52957;52968;56033;56043;56046;56065;56067;56074;56079;56173];
[st_num,~] = size(SiteNo);
% ***
% Load the input *.mat files (from preprocessor.m)
matindir = [root_dir 'mat_input\'];
out_dir  = [root_dir 'results\'];
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
%	- mtsfg_obs : obs mtsfg at each site
load([matindir 'Site_obs_mtsfg.mat']);
% Load the initial parameters
ParaIni = load('ParaIni.txt');   % add Ds, mark-180210
% ***
% Load the initial parameters
%	- initial values of parameters
ini_alp    = ParaIni(1,1);
ini_p      = ParaIni(4,1);
ini_Ds     = ParaIni(7,1);
%	- uniform random values of parameters
rndmat_a   = unifrnd(ParaIni(1,2),ParaIni(1,3), st_num,LoopNum);
rndmat_p   = unifrnd(ParaIni(4,2),ParaIni(4,3), st_num,LoopNum);
rndmat_d   = unifrnd(ParaIni(7,2),ParaIni(7,3), st_num,LoopNum);
%	- record values of parameters in each loop
Para_a     = zeros(LoopNum, st_num);
Para_p     = zeros(LoopNum, st_num);
Para_d     = zeros(LoopNum, st_num);
%	- record RMSE at each station for each k-fold
rmse       = zeros(st_num, fold_n);
% ***
tic
% Loop of stations (SiteNo)
for st = 1:st_num
% GET station row(strow) in *.mat
	[strow,~] = find(stn_list==SiteNo(st));
% GET station input (_st) forcing data or parameters in *.mat
	AP_st   = AP((strow-1)*AP_mnum+1:strow*AP_mnum, 1+yrgp:end-1);
	DDF_st  = DDF(strow, 1+yrgp:end-1);
	nf_st   = nf(strow);
	rhod_st = rhod(strow);
	sand_st = sand(strow);
	ths_st  = ths(strow);
% GET station input (_st) obs data in *.mat
	Zobs_st = mtsfg_obs(strow, :);
% GET k-fold cross validation index
	indices = crossvalind('Kfold', k_st, fold_n);
% ***
% Loop of train folds (k-fold validation)
for ifold = 1:fold_n
	valid = (indices == ifold); % test subsamples
	train = ~valid; % subsamples for training
% select the training data (*_train)
	AP_train   = AP_st(:,train);
	DDF_train  = DDF_st(:,train);
	Zobs_train = Zobs_st(:,train);
	k_train    = size(Zobs_train,2);
% select the testing data (*_valid)
	AP_valid   = AP_st(:,valid);
	DDF_valid  = DDF_st(:,valid);
	Zobs_valid = Zobs_st(:,valid);
	k_valid    = size(Zobs_valid,2);	
	Zsim_valid = zeros(1,k_valid);
% Loop of MCMC-loops
	for lp = 1:LoopNum
		if lp == 1
% Initial: GET "current" data (_crt) from ini_* (ParaIni)
			alp_crt = ini_alp;
			p_crt   = ini_p;
			Ds_crt  = ini_Ds;
		end
% Solo (1-term) API model
		bet_crt=0; q_crt=0; gam_crt=0; w_crt=0;
% Updating-Step 1: alpha (alp_)
		alp_rnd = rndmat_a(st,lp);
		r_A = LikelihoodFunc2(k_train,sig*Ds_crt,Zobs_train,...
			alp_rnd,bet_crt,gam_crt,p_crt,q_crt,w_crt,...
			AP_train,DDF_train,Ds_crt,nf_st,rhod_st,sand_st,ths_st);
		r_B = LikelihoodFunc2(k_train,sig*Ds_crt,Zobs_train,...
			alp_crt,bet_crt,gam_crt,p_crt,q_crt,w_crt,...
			AP_train,DDF_train,Ds_crt,nf_st,rhod_st,sand_st,ths_st);
		ratio = min(r_A / r_B, 1);
		U = unifrnd(0,1);
		if ratio > U
			alp_crt = alp_rnd;
		end
% Updating-Step 2: p
		p_rnd = rndmat_p(st,lp);
		r_A = LikelihoodFunc2(k_train,sig*Ds_crt,Zobs_train,...
			alp_crt,bet_crt,gam_crt,p_rnd,q_crt,w_crt,...
			AP_train,DDF_train,Ds_crt,nf_st,rhod_st,sand_st,ths_st);
		r_B = LikelihoodFunc2(k_train,sig*Ds_crt,Zobs_train,...
			alp_crt,bet_crt,gam_crt,p_crt,q_crt,w_crt,...
			AP_train,DDF_train,Ds_crt,nf_st,rhod_st,sand_st,ths_st);
		ratio = min(r_A / r_B, 1);
		U = unifrnd(0,1);
		if ratio > U
			p_crt = p_rnd;
		end
% Updating-Step 3: Ds
		Ds_rnd = rndmat_d(st,lp); % mark-180210
		r_A = LikelihoodFunc2(k_train,sig*Ds_rnd,Zobs_train,...
			alp_crt,bet_crt,gam_crt,p_crt,q_crt,w_crt,...
			AP_train,DDF_train,Ds_rnd,nf_st,rhod_st,sand_st,ths_st);
		r_B = LikelihoodFunc2(k_train,sig*Ds_crt,Zobs_train,...
			alp_crt,bet_crt,gam_crt,p_crt,q_crt,w_crt,...
			AP_train,DDF_train,Ds_crt,nf_st,rhod_st,sand_st,ths_st);
		ratio = min(r_A / r_B, 1);
		U = unifrnd(0,1);
		if ratio > U
			Ds_crt = Ds_rnd;
		end
% ***
% SAVE the current para-set of the loop 'lp'
		Para_a(lp,st) = alp_crt;
		Para_p(lp,st) = p_crt;
		Para_d(lp,st) = Ds_crt;
% Process display
		if mod(lp, LoopNum/5) == 0
			disp(['Finished: ' num2str(lp) ' of ' num2str(LoopNum) ' Loops'...
			' | ' num2str(ifold) ' of ' num2str(fold_n) ' Folds'...
			' | st-' num2str(st) ' of ' num2str(st_num) ' Sites'])
		end
	end % end of lp (i-loop)
% ***
% VALIDATE of k-fold cross validation
	m_a = median(Para_a(lhos,st));
	m_p = median(Para_p(lhos,st));
	m_d = median(Para_d(lhos,st));
% Solo (1-term) API model
	m_b=0; m_g=0; m_q=0; m_w=0;
% Loop of testing samples
	for iv = 1:k_valid
		if Zobs_valid(iv) > 0 && DDF_valid(iv) > 0
			Zsim_valid(iv) = StefanFunc( ...
				m_a, m_b, m_g, m_p, m_q, m_w, ...
				AP_valid(:,iv),DDF_valid(iv), ...
				m_d,nf_st,rhod_st,sand_st,ths_st);
		else
			Zobs_valid(iv) = NaN;
			Zsim_valid(iv) = NaN;
		end
	end
	ind_nonan = ~isnan(Zobs_valid);
	rmse(st,ifold) = ...
		sqrt(mean((Zsim_valid(ind_nonan)-Zobs_valid(ind_nonan)).^2));
end % end of k-fold (ifold)
disp(['---- Finish of st-' num2str(st) ' ----'])
toc
end % end of station (st)
% ***
save([out_dir 'kfold_rmse.mat'],'rmse');
