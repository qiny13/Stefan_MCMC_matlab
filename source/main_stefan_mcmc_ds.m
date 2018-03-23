% Date: 2017-03-09 | Author: Y.Qin
% update: 17-03-27 | ~,alpha,beta ->> alpha,beta,gamma
% update: 17-04-13 | one-thread ->> multi-thread (parfor)
% update: 17-04-19 | one-chain ->> multi-chain (PSRF cal)
% update: 17-09-03 | 3-term model ->> 1-term model
% update: 18-02-10 | API(alpha,p) ->> API(alpha,p, Ds)
% main_stefan_mcmc_ds.m:
%	The main process of MCMC-based Stefan Solution with 1-term API and Ds
%   (adapted from main_stefan_mcmc.m)
% Function(s) Needed:
%   - StefanFunc.m
%   - LikelihoodFunc.m
% Input Data Needed:
%   - .\mat_input\Site_data_input.mat
%    >> climatic forcing (DDF & AP) and land surface parameters
%   - .\mat_input\Site_obs_mtsfg.mat
%    >> prior observation data (Zobs & k)
clc
clear
% MULTI-THREAD CALCULATION (open # workers)
matlabpool local 6;
% ***
% Set workspace direction
root_dir = '..\';
% Input Constant
%	- SET: number of chains: '1' for solo cal; others for PSRF cal.
mulchain   = 5;
%	- SET: random seeds for initial condition
seed       = [3;7;13;17;23;31;37];
%	- SET: number of loops in MCMC
%   - runtime: 68 sec per 10,000 loop per station (solo-thread)
LoopNum    = 10 *10^4;
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
%	- num_obs   : number of valid data in each site
load([matindir 'Site_obs_mtsfg.mat']);
% Load the initial parameters
ParaIni = load('ParaIni.txt');   % add Ds, mark-180210
% ***
% Start of multi-chain
for ichain = 1:mulchain
% Load the initial parameters
%	- initial values of parameters
if ichain == 1 % the first chain (ini = mid_boundary)
	ini_alp    = ParaIni(1,1);
	ini_p      = ParaIni(4,1);
	ini_Ds     = ParaIni(7,1);
else % multiple chains
	rng(seed(1));
	inirnd_a   = unifrnd(ParaIni(1,2),ParaIni(1,3), 1,mulchain);
	ini_alp    = inirnd_a(ichain);
	rng(seed(4));
	inirnd_p   = unifrnd(ParaIni(4,2),ParaIni(4,3), 1,mulchain);
	ini_p      = inirnd_p(ichain);
	rng(seed(7));
	inirnd_d   = unifrnd(ParaIni(7,2),ParaIni(7,3), 1,mulchain);
	ini_Ds     = inirnd_d(ichain);
end
%	- uniform random values of parameters
rndmat_a   = unifrnd(ParaIni(1,2),ParaIni(1,3), st_num,LoopNum);
rndmat_p   = unifrnd(ParaIni(4,2),ParaIni(4,3), st_num,LoopNum);
rndmat_d   = unifrnd(ParaIni(7,2),ParaIni(7,3), st_num,LoopNum);
[pa_num,~] = size(ParaIni);
%	- record values of parameters in each loop
Para       = zeros(LoopNum, pa_num*st_num);
Para_a     = zeros(LoopNum, st_num);
Para_p     = zeros(LoopNum, st_num);
Para_d     = zeros(LoopNum, st_num);
% ***
tic  % display-timer start
% Parallel-Loop of stations (SiteNo)
parfor st = 1:st_num
% GET station row(strow) in *.mat
	[strow,~] = find(stn_list==SiteNo(st));
% GET station input (_st) forcing data or parameters in *.mat
	AP_st   = AP((strow-1)*AP_mnum+1 : strow*AP_mnum, :);
	DDF_st  = DDF(strow, :);
	nf_st   = nf(strow);
	rhod_st = rhod(strow);
	sand_st = sand(strow);
	ths_st  = ths(strow);
% GET station input (_st) obs data in *.mat
	Zobs_st = mtsfg_obs(strow, :);
	k_st    = num_obs(strow);
% ***
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
		r_A = LikelihoodFunc(k_st,sig*Ds_crt,Zobs_st,...
			alp_rnd,bet_crt,gam_crt,p_crt,q_crt,w_crt,...
			AP_st,DDF_st,Ds_crt,nf_st,rhod_st,sand_st,ths_st);
		r_B = LikelihoodFunc(k_st,sig*Ds_crt,Zobs_st,...
			alp_crt,bet_crt,gam_crt,p_crt,q_crt,w_crt,...
			AP_st,DDF_st,Ds_crt,nf_st,rhod_st,sand_st,ths_st);
		ratio = min(r_A / r_B, 1);
		U = unifrnd(0,1);
		if ratio > U
			alp_crt = alp_rnd;
		end
% Updating-Step 2: p
		p_rnd = rndmat_p(st,lp);
		r_A = LikelihoodFunc(k_st,sig*Ds_crt,Zobs_st,...
			alp_crt,bet_crt,gam_crt,p_rnd,q_crt,w_crt,...
			AP_st,DDF_st,Ds_crt,nf_st,rhod_st,sand_st,ths_st);
		r_B = LikelihoodFunc(k_st,sig*Ds_crt,Zobs_st,...
			alp_crt,bet_crt,gam_crt,p_crt,q_crt,w_crt,...
			AP_st,DDF_st,Ds_crt,nf_st,rhod_st,sand_st,ths_st);
		ratio = min(r_A / r_B, 1);
		U = unifrnd(0,1);
		if ratio > U
			p_crt = p_rnd;
		end
% Updating-Step 3: Ds
		Ds_rnd = rndmat_d(st,lp);  % mark-180210
		r_A = LikelihoodFunc(k_st,sig*Ds_rnd,Zobs_st,...
			alp_crt,bet_crt,gam_crt,p_crt,q_crt,w_crt,...
			AP_st,DDF_st,Ds_rnd,nf_st,rhod_st,sand_st,ths_st);
		r_B = LikelihoodFunc(k_st,sig*Ds_crt,Zobs_st,...
			alp_crt,bet_crt,gam_crt,p_crt,q_crt,w_crt,...
			AP_st,DDF_st,Ds_crt,nf_st,rhod_st,sand_st,ths_st);
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
			disp(['Finished: ' num2str(lp) ' of ' num2str(LoopNum) ...
			' Loops | st-' num2str(st) ' of ' num2str(st_num) ' Sites'])
		end
	end
	disp(['----- Finish loops of st-' num2str(st)])
end % end of parfor
toc
% save the output sliced variable
for st = 1:st_num
	Para(:,pa_num*(st-1)+1) = Para_a(:,st);
	Para(:,pa_num*(st-1)+4) = Para_p(:,st);
	Para(:,pa_num*(st-1)+7) = Para_d(:,st);
end
% ***
if mulchain == 1 % solo chain to save
	savename =  'main_stefan_mcmc';
else % multiple chains to save
	savename = ['main_stefan_mcmc_c' num2str(ichain)];
end
save([out_dir savename '.mat'],'Para','SiteNo','ParaIni');
disp(['File Saved: ' savename '.mat'])
end % end of ichain
% ***
matlabpool close; % close multi-thread
