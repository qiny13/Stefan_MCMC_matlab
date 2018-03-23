% Date: 2017-06-29 | Author: Y.Qin
% update: 18-02-22 | API(alpha,p) ->> API(alpha,p, Ds)
%                    num_prior = 5 ->> 7
% sens_prior_selct_ds.m:
%	Sensitive analysis of the priors of parameters
%   (adapted from: 'main_stefan_mcmc_ds.m')
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
%	- SET: number of priors
pr_num     = 7;
%	- SET: number of loops in MCMC
LoopNum    = 4 *10^4;
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
%	- num_obs   : number of valid data in each site
load([matindir 'Site_obs_mtsfg.mat']);
% Load the initial parameters
ParaIni_0 = load('ParaIni.txt');
[pa_num,~] = size(ParaIni_0);
% Set the output matrix
Para = zeros(LoopNum, pa_num*st_num*pr_num);
% ***
% Prior selection: para-set 1 to pr_num
for ipr = 1:pr_num
switch ipr
	case 1
		f_rec = 1.0; f_exp = 1.0; f_ds = 1.0;
	case 2
		f_rec = 0.5; f_exp = 1.0; f_ds = 1.0;
	case 3
		f_rec = 2.0; f_exp = 1.0; f_ds = 1.0;
	case 4
		f_rec = 1.0; f_exp = 0.5; f_ds = 1.0;
	case 5
		f_rec = 1.0; f_exp = 2.0; f_ds = 1.0;
	case 6
		f_rec = 1.0; f_exp = 1.0; f_ds = 0.5;
	case 7
		f_rec = 1.0; f_exp = 1.0; f_ds = 2.0;
end
% Load the initial parameters
%   - set para-set senarios
ParaIni(1:3,:) = f_rec * ParaIni_0(1:3,:);
ParaIni(4:6,:) = f_exp * ParaIni_0(4:6,:);
ParaIni(7,  :) = f_ds  * ParaIni_0(7,  :);
%	- initial values of parameters
ini_alp    = ParaIni(1,1);
ini_bet    = ParaIni(2,1);
ini_gam    = ParaIni(3,1);
ini_p      = ParaIni(4,1);
ini_q      = ParaIni(5,1);
ini_w      = ParaIni(6,1);
ini_Ds     = ParaIni(7,1);
%	- uniform random values of parameters
rndmat_a   = unifrnd(ParaIni(1,2),ParaIni(1,3), st_num,LoopNum);
rndmat_b   = unifrnd(ParaIni(2,2),ParaIni(2,3), st_num,LoopNum);
rndmat_g   = unifrnd(ParaIni(3,2),ParaIni(3,3), st_num,LoopNum);
rndmat_p   = unifrnd(ParaIni(4,2),ParaIni(4,3), st_num,LoopNum);
rndmat_q   = unifrnd(ParaIni(5,2),ParaIni(5,3), st_num,LoopNum);
rndmat_w   = unifrnd(ParaIni(6,2),ParaIni(6,3), st_num,LoopNum);
rndmat_d   = unifrnd(ParaIni(7,2),ParaIni(7,3), st_num,LoopNum);
%	- record values of parameters in each loop
Para_a     = zeros(LoopNum, st_num);
Para_b     = zeros(LoopNum, st_num);
Para_g     = zeros(LoopNum, st_num);
Para_p     = zeros(LoopNum, st_num);
Para_q     = zeros(LoopNum, st_num);
Para_w     = zeros(LoopNum, st_num);
Para_d     = zeros(LoopNum, st_num);
% ***
tic
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
			bet_crt = ini_bet;
			gam_crt = ini_gam;
			p_crt   = ini_p;
			q_crt   = ini_q;
			w_crt   = ini_w;
			Ds_crt  = ini_Ds;
		end
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
% Updating-Step 3: beta (bet_)
		bet_rnd = rndmat_b(st,lp);
		r_A = LikelihoodFunc(k_st,sig*Ds_crt,Zobs_st,...
			alp_crt,bet_rnd,gam_crt,p_crt,q_crt,w_crt,...
			AP_st,DDF_st,Ds_crt,nf_st,rhod_st,sand_st,ths_st);
		r_B = LikelihoodFunc(k_st,sig*Ds_crt,Zobs_st,...
			alp_crt,bet_crt,gam_crt,p_crt,q_crt,w_crt,...
			AP_st,DDF_st,Ds_crt,nf_st,rhod_st,sand_st,ths_st);
		ratio = min(r_A / r_B, 1);
		U = unifrnd(0,1);
		if ratio > U
			bet_crt = bet_rnd;
		end
% Updating-Step 4: q
		q_rnd = rndmat_q(st,lp);
		r_A = LikelihoodFunc(k_st,sig*Ds_crt,Zobs_st,...
			alp_crt,bet_crt,gam_crt,p_crt,q_rnd,w_crt,...
			AP_st,DDF_st,Ds_crt,nf_st,rhod_st,sand_st,ths_st);
		r_B = LikelihoodFunc(k_st,sig*Ds_crt,Zobs_st,...
			alp_crt,bet_crt,gam_crt,p_crt,q_crt,w_crt,...
			AP_st,DDF_st,Ds_crt,nf_st,rhod_st,sand_st,ths_st);
		ratio = min(r_A / r_B, 1);
		U = unifrnd(0,1);
		if ratio > U
			q_crt = q_rnd;
		end
% Updating-Step 5: gamma (gam_)
		gam_rnd = rndmat_g(st,lp);
		r_A = LikelihoodFunc(k_st,sig*Ds_crt,Zobs_st,...
			alp_crt,bet_crt,gam_rnd,p_crt,q_crt,w_crt,...
			AP_st,DDF_st,Ds_crt,nf_st,rhod_st,sand_st,ths_st);
		r_B = LikelihoodFunc(k_st,sig*Ds_crt,Zobs_st,...
			alp_crt,bet_crt,gam_crt,p_crt,q_crt,w_crt,...
			AP_st,DDF_st,Ds_crt,nf_st,rhod_st,sand_st,ths_st);
		ratio = min(r_A / r_B, 1);
		U = unifrnd(0,1);
		if ratio > U
			gam_crt = gam_rnd;
		end
% Updating-Step 6: w
		w_rnd = rndmat_w(st,lp);
		r_A = LikelihoodFunc(k_st,sig*Ds_crt,Zobs_st,...
			alp_crt,bet_crt,gam_crt,p_crt,q_crt,w_rnd,...
			AP_st,DDF_st,Ds_crt,nf_st,rhod_st,sand_st,ths_st);
		r_B = LikelihoodFunc(k_st,sig*Ds_crt,Zobs_st,...
			alp_crt,bet_crt,gam_crt,p_crt,q_crt,w_crt,...
			AP_st,DDF_st,Ds_crt,nf_st,rhod_st,sand_st,ths_st);
		ratio = min(r_A / r_B, 1);
		U = unifrnd(0,1);
		if ratio > U
			w_crt = w_rnd;
		end
% Updating-Step 7: Ds
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
		Para_b(lp,st) = bet_crt;
		Para_g(lp,st) = gam_crt;
		Para_p(lp,st) = p_crt;
		Para_q(lp,st) = q_crt;
		Para_w(lp,st) = w_crt;
		Para_d(lp,st) = Ds_crt;
% Process display
		if mod(lp, LoopNum/5) == 0
			disp(['Finished: ' num2str(lp) ' of ' num2str(LoopNum) ' Loops'...
			' | st-' num2str(st) ' of ' num2str(st_num) ' Sites' ...
			' | Prior-' num2str(ipr)])
		end
	end
	disp(['----- Finish loops of st-' num2str(st) ...
		  ' | Prior-' num2str(ipr) ' -----'])
end % end of parfor
toc
% save the output sliced variable
for st = 1:st_num
	Para(:,pa_num*st_num*(ipr-1)+pa_num*(st-1)+1) = Para_a(:,st);
	Para(:,pa_num*st_num*(ipr-1)+pa_num*(st-1)+2) = Para_b(:,st);
	Para(:,pa_num*st_num*(ipr-1)+pa_num*(st-1)+3) = Para_g(:,st);
	Para(:,pa_num*st_num*(ipr-1)+pa_num*(st-1)+4) = Para_p(:,st);
	Para(:,pa_num*st_num*(ipr-1)+pa_num*(st-1)+5) = Para_q(:,st);
	Para(:,pa_num*st_num*(ipr-1)+pa_num*(st-1)+6) = Para_w(:,st);
	Para(:,pa_num*st_num*(ipr-1)+pa_num*(st-1)+7) = Para_d(:,st);
end
end % end of priors (ipr)
% ***
save([out_dir 'sens_prior_para.mat'],'Para',...
	'st_num','pa_num','pr_num','lhos','SiteNo');
% ***
matlabpool close; % close multi-thread
