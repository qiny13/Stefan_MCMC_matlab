% Date: 2017-08-09 | Author: Y.Qin
% update: 18-02-22 | API(alpha,p) ->> API(alpha,p, Ds)
%                    include sens_model_*.m 2 in 1
% sens_model_selct_ds.m:
%	Sensitive analysis of the Stefan model selection
%   i.e., one, two, or three term(s)
%   using the AIC & BIC evaluation metrics
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
%	- num_obs   : number of valid data in each site
load([matindir 'Site_obs_mtsfg.mat']);
% Load the initial parameters
ParaIni = load('ParaIni.txt');
% ***
% Load the initial parameters
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
[pa_num,~] = size(ParaIni);
%	- record values of parameters in each loop
Para_a     = zeros(LoopNum, st_num);
Para_b     = zeros(LoopNum, st_num);
Para_g     = zeros(LoopNum, st_num);
Para_p     = zeros(LoopNum, st_num);
Para_q     = zeros(LoopNum, st_num);
Para_w     = zeros(LoopNum, st_num);
Para_d     = zeros(LoopNum, st_num);
%	- record AIC/BIC/RMSE at each station for each k-fold
rmse       = zeros(st_num, 3);
aicbic     = zeros(st_num, 3*2);
% ***
% Model selection: term 1 to 3
for term = 1:3
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
		if term == 1
			bet_crt=0; q_crt=0; gam_crt=0; w_crt=0;
		elseif term == 2
			gam_crt=0; w_crt=0;
		end
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
	if term >= 2
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
	if term >= 3
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
	end % term >= 3
	end % term >= 2
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
			' | Term-' num2str(term)])
		end
	end % end of loop
	disp(['----- Finish loops of st-' num2str(st) ...
		  ' | Term-' num2str(term) ' -----'])
end % end of parfor
toc
% ***
% VALIDATE process-1: AIC/BIC
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
% GET station input (_st) obs data in *.mat
	Zobs_st = mtsfg_obs(strow, :);
	k_st    = num_obs(strow);
% Maximum Likelihood: iterative loop of para_a,b,g,p,q,w
	Lkhd = zeros(LoopNum/2, 1);
	for lh = LoopNum/2+1 : LoopNum
		p_a = Para_a(lh,st);
		p_b = Para_b(lh,st);
		p_g = Para_g(lh,st);
		p_p = Para_p(lh,st);
		p_q = Para_q(lh,st);
		p_w = Para_w(lh,st);
		p_d = Para_d(lh,st);
		Lkhd(lh-LoopNum/2, 1) = LikelihoodFunc(...
			k_st,sig*p_d,Zobs_st,...
			p_a, p_b, p_g, p_p, p_q, p_w,...
			AP_st,DDF_st,p_d,nf_st,rhod_st,sand_st,ths_st);
	end
% Calculate AIC and BIC, and save the results
	AIC = -2*log(max(Lkhd)) + (term*2)*2;
	BIC = -2*log(max(Lkhd)) + (term*2)*log(k_st);
	aicbic(st, 2*term-1) = AIC;
	aicbic(st, 2*term)   = BIC;
end
% ***
% VALIDATE process-2: RMSE
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
% Optimized: para_a,b,g,p,q,w
	m_a = median(Para_a(lhos,st));
	m_b = median(Para_b(lhos,st));
	m_g = median(Para_g(lhos,st));
	m_p = median(Para_p(lhos,st));
	m_q = median(Para_q(lhos,st));
	m_w = median(Para_w(lhos,st));
	m_d = median(Para_d(lhos,st));
% Loop of testing samples
	Zsim_st = zeros(1,length(Zobs_st));
	for iv = 1:length(Zobs_st)
		if Zobs_st(iv) > 0 && DDF_st(iv) > 0
			Zsim_st(iv) = StefanFunc( ...
				m_a, m_b, m_g, m_p, m_q, m_w, ...
				AP_st(:,iv),DDF_st(iv), ...
				m_d,nf_st,rhod_st,sand_st,ths_st);
		else
			Zobs_st(iv) = NaN;
			Zsim_st(iv) = NaN;
		end
	end
	ind_nonan = ~isnan(Zobs_st);
	rmse(st,term) = ...
		sqrt(mean((Zsim_st(ind_nonan)-Zobs_st(ind_nonan)).^2));
end
end % end of term (1:3)
% ***
matlabpool close; % close multi-thread
save([out_dir 'sens_model_aicbic.mat'],'aicbic');
save([out_dir 'sens_model_rmse.mat'],'rmse');
