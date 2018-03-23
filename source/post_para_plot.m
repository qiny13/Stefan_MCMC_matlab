% Date: 2017-03-10 | Author: Y.Qin
% update: 17-03-27 | ~,alpha,beta ->> alpha,beta,gamma
% update: 17-04-13 | optmiz: mode ->> median
% update: 17-04-17 | save both: mode + median (delete optmiz)
% update: 17-05-26 | solo chain ->> ensemble para dataset
% update: 17-06-05 | xlim(xl) ->> read 'ParaIni' to set xlim
% update: 18-02-11 | pa_num = 6->7, API(alpha,p) ->> API(alpha,p,Ds)
% post_para_plot.m:
%	Post-process main_stefan_mcmc.mat and plot parameter pdf figures
%clear
% ***
% Set workspace direction
root_dir  = '..\';
out_dir   = [root_dir 'results\'];
chain_dir = [out_dir  'multi_chain\'];
% Go through the mat files
dirs = dir([chain_dir '*.mat']);
% Struct to cell | Transpose to column
dircell = struct2cell(dirs)' ;
N = dircell(:,1);
Chain_num = size(N,1);
% ***
% [*Chain-by-Chain Plot*]
tic
for ichain = 1:Chain_num
% LOOP of chain.mat files
	filename = char(N(ichain));
% Load the main_stefan_mcmc result
%	- Para : row = LoopNum | col = (alp,bet,gam,p,q,w,Ds) x stations
%	- SiteNo : row = st_num
%   - ParaIni: row = pa_num, col(1,2,3) = ini_value,low_limit,high_limit, 
	load([chain_dir filename]);
	[LoopNum, pa_col] = size(Para);
	[st_num, ~]       = size(SiteNo);
	pa_num = pa_col/st_num;
	rank   = 1 : LoopNum; % rank: full sample
	lhos   = LoopNum/2+1 : LoopNum; % lhos: last half of samples
	ParaNameList = {'alpha','beta','gamma','p','q','w','Ds'};
% Initialize of metrics matrix
	if ichain == 1
		median = zeros(pa_num,st_num);  % sample median
		pdmode = zeros(pa_num,st_num);  % pdf mode (peak of density)
		standd = zeros(pa_num,st_num);  % standard deviation
		Ensemb = zeros(LoopNum/2*Chain_num, pa_col); % ensembled para
	end
% Save to the ensemble para dataset
	Ensemb(LoopNum/2*(ichain-2)+lhos, :) = Para(lhos, :);
% LOOP of stations and ipara
	for st = 1:st_num
	for ipa = 1:pa_num
% Set : get data samples for station:st, parameter:ipa
		sample = Para(rank, pa_num*(st-1)+ipa); % full sample
		smppdf = Para(lhos, pa_num*(st-1)+ipa); % last half of samples (lhos)
		figure('visible','off')
% plot-1 para_value-rank time series
		subplot(2,2,[1 2]);
		plot(rank, sample)
		xlabel('Number of iteration')
		ylabel(ParaNameList{ipa})
% plot-2 hist graph
		subplot(2,2,3);
		hist(smppdf, 50)
		[counts,center] = hist(smppdf, 50);% get hist parameters
		[max_h,index_h] = max(counts);	   % max of density (mode)
		line([center(index_h) center(index_h)],...
			 [0 max_h],'Color','r')
		xlabel(['Parameter: ' ParaNameList{ipa}])
		ylabel('Frequency')
% plot-3 pdf graph	
		subplot(2,2,4);
		[ksd_f,ksd_x]   = ksdensity(smppdf);% kernel-smooth
		[max_d,index_d] = max(ksd_f);	    % get max of density
		plot(ksd_x,ksd_f,'LineWidth',1)
		medi = quantile(smppdf,0.50); % 50-percentile of samples (median)
		mode = ksd_x(index_d); % max of density (mode)
		optm = medi; % optmiz for plotting [median .or. mode]
		line([optm optm],[0 max_d],'Color','r')
		xlim([ParaIni(ipa,2),ParaIni(ipa,3)]);	% set the x-axis min/max
		xlabel(['Parameter: ' ParaNameList{ipa}])
		ylabel('Density')
% OUTPUT of figures
		figname = [out_dir 'st' num2str(SiteNo(st)) ...
			'_para_' ParaNameList{ipa}(1) ...
			'_c' num2str(ichain) '.png'];
		saveas(gcf, figname);
	end
	disp(['ParaPlot: chain-' num2str(ichain) ' of ' num2str(Chain_num) ' chains'...
		' st-' num2str(st) ' of ' num2str(st_num) ' sites'])
	end
	toc
end
% ***
% [*Ensemble*]
% LOOP of stations and ipara
disp('Plotting the ensemble figure ...')
for st = 1:st_num
	for ipa = 1:pa_num
		esmpdf = Ensemb(:, pa_num*(st-1)+ipa); % last half of samples (lhos)
		[ksd_f,ksd_x]   = ksdensity(esmpdf);% kernel-smooth
		[max_d,index_d] = max(ksd_f);	    % get max of density
% plot-0 ensemble figure plot
		figure('visible','off')
		plot(ksd_x,ksd_f,'LineWidth',3)
		medi = quantile(esmpdf,0.50); % 50-percentile of samples (median)
		mode = ksd_x(index_d); % max of density (mode)
		optm = medi; % optmiz for plotting [median .or. mode]
		line([optm optm],[0 max_d],'Color','r','LineWidth',3)
		xlim([ParaIni(ipa,2),ParaIni(ipa,3)]);	% set the x-axis min/max
		xlabel(['Parameter: ' ParaNameList{ipa}],'FontSize',28)
		ylabel('Density','FontSize',28)
        set(gca,'FontSize',28);
% OUTPUT of figures
		figname = [out_dir 'st' num2str(SiteNo(st)) ...
			'_para_' ParaNameList{ipa}(1) '_esm.png'];
		saveas(gcf, figname);
% OUTPUT of statistics of parameters
		median(ipa,st)  = medi;
		pdmode(ipa,st)  = mode;
		standd(ipa,st)  = std(esmpdf); % std deviation of sample
	end
end
toc
% OUTPUT of xls: save mode, mid, and std of each parameter
xlswrite([out_dir 'midstd.xlsx'],pdmode,1, 'B2')
xlswrite([out_dir 'midstd.xlsx'],median,1,'B10')
xlswrite([out_dir 'midstd.xlsx'],standd,1,'B18')
% OUTPUT of *.mat: for post_valid.m
save([out_dir 'post_para.mat'],'Ensemb','SiteNo');
