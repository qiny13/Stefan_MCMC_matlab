% Date: 2017-09-15 | Author: Y.Qin
% Update: 18-03-19 | Ds = 2.0 ->> 'krg_ds'
% post_valid_kriging.m:
%	Post-process load kriging_data_input.mat and generate the
%   "kriginged" spatial MTSFG from MCMC results
%   (adpated from 'post_valid.m')
% SELECTION: *** idw_* or krg_* (at 4 places) ***
clear
% ***
% Set workspace direction
root_dir = '..\';
mat_dir = [root_dir 'kriging\'];
out_dir = [root_dir 'results\'];
asc_dir = [root_dir 'kriging\input\'];
% Input Constant
%	- SET: average soil depth (m)
IniYear  = 1981;
EndYear  = 2015;
yr_num   = EndYear-IniYear+1;
% Load the stations climatic forcing and land surface parameters
%	'IniYear','EndYear','Xnum_grids','Ynum_grids', ...
%	'basin','rhod','sand','ths', ...
%	'idw_alp','idw_nf','idw_p','idw_pexp', ...
%	'krg_alp','krg_nf','krg_p','krg_pexp','krg_ds', ...
%	'DDF','AP1','AP2','AP3'
load([mat_dir 'kriging_data_input.mat']);
% Load the annual max frozensoil depth from gbehm simul
load([mat_dir 'FRS_sim_YRSR_minmax.mat'], 'frs_ymax');
% Load the permafrost zone (1-permafrost, 0-non)
pmfr   = arcgridread([asc_dir 'frs_pmfr_d2.asc']);
% Initialize of metrics matrix
%	- simulation of each year in each loop
Zmc_yr = zeros(Ynum_grids, Xnum_grids, EndYear-IniYear);
Zgb_yr = frs_ymax;
Mat_0  = zeros(Ynum_grids, Xnum_grids);
plot_zmc = Mat_0;
plot_zgb = Mat_0;
% ***
% Model selection: idw_* or krg_* (at 4 places)
alp_cal = krg_alp;  p_cal = krg_p;      ds_cal = krg_ds;
bet_cal = Mat_0;	q_cal = Mat_0;
gam_cal = Mat_0;	w_cal = Mat_0;
nf_cal = krg_nf;	ap_cal = AP1;		ddf_cal = DDF;
rhod_cal = rhod;	sand_cal = sand;	ths_cal = ths;
% ***
% Loop of Years
for stryr = IniYear+1 : EndYear % from 1981+1 to 2015
iyr = stryr - IniYear; % from 1 to 34
% Loop of grids (x by y)
for m = 1:Ynum_grids
    for n = 1:Xnum_grids
        if basin(m,n) > 0
            Zmc_yr(m,n,iyr) = StefanFunc2( ...
                alp_cal(m,n),bet_cal(m,n),gam_cal(m,n), ...
                p_cal(m,n),q_cal(m,n),w_cal(m,n), ...
                ap_cal(m,n,iyr),ddf_cal(m,n,iyr), ...
                ds_cal(m,n),nf_cal(m,n),rhod_cal(m,n), ...
                sand_cal(m,n),ths_cal(m,n));
            % Calculate average frozen depth
            plot_zmc(m,n) = plot_zmc(m,n) + ...
                Zmc_yr(m,n,iyr)  /(EndYear-IniYear);
            plot_zgb(m,n) = plot_zgb(m,n) + ...
                Zgb_yr(m,n,iyr+1)/(EndYear-IniYear);
        else % Out of the basin
            Zmc_yr(m,n,iyr) = -9999;
        end
    end
end
disp(['Finish: ' num2str(stryr)])
end % END of year loop
% Mask of permafrost zone
plot_zmc(pmfr>0) = 100;	plot_zmc(isnan(basin)) = -9999;
plot_zgb(pmfr>0) = 100;	plot_zgb(isnan(basin)) = -9999;
% ***
% PLOT ASC Files
MatToAsc(plot_zmc, [out_dir 'Zmc_avg_map'], ...
	270, 170, -440000, -1460000, 3000, -9999);
MatToAsc(plot_zgb, [out_dir 'Zgb_avg_map'], ...
	270, 170, -440000, -1460000, 3000, -9999);
% OUTPUT of *.mat
save([out_dir 'post_valid_kriging.mat'], ...
	'Zmc_yr','Zgb_yr','plot_zmc','plot_zgb');
