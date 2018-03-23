% Date: 2017-03-06 | Author: Y.Qin
% update: 17-03-27 | ~,alpha,beta -> alpha,beta,gamma
% update: 17-04-04 | th_liq - (0.00, 0.99*ths)
% StefanFunc.m:
%	calculate the maximum thickness of seasonally frozen ground
%	(MTSFG) by using the Stefan Solution with physical parameters
% References:
%	Qin et al.2016| DOI:10.1016/j.jhydrol.2016.09.008 | Eq(1)
%	Li and Koike.2003| DOI:10.1016/S0165-232X(03)00009-0 | Eq(20)
function [Z,th_liq] = StefanFunc(alpha,beta,gamma,p,q,w,AP,DDF,Ds,nf,rhod,sand,ths)
% ***
% Initialization of Constant
lamd_q  = 7.70;   % quartz (Johansen, 1975)
lamd_o  = 2.00;   % others (Johansen, 1975)
lamd_w  = 0.57;   % Ming-ko Woo BOOK
lamd_i  = 2.20;   % Ming-ko Woo BOOK
rho_w   = 1000;   % water density (kg/m3)
L       = 334000; % phase change heat (J/kg)
tao     = 86400;  % seconds per day
% ***
% 0. Calc the annual liquid soil water content (Qin,2016)
AP      = AP/1000;   % unit: from mm to m
th_liq  = alpha*((AP(1)/Ds)^p) + beta*((AP(2)/Ds)^q) + gamma*((AP(3)/Ds)^w);
%       - ADD boundary of th_liq (0.00, 0.99*ths)
if th_liq > 0.999*ths
	th_liq = 0.999*ths; % SET upper boundary
elseif th_liq < 0.001*ths
	th_liq = 0.001*ths; % SET lower boundary
end
% 1. Kf - thermal conductivity of frozen soil (Fakoui,1981)
% delta -percentage of quartz (50% of sand content)
delta   = sand/100 *0.5;
% solids thermal conductivity (Johansen, 1975)
lamd_s  = (lamd_q^(delta))*(lamd_o^(1-delta));
% K of saturated frozen soil (Fakoui,1981)
Ksat    = (lamd_s^(1-ths)) * (lamd_w^(th_liq)) ...
	* (lamd_i^(ths-th_liq));
% K of dry frozen soil (Fakoui,1981)
Kdry    = (0.135 * 1000*rhod + 64.7) / (2700 - 0.947 * 1000*rhod);
% Kersten number (Fakoui,1981)
Ke      = th_liq / ths;
% Kf - frozen soil thermal conductivity (Fakoui,1981)
Kf      = Ke*(Ksat-Kdry) + Kdry;
% 2. Calculate the E Value: edaphic factor
% 	- based on Qin et al.,2016, Formula(4-9)
A = 2 * Kf * nf * tao;
B = rho_w * th_liq * L;
% StefanE (Qin,2016; Li and Koike,2003)
StefanE = sqrt(A/B);
% 3. Calculate annual MTSFG(Z) : based on StefanE & DDF
Z = StefanE * sqrt(DDF);
% ***
end  % end of function
