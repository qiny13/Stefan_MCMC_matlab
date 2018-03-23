% Date: 2017-03-09 | Author: Y.Qin
% update: 17-03-27 | ~,alpha,beta -> alpha,beta,gamma
% LikelihoodFunc.m:
%	calculate the likelihood of maximum thickness of seasonally frozen 
%	ground (MTSFG) by using the Stefan Solution and observation
% Function(s) Needed:
%   - StefanFunc.m
function L = LikelihoodFunc(k,sig,Zobs, ...
							alpha,beta,gamma,p,q,w,AP,DDF,Ds,nf,rhod,sand,ths)
% ***
% Initialization of Constant
pi    = 3.141593;		% pi constant
V     = 0;				% reset var: variance of sim and obs
yrgp  = 1967 - 1961;	% year gap of obs(1967-2015) and sim(1961-2016)
[~,n_sample] = size(Zobs);
% ***
% 1. Calculate the variance: sigma(Zobs(i)-Zsim(i))^2
for i_sam = 1:n_sample
	if Zobs(i_sam) > 0 && DDF(i_sam + yrgp) > 0
		Zsim = StefanFunc(alpha,beta,gamma,p,q,w, ...
						  AP(:,i_sam + yrgp),DDF(i_sam + yrgp), ...
						  Ds,nf,rhod,sand,ths);
		V = V + (Zobs(i_sam) - Zsim)^2;
	end
end
% 2. Calculate the likelihood
L = (2*pi)^(-k/2) * sig^(-k) * exp(-1/(2*sig^2) * V);
% ***
end  % end of function
