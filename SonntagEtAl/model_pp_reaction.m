function dcdt = model_pp_reaction(t,c,pars)
%MODEL_PP_REACTION Summary of this function goes here
%   input c is in [g/l]
%   M is in [g/umol]
cPP=c(1)./pars.M_PP; % umol/l
cP=c(2)./pars.M_P;
cOH=1; %c(3)./pars.M_LiOH;

r1=pars.k1*cOH*cPP; % k1 units: [l/umol/day]
rdeg=pars.kdeg*cP*cOH;
% Define Equations  
dcPPdt  = -r1; 
dcPdt = r1-rdeg;
dcOH = -r1-rdeg;
dcdt=[dcPPdt.*pars.M_PP;dcPdt.*pars.M_P;dcOH.*pars.M_LiOH]; %output in g/l
end

