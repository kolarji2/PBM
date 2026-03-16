%
%   Example of usage of modified PBM class for 
%		modelling of pH dependent dissolution of telmisartan.
%		
%		Modified @PBM class can be found in directory: @PBMTelmi 
%
% Show initial Particle size distribution and fit using lognormal PSD used
% for Telmisartan particles
%

m=PBMTelmi();
if ~exist('TimeEnd45','var')
    [exps,TimeEnd45,TimeEnd55,TimeEnd65,TimeEndData,prop] = load_exp_data_telmisartan(m);
    [phi0expdata] = load_psd_telmisartan();
end

m=PBMTelmi();

w=[0,0,1];
m.make_loggrid(-8,-2,300);
phi0exp_data=cellfun(@(x) x.phi0exp,phi0expdata,'UniformOutput',false);
dexp_data=cellfun(@(x) x.dexp,phi0expdata,'UniformOutput',false);
% Plot all three PSD and corresponding fit
% Final combination of the PSDs given by w is not displayed as is simply a
% linear combination of the individual PSDs
dispfit=true
set_init_distFromExpData(m,dexp_data,phi0exp_data,w,[],dispfit);


