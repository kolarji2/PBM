function [exps,TimeEnd45,TimeEnd55,TimeEnd65,timeEndData,prop] = load_exp_data_telmisartan(par)
%load_exp_data_telmisartan Summary of this function goes here
%   Detailed explanation goes here

basedir='MuzikEtAl/expdata/'; % default path to files

Nexp=6;
exps=cell(Nexp,1);



for i=1:Nexp
    data=readtable(fullfile(basedir, 'telmiexps.xlsx'), ...
        'Sheet',['data' num2str(i)],'VariableNamingRule','modify');
    exps{i}={};
    exps{i}.cT=data.cA_g_l_/par.MtelH;% mol/l
    exps{i}.cTstd=data.cstd/par.MtelH;% mol/l
    exps{i}.time=data.t_min_;
    exps{i}.cS0=data.cT0(1)/par.MtelH;% mol/l
    exps{i}.cB0=data.cB0(1)/par.MNa2CO3;% mol/l
    exps{i}.T=data.T(1)+273.15;
    sel=~isnan(data.tpH_min_);
    exps{i}.time_pH=data.tpH_min_(sel);
    exps{i}.pH=data.pH(sel);
    exps{i}.conc_units='moll';
end

TimeEnd45= @(cB) 5820.3384./(cB*par.MNa2CO3).^1.28963; %cB % mol/l; time [minutes]
TimeEnd55= @(cB) 180.08./(cB*par.MNa2CO3).^0.5926; %cB % mol/l; time [minutes]
TimeEnd65= @(cB) -7.6195*log(cB*par.MNa2CO3) + 35.4875; %cB % mol/l; time [minutes]
timeEndData={};
timeEndData.cBmin=12/par.MNa2CO3;% mol/l
timeEndData.cBmax=48/par.MNa2CO3;% mol/l
timeEndData.cTel=30/par.MtelH; % mol/l 0.6g Telmi do 20ml

prop = load_viscosity_telmisartan();
end

