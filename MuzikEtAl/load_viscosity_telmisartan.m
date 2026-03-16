function prop = load_viscosity_telmisartan()
%load_viscosity_telmisartan Summary of this function goes here
%   Detailed explanation goes here
basedir='MuzikEtAl/expdata/'; % default path to files
data_visc=readtable(fullfile(basedir, 'telmiexps.xlsx'),'Sheet','viscosity','VariableNamingRule','modify');
%data_visc.c_g_l_;
%data_visc.T_C_;
%data_visc.etaf_Pas_;
p = polyfit(data_visc.T_C_+273.15,log10(data_visc.etaf_Pas_),2);
prop={};
prop.etafTelmi230=@(x) 10.^(polyval(p,x));

% figure
% plot(data_visc.T_C_+273.15,data_visc.etaf_Pas_,'or')
% hold on 
% Tarr=linspace(25,66,100)+273.15;
% plot(Tarr,prop.etafTelmi230(Tarr))

end

