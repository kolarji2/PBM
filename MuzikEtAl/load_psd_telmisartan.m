function [phi0expdata] = load_psd_telmisartan()
%LOAD_PSD_TELMISARTAN Summary of this function goes here
%   Detailed explanation goes here
fname='MuzikEtAl/expdata/telmisartan_particlesize.xlsx';
data=readtable(fname,'sheet','phi');
xbounds=[0;data.diam_m_*1e-6]; %m

utils=Utils;
xcexp=utils.get_xclass_fromBounds(xbounds);
phiexp=(data.tween1+data.tween2)*0.5;
phiexp=phiexp/sum(phiexp);
phiexpAgg=data.h2o;
phiexpAgg=phiexpAgg/sum(phiexpAgg);
phiexpAgg2=data.agg;
phiexpAgg2=phiexpAgg2/sum(phiexpAgg2);


phi0expdata=cell(3,1);

phi0expdata{1}.phi0exp=phiexp;
phi0expdata{2}.phi0exp=phiexpAgg;
phi0expdata{3}.phi0exp=phiexpAgg2;

for i=1:length(phi0expdata)
   phi0expdata{i}.dexp=data.diam_m_*1e-6; 
   phi0expdata{i}.xcexp=xcexp;
   phi0expdata{i}.dxexp=utils.get_dx_fromBounds(xbounds);   
end

end

