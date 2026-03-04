function [Vmix,rho] = vol_contraction_ipa(Vipa,Vw,rho_ipa,rho_water)
%VOL_CONTRACTION_IPA Calculate real volume and density of water isopropanol mixture
% Vipa,Vw [l]
% rho_ipa,rho_wate [g/l]
% Alcohol volume contraction
% data from Properties of mixtures of isopropyl alcohol and water Robert B. Lebo 1921
wIpa=[100	90.35	85.09	74.35	65.22	53.07	43.02	33.17	21.39	9.58	0]./100; % [-]
rhoMix=[0.78556	0.80866	0.82282	0.84828	0.87003	0.89868	0.92418	0.9459	0.96847	0.98293	0.998].*1e3; % [g/l] [kg/m3] 

%Default densities
if rho_ipa==0
    rho_ipa=785.56; %[g/l] [kg/m3]
end
if rho_water==0
    rho_water=998.0; %[g/l] [kg/m3]
end
% Calculation
mipa=rho_ipa*Vipa;
mw=rho_water*Vw;
mtot=(mipa+mw);
w=mipa/mtot;
rho=interp1(wIpa,rhoMix,w);
Vmix=mtot/rho;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% water:IPA 3:10
% Water IPA mixture     Robert, Lebo Properties of mixture of isopropyl
% alcohol and water 1921
% Wt.   Vol. Density
% 70    76.3 0.8585
% 71    77.2 0.8561
% 72    77.9 0.8538
% 73    78.8 0.8514
%rhof=interp1([76.3,77.2,77.9],[858.5,856.1,853.8],phit_IPA*100); % 76.9231
%w_ipa=phit_IPA*rho_ipa/(phit_IPA*rho_ipa+phi_water*rho_water);

end

