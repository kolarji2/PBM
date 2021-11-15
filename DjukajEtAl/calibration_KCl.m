
%% Create calibration of KCl
rho_water=998; %kg/m3 == g/l
rho_ipa=785; %kg/m3  == g/l https://pubchem.ncbi.nlm.nih.gov/compound/Isopropyl-alcohol#section=Density
rhos = 1980;%kg/m3 drug=1200;
[cond2conc,conc2cond] = build_cond2conc_interp(rho_ipa,rho_water,rhos,'DjukajEtAl/expdata/calibrationKCl.xlsx');
save('DjukajEtAl/expdata/conductivity.mat','cond2conc','conc2cond')

function [cond2conc,conc2cond] = build_cond2conc_interp(rho_ipa,rho_water,rhos,calib_file)
%COND2CONC Return interpolator for calculateion of concentration
%           from conductivity
%
%   Use conductivity model or interpolation of experimental points.
%   Assumes conductivity in same units as experimental data.
%   Returnes concnetratin in same units as experimental data.
%   
%   
%    Model is based on Fuoss-Onsager equation
%   model_associated: finished
%   model_unassociated: TODO
%
%
%   
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conductivity_method='model_associated'; % 'interp' 'model'
cmax = 20; % for model to create inverse function easily
Ninterppoints=10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Suada's old calib data
% Note: check calibration (with 12 g/l max)
calib_a=186.852173913044;
calib_b=172.956521739131;
% Note: check calibration (with 16 g/l max)
% calib_a=145.51;
% calib_b=292.96;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%conc=(conductivity-calib_b)./calib_a; %Concentration [g/l] == [kg/m3]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data preparation
calib_data=readmatrix(calib_file);
ms=calib_data(:,1)*1e-3; %g
Vipa=100e-3; %l
Vw=30e-3; %l
[Vmix,rhomix] = vol_contraction_ipa(Vipa,Vw,rho_ipa,rho_water);
V=Vmix+ms./rhos; % l
cexp=ms./V; %g/l
condexp=calib_data(:,3); % uS/cm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if contains(conductivity_method,'interp')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation
    cond2conc=@(cond) interp1([0;condexp],[0;cexp],cond,'makima','extrap');
    conc2cond=@(conc) interp1([0;cexp],[0;condexp],conc,'makima','extrap');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif contains(conductivity_method,'model_associated')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Classic conductivity model
    eps=1e-20;
    % Fuoss-Onsager equation for associated electrolytes
    K=@(a,c) (a(1)-a(2).*sqrt(c)+a(3).*c.*log(a(4).*c+eps)+a(5).*c).*c./(1+a(6).*c);
    % Fit to experimental data
    objFun= @(x) sum((K(x,cexp)-condexp).^2);
    eqpars=fmincon(objFun,[100,1,1,1,1,1]);
    % Create functions
    conc2cond=@(conc) K(eqpars,conc);
    % Inverse function created by inverting interpolation ranges
    c=linspace(0,cmax,Ninterppoints);    
    cond2conc=@(cond) interp1(conc2cond(c),c,cond,'makima','extrap');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end

