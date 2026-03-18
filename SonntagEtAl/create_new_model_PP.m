function m = create_new_model_PP(T,M_PP)
    %% Create new model instance
    m=PBM;
    
    m.rhof    = 977.8;       % water density @70 °C [kg/m^3] (Holecek etabulky)
    m.rhos    = 1200;        % PP density [kg/m^3] (http://www.chemspider.com/Chemical-Structure.8028457.html)
    etaf_water_fun= @(T) 0.1/(2.20065*(T-282.92341+(8761.27+(T-282.92341)^2)^.5)-129.908);

%%   Diffusivity estimation (Wilke Chang method)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Molar volume estimation based on:
    %               Technicka termodynamika Basarova, Haidl)
    %           elektronický učební materiál pro výuku předmětu
    %               Technická termodynamika N409023
    %
    % Paliperidone palmitate formula: C39H57FN4O4
    %   C = +39*14.8
    %   H = +57*3.7
    %   F = +1*8.7
    % N4:
    %   2x double bond N = +2*15.6
    %   2x ternary amin  =  +2*12   (taken as secondary, unkwon value for ternary)
    % O4:
    %   1x higher ester = +11
    %   2x acid =0     =  +2*12
    %   1x bonded with N = +8.3
    % 1x cycle-naftalen = -30
    % 2x cycle-6 = -2*15
    % 1x cycle-5 = -1*11.5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Vmol_solute=39*14.8+57*3.7+1*8.7+2*15.6 +2*12+11+2*12+8.3-30 -2*15-1*11.5; % cm3/mol    
    %% Update class
    m.T=T;
    m.Mf=18.015;%g/mol Water
    m.Ms=M_PP*1e6;
    m.assoc_par=2.6; %for Water
    m.etaf = etaf_water_fun(T);
    m.Vmol_solute=Vmol_solute;
    
    m.impeler_rpm=180;
    m.impeler_diam=26.2e-3;
    m.tank_diam=70e-3;
    % Levins pars    
    m.DsDt=m.impeler_diam./m.tank_diam;
    m.diff=m.calc_diff_wilkechang();
end

