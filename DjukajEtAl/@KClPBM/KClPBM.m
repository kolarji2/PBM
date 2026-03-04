classdef KClPBM < PBM
    %% Population Balance Model (PBM) for Dissolution of KCl
    %   This is modification of PBM class, which was used for analysis
    %   of KCl dissolution. Results published in:
    %   S. Djukaj, J. Kolář, R. Lehocký, A. Zadražil, F. Štěpánek,
    %   Design of particle size distribution for custom dissolution
    %   profiles by solving the inverse problem, Powder Technology, 2021
    %   https://doi.org/10.1016/j.powtec.2021.10.023
    %% About KClPBM
    %  
    %   This class is used just for initialization of the original PBM
    %   class with compound specific and experiment set-up specific
    %   properties.
    %
    %   Namely: 
    %
    %   KCl: 
    %       density
    %       relation between conductivity and KCl concentration 
    %
    %   Solvent: Water + isopropanol
    %           Takes into account:
    %               - alcohol volumetric contraction
    %               - nonideal viscosity of the mixture
    %
    %   Functions for loading experimental data from excel files
    %    
    
    %% KClPBM Properties
    properties
        
    end
    
    methods
        function obj=KClPBM(infoi,datai)
            % KClPBM(obj,infoi,datai) - initialize PBM class for KCl
            % dissolution
            % infoi - info about experiment
            % datai - contains initial particle size distribution
            %
            T=25+273.15; 
            rho_water=998; %kg/m3 == g/l
            rho_ipa=785; %kg/m3  == g/l https://pubcheobj.ncbi.nlobj.nih.gov/compound/Isopropyl-alcohol#section=Density
            rhos = 1980;%kg/m3 drug=1200;

            %% Volume fraction of IPA and water
            phi0_water=infoi.phiw;
            phi0_ipa=infoi.phiipa;
            %% Viscosity
            etafit=@(phi_ipa) (-6.2413*phi_ipa.^2+7.2342*phi_ipa+1.0607)*1e-3;
            etaf = etafit(phi0_ipa);

            %% Correction of initial concentration (alcohol contraction)
            V0=infoi.V*1e-3; %l
            [Vmix,rhomix]=vol_contraction_ipa(phi0_ipa*V0,phi0_water*V0,rho_ipa,rho_water); %l, g/l
            V0=(Vmix+infoi.mtot./rhos);
            cs0=infoi.mtot/V0; %[g/l]

            %% For Diffusivity estimation (Wilke-Chang method)
            M_water=18.015;%g/mol
            M_ipa=60.096; %g/mol    
            x_water=2.6; %water association parameter
            x_ipa=1.5; % propanol association parameter (unknown) can be 1
            MKCl=74.55; %g/mol
            Vmol_solute=MKCl/rhos*1e3;
            n_water=phi0_water*V0*rho_water/M_water;
            n_ipa=phi0_ipa*V0*rho_ipa/M_ipa;
            y_water=n_water/(n_water+n_ipa);
            y_ipa=1-y_water;
            M_solvent=y_water*M_water+y_ipa*M_ipa;
            x=(x_water*y_water*M_water+x_ipa*y_ipa*M_ipa)/M_solvent;

            %% Update class
            obj.T=T;
            obj.Mf=M_solvent;
            obj.Ms=MKCl;
            obj.assoc_par=x;
            obj.Vmol_solute=Vmol_solute; % cm3/mol
            obj.csat=13.79; % cond2conc(2030)
            obj.rhos=rhos;
            obj.rhof=rhomix;
            obj.etaf=etaf;
            obj.diff=obj.calc_diff_wilkechang();
            
            % Method for estimation of initial PSD from experimental data
            obj.dist_est_method='2lognormal';
            %% Init model
            obj.make_loggrid(-7,-2,400);
            obj.set_init_distFromExpData(datai.all.dexp,datai.phi0exp,infoi.w,[],false);
            obj.set_initial_conditions(0,cs0,'gl');
        end
        
        function set_optimized_pars(obj)
            %% Optimized model pars
            obj.epsilon=8.50e-2;
            a1=0.5;
            a2=2.71e-01;
            obj.fAfun=@(r,x) pi*(1.0+a1*exp(-r/a2).*ones(size(x)));
            obj.Sh_method='kolmogorov';            
        end
        
    end
end