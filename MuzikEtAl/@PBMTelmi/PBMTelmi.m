classdef PBMTelmi < PBM
    %% Population Balance Model (PBM) for Dissolution of Telmisartan
    %   This is modification of PBM class suitable for description of 
    %   pH dependent dissolution
    %% About PBM
    %   Model assumes acid-base equilibrium between 
    %   * all 8 possible protonated/deprotonated states of Telmisartan.
    %   * HCO3(1-), CO3(2-), CO2
    %   * CO2(aq.) and CO2(l)
    %    
    %% PBMTelmi Properties
    properties
        
        % Acid dissociation constants of 3 telmi acidic groups as func of T
        Kafun_T
        
        % Return composition of Telmi forms based on Ka and cH
        compfun_KacH
        
        % Indices and charges of Telmisartan forms
        ind
        
        %% Equillibrium constants
        % Acid dissociation constants of 3 telmi acidic groups
        % In order: T-COOH, N1-H, N2-H
        Ka
        
        % Water eq. constant
        Kw
        
        % activity of H ions
        gH
        
        % Acid dissociation constant of HCO3(1-)
        KaHCO3
        
        % Acid dissociation constant of H2CO3
        KaH2CO3
        
        % CO2 Liquid-gas equillibrium constant
        KgCO2
        
        % CO2 Liquid-gas mass transfer coefficient kM*a
        kMaCO2
        
        % Atmosferic concentration of CO2 [Pa]
        pCO2atm=40; %Pa
        
        % Atmosferic pressure [Pa]
        p0atm=101325;
        
        %% Molar mass
        
        % Molar mass of Telmisartan Na salt [g/mol]         
        MtelNa= 536.6;
        
        % Molar mass of Telmisartan free acid [g/mol]
        MtelH= 514.6;
        
        % Molar mass of Na2CO3 [g/mol]
        MNa2CO3= 23*2 + 12 + 3*16;
        
        % Molar mass of Water [g/mol]
        MH2O=18.015;
        
        %% Density
        rhoTel
        rhoWater
        rhoNa2CO3        
        
        %% Concentration
        % Saturated concentration of neutral forms of telmisartan [mol/l]
        cTintrinsic
        
        %Total concentration of cCO3 in solution [mol/l]
        cCO3tot
        
        %Total concentration of Na ionts in solution [mol/l]
        cNa
    end
    
    %% PBM methods
    
    methods 
        function obj=PBMTelmi(obj)  
            % Constructor of PBMTelmi class
        end
        
        function set_props(obj,T,cNa2CO3,useactivity,etafTelmi230_h)
        % Set properties
            %% Concentrations
            obj.cNa=2*cNa2CO3;
            obj.cCO3tot=cNa2CO3;
            
            %%
            obj.T=T;
            obj.Ms=obj.MtelNa;
            
            %% Density            
            obj.rhoTel=1200; %kg/m3
            obj.rhoWater=999.65+0.20438*(T-273.15)-0.061744*(T-273.15)^1.5; %
            obj.rhoNa2CO3=2530;
            cH2O=obj.rhoWater/obj.MH2O;            
            obj.rhof_method='rhof_fun';
            %obj.rhof_fun_h= @(c) c(1)*obj.MtelNa+c(2)*obj.MNa2CO3+cH2O*obj.MH2O;
            
            obj.rhof_fun_h= @(c) c(1)*obj.MtelNa*(1-obj.rhoWater/obj.rhoTel)+c(2)*obj.MNa2CO3*(1-obj.rhoWater/obj.rhoNa2CO3)+obj.rhoWater;
            obj.rhos=obj.rhoTel;
            
            %% Viscosity
            etafTelmi230=etafTelmi230_h(T); %Pas %T in K            
            etaH2O=0.1/(2.20065*(T-282.92341+(8761.27+(T-282.92341)^2)^.5)-129.908); %Pas T in K
            cT230=230./obj.MtelH;
            obj.etaf_method='etaf_fun';
            obj.etaf_fun_h=@(c) 10.^((log10(etafTelmi230/etaH2O)/cT230).*c(1)+log10(etaH2O));
            
            %% Diffusivity
            Vmol_TH=33*14.8+30*3.7+2*15.6+2*12+2*12-4*15-2*11.5;
            obj.Vmol_solute=Vmol_TH;
            obj.assoc_par=2.6; %water
            obj.Mf=obj.MH2O;
            obj.diff_method='diff_wilkechang';
            
            %% Reaction
            obj.solver_name='dissolve_qamar_fxfvm';
            obj.model_reaction_h=@(t,c) obj.model_reaction(t,c);
            obj.ode_use_reaction=true;
            %% Equillibrium constants
            if useactivity
                IS = 0.5*(cNa2CO3*2 + cNa2CO3*4);
            else
                IS=0;
            end
            gH = obj.gammaX(0.9, 1,-1, IS);
            gOH = obj.gammaX(0.35, 1,-1, IS);
            gCO3 = obj.gammaX(0.45, 1,-2, IS);
            gHCO3 = obj.gammaX(0.45, 1,-1, IS);
            gT = obj.gammaX(0.4, 1,-1, IS);
            %gNa = obj.gammaX(0.45, 1,-2, IS);
            

            %%%%%%%%%%%%%%%%%%%%%%
            %pKaCO2 = -2622.38/T - 0.0178471*T + 15.5873; %(Harned and Davis, 1943)
            %pKaH2CO3 = 3404.71/T + 0.032786*T - 14.8435; %(Harned and Davis, 1943)
            %pKaHCO3 = 2902.39/T + 0.02379*T - 6.4980; %(Harned and Scholes, 1941)
            
            % Plummer and Busenberg, 1982, The solubilities of calcite, ...
            % https://doi.org/10.1016/0016-6817037(82)90056-4
            % Valid for T=0-100C, p=1atm
            pKaCO2=-108.3865-0.01985076*T+6919.53./T+40.45154*log10(T)-669365./T.^2;
            pKaH2CO3=356.3094+0.06091964*T-21834.37./T-126.8339*log10(T)+1684915./T.^2;
            pKaHCO3=107.8871+0.03252849*T-5151.79./T-38.92561.*log10(T)+563713.9./T.^2;
            
            pKw = -61.706 +5.864e3./T +22.645*log10(T); % Bandura, Lvov 2005, fit to data https://srd.nist.gov/JPCRD/jpcrd696.pdf,
            %pKw -log10(exp(148.9802 - 13847.26/T - 23.6521*log(T))); %(Dickson and Riley, 1979)
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.gH=gH;
            obj.Ka = obj.Kafun_T(T)/gH/gT;
            obj.KgCO2=10^-pKaCO2;
            obj.KaH2CO3= (10^-pKaH2CO3)/gHCO3/gH;
            obj.KaHCO3 = (10^-pKaHCO3)*gHCO3/(gH*gCO3);
            obj.Kw = (10^-pKw)/(gH*gOH);
        end % end of set_props
        
        function pH=cH2pH(obj,cH)
            % converts cH to pH, consider activity
            pH=-log10(cH.*obj.gH);
        end
        
        function cH=pH2cH(obj,pH)
            % converts pH to cH, consider activity
            cH=(10.^-pH)./obj.gH;
        end
        
        function cCO3=get_cCO3(obj,cH,cCO3tot)
            cCO3=cCO3tot.*obj.KaH2CO3.*obj.KaHCO3./(cH.^2+cH.*obj.KaH2CO3+obj.KaH2CO3.*obj.KaHCO3);
        end
        
        function cHCO3=get_cHCO3(obj,cH,cCO3tot)
            cHCO3=cCO3tot.*cH.*obj.KaH2CO3./(cH.^2+cH.*obj.KaH2CO3+obj.KaH2CO3.*obj.KaHCO3);
        end
        
        function cH2CO3=get_cH2CO3(obj,cH,cCO3tot)
            cH2CO3=cCO3tot.*cH.^2./(cH.^2+cH.*obj.KaH2CO3+obj.KaH2CO3.*obj.KaHCO3);
        end
        
        function psatCO2=get_psatCO2(obj,cH,cCO3tot)
            cH2CO3=get_cH2CO3(obj,cH,cCO3tot);
            psatCO2=cH2CO3./obj.KgCO2*obj.p0atm; %Pa
        end
        
        function cOH=get_cOH(obj,cH)
            cOH=obj.Kw./cH;
        end
        
        function charge=charge_balance(obj,cH,cTtot,cCO3tot,cNa)
           % charge balance equation in solution                
           cCO3=get_cCO3(obj,cH,cCO3tot);
           cHCO3=get_cHCO3(obj,cH,cCO3tot);
           cH2CO3=get_cH2CO3(obj,cH,cCO3tot);
           cOH=get_cOH(obj,cH);
           % total charge of Telmisartan
           cTcharge=sum(cTtot*obj.compfun_KacH(obj.Ka,reshape(cH,length(cH),1)).*obj.ind.charge,2);
           charge=cTcharge - cHCO3 - 2*cCO3 - cOH + cH + cNa;
        end
        
        function cH=solve_charge_balance(obj,cTtot,cCO3tot,cNa)
            % solve_charge_balance: (obj,cTtot,cCO3tot,cNa)
            cH=obj.fzerobnd(@(cH) obj.charge_balance(cH,cTtot,cCO3tot,cNa),1e-20,1,1e-15);
        end
        
        function [gX] = gammaX(obj,aX, zX1, zX2, IS)
            %gammaX(obj,aX, zX1, zX2, IX) activity of ion based on Davies eq
            %   %Debey Huckel, Davies eq. IS<0.5
            A=0.51;
            %IS=min([IS,0.5]); % IS<0.5
            %log10g = -A*abs(zX1*zX2)*(sqrt(IS)/(1 + 3.3*aX*sqrt(IS)) - 0.3*IS);
            log10g = -A*abs(zX1*zX2)*(sqrt(IS)/(1 + sqrt(IS)) - 0.3*IS);
            gX = 10^log10g;
        end
        
        function ceqred=get_ceqred(obj,cH)
            %Calculates reduced solubility c/c_intrinsic at given cH.
            selneutral=zeros(1,length(obj.ind.charge));
            selneutral(obj.ind.ineutral)=1;
            ceqred=1./sum(obj.compfun_KacH(obj.Ka,reshape(cH,length(cH),1)).*selneutral,2);
        end
        
        function composition=get_composition(obj,cH,cTtot,cCO3tot)
           % Returns current ionic composition of the solution
           cCO3=get_cCO3(obj,cH,cCO3tot);
           cHCO3=get_cHCO3(obj,cH,cCO3tot);
           cH2CO3=get_cH2CO3(obj,cH,cCO3tot);
           cOH=get_cOH(obj,cH);
           xTneutral=1./get_ceqred(obj,cH);
           cTneutral=cTtot.*xTneutral;           
           cTother=cTtot.*(1-xTneutral);
           composition=[cTneutral,cTother,cOH,cCO3,cHCO3,cH2CO3];
        end
        
        function [csat_fun,cHsat_fun] = get_csat_interp(obj)
        %   Create interpolant for prediction of saturated concentration based on
        %   current base concentration
        %       - Used to model impact of CO2 desorption from the solution
        %       - Valid only for defined initial concentration of base 
        %         (all Na+ ionts are originally from the base)
        %   Change in viscosity and diffusion coefficient seems that has no effect 
        %   on resulting saturated concentration, now is neglected
        cB0max=obj.cCO3tot;
        wB0=linspace(0.3,1,50);
        cB0arr=cB0max.*wB0;
        csatarr=zeros(size(cB0arr));
        cHsatarr=zeros(size(cB0arr));
            for i=1:length(cB0arr)
                [csatarr(i),cHsatarr(i)] = calc_csat(obj,obj.cTintrinsic,cB0arr(i),obj.cNa);
            end
            csat_fun=@(c) interp1(cB0arr,csatarr,c(2),'linear','extrap');
            cHsat_fun=@(c) interp1(cB0arr,cHsatarr,c(2),'linear','extrap');
        end
        
        function dcdt=model_reaction(obj,t,c)   
            % Model describing escape of CO2 to the atmosphere
            cTtot=c(1);
            cCO3tot=c(2);
            cH=obj.solve_charge_balance(cTtot,cCO3tot,obj.cNa);
            psatCO2=obj.get_psatCO2(cH,cCO3tot);           

            if psatCO2>obj.pCO2atm
                dcCO2totdt=-obj.kMaCO2.*(psatCO2-obj.pCO2atm)./obj.p0atm;
            else
                dcCO2totdt=0;
            end      
            dcdt=[0;dcCO2totdt];
            if abs(cTtot-obj.csat_fun_h(c))<1e-6
                dcdt=zeros(size(dcdt));
            end
        end
        
    end % end of methods
end % end of class PBMTelmi