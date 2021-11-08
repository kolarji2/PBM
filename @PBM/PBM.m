classdef PBM < Utils & matlab.mixin.Copyable & matlab.mixin.SetGet 
    %% Population Balance Model (PBM) for Dissolution
    %  
    %% About PBM
    % 
    % PBM is a collection of matlab code for solving the population balance 
    % equation (PBE) with size-dependent negative growth term, the case of
    % dissolution.
    % 
    % PBE is solved via finite volume method with flux limiter.
    % 
    % TODO: 
    %   dissolve_fxFVM
    %   dissolve_fxxxFVM
    %   dissolve_fnFD
    %   dissolve_Naive
    %
    % TODO write and check documentation:
    % Particle size distribution from experimental data
    % Dissolution rate description
    % Optimization, finding correct model parameters
    % Reaction kinetics
    %
    %%  Usage
    %
    % PBM class contains all functionality
    %
    %   >> m = PBM; % create a new instance of the model
    %
    %   View list of parameters that you have to set-up
    %   >> m.check_settings(); 
    %
    %   Some parameters must be setup manualy, others have initialization
    %   functions, check bellow
    %   >> m.... % 
    %   >> m.... % set physico-chemical properties
    %   >> m.... %
    %
    %   Creates discretization, values are log10 of the size, units in [m]
    %   assigns x0, xb, dx, dxii
    %   >> m.make_loggrid(-7,-2,400); 
    %
    %   create initial particle size distribution
    %   assigns: 
    %       phi0_distestAll - any number of normalized volume distributions
    %       w               - required weight/volume fraction how to
    %                         combine distributions from phi0_distestAll 
    %       phi0_distest    - initial normalized volume distribution of 
    %                         particle sizes.
    %   >> m.set_init_dist([5e-4],[1e-3],[1],[],'lognormal');
    %
    %   Setups the concentration and model initial conditions. Prepares
    %   initial number based distribution from phi0_distest and cs0.
    %   assigns: 
    %       conc_units, vol2c, c2vol, cs0, cstart, fnstart, fxstart, y0fx,
    %       y0fn
    %   >> m.set_initial_conditions(0,cs0,'gl');
    %
    %   % Solves the model
    %   >> time=[0:1:60];
    %   >> m.solve(time);
    %
    
    %% PBM Properties
        
    properties
        
        %% General info
        
        % Name/identificator
        title
        
        %% Default methods for model calculation
        
        % Name of the solver used for PBE solution        
        % Defines numerical method:
        % * dissolve_fxfvm - finite volume method on fn*x
        % * dissolve_fvfvm - finite volume method on fn*x^3        
        solver_name='dissolve_fxfvm';
                
        % Use flux limiter
        ode_flux_limiter=1;
        
        % ODE scaling parameter to improve performance
        ode_scalepar=1e3;
        
        % Options for matlab ode solver
        ode_options
        
        % Method for estimation of size distribution
        dist_est_method='3lognormal'; 
        
        % Calculate mass trasnfer and growth rate
        Gfun=@(obj,t,c,xb) calc_G(obj,t,c,xb);
        
        % Method for mass transfer estimation (used by calc_G)
        Sh_method='kolmogorov';
        
        % Method for estimating solute saturated concentration (used by calc_G)
        % Possible methods
        % * csat_const
        % * csat_fun
        csat_method='csat_const'
        
        % Method for density calculation (used by calc_G)
        % Possible methods
        % * rhof_const
        % * rhof_fun
        rhof_method='rhof_const';
        
        % Method for viscosity calculation (used by calc_G)
        % Possible methods
        % * etaf_const
        % * etaf_fun
        etaf_method='etaf_const';
        
        % Method for difusivity calculation (used by calc_G)
        % Possible methods
        % * diff_const
        % * diff_wilkechang
        diff_method='diff_const';
        
        % Number of compounds (used for reactive dissolution)
        Nconc=1;
        
        % Initial concentration of solid particles
        cs0
        
        % Initial concentration of dissolved particles (solute)
        cstart
        
        %% Calculation results
        
        % Time used by model
        % Column vector
        tm
        
        % Calculated concentration by model
        % Columns correspond to each concentration in cstart
        % Rows to timesteps
        cm
        
        % Calculated fn (number distribution function) for model
        % Columns correspond to timesteps
        % matrix(length(x0),length(tm))
        fnm
        
        % Calculated volume distribution by model
        % Columns correspond to timesteps
        % matrix(length(x0),length(tm))
        phim_dist
        
        %% FVM Discretization
        
        % number of cells
        N
        
        % cell centers
        x0
        
        % Boundaries positions
        xb
        
        % cell sizes
        dx
        
        % distance between cell centers
        dxii
        
        % Volume fraction distribution of initial particle sizes
        phi0_distest
        
        % For multimodal particle size distribution contains dist of each fraction        
        phi0_distestAll
        
        % Weight/volume fraction for combination of multimodal particle size dist
        w
        
        % Average mean values for each distribution
        mu_arr
        
        % Average std values for each distribution
        sigma_arr
               
        % Number based PB of particles [#particles/volume]
        fnstart
        
        % Number based PB of particles multiplied by x [#particles.m/volume]
        fxstart
        
        % fnfvm initial model condition
        % [cstart; reshape(fnstart,N,1)];
        y0fn
        
        % fxfvm initial model condition
        % [cstart; reshape(fxstart,N,1)];
        y0fx        
        
        %% Matter properties (density, viscosity, diffusivity ...)
        
        % Temperature [K]
        T
        
        % Fluid density [kg/m^3]
        rhof
        
        % Function handle to obtain fluid density [kg/m^3]
        % e.g. @(c) my_etaf(c,pars);
        rhof_fun_h
        
        % Solute solid density [kg/m^3]
        rhos
        
        % Fluid viscosity [Pa.s]
        etaf
        
        % Function handle to obtain fluid viscosity [Pa.s]
        % e.g. @(c) my_etaf(c,pars);
        etaf_fun_h;
        
        % Molar mass of fluid [g/mol]
        Mf
        
        % Molar mass of solute (solid) [g/mol]
        Ms  
        
        % Diffusion constant
        diff
        
        % Association parameter between liquid (solvent) and solute
        assoc_par
        
        % Molar volume of solute [cm3/mol]
        Vmol_solute
        
        % Solute saturated concentration [kg/m3]
        % Used if method is csat_const
        csat
        
        % Function handle to calculate solute saturated concentration [kg/m3]
        % Used if method is csat_fun
        % @(c) mycsat_fun(c,pars)
        csat_fun_h
        
        %% Units conversion
        
        % Conversion constant to convert volume to concentration
        vol2c
        
        % Conversion constant to convert concentration to volume
        c2vol
        
        % Units used for concentration
        % Possible values:
        % 'moll'
        % 'gl'
        % 'kgm3'
        conc_units
        
        %% Parameters for mass transfer estimation
        
        % Molar gass constant
        R=8.314;
        
        % Gravity accelaration constant
        g=9.81; %m2/s
        
        % Accelaration for settling velocity calculation
        acc
        
        % Mass transfer coefficient
        kM
        
        % Volume shape factor
        fV=pi/6;
        
        % Area shape factor
        fA=pi;
        
        % Modification of surface area
        % e.g. @(t,crel) obj.fA
        fAfun;
        
        % mixing energy dissipation per unit mass
        epsilon=0.01;
        
        %Velocity Correction factors, (1 = no correction)
        vcorr=1;
        
        %Settling velocity Correction factors, (1 = no correction)
        vtcorr=1;
        
        % kM Correction factors, (1 = no correction)
        kmcorr=1;
        
        % diffusion Correction factors (1 = no correction)
        diffcorr=1;
        
        % Additional scalar velocity (0 = off)
        v=0;
        
        % Impeler RPM (Sh est by archimedes)
        impeler_rpm
        
        % Impeler diameter [m] (Sh est by archimedes,levins)
        impeler_diam
        
        % Tank diameter [m] (Sh est by archimedes,levins)
        tank_diam
        
        %Additional correction to slip velocity (Sh est by levins)
        levinsSlipCorr= 0;
        
        % Ratio: impeler_diam/tank_diam (Sh est by levins) (default 1)
        DsDt=1;
        
        % Kolmogorov parameterA (Sh est by kolgomorov)
        kolmogorovA=0.52;
        
        % Kolmogorov parameterB (Sh est by kolgomorov)
        kolmogorovB=0.52;
        
        % Kolmogorov parameterC (Sh est by kolgomorov)
        kolmogorovC=1/3;
        
        %% Reactive dissolution
        
        % Function handle that calculates the reaction kinetics
        % @(t,c) mymodel_reaction(t,c,pars)
        model_reaction_h
    end
    
    properties (SetAccess = private)
        
        % solver_name options
        solver_name_list={'dissolve_fxfvm','dissolve_fvfvm','reactive_dissolve_fxfvm'};
        
        % properties needed for Wilke-Chang method
        check_diff_wilkechang_list={'etaf','T','assoc_par','Mf','Vmol_solute'};
        
        % properties needed for model solution
        check_model_list={'cs0','cstart','fxstart','fnstart','x0','xb','dx', ...
            'dxii','y0fn','y0fx','vol2c','c2vol'};
        
        % properties that must be set manually, compound specific
        % properties
        check_matter_list={'rhof','rhos','etaf','diff','csat','Ms'}
        
        % handles needed for specific options
        % triplets: option,value, handle_needed
        check_handles_list={'csat_method','csat_fun','csat_fun_h'; ...
         'etaf_method','etaf_fun','etaf_fun_h'; ...
         'rhof_method','rhof_fun','rhof_fun_h'; ...
         'solver_name','reactive_dissolve_fxfvm','model_reaction_h'};
            
    end
    
    %% PBM methods
    
    methods       
        
        function obj=PBM(obj)
            % Constructor of the PBM
            obj.acc=obj.g;            
            obj.fAfun=@(t,crel) obj.fA;
        end
        
        function make_loggrid(obj,log10xmin,log10xmax,Ncells)
            % make_loggrid (obj,log10xmin,log10xmax,Ncells): create log scale grid
            % Create log scale discretization grid and includes 0
            % log10xmin - log10 of minimum particle size [m] (do not set this to 0)
            % log10xmax - log10 of maximum particle size [m]
            % Ncells - number of cells in FVM method
            obj.N = Ncells;
            obj.xb=[0,logspace(log10xmin,log10xmax,Ncells)]';
            obj.x0=obj.get_xclass_fromBounds(obj.xb);
            obj.dx=obj.get_dx_fromBounds(obj.xb);
            obj.dxii=obj.x0(2:end)-obj.x0(1:end-1); 
        end
        
        function set_init_dist(obj,mu_arr,sigma_arr,w,sieve,pdf_type)
            % set_init_dist(obj,mu_arr,sigma_arr,sieve,pdf_type):
            % Assign initial volume distributions to phi0_distestAll
            % and phi0_distest as well.
            %
            % Inputs:
            %   mu_arr - mean particle size [m]
            %       array(1,Ndists)
            %   sigma_arr - std of particle size [m]
            %       array(1,Ndists)
            %   w - weight/volume ratio of fractions
            %       array(1,Ndists)
            %   sieve - array of min/max boundaries for each phi0 distribution
            %       array(Ndists,2) or [] to turn off
            %   pdf_type - type of probability density function
            %       "lognormal" or "normal"
            
            if length(mu_arr)~=length(sigma_arr)
                error('Size of mu_arr and sigma_arr is not same!')
            end           
            
            if contains(pdf_type,'lognormal')
                    dist0 = @obj.distfunLogNormal;
            else
                dist0 = @obj.distfunNormal;
            end
            
            Ndist=length(mu_arr);
            obj.phi0_distestAll=zeros(obj.N,Ndist);
            for i=1:Ndist
                if isempty(sieve)
                    sievei=[0,1e9]; % do no cut dist0
                else
                    sievei=sieve(i,:);
                end               
                obj.phi0_distestAll(:,i)=dist0(obj.x0,obj.dx,mu_arr(i), ...
                                          sigma_arr(i),sievei(1),sievei(2));
            end
            % Setting w, assigns distribution phi0_distest as well:
            obj.w=w;
            
        end
        
        function set_init_distFromExpData(obj,dexp_data,phi0exp_data,w,sieve,dispfit)
            %get_init_distFromExpData (obj,dexp_data,phi0exp_data,w,sieve,dispfit):
            % Create distribution estimation using lognormal or normal
            % distribution.
            %
            % Inputs:
            %   dexp_data - contains dexp for each phi0 distribution
            %       cell(Ndists,1)
            %   phi0exp_data
            %       cell(Ndists,1)
            %   w - weight/volume ratio of fractions
            %       array(1,Ndists)
            %   sieve - array of min/max boundaries for each phi0 distribution
            %       array(Ndists,2) or [] to turn off
            %   dispfit - if true, resulting fit will be displayed
            %           true/false
            % return:
            %   fnstart f(x) - number density function
            %   fxstart: f(x)*x

            if ~iscell(phi0exp_data)
               error('phi0exp must be cell containing phi for each fraction!')
            end
            phi0_distestAll=zeros(length(obj.x0),length(phi0exp_data));
            mu_arr=zeros(length(phi0exp_data),1);
            sigma_arr=zeros(length(phi0exp_data),1);
            for i=1:length(phi0exp_data)

                if isempty(sieve)
                    sievei=[0,1e9]; % do no cut phi0
                else
                    sievei=sieve(i,:);
                end    
                xbexp=dexp_data{i};   
                if length(xbexp)==length(phi0exp_data{i}) %prepend zero if needed
                    xexp=obj.get_xclass_include0(xbexp);  
                    dxexp=obj.get_dx_include0(xbexp);  
                else
                    xexp=obj.get_xclass_fromBounds(xbexp);
                    dxexp=obj.get_dx_fromBounds(xbexp);   
                end
                             
                
                [phi0_distest,dist0fun,mu1,sigma1]= obj.est_dist(obj.x0,obj.dx,xexp,dxexp,phi0exp_data{i},obj.dist_est_method,sievei,dispfit);    
                phi0_distestAll(:,i) =phi0_distest;
                mu_arr(i)=mu1;
                sigma_arr(i)=sigma1;
            end
            obj.phi0_distestAll=phi0_distestAll;            
            obj.mu_arr=mu_arr;
            obj.sigma_arr=sigma_arr;
            
            % Setting w, assigns distribution as well:
            obj.w=w;            
            % phi0_distest=sum(w.*phi0_distestAll,2);
            % obj.phi0_distest=phi0_distest./sum(phi0_distest.*obj.dx);
            
        end
        
        function set.w(obj,w)
           obj.w=w;
           if ~isempty(obj.phi0_distestAll)
               if size(w,2)~=size(obj.phi0_distestAll,2)                   
                warning('Size of w=(%d) and phi0_distestAll=(%d) not consistent. w must be a row vector.', ...
                    size(w,2),size(obj.phi0_distestAll,2))
               else
                phi0_distest=sum(w.*obj.phi0_distestAll,2);
                obj.phi0_distest=phi0_distest./sum(phi0_distest.*obj.dx);
               end
           else
               warning('Setting w, but no phi0_distestAll are set')
           end           
        end
                
        function set_initial_conditions(obj,cstart,cs0,conc_units)
            % set_initial_conditions(obj,cstart,cs0,conc_units):
            % Set initial model conditions:
            %   concetrations: cstart,cs0
            %
            %   conc_units: concentration units (must be same for both
            %               cstart and cs0)
            % Set fnstart, fxstart from phi0_distest
            obj.conc_units=conc_units;
            obj.vol2c=obj.units_conv(1,'m3m3',obj.conc_units,obj.Ms,obj.rhos);        
            obj.c2vol=obj.units_conv(1,obj.conc_units,'m3m3',obj.Ms,obj.rhos);
            
            obj.cs0=cs0;
            obj.cstart=cstart;
            
            obj.fnstart=obj.cs0.*obj.c2vol.*obj.phi0_distest./(obj.x0.^3*obj.fV);
            obj.fxstart=obj.fnstart.*obj.x0; %[particles/m3]
            obj.y0fx = [cstart; reshape(obj.fxstart,obj.N,1)];
            obj.y0fn = [cstart; reshape(obj.fnstart,obj.N,1)];
        end

        function diff=calc_diff_wilkechang(obj)
            % Calculate diffusion using Wilke-Chang method
            diff=obj.wilkechang(obj.etaf,obj.T,obj.assoc_par,obj.Mf,obj.Vmol_solute);
        end
        
    end
    
end