# PBM
PBM is a collection of matlab code files for solving the population balance equation (PBE) with size-dependent negative growth term, the case of dissolution.

PBE is solved via finite volume method with flux limiter.
    
##  Usage

Create a new instance of the model

	m = PBM; 

    
View list of parameters that you have to set-up

	m.check_settings(); 


Some parameters must be setup manualy, others have initialization
functions, check bellow.

Set physico-chemical properties

	m.rhof=1000; % Density of fluid [kg/m3 ~ g/l]
	m.rhos=1200; % Density of solids [kg/m3 ~ g/l]
	m.etaf=0.001; % Viscosity of fluid [Pa.s]
	m.diff=1e-9; % Diffusivity of the dissolved solids [m^2/s]
	m.csat=20;% Saturate concentration [kg/m3 ~ g/l]
	m.Ms=58.5; % Molar mass of solid [g/mol]


Creates discretization, values are log10 of the size, units in [m]
assigns x0, xb, dx, dxii

	m.make_loggrid(-7,-2,400); 


Create initial particle size distribution

Assigns: 

- phi0_distestAll - any number of normalized volume distributions
- w               - required weight/volume fraction how to combine distributions from phi0_distestAll
- phi0_distest    - initial normalized volume distribution of particle sizes.

```
m.set_init_dist([5e-4],[1e-3],[1],[],'lognormal');
```

Setups the concentration and model initial conditions. Prepares
initial number based distribution from phi0_distest and cs0.

Assigns: conc_units, vol2c, c2vol, cs0, cstart, fnstart, fxstart, y0fx, y0fn

	cs0=6; % Initial concentration of solids [kg/m3 ~ g/l]
	m.set_initial_conditions(0,cs0,'gl');


Solves the model

	time=[0:1:60];
	m.solve(time);

Plot results
	
	figure
	plot(m.tm,m.cm)

    

