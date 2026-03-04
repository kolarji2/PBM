%
%   Example of usage of PBM class for modelling of KCl dissolution
%   with elimination kinetics
%

%% Load experimental data
load_datasetKCl

%% Init model
infoi=infoKCl{1};

%% Create new model instance
m=KClPBM(infoKCl{1},dataKCl{1});
m.set_optimized_pars();
m.title='Elimination model';
        
%% Init grid & particle size bases
m.make_loggrid(-7,-1,400);
Ndist=40;
mu_arr=logspace(log10(100),log10(5000),Ndist)*1e-6;
sigma_arr=mu_arr*0.1;
m.set_init_dist(mu_arr,sigma_arr,ones(size(mu_arr)),[],'lognormal');
m.solver_name='reactive_dissolve_fxfvm';

%% Elimination
kelim=0.01;
m.model_reaction_h=@(t,c) -kelim*c;

%% Set initial conditions
w=rand(1,Ndist);
m.w=w./sum(w);
cs0=20;
m.set_initial_conditions(0,cs0,'gl');
%% Solve model
m.ode_flux_limiter='none'; % turn off flux limiter (faster solution)
time=linspace(0,1200,200);
m.solve(time); 

%% Visualize
figure('Position',[100,100,800,600])
subplot(1,2,1)
semilogx(m.x0*1e6,m.phi0_distest.*m.dx);
title('fV')
subplot(1,2,2)
plot(m.tm,m.cm);
title('Concentration')

