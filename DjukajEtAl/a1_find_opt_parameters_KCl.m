%
%   Example of optimization of model parameters for KCl dissolution
%   
%

%% Load experimental data
load_datasetKCl

%% Init model
Nexp=length(dataKCl);
model_arr=cell(length(dataKCl),1);
names={'monomodal','bimodal','trimodal'};

for iexp=1:Nexp
    %% Create new model instance
    m=KClPBM(infoKCl{iexp},dataKCl{iexp});
    % m.set_optimized_pars();
    m.ode_flux_limiter='Koren';
    m.solver_name='dissolve_qamar_fxfvm';
    m.ode_solver_h=@ode15s;
    m.title=names{iexp};    
    m.check_settings;
    model_arr{iexp}=m;
end

%% Optimization to find parameters (epsilon, a2) from first experiment
opt=OptimizePBM;
opt.models=model_arr(1);
opt.expdata=dataKCl(1);
opt.solver_method='fminsearch';
opt.lb=[0,0.0001];
opt.ub=[0.2,0.5];
opt.x0=[0.1,0.1];
opt.target_list={'cm','cexp'}; %, ,

opt.prop_list={'epsilon'};
opt.handle_list={'fAfun'};
opt.handle_gen_list={@(xopt) @(r,x) pi*(1.0+0.5*exp(-r/xopt(2)).*ones(size(x)))};
opt.time_max=615;
opt.objfun_generic([8.5e-02, 2.71e-01])

xres=opt.optimize('test (epsilon, a2)','DjukajEtAl/opt_results.log');
fprintf('epsilon: %.2e a2: %.2e\n',xres)
%% Test on all 3 experiments
opt.models=model_arr;
opt.expdata=dataKCl;
opt.objfun_generic(xres)
opt.objfun_generic([8.5e-02, 2.71e-01])