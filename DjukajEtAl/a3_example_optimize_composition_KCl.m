%
%   Example of optimization of composition
%   Inverse problem solution
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
    m.set_optimized_pars();
    m.title=names{iexp};    
    %% Save
    model_arr{iexp}=m;
end

%% Optimization to find optimal composition for Trimodal experiment
opt=OptimizePBM;
opt.models=model_arr(3);
opt.expdata=dataKCl(3);
opt.solver_method='fminsearch';
opt.x0=[0.7,0.2,0.2];
opt.target_list={'cm','cexp'};
opt.handle_list={'w'};
opt.handle_gen_list={@(x) x./sum(x)};
opt.log_to_file=false;

xres=opt.optimize('test (composition)','suada/opt_results.log')
xres=xres./sum(xres);
fprintf('found opt w: %f %f %f\n',xres)
fprintf('experimental w: %f %f %f\n',info{iexp}.w)