%
%   Finds optimal value of kM for normal dissolution
%

%% Load experimental data
readdata=true;
load_dataset_paliperidone
data=data_nd;

%% Init model
Nexp=length(data);
model_arr=cell(Nexp,1);
data_arr=cell(Nexp,1);

M_P  = 426.5e-6; %g/umol Paliperidone
M_PP = 664.9e-6; %g/umol Paliperidone Palmitate
M_SDS = 288.372e-6; %g/umol
M_LiOH = 23.95e-6; % g/umol

for iexp=1:Nexp
   %% Create new model instance
    T=70+273.15;
    m=create_new_model_PP(T,M_PP);
    m.title=data{iexp}.expID;
    
    %% Optimized model pars
    m.csat =1540.*M_PP; % [g/l] = [kg/m3] %1540
    kMsol=1.33e-3;
    m.kM_method='kM_fun';
    m.kM_fun_h=@(kM0) kMsol;
    
    %% Init model
    msusp=[0.156,0.125,0.1,0.1];
    V0=250*1e-3; %l
    w=data{iexp}.w;
    cs0=sum(msusp.*w)/V0; %[g/l] = [kg/m3]
    w=msusp.*w./sum(msusp.*w); % correct w    
    m.make_loggrid(-8,-3,200);
    m.dist_est_method='3lognormal';
    m.set_init_distFromExpData(data{iexp}.all.dexp,data{iexp}.phi0exp,w,[],false);
    m.set_initial_conditions(0,cs0,'gl');
    
    %% Solve model
    tmax=max(data{iexp}.texp);
    time=linspace(0,tmax,300);
    %[time,c,fn,phi_dist] = m.solve(time);    
    %% Save
    model_arr{iexp}=m;    
    data_arr{iexp}=struct;
    data_arr{iexp}.texp=data{iexp}.tscatter;
    data_arr{iexp}.cscatterPP=data{iexp}.cscatterPP.*M_PP;
end

%% Optimization to find parametersfrom first experiment
opt=OptimizePBM;
opt.models=model_arr(1:5);
opt.expdata=data_arr(1:5);
opt.solver_method='fminsearch';
%kmconst
opt.x0=[0.0013];
opt.TolX=1e-6;
opt.target_list={'cm',[1],'cscatterPP',[1]}; % 'cscatterPP' relative error
opt.handle_list={'kM_fun_h'};
opt.handle_gen_list={@(x) @(kM0) x(1)};
opt.log_to_file=false;
xres=opt.optimize([],[]); % optimize with no log (idname,logFileName)
fprintf('optimal kM: %.5e m/day\n',xres); % 0.00133
%% Test on all experiments
opt.models=model_arr;
opt.expdata=data_arr;
opt.objfun_generic(xres)


