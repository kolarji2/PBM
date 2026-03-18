%
%   Find optimal saturated concentration for reactive dissolution
%

%% Load experimental data
readdata=true;
load_dataset_paliperidone
data=data_rd;

%% Init model
Nexp=length(data);
model_arr=cell(length(data),1);
data_arr=cell(length(data),1);

%% Molar mass
M_P  = 426.5e-6; %g/umol Paliperidone
M_PP = 664.9e-6; %g/umol Paliperidone Palmitate
M_SDS = 288.372e-6; %g/umol
M_LiOH = 23.95e-6; % g/umol   

for iexp=1:Nexp
    %% Create new model instance
     T=70+273.15;
    m=create_new_model_PP(T,M_PP);
    m.title=data{iexp}.expID;
    m.ode_options=odeset('RelTol',1e-6,'AbsTol',1e-8,'NonNegative',1);
    
    %% Optimized model pars   
    % estimated csat at basic pH
    m.csat =0.525; % [g/l] = [kg/m3] %1540
    kMsol=1.33e-3; % m/day
    m.kM_method='kM_fun';
    m.kM_fun_h=@(kM0) kMsol;  
    
    %% Reaction submodel
    pars={}; 
    pars.k1=6.873; % m/day
    pars.kdeg=0.02; % m/day
    pars.M_P  = M_P; %g/umol Paliperidone
    pars.M_PP = M_PP; %g/umol Paliperidone Palmitate
    pars.M_LiOH = M_LiOH;
    m.ode_use_reaction=true;
    m.model_reaction_h=@(t,c) model_pp_reaction(t,c,pars);
    
    %% Init model
    msusp=[0.156,0.125,0.1,0.1];
    w=data{iexp}.w;
    V0=250*1e-3; %l
    cs0=sum(msusp.*w)/V0; %[g/l] = [kg/m3]    
    w=msusp.*w./sum(msusp.*w); % correct w
    m.make_loggrid(-8,-3,200);
    m.set_init_distFromExpData(data{iexp}.all.dexp,data{iexp}.phi0exp,w,[],false);        
    m.Nconc=3;        
    cB0=1.6701e+04.*M_LiOH; %  [g/l]  % LiOH concentration  pH ~ 11 
    cstart=[0;0;cB0]; % initiate concentration [cPP,cP,cLiOH];
    m.set_initial_conditions(cstart,cs0,'gl');
    
    %% Solve model
    tmax=max(data{iexp}.texp);
    time=linspace(0,tmax,300);
    %[time,c,fn,phi_dist] = m.solve(time);    
    %% Save
    model_arr{iexp}=m;    
    data_arr{iexp}=struct;
    data_arr{iexp}.texp=data{iexp}.texp;
    data_arr{iexp}.cscatterPP=data{iexp}.cscatterPP.*M_PP;
    data_arr{iexp}.cscatterP=data{iexp}.cscatterP.*M_P;
    
end

%% Optimization to find opt csat for experiment with reaction
opt=OptimizePBM;
opt.models=model_arr(1:5);
opt.expdata=data_arr(1:5);
opt.solver_method='fminsearch';

%csat
opt.x0=[0.5];
opt.target_list={'cm',[1],'cscatterPP',[1]; ...
                 'cm',[2],'cscatterP',[1]}; %,
opt.prop_list={'csat'};
opt.log_to_file=false;
xres=opt.optimize([],[]);
fprintf('optimal csat: %.5f g/l\n',xres); % 0.564 g/l
%% Test on all experiments
opt.models=model_arr;
opt.expdata=data_arr;
opt.objfun_generic(xres)


