%
%   Normal dissolution of Paliperidone-Palmitate dissolution
%

%% Load experimental data
readdata=true;
load_dataset_paliperidone
data=data_nd(1:5);

%% Init model
Nexp=length(data);
model_arr=cell(length(data),1);

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
    kMsol=1.33e-3; % m/day
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
    time=linspace(0,tmax,100);
    [time,c,fn,phi_dist] = m.solve(time);    
    %% Save
    model_arr{iexp}=m;
end

%% Visualize
fig=figure;
ax=gca(fig);
cols=lines(Nexp);
for iexp=1:Nexp
    m=model_arr{iexp};
    plot(ax,m.tm,m.cm./m.cs0,'DisplayName',m.title,'Color',cols(iexp,:))
    hold(ax,'on')
    cm=interp1(m.tm,m.cm,data{iexp}.tscatter);
    plot(ax,data{iexp}.tscatter,data{iexp}.cscatterPP.*M_PP./m.cs0,'o','DisplayName',sprintf('exp - %s',m.title),'Color',cols(iexp,:))
end
xlabel(ax,'Time [day]')
ylabel(ax,'Released [%]')
legend(ax,'Location','best')
