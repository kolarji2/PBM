%
%   Plot initial PSD
%

%% Load experimental data
readdata=true;
load_dataset_paliperidone
data=data_nd(1:5);

%% Init model with just initial PSD
Nexp=length(data);
model_arr=cell(length(data),1);

for iexp=1:Nexp
    %% Create new model instance
    m=PBM;
    m.title=data{iexp}.expID;

    %% Init model
    msusp=[0.156,0.125,0.1,0.1];
    w=data{iexp}.w;
    V0=250*1e-3; %l
    cs0=sum(msusp.*w)/V0; %[g/l] = [kg/m3]
    w=msusp.*w./sum(msusp.*w); % correct w
    m.make_loggrid(-8,-3,200);
    m.dist_est_method='3lognormal';
    m.set_init_distFromExpData(data{iexp}.all.dexp,data{iexp}.phi0exp,w,[],false);
    %% Save
    model_arr{iexp}=m;
end


%% Visualize
fig=figure('Position',[100,100,1000,800]);
ax=gca(fig);


for iexp=1:Nexp    
    m=model_arr{iexp};
    % m.phi0_distest - vol. distribution [m^{-1}]    
    phim=m.phi0_distest.*m.dx.*100; % [vol. %]
    %sum(phim)==100
    semilogx(ax,m.x0,phim,'DisplayName',m.title)
    hold(ax,'on');
end
legend(ax)

return

%% Example of coarser discretization
m=model_arr{1};
xb_coarse=logspace(-8,-3,50)';
dx_coarse=xb_coarse(2:end)-xb_coarse(1:end-1);
[phim_coarse,x0_coarse]=get_mean_integral_values(m,m.x0,m.phi0_distest,xb_coarse);
phim_coarse=phim_coarse./sum(phim_coarse.*dx_coarse);

x0_coarse_um=x0_coarse*1e6;
ax=gca(figure);
bar(ax,log10(x0_coarse_um),phim_coarse.*dx_coarse*100,0.7);
%set(ax, 'XScale', 'log')
xticks_=[-1:3];
xticks(ax,xticks_);
xticklabels(ax,arrayfun(@(a) sprintf('%.1f',10.^a),xticks_,'UniformOutput',false));

