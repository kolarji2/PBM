%
%   Ideal particle size for steady state dissolution
%

%% Init model

M_P  = 426.5e-6; %g/umol Paliperidone
M_PP = 664.9e-6; %g/umol Paliperidone Palmitate
M_SDS = 288.372e-6; %g/umol
M_LiOH = 23.95e-6; % g/umol



mu_arr0=[1,5,10,100,100].*1e-6;
cs0_arr=[0.4752,2.3762,4.7525];
kmultiplier=[1,0.1,0.02];
tmax_arr=[12,30,100];
area_init=600;
Nmu=length(mu_arr0);
Nk=length(kmultiplier);
Nc=length(cs0_arr);
model_arr=cell(Nmu,Nc,Nk);

dts_1_100um=load('SonntagEtAl/expdata/opt1_100um_fine.mat');
%dts_1_100um=load('SonntagEtAl/expdata/opt1_100um_e100.mat');
w100um=dts_1_100um.w100um;
tstable_arr=dts_1_100um.tstable_arr;
    
for imu=1:Nmu    
    imu
    for ik=1:Nk
        parfor ic=1:Nc
        
        %% Create new model instance
        T=70+273.15;
        m=create_new_model_PP(T,M_PP);
        m.title=sprintf('%.0f um',mu_arr0(imu)*1e6);
        m.ode_options=odeset('RelTol',1e-6,'AbsTol',1e-8,'NonNegative',1);
        m.ode_flux_limiter='none'
        %% Optimized model pars
        m.csat =0.525; % [g/l] = [kg/m3]
        kMsol=1.33e-3; % m/day 
        m.kM_method='kM_fun';
        m.kM_fun_h=@(kM0) kMsol;
        
        %% Reaction submodel
        pars={};
        pars.k1=6.873*kmultiplier(ik); % m/day
        pars.kdeg=0.02*kmultiplier(ik); % m/day
        pars.M_P  = M_P; %g/umol Paliperidone
        pars.M_PP = M_PP; %g/umol Paliperidone Palmitate
        pars.M_LiOH = M_LiOH;
        m.ode_use_reaction=true;
        m.model_reaction_h=@(t,c) model_pp_reaction(t,c,pars);
        
        %% Init model
        
        V0=250*1e-3; %l
        
        m.make_loggrid(-8,-3,300);    
        if imu==length(mu_arr0)
            mu_arr=[1e-6 mu_arr0(imu)];
            w=[1-w100um(ic,ik),w100um(ic,ik)];
        else
            mu_arr=[mu_arr0(imu)];        
            w=[1];
        end
        sigma_arr=mu_arr.*0.1;
        m.set_init_dist(mu_arr,sigma_arr,w,[],'lognormal');
        area_perkg=sum(m.phi0_distest.*m.dx./(m.x0.^3.*m.fV).*(m.x0.^2.*m.fA))./m.rhos;
        msusp=area_init/area_perkg;    
        cs0=cs0_arr(ic); %msusp/V0; %[g/l] = [kg/m3]    
        
        m.Nconc=3;
        cB0=1e9;%  place holder
        cstart=[0;0;cB0]; % initiate concentration [cPP,cP,cLiOH];
        m.set_initial_conditions(cstart,cs0,'gl');
        
        %% Solve model
        tmax=tmax_arr(ik);
        time=linspace(0,tmax,400);
        [time,c,fn,phi_dist] = m.solve(time);    
        %% Save
        model_arr{imu,ic,ik}=m;    
        end
    end
end

%% Visualize
utils=Utils();
fig=figure('Position',[100 100 1200 600]);
ax_all={};
for ik=1:Nk
    ax_all{ik}=subplot(1,3,ik,'Position',[(0.07+(ik-1)*0.33) 0.35 0.25 0.55]);
end


for ik=1:Nk
ax=ax_all{ik};
cfg={'LineWidth',2};
fs=18;

%ax2=subplot(1,2,2);
Ncols=round(Nmu*1.5);
cols_all={winter(Ncols),autumn(Ncols),summer(Ncols)};
conc_umol_l=[M_PP,M_P,M_LiOH];
%cols=lines(Nexp);    
ls={'--',':','-.','--','-'};
larr=[];
for ic=1:Nc
    cols=cols_all{ic};
    for imu=1:Nmu        
    m=model_arr{imu,ic,ik};    
    if imu==length(mu_arr0)
        m.title=sprintf('1 + %.0f $\\mu$m (%.1f g/l)',mu_arr0(imu)*1e6,m.cs0);
    else
        m.title=sprintf('%.0f $\\mu$m (%.1f g/l)',mu_arr0(imu)*1e6,m.cs0);
    end
    csm=sum(m.fnm.*m.dx.*(m.x0.^3.*m.fV),1).*m.rhos;
    title0=m.title;
    
    l=plot(ax,m.tm,m.cm(:,1),ls{mod(imu-1,length(ls))+1},'DisplayName',title0,'Color',cols(imu,:),cfg{:});
    larr=[larr l];
    % Threshold
    cs_threshold=0.05;
    cs_target=max(m.cm(:,1))./(1+cs_threshold);
    cmax_target=cs_target*(1+cs_threshold);
    cmin_target=cs_target*(1-cs_threshold); 
    if imu==length(mu_arr0)
        plot(ax,[0 max(m.tm)],[cmax_target cmax_target],'--','Color',0.5*ones(1,3),cfg{:})
        lth=plot(ax,[0 max(m.tm)],[cmin_target cmin_target],'--','Color',0.5*ones(1,3),'DisplayName','$\pm 5$\% (1 + 100 $\mu$m)',cfg{:});
    end    
    %plot(ax2,m.tm,m.cm(:,2)./conc_umol_l(2),'DisplayName',title0,'Color',cols(iexp,:))
    hold(ax,'on')       
    %hold(ax2,'on')       
    end
end


xlim(ax,[0,max(m.tm)])
grid(ax,'on')
xlabel(ax,'Time [days]','Interpreter','latex')
ylabel(ax,'PP concentration [g/l]','Interpreter','latex')

if kmultiplier(ik)==1
    title(ax,'Hydrolysis rate 1x','Interpreter','latex')
else
    title(ax,sprintf('Hydrolysis slowed %dx',1/kmultiplier(ik)),'Interpreter','latex')
end
set(ax,'FontSize',fs,'TickLabelInterpreter','latex')

end
for i=1:4
    l0=plot(ax,0,0,'DisplayName','','Color',ones(1,3));
    larr=[larr l0];
end
larr=[larr lth];
legend(ax,larr,'NumColumns',5,'Orientation','horizontal', ...
    'FontSize',fs,'Interpreter','latex','Position',[0.1 0.03 0.8 0.15])

utils.savefig_article(fig,'SonntagEtAl/results/quasi_steady_state_article');
