%
%   Example of usage of PBM class for modelling of KCl dissolution
%   with elimination kinetics
%

%% Load experimental data
load_datasetKCl
% Directories for results
results_dir='DjukajEtAl/results';

%% Control variales
calculate_elimination=false; % prepare data and save them, can take long time
plot_elimination=true;

%% Init model
infoi=infoKCl{1};

%% Create new model instance
m=KClPBM(infoKCl{1},dataKCl{1});
m.set_optimized_pars();
m.title='Elimination model';
        
%% Init grid & particle size bases
m.make_loggrid(-7,-1,400);
Ndist=50;
mu_arr=logspace(log10(10),log10(5000),Ndist)*1e-6;
sigma_arr=mu_arr*0.1;
m.set_init_dist(mu_arr,sigma_arr,ones(size(mu_arr)),[],'lognormal');
m.solver_name='reactive_dissolve_fxfvm';

%% Elimination
kelim=0.01;
cs_target=0.5;
cs_threshold=0.05;
m.model_reaction_h=@(t,c) -kelim*c;

%% Set initial conditions
w=rand(1,Ndist);
m.w=w./sum(w);
cs0=20;
m.set_initial_conditions(0,cs0,'gl');
m.ode_flux_limiter='none'; % turn off flux limiter (faster solution)

%% Calculate elimination
if calculate_elimination
    csi_arr=zeros(Ndist-1,1);
    time_max=zeros(Ndist-1,1);
    m_res=cell(Ndist-1,1);    
    time=linspace(0,2500,2500);

    for i=2:Ndist
        m0=m.copy();
        fun=@(csi) objfun_elimkin_simple(csi,i,time,m0,cs_target,cs_threshold);
        cres=m0.fzerobnd(fun,0,50,1e-6);
        csi_arr(i-1)=cres;
        [~,time_max(i-1),m_res{i-1}]=fun(cres);
    end

    save('DjukajEtAl/elimination_results.mat','csi_arr','time_max','m_res','time');
else
    dts=load('DjukajEtAl/elimination_results.mat');
    csi_arr=dts.csi_arr;
    time_max=dts.time_max;
    m_res=dts.m_res;
    time=dts.time;
end

if plot_elimination
    davg_arr=zeros(Ndist-1,1);
    for i=1:Ndist-1
        m=m_res{i};
        davg_arr(i)=sum(m.phi0_distestAll(:,i+1).*m.dx.*m.x0);
    end
    init_colors_articleDjukajEtAl;
    plotcfg={'LineWidth',2.5};
    fs={'FontSize',28};
    fstics={'FontSize',24};

    % Stable release
    fig=figure('Position',[100 100 800 600]);
    ax=gca(fig);
    colororder(ax,{blue,green})
    yyaxis(ax,'left');
    semilogx(ax,davg_arr*1e6,time_max,'-','color',blue,'DisplayName','Release time',plotcfg{:})
    xlabel(ax,'Initial particle diameter [\mum]',fs{:})
    ylabel(ax,'Stable release duration [s]',fs{:})
    yyaxis(ax,'right');
    semilogx(ax,davg_arr*1e6,csi_arr,'-','color',green,'DisplayName','Required {\it C}_{s,0}',plotcfg{:})
    ylabel(ax,'{\it C}_{s,0} [kg m^{-3}]',fs{:})
    legend(ax,'Location','northwest',fs{:})
    xlim(ax,[10,5000])
    xticks_=[1 10 100 1000 10000];
    xticks(ax,xticks_)
    xticklabels(ax,arrayfun(@(a) sprintf('%.0f',a),xticks_,'UniformOutput',false))    
    set(ax,'LineWidth',1.5,fstics{:});
    grid(ax,'on')
    Utils.savefig_article(fig,fullfile(results_dir,'elim_tmax'));
    
    % Example dissolution
    sel_dists=[29,44];
    cs_upper=cs_target.*(1+cs_threshold);
    cs_target=cs_target;
    cs_lower=cs_target.*(1-cs_threshold);

    fig=figure('Position',[100 100 800 600]);
    leg_={'Short release';'Long release'};
    ax=gca(fig);
    l0=cell(1,length(sel_dists));

    plot(ax,time,ones(size(time))*cs_lower,'--',plotcfg{:},'color',colgray(1,:));
    hold(ax,'on')
    l1=plot(ax,time,ones(size(time))*cs_target,'-.','color',yellow,'DisplayName','Target concentration',plotcfg{:});
    l2=plot(ax,time,ones(size(time))*cs_upper,'--',plotcfg{:},'color',colgray(1,:),'DisplayName','5% Threshold');

    h=patch(ax,[0,800,800,0],[cs_lower,cs_lower,cs_upper,cs_upper],0,'FaceColor',colgray(3,:),'EdgeColor','none');
    h.FaceAlpha=0.3;
    h.DisplayName='5% Threshold conc.';

    for i=1:length(sel_dists)
        idist=sel_dists(i);
        tmax=time_max(idist);
        m=m_res{idist};        
        l0{i}=plot(ax,m.tm,m.cm,'-','color',colm(i,:),plotcfg{:},'DisplayName',leg_{i});
        l3=plot(ax,[tmax tmax],[cs_upper,0],'--','color',green,plotcfg{:}, ...
                                    'DisplayName','End of stable release');
    end

    legend(ax,[l0{:},l3,h],'Location','south','FontSize',24);
    grid(ax,'on')
    xlim(ax,[-5,800])
    ylim(ax,[0,0.6])
    xlabel('Time [s]',fs{:});
    ylabel('{\it C}_b [kg m^{-3}]',fs{:});
    set(ax,'LineWidth',1.5,fstics{:});
    %grid(ax,'on')
    Utils.savefig_article(fig,fullfile(results_dir,'elim_cprofile'));
    
    % Example distribution
    fig=figure('Position',[100 100 800 600]);
    ax=gca(fig);
    l0=cell(1,length(sel_dists));
    for i=1:length(sel_dists)
        idist=sel_dists(i);
        m=m_res{idist};
        %leg_=sprintf('%.0f \\mum',davg_arr(idist)*1e6); 
        phim=m.phim_dist(:,1).*m.dx;        
        semilogx(m.x0*1e6,phim*100,'-','color',colm(i,:),plotcfg{:},'DisplayName',leg_{i});
        hold(ax,'on')
    end

    legend(ax,'Location','northwest',fs{:});
    xlim(ax,[1,5000]);
    xticks_=[1 10 100 1000];
    xticks(ax,xticks_)
    xticklabels(ax,arrayfun(@(a) sprintf('%.0f',a),xticks_,'UniformOutput',false))    
    xlabel('Diameter [\mum]',fs{:});
    ylabel('dQ_3 [%]',fs{:});
    set(ax,'LineWidth',1.5,fstics{:});
    grid(ax,'on')
    Utils.savefig_article(fig,fullfile(results_dir,'elim_dist'));
end

function [err,tmax,m]=objfun_elimkin_simple(csi,i,time,m,cs_target,cs_threshold)
    Ndist=length(m.w);
    cmax_target=cs_target*(1+cs_threshold);   
    cmin_target=cs_target*(1-cs_threshold);   
    mfraci=zeros(1,Ndist);
    mfraci(1)=cmin_target;
    mfraci(i)=csi;
    cs0=sum(mfraci);
    w=mfraci./cs0;    
    m.w=w;
    m.set_initial_conditions(m.cstart,cs0,m.conc_units);
    [tt,c] = m.solve(time);
    [cmax,cmaxi]=max(c);
    if cmax>cmax_target
        penalty=(cmax-cmax_target)*10;
    else
        penalty=0;
    end
    err=cmax-cmax_target;    
    sel=tt>tt(cmaxi);
    tmax=interp1(c(sel),tt(sel),cmin_target);    
end

