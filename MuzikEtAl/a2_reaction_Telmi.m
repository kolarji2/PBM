%
%   Example of usage of modified PBM class for 
%		modelling of pH dependent dissolution of telmisartan.
%		
%		Modified @PBM class can be found in directory: @PBMTelmi 
%

m=PBMTelmi();
loadData=0;
if ~exist('TimeEnd45','var')
    [exps,TimeEnd45,TimeEnd55,TimeEnd65,TimeEndData,prop] = load_exp_data_telmisartan(m);
    [phi0expdata] = load_psd_telmisartan();
end

wpure_Na2CO3=0.928;
wpure_Telmi=0.996;

Nexp=length(exps);
model_arr=cell(Nexp,1);

fname='MuzikEtAl/results/reaction_results.xlsx';

exps_sel=[1:6]
for iexp=exps_sel
    tic
        
    cNa2CO3=exps{iexp}.cB0*wpure_Na2CO3;
    cS0=exps{iexp}.cS0*wpure_Telmi;
    
    % PSD (i=1) - sonicated, true particle size
    % PSD (i=2) - low concentration, almost no agglomeration
    % PSD (i=3) - higher concentrations, mild agglomeration

    if cS0<0.01 % almost no agglomeration, smaller PSD for low conc
        w=[0,1,0]; % for exps 2,3
    else % mild agglomeration, bigger PSD for concentrated solutions
        w=[0,0,1];        
    end
    S0funT=@(T) exp(13.70760 + (-8707.17148./T));
    cTintrinsic=S0funT(exps{iexp}.T);
    
    m=PBMTelmi();
    m.init_telmi_comp_solver('MuzikEtAl/expdata/pka_dist.xlsx');
    useactivity=false;
    m.set_props(exps{iexp}.T,cNa2CO3,useactivity,prop.etafTelmi230);
    m.title=sprintf('cB = %.2f g/l cS= %.2f g/l',cNa2CO3*m.MNa2CO3,cS0*m.MtelH);
    sheet=sprintf('e%d_T=%.0fC_cB=%.2f_cS=%.2f',iexp,exps{iexp}.T-273.15,cNa2CO3*m.MNa2CO3,cS0*m.MtelH);
    
    m.solver_name='dissolve_qamar_fxfvm';
    m.ode_solver_h=@ode23t;
    m.ode_use_reaction=true;
    
    %% csat            
    m.cTintrinsic=cTintrinsic;
    m.csat_method='csat_fun';
    m.csat_fun_h=m.get_csat_interp();
    m.csat=calc_csat(m,m.cTintrinsic,m.cCO3tot,m.cNa);

    fprintf('(Exp %d) %s\n', iexp, m.title);
    fprintf('logS0(%.0f) = %.2f mol/l\n', m.T-273.15, log10(cTintrinsic));
    %% Optimized model pars
    m.epsilon=0.01; % m2/s3
    m.ode_scalepar=1e6;
    m.ode_flux_limiter='Koren';        % Koren
    
    m.kMaCO2=5e-4;
    if iexp==5
        m.kMaCO2=2.5e-4;        
    elseif iexp==4
        m.kMaCO2=8e-4;        
    end
    
    %% Init model
    m.make_loggrid(-8,-2,300);    
    phi0exp_data=cellfun(@(x) x.phi0exp,phi0expdata,'UniformOutput',false);
    dexp_data=cellfun(@(x) x.dexp,phi0expdata,'UniformOutput',false);    
    set_init_distFromExpData(m,dexp_data,phi0exp_data,w,[],false);        
    m.Nconc=2;            
    cstart=[0;cNa2CO3];
    m.set_initial_conditions(cstart,cS0,'moll');
    
    %% Solve model
    tmax=max(exps{iexp}.time)*60;
    time=[0, logspace(-1,log10(tmax),300)];
    m.ode_options=odeset('RelTol',1e-4,'AbsTol',1e-6);
    [time,c,fn,phi_dist] = m.solve(time);    
    %% Save
    model_arr{iexp}=m;
    toc
    
    %% Save to excel
    pHest=zeros(length(m.cm),1);
    for i=1:length(pHest)
       pHest(i)=m.cH2pH(solve_charge_balance(m,m.cm(i,1),m.cm(i,2),m.cNa));
       %pHest(i)=cH2pH(cHfun_cT(cT(i)));        
    end
    mdata=struct;
    mdata.tmodel_min=m.tm./60;
    mdata.cTmodel_gl=m.cm(:,1)*m.MtelH;
    mdata.cCO3model_gl=m.cm(:,2)*m.MNa2CO3;    
    mdata.pHmodel=pHest;
    edata=struct;
    edata.texp_min=exps{iexp}.time;
    edata.cTexp_gl=exps{iexp}.cT.*m.MtelH;
    edata.cTexp_std_gl=exps{iexp}.cTstd.*m.MtelH;
    edata.time_pH=exps{iexp}.time_pH;
    edata.pHexp=exps{iexp}.pH;
    edata.cS0_gl=[cS0.*m.MtelH];
    edata.cN2CO3_gl=[cNa2CO3.*m.MNa2CO3];
    fields={{'tmodel_min','cTmodel_gl','cCO3model_gl','pHmodel'},...
     {'texp_min','cTexp_gl','cTexp_std_gl','time_pH','pHexp','cS0_gl','cN2CO3_gl'}};
    m.save_to_excel({mdata,edata},fields,fname,sheet);    
end

%% Visualize
fig=figure('Position',[100,100,1200,800]);
cols=lines(2*Nexp);
for iexp=exps_sel
    ax=subplot(2,3,iexp);
    %fig=figure;ax=gca(fig);
    m=model_arr{iexp};
    yyaxis(ax,'left')
    if max(m.tm)/60>80
       tohours=60;
    else
        tohours=1;
    end
    plot(ax,m.tm./60./tohours,m.cm(:,1)*m.MtelH,'DisplayName','model','Color',cols(iexp,:))
    hold(ax,'on')
    title(ax,m.title)
    errorbar(ax,exps{iexp}.time./tohours,exps{iexp}.cT*m.MtelH,exps{iexp}.cTstd*m.MtelH,'o','DisplayName','exp','Color',cols(iexp,:))
    ylabel('c_{Telmi} [g/l]')
    legend(ax,'Location','south') 
    yyaxis(ax,'right')
    
    pHest=zeros(length(m.cm),1);
    for i=1:length(pHest)
       pHest(i)=m.cH2pH(solve_charge_balance(m,m.cm(i,1),m.cm(i,2),m.cNa));
       %pHest(i)=cH2pH(cHfun_cT(cT(i)));        
    end
    plot(ax,m.tm./60./tohours,pHest,'DisplayName','pH','Color',cols(iexp+1,:))    
    plot(ax,exps{iexp}.time_pH./tohours,exps{iexp}.pH,'DisplayName','pH - exp','Color',cols(iexp+1,:))
    ylim(ax,[8,12])
    if max(m.tm)/60>80
        xlabel('Time [h]')
    else
        xlabel('Time [min]')
    end
    
    %m.savefig_article(fig,sprintf('MuzikEtAl/results/reaction%d',iexp))
    %plot(ax,m.time./60,m.c(:,2),'DisplayName','CO3','Color',cols(iexp+1,:))
end
m.savefig_article(fig,'MuzikEtAl/results/reaction_variableKMaCO2')