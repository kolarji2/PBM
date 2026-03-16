%
%   Example of usage of modified PBM class for 
%		modelling of pH dependent dissolution of telmisartan.
%		
%		Modified class can be found in directory: @PBMTelmi 
%   
% Prediction of reaction end for 4 theoretical scenarios

m=PBMTelmi();
loadData=0;
if ~exist('TimeEnd45','var')
    [exps,TimeEnd45,TimeEnd55,TimeEnd65,TimeEndData,prop] = load_exp_data_telmisartan(m);
    [phi0expdata] = load_psd_telmisartan();
end

% Scenarios (initial conditions):
cexpTelmi=[50,100,50,100]; %g/l
excessB=[1.04, 1.04, 1.37, 1.37]; % molar excess
tend0=[6,16,1,6]; %h 

wpure_Na2CO3=0.928;
wpure_Telmi=0.996;

Nexp=length(exps);
model_arr=cell(Nexp,1);

exps_sel=[1:4]; % [1,2,3] % 1:6;
for iexp=exps_sel %[3:6] %[4,6] %[1:5] %[1:5] %1:5 %[1:4]
    iexp
    tic

    m=PBMTelmi();
        
    Texp=65+273.15;
    cS0=cexpTelmi(iexp)./m.MtelH;
    cNa2CO3=excessB(iexp).*cS0.*wpure_Telmi;
    cS0=cS0*wpure_Telmi;

    if cS0<0.01 % almost no agglomeration, smaller PSD for low conc
        w=[0,1,0]; % for exps 2,3
    else % mild agglomeration, bigger PSD for concentrated solutions
        w=[0,0,1];        
    end
    S0funT=@(T) exp(13.70760 + (-8707.17148./T));
    cTintrinsic=S0funT(Texp);
    
    m.init_telmi_comp_solver('MuzikEtAl/expdata/pka_dist.xlsx');
    useactivity=false;
    m.set_props(Texp,cNa2CO3,useactivity,prop.etafTelmi230);
    m.title=sprintf('cB = %.2f g/l cS= %.2f g/l',cNa2CO3*m.MNa2CO3,cS0*m.MtelH);
    sheet=sprintf('e%d_T=%.0fC_cB=%.2f_cS=%.2f',iexp,Texp-273.15,cNa2CO3*m.MNa2CO3,cS0*m.MtelH);
    
    m.solver_name='dissolve_qamar_fxfvm';
    m.ode_solver_h=@ode23t;
    m.ode_use_reaction=true;
    
    %% csat            
    m.cTintrinsic=cTintrinsic;
    m.csat_method='csat_fun';
    m.csat_fun_h=m.get_csat_interp();
    m.csat=calc_csat(m,m.cTintrinsic,m.cCO3tot,m.cNa);
    
    %% Optimized model pars
    m.epsilon=0.01;
    m.ode_scalepar=1e6;
    m.ode_flux_limiter='Koren';    % Koren | None      
    
    m.kMaCO2=5e-4;
    
    %% Init model
    m.make_loggrid(-8,-2,300);    
    phi0exp_data=cellfun(@(x) x.phi0exp,phi0expdata,'UniformOutput',false);
    dexp_data=cellfun(@(x) x.dexp,phi0expdata,'UniformOutput',false);    
    set_init_distFromExpData(m,dexp_data,phi0exp_data,w,[],false);        
    m.Nconc=2;            
    cstart=[0;cNa2CO3];
    m.set_initial_conditions(cstart,cS0,'moll');
    
    %% Solve model
    tmax=3600*tend0(iexp);
    time=[0, logspace(-1,log10(tmax),300)];
    m.ode_options=odeset('RelTol',1e-4,'AbsTol',1e-6);
    [time,c,fn,phi_dist] = m.solve(time);    
    %% Save
    model_arr{iexp}=m;
    toc  
end

fig=figure('Position',[100,100,1200,800]);
cols=lines(2*Nexp);
for iexp=exps_sel
    ax=subplot(2,3,iexp);
    m=model_arr{iexp};
    1
    if max(m.tm)/60>80
       tohours=60;
    else
        tohours=1;
    end
    sol_threshold = 0.975; % Threshold to determine end of reaction
    plot(ax,m.tm./60./tohours,m.cm(:,1)./m.cs0,'DisplayName',m.title,'Color',cols(iexp,:))
    hold(ax,'on')
    plot(ax,[0,max(m.tm)./60./tohours],[sol_threshold,sol_threshold],'--r','DisplayName','Threshold (end of reaction)')
    ylabel('c_{Telmi}/ c_{Telmi, max}  [-]')
    legend(ax,'Location','south')
    if max(m.tm)/60>80
        xlabel('Time [h]')
        time_un='h'
    else
        xlabel('Time [min]')
        time_un='min'
    end
    cr=m.cm(:,1)./m.cs0;
    th=m.tm./60./tohours;
    tend_est=interp1(cr(cr<0.9999),th(cr<0.9999),sol_threshold);
    plot(ax,tend_est,sol_threshold,'or', 'DisplayName','t_{end}')
    tend_exp = tend0(iexp);
    tend_str=sprintf('t_{end} = %.2f %s,  c_{sat} = %.2f g/l',tend_est(1),time_un,m.csat.*m.MtelH)
    title(ax,tend_str)
end
m.savefig_article(fig,'MuzikEtAl/results/time_end_reaction')