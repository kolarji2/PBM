%
%   Reactive dissolution of Paliperidone palmitate
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Reaction has no effect on mass transfer in film
%   Damkohler~0; enhancement factor~1%
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d0 = 1e-4;
    D = 7.6599e-10; % m2/s
    k1 = 6.863/3600/24; % 1/s
    Da = k1*d0^2/D;% == 0.001 ~ 0
    kM = 2*D/d0; % 1/s wihtout barier Sh*D/x
    ef=sqrt(D*k1/kM^2)*coth(sqrt(D*k1/kM^2));% reaction enhancment factor =1, reaction no influence;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Molar mass
M_P  = 426.5e-6; %g/umol Paliperidone
M_PP = 664.9e-6; %g/umol Paliperidone Palmitate
M_SDS = 288.372e-6; %g/umol
M_LiOH = 23.95e-6; % g/umol    

%% Load experimental data
readdata=true;
load_dataset_paliperidone
data=data_rd;

%% Init model
Nexp=length(data);
model_arr=cell(length(data),1);

for iexp=1:Nexp
    %% Create new model instance
    T=70+273.15;
    m=create_new_model_PP(T,M_PP);
    m.title=data{iexp}.expID;
    m.ode_options=odeset('RelTol',1e-6,'AbsTol',1e-8,'NonNegative',1);
    
    %% Optimized model pars   
    % estimated csat at basic pH
    m.csat =0.525; % [g/l] = [kg/m3]
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
    
%     m.csat=0.51731;
%     pars.k1=8.25474;
%     pars.kdeg=0.02036;
    
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
    tmax=18; %max(data{iexp}.texp);
    time=linspace(0,tmax,500);
    [time,c,fn,phi_dist] = m.solve(time);    
    %% Save
    model_arr{iexp}=m;
    
end

%% Visualize
% fig=figure('Position',[0 0 1500 900]);
% cols=lines(2*Nexp);
% conc_umol_l=[M_PP,M_P,M_LiOH];
% for iexp=1:4
%     ax=subplot(2,2,iexp);
%     m=model_arr{iexp};
%     
%     plot(ax,m.tm,m.cm(:,1)./conc_umol_l(1),'DisplayName',sprintf('PP - %s',m.title),'Color',cols(iexp,:))
%         hold(ax,'on')
% 	plot(ax,m.tm,m.cm(:,2)./conc_umol_l(2),'DisplayName',sprintf('P - %s',m.title),'Color',cols(iexp+1,:))
%     
%     plot(ax,data{iexp}.tscatter,data{iexp}.cscatterPP,'o','DisplayName',sprintf('PP exp - %s',m.title),'Color',cols(iexp,:))
%     plot(ax,data{iexp}.tscatter,data{iexp}.cscatterP,'o','DisplayName',sprintf('P exp - %s',m.title),'Color',cols(iexp+1,:))
%     legend(ax,'Location','best')    
%     xlabel('Time [days]')
%     ylabel('c [\mumol/l]')
% end

% Relative
fig=figure('Position',[0 0 1500 900]);
ax=cell(4,1);
for i=1:4
    ax{i}=subplot(2,2,i);
end
cols=lines(Nexp);
conc_umol_l=[M_PP,M_P,M_LiOH];
for iexp=1:Nexp    
    m=model_arr{iexp};   
    plot(ax{1},m.tm,m.cm(:,1)./m.cs0*100,'DisplayName',sprintf('PP - %s',m.title),'Color',cols(iexp,:))
    
	plot(ax{2},m.tm,m.cm(:,2)./m.cs0*100,'DisplayName',sprintf('P - %s',m.title),'Color',cols(iexp,:))
    
    plot(ax{3},data{iexp}.tscatter,data{iexp}.cscatterPP.*conc_umol_l(1)./m.cs0*100,'--o','DisplayName',sprintf('PP exp - %s',m.title),'Color',cols(iexp,:))
    plot(ax{4},data{iexp}.tscatter,data{iexp}.cscatterP.*conc_umol_l(2)./m.cs0*100,'--o','DisplayName',sprintf('P exp - %s',m.title),'Color',cols(iexp,:))
    for i=1:4
        hold(ax{i},'on')
    end
end
for i=1:4
    legend(ax{i},'Location','best')
    xlabel(ax{i},'Time [days]')
    ylabel(ax{i},'Released PP [%]')
    xlim(ax{i},[0,14])
end

