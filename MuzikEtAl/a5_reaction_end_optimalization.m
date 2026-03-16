%
%   Example of usage of modified PBM class for 
%		modelling of pH dependent dissolution of telmisartan.
%		
%		Modified class can be found in directory: @PBMTelmi 
%   
% Optimalization of reaction time
%  - compute time of reaction for increasing initial telmisartan
%  concentration
%  - Assume constant 1.05 molar excess of Na2CO3

m=PBMTelmi();
loadData=0;
if ~exist('TimeEnd45','var')
    [exps,TimeEnd45,TimeEnd55,TimeEnd65,TimeEndData,prop] = load_exp_data_telmisartan(m);
    [phi0expdata] = load_psd_telmisartan();
end

cexpTelmi=logspace(1,log10(300),200); %g/l
excessB=1.05*ones(size(cexpTelmi)); % molar excess
tend0=120*ones(size(cexpTelmi)); %h

wpure_Na2CO3=0.928;
wpure_Telmi=0.996;

Nexp=length(cexpTelmi);
model_arr=cell(Nexp,1);
tend_arr=zeros(Nexp,1);
cn_arr=zeros(Nexp,2);

% Track progress
global n_completed
n_completed = 0;
progress_q = parallel.pool.DataQueue;
afterEach(progress_q, @(~) update_progress(Nexp));


% Loop
parfor iexp=1:Nexp
    m=PBMTelmi();
        
    Texp=65+273.15;
    cS0=cexpTelmi(iexp)./m.MtelH;
    cNa2CO3=excessB(iexp).*cS0.*wpure_Telmi;
    cS0=cS0*wpure_Telmi;
    cn_arr(iexp,:)=[cS0,cNa2CO3];
    % mild agglomeration, bigger PSD for concentrated solutions
    w=[0,0,1];

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
    m.epsilon=0.01; % m2/s3
    m.ode_scalepar=1e6;
    m.ode_flux_limiter='Koren';   % Koren | None       
    
    % Case 1: Normal CO2 release
    m.kMaCO2=5e-4; % 1/s 
    % Case 2: Enhanced CO2 release
    %m.kMaCO2=20e-4; % 1/s 
    
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
    time=[0, logspace(-1,log10(tmax),500)];
    m.ode_options=odeset('RelTol',1e-4,'AbsTol',1e-6,'Events',@(t,y) stop_dissolution(t,y,m,0.9999));
    [time,c,fn,phi_dist] = m.solve(time);    
    %% Save
    model_arr{iexp}=m;

    %% Calc tend
    cr=m.cm(:,1)./m.cs0;
    th=m.tm;
    tend_est=interp1(cr(cr<0.9999),th(cr<0.9999),0.975);
    tend_arr(iexp)=tend_est
    send(progress_q, iexp);
end

cm_arr=cn_arr.*[m.MtelH,m.MNa2CO3];
matfile_path=sprintf('MuzikEtAl/results/time_end_ka_%d.mat',m.kMaCO2*1e4);
save(matfile_path,'tend_arr','cn_arr','cm_arr');

%
fig=figure('Position',[100,100,800,800]);
plot(cn_arr(:,1).*m.MtelH,tend_arr./3600)

function [position,isterminal,direction] = stop_dissolution(t,y,m,eps)
  position = y(1)-m.cs0*eps; % The value that we want to be zero
  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end

function update_progress(Nexp)
    global n_completed;
    n_completed = n_completed + 1;
    fprintf('Progress: %.1f %%\n', n_completed/Nexp*100);
end