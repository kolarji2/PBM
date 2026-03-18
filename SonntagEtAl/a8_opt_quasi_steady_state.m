%
%   Ideal particle size for steady state dissolution (optimization)
%

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
opt_model_arr=cell(Nc,Nk); 
cs_target_arr=zeros(Nc,Nk);

for ik=1:Nk
    parfor ic=1:Nc
    
    %% Create new model instance
    T=70+273.15;
    m=create_new_model_PP(T,M_PP);
    m.title='opt';
    m.ode_options=odeset('RelTol',1e-4,'AbsTol',1e-6,'NonNegative',1);
    m.ode_flux_limiter='none';
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
    mu_arr=[1e-6 100e-6];
    w=[0,1];
    sigma_arr=mu_arr.*0.1;
    m.set_init_dist(mu_arr,sigma_arr,w,[],'lognormal'); 
    cs0=cs0_arr(ic); %msusp/V0; %[g/l] = [kg/m3]    
    
    m.Nconc=3;        
    cB0=1e9;%  place holder
    cstart=[0;0;cB0]; % initiate concentration [cPP,cP,cLiOH];
    m.set_initial_conditions(cstart,cs0,'gl');
    
    %% Solve model
    tmax=tmax_arr(ik);
    time=linspace(0,tmax,400);
    [time,c,fn,phi_dist] = m.solve(time);    
    cs_target_arr(ic,ik)=max(m.cm(:,1));
    %% Save
    opt_model_arr{ic,ik}=m;
    opt_data_arr{ic,ik}=struct;
    
    end
end


%% Optimization to find opt csat for experiment with reaction
w100um=zeros(Nc,Nk);
tstable_arr=zeros(Nc,Nk);
for ik=1:Nk
    for ic=1:Nc
        ik
        ic
    opt=OptimizePBM;
    %ic=1;
    %ik=3;
    opt.models=opt_model_arr(ic,ik);
    opt.solver_method='fminbnd';
    opt.TolX=1e-5;
    opt.swarm_size=6;
    opt.objfun_method='handle';
    cs_threshold=0.05;
    time=linspace(0,100,400);
    cmax_target=cs_target_arr(ic,ik);
    opt.objfun_h=@(x,models,expdata) objfun_elimkin_simple0(x,time,models{1},cmax_target,cs_threshold);
    myobjfun=@(x) objfun_elimkin_simple0(x,time,opt.models{1},cmax_target,cs_threshold);

    %w1
    xarr=linspace(0,1,21);
    errval=zeros(size(xarr));
    parfor i=1:length(xarr)
        errval(i)=myobjfun(xarr(i));
    end
    figure
    plot(xarr,errval)
    [~,imin]=min(errval);   
    opt.x0=[0.9];
    opt.lb=[xarr(max(1,imin-2))];
    opt.ub=[1];

    xres=opt.optimize([],[]);
    w100um(ic,ik)=xres;
    [err,tmax,tmin,m]=myobjfun(xres);
    tstable_arr(ic,ik)=tmax-tmin;
    
    %plot
    figure
    plot(m.tm,m.cm(:,1))
    hold on
    [~,~,~,m2]=myobjfun(1);
    plot(m.tm,m2.cm(:,1),'g')
    cs_target=max(m.cm(:,1))./(1+cs_threshold);
    cmax_target=cs_target*(1+cs_threshold);   
    cmin_target=cs_target*(1-cs_threshold); 
    hold on
    plot([0 max(m.tm)],[cmax_target cmax_target],'--r')
    plot([0 max(m.tm)],[cmin_target cmin_target],'--r')
    %drawnow();
    end
end
m=opt_model_arr{1,1};
mu_arr=m.mu_arr;
sigma_arr=m.sigma_arr;
save('SonntagEtAl/expdata/opt1_100um_fine.mat','w100um','tstable_arr','mu_arr','sigma_arr','opt_model_arr');



function [err,tmax,tmin,m]=objfun_elimkin_simple0(x,time,m,cmax_target0,cs_threshold)
    w=[1-sum(x) x];
    if w(1)<0
        w(1)=0;
    end    
    w=w./sum(w);    
    m.w=w;
    m.set_initial_conditions(m.cstart,m.cs0,m.conc_units);
    [tt,c] = m.solve(time);
    cs_target=max(m.cm(:,1))./(1+cs_threshold);
    cmax_target=cs_target*(1+cs_threshold);   
    cmin_target=cs_target*(1-cs_threshold);   
    sel=time>0;
    cm=c(:,1);
    [cmax,cmaxi]=max(cm); 
    sel=tt>=tt(cmaxi) & cm>0;    
    sel2=tt<=tt(cmaxi);
    if cmin_target<cm(end)
        tmax=max(time);
    else
        tmax=interp1(cm(sel),tt(sel),cmin_target);
    end    
    penalty=1;
    if abs(cmax_target-m.csat)<0.02
        penalty=0;
    end
%     if cmax>cmax_target0*1.1
%         penalty2=(cmax-cmax_target0*1.1)./cmax_target0/1.1*10;
%     else
%         penalty2=0;
%     end
    
    if sum(sel2)>1
        tmin=interp1(cm(sel2),tt(sel2),cmin_target);
    else
        tmin=0;
    end
    err=-tmax*penalty+tmin; %*10+penalty2;    
end