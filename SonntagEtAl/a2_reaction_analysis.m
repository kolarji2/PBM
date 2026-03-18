%
%   Determine reaction rate from experiments by fitting a reaction kinetics
%   model
%

load_dataset_paliperidone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=data_kin{1};
data.cexpPP_mean(1)=data.cPP0; %umol/l

x0=[6, 0.02];
fun=@(x) objfun_pp_reaction(x,data);
res=fminsearch(fun,x0);
fprintf('fminsearch res:\n');fprintf('%e, ',res);fprintf('\n');
res=[6.873,0.02];
sel=data.texp<1e9;
time=data.texp(sel);
t_integ=linspace(0,max(time),200);
[cPP,cP]= solve_pp_reaction(res, t_integ, data);
figure;
plot(t_integ,cPP,'-r')
hold on
errorbar(time,data.cexpPP_mean(sel),data.cexpPP_std(sel),'or')
plot(t_integ,cP,'-b')
errorbar(time,data.cexpP_mean(sel),data.cexpP_std(sel),'ob')

function [cPP,cP]= solve_pp_reaction(x, t_eval, data)
    %solve_pp_reaction Solve model for given kinetics parameters

    pars={};
    pars.k1=x(1);
    pars.kdeg=x(2);
    pars.M_P  = 426.5e-6; %g/umol Paliperidone
    pars.M_PP = 664.9e-6; %g/umol Paliperidone Palmitate
    pars.M_LiOH = 23.95e-6; % g/umol  
    c0=[data.cPP0.*pars.M_PP;0;data.cB0.*pars.M_LiOH]; %g/l
    options = odeset('RelTol',1e-6,'AbsTol',1e-8,'NonNegative',1); %
    [tt,yy] = ode15s(@(t,c) model_pp_reaction(t,c,pars),t_eval,c0,options);
    cPP=yy(:,1)./pars.M_PP;
    cP=yy(:,2)./pars.M_P;
end

function [err,cPP,cP]= objfun_pp_reaction(x,data)
    %OBJFUN_PP_REACTION Objective function to determine reaction rate
    
    sel=data.texp<1e9;
    t_eval=data.texp(sel);
    [cPP, cP] =solve_pp_reaction(x, t_eval, data);
    err=100*sum((cPP-data.cexpPP_mean(sel)).^2./data.cexpPP_std(sel)./data.cexpPP_mean(sel));
    err=err+sum((cP-data.cexpP_mean(sel)).^2./data.cexpP_std(sel)./data.cexpP_mean(sel));
end