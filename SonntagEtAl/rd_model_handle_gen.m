function fun_h = rd_model_handle_gen(x)
%RD_MODEL_HANDLE_GEN For optimization by generic optimizer
%   Detailed explanation goes here
    pars={};
    pars.k1=x(2);
    pars.kdeg=x(3);
    pars.M_P  = 426.5e-6; %g/umol Paliperidone
    pars.M_PP = 664.9e-6; %g/umol Paliperidone Palmitate
    pars.M_LiOH = 23.95e-6; % g/umol    
    fun_h=@(t,c) model_pp_reaction(t,c,pars);
end