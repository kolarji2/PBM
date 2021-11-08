function dydt = reactive_dissolve_fxfvm(obj,t,y)
% reactive_dissolve_fxfvm - 
%   solve PBE coupled with reaction kinetics
%   reuse dissolve_fxfvm to solve PBE with flux limited finite volume method
%   
%   Can not handle reaction in film ( Da << 1 )
%   Reaction rate and growth rate must have same time units!

% Separate conc and particle size dist
Nconc=obj.Nconc;
c=y(1:Nconc);

% Calculate reaction (any reaction model can be used here)
dcrdt=obj.model_reaction_h(t,c);

% Reuse PBE solver
dydt = obj.dissolve_fxfvm(t,y);

% Update derivatives
dydt(1:obj.Nconc)=dydt(1:obj.Nconc)+dcrdt;
end

