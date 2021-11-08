function diff_WC = wilkechang(obj,etaf,T,assoc_par,Mf,Vmol_solute)
%wilkechang(obj,etaf,T,assoc_par,Mf,Vmol_solute): calc difusivity
% Calculate diffusivity using Wilke-Chang method:
% D=7.4e-8*(x*Msolvent)^0.5*T/eta_solvent/Vmol_solute^0.6/1e4 ; eta [mPas]


%Estimation of molecular diffusivity in liquid phase systems by the
% Wilke–Chang equation
diff_WC=7.4e-8*(assoc_par*Mf)^0.5*T./etaf/1e3/(Vmol_solute)^0.6/1e4;
%modified Wilke–Chang equation assoc_par is used as well for Vmol (Miyabe Journal
%of Chromatography A, 1218 (2011) 6639– 6645)
%diff_WCmod=7.4e-8*(obj.assoc_par*obj.Mf)^0.5*T/etaf/1e3/(obj.assoc_par*obj.Vmol_solute)^0.6/1e4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

