function [csat,cHsat] = calc_csat(obj,cTintrinsic,cCO3tot,cNa)
%calc_csat(obj,cTintrinsic,cCO3tot,cNa): Calculate csat, slow!
%   Estimates saturated concentration
%   Calculate maximum possible concentration at solid-liquid boundary
%   including consideration of pH change necessary to obtain this
%   concentration
%
%   Slow calculation!
%   
%   Inputs:
%       intrinsicS0: intrinsic solubility of solute
%       cCO3tot
%       cNa
%   Returns:
%       csat: maximum possible concentration at solid boundary
%       cHsat: pH at solid boundary

cTmax=1;
cTmin=0;


cH_cTmax=obj.solve_charge_balance(cTmax,cCO3tot,cNa);
cH_cTmin=obj.solve_charge_balance(cTmin,cCO3tot,cNa);


cTsolubility=@(cH) cTintrinsic.*obj.get_ceqred(cH);

eqcT=@(cT) cT-cTsolubility(obj.solve_charge_balance(cT,cCO3tot,cNa));

amax=eqcT(cTmax);
amin=eqcT(cTmin);

if amax>0 && amin>=0
    csat=cTmin;
    cHsat=cH_cTmin;
    return
elseif amax<0 && amin<0
   csat=cTmax;
   cHsat=cH_cTmax;
else
    tol=max(1e-15,cTintrinsic/1e3);
    csat=obj.fzerobnd(eqcT,cTmin,cTmax,tol);
    cHsat=obj.solve_charge_balance(csat,cCO3tot,cNa);
    if csat>cTmax
       csat=cTmax;
       cHsat=cH_cTmax;      
    elseif csat<cTmin
        csat=cTmin;
        cHsat=cH_cTmin;
    end
end

end

