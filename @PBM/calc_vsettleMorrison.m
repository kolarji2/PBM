function  vsettle = calc_vsettleMorrison(obj,xb)
%calc_vsettleMorrison (obj,xb): Calculation of settling velocity
% Calculation of settling velocity according to Solid-Liquid Mixing Chapter
% 10 V.A. Atiemo-Obeng, W.R. Penney, P. Armenante, Handbook of Industrial
% Mixing: Science and Practice, Wiley 2004
acc=obj.acc;
rhos=obj.rhos;
rhof=obj.rhof;
etaf=obj.etaf;


vsettle=zeros(length(xb),1);
if acc>0
    vtLam=acc*(rhos-rhof)*xb.^2/18/etaf; % settling velocity

    Reu=@(v,d) v.*d*rhof/etaf;
    %https://pages.mtu.edu/~fmorriso/DataCorrelationForSphereDrag2016.pdf
    Cd_morrison=@(v,d) 24./Reu(v,d)+ 2.6*(Reu(v,d)./5)./(1+(Reu(v,d)./5).^1.52) ...
                        + 0.411*(Reu(v,d)/2.63e5).^-7.94./(1+(Reu(v,d)/2.63e5).^-8) ...
                        + 0.25*Reu(v,d)/1e6./(1+Reu(v,d)/1e6); % Morrison correlation for whole range
    vtInter=zeros(length(xb),1);
    vt_res=@(v,d) acc*(rhos-rhof)/rhof-Cd_morrison(v,d).*3/4./d.*v.^2;
    for i=1:length(xb)  
        vsettle(i)=fzerobnd(@(v) vt_res(v,xb(i)),1e-20,2,1e-20);
    end
end
end

