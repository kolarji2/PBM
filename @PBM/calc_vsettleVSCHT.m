function  vsettle = calc_vsettleVSCHT(obj,xb)
%calc_vsettleVSCHT (obj,xb): Calculation of settling velocity according to lecture notes from VSCHT.
acc=obj.acc;
rhos=obj.rhos;
rhof=obj.rhof;
etaf=obj.etaf;

if acc==0
    vsettle=zeros(length(xb),1);
else
    vtLam=acc*(rhos-rhof)*xb.^2/18/etaf; % settling velocity

    Reu=@(v,d) v.*d*rhof/etaf;
    Cd_inter=@(v,d) 24./Reu(v,d).*(1+0.125*Reu(v,d).^0.72);

    vtInter=zeros(length(xb),1);
    vtInter_res=@(v,d) acc*(rhos-rhof)/rhof-Cd_inter(v,d).*3/4./d.*v.^2;
    for i=1:length(xb)
        if xb(i)<1e-9
            continue
        end
        vtInter(i)=fzerobnd(@(v) vtInter_res(v,xb(i)),1e-10,1,1e-9);
    end
    vtTurb=1.74*sqrt(acc*xb*(rhos-rhof)/rhof);

    ReuLam=Reu(vtLam,xb);
    ReuInter=Reu(vtInter,xb);
    ReuTurb=Reu(vtTurb,xb);
%     if any(ReuInter>1e3)
%         error('Reynolds number too high, settling in turbulent regime not implemented.')
%     end
    sellam=ReuLam<0.2;
    sellinter=ReuInter<1000;
    vsettle=zeros(length(xb),1);
    vsettle(sellam)=vtLam(sellam);
    vsettle(~sellam&sellinter)=vtInter(~sellam&sellinter);
    vsettle(~sellinter)=vtTurb(~sellinter);
end
end

