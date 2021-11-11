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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BSD 3-Clause License
% 
% Copyright (c) 2021, Jiri Kolar
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

