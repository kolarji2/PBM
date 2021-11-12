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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BSD 3-Clause License
% 
% Copyright (c) 2021, Jiri Kolar, Suada Djukaj
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

