function [G,kM,kM0] = calc_G(obj,t,c,xb)
% calc_G(obj,t,c,xb): calculate mass transfer coefficient from selected correlation
% for each diameter in xb, include Ostwald effect to increase solubility of
% small particles.
%   
% Input:
%   t: time - (usually not used / can be used for area modification)
%   c: float or array - current concentration
%                       - c(1) - solute concentration
%                       - function handles get whole vector c
%   xb: vector of float - boundary particle diameter
%   obj: model parameters
% Output:
%   G: growth rate
%   kM: mass transfer
%
%
%   Note: several correction factors can be used, ussually are all set to 1. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%phim=0.63;
%cv=(obj.cS0-c).*obj.c2vol;
%etaSlurry=obj.etaf_fun(c).*(1+1.25*cv./(phim-cv)).^2;

if strcmp(obj.etaf_method,'etaf_fun')
    obj.etaf=obj.etaf_fun_h(c);
end

if strcmp(obj.rhof_method,'rhof_fun')
    obj.rhof=obj.rhof_fun_h(c);
end

if strcmp(obj.csat_method,'etaf_fun')
    obj.csat=obj.csat_fun_h(c);
end

if strcmp(obj.diff_method,'diff_wilkechang')
    obj.diff=obj.wilkechang(obj.etaf,obj.T,obj.assoc_par,obj.Mf,obj.Vmol_solute);
end

Sc = obj.etaf./(obj.diffcorr.*obj.diff.*obj.rhof);
kM0=ones(size(xb));

if (strcmp(obj.Sh_method,'settling') || strcmp(obj.Sh_method,'levins')) & ~isfield(obj,'vsettle')
    warning('vsettle not set in obj, calculating vsettle...');
    vsettle = calc_vsettle(xb,obj);
    obj.vsettle=vsettle;
end

if strcmp(obj.Sh_method,'settling')
    vsettle=obj.vsettle;
    vt=obj.vtcorr*obj.vsettle;
    velo=obj.vcorr*sqrt(obj.v.^2+vt.^2);
    if length(xb)~=length(velo)
       velo=interp1(obj.xb,velo,xb); 
    end
    Rey= velo.*xb.*obj.rhof./obj.etaf;
    %Sh = 2+0.6*Rey.^0.5*Sc.^(1/3); % Ranz-Marshall Correlation for sphere
    Sh = 2+0.44.*Rey.^0.5.*Sc.^(0.38); 
    % Maybe Nienow and Miles (1978) developed the Froessling equation based on the slip velocity theory  - Could find the equation in the original paper
    % Levins and Glastonbury equation  ( mentioned in Pangarkar2002) - original paper is unavailable
elseif strcmp(obj.Sh_method,'levins')
    % Mass transfer defined by Levins and Glastonbury 1972
    Reeps= obj.epsilon^(1/3)*xb.^(4/3)*obj.rhof/obj.etaf;
    Shp=2+0.47*(Reeps.*obj.DsDt.^0.28).^0.62.*Sc.^0.36;
    kM0= obj.kmcorr*Shp.*obj.diff./xb; % for particles with zero density difference
    ve=0.93*obj.etaf/obj.rhof./xb.*Reeps.^1.23*obj.DsDt.^0.35; % effective velocity
    vt=obj.vtcorr*obj.vsettle;
    vs=obj.levinsSlipCorr.*vt; % slip velocity
    velo=obj.vcorr*sqrt(ve.^2+vt.^2+vs.^2);
    Rep= velo.*xb*obj.rhof/obj.etaf;
    Sh = 2+0.44*Rep.^0.5*Sc.^(0.38);
elseif strcmp(obj.Sh_method,'kolmogorov')
    % Based on Kolmogorov theory of isotropic turbulence
    % Levins Glastonbury 1972 Application of Kolgomorov theory to mass transfer
    % Default obj from  Piero Armenante Donald Kirwan Mass transfer to
    % microparticles in agitated systems 1989
%     obj.kolmogorov=0.52;%0.52;
%     obj.kolmogorov=0.52; %0.52;
%     obj.kolmogorov=1/3; %%1/3;
    Reeps= obj.epsilon^(1/3)*xb.^(4/3)*obj.rhof/obj.etaf;
    Sh=2+obj.kolmogorovA*Reeps.^(obj.kolmogorovB).*Sc.^(obj.kolmogorovC);    
elseif  strcmp(obj.Sh_method,'archimedes')
    Ar=obj.g*xb.^3*obj.rhof*(obj.rhos-obj.rhof)/obj.etaf^2;
    N=obj.impeler_rpm/60;% impeler speed 1/s (frequency) RPS=RPM/60;
    D=obj.impeler_diam; % impeler diameter [m]
    Rep=obj.rhof*pi*xb./obj.etaf*N*D;    
    Sh = 2+0.0096*Rep.^(2/3).*Sc.^(1/3).*Ar.^(0.16);
end

kM = obj.kmcorr*obj.diffcorr*Sh.*obj.diff./xb;

% G [m/s]
% dc -> m3/m3
if strcmp(obj.csat_method,'csat_fun')
    csat=obj.csat_fun_h(c);
elseif strcmp(obj.csat_method,'csat_const')
    csat=obj.csat;
end

crel=c(1)/obj.cs0;
G = obj.fAfun(crel,xb)./(3*obj.fV).*obj.c2vol.*kM.*(c(1)-csat);
    
if isnan(G(1)) || isinf(G(1)) % correction if xb(1)=0
    G(1)=G(2);
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

