function dydt = dissolve_qamar_fxfvm(obj,t,y)
%dissolve_fxfvm(obj,t,y) Solve population balance with finite volume method
%   Use flux limiter to reduce numerical diffusion
%
%   use modified population balance: fn*x
%   (Number based particle size distribution multiplied by x)
%   This method is taken from article:
%   S. Qamar et al.: On the solution of population balances for nucleation,
%   growth, aggregation and  breakage processes (Chemical Engineering
%   Science 2009)
%
%
%   In original paper x should be volume. However works as well for x as diamerer.
%   If x as volume would be supplied, calculation of concentration must be modified!
%
%   It seems to be faster than regular FVM where Number based particle size
%   distribution is directly solved.
%
%   Note: Several parameters must be set in PBM class (object obj), see below!
%
%
% x0=obj.x0*obj.ode_scalepar; % cells center position (N,1)
% xb=obj.xb*obj.ode_scalepar; % cells boundary position (N+1)
% dx=obj.dx*obj.ode_scalepar;

% ode_x0=obj.x0.^obj.ode_xorder*obj.ode_scalepar;
% ode_xb=obj.xb.^obj.ode_xorder*obj.ode_scalepar;
% ode_dx=xb(2:end)-xb(1:end-1); % cell width
% ode_dxii=x0m(2:end)-x0m(1:end-1);
% ode_Gscale=obj.xb.^(obj.ode_xorder-1).*obj.ode_scalepar;
% ode_Vscale=x0.^(3-obj.ode_xorder).*dx./obj.ode_scalepar^3;

x0=obj.ode_x0; %obj.x0.^obj.ode_xorder*obj.ode_scalepar; % cells center position (N,1)
xb=obj.ode_xb;%obj.xb.^obj.ode_xorder*obj.ode_scalepar; % cells boundary position (N+1)
dx=obj.ode_dx;%xb(2:end)-xb(1:end-1); % cell width
dxii=obj.ode_dxii; %x0(2:end)-x0(1:end-1); %distance between cell centers
N=obj.N; % number of cell

c=y(1:obj.Nconc); % concentration, must be in kg/m3 or g/l (same value)!
fx=y(obj.Nconc+1:end);
ode_Gscale=(obj.x0.*obj.ode_scalepar).^(obj.ode_xorder-1).*obj.ode_scalepar;
G=obj.Gfun(obj,t,c,obj.x0).*ode_Gscale; %obj.xb.^(obj.ode_xorder-1).*obj.ode_scalepar; % calculate growth rate (dxdt) at boundaries xb
Gf=zeros(N+1,1); %G*fx
Gfc=G.*fx; %value in cell center

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Different flux limiters can be used:
%
%flux limiting function
phir=obj.ode_flux_limiter_h;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Different formulation for positive and negative growth rate
%   Properly tested and validated only for negative growth rate
%   (dissolution)
%
    eps=1e-25; % must be very small number
    if G(1)>0
        Gf(1)=0; %no nucleation
        %Gf(2)=G(1)*(fx(1));
        ri=(Gfc(2)-Gfc(1)+eps)/((Gfc(1)-0)+eps);    
        Gf(2)=(Gfc(1)+0.5*phir(ri)*dx(1)/dxii(1)*(Gfc(2)-Gfc(1)));   
        Gf(N+1)=(Gfc(N)+0.5*dx(N)/dxii(N-1)*(Gfc(N)-Gfc(N-1)));
    else
        Gf(1)=(Gfc(1)-0.5*(Gfc(2)-Gfc(1))*dx(1)/dxii(1));
        Gf(N)=(Gfc(N)-0.5*(Gfc(N)-Gfc(N-1))*dx(N)/dxii(N-1));
        Gf(N+1)=0; %no addition of large particles
    end
    if G(1)>0
        for i=3:N
            ic=i-1;
            %ri upwind ratio of two consecutive solution gradients
            ri=((Gfc(ic+1)-Gfc(ic))+eps)/((Gfc(ic)-Gfc(ic-1))+eps)*dxii(ic-1)/dxii(ic); 
            gfxbound=(Gfc(ic)+0.5*phir(ri)*dx(ic)/dxii(ic)*(Gfc(ic)-Gfc(ic-1)));
            if gfxbound<0
                gfxbound=0;
            end
            Gf(i)=gfxbound;    
        end
    else %negative growth
        for i=2:N-1
            ic=i-1;
            %ri upwind ratio of two consecutive solution gradients
            ri=((Gfc(ic)-Gfc(ic+1))+eps)/((Gfc(ic+1)-Gfc(ic+2))+eps)*dxii(ic+1)/dxii(ic);            
            Gf(i)=(Gfc(ic+1)+0.5*phir(ri)*dx(ic+1)/dxii(ic+1)*(Gfc(ic+1)-Gfc(ic+2)));
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Final mass balance for particles
%   Last term arise from initial multiplication of diff equation by x
%
dfxdt=(Gf(1:end-1)-Gf(2:end))./dx+0.5*(Gf(1:end-1)+Gf(2:end))./x0;
dVdt=sum(-obj.fV.*dfxdt.*dx.*obj.x0.^(3-obj.ode_xorder).*obj.ode_Vscale);    
dcdt_particles = obj.vol2c.*dVdt;
dcdt=zeros(size(c));
dcdt(1)=dcdt_particles;
dydt=[dcdt;dfxdt];
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

