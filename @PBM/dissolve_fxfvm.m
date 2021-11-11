

function dydt = dissolve_fxfvm(obj,t,y)
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


x0=obj.x0*obj.ode_scalepar; % cells center position (N,1)
xb=obj.xb*obj.ode_scalepar; % cells boundary position (N+1)
dx=obj.dx*obj.ode_scalepar; % cell width
dxii=x0(2:end)-x0(1:end-1); %distance between cell centers
N=obj.N; % number of cell

c=y(1:obj.Nconc); % concentration, must be in kg/m3 or g/l (same value)!
fx=y(obj.Nconc+1:end);
G=obj.Gfun(obj,t,c,obj.xb).*obj.ode_scalepar; % calculate growth rate (dxdt) at boundaries xb
Gf=zeros(N+1,1); %G*fx


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Different flux limiters can be used:
%
%flux limiting function
if obj.ode_flux_limiter==0
    phir= @(r) 0; %no limiter
else
%phir= @(r) (abs(r)+r)/(1+abs(r)); %van Leer
    phir= @(r) max(0,min(2*r,min(1/3+2/3*r,2))); % Best performance for normal dist
end
%phir= @(r) max([0,min(2*r,1),min(r,2)]); %super bee
%phir= @(r) max(0,min(1,r)); % minmod
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
        ri=(fx(2)-fx(1)+eps)/((fx(1)-0)+eps);    
        Gf(2)=G(2)*(fx(1)+0.5*phir(ri)*dx(1)/dxii(1)*(fx(2)-fx(1)));   
        Gf(N+1)=G(N+1)*(fx(N)+0.5*dx(N)/dxii(N-1)*(fx(N)-fx(N-1)));
    else
        Gf(1)=G(1)*(fx(1)-0.5*(fx(2)-fx(1))*dx(1)/dxii(1));
        Gf(N)=G(N)*(fx(N)-0.5*(fx(N)-fx(N-1))*dx(N)/dxii(N-1));
        Gf(N+1)=0; %no addition of large particles
    end
    if G(1)>0
        for i=3:N
            ic=i-1;
            %ri upwind ratio of two consecutive solution gradients
            ri=((fx(ic+1)-fx(ic))+eps)/((fx(ic)-fx(ic-1))+eps)*dxii(ic-1)/dxii(ic); 
            fxbound=(fx(ic)+0.5*phir(ri)*dx(ic)/dxii(ic)*(fx(ic)-fx(ic-1)));
            if fxbound<0
                fxbound=0;
            end
            Gf(i)=G(i)*fxbound;    
        end
    else %negative growth
        for i=2:N-1
            ic=i-1;
            %ri upwind ratio of two consecutive solution gradients
            ri=((fx(ic)-fx(ic+1))+eps)/((fx(ic+1)-fx(ic+2))+eps)*dxii(ic+1)/dxii(ic);
            fxbound=(fx(ic+1)+0.5*phir(ri)*dx(ic+1)/dxii(ic+1)*(fx(ic+1)-fx(ic+2)));    
            Gf(i)=G(i)*fxbound;    
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Final mass balance for particles
%   Last term arise from initial multiplication of diff equation by x
%
dfxdt=(Gf(1:end-1)-Gf(2:end))./dx+0.5*(G(1:end-1)+G(2:end)).*fx./x0;

dVdt=sum(-obj.fV.*dfxdt.*x0.^2.*dx);    
dcdt_particles = obj.vol2c.*dVdt./obj.ode_scalepar^3;
dcdt=zeros(size(c));
dcdt(1)=dcdt_particles;
dydt=[dcdt;dfxdt];
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

