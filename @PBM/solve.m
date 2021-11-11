function [time,c,fn,phi_dist] = solve(obj,time)
% solve(obj,time): Call any solver specified in obj.solver_name and return results.

if isempty(obj.ode_options)
    options = odeset('RelTol',1e-6,'AbsTol',1e-8);
else
    options =obj.ode_options;
end

%options = odeset('RelTol',1e-8,'AbsTol',1e-10);
solverh=@(t,y) obj.(obj.solver_name)(t,y); %str2func

if contains(obj.solver_name,'fx') || contains(obj.solver_name,'fvx')
    y0=obj.y0fx;
else
    y0=obj.y0fn;
end

[tt,yy] = ode15s(solverh,time, y0,options); % %data.times

time=tt;
c=yy(:,1:obj.Nconc);
fn=yy(:,obj.Nconc+1:end)'; % To have same ordering as expdata, row=diams, colums = time steps
if contains(obj.solver_name,'fx')
    fn=fn./obj.x0; % fx convert to fn
end

Vpar=obj.x0.^3*obj.fV;
phi_dist=fn.*Vpar;
phi_dist=phi_dist./sum(phi_dist.*obj.dx,1);

obj.tm=time;
obj.cm=c;
obj.fnm=fn;
obj.phim_dist=phi_dist;
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

