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

