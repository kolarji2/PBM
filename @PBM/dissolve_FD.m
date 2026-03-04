function dydt = dissolveFD(obj,t,y)
%dissolveFD(obj,t,y) Solve population balance with finite difference scheme
% fast, but very diffusive.
%
% tau - time
% states - state variables (c(t) and n(t,x)
% pars - model parameters
% x - position of discretized points for n(t,x)

% unpack state variables
c=y(1:obj.Nconc); % concentration, must be in kg/m3 or g/l (same value)!
n=y(obj.Nconc+1:end); % discretized population density

x=obj.x0;
dx=obj.dx;

if(size(n) ~= size(x))
    error('Size of vector "x" does not correspond to number of states')
end


% calculate mesh size
N = max(size(x));   % max returns the max element of an array
G=obj.Gfun(obj,t,c,x);

dndt = zeros(N,1);
for i=1:N-1
    dndt(i)=(G(i)*n(i)-G(i+1)*n(i+1))/dx(i);
end
dndt(N) = 0; %(G(N)*n(N)-0)/dx(N); % no particles from larger class

dcdt_particles = sum(-obj.vol2c.*obj.fV*x.^3.*dx.*dndt);

% pack derivatives of state variables
dcdt=zeros(size(c));
dcdt(1)=dcdt_particles;
dydt = [dcdt; dndt];
end