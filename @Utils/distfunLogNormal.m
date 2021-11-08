function [pdf0] = distfunLogNormal(obj,x,dx,mu,sigma,minx,maxx)
%distfunLogNormal (obj,x,dx,mu,sigma,minx,maxx): Generate pdf function.
% Input:
%   x - cell centers
%   dx - width of cells, used for estimation of integral as sum(pdf.*dx)
%   mu - mean size
%   sigma - standard deviation
%   minx - minimum size (anything below is set to 0)
%   maxx - maximum size (anything above is set to 0)
%
%   Output
%   pdf0 - discretized probability density function sum(pdf0.*dx)=1
%
if mu<=0
    pdf0=zeros(size(x));
else
    mux=log(mu^2/sqrt(sigma.^2+mu^2));
    sigmax=sqrt(log(1+sigma^2/mu^2));
    pdf0 = pdf('lognormal',x,mux,sigmax);
    pdf0(x<minx)=0;
    pdf0(x>maxx)=0;
    if sum(pdf0.*dx)>0
        pdf0=pdf0./sum(pdf0.*dx);
    end    
end
end

