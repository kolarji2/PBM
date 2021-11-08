function [pdf0] = distfunNormal(obj,x,dx,mu,sigma,minx,maxx)
%distfunNormal (obj,x,dx,mu,sigma,minx,maxx): Generate pdf normal dist. function.
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
    pdf0 = pdf('normal',x,mu,sigma);
    pdf0(x<minx)=0;
    pdf0(x>maxx)=0;
    if sum(pdf0.*dx)>0
        pdf0=pdf0./sum(pdf0.*dx);
    end    
end
end

