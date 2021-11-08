function [pdf0] = distfunMultiModal(obj,x,dx,w,mu,sigma,dist0)
%distfunMultiModal (obj,x,dx,w,mu,sigma,dist0): calculate multimodal distribution from dist functions
%arrays of values for each dist: w,mu,sigma
% Input:
%   x - cell centers
%   dx - width of cells, used for estimation of integral as sum(pdf.*dx)
%   w - array of fractions, weights how distributions are combined
%   mu - array of mean sizes
%   sigma - array of standard deviations
%   dist0 - probability density function of @(x,dx,mui,sigmai)
%
%   Output
%   pdf0 - discretized probability density function sum(pdf0.*dx)=1
%

Ndist=length(w);
pdf0=zeros(size(x));
for i=1:Ndist
    pdf0=pdf0+w(i)*dist0(x,dx,mu(i),sigma(i));
end
if sum(pdf0.*dx)>0
        pdf0=pdf0./sum(pdf0.*dx);
end    

end

