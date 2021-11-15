classdef Utils < handle
%% Utils
% Collection of functions used by the model

%% methods of utils
methods       
    
function obj=Utils(obj)
end

function [ dx ] = get_dx_fromBounds(obj, boundaries )
%get_dx_fromBounds (obj, boundaries ): Get meshspacing of size classes
    dx = boundaries(2:end)-boundaries(1:end-1);
end

function [ dx ] = get_dx_include0(obj, x )
%get_dx_include0 ( x ) Get meshspacing of size classes including point 0
    if size(x,1)<=1
       error('use_column_vectors');
    end
    dx=x-[0;x(1:end-1)];
end

function [xclass] = get_xclass_fromBounds(obj,xbound)
%get_xclass_fromBounds (obj,xbound): return class centers from bounds
    xclass=0.5*(xbound(1:end-1)+xbound(2:end));
end

function [xclass] = get_xclass_include0(obj,xbound)
%get_xclass_include0 (obj,xbound): return class centers from bounds, including point 0
    if size(xbound,1)<=1
       error('use_column_vectors');
    end
    xbound=[0;xbound];
    xclass=0.5*(xbound(1:end-1)+xbound(2:end));
end

function fc=get_mean_integral_values(obj,xsmooth,fsmooth,xb)
%mean_integral_values(obj,xsmooth,fsmooth,xb): Calcalutes mean integral value
%   Evaluate function in new points (coarse discretization)
    fc=zeros(length(xb)-1,1);
    Np=1000;
    for i=1:length(xb)-1
        x=linspace(xb(i),xb(i+1),Np);
        f=interp1(xsmooth,fsmooth,x);
        fc(i)=trapz(x,f)./(x(end)-x(1));
    end

end

function phim_scaled=get_scaled_dist(obj,x,dx,phim_dist,xb_coarse)
   %get_scaled_dist(obj,x,dx,phim_dist,xb_coarse)
   %Scale fine distribution so it look comparably to a histogram on a coarse mesh
   %    Warning - works best if each coarse interval contains approximately
   %              same number of points from the fine discretization.
   %
   % Input:
   %    x - cell centers
   %    dx - width of cells
   %    phim_dist - distribution function (integral or sum(phim_dist.*dx)=1)
   %    xb_coarse - boundaries for coarse mesh (excluding 0)
   %
   % Output:
   %    phim_scaled - scaled histogram-like distribution 
   %            sum(phim_scaled)!=1 however it looks comparably to a histogram
   %            on a given coarse mesh
   %
    xb_coarse=[0;reshape(xb_coarse,length(xb_coarse),1)];   
    xexp=obj.get_xclass_fromBounds(xb_coarse);
    [Nvals_in_bin,~,bins]=histcounts(x,xb_coarse);    
    scale=interp1(xexp,Nvals_in_bin,x,'linear',mean(Nvals_in_bin));
    phim_scaled=phim_dist.*dx.*mean(scale);
end


end % end of methods

methods(Static)
    function savefig_article(fig,fpath)
    %savefig_article (obj,fig,fpath): save any figure in .fig and .png format
        savefig(fig,[fpath '.fig']);
        print(fig,'-dpng','-r150',[fpath '.png']);
    end
end

end % end of class Utils


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

