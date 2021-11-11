function valcaled = units_conv(obj,val,inUnits,outUnits,M,rho)
%units_conv (obj,val,inUnits,outUnits,M,rho): conversion between units
%   units: 
%          moll [mol/l]
%          gl [g/l]
%           kgm3 [kg/m3]
%           m3m3 [m3/m3] Volume of compound
%           ll [l/l] Volume of compound
%   M [g/mol]
%   rho [g/l] | [kg/m3] 

mol2g=M;
g2mol=1./M;
g2l= 1./rho; % == kg2m3
l2g= rho; % == m32kg

scale=1;
if strcmp(inUnits,'moll')    
    scale=scale.*mol2g;            
elseif strcmp(inUnits,'gl') || strcmp(inUnits,'kgm3')
    scale=scale.*1;        
elseif strcmp(inUnits,'m3m3') || strcmp(outUnits,'ll')
    scale=scale.*l2g;
else
     ME = MException('units_conv:UnknownInputUnits', ...
        'input units %s not supported',inUnits);
    throw(ME) 
end

% now units are all in g/l or kg/m3
if strcmp(outUnits,'moll')    
    scale=scale.*g2mol;            
elseif strcmp(outUnits,'gl') || strcmp(outUnits,'kgm3')
    scale=scale.*1;        
elseif strcmp(outUnits,'m3m3') || strcmp(outUnits,'ll')
    scale=scale.*g2l;
else
     ME = MException('units_conv:UnknownOutputUnits', ...
        'output units %s not supported',outUnits);
    throw(ME) 
end

valcaled=val*scale;
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

