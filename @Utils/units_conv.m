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

