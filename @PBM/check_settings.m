function warn = check_settings(obj)
%CHECK_SETTINGS - run diagnostic check if all params are set
warn=false;
%% Check if known solver is used
warn=warn|check_list_options(obj,'solver_name',obj.solver_name_list);

warn=warn|check_property_defined(obj,obj.check_matter_list, ...
    'Missing physico-chemical properties:');

%% Generic model properties
warn=warn|check_property_defined(obj,obj.check_model_list, ...
    'Missing properties for model solution:');

%% Wilke-Chang Method
if contains(obj.diff_method,'wilkechang')
    warn=warn|check_property_defined(obj,obj.check_diff_wilkechang_list, ...
        'Checking properties for Wilke-Chang method...');
end

%% Check handles
for i=1:size(obj.check_handles_list,1)
    opt=obj.check_handles_list{i,1};
    value=obj.check_handles_list{i,2};
    prop=obj.check_handles_list{i,3};
    if strcmp(obj.(opt),value)
         warn=warn|scheck_property_defined(obj,{prop}, ...
            sprintf('%s: "%s" needs following handle:',opt,value));
    end
end

if ~warn
    warning('off','backtrace');
    warning(sprintf('\n\nModel check completed:\n\tAll seems OK\n\n'));
    warning('on','backtrace')
end
end

function warn=check_list_options(obj,opt,list)
% Warning if property does not have suitable value from list
    warn=~any(strcmp(obj.(opt),list));
    if warn
        warning('off','backtrace');
        values=sprintf('\t%s\n',list{:})
        warning('Invalid value for property: "%s"\nPossible values:\n%s',opt,values)
        warning('on','backtrace')
    end
end

function warn=check_property_defined(obj,list_properties,info)
warning('off','backtrace')
    warn=false;
    warn_string=sprintf('%s\n',info);
    for i=1:length(list_properties)
       if isempty(obj.(list_properties{i}))           
           warn_string=[warn_string sprintf('Missing property: %s\n',list_properties{i})];
           warn=true;
       end
    end
if warn
    warning(warn_string);
end
warning('on','backtrace')
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

