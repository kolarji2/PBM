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
end