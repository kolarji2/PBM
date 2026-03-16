classdef OptimizePBM < matlab.mixin.Copyable
%% OptimizePBM - Class for finding optimal parameters of the model
    %  
    %% About OptimizePBM
    % 
    % OptimizePBM is used to simplify optimization on multiple experimental data
    % 
    %   
    %%  Usage
    %   Important note, Matlab always pass class as a reference. So all PBM
    %   classes passed to variable "models" are still references to the same
    %   object. 
    %
    %   This means that evaluation of the objective function change model
    %   parameters!
    %   
    
    %% PBM Properties
        
    properties
        
        % models - cell of initialized PBM model for each experiment
        % remains connected, matlab pass class always as a reference!
        % Pass a copy of the PBM class if you do not want to share reference
        % e.g. models={m.copy}
        models
        
        % cell of experimental data for each experiment
        expdata
        
        % lower bounds of variables, if neeeded
        lb
        % upper bounds of variables, if neeeded
        ub
        
        %
        
        % Method for objective function
        %   'generic' - use generic objfun
        %   'handle' - use objfun handle: err=@(x,models,expdata)
        objfun_method='generic';
        
        % Numerical method for solution of the optimization problem
        % particleswarm - needs lb,ub
        % patternsearch - needs x0,lb,ub
        % fminsearch - needs only x0
        % fmincon - needs x0,lb,ub
        % fminbnd - only 1D, needs lb,ub
        solver_method='particleswarm';
        
        % Method for error (performance) calculation
        % mse - mean squared error
        % mae - mean absolute error
        % msre - mean square root error
        error_method='mse';
        
        % objfun handle: err=@(x,models,expdata)
        % must accept (x,models,expdata) and return error 
        objfun_h
        
        % Number of points evaluated in each iteration (particleswarm method)
        swarm_size='100';
        
        % Initial guess for unconstrained optimization
        x0
        
        % Results of the solver
        xres
        
        % Iterations of the solver
        xiter
        
        % Iteration error of the solver
        erriter
        
        %% Settings for solver accuracy
        
        % Termination tolerance on x, a positive scalar.
        % The default value is 1e-4
        TolX=1e-4;
        
        % Termination tolerance on the function value, a positive scalar.
        % The default is 1e-4.
        TolFun=1e-4;
        
        % Maximum number of iterations allowed, a positive integer.
        % The default value is 200*numberOfVariables (assigned if empty)
        Maxiter=[]
        numberOfVariables
        
        %% Settings for generic objfun
        
        %prop_list - cell of property names to set in model
        %          - e.g. There are 3 properties in prop_list,
        %            first 3 values in x are used to set properties
        %              
        prop_list
        
        %handle_list - cell of handles names/properties to set in model
        %              handle will be called with all values of x
        %              can be used to set vector properties
        handle_list
        
        % handle_gen_list - cell of handle generators from x
        %
        handle_gen_list
        
        %target_list - target fields in PBM model and expdata
        %              format: 
        %              {model_field, [column subset], expdata_field, [column subset]}
        %              e.g. {'cm',[1],'cexp',[1]};
        %              e.g. {'cm',[1],'cexp',[1],'cstd',[1]};
        %              e.g. {'cm','cexp'} if both cm,cexp have same number of columns;
        %              multiple targets can be given in rows
        target_list
        
        % maximum time for integration
        time_max=1e9;
        
        % if display results
        log_to_screen=true;
        
        % if log results to file
        log_to_file=true;
        
    end
    
    methods
        function obj=OptimizePBM(obj)
            % Constructor of OptimizePBm class
        end
        
        function err=objfun_generic(obj,x)
            % Generic objective function
            %   x - optimized variables
            %      array(Nx,1)
            %   Needs following:
            %   models - cells of PBM models
            %          - cell(N,1)
            %
            %   expdata - cells of experimenta data
            %           each cell is a struct with fields, must contain texp            
            %           -cell(N,1)
            %   prop_list - cell of property names to set in model
            %             - e.g. There are 3 properties in prop_list,
            %               first 3 values in x are used to set properties
            %              
            %   handle_list - cell of handles names/properties to set in model
            %                 handle will be called with all values of x
            %                 can be used to set vector properties
            %   handle_gen_list - cell of handle generators from x
            %
            %   target_list - target fields in PBM model and expdata
            %               format: 
            %               {model_field, [column subset], expdata_field, [column subset]}
            %               e.g. {'cm',[1],'cexp',[1]};
            %               e.g. {'cm',[1],'cexp',[1],'cstd',[1]};
            %               e.g. {'cm','cexp'} if both cm,cexp have same number of columns;
            %               multiple targets can be given in rows
            
            N=length(obj.models);
            Nvals=size(obj.target_list,2);            
            if Nvals==1 || Nvals==5 || Nvals>6
                error('Property target_list is in wrong format.')
            end
            selcolumns=mod(Nvals,2)==0 & Nvals>3;
            k0=1;
            if selcolumns dk=2; else dk=1; end;
            
            err=0;
            for i=1:N
                m=obj.models{i};
                expi=obj.expdata{i};
               %% Set properties
                for j=1:length(obj.prop_list)
                   m.(obj.prop_list{j})=x(j);
                end
               %% Set handles or vector properties
                for j=1:length(obj.handle_list)
                    hgen=obj.handle_gen_list{j};
                    m.(obj.handle_list{j})=hgen(x);
                end
                m.w=m.w; % to update phi0dist if it would change
                m.set_initial_conditions(m.cstart,m.cs0,m.conc_units);
                sel=expi.texp<=obj.time_max;
                time=expi.texp(sel);
                time_unique=unique(time);                
                m.solve(time_unique);                    
                %% Calculate error                
                for j=1:size(obj.target_list,1)
                    k0=1;
                    ym_unique=obj.get_target(m,obj.target_list(j,k0:k0+dk-1));
                    k0=k0+dk;
                    yexp=obj.get_target(expi,obj.target_list(j,k0:k0+dk-1));
                    yexp=yexp(sel,:);
                    k0=k0+dk;                    
                    ym=interp1(time_unique,ym_unique,time);
                    ym=reshape(ym,size(ym));
                    if Nvals==3 || Nvals==6
                        errweight=obj.get_target(expi,obj.target_list(j,k0:k0+dk-1));                        
                        errweight=errweight(sel,:);
                        err=err+obj.(obj.error_method)(ym,yexp,errweight);
                    else
                        err=err+obj.(obj.error_method)(ym,yexp,1);
                    end                    
                end
            end
        end
        
        function x=get_target(obj,data,target)
            x=data.(target{1});
            if length(target)>1
                x=x(:,target{2});    
            end
        end
        
        function err=mse(obj,x,y,w)
           % Calculate error Mean Squered error  
           w(w==0)=1;
           err=sum((x-y).^2./w,'all')./prod(size(x));
        end
        
        function err=msre(obj,x,y,w)
           % Calculate error Mean Squered error  
           w(w==0)=1;
           err=sqrt(sum((x-y).^2./w,'all')./prod(size(x)));
        end
        
        function err=mae(obj,x,y,w)
           % Calculate error Mean absolute error (MAE)
           w(w==0)=1;
           err=sum(abs(x-y)./w,'all')./prod(size(x));
        end
                
        function res=optimize(obj,idname,logFileName)
            % idname - identificator string            
            % logFileName - path to file where opt results are written
            %        - use logFileName=[] to turn off logging to file
            
            obj.xiter=[];
            obj.erriter=[];
            if isempty(logFileName)
                obj.log_to_file=false;
            end
            % Find optimal model parameters
            if strcmp(obj.objfun_method,'generic')
                fun=@(x) obj.objfun_generic(x);
            elseif strcmp(obj.objfun_method,'handle')
                fun=@(x) obj.objfun_h(x,obj.models,obj.expdata);
            else
                error('Unknown value in objfun_method');
            end
            
            if strcmp(obj.solver_method,'particleswarm')
                res=obj.solver_particleswarm(fun,obj.lb,obj.ub,obj.swarm_size,idname,logFileName);
            elseif strcmp(obj.solver_method,'patternsearch')
                res=obj.solver_patternsearch(fun,obj.x0,obj.lb,obj.ub,idname,logFileName);
            elseif strcmp(obj.solver_method,'fminsearch')
                res=obj.solver_fminsearch(fun,obj.x0,idname,logFileName);
            elseif strcmp(obj.solver_method,'fmincon')
                res=obj.solver_fmincon(fun,obj.x0,obj.lb,obj.ub,idname,logFileName);
            elseif strcmp(obj.solver_method,'fminbnd')
                res=obj.solver_fminbnd(fun,obj.lb,obj.ub,idname,logFileName);
            else
                error('Unknown value in solver_method');
            end
            obj.xres=res;
        end
        
        function err=eval(obj,x)
            % Find optimal model parameters
            if strcmp(obj.objfun_method,'generic')
                fun=@(x) obj.objfun_generic(x);
            elseif strcmp(obj.objfun_method,'handle')
                fun=@(x) obj.objfun_h(x,obj.models,obj.expdata);
            else
                error('Unknown value in objfun_method');
            end
            err=fun(x);
        end
        
        %% Solvers
        function [ x ] = solver_fminsearch(obj,fun,x0,id,logFileName)
        %solver_fminsearch (obj,fun,x0,id,logFileName)
        %Find optimal values of x with Nelder-Mead method (do not use derivatives)
        %   
            obj.check_empty(fun,'fun')
            obj.check_empty(x0,'x0')
            obj.numberOfVariables=length(x0);
            if isempty(obj.Maxiter)
                obj.Maxiter=200*obj.numberOfVariables;
            end
            if obj.log_to_screen display_opt='iter'; else display_opt='none'; end;
            options = optimset('Display',display_opt,'OutputFcn',@obj.report_fminsearch, ...
                                    'TolX',obj.TolX,'TolFun',obj.TolFun,'MaxIter',obj.Maxiter);

            [x,fval,exitflag,output] = fminsearch(fun,obj.x0,options);
            obj.append_to_file(x,fval,id,logFileName);
        end
        
        function [ x ] = solver_particleswarm(obj,fun,lb,ub,swarm_size,id,logFileName)
            %solver_particleswarm 
            %Find optimal values of x wiht particleswarm algorith, append ressults to
            %the file
            obj.check_empty(fun,'fun')            
            obj.check_empty(lb,'lb')
            obj.check_empty(ub,'ub')                        
            obj.check_empty(swarm_size,'swarm_size')
            obj.numberOfVariables=length(lb);
            if isempty(obj.Maxiter)
                obj.Maxiter=200*obj.numberOfVariables;
            end
            
            if obj.log_to_screen display_opt='iter'; else display_opt='none'; end;
            options = optimoptions('particleswarm','UseParallel',true, ...
                    'SwarmSize',swarm_size,'HybridFcn',@patternsearch, ...
                    'Display',display_opt,'OutputFcn',@obj.reportx_ps);

            [x,fval,exitflag,output] = particleswarm(fun,length(lb),lb,ub,options);
            obj.append_to_file(x,fval,id,logFileName);
        end

        function [ x ] = solver_patternsearch(obj,fun,x0,lb,ub,id,logFileName)
        %solver_patternsearch 
        %Find optimal values of x wiht patternsearch algorith
        %
            obj.check_empty(fun,'fun')
            obj.check_empty(x0,'x0')
            obj.check_empty(lb,'lb')
            obj.check_empty(ub,'ub')    
            obj.numberOfVariables=length(x0);
            if isempty(obj.Maxiter)
                obj.Maxiter=200*obj.numberOfVariables;
            end
            
            if obj.log_to_screen display_opt='iter'; else display_opt='none'; end;
            options = optimoptions('patternsearch','UseParallel',true, ...
                    'Display',display_opt,'OutputFcn',@obj.report_patternsearch, ...
                    'TolX',obj.TolX,'TolFun',obj.TolFun,'MaxIter',obj.Maxiter);
            
            [x,fval] = patternsearch(fun,x0,[],[],[],[],lb,ub,[],options);
            obj.append_to_file(x,fval,id,logFileName);
        end
        
        function [ x ] = solver_fmincon(obj,fun,x0,lb,ub,id,logFileName)
        %solver_fmincon (obj,fun,x0,lb,ub,id,logFileName)
        % Find minimum of constrained nonlinear multivariable function (estimate derivatives)
        %   
            obj.check_empty(fun,'fun')
            obj.check_empty(x0,'x0')
            obj.check_empty(lb,'lb')
            obj.check_empty(ub,'ub')    
            obj.numberOfVariables=length(x0);
            if isempty(obj.Maxiter)
                obj.Maxiter=200*obj.numberOfVariables;
            end
            
            if obj.log_to_screen display_opt='iter'; else display_opt='none'; end;
            options = optimoptions('fmincon','UseParallel',true, ...
                    'Display',display_opt,'OutputFcn',@obj.report_fminsearch, ...
                    'StepTolerance',obj.TolX,'MaxIterations',obj.Maxiter);
            
            [x,fval] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
            obj.append_to_file(x,fval,id,logFileName);
        end
        
        function [ x ] = solver_fminbnd(obj,fun,lb,ub,id,logFileName)
        %solver_fminbnd (obj,fun,lb,ub,id,logFileName)
        % Find minimum in bounds, only for 1d problems
        %   
            obj.check_empty(fun,'fun')
            obj.check_empty(lb,'lb')
            obj.check_empty(ub,'ub')    
            obj.numberOfVariables=length(lb);
            if obj.numberOfVariables>1
                error('\n<strong> Solver fminbnd supports only 1D problems! </strong>\n')
            end
            if isempty(obj.Maxiter)
                obj.Maxiter=200
            end
            
            if obj.log_to_screen display_opt='iter'; else display_opt='none'; end;
            
            xa0=lb;
            a0=0;
            b0=ub-lb;
            seta1=true;
            setb1=true;
            p=1; % TODO: directional search 
            r=(3-sqrt(5))/2;
            tol=obj.TolX;
            for i=1:obj.Maxiter         
                da=b0-a0;
                if seta1
                    seta1=false;
                    a1=r*da+a0;
                    x1=xa0+a1.*p;
                    fx1=fun(x1);
                    if obj.log_to_screen
                        obj.print_res(x1,fx1);
                    end
                end
                if setb1
                    setb1=false;
                    b1=(1-r)*da+a0;
                    x2=xa0+b1.*p;
                    fx2=fun(x2);
                    if obj.log_to_screen
                        obj.print_res(x2,fx2);
                    end
                end
                if (b0-a0)<tol
                    break
                end
                if fx1<fx2
                    b0=b1;
                    b1=a1;
                    fx2=fx1;
                    seta1=true;
                else
                    a0=a1;
                    a1=b1;
                    fx1=fx2;
                    setb1=true;
                end                
            end
                               
            x=xa0+(a1+b1)/2.*p;
            fval=fun(x);
            if obj.log_to_screen
                 obj.print_res(x,fval);
            end
            obj.append_to_file(x,fval,id,logFileName);
        end
        
        %%Report function
        function stop = report_fminsearch(obj,x,optimvalues,state)
            % report function for patternsearch
            stop=false;
            print_res(obj,x,optimvalues.fval)
            obj.append_to_file(x,optimvalues.fval,[],'fminsearch_progress.log');
        end
        
        function [stop,options,optchanged] = report_patternsearch(obj,optimvalues,options,flag)
            % report function for patternsearch
            stop=false;
            optchanged=false;
            print_res(obj,optimvalues.x,optimvalues.fval)
            obj.append_to_file(optimvalues.x,optimvalues.fval,[],'patternsearch_progress.log');
        end
        
        function stop=reportx_ps(obj,optimValues,state)
            % report function for particleswarm
            stop=false;
            print_res(obj,optimValues.bestx,optimValues.bestfval)          
            obj.append_to_file(optimValues.bestx,optimValues.bestfval,[],'ps_progress.log');
        end
        
        function print_res(obj,x,fval)
            % print to output
            obj.xiter=[obj.xiter;x];
            obj.erriter=[obj.erriter;fval];
            if obj.log_to_screen
                fprintf('x: ')
                fprintf('%e, ',x)
                fprintf('fval: %e \n',fval)
            end
        end
        
        function append_to_file(obj,x,fval,id,logFileName)
            % append x,fval and if to a file
            if obj.log_to_file
                fileID = fopen(logFileName,'a');
                fprintf(fileID,'%e, ',x(:));
                fprintf(fileID,'\t fval: %e, ',fval);
                if ~isempty(id)
                    fprintf(fileID,'\t %s \n',id);
                end
                fclose(fileID);
            end
        end
        
        %% Checks
        function check_empty(obj,x,property_name)
            % Check if variable is empty
            if isempty(x)
                error(sprintf('\n<strong> Empty input variable/property %s! </strong>\n',property_name))
            end
        end        
    end
end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BSD 3-Clause License
% 
% Copyright (c) 2026, Jiri Kolar, Suada Djukaj
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

