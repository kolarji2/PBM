function err = objfun_composition(w,models,expdata)
    [penalty,w] = rosenbrock_penalty(w);
    %OBJFUN_CONC Calculate error based on composition
    time_sol=linspace(0,600,201);    
    err=0;
    for iexp=1:length(models)
        m=models{iexp};
        datai=expdata{iexp};
        m.w=w;
        m.set_initial_conditions(m.cstart,m.cs0,m.conc_units);
        [tt,c]=m.solve(time_sol);        
        sel=datai.texp<600;
        cmodel=interp1(tt,c,datai.texp(sel)');
        dc=(cmodel-datai.cexp(sel));
        %abs error
        err=err+sum(abs(dc))./length(dc)*penalty;
    end
end

function [penalty,w] = rosenbrock_penalty(w)
%ROSENBROCK_PENALTY rosenbrock penalty function 
% For condition all(w>=0) && all(w<=1)

A=1e-9;
scale=10;
penalty=1;

for i=1:length(w)
    penalty=penalty*psi_penalty(w(i),0-A,1+A,A,scale);
end
w(w<0)=0;
w(w>1)=1;

end

function penalty = psi_penalty(x,a,b,A,scale)
%PSI_PENALTY Rosenbrock penalty function for minimizing problems

if x<a
    penalty=2*scale;
elseif x>a && x<a+A
    penalty=1 + ((a+A-x)/A)^2*scale;
elseif a+A<=x && x<=b-A
    penalty=1;    
elseif x>b-A && x<b
    penalty=1 + ((x-b+A)/A)^2*scale;
else
    penalty=2*scale;
end
end
