function [ cb ] = fzerobnd(obj,fun,lb,rb,tol)
%fzerobnd (obj,fun,lb,rb,tol): find root of 1-D function by interval halving
%   fun - function handle
%   lb - left bound
%   rb - right bound
%   tol - tolerance
    maxiter=1000;
    for i=1:maxiter
        flb=fun(lb);
        frb=fun(rb);
        if flb*frb>0
            error(sprintf('No root in the interval. Both values have same sign!\n f(%e)=%e\n f(%e)=%e',lb,flb,rb,frb));
        end
        cb=(lb+rb)/2;
        if abs(lb-rb)<tol
           break
        end
        fcb=fun(cb);
        if abs(fcb)<tol; return; end
        if fcb*flb<0
            rb=cb;
        else
            lb=cb;
        end    
    end
end