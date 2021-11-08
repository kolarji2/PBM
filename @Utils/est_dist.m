function [phi0_distest,dist0fun,mu1,sigma1] = est_dist(obj,x,dx,xexp,dxexp,phi0exp,method,sieve,dispfit)
%EST_DIST Estimate model distribution from experimental distribution using
%   selected method
%
%   Input:
%       x - cell centers (output dicretization)
%       dx - width of each cell
%       xexp - experimental cell centers
%       dxexp - width of each experimental cell
%       phi0exp - experimental distribution function (histogram)
%               - (integrated form sum(phi0exp)=1)
%       method - method used for estimation
%               * 3lognormal - fit using 3 lognormal distributions
%               * lognormal - fit using single lognormal distribution
%               * 3normal - fit using 3 normal distribution
%               * pchip - fit using Piecewise Cubic Hermite Interpolating Polynomial
%       sieve - min/max bounds for distribution, outside bounds function is set to 0
%       dispfit - display resulting fit
%               true/false
%   Output:
%       phi0_distest - resulting distribution function evaluated at cell centers
%       dist0fun - universal dist fun @(x,dx) can be evaluated at any points
%       mu1 - mean size of resulting dist
%       sigma1 - standard deviation of resulting dist
%
Np=length(x);
Npexp=length(xexp);

%resize
x=reshape(x,Np,1);
dx=reshape(dx,Np,1);
xexp=reshape(xexp,Npexp,1);
dxexp=reshape(dxexp,Npexp,1);
phi0exp=reshape(phi0exp,Npexp,1);

phi0_dist=phi0exp./dxexp;

% Bounds for distribution, outside bounds function is set to 0
minx=sieve(1);
maxx=sieve(2);

if ~strcmp(method,'pchip')
    npeaks=str2num(method(1));
    if isempty(npeaks)
        if strcmp(method,'lognormal')
            dist0 = @(x,dx,mu,sigma) obj.distfunLogNormal(x,dx,mu,sigma,minx,maxx); %Use cut normal dist
        elseif strcmp(method,'normal')
            dist0 = @(x,dx,mu,sigma) obj.distfunNormal(x,dx,mu,sigma,minx,maxx); %Use cut normal dist
        end
        objFun=@(a) sum((phi0_dist-dist0(xexp,dxexp,a(1),a(2))).^2);
        mu0=sum(phi0_dist.*dxexp.*xexp);
        sigma0=sqrt(sum((xexp-mu0).^2.*phi0_dist.*dxexp));
        a0=[mu0,sigma0];    
%         options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','none');
%         [res]=fsolve(objFun,a0,options);
        options= optimset('Display','none');
        [res]=fminsearch(objFun,a0,options);
        mu1=res(1);
        sigma1=res(2);
        dist0fun= @(x,dx) dist0(x,dx,mu1,sigma1);
        phi0_distest=dist0(x,dx,mu1,sigma1);
    else
        %multimodal dist
        method=method(2:end);
        if strcmp(method,'lognormal')
            dist0 = @(x,dx,mu,sigma) obj.distfunLogNormal(x,dx,mu,sigma,minx,maxx); %Use cut normal dist
        elseif strcmp(method,'normal')
            dist0 = @(x,dx,mu,sigma) obj.distfunNormal(x,dx,mu,sigma,minx,maxx); %Use cut normal dist
        end
        %distfunMultiModal(x,dx,w,mu,sigma,dist0)
        distMulti=@(x,dx,a) obj.distfunMultiModal(x,dx,a(1:npeaks),a((npeaks+1):2*npeaks),a((2*npeaks+1):3*npeaks),dist0);
        
        objFun=@(a) 1e-12*sum((phi0_dist-distMulti(xexp,dxexp,a)).^2)+sum((phi0_dist.*dxexp-distMulti(xexp,dxexp,a).*dxexp).^2);
        a0=zeros(1,3*npeaks);
        %w0=1/npeaks;
        %mu0=sum(phi0_dist.*dxexp.*xexp);
        %sigma0=sqrt(sum((xexp-mu0).^2.*phi0_dist.*dxexp));
        [w0,locs,widths]=findpeaks(phi0_dist.*dxexp,xexp);
        
        [w0,sorti]=sort(w0','descend');
        mu0=locs(sorti)';        
        sigma0=widths(sorti)';        
        if length(mu0)<npeaks
            mu0=repmat(mu0,1,npeaks);
            w0=repmat(w0,1,npeaks);
            sigma0=repmat(sigma0,1,npeaks);
        end
        sigma0=sigma0(1:npeaks);
        w0=w0(1:npeaks);
        a0(1:npeaks)=w0./sum(w0);
        a0((npeaks+1):2*npeaks)=mu0(1:npeaks); % mu0;
        a0((2*npeaks+1):3*npeaks)=sigma0./2;
        %a0((npeaks+1):2*npeaks)=a0((npeaks+1):2*npeaks).*linspace(0.7,1.2,npeaks);        
        %options = optimoptions('fsolve','Display','iter'); %,'Algorithm','levenberg-marquardt'
        options= optimset('Display','none');
        [res]=fminsearch(objFun,a0,options);        
        dist0fun= @(x,dx) distMulti(x,dx,res);        
        phi0_distest=dist0fun(x,dx);
    end
else
    dist0fun= @(x,dx) pchip(xexp,phi0_dist,x)./sum(pchip(xexp,phi0_dist,x).*dx);
    phi0_distest=dist0fun(x,dx);    
end

mu1=sum(phi0_distest.*dx.*x);
sigma1=sum((x-mu1).^2.*phi0_distest.*dx);

if dispfit
    figure
    subplot(2,1,1)
    semilogx(xexp,phi0_dist,'o--r')
    hold on
    %semilogx(xexp,dist0(xexp,dxexp,mu0,sigma0),'-g')
    txt_=sprintf('mu= %.3e \nsigma= %.3e\n',mu1,sigma1);
    text(1e-6,mean(phi0_distest),txt_);
    semilogx(x,phi0_distest,'-b')    
    legend({'exp';'final fit'}) %'guess';
    
    subplot(2,1,2)
    phi0hist=interp1(x,phi0_distest,xexp).*dxexp;
    semilogx(xexp,phi0exp,'o--r')
    hold on
    semilogx(xexp,phi0hist,'-b')  
    legend({'exp';'final fit'}) %'guess';
end
end

