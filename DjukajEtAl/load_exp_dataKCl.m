function [data,info] = load_exp_dataKCl(input_file,cond2conc)
%LOAD_EXP_DATA Load dissolution data stored in excel file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Excel must contain following sheets:
%   info (exclude,nexp,nphi,ns,V,phiipa,phiw,mtot,w1..wn )
%   phi1 ... phin
%   conductivity
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
condsheet='conductivity';
distbasesheet='phi';
infosheet='info';

inforaw=readcell(input_file,'Sheet',infosheet);
%Parse info
info=struct;
w=[];
sieve=[];
for i=1:size(inforaw,1)
    key=inforaw{i,1};
    if ismissing(key)
        continue
    end
    value=inforaw{i,2};
    if strcmp(key(1),'w')
       w=[w inforaw{i,2}];
       sieve=[sieve; inforaw{i,3},inforaw{i,4}];
    elseif contains(key,'exclude')
        info.(key)=str2num(value);
    else
        info.(key)=value;
    end
end
info.w=w;
info.sieve=sieve;

Nphi=info.nphi;
Nexp=info.nexp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conductivityData=readmatrix(input_file,'Sheet',condsheet);
time_cond=conductivityData(:,1);
conductivityAll=conductivityData(:,2:end);
conductivity=conductivityAll(:,~ismember(1:Nexp,info.exclude));
conc = cond2conc(conductivity);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   TODO: find better solution !!!
%   Maybe substract it instead? For monomodal negative values are for 0 time
%
if ~ismember(0,time_cond)
   error('Time zero not in conductivity data!\n') ;
end
conc(1,:)=0; % force zero concentration at time 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmean=mean(conc,2);
cstd=std(conc,[],2);
%pars = getpars_KCl();
%cmean(cmean>pars.csat)=pars.csat; % Concentration can not be above saturated conc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %For Horiba data
% % dexp(0) = 0
% % phi(i) is for this class dexp(i-1) < d < dexp(i)
% % %each volume fraction i means it is for class less than dexp(i) and
% % bigger than dexp(i-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=struct;
data.all=struct;
data.all.phiexp=cell(1,1);
data.all.texp=cell(1,1);
data.all.cexp=cell(1,1);
data.all.dexp=cell(1,1);

k=1;
for i=1:Nphi
    if ismember(i,info.exclude)
        continue
    end
    datai=readmatrix(input_file,'Sheet',[distbasesheet sprintf('%d',i)]);
    diamsi=datai(2:end,1);
    phisAll=datai(2:end,2:end);
    timeAll=datai(1,2:end);
    timei=unique(timeAll);
    Ntimes=length(timei);
    Ndiams=length(diamsi);
    phii=zeros(Ndiams,Ntimes);
    for j=1:Ntimes
        sel=ismember(timeAll,timei(j));
        phii(:,j)=mean(phisAll(:,sel),2);
    end
    phii=phii./sum(phii,1);  
    if i==1 %first experiment, might be used average of all, but interpolation is necessary (different time for each exp)
        sel=timei<min(1200,max(time_cond))+11; %601;
        data.dexp=diamsi*1e-6;% [m]
        data.phiexp=phii(:,sel); %[-]
        data.texp=timei(sel);% [s]
        
        data.cond=interp1(time_cond,mean(conductivity,2),data.texp,'makima','extrap')';
        data.condstd=interp1(time_cond,std(conductivity,[],2),data.texp,'makima','extrap')';
        data.cexp=interp1(time_cond,cmean,data.texp,'makima','extrap')';
        data.cstd=interp1(time_cond,cstd,data.texp,'makima','extrap')';
        data.cexpAll=cmean';
        data.tcexpAll=time_cond';
        data.cstdAll=cstd';
    end
    data.all.phiexp{k}=phii;
    data.all.texp{k}=timei;
    data.all.cexp{k}=interp1(time_cond,conc(:,k),timei)';
    
    k=k+1;
end

%Init dist
phi0raw=readcell(input_file,'Sheet','phi0');
ids=phi0raw(1,2:end);
phi0all=cell2mat(phi0raw(2:end,2:end));
diamsi=cell2mat(phi0raw(2:end,1));
data.phi0exp=cell(info.ns,1);
for i=1:info.ns
    data.all.dexp{i}=diamsi*1e-6;
    phi0=phi0all(:,ismember(ids,sprintf('w%d',i)));
    %phi0=phi0(:,~ismember(1:info.nexp,info.exclude)); %exclude
    data.phi0exp{i}=mean(phi0,2)./sum(mean(phi0,2));    
end

utils=Utils;
data.dxexp=utils.get_dx_include0(data.dexp);
data.xexp=utils.get_xclass_fromBounds([0;data.dexp]);
data.bounds=[0; data.dexp];
end

