function alldata = load_exp_fromfile(fname,fname_dist)
%load_exp_fromfile Load dissolution data stored in excel file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Inputs:
%      fname - excel file with experimental concentration, additional info
%      fname_dist - excel file with particle size distribution data
%   
% Output:
%      data - loaded experimental data

%   Expect the excel file 'fname' to contain any or multiple sheets with prefix:
%   rd_ for reactive dissolution
%   nd_ for normal dissolution
%   kin_ for experiments for reaction kinetics evaluation
%
%   Each sheet name of rd_ and nd_ prefix contains info of composition
%   e.g. C70M30 -> 70 wt. % of crude fraction and 30 wt. % of milled
%   fraction
%   for kin_ sheets  these are ignored

%   Optional is additional "info" sheet where other info can be supplied
%   if any field is not specified a default values is used.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distbasesheet='phi0';
infosheet='info';
sheets= sheetnames(fname);
info=struct;
if any(contains(sheets,'info'))
    inforaw=readcell(fname,'Sheet',infosheet);
    %Parse info    
    for i=1:size(inforaw,1)
        key=matlab.lang.makeValidName(inforaw{i,1});
        value=inforaw{i,2};
        info.(key)=value;
    end
end
info.w=[];
info.ns=4;
alldata=cell(0);
for i=1:length(sheets)
    expid=sheets(i);
    if contains(expid,'info')
        continue
    end
    ss=split(expid,'_');
    prefix=ss(1);
    expinfo=ss(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load concentration
    
    concData=readtable(fname,'Sheet',expid,'HeaderLines',1);
    keys=concData.Properties.VariableNames;
    if sum(contains(keys,'Time'))==1
        time_usedefault=true; %only one time column is supplied
    else
        time_usedefault=false; %time for each experiment is supplied
    end
    
    
    conc_PP=cell(0);
    conc_P=cell(0);
    time_all=cell(0);
    
    for j=1:length(keys)
        key=keys{j};     
        values=concData.(key);
        values=values(~isnan(values));
        if contains(key,'Time')
            time_all{length(time_all)+1}=values;
        elseif contains(key,'cPP')
            conc_PP{length(conc_PP)+1}=values;
        elseif contains(key,'cP')
            conc_P{length(conc_P)+1}=values;
        end        
    end
    
    if time_usedefault
        % use one time for all exps, all exps must be alligned in time
        texp=time_all{1};
        cexpPP=cell2mat(conc_PP);
        cexpP=cell2mat(conc_P);        
        cscatterPP=cexpPP;
        cscatterP=cexpP;
        tscatter=texp.*ones(size(cexpPP));
    else
        % different time for each experiment
        % resulting cexp and std is calculated by linear interpolation
        cscatterPP=[];
        cscatterP=[];
        tscatter=[];
        for j=1:length(conc_PP)
            tscatter=[tscatter; time_all{j}];
            cscatterPP=[cscatterPP; conc_PP{j}];
        end
        
        for j=1:length(conc_P)
            cscatterP=[cscatterP; conc_P{j}];
        end
        
        texp=unique(tscatter);
        Np=length(texp);
        cexpPP=zeros(Np,length(conc_PP));
        cexpP=zeros(Np,length(conc_P));
        for j=1:length(conc_PP);
            cexpPP(:,j)=interp1(time_all{j},conc_PP{j},texp,'linear');
        end
        
        for j=1:length(conc_P);
            cexpP(:,j)=interp1(time_all{j},conc_P{j},texp,'linear');
        end
    end
       


    cexpPP_std=std(cexpPP,[],2,'omitnan');
    cexpP_std=std(cexpP,[],2,'omitnan');
    cexpPP_mean=mean(cexpPP,2,'omitnan');
    cexpP_mean=mean(cexpP,2,'omitnan');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data=struct;
    data.all=struct;
    data.texp=texp;
    data.cexpPP_mean=cexpPP_mean;
    data.cexpP_mean=cexpP_mean;
    data.cexpPP_std=cexpPP_std;
    data.cexpP_std=cexpP_std;

    data.tscatter=tscatter;
    data.cscatterPP=cscatterPP;
    data.cscatterP=cscatterP;

    data.all.conc_P=conc_P;
    data.all.conc_PP=conc_PP;
    data.all.texp=time_all;

    % Mapping
    charids=['C';'M';'S';'P'];
    wids=[1,2,3,4];

    if ~contains(prefix,'kin')
        for j=1:length(wids)
            expinfo=replace(expinfo,charids(j),['_' charids(j)]);
        end
        expinfo=expinfo.split('_');
        w=zeros(1,length(wids));
        for j=1:length(expinfo)
           if isempty(expinfo{j})
               continue
           end
           expinfo0=expinfo{j};
           cid=expinfo0(1);
           w0=str2num(expinfo0(2:end));
           w(wids(charids==cid))=w0;
        end
        data.w=w./sum(w);
        data.sieve=[[0,1e9];[0,1e9];[0,1e9];[0,1e9]];

        % %For Horiba data
        % % dexp(0) = 0
        % % phi(i) is for this class dexp(i-1) < d < dexp(i)
        % % %each volume fraction i means it is for class less than dexp(i) and
        % % bigger than dexp(i-1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % data.all.texp=cell(1,1);
        % data.all.cexp=cell(1,1);
        % 
        % %Init dist
        data.phi0exp=cell(info.ns,1);
        data.all.dexp=cell(info.ns,1);

        for j=1:info.ns
            phi0raw=readcell(fname_dist,'Sheet',sprintf('phiw%d',j));
            ids=phi0raw(1,2:end);
            phi0all=cell2mat(phi0raw(2:end,2:end));    
            diamsi=cell2mat(phi0raw(2:end,1));
            phi0=phi0all(:,ismember(ids,sprintf('w%d',j)));
            %phi0=phi0(:,~ismember(1:info.nexp,info.exclude)); %exclude
            data.phi0exp{j}=mean(phi0,2)./sum(mean(phi0,2));   
            data.all.dexp{j}=diamsi*1e-6;
        end
        % %%%%%%%%%%%%%%%
        data.dexp=diamsi*1e-6;% [m]
        utils=Utils;
        data.dxexp=utils.get_dx_include0(data.dexp);
        data.xexp=utils.get_xclass_fromBounds([0;data.dexp]);
    end

    data.cB0=0;
    data.cPP0=0;
    data.cP0=0;
    data.isPP=true;
    data.isP=false;
    if ~contains(prefix,'nd')    
        data.isP=true;
    end
    if contains(prefix,'kin')
       data.cB0=info.cB0;
       data.cPP0=info.cPP0;
       data.cP0=info.cP0;
       data.cexpPP_mean(1)=data.cPP0;
    end
    data.expID=replace(expid,'_',' ');
    data.cunits='[umol/l]';
    data.tunits='[days]';
    alldata{length(alldata)+1}=data;
end
end

