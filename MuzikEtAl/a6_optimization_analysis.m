%
%   Example of usage of modified PBM class for 
%		modelling of pH dependent dissolution of telmisartan.
%		
%		Modified class can be found in directory: @PBMTelmi 
%   
% Compares reaction time and time required to deposit solution onto particle
% carriers:
%  - constant spray rate 5, 10, 20 or 100 ml/min
%  - constant required mass of deposited Telmisartan - 1 kg
%  - variable concentration of Telmisartan solution
%  
%  Concentrated solutions are sprayed faster, but it takes longer to
%  prepare them



%% Calculation
% ka = 5e-4
data_normal=load('MuzikEtAl/results/time_end_ka_5.mat');

% ka = 20e-4
data_enhanced=load('MuzikEtAl/results/time_end_ka_20.mat');
data_arr={data_normal,data_enhanced};
flabel_arr={'ka=5e-4','ka=20e-4'};

fs=12;
cfg={'FontSize',fs};
cfgline={'LineWidth',1.5};

fig=figure('Position',[100,100,1200,800]);
Vdot_arr=[5,10,20,100]; %ml/min
for i=1:2
    ax=subplot(1,2,i);

    data=data_arr{i};
    flabel=flabel_arr{i};
    mspray=1; %kg of product sprayed
    for j=1:4
        Vdot=Vdot_arr(j); %ml/min
       
        
        cm=data.cm_arr(:,1);
        cn=data.cn_arr(:,1);
        tend=data.tend_arr;
        
        
        cm_prod=sum(data.cm_arr,2);
        Vspray=mspray./cm; %m3
        tspray=Vspray./(Vdot./60*1e-6);
        
        tend_h=tend./3600;
        tspray_h=tspray./3600;
        
        
        ttotal=tend_h+tspray_h;
        
        
        tr_fun=@(c) interp1(cm,tend_h,c);
        tsp_fun=@(c) interp1(cm,tspray_h,c);
        
        tmin_fun=@(c) interp1(cm,ttotal,c);
        
        cmin=fminbnd(tmin_fun,1,300);
        
        %% Plot
        
        %Total
        plot(ax,cm,ttotal,cfgline{:},'DisplayName',sprintf('%d ml/min',Vdot))
        hold on
        tmin_arr(j)=tmin_fun(cmin);
        cmin_arr(j)=cmin;
    end
    plot(ax,cmin_arr,tmin_arr,'-or',cfgline{:},'DisplayName','Min')
    grid on
    legend(ax,cfg{:})
    xlabel(ax,'cTelmi [g/l]',cfg{:})
    ylabel(ax,'Process time [h]',cfg{:})
    ylim(ax,[0,250])
    title(ax,sprintf('Time comparison (%s)',flabel))

    
    % mdata=struct;
    % mdata.cnT_moll=cn;
    % mdata.cmT_gl=cm;
    % mdata.cm_prod_gl=cm_prod;
    % mdata.t_reaction_h=tend_h;
    % mdata.t_spray_h=tspray_h;
    % mdata.t_total_h=ttotal;
    % mdata.Vtotal_l=Vspray*1e3;
    % mdata.cost_fun=ttotal_w;
    % 
    % fields={{'cnT_moll','cmT_gl','cm_prod_gl','t_reaction_h','t_spray_h','t_total_h','Vtotal_l','cost_fun'}};
    % fname='MuzikEtAl/results/data_opt_time_reaction.xlsx';
    % sheet=sprintf('opt_%s',flabel);
    % m.save_to_excel({mdata},fields,fname,sheet);   
end

m=PBMTelmi();
m.savefig_article(fig,sprintf('MuzikEtAl/results/opt_time_reaction_%s',flabel))