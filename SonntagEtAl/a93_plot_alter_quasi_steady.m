%
% Concentration profiles for various PSDs at different hydrolysis rate
%
%       Run a9_quasi_steady_state.m before this script
%       or load SonntagEtAl/expdata/quasi_steady_data.mat

%load('SonntagEtAl/expdata/quasi_steady_data.mat')

%% Visualize
utils=Utils();

ax_all={}

for ik=1:Nk
fig=figure('Position',[100 100 800 600]);
ax=gca(fig);
cfg={'LineWidth',2};
fs=18;

%ax2=subplot(1,2,2);
Ncols=round(Nmu*1.5);
cols_all={winter(Ncols),autumn(Ncols),summer(Ncols)};
conc_umol_l=[M_PP,M_P,M_LiOH];
%cols=lines(Nexp);    
ls={'--',':','-.','--','-'};
larr=[];
for ic=1:Nc
    cols=cols_all{ic};
    for imu=1:Nmu        
    m=model_arr{imu,ic,ik};    
    if imu==length(mu_arr0)
        m.title=sprintf('1 + %.0f $\\mu$m (%.1f g/l)',mu_arr0(imu)*1e6,m.cs0);
    else
        m.title=sprintf('%.0f $\\mu$m (%.1f g/l)',mu_arr0(imu)*1e6,m.cs0);
    end
    csm=sum(m.fnm.*m.dx.*(m.x0.^3.*m.fV),1).*m.rhos;
    title0=m.title;
    
    l=plot(ax,m.tm,m.cm(:,1),ls{mod(imu-1,length(ls))+1},'DisplayName',title0,'Color',cols(imu,:),cfg{:});
    larr=[larr l];
    % Threshold
    cs_threshold=0.05;
    cs_target=max(m.cm(:,1))./(1+cs_threshold);
    cmax_target=cs_target*(1+cs_threshold);
    cmin_target=cs_target*(1-cs_threshold); 
    if imu==length(mu_arr0)
        plot(ax,[0 max(m.tm)],[cmax_target cmax_target],'--','Color',0.5*ones(1,3),cfg{:})
        lth=plot(ax,[0 max(m.tm)],[cmin_target cmin_target],'--','Color',0.5*ones(1,3),'DisplayName','$\pm 5$\% (1 + 100 $\mu$m)',cfg{:});
    end    
    %plot(ax2,m.tm,m.cm(:,2)./conc_umol_l(2),'DisplayName',title0,'Color',cols(iexp,:))
    hold(ax,'on')       
    %hold(ax2,'on')       
    end
end


xlim(ax,[0,max(m.tm)])
grid(ax,'on')
xlabel(ax,'Time [days]','Interpreter','latex')
ylabel(ax,'PP concentration [g/l]','Interpreter','latex')

if kmultiplier(ik)==1
    title(ax,'Hydrolysis rate 1x','Interpreter','latex')
else
    title(ax,sprintf('Hydrolysis slowed %dx',1/kmultiplier(ik)),'Interpreter','latex')
end
set(ax,'FontSize',fs,'TickLabelInterpreter','latex')
larr=[larr lth];
legend(ax,larr, 'Location','eastoutside', ... %'NumColumns',5,'Orientation','horizontal',
    'FontSize',fs,'Interpreter','latex')
utils.savefig_article(fig,sprintf('SonntagEtAl/results/quasi_steady_state_article_alt%d',ik));
end



