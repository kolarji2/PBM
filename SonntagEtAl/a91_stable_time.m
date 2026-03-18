%
%   Plot duration reached PP concentration and duration of stable release 
%   at different hydrolysis rate
%
%       Run a9_quasi_steady_state.m before this script
%       or load SonntagEtAl/expdata/quasi_steady_data.mat

%load('SonntagEtAl/expdata/quasi_steady_data.mat')

tstable_all=zeros(Nmu,Nc,Nk);
cstable_all=zeros(Nmu,Nc,Nk);
cbar_pos=zeros(Nmu,Nc,Nk);
cbar_neg=zeros(Nmu,Nc,Nk);
for imu=1:Nmu        
    for ik=1:Nk
    for ic=1:Nc
    m=model_arr{imu,ic,ik};
    cs_threshold=0.05;
    [tstable_all(imu,ic,ik),cstable_all(imu,ic,ik),cbar_pos(imu,ic,ik),cbar_neg(imu,ic,ik)]=time_stable(m,cs_threshold);
    end
    end
end


tmax_fig=[8 14 65];
leg_loc={'ne','e','se'};
icdraw={3,[1],[1]};
for ik=1:Nk
    cols=lines(4);
fig=figure('Position',[100 100 600 600]);
cfg={'LineWidth',2};
fs=24;
ax=gca(fig);
larr=[];
for ic=1:Nc
%         if imu==length(mu_arr0)
%             m.title=sprintf('1 + %.0f $\\mu$m (%.1f g/l)',mu_arr0(imu)*1e6,m.cs0);
%         else
%             m.title=sprintf('%.0f $\\mu$m (%.1f g/l)',mu_arr0(imu)*1e6,m.cs0);
%         end
        m=model_arr{1,ic,ik};
        title0=sprintf('%.1f g/l',m.cs0);
        pos=cbar_pos(:,ic,ik);
        neg=cbar_neg(:,ic,ik);
        l=errorbar(ax,tstable_all(:,ic,ik),cstable_all(:,ic,ik),neg,pos,'-o','Color',cols(ic,:),'DisplayName',title0,cfg{:});
        larr=[larr l];
        hold on
        if any(ic==icdraw{ik})
        for imu=1:Nmu
            m=model_arr{imu,ic,ik};            
            if imu==length(mu_arr0) & (ik>1 || (ik==1 && ic==3))
                txt0=sprintf('1 + %.0f $\\mu$m',mu_arr0(imu)*1e6);
                text(ax,tstable_all(imu,ic,ik),cstable_all(imu,ic,ik),txt0,'Interpreter','latex', ...
                'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',fs-2)
            elseif imu<length(mu_arr0)-1
                txt0=sprintf('%.0f $\\mu$m',mu_arr0(imu)*1e6);
                text(ax,tstable_all(imu,ic,ik),cstable_all(imu,ic,ik),txt0,'Interpreter','latex', ...
                'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',fs-2)
            end
            
            
        end
        end
end
xlabel(ax,'Stable release [days]','Interpreter','latex')
ylabel(ax,'Mean PP concentration [g/l]','Interpreter','latex')
%ylim([0,0.55]);
xlim([0,tmax_fig(ik)]);

%legend(ax,larr,'Location','eastoutside','FontSize',fs,'Interpreter','latex')
grid(ax,'on')
legend(ax,larr,'Location',leg_loc{ik},'FontSize',fs,'Interpreter','latex')
if kmultiplier(ik)==1
    title(ax,'Hydrolysis rate 1x','Interpreter','latex')
else
    title(ax,sprintf('Hydrolysis slowed %dx',1/kmultiplier(ik)),'Interpreter','latex')
end
set(ax,'FontSize',fs,'TickLabelInterpreter','latex')
utils.savefig_article(fig,sprintf('SonntagEtAl/results/stable_release_article_%d',ik));
end

%% Visualize
% utils=Utils();
% for ik=1:Nk

% ax=gca(fig);%subplot(1,2,1);
% %ax2=subplot(1,2,2);
% Ncols=round(Nmu*1.5);
% cols_all={winter(Ncols),autumn(Ncols),summer(Ncols)};
% conc_umol_l=[M_PP,M_P,M_LiOH];
% %cols=lines(Nexp);    
% ls={'--',':','-.','--','-'};
% larr=[];
% for ic=1:Nc
%     cols=cols_all{ic};
%     for imu=1:Nmu        
%     m=model_arr{imu,ic,ik};    
%     if imu==length(mu_arr0)
%         m.title=sprintf('1 + %.0f $\\mu$m (%.1f g/l)',mu_arr0(imu)*1e6,m.cs0);
%     else
%         m.title=sprintf('%.0f $\\mu$m (%.1f g/l)',mu_arr0(imu)*1e6,m.cs0);
%     end
%     csm=sum(m.fnm.*m.dx.*(m.x0.^3.*m.fV),1).*m.rhos;
%     title0=m.title;
%     
%     l=plot(ax,m.tm,m.cm(:,1),ls{mod(imu-1,length(ls))+1},'DisplayName',title0,'Color',cols(imu,:),cfg{:});
%     larr=[larr l];
%     % Threshold
%     cs_threshold=0.05;
%     cs_target=max(m.cm(:,1))./(1+cs_threshold);
%     cmax_target=cs_target*(1+cs_threshold);
%     cmin_target=cs_target*(1-cs_threshold); 
%     if imu==length(mu_arr0)
%         plot(ax,[0 max(m.tm)],[cmax_target cmax_target],'--','Color',0.5*ones(1,3),cfg{:})
%         lth=plot(ax,[0 max(m.tm)],[cmin_target cmin_target],'--','Color',0.5*ones(1,3),'DisplayName','$\pm 5$\%',cfg{:})
%     end    
%     %plot(ax2,m.tm,m.cm(:,2)./conc_umol_l(2),'DisplayName',title0,'Color',cols(iexp,:))
%     hold(ax,'on')       
%     %hold(ax2,'on')       
%     end
% end
% larr=[larr lth];



function [tstable,cstable,cbar_pos,cbar_neg]=time_stable(m,cs_threshold)
    tt=m.tm;
    cs_target=max(m.cm(:,1))./(1+cs_threshold);
    cmax_target=cs_target*(1+cs_threshold);   
    cmin_target=cs_target*(1-cs_threshold);   
    sel=tt>0;
    cm=m.cm(:,1);
    [cmax,cmaxi]=max(cm); 
    sel=tt>=tt(cmaxi) & cm>0;    
    sel2=tt<=tt(cmaxi);
    if cmin_target<cm(end)
        tmax=max(tt);
    else
        tmax=interp1(cm(sel),tt(sel),cmin_target);
    end        
    if sum(sel2)>1
        tmin=interp1(cm(sel2),tt(sel2),cmin_target);
    else
        tmin=0;
    end
    tstable=tmax-tmin;    
    sel=tt<=tmax & tt>=tmin;
    cstable=mean(cm(sel));
    cbar_pos=cmax_target-cstable;
    cbar_neg=cstable-cmin_target;
end