%
%   Full version for optimization of composition
%   Generate graphs as in the article
%   Inverse problem solution
%
%
% Directories for results
results_dir='DjukajEtAl/results';
%% Control variables
calculate_error=false; % calculate error and save results, can take long time
plot_bimodal=true;
plot_trimodal=true;
%% Load experimental data
load_datasetKCl

%% Init model
Nexp=length(dataKCl);
model_arr=cell(length(dataKCl),1);
names={'monomodal','bimodal','trimodal'};

for iexp=1:Nexp
    %% Create new model instance
    m=KClPBM(infoKCl{iexp},dataKCl{iexp});
    m.set_optimized_pars();
    m.ode_flux_limiter='Koren';
    m.solver_name='dissolve_qamar_fxfvm';
    m.ode_solver_h=@ode23t;
    m.title=names{iexp};  
    %% Save
    model_arr{iexp}=m;
end

%% Screening for error
if calculate_error
    data_bimodal=struct;
    data_trimodal=struct;
    search_bimodal=struct;
    search_trimodal=struct;
    
    
    Nwbimodal=0 %101
    Nwtrimodal=0 %101 %101

    for iexp=2:3
        % Initialize optimization class
        opt=OptimizePBM;
        m=model_arr{iexp};
        opt.models={m.copy()};
        opt.expdata=dataKCl(iexp);
        opt.solver_method='fminsearch';  
        opt.objfun_method='handle';        
        opt.objfun_h=@(x,models,expdata) objfun_composition([x,1-sum(x)],models,expdata);
        opt.log_to_file=false;
        
        % Bimodal
        if iexp==2        
            Nw=Nwbimodal;
            warr=linspace(0,1,Nw);    
            wdata=zeros(Nw,2);
            errdata=zeros(Nw,1);
            cdata=cell(Nw,1);
            tdata=cell(Nw,1);  
            tic
            parfor i=1:Nw
                w=[warr(i),1-warr(i)];
                opt0=opt.copy(); % to make sure that each thread has own copy
                err = opt0.eval(w(1:end-1));
                errdata(i)=err;
                wdata(i,:)=w;
                cdata{i}=opt0.models{1}.cm;
                tdata{i}=opt0.models{1}.tm;      
            end
            fprintf('Bimodal completed in %f min\n',toc/60);
            % Save results
            data_bimodal.w=wdata;
            data_bimodal.err=errdata;
            data_bimodal.c=cdata;
            data_bimodal.t=tdata;    
            [errmin,imin]=min(errdata);
            data_bimodal.errmin=errmin;
            data_bimodal.wmin=warr(imin);    
            
            % Calculate optimal error
            opt.x0=[0.8];
            opt.objfun_h=@(x,models,expdata) objfun_composition([x,1-sum(x)],models,expdata);  
            xres=opt.optimize('','');
            search_bimodal.xres=xres;
            search_bimodal.warr=opt.xiter;
            search_bimodal.err=opt.erriter;
            
        else
        %%%%%%%%%%%%%%%%%%%%%%%% Trimodal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            Nw=Nwtrimodal;
            [w1_m,w2_m]=meshgrid(linspace(0,1,Nw),linspace(0,1,Nw));
            errdata=zeros(Nw,Nw);
            wdata=zeros(Nw,Nw,3);
            cdata=cell(Nw,Nw);
            tdata=cell(Nw,Nw);  
            ttotal=0;
            for i=1:Nw
                fprintf('Trimodal completed: %.2f %%\n',(i-1)/Nw*100)
                tic;
                parfor j=1:Nw                    
                    w=[w1_m(i,j),w2_m(i,j),0];
                    w(3)=1-sum(w);                    
                    if w(3)<0
                        continue
                    end
                    opt0=opt.copy(); % to make sure that each thread has own copy
                    err = opt0.eval(w(1:end-1));                    
                    wdata(i,j,:)=w;
                    errdata(i,j)=err;
                    cdata{i,j}=opt0.models{1}.cm;
                    tdata{i,j}=opt0.models{1}.tm;
                end 
                ttotal=ttotal+toc;
            end
            fprintf('Trimodal completed in %f min\n',ttotal/60);
            % Save results
            data_trimodal.err=errdata;
            data_trimodal.w=wdata;
            data_trimodal.c=cdata;
            data_trimodal.t=tdata;
            minVal=min(errdata,[],'all');
            [imin,jmin]=find(errdata==minVal);
            data_trimodal.errmin=minVal;
            data_trimodal.wmin=[w1_m(imin,jmin),w2_m(imin,jmin),1-w1_m(imin,jmin)-w2_m(imin,jmin)];
            
            % Calculate optimal error
            opt.x0=[0.3,0.7];              
            xres=opt.optimize('','');                   
            search_trimodal.xres=xres;
            search_trimodal.warr=opt.xiter;
            search_trimodal.err=opt.erriter;
        end
    end
    save('DjukajEtAl/results_composition.mat','data_bimodal','data_trimodal');
    save('DjukajEtAl/results_composition_search.mat','search_bimodal','search_trimodal');
else
    dts_comp=load('DjukajEtAl/results_composition.mat');
    data_bimodal=dts_comp.data_bimodal;
    data_trimodal=dts_comp.data_trimodal;

    dts_search=load('DjukajEtAl/results_composition_search.mat');
    search_bimodal=dts_search.search_bimodal;
    search_trimodal=dts_search.search_trimodal;
end

%% Visualize
init_colors_articleDjukajEtAl;



plotcfg={'LineWidth',2.0};
fs={'FontSize',20};
fstics={'FontSize',18};

if plot_bimodal
    infoi=infoKCl{2};
    fig=figure('Position',[100 100 800 600]);
    ax=gca(fig);
    
    Nw=size(data_bimodal.w,1);    
    
    for i=1:5:Nw
        w=data_bimodal.w(i,:);
            
       tt=data_bimodal.t{i};
       c=data_bimodal.c{i};
       
       if w(1)==0 || w(1)==1
            ls='-';
            coli=red;
            ilimit=i;
            continue
        elseif w(1)==0.5
            ls='-';
            coli=green;
            continue
        else
            ls='--';
            coli=colgray(2,:);
        end    
       
       l2=plot(ax,tt,c,ls,'color',coli,'LineWidth',1.5,'DisplayName','Combinations');    
       hold(ax,'on')       
    end
    
    tt=data_bimodal.t{ilimit};
    c=data_bimodal.c{ilimit};
    plot(ax,tt,c,'-','color',red,'LineWidth',1.5,'DisplayName','Limiting cases');
    tt=data_bimodal.t{1};
    c=data_bimodal.c{1};
    l0=plot(ax,tt,c,'-','color',red,'LineWidth',1.5,'DisplayName','Limiting cases');
    l1=plot(ax,dataKCl{2}.texp,dataKCl{2}.cexp,'o','color',green,'LineWidth',2.5,'DisplayName','Exp.');
    xlabel(ax,'Time [s]',fs{:})
    ylabel(ax,'{\it C}_b [kg m^{-3}]',fs{:})   
    xlim([0,600]);
    set(ax,'LineWidth',1.5,fstics{:});  
    legend(ax,[l0,l1,l2],'Location','southeast',fs{:})
    grid(ax,'on'); 
    m.savefig_article(fig,fullfile(results_dir,'comp_bimodal_limits'))
    
    % Error composition
    fig=figure('Position',[100 100 800 600]);
    w=linspace(0.4,0.6,101);
    fun =@(w) interp1(data_bimodal.w(:,1),data_bimodal.err,w,'cubic');
    ax=gca(fig);
    l1=plot(ax,w,fun(w),'-','color',blue,'DisplayName','Bimodal fit error',plotcfg{:});
    hold(ax,'on');
    errmax=max(data_bimodal.err);
    l2=plot(ax,[0.5 0.5],[0 errmax],'DisplayName','Target composition','color',green,plotcfg{:});
    wmin=search_bimodal.xres; %fminsearch(fun,0.5); 
    wt=[0.5, 0.5];
    wmin_r=round(wmin,4);
    wp=[wmin_r,1-wmin_r];    
    fprintf('%f %f \n',wp);
    fprintf('Aprox error: %f %f [%%]\n',(wp-wt)./wt*100)
    %search way
    l3=plot(ax,search_bimodal.warr,search_bimodal.err,'--o','DisplayName','Search progress','color',red,plotcfg{:});
    plot(ax,search_bimodal.warr(end),search_bimodal.err(end),'x','color',red,'LineWidth',3.5,'MarkerSize',12);
    l4=plot(ax,[wmin wmin],[0 search_bimodal.err(end)],'-','DisplayName','Predicted composition','color',red,plotcfg{:});
    %plot(ax,[0.4 0.6],[stderror stderror],'--','DisplayName','Threshold','color',colgray(1,:),plotcfg{:})
    xlabel(ax,'w_1 [-]',fs{:});
    ylabel(ax,'MAE',fs{:});
    set(ax,'LineWidth',1.5,fstics{:},'GridLineStyle','-'); 
    legend(ax,[l1,l2,l3,l4],'Location','north',fs{:});
    xlim([0.4,0.6])
    ylim([0,0.5])
    xarrow=[0.6 0.554];
    yarrow=[0.23 0.145];
    annotation('textarrow',xarrow,yarrow,'String',sprintf('%.4f',wmin),'FontSize',20,'Units','normalized');
    %xticklabels(ax,arrayfun(@(a) sprintf('%.2f',a),xticks_,'UniformOutput',false))  
    grid(ax,'on');
    m.savefig_article(fig,fullfile(results_dir,'comp_bimodal'))
end

if plot_trimodal
    %%%%%%%%%%%%%%%%%%%%%%%% Trimodal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    infoi=infoKCl{3};
    %%%%%%%% limits
    fig=figure('Position',[100 100 800 600]);
    ax=gca(fig);
    
    Nw=size(data_trimodal.w,1);    
    
    
    for i=1:10:Nw
        for j=1:10:Nw
            if isempty(data_trimodal.c{i,j})
                continue
            end
            w=data_trimodal.w(i,j,:);            

           tt=data_trimodal.t{i,j};
           c=data_trimodal.c{i,j};

           if w(1)==1 || w(3)==1
                ls='-';
                coli=red;
                ilimit=[i,j];
                continue
            else
                ls='--';
                coli=colgray(2,:);
            end   

           l2=plot(ax,tt,c,ls,'color',coli,'LineWidth',1.5,'DisplayName','Combinations');  
           hold(ax,'on')       
        end
    end
    
    tt=data_trimodal.t{ilimit(1),ilimit(2)};
    c=data_trimodal.c{ilimit(1),ilimit(2)};
    l0=plot(ax,tt,c,'-','color',red,plotcfg{:},'DisplayName','Limiting cases');
    tt=data_trimodal.t{1,1};
    c=data_trimodal.c{1,1};
    l0=plot(ax,tt,c,'-','color',red,plotcfg{:},'DisplayName','Limiting cases');
    l1=plot(ax,dataKCl{3}.texp,dataKCl{3}.cexp,'o','color',green,'LineWidth',2.5,'DisplayName','Exp.');
    xlabel(ax,'Time [s]',fs{:})
    ylabel(ax,'{\it C}_b [kg m^{-3}]',fs{:})   
    xlim([0,600]);
    set(ax,'LineWidth',1.5,fstics{:});  
    legend(ax,[l0,l1,l2],'Location','southeast',fs{:})
    grid(ax,'on');
    m.savefig_article(fig,fullfile(results_dir,'comp_trimodal_limits'))

    %%%%%%%% composition diagram
    Nw=size(data_trimodal.w,1);
    [w1_m,w2_m]=meshgrid(linspace(0,1,Nw),linspace(0,1,Nw));
    
    fig=figure('Position',[100 100 800 600]);
    ax=gca(fig);
    levels=[0.07,0.1,0.15,0.2,0.25,0.3,0.4]; %linspace(0.1,0.5,17);
    
    %[cb2,h2]=contourf(ax,w1_m,w2_m,data_trimodal.err,[stderror],'DisplayName','Trimodal fit error',plotcfg{:});
    [cb,h]=contour(ax,w1_m,w2_m,data_trimodal.err,levels,'DisplayName','Trimodal fit error',plotcfg{:});
    clabel(cb,h);
    hold(ax,'on');
    %plot(infoi.w(1),infoi.w(2),'or');
    colormap('parula');
    hbar=colorbar(ax);
    set(get(hbar,'label'),'string','MAE');
    
    plot(ax,infoi.w(1),infoi.w(2),'x','DisplayName','Target composition','color',green,'LineWidth',5,'MarkerSize',12);
    fun=@(w) interp2(w1_m,w2_m,data_trimodal.err,w(1),w(2),'cubic');
    %wmin=fminsearch(fun,[0.16,0.3]) 
    wmin=search_trimodal.xres;
    wt=round([infoi.w(1),infoi.w(2),1-infoi.w(1)-infoi.w(2)],4);
    wmin_r=round(wmin,4);
    wp=[wmin_r(1),wmin_r(2),1-sum(wmin_r)];
    fprintf('%f %f %f \n',wp);
    fprintf('Aprox error: %f %f %f [%%]\n',(wp-wt)./wt*100)
    
    plot(ax,search_trimodal.warr(:,1),search_trimodal.warr(:,2),'--o','DisplayName','Search progress','color',red,plotcfg{:})
    
    plot(ax,wmin(1),wmin(2),'x','DisplayName','Predicted composition','color',red,'LineWidth',3.5,'MarkerSize',12);
    ylim([0.2,0.6]);
    xlim([0,0.4]);    
    xlabel(ax,'w_1 [-]',fs{:});
    ylabel(ax,'w_2 [-]',fs{:});
    set(ax,'LineWidth',1.5,fstics{:});  
    legend(ax,'Location','southwest',fs{:});    
    grid(ax,'on');
    m.savefig_article(fig,fullfile(results_dir,'comp_trimodal'))


end
