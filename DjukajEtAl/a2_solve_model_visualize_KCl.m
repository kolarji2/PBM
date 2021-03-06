%
%   Example of usage of PBM class for modelling of KCl dissolution
%   
%% Description
%   Solve model and plot model results vs experimental data

%% Load experimental data
load_datasetKCl

%% Init model
Nexp=length(dataKCl);
model_arr=cell(length(dataKCl),1);
names={'monomodal','bimodal','trimodal'};
tic
for iexp=1:Nexp
    %% Create new model instance
    m=KClPBM(infoKCl{iexp},dataKCl{iexp});
    m.set_optimized_pars();
    m.title=names{iexp};
    m.check_settings;
    %% Solve model
    time=[0:1:600];
    [time,c,fn,phi_dist] = m.solve(time);    
    %% Save
    model_arr{iexp}=m;
end
toc
%% Visualize
%%% What to plot %%%%
plot_dists=true;
plot_conc=true;
plot_init_dists=false;
plot_exp_dists=false;
%%%%%%%%%%%%%%%%%%%%%%

init_colors_articleDjukajEtAl

% Directories for results
results_dir='DjukajEtAl/results';

% Plot concentration
if plot_conc
    % Font sizes, linewidth for conc plot
    plotcfg={'LineWidth',2.5};
    fs={'FontSize',20};
    fstics={'FontSize',18};    
    
    % plot t conc
    fig=figure('Position',[100 100 1000 600]);
    ax=gca(fig);

    names_fig={'Monomodal','Bimodal','Trimodal'};

    for iexp=1:3
        datai=dataKCl{iexp};
        leg1=sprintf('%s %s',names_fig{iexp},'Model');
        leg2=sprintf('%s %s',names_fig{iexp},'Exp.');
        errorbar(ax,datai.tcexpAll,datai.cexpAll,datai.cstdAll,'o','Color',colm(iexp,:),'DisplayName',leg2,plotcfg{:});
        hold(ax,'on')
        plot(ax,model_arr{iexp}.tm,model_arr{iexp}.cm,'-','Color',colm(iexp,:),'DisplayName',leg1,plotcfg{:}) 
    end
    xlabel(ax,'Time [s]',fs{:})
    ylabel(ax,'{\it C}_b [kg m^{-3}]',fs{:})
    xlim([0,600])
    legend(ax,'Location','southeast',fs{:})
    set(ax,'LineWidth',1.5,fstics{:}); 
    grid(ax,'on');
    model_arr{1}.savefig_article(fig,fullfile(results_dir,'conc'));
end

% Plot distribution
if plot_dists
    lw=1.0;
    plotcfg={'LineWidth',1.5};
    fs={'FontSize',12};
    fstics={'FontSize',10};    
    names={'monomodal','bimodal','trimodal'};
    % Plot distribution only at these times
    time_sel={};
    time_sel{1}=[0,30,50,150,250,350,450,550];
    time_sel{2}=[0,10,30,50,70,90,150,550];
    time_sel{3}=[0,10,30,50,70,150,250,530]; % 550 was only for 2 dist

        
    for iexp=1:3
        % Load model results
        m=model_arr{iexp};
        
        fig=figure('Position',[100 100 2000 1000]);
        time_seli=time_sel{iexp};
        for i=1:8
            ax=subplot(2,4,i);
            xexp=dataKCl{iexp}.xexp;
            ttarget=time_seli(i);
            l=[];
            k=1;            
            for j=1:length(dataKCl{iexp}.all.phiexp)
               phiexp=dataKCl{iexp}.all.phiexp{j};
               texp=dataKCl{iexp}.all.texp{j};
               sel=texp==ttarget;
               phiexp=phiexp(:,sel)*100;
               lege=sprintf('Exp %d (%d s)',j,ttarget);
               l(k)=semilogx(ax,xexp*1e6,phiexp,'--o','Color',colgray(j,:), ...
                                            'DisplayName',lege,plotcfg{:});
               hold(ax,'on')
               k=k+1;
               % Peaks
               [pkse,locse,w2,p2] = findpeaks(phiexp,xexp.*1e6);               
               plot(ax,locse,pkse+0.3,'v','Color',colgray(j,:),'MarkerFaceColor',colgray(j,:))
            end
            legm=sprintf('Model (%d s)',ttarget);
            sel=m.tm==ttarget;    
            phim_scaled=m.get_scaled_dist(m.x0,m.dx,m.phim_dist(:,sel),dataKCl{iexp}.dexp);
            l(length(l)+1)=semilogx(ax,m.x0*1e6,phim_scaled*100,'-','Color',colm(iexp,:),'DisplayName',legm,plotcfg{:});
            % Peaks
            [pks0,locs0,w2,p2] = findpeaks(phim_scaled*100,m.x0.*1e6);
            plot(ax,locs0,pks0+0.3,'v','Color',colm(iexp,:),'MarkerFaceColor',colm(iexp,:))
            
            axlim=[10,1500];
            xlim(ax,axlim);
            xticks_=[1 10 100 1000];
            xticks(ax,[1 10 100 1000])
            xticklabels(ax,arrayfun(@(a) sprintf('%.0f',a),xticks_,'UniformOutput',false))    
            %New_XTickLabel = get(ax,'xtick');
            %New_XTickLabel=[1,10,1000];
            %set(ax,'XTickLabel',New_XTickLabel)
            set(ax,'LineWidth',lw,fstics{:});
            xlabel(ax,'Diameter [\mum]',fs{:});
            ylabel(ax,'dQ_3 [%]',fs{:});
            legend(ax,l,'Location','northwest',fs{:});%[l0,l1],


            grid(ax,'on')
        end
        m.savefig_article(fig,fullfile(results_dir,names{iexp}));    
    end
end


names={'monomodal','bimodal','trimodal'};
names_fig={'Monomodal','Bimodal','Trimodal'};
plotcfgexp={'LineWidth',1.5};
plotcfg={'LineWidth',2};
fs={'FontSize',14};
fstics={'FontSize',12};

% Plot initial particle size distribution
if plot_init_dists
    
    coldist0=darkgreen;
    fig=figure('Position',[100 100 1800 600]);
    for iexp=1:3
        m=model_arr{iexp};
        datai=dataKCl{iexp};
        xexp=dataKCl{iexp}.xexp;
        phim_scaled=m.get_scaled_dist(m.x0,m.dx,m.phim_dist(:,1),dataKCl{iexp}.dexp);
        % exp dist
        ax=subplot(1,3,iexp);
        lege0=sprintf('Measured %s',names_fig{iexp});

        semilogx(ax,xexp*1e6,datai.phiexp(:,1)*100,'-o','Color',coldist0,'DisplayName',lege0,plotcfgexp{:});

        hold(ax,'on')
        %semilogx(ax,xexp*1e6,phim*100,'-','Color',colm(1,:),'DisplayName','hist from model',plotcfgexp{:});
        legm=sprintf('Model fit'); %names_fig{iexp}
        semilogx(ax,m.x0*1e6,phim_scaled*100,'-','Color',colm(iexp,:),'DisplayName',legm,plotcfg{:});

        sieve=infoKCl{iexp}.sieve;
        for j=1:length(datai.phi0exp)        
            lege=sprintf('%d-%d \\mum',sieve(j,1),sieve(j,2));
            semilogx(ax,xexp*1e6,m.w(j)*datai.phi0exp{j}*100,'--o','Color',colgray(j,:),'DisplayName',lege,plotcfgexp{:});
        end

        axlim=[1,1500];
        xlim(ax,axlim);
        xticks_=[1 10 100 1000];
        xticks(ax,[1 10 100 1000])
        xticklabels(ax,arrayfun(@(a) sprintf('%.0f',a),xticks_,'UniformOutput',false))    
        set(ax,'LineWidth',1.5,fstics{:});
        xlabel(ax,'Diameter [\mum]',fs{:});
        ylabel(ax,'dQ_3 [%]',fs{:});
        legend(ax,'Location','northwest',fs{:});
        grid(ax,'on')

    end
    m.savefig_article(fig,fullfile(results_dir,'initdist'));
end

% Plot initial experimental particle size distribution
if plot_exp_dists
    fig=figure('Position',[100 100 1800 600]);
    for iexp=1:3
        infoi=infoKCl{iexp};
        datai=dataKCl{iexp};
        xexp=dataKCl{iexp}.xexp;
        dxexp=dataKCl{iexp}.dxexp;
        
        
        sieve=infoKCl{iexp}.sieve;
        phiexpdist=zeros(length(dxexp),1);
        for j=1:length(datai.phi0exp)
            phi0exp_dist=datai.phi0exp{j}./dxexp;
            phiexpdist=phiexpdist+infoKCl{iexp}.w(j).*phi0exp_dist;        
        end    
        lege=sprintf('%s',names_fig{iexp});
        phiexpdist=phiexpdist./sum(phiexpdist.*dxexp);
        ax=subplot(1,3,iexp);
        semilogx(ax,xexp*1e6,phiexpdist.*dxexp*100,'-o','Color',colm(iexp,:),'DisplayName',lege,plotcfgexp{:});
        
        axlim=[10,1500];
        xlim(ax,axlim);
        xticks_=[1 10 100 1000];
        xticks(ax,[1 10 100 1000])
        xticklabels(ax,arrayfun(@(a) sprintf('%.0f',a),xticks_,'UniformOutput',false))    
        set(ax,'LineWidth',1.5,fstics{:});
        xlabel(ax,'Diameter [\mum]',fs{:});
        ylabel(ax,'dQ_3 [%]',fs{:});
        legend(ax,'Location','northwest',fs{:});
        grid(ax,'on')
    end
    m.savefig_article(fig,fullfile(results_dir,'expdist'));
end
