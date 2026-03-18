%
%   Optimal content of 1um particle fraction for stable release at
%   different hydrolysis rate
%
%       Run a9_quasi_steady_state.m before this script
%       or load SonntagEtAl/expdata/quasi_steady_data.mat

%load('SonntagEtAl/expdata/quasi_steady_data.mat')

utils=Utils();
fs=21;
dts_1_100um=load('SonntagEtAl/expdata/opt1_100um_fine.mat');
%dts_1_100um=load('SonntagEtAl/expdata/opt1_100um_e100.mat');
opt_model_arr=dts_1_100um.opt_model_arr;
w100um=dts_1_100um.w100um;
tstable_arr=dts_1_100um.tstable_arr;
cs0_arr
cfg={'LineWidth',2};
txtcfg={'FontSize',fs,'Interpreter','latex'};
%fig=figure('Position',[100,100,900,600]);
%ax1=subplot(1,2,1);
%ax2=subplot(1,2,2);
fig1=figure('Position',[100,100,600,600]);
ax1=gca(fig1);


for ik=1:Nk
    if kmultiplier(ik)==1
        title0='Hydrolysis rate 1x';
    else
        title0=sprintf('Hydrolysis slowed %dx',1/kmultiplier(ik));
    end
        
    plot(ax1,cs0_arr,1-w100um(:,ik),'-o','DisplayName',title0,cfg{:})
    for ic=1:3
       text(ax1,cs0_arr(ic),1-w100um(ic,ik),sprintf('%.0f mg/l',(1-w100um(ic,ik)).*cs0_arr(ic)*1e3), ...
           'HorizontalAlignment','left','VerticalAlignment','bottom',txtcfg{:},'FontSize',fs-2);
    end
    %plot(cs0_arr,(1-w100um(:,ik)).*cs0_arr','-o','DisplayName',title0,cfg{:})
    hold(ax1,'on')
end

xlabel(ax1,'Initial solids conc. [g/l]',txtcfg{:})
ylabel(ax1,'Weight frac. of 1$\mu$m part. [-]',txtcfg{:})
legend(ax1,'Location','ne',txtcfg{:})
grid(ax1,'on')
set(ax1,'FontSize',fs,'TickLabelInterpreter','latex')
%ylim([0.86,1.01])
ylim(ax1,[0,0.15])
xlim(ax1,[0 7])
utils.savefig_article(fig1,'SonntagEtAl/results/opt_composition_article');

% Particle size example
fig2=figure('Position',[100,100,600,600]);
ax2=gca(fig2);
m=opt_model_arr{1,3};
w=[1-0.87 0.87];
semilogx(ax2,m.x0*1e6,w.*m.phi0_distestAll.*m.dx*100,cfg{:});
xlim(ax2,[5e-1,3e2])
%ylim(ax2,[0,0.17])
grid(ax2,'on')

legend(ax2,{'1 $\mu$m (13\%)','100 $\mu$m (87\%)'},'Location','nw',txtcfg{:})
xlabel(ax2,'Diameter [$\mu$m]',txtcfg{:})
ylabel(ax2,'Volume percent',txtcfg{:})
%title('Particle size example (w=0.87)')

set(ax2,'FontSize',fs,'TickLabelInterpreter','latex')
utils.savefig_article(fig2,'SonntagEtAl/results/opt_particle_size_article');


