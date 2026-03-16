
% Plot solubility as pH function


marr={};

S0_25C=10^(-6.73); % mol/l
S0_55C=10^(-5.57); %5e-06; % mol/l
S0_65C=10^(-5.23); % mol/l
Tarr=[25,55,65];
S0arr=[S0_25C,S0_55C,S0_65C];
leg={'25°C','55°C','65°C'};

logS=cell(3,1);
pHarr=linspace(0,14,100)';
fig=figure;
ax=gca;
for i=1:length(Tarr)
    m=PBMTelmi();
    useactivity=false;
    cNa2CO3=0;
    m.init_telmi_comp_solver('MuzikEtAl/expdata/pka_dist.xlsx');
    m.set_props(Tarr(i)+273.15,0,useactivity,@(x)0);
    Scalc=S0arr(i)*m.MtelH*m.get_ceqred(m.pH2cH(pHarr));
    Scalc(Scalc>m.MtelH)=m.MtelH;
    S{i}=Scalc;
    
    semilogy(ax,pHarr,S{i},'DisplayName',leg{i},'LineWidth',1.2)
    hold(ax,'on')
end
grid on
ylabel(ax,'Solubility [g/l]')  % g/l = mg/ml
xlabel(ax,'pH')
legend(ax,'Location','n')
u=Utils();
u.savefig_article(fig,'MuzikEtAl/results/visu_solubility')

% Save to excel
mdata=struct;
mdata.pHarr=pHarr;
mdata.S25_gl=S{1};
mdata.S55_gl=S{2};
mdata.S65_gl=S{3};
mdata.logS25_gl=log10(S{1});
mdata.logS55_gl=log10(S{2});
mdata.logS65_gl=log10(S{3});
mdata.logS25_moll=log10(S{1}./m.MtelH);
mdata.logS55_moll=log10(S{2}./m.MtelH);
mdata.logS65_moll=log10(S{3}./m.MtelH);

fields={{'pHarr','S25_gl','S55_gl','S65_gl','logS25_gl','logS55_gl','logS65_gl','logS25_moll','logS55_moll','logS65_moll'}};
fname='MuzikEtAl/results/predicted_solubilityTH.xlsx';
m.save_to_excel({mdata},fields,fname,'solubility');  