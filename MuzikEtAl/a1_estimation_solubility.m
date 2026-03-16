%
%   Example of usage of modified PBM class for 
%		modelling of pH dependent dissolution of telmisartan.
%		
%		Modified @PBM class can be found in directory: @PBMTelmi 
%
% Estimates intrinsic solulbility of Telmisartan
% Uses simplified Apelblat equation to interpolate between 25C and 65C
%

m=PBMTelmi();
loadData=0;
if ~exist('TimeEnd45','var')
    [exps,TimeEnd45,TimeEnd55,TimeEnd65,TimeEndData,prop] = load_exp_data_telmisartan(m);
    [phi0expdata] = load_psd_telmisartan();
end

wpure_Na2CO3=0.928;
wpure_Telmi=0.996;

% Experimental data
dataSol=[65, 8.21, 2.21331260483574; ...
65, 8.21, 2.27672544286708; ...
65, 10.24, 246.123634578407];




Tarr=dataSol(:,1);
pHarr=dataSol(:,2);
csol=dataSol(:,3);
cTintrinsic_arr=zeros(size(csol));

m=PBMTelmi();
useactivity=false;
cNa2CO3=0.2;
m.init_telmi_comp_solver('MuzikEtAl/expdata/pka_dist.xlsx');


for i=1:length(csol)
    m.set_props(Tarr(i)+273.15,cNa2CO3,useactivity,prop.etafTelmi230);
    cTintrinsic_arr(i)=csol(i)/m.MtelH/m.get_ceqred(m.pH2cH(pHarr(i)));
end


%
S0_25C=10^(-6.73); % mol/l https://pub.iapchem.org/ojs/index.php/admet/article/view/686/pdf
S0_65C=mean(cTintrinsic_arr); % mol/l

%% Fit solubility to function
% Solubility Apelblat equation
% ln x = A+B/T +c*ln(T)
% Simplified Apelblat for ideal solution
% ln x = A+B/T
S0arr=[S0_25C,S0_65C];

logS0=log(S0arr);
Tarr=[25,65]+273.15;
logSfun=@(T,a) a(1)+a(2)./T;
objfun=@(a) mean((logSfun(Tarr,a)-logS0).^2);
res=fminsearch(objfun,[1,1]);


fprintf('\n')
% S0funT=@(T) exp(13.70760 + (-8707.17148./T));
fprintf('Obtained fit: S0 = exp(%.5f + (%.5f./T))\n',res(1),res(2));
fprintf('ln(S0) = %.5f + (%.5f/T)\n',res(1),res(2))
fprintf('Fit MSE: %e\n', objfun(res))
fprintf('\n')

figure
subplot(1,2,1)
plot(Tarr-273.15,logS0,'o','DisplayName','Experiment')
Tarr_fine=linspace(273.15+25,273.15+66,100);
hold on
plot(Tarr_fine-273.15,logSfun(Tarr_fine,res),'DisplayName','Apelblat (fit)')
ylabel('logS0')
xlabel('T [C]')
legend()

subplot(1,2,2)
plot(Tarr-273.15,S0arr,'o','DisplayName','Experiment')
Tarr_fine=linspace(273.15+25,273.15+66,100);
hold on
plot(Tarr_fine-273.15,exp(logSfun(Tarr_fine,res)),'DisplayName','Apelblat (fit)')
ylabel('S0')
xlabel('T [C]')
legend()

S0_55C=exp(logSfun([55+273.15],res));  % interpolated value

fprintf('Estimated solubility:\n')
fprintf('S0(25C) = %.4e mol/l\n',S0_25C)
fprintf('S0(55C) = %.4e mol/l\n',S0_55C)
fprintf('S0(65C) = %.4e mol/l\n',S0_65C)
fprintf('\n')
% logS0(25C) = -6.73 mol/l
% logS0(55C) = -5.57 mol/l
% logS0(65C) = -5.23 mol/l
fprintf('log10[S0(25C)] = %.2f mol/l\n',log10(S0_25C))
fprintf('log10[S0(55C)] = %.2f mol/l\n',log10(S0_55C))
fprintf('log10[S0(65C)] = %.2f mol/l\n',log10(S0_65C))
fprintf('\n')
fprintf('ln[S0(25C)] = %.2f mol/l\n',log(S0_25C))
fprintf('ln[S0(55C)] = %.2f mol/l\n',log(S0_55C))
fprintf('ln[S0(65C)] = %.2f mol/l\n',log(S0_65C))


