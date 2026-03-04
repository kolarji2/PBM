%
%   Simple show-case of the PBM class usage
%

m = PBM; % create a new instance of the model

%View list of parameters that you have to set-up
m.check_settings(); 

%Some parameters must be setup manualy, others have initialization
%functions, check bellow
% set physico-chemical properties
m.rhof=1000; % 
m.rhos=1200;
m.etaf=0.001;
m.diff=1e-9;
m.csat=20;%g/l
m.Ms=58.5;

%Creates discretization, values are log10 of the size, units in [m]
%assigns x0, xb, dx, dxii
m.make_loggrid(-7,-2,400); 

%create initial particle size distribution
m.set_init_dist([5e-4],[5e-4],[1],[],'lognormal');

%Setups the concentration and model initial conditions.
cs0=6; %g/l
m.set_initial_conditions(0,cs0,'gl');

% Solves the model
time=[0:1:180];
m.solve(time);

% Visualize
figure
subplot(1,2,1)
for i=1:60:181
    semilogx(m.x0,m.phim_dist(:,i),'DisplayName',sprintf('%d s',time(i)))
    hold on
    legend('Location','best')
end
xlabel('size [m]')
ylabel('Volume distribution [m^{-1}]')
subplot(1,2,2)
plot(m.tm,m.cm)
xlabel('time [s]')
ylabel('c [g/l]')
