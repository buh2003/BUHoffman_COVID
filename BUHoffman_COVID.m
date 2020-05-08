%% ################# SEIAQHRRRDP Model ######################
set(0,'DefaultFigureWindowStyle','docked')

% S:   number of succeptibles
% E:   number of exposed
% I:   number of infectious (not quarantined)
% A:   number of infectious Undocumented
% Q:   number of quarantined, active cases not requiring hospitalization
% H:   number of hospitalized, active cases requiring hospitalization
% Ra:  number of recovered Undocumented cases
% Rq:  number of recovered cases not requiring hospitalization
% Rh:  number of recovered cases requiring hospitalization
% D:   number of dead
% P:   number of protected (by public health measures, ie social distancing)

%% Import Data
state = 'NY';
[Confirmed,Hospitalized_tot,Hospitalized_active,Discharged,Deaths,Time] = getUS_Covid_data(state);

% fix Hospt_total data, as they are currently incomplete
% from Cuomo press videos (updated 4/29)
date_updated = 56;
newhosp = zeros(size(Hospitalized_tot));
newhosp(16:date_updated,1) = [496;605;696;836;923;1263;1564;1777;1925;1959;2210;2411;2736;3038;3177;3165;2721;2486;2563;...
    2825;2923;2772;2689;2389;2119;2016;2039;2156;2045;1833;1616;1389;1408;1404;1367;1224;...
    1119;1076;973;970;933];
Hospitalized_tot = [Hospitalized_tot(1:15,1); cumtrapz(newhosp(16:date_updated,1))];

% estimate incomplete Disharge data
Discharged(1:date_updated,1) = Hospitalized_tot(1:date_updated,1)-Deaths(1:date_updated,1)-Hospitalized_active(1:date_updated,1);
Discharged(Discharged < 0) = 0;

% Data is incomplete as of 4/29/20 stop downloading after then
Confirmed = Confirmed(1:date_updated,:);
Hospitalized_tot = Hospitalized_tot(1:date_updated,:);
Hospitalized_active = Hospitalized_active(1:date_updated,:);
Discharged = Discharged(1:date_updated,:);
Deaths = Deaths(1:date_updated,:);
Time = Time(1:date_updated,:);

figure;plot(Time, Confirmed, Time, Hospitalized_tot,Time,Hospitalized_active, Time, Discharged, Time, Deaths); legend('Confirmed','Hospitalized Total','Hospitalized Active','Discharged','Deaths');

%% Population estimate for NYS data (selecting only major affected counties as this is where the data is primarily from)
% Dat: https://www.census.gov/quickfacts/fact/table/
% Major affected counties: 
Queens = 2.25e6; Kings= 2.56e6; Nassau = 1.36e6; Bronx = 1.42e6; Suffolk = 1.42e6; Westchester = 0.97e6;
Manhattan = 1.63e6; Rockland = 0.33e6; Richmond = 0.48e6; Orange = 0.38e6; Dutchess = 0.29e6;  Erie = 0.92e6; Monroe = 0.74e6;
Npop = Queens + Kings + Nassau + Bronx + Suffolk + Westchester + Manhattan + Rockland + Richmond + Orange + Dutchess + Erie + Monroe;

% Land area (miles^2)
Queens_a = 108.53; Kings_a = 70.82; Nassau_a = 284.72; Bronx_a = 42.1; Suffolk_a = 912.05; Westchester_a = 430.5;
Manhattan_a = 22.83; Rockland_a = 173.55; Richmond_a = 58.37; Orange_a = 811.69; Dutchess_a = 795.63;  Erie_a = 1042.69; Monroe_a = 657.21;
Total_a = Queens_a + Kings_a + Nassau_a + Bronx_a + Suffolk_a + Westchester_a + Manhattan_a + Rockland_a + Richmond_a + Orange_a + Dutchess_a + Erie_a + Monroe_a;

% Population density (#1000s of people per mile^2)
Pop_density = (Npop/Total_a)/1000;

clearvars Queens Kings Nassau Bronx Suffolk Westchester Manhattan Rockland Richmond Orange Dutchess Erie Monroe...
    Queens_a Kings_a Nassau_a Bronx_a Suffolk_a Westchester_a Manhattan_a Rockland_a Richmond_a Orange_a Dutchess_a Erie_a Monroe_a Total_a

%% ###### Estimate parameters #####
%  ###### For parameter estimation a simplified model is used where Q and Rq are combined to fit the reported data from NYC

% S:   to be modeled
% E:   to be modeled
% I:   to be modeled
% A:   to be modeled
% Q:   Not_Hospitalized = Confirmed - Hospitalized_tot (number of quarantined, all cases not
%       requiring hospitalization)
% H:   Hospitalized_active (less confident about these data)
% Ra:  to be modeled
% Rh:  Discharged (number of recovered cases requiring hospitalization)
% D:   Deaths 
% P:   to be modeled

% ## Set Asympotomtic Rate
Asymp_rate = 0.75;

% ## Set initial conditions for states
E0 =  (Confirmed(1)+(Confirmed(1)*Asymp_rate)/(1-Asymp_rate))*20;
I0 =  0;%(Confirmed(1)+(Confirmed(1)*Asymp_rate)/(1-Asymp_rate));
A0 =  0;%(Confirmed(1)*Asymp_rate)/(1-Asymp_rate);
Q0 =  0;
H0 =  0;
Ra0 = 0;
Rh0 = 0;
D0 =  0;
P0 =  0;
S0 =  Npop - E0 - I0 - A0 - Q0 - H0 - Ra0 - Rh0 - D0 - P0;
ic =  [S0; E0; I0; A0; Q0; H0; Ra0; Rh0; D0; P0];

% ## Set initial guess of paramters
Coef(1) = 0.1046;  % alpha: effective protection rate
Coef(2) = 0.0162;  % zeta: effective protection leak rate
Coef(3) = 0.2809;  % beta: effective contact rate
Coef(4) = 0.1752;  % lambda: basic intervention recovery rate
Coef(5) = 0.2224;  % epsilon: intervention recovery development rate constant;
Coef(6) = 0.5106;   % sigma: intervention recovery efficacy;
Coef(7) = 0.0957;  % kappa: basic intervention death rate
Coef(8) = 0.4069;  % psi: death intervention development rate constant
Coef(9) = 2.6891;   % mu: intervention death efficacy;

% #### Fixed Coefficients #####
rho = Pop_density; 
gamma = 1/3;   % gamma: effective exposed to infected rate
delta = 1/5;   % delta: effective hospitalization rate

% theta is constrained to delta to represent (1 - hospitalization rate), with a delay of hospitalization by 1 days
hosp_rate = (Hospitalized_tot(1:end,1)./(Confirmed(1:end)));
non_hosp_rate_ratio = ((1-hosp_rate)./hosp_rate);
theta = delta*mean(non_hosp_rate_ratio(46:56,1),1); %effective infected to quarantined rate 

% ## Aymptomatic rate calculation
nu = (Asymp_rate*(theta+delta))/(1-Asymp_rate);
phi = 1/9;
clearvars hosp_rate non_hosp_rate_ratio; 
SD_delay = 18;  % Days to start Social Distancing measure from first case
SD_remove = 1;
zeta_factor = 1;

% Immunity
Immunity_fraction = 0.7;% fraction of recovered who develop immunity
tau = 1-Immunity_fraction; % fraction of recovered who DO NOT develop immunity
pi = 0; % 1/(Duration of immunity)

% ##### Initialize data ####
Not_Hospitalized = Confirmed - Hospitalized_tot;
yvals = [Not_Hospitalized'; Hospitalized_active'; Discharged'; Deaths'];
yvals(isnan(yvals))=0;

bounds_factor = 1/1.5;
Coef1_bounds = [0.030,0.170];
Coef2_bounds = [0.014,0.0162];
Coef3_bounds = [Coef(3)-(Coef(3)*bounds_factor),Coef(3)+(Coef(3)*bounds_factor)];
Coef4_bounds = [Coef(4)-(Coef(4)*bounds_factor),Coef(4)+(Coef(4)*bounds_factor)];
Coef5_bounds = [Coef(5)-(Coef(5)*1.5*bounds_factor),Coef(5)+(Coef(5)*2*bounds_factor)];
Coef6_bounds = [Coef(6)-(Coef(6)*1.5*bounds_factor),Coef(6)+(Coef(6)*2*bounds_factor)];
Coef7_bounds = [Coef(7)-(Coef(7)*bounds_factor),Coef(7)+(Coef(7)*bounds_factor)];
Coef8_bounds = [Coef(8)-(Coef(8)*bounds_factor),Coef(8)+(Coef(8)*bounds_factor)];
Coef9_bounds = [Coef(9)-(Coef(9)*bounds_factor),Coef(9)+(Coef(9)*bounds_factor)];

LB = [Coef1_bounds(1),Coef2_bounds(1),Coef3_bounds(1),Coef4_bounds(1),Coef5_bounds(1),Coef6_bounds(1),Coef7_bounds(1),Coef8_bounds(1),Coef9_bounds(1)];
UB = [Coef1_bounds(2),Coef2_bounds(2),Coef3_bounds(2),Coef4_bounds(2),Coef5_bounds(2),Coef6_bounds(2),Coef7_bounds(2),Coef8_bounds(2),Coef9_bounds(2)];

Coef_r = optimvar('Coef_r',9,"LowerBound",LB,"UpperBound",UB);
tspan = 1:length(Time);

% Import Flu/HCov data to estimate seasonal variance in infectivitiy
% From CDC FluView Website
% https://gis.cdc.gov/grasp/fluview/fluportaldashboard.html
% Averaged "WeightedILI" column over 1997-2020
% 2014-2017 From https://www.sciencedirect.com/science/article/pii/S1386653218300325#bib0035
% 2018-2020 From https://www.cdc.gov/surveillance/nrevss/coronavirus/index.html
flufactor = 5;
chi = importflu(Time,[tspan(1) tspan(end)],'factor',flufactor);

myfcn = fcn2optimexpr(@SEIAQHRRDP_ConvODE2, Coef_r, tspan, Npop, ic, rho, gamma, nu, phi, delta, theta, tau, pi, SD_delay, SD_remove, zeta_factor, chi);
obj = sum(sum((myfcn - yvals).^2));
prob = optimproblem("Objective",obj);
constraint1 = Coef_r(7) <= Coef_r(4);
prob.Constraints.constraint1 = constraint1;
Coef_r0.Coef_r = Coef;
opts = optimoptions(@fmincon,'MaxIterations', 6e2,'MaxFunctionEvaluations',2e3,...
    'Display','iter','Diagnostics','on','StepTolerance',1e-16,'OptimalityTolerance',1e-16,...
    'Algorithm','sqp');

tic
[rsol,fval,flag,output,lagrange] = solve(prob, Coef_r0,'Options',opts);
toc

%% Figures S2-3: Sensitivity analysis, bivariate, effect on Deaths, and global sensitivity analysis

% Immunity
Immunity_fraction = 0.70;% fraction of recovered who develop immunity
tau = 1-Immunity_fraction; % fraction of recovered who DO NOT develop immunity
pi = 0/1825; % 1/(Duration of immunity)

Coef1s = Coef1_bounds(1):(Coef1_bounds(2)-Coef1_bounds(1))/19:Coef1_bounds(2);
Coef2s = 0.0012:(0.03-0.0012)/19:0.03;
Coef3s = Coef3_bounds(1):(Coef3_bounds(2)-Coef3_bounds(1))/19:Coef3_bounds(2);
Coef4s = Coef4_bounds(1):(Coef4_bounds(2)-Coef4_bounds(1))/19:Coef4_bounds(2);
Coef5s = Coef5_bounds(1):(Coef5_bounds(2)-Coef5_bounds(1))/19:Coef5_bounds(2);
Coef6s = Coef6_bounds(1):(Coef6_bounds(2)-Coef6_bounds(1))/19:Coef6_bounds(2);
Coef7s = Coef7_bounds(1):(Coef7_bounds(2)-Coef7_bounds(1))/19:Coef7_bounds(2);
Coef8s = Coef8_bounds(1):(Coef8_bounds(2)-Coef8_bounds(1))/19:Coef8_bounds(2);
Coef9s = Coef9_bounds(1):(Coef9_bounds(2)-Coef9_bounds(1))/19:Coef9_bounds(2);

Coefs_cells = cell(1,9);
Coefs_cells{1,1} = Coef1s;
Coefs_cells{1,2} = Coef2s;
Coefs_cells{1,3} = Coef3s;
Coefs_cells{1,4} = Coef4s;
Coefs_cells{1,5} = Coef5s;
Coefs_cells{1,6} = Coef6s;
Coefs_cells{1,7} = Coef7s;
Coefs_cells{1,8} = Coef8s;
Coefs_cells{1,9} = Coef9s;

% WINDOW FOR SIMULATION
window_sim = [1, 182];
chi = importflu(Time,window_sim,'factor',flufactor);
megayout = cell(1,45);
Deathsout = zeros(20,20);

tic
count = 0;
for g = 1:9
    for h = g:9
        Coefs_combo = Coef;
        Coefs_cells_1_g = Coefs_cells{1,g};
        Coefs_cells_1_h = Coefs_cells{1,h};   
        parfor i = 1:20
            Coefs_combo_parfor = Coefs_combo;
            Coefs_combo_parfor(g) = Coefs_cells_1_g(i);  
            for j = 1:20             
                if g == h
                    if i == j
                        Coefs_combo_parfor(h) = Coefs_cells_1_h(j);
                        [~,y] = ode45(@(t,y) SEIAQHRRRDP_deqs(t, y, Coefs_combo_parfor, Npop, rho, gamma,...
                            nu, phi, delta, theta, omega,tau, pi, SD_delay, ...
                            SD_remove, zeta_factor, window_sim, chi), window_sim(1):window_sim(end), ic);
                        Deathsout(i,j) = y(end,10);
                    else
                        Deathsout(i,j) = NaN;
                    end
                else
                    Coefs_combo_parfor(h) = Coefs_cells_1_h(j);
                    [~,y] = ode45(@(t,y) SEIAQHRRRDP_deqs(t, y, Coefs_combo_parfor, Npop, rho, gamma,...
                        nu, phi, delta, theta, omega,tau, pi, SD_delay, ...
                        SD_remove, zeta_factor, window_sim, chi), window_sim(1):window_sim(end), ic);
                    Deathsout(i,j) = y(end,10); 
                end
            end
        end
        disp([g h]);
        count = count +1;
        megayout{1,count} = Deathsout;
    end
end
toc

decoder = [9 17 24 30 35 39 42 44 45;...
           8 16 23 29 34 38 41 43 44;...
           7 15 22 28 33 37 40 41 42;...
           6 14 21 27 32 36 37 38 39;...
           5 13 20 26 31 32 33 34 35;...
           4 12 19 25 26 27 28 29 30;...
           3 11 18 19 20 21 22 23 24;...
           2 10 11 12 13 14 15 16 17;...
           1 2 3 4 5 6 7 8 9];

maxmax = log10(max(max(cell2mat(megayout))));
minmin = log10(min(min(cell2mat(megayout))));
hand = figure;
count = 1;
for i = 1:9 % rows
    ii = 10-i;
    for j =1:9 % columns
        if ii >= j
            maxdeaths = rot90(rot90(megayout{1,decoder(i,j)},2),-1);
        else
            maxdeaths = flipud(megayout{1,decoder(i,j)});
        end
            
        subplot(9,9,count); heatmap(Coefs_cells{1,j},Coefs_cells{1,ii},log10(maxdeaths),...
            'GridVisible','off','CellLabelColor','none','Colormap',parula,...
            'ColorbarVisible','off',...
            'ColorLimits',[minmin maxmax],'XLabel',[],'YLabel',[],...
            'XDisplayLabels',{'','','','','','','','','','','','','','','','','','','',''},...
            'YDisplayLabels',{'','','','','','','','','','','','','','','','','','','',''})
        count = count +1;
    end
end

% For colorbar scale
figure;
heatmap(Coefs_cells{1,1},Coefs_cells{1,1},log10(maxdeaths),...
            'GridVisible','off','CellLabelColor','none','Colormap',parula,...
            'ColorLimits',[minmin maxmax])
        
axs = struct(gca);
cb = axs.Colorbar;
cb.Ticks = [...
    log10(5e1) log10(6e1)...
        log10(7e1) log10(8e1) log10(9e1)...
    log10(1e2) log10(2e2) log10(3e2) log10(4e2) log10(5e2) log10(6e2)...
        log10(7e2) log10(8e2) log10(9e2)...
    log10(1e3) log10(2e3) log10(3e3) log10(4e3) log10(5e3) log10(6e3)...
        log10(7e3) log10(8e3) log10(9e3)...
    log10(1e4) log10(2e4) log10(3e4) log10(4e4) log10(5e4) log10(6e4)...
        log10(7e4) log10(8e4) log10(9e4)...
    log10(1e5) log10(2e5) log10(3e5) log10(4e5) log10(5e5) log10(6e5)...
    ]; 

% Sobol sensitivity analysis
% Using the "Global Sensitivity Analysis Toolbox"
% Citation :flax (2020). Global Sensitivity Analysis Toolbox (https://www.mathworks.com/matlabcentral/fileexchange/40759-global-sensitivity-analysis-toolbox), MATLAB Central File Exchange. Retrieved May 8, 2020.

% create a new project 
pro = pro_Create();

% add 5 input variables to the project, named param*, distributed in the 
% range [0 1] and indicate that the variables will be sampled following a 
% Sobol set
pro = pro_AddInput(pro, @()pdf_Sobol(Coef1_bounds), 'param1');
pro = pro_AddInput(pro, @()pdf_Sobol([0.01 0.02]), 'param2');
pro = pro_AddInput(pro, @()pdf_Sobol(Coef3_bounds), 'param3');
pro = pro_AddInput(pro, @()pdf_Sobol(Coef4_bounds), 'param4');
pro = pro_AddInput(pro, @()pdf_Sobol(Coef5_bounds), 'param5');
pro = pro_AddInput(pro, @()pdf_Sobol(Coef6_bounds), 'param6');
pro = pro_AddInput(pro, @()pdf_Sobol(Coef7_bounds), 'param7');
pro = pro_AddInput(pro, @()pdf_Sobol(Coef8_bounds), 'param8');
pro = pro_AddInput(pro, @()pdf_Sobol(Coef9_bounds), 'param9');

Rq0 = 0;
ic = [S0; E0; I0; A0; Q0; H0; Ra0; Rq0; Rh0; D0; P0];
omega = 1/14;

% Immunity
Immunity_fraction = 0.70;% fraction of recovered who develop immunity
tau = 1-Immunity_fraction; % fraction of recovered who DO NOT develop immunity
pi = 0/1825; % 1/(Duration of immunity)

% WINDOW FOR SIMULATION
window_sim = [1, 182];

chi = importflu(Time,window_sim,'factor',flufactor);

% Simulate based on parameters
pro = pro_SetModel(pro, @(x) ODEsforSensitivity(x, Npop, rho, gamma, nu, phi,...
    delta, theta, omega,tau, pi, SD_delay, SD_remove, zeta_factor,chi, ic, window_sim), 'ODEs');

% set the number of samples for the quasi-random Monte Carlo simulation
pro.N = 10000;

tic
% initialize the project by calculating the model at the sample points
pro = GSA_Init_MultiOut_MultiSI(pro);
toc

sound_end=load('handel.mat');
sound(sound_end.y,sound_end.Fs);

%
% Calculate the first order global sensitivity coefficient for the set of
% variables comprising only the first FUNDAMENTAL, DISTRIBUTED parameter
[S1, eS1, pro] = GSA_GetSy_MultiOut_MultiSI(pro, {1});
disp(1)
% Calculate the first order global sensitivity coefficient for the set of
% variables comprising only the second FUNDAMENTAL, DISTRIBUTED parameter
[S2, eS2, pro] = GSA_GetSy_MultiOut_MultiSI(pro, {2});
disp(2)
% Calculate the first order global sensitivity coefficient for the set of
% variables comprising only the third FUNDAMENTAL, DISTRIBUTED parameter
[S3, eS3, pro] = GSA_GetSy_MultiOut_MultiSI(pro, {3});
disp(3)
% Calculate the first order global sensitivity coefficient for the set of
% variables comprising only the fourth FUNDAMENTAL, DISTRIBUTED parameter
[S4, eS4, pro] = GSA_GetSy_MultiOut_MultiSI(pro, {4});
disp(4)
% Calculate the first order global sensitivity coefficient for the set of
% variables comprising only the fifth FUNDAMENTAL, DISTRIBUTED parameter
[S5, eS5, pro] = GSA_GetSy_MultiOut_MultiSI(pro, {5});
disp(5)
% Calculate the first order global sensitivity coefficient for the set of
% variables comprising only the fifth FUNDAMENTAL, DISTRIBUTED parameter
[S6, eS6, pro] = GSA_GetSy_MultiOut_MultiSI(pro, {6});
disp(6)
% Calculate the first order global sensitivity coefficient for the set of
% variables comprising only the fifth FUNDAMENTAL, DISTRIBUTED parameter
[S7, eS7, pro] = GSA_GetSy_MultiOut_MultiSI(pro, {7});
disp(7)
% Calculate the first order global sensitivity coefficient for the set of
% variables comprising only the fifth FUNDAMENTAL, DISTRIBUTED parameter
[S8, eS8, pro] = GSA_GetSy_MultiOut_MultiSI(pro, {8});
disp(8)
% Calculate the first order global sensitivity coefficient for the set of
% variables comprising only the fifth FUNDAMENTAL, DISTRIBUTED parameter
[S9, eS9, pro] = GSA_GetSy_MultiOut_MultiSI(pro, {9});
disp(9)
toc

Si_sobol = [S1 S2 S3 S4 S5 S6 S7 S8 S9];

figure; 
subplot(4,3,1); bar(Si_sobol(1,:)); title('Susceptibles'); ylim([0 1])
subplot(4,3,2); bar(Si_sobol(2,:)); title('Exposed'); ylim([0 1])
subplot(4,3,3); bar(Si_sobol(3,:)); title('Infected'); ylim([0 1])
subplot(4,3,4); bar(Si_sobol(4,:)); title('Undocumented Infected'); ylim([0 1])
subplot(4,3,5); bar(Si_sobol(5,:)); title('Quarantined'); ylim([0 1])
subplot(4,3,6); bar(Si_sobol(6,:)); title('Hospitalized'); ylim([0 1])
subplot(4,3,7); bar(Si_sobol(7,:)); title('Recovered Undocumented'); ylim([0 1])
subplot(4,3,8); bar(Si_sobol(8,:)); title('Recovered Quarantined'); ylim([0 1])
subplot(4,3,9); bar(Si_sobol(9,:)); title('Recovered Hospitalized'); ylim([0 1])
subplot(4,3,10); bar(Si_sobol(10,:)); title('Dead'); ylim([0 1])
subplot(4,3,11); bar(Si_sobol(11,:)); title('Protected'); ylim([0 1])

%% Figure S4 -- Residual Analysis

Coef = rsol.Coef_r';

Rq0 = 0;
ic = [S0; E0; I0; A0; Q0; H0; Ra0; Rq0; Rh0; D0; P0];

omega = 1/14;

% Immunity
Immunity_fraction = 0.70;% fraction of recovered who develop immunity
tau = 1-Immunity_fraction; % fraction of recovered who DO NOT develop immunity
pi = 0/1825; % 1/(Duration of immunity)

% WINDOW FOR SIMULATION
window_sim = [1, size(Confirmed,1)];

chi = importflu(Time,window_sim,'factor',flufactor);

% Simulate based on parameters
[t,y] = ode45(@(t,y) SEIAQHRRRDP_deqs(t, y, Coef, Npop, rho, gamma, nu, phi, delta, theta, omega,tau, pi, SD_delay, SD_remove, zeta_factor, window_sim,chi), window_sim(1):window_sim(2), ic);

residual_times = Time(1:size(Confirmed,1));

figure;
subplot(5,3,1); plot(residual_times,(sum(y(1:size(Confirmed,1),[5,6,8:10]),2)),'LineWidth',2); title({'# Total Confirmed Infections'});  xlim([residual_times(1) residual_times(end)]);
    hold on; scatter(residual_times,(Confirmed),'LineWidth',1,'MarkerEdgeColor','#D95319')    
subplot(5,3,2); scatter(residual_times, (Confirmed)-(sum(y(1:size(Confirmed,1),[5,6,8:10]),2)),'LineWidth',1); ylim([-10000 10000]);
    hold on; plot(residual_times,zeros(size(residual_times,1)),'r--')
    title({'Residuals';'# Total Confirmed Infections'});
subplot(5,3,3); qqplot((Confirmed)-(sum(y(1:size(Confirmed,1),[5,6,8:10]),2)));
    title({'QQ-plot Residuals';'# Total Confirmed Infections'});

subplot(5,3,4); plot(residual_times,(sum(y(1:size(Confirmed,1),[6,9,10]),2)),'LineWidth',2); title({'# Total Confirmed Hospitalized'}); xlim([residual_times(1) residual_times(end)]);
    hold on; scatter(residual_times,(Hospitalized_active+Deaths+Discharged),'LineWidth',1,'MarkerEdgeColor','#D95319')
subplot(5,3,5); scatter(residual_times, (Hospitalized_active+Deaths+Discharged)-(sum(y(1:size(Confirmed,1),[6,9,10]),2)),'LineWidth',1); 
    xlim([residual_times(1) residual_times(end)]); ylim([-1000 1000]);
    title({'Residuals';'# Total Confirmed Hospitalized'});
    hold on; plot(residual_times,zeros(size(residual_times,1)),'r--')
subplot(5,3,6); qqplot((Hospitalized_active+Deaths+Discharged)-(sum(y(1:size(Confirmed,1),[6,9,10]),2)));  ylim([-1000 1000])
    title({'QQ-plot Residuals';'# Total Confirmed Hospitalized'});
    
subplot(5,3,7); plot(residual_times,y(1:size(Confirmed,1),6),'LineWidth',2); title({'# Active Hospitalized'}); xlim([residual_times(1) residual_times(end)]);
    hold on; scatter(residual_times,(Hospitalized_active),'LineWidth',1,'MarkerEdgeColor','#D95319')
subplot(5,3,8); scatter(residual_times, (Hospitalized_active)-y(1:size(Confirmed,1),6),'LineWidth',1); 
    xlim([residual_times(1) residual_times(end)]); ylim([-1000 1000]);
    hold on; plot(residual_times,zeros(size(residual_times,1)),'r--')
    title({'Residuals';'# Active Hospitalized'});
subplot(5,3,9); qqplot((Hospitalized_active)-y(:,6)); ylim([-1000 1000])
    title({'QQ-plot Residuals';'# Active Hospitalized'});

subplot(5,3,10); plot(residual_times,(y(1:size(Confirmed,1),9)),'LineWidth',2); title({'# Total Hospital Discharges'}); xlim([residual_times(1) residual_times(end)]);
    hold on; scatter(residual_times,(Discharged),'LineWidth',1,'MarkerEdgeColor','#D95319')
subplot(5,3,11); scatter(residual_times, (Discharged)-(y(1:size(Confirmed,1),9)),'LineWidth',1); 
    xlim([residual_times(1) residual_times(end)]); ylim([-1000 1000]);
    hold on; plot(residual_times,zeros(size(residual_times,1)),'r--')
    title({'Residuals';'# Total Hospital Discharges'});
subplot(5,3,12); qqplot((Discharged)-(y(1:size(Confirmed,1),9)));
    title({'QQ-plot Residuals';'# Total Hospital Discharges'});

subplot(5,3,13); plot(residual_times,y(1:size(Confirmed,1),10),'LineWidth',2); title('# Deaths'); xlim([residual_times(1) residual_times(end)]);
    hold on; scatter(residual_times,(Deaths),'LineWidth',1,'MarkerEdgeColor','#D95319')
subplot(5,3,14); scatter(residual_times, (Deaths)-(y(1:size(Confirmed,1),10)),'LineWidth',1); 
    xlim([residual_times(1) residual_times(end)]); ylim([-1000 1000]);
    hold on; plot(residual_times,zeros(size(residual_times,1)),'r--')
    title({'Residuals';'# Deaths'});
subplot(5,3,15); qqplot((Deaths)-(y(1:size(Confirmed,1),10)))
    title({'QQ-plot Residuals';'# Deaths'});



%%
% Estimate CI of parameters
btsrpNO = 50;
optboot = statset('UseParallel',true);

LB(2) = 0.0012;
mean_temp = Coef(2)+ 0.0005; % mean
var_temp =  Coef(2)*0.025; % variance
mu_temp = log((mean_temp^2)/sqrt(var_temp+mean_temp^2));
sigma_temp = sqrt(log(var_temp/(mean_temp^2)+1));
Coef2_UB = abs(lognrnd(mu_temp,sigma_temp,500,1));
mean(Coef2_UB)

tic
% parpool
bootCoef = bootstrp(btsrpNO,@(bootr) bootSolveODE_v2(bootr, tspan, Npop,...
    ic, rho, gamma, nu, phi, delta, theta, tau, pi, SD_delay, SD_remove, zeta_factor,...
    chi, UB, LB, yvals), Coef2_UB,'Options', optboot);
% delete(gcp)
toc

bootCoefCI =cell(1,size(bootCoef,2));
for i = 1:size(bootCoef,2)
    SEM = std(bootCoef(:,i))/sqrt(size(bootCoef,1));               % Standard Error
    ts = tinv([0.025  0.975],size(bootCoef,1)-1);      % T-Score
    bootCoefCI{1,i} = mean(bootCoef(:,i))+ ts*SEM; 
end
%
sound_end=load('handel.mat');
sound(sound_end.y,sound_end.Fs);

%% Figure S5

disp(rsol.Coef_r);

Coef = rsol.Coef_r';

Rq0 = 0;
ic = [S0; E0; I0; A0; Q0; H0; Ra0; Rq0; Rh0; D0; P0];

omega = 1/14;

% Immunity
Immunity_fraction = 0.70;% fraction of recovered who develop immunity
tau = 1-Immunity_fraction; % fraction of recovered who DO NOT develop immunity
pi = 0/1825; % 1/(Duration of immunity)

% WINDOW FOR SIMULATION
window_sim = [1, 182];

chi = importflu(Time,window_sim,'factor',flufactor);

% Simulate based on parameters
[t,y] = ode45(@(t,y) SEIAQHRRRDP_deqs(t, y, Coef, Npop, rho, gamma, nu, phi, delta, theta, omega,tau, pi, SD_delay, SD_remove, zeta_factor, window_sim,chi), window_sim(1):window_sim(2), ic);

% Estimate CI of simulation
bootyCI = simCI(bootCoef,Npop, rho, gamma, nu, phi, delta, theta, omega,tau, pi, SD_delay, SD_remove, zeta_factor, window_sim,chi, ic);

% Calculate Ro
R_naught = ((chi*rho*Coef(3))/(nu+theta+delta)) + ((chi*rho*Coef(3)*nu)/(phi*(nu+theta+delta)));
R_naught_peak = max(((chi*rho*Coef(3))/(nu+theta+delta)) + ((chi*rho*Coef(3)*nu)/(phi*(nu+theta+delta))));

% Figure S5
% Comprehensive Outputs
date_times = Time(1):(Time(1)+t(end)-1);
real_data_length = length(Confirmed);

fig=figure;
h=subplot(4,5,1); plotwithCI(date_times,y(:,2),bootyCI(:,2),real_data_length,h); title('# Exposed'); xlim([date_times(1) date_times(end)]); %ylim([0 y_limit2]);
h=subplot(4,5,2); plotwithCI(date_times,sum(y(:,3:4),2),bootyCI(:,3:4),real_data_length,h); title({'# Infected';'Not Quarantined'}); xlim([date_times(1) date_times(end)]); %ylim([0 y_limit2]);
h=subplot(4,5,3); plotwithCI(date_times,sum(y(:,4:10),2),bootyCI(:,4:10),real_data_length,h); title({'# Total Infections';'(Undocumented + confirmed)'}); xlim([date_times(1) date_times(end)]);
h=subplot(4,5,4); plotwithCI(date_times,sum(y(:,[5,6,8:10]),2),bootyCI(:,[5,6,8:10]),real_data_length,h); title({'# Total';'Confirmed Infections'}); xlim([date_times(1) date_times(end)]);
    hold on; scatter(date_times(1:length(Confirmed)),Confirmed,'LineWidth',1)
h=subplot(4,5,5); plotwithCI(date_times,sum(y(:,5:6),2),bootyCI(:,5:6),real_data_length,h); title({'# Active';'Confirmed Infections'}); xlim([date_times(1) date_times(end)]);
    hold on; scatter(date_times(1:length(Confirmed)),Hospitalized_active+Not_Hospitalized-y(1:length(Confirmed),8),'LineWidth',1)
h=subplot(4,5,6); plotwithCI(date_times(2:end),diff(sum(y(:,[5,6,8,9,10]),2)),bootyCI(:,[5,6,8,9,10]),real_data_length,h,'diff'); title({'# New';'Confirmed Infections'}); xlim([date_times(1) date_times(end)]);
    hold on; scatter(date_times(2:length(Not_Hospitalized)),diff(Confirmed),'LineWidth',1)
h=subplot(4,5,7); plotwithCI(date_times,sum(y(:,[4,7]),2),bootyCI(:,[4,7]),real_data_length,h); title({'# Total';'Undocumented Infections'}); xlim([date_times(1) date_times(end)]);
h=subplot(4,5,8); plotwithCI(date_times,y(:,4),bootyCI(:,4),real_data_length,h); title({'# Active';'Undocumented Infections'}); xlim([date_times(1) date_times(end)]);
h=subplot(4,5,9); plotwithCI(date_times,sum(y(:,[5,8]),2),bootyCI(:,[5,8]),real_data_length,h); title({'# Total Confirmed';'Non-Hospitalized'}); xlim([date_times(1) date_times(end)]);
    hold on; scatter(date_times(1:length(Not_Hospitalized)),Not_Hospitalized,'LineWidth',1)
h=subplot(4,5,10); plotwithCI(date_times,y(:,5),bootyCI(:,5),real_data_length,h); title({'# Confirmed Non-Hospitalized';'Active Infections'}); xlim([date_times(1) date_times(end)]);
h=subplot(4,5,11); plotwithCI(date_times,y(:,8),bootyCI(:,8),real_data_length,h); title({'# Confirmed Non-Hospitalized';'Recoveries'}); xlim([date_times(1) date_times(end)]);
h=subplot(4,5,12); plotwithCI(date_times,sum(y(:,[6,9,10]),2),bootyCI(:,[6,9,10]),real_data_length,h); title({'# Total Confirmed';'Hospitalized'}); xlim([date_times(1) date_times(end)]);
    hold on; scatter(date_times(1:length(Hospitalized_tot)),Hospitalized_active+Deaths+Discharged,'LineWidth',1)
h=subplot(4,5,13); plotwithCI(date_times,y(:,6),bootyCI(:,6),real_data_length,h); title('# Active Hospitalized'); xlim([date_times(1) date_times(end)]);
    hold on; scatter(date_times(1:length(Hospitalized_active)),Hospitalized_active,'LineWidth',1)
h=subplot(4,5,14); plotwithCI(date_times,y(:,9),bootyCI(:,9),real_data_length,h); title('# Hospital Discharges'); xlim([date_times(1) date_times(end)]);
    hold on; scatter(date_times(1:length(Discharged)),Discharged,'LineWidth',1)
h=subplot(4,5,15); rd2 = plotwithCI(date_times(2:end),(diff(y(:,9))),bootyCI(:,9),real_data_length,h,'diff'); title('# New Discharges'); xlim([date_times(1) date_times(end)]);
    hold on; rd1 = scatter(date_times(2:length(Discharged)),diff(Discharged),'LineWidth',1);
h=subplot(4,5,16); plotwithCI(date_times,y(:,10),bootyCI(:,10),real_data_length,h); title('# Deaths'); xlim([date_times(1) date_times(end)]);
    hold on; scatter(date_times(1:length(Deaths)),Deaths,'LineWidth',1,'MarkerEdgeColor','#D95319')

h=subplot(4,5,17); plotwithCI(date_times(2:end),(diff(y(:,10))),bootyCI(:,10),real_data_length,h,'diff'); title('# New Deaths'); xlim([date_times(1) date_times(end)]);
    hold on; scatter(date_times(2:length(Deaths)),diff(Deaths),'LineWidth',1)
ax = subplot(4,5,18); plot(date_times,y(:,10)./(sum([y(:,9) y(:,10)],2)),'LineWidth',2);
    hold on; plot(date_times,y(:,9)./(sum([y(:,10) y(:,9)],2)),'LineWidth',2); title({'Instantaneous Hospital';'Recovery and Death Rate'}); xlim([date_times(1) date_times(end)]);
    legend(ax,{'Death Rate', 'Recovery Rate'});
subplot(4,5,19); plot(date_times,((y(:,10)./(sum([y(:,10) y(:,9)],2)))./(y(:,9)./(sum([y(:,10) y(:,9)],2)))),'LineWidth',2);
    title('Ratio Death Rate:Recovery Rate'); xlim([date_times(1) date_times(end)]); 
temp_ax = subplot(4,5,20);set(temp_ax,'Color','none','Visible','off');
    legend(temp_ax,[rd2 rd1],{'Simulated','Real'},'Location','best','FontSize',12)
    
DR = y(end,10)/(sum(y(end,7:10)));

disp(['Death rate is: ', num2str(DR)]);
disp(['Peak R0 is: ', num2str(R_naught_peak)]);
% disp(['Recovery rate is: ', num2str(RR)]);

% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'# of people', 'FontSize', 16, 'FontWeight','bold');
xlabel(han,'time (d)','FontSize', 16, 'FontWeight','bold');

%% ##### SEIAQHRRRDP Model ###### 
% Figure 1
% Summary Outputs
date_times = Time(1):(Time(1)+t(end)-1);
fig=figure;
h=subplot(3,5,1); plotwithCI(date_times,sum(y(:,4:10),2),bootyCI(:,4:10),0,h); title({'# Total Infections';'(Undocumented + confirmed)'}); ylim([0 20e5]); xlim([date_times(1) date_times(end)]);
h=subplot(3,5,2); plotwithCI(date_times,sum(y(:,[4,7]),2),bootyCI(:,[4,7]),0,h);  title({'# Total';'Undocumented Infections'}); ylim([0 20e5]);xlim([date_times(1) date_times(end)]);
h=subplot(3,5,3); plotwithCI(date_times,sum(y(:,[5,6,8:10]),2),bootyCI(:,[5,6,8:10]),real_data_length,h); title({'# Total';'Confirmed Infections'}); ylim([0 20e5]); xlim([date_times(1) date_times(end)]);
    hold on; scatter(date_times(1:length(Confirmed)),Confirmed,'LineWidth',1, 'MarkerEdgeColor', '#D95319')
h=subplot(3,5,4); plotwithCI(date_times,sum(y(:,5:6),2),bootyCI(:,5:6),real_data_length,h); title({'# Active';'Confirmed Infections'}); xlim([date_times(1) date_times(end)]); ylim([0 10e4])
    hold on; scatter(date_times(1:length(Confirmed)),Hospitalized_active+Not_Hospitalized-y(1:length(Confirmed),8),'LineWidth',1, 'MarkerEdgeColor', '#D95319')
h=subplot(3,5,5); plotwithCI(date_times(2:end),diff(sum(y(:,[5,6,8,9,10]),2)),bootyCI(:,[5,6,8,9,10]),real_data_length,h,'diff'); title({'# New';'Confirmed Infections'}); xlim([date_times(1) date_times(end)]);ylim([0 12e3])
    hold on; scatter(date_times(2:length(Not_Hospitalized)),diff(Confirmed),'LineWidth',1, 'MarkerEdgeColor', '#D95319')
h=subplot(3,5,6); plotwithCI(date_times,sum(y(:,[6,9,10]),2),bootyCI(:,[6,9,10]),real_data_length,h); title({'# Total Confirmed';'Hospitalized'}); xlim([date_times(1) date_times(end)]);
    hold on; scatter(date_times(1:length(Hospitalized_tot)),Hospitalized_active+Deaths+Discharged,'LineWidth',1, 'MarkerEdgeColor', '#D95319')
h=subplot(3,5,7); plotwithCI(date_times,y(:,6),bootyCI(:,6),real_data_length,h); title('# Active Hospitalized'); xlim([date_times(1) date_times(end)]);
    hold on; scatter(date_times(1:length(Hospitalized_active)),Hospitalized_active,'LineWidth',1, 'MarkerEdgeColor', '#D95319')
h=subplot(3,5,8); plotwithCI(date_times,y(:,9),bootyCI(:,9),real_data_length,h); title('# Hospital Discharges'); xlim([date_times(1) date_times(end)]);
    hold on; scatter(date_times(1:length(Discharged)),Discharged,'LineWidth',1, 'MarkerEdgeColor', '#D95319')
h=subplot(3,5,9); rd2 = plotwithCI(date_times(2:end),(diff(y(:,9))),bootyCI(:,9),real_data_length,h,'diff'); title('# New Discharges'); xlim([date_times(1) date_times(end)]);ylim([0 2500])
    hold on; rd1 = scatter(date_times(2:length(Discharged)),diff(Discharged),'LineWidth',1, 'MarkerEdgeColor', '#D95319');
h=subplot(3,5,10); plotwithCI(date_times,y(:,10),bootyCI(:,10),real_data_length,h); title('# Deaths'); xlim([date_times(1) date_times(end)]);
    hold on; scatter(date_times(1:length(Deaths)),Deaths,'LineWidth',1, 'MarkerEdgeColor', '#D95319')
h=subplot(3,5,11); plotwithCI(date_times(2:end),(diff(y(:,10))),bootyCI(:,10),real_data_length,h,'diff'); title('# New Deaths'); xlim([date_times(1) date_times(end)]); ylim([0 800])
    hold on; scatter(date_times(2:length(Deaths)),diff(Deaths),'LineWidth',1, 'MarkerEdgeColor', '#D95319')
subplot(3,5,12); plot(date_times(1:end),R_naught,'LineWidth',2); title({'Basic';'Reproductive #'}); xlim([date_times(1) date_times(end)]);
subplot(3,5,13); plot(date_times(1:end),R_naught'.*y(:,1)/Npop,'LineWidth',2); title({'Effective';'Reproductive #'}); xlim([date_times(1) date_times(end)]);

temp_ax = subplot(3,5,14); set(temp_ax,'xcolor','none'); set(temp_ax,'ycolor','none')
     text(0.1,0.75,{['Peak Reproductive Number:    ',num2str(round(R_naught_peak,2))];...
                    ['Undocumented Rate:        ',num2str(Asymp_rate*100),'%'];...
                    ['Total Infected:                 ',num2str(round(max(sum(y(:,4:10),2))))];...
                    ['Total Dead:                      ',num2str(round(max(y(:,10))))];...
                    ['Mortality Rate:                 ',num2str(round(DR*100,2)),'%'];...
         },'Units','Normalized','FontSize', 12);
temp_ax = subplot(3,5,15);set(temp_ax,'Color','none','Visible','off');
    legend(temp_ax,[rd2 rd1],{'Simulated','Real'},'Location','best','FontSize',12)
  
DR = y(end,10)/(sum(y(end,7:10)));
RR = y(end,7:9)/(sum(y(end,7:10)));

disp(['Death rate is: ', num2str(DR)]);
disp(['R0 is: ', num2str(R_naught_peak)]);

% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'# of people', 'FontSize', 16, 'FontWeight','bold');
sgtitle({'BUHoffman-COVID-Model, New York:';'100% Sustained Immunity, Permanent Social Distancing'},'FontSize', 16, 'FontWeight','bold')

%% Figure 2, Figure S6-7
% Effect on Time to steady state social distancing and steady state level
% total deaths by 9/1/2020

% Immunity
Immunity_fraction = 0.70;% fraction of recovered who develop immunity
tau = 1-Immunity_fraction; % fraction of recovered who DO NOT develop immunity
pi = 1/1825; % 1/(Duration of immunity)

SD_remove = cell(1,20);
for i = 1:20
    time_factor = 5*i;
    SDout = ones(1,3);
    SDout(1,1) = 40+time_factor;
    SDout(1,2:3) = 1;
    SD_remove{1,i} = SDout;
end

zeta_factor = cell(1,20);
for i = 1:20
    relax_factor = (1/20)*i;
    zetaout = ones(1,3);
    for j = 1:size(zetaout,2)
        zetaout(1,j) = 1/(1.05-relax_factor);
    end
    zeta_factor{1,i} = zetaout;
end

% WINDOW FOR SIMULATION
window_sim = [1, 182];
chi = importflu(Time,window_sim,'factor',flufactor);
yout = cell(20,20);
tic
for i = 1:20
    for j = 1:20
        [t,y] = ode45(@(t,y) SEIAQHRRRDP_deqs(t, y, Coef, Npop, rho, gamma,...
            nu, phi, delta, theta, omega,tau, pi, SD_delay, ...
            SD_remove{1,j}, zeta_factor{1,i}, window_sim, chi), window_sim(1):window_sim(end), ic);
        yout{i,j} = y;
    end
end
toc

time2SS_SD = zeros(20,1);
for i = 1:20
    time2SS_SD(i,1) = sum(SD_remove{1,i}(:,:));
end
SS_SD = zeros(20,1);
for i = 1:20
    SS_SD(i,1) = 100-(100/(zeta_factor{1,21-i}(3)));
end

maxdeaths = zeros(20,20);
for i = 1:20
    for j = 1:20
        maxdeaths(i,j) = yout{i,j}(end,10);
    end
end
maxdeaths = flipud(maxdeaths);

R_naught = ((chi*rho*Coef(3))/(nu+theta+delta)) + ((chi*rho*Coef(3)*nu)/(phi*(nu+theta+delta)));
R_naught_peak = max(((chi*rho*Coef(3))/(nu+theta+delta)) + ((chi*rho*Coef(3)*nu)/(phi*(nu+theta+delta))));

ReThresh = zeros(20,20);
for i = 1:20
    for j = 1:20
        Reff = (R_naught'.*yout{i,j}(:,1)/Npop);
        above = find( Reff > 1);
        above2 = find( above > (sum(SD_remove{1,j})+ SD_delay));
        if ~isempty(above2)
%             ReThresh(i,j) = above(above2(1));
            ReThresh(i,j) = 1;
        end
    end
end
ReThresh = flipud(ReThresh);

figure;
subplot(2,2,1);heatmap(date_times(1)+SD_delay+time2SS_SD,SS_SD,maxdeaths,...
    'GridVisible','off','CellLabelColor','none','Colormap',parula,...
    'XLabel','Date of Social Distancing Relaxation',...
    'YLabel','% Reduction in Social Distancing')

subplot(2,2,2);heatmap(date_times(1)+SD_delay+time2SS_SD,SS_SD,ReThresh,...
    'GridVisible','off','CellLabelColor','none','Colormap',[0.2422,0.1504,0.6603;0.9769,0.9839,0.0805],...
    'XLabel','Date of Social Distancing Relaxation',...
    'YLabel','% Reduction in Social Distancing')

%% Fig 2  A B C E
y1 = yout{4,6};
y2 = yout{7,6};
y3 = yout{11,6};

zeta_factor1 = zeta_factor{:,4};
zeta_factor2 = zeta_factor{:,7};
zeta_factor3 = zeta_factor{:,11};

SD_remove1 = SD_remove{:,6};
SD_remove2 = SD_remove{:,6};
SD_remove3 = SD_remove{:,6};


[lambda_vector,~] = lambdakappa(t,Coef);
Total_infxns1 = (cumtrapz(y1(:,3))*nu) + (cumtrapz(y1(:,3))*theta) + (cumtrapz(y1(:,3))*delta);
Total_infxns2 = (cumtrapz(y2(:,3))*nu) + (cumtrapz(y2(:,3))*theta) + (cumtrapz(y2(:,3))*delta);
Total_infxns3 = (cumtrapz(y3(:,3))*nu) + (cumtrapz(y3(:,3))*theta) + (cumtrapz(y3(:,3))*delta);

Total_confirmed1 = (cumtrapz(y1(:,3))*theta) + (cumtrapz(y1(:,3))*delta);
Total_confirmed2 = (cumtrapz(y2(:,3))*theta) + (cumtrapz(y2(:,3))*delta);
Total_confirmed3 = (cumtrapz(y3(:,3))*theta) + (cumtrapz(y3(:,3))*delta);

% Summary Outputs
All_recoveries1 = (cumtrapz(y1(:,4))*phi) + (cumtrapz(y1(:,5))*omega) + (cumtrapz(y1(:,6)).*lambda_vector);

DR1 = y(end,10)/(All_recoveries(end)+y1(end,10));

date_times = Time(1):(Time(1)+t(end)-1);
fig=figure;
subplot(3,5,1); plot(date_times(1:end),Total_infxns1,'LineWidth',2);
    hold on; title({'# Total Infections';'(Undocumented + confirmed)'}); ylim([0 max(max([Total_infxns1,Total_infxns2,Total_infxns3]))]); xlim([date_times(1) date_times(end)]);
    hold on;plot(date_times(1:end),Total_infxns2,'LineWidth',2)
    plot(date_times(1:end),Total_infxns3,'LineWidth',2)
    
subplot(3,5,2); plot(date_times,(cumtrapz(y1(:,3))*nu),'LineWidth',2); title({'# Total';'Undocumented Infections'}); ylim([0 max(Total_infxns)]); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,(cumtrapz(y2(:,3))*nu),'LineWidth',2);
    plot(date_times,(cumtrapz(y3(:,3))*nu),'LineWidth',2);

subplot(3,5,3); plot(date_times,Total_confirmed1,'LineWidth',2); title({'# Total';'Confirmed Infections'}); title({'# Total';'Confirmed Infections'}); ylim([0 max(Total_infxns)]); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,Total_confirmed2,'LineWidth',2);
    plot(date_times,Total_confirmed3,'LineWidth',2);
    scatter(date_times(1:length(Confirmed)),Confirmed,'LineWidth',1,'MarkerEdgeColor','#D95319')

subplot(3,5,4); plot(date_times,sum(y1(:,5:6),2),'LineWidth',2); title({'# Active';'Confirmed Infections'}); title({'# Active';'Confirmed Infections'}); xlim([date_times(1) date_times(end)]);
    hold on;plot(date_times,sum(y2(:,5:6),2),'LineWidth',2);
    plot(date_times,sum(y3(:,5:6),2),'LineWidth',2);
    scatter(date_times(1:length(Confirmed)),Hospitalized_active+Not_Hospitalized-y(1:length(Confirmed),8),'LineWidth',1, 'MarkerEdgeColor','#D95319')
    
subplot(3,5,5); plot(date_times(2:end),diff(Total_confirmed1),'LineWidth',2); title({'# New';'Confirmed Infections'}); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times(2:end),diff(Total_confirmed2),'LineWidth',2);
    plot(date_times(2:end),diff(Total_confirmed3),'LineWidth',2);
    scatter(date_times(2:length(Not_Hospitalized)),diff(Confirmed),'LineWidth',1, 'MarkerEdgeColor','#D95319')
    
subplot(3,5,6); plot(date_times,(cumtrapz(y1(:,3))*delta),'LineWidth',2); title({'# Total Confirmed';'Hospitalized'}); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,(cumtrapz(y2(:,3))*delta),'LineWidth',2);
    plot(date_times,(cumtrapz(y3(:,3))*delta),'LineWidth',2);
    scatter(date_times(1:length(Hospitalized_tot)),Hospitalized_active+Deaths+Discharged,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,7); plot(date_times,y1(:,6),'LineWidth',2); title('# Active Hospitalized'); xlim([date_times(1) date_times(end)]);
    hold on;  plot(date_times,y2(:,6),'LineWidth',2);
    plot(date_times,y3(:,6),'LineWidth',2);
    scatter(date_times(1:length(Hospitalized_active)),Hospitalized_active,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,8); plot(date_times,(cumtrapz(y1(:,6)).*lambda_vector),'LineWidth',2); title({'# Total';'Hospital Discharges'}); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,(cumtrapz(y2(:,6)).*lambda_vector),'LineWidth',2);
    plot(date_times,(cumtrapz(y3(:,6)).*lambda_vector),'LineWidth',2);
    scatter(date_times(1:length(Discharged)),Discharged,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,9); rd1 = plot(date_times(2:end),diff((cumtrapz(y1(:,6)).*lambda_vector)),'LineWidth',2); title('# New Discharges'); xlim([date_times(1) date_times(end)]);
    hold on; rd2 = plot(date_times(2:end),diff((cumtrapz(y2(:,6)).*lambda_vector)),'LineWidth',2);
    rd3 = plot(date_times(2:end),diff((cumtrapz(y3(:,6)).*lambda_vector)),'LineWidth',2);
    rd4 = scatter(date_times(2:length(Discharged)),diff(Discharged),'LineWidth',1,'MarkerEdgeColor','#D95319');
    
subplot(3,5,10); plot(date_times,y1(:,10),'LineWidth',2); title('# Deaths'); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,y2(:,10),'LineWidth',2);
    plot(date_times,y3(:,10),'LineWidth',2);
    scatter(date_times(1:length(Deaths)),Deaths,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,11); plot(date_times(2:end),(diff(y1(:,10))),'LineWidth',2); title('# New Deaths'); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times(2:end),(diff(y2(:,10))),'LineWidth',2);
    plot(date_times(2:end),(diff(y3(:,10))),'LineWidth',2);
    scatter(date_times(2:length(Deaths)),diff(Deaths),'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,12); plot(date_times(1:end),R_naught'.*y1(:,1)/Npop,'LineWidth',2); title({'Effective';'Reproductive #'}); xlim([date_times(1) date_times(end)]); ylim([0 2]);
    hold on; plot(date_times(1:end),R_naught'.*y2(:,1)/Npop,'LineWidth',2);
    plot(date_times(1:end),R_naught'.*y3(:,1)/Npop,'LineWidth',2);

temp_ax = subplot(3,5,13); set(temp_ax,'xcolor','none'); set(temp_ax,'ycolor','none')
     text(-0.10,0.6,{...
                    'Summary Values';...
                    ['  Peak Reproductive Number:     ',num2str(round(R_naught_peak,2))];...
                    ['  Undocumented Rate:              ',num2str(Asymp_rate*100),'%'];...
                    ['  SD Started:                            ',datestr(date_times(SD_delay+1),'mm/dd/yy')];...
                    ['  Mortality Rate:                       ',num2str(round(DR1*100,2)),'%'];...
                    ' ';...
                    'Return To Society Model 1';
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor1(1))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove1(1)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor1(2))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove1(1:2)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor1(3))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove1(1:3)),'mm/dd/yy')];...
                    ['  Total Infected:                       ',num2str(round(max(All_recoveries1+y1(:,10))))];...
                    ['  Total Dead:                            ',num2str(round(max(y1(:,10))))];...
                    ['  % Increase in Deaths:           ',num2str(round(100*((max(y1(:,10))) - (max(y1(:,10))))/((max(y1(:,10)))),1)),'%'];...
                     },'Units','Normalized','FontSize', 12);
 
temp_ax = subplot(3,5,14); set(temp_ax,'xcolor','none'); set(temp_ax,'ycolor','none')
     text(-0.10,0.5,{...     
                    'Return To Society Model 2';
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor2(1))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove2(1)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor2(2))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove2(1:2)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor2(3))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove2(1:3)),'mm/dd/yy')];...
                    ['  Total Infected:                       ',num2str(round(max(All_recoveries1+y2(:,10))))];...
                    ['  Total Dead:                            ',num2str(round(max(y2(:,10))))];...
                    ['  % Increase in Deaths:           ',num2str(round(100*((max(y2(:,10))) - (max(y1(:,10))))/((max(y2(:,10)))),1)),'%'];...
                    '';
                    'Return To Society Model 3';
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor3(1))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove3(1)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor3(2))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove3(1:2)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor3(3))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove3(1:3)),'mm/dd/yy')];...
                    ['  Total Infected:                       ',num2str(round(max(All_recoveries1+y3(:,10))))];...
                    ['  Total Dead:                            ',num2str(round(max(y3(:,10))))];...
                    ['  % Increase in Deaths:           ',num2str(round(100*((max(y3(:,10))) - (max(y1(:,10))))/((max(y3(:,10)))),1)),'%'];...
                     },'Units','Normalized','FontSize', 12);     
                
temp_ax = subplot(3,5,15);set(temp_ax,'Color','none','Visible','off');
    legend(temp_ax,[rd1 rd2 rd3 rd4],...
        {'Return to Society Model 1',...
        'Return to Society Model 2',...
        'Return to Society Model 3',...
        'Real Data'},'Location','best','FontSize',12)
    
disp(['Death rate is: ', num2str(DR)]);
disp(['Peak R0 is: ', num2str(R_naught_peak)]);

% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'# of people', 'FontSize', 16, 'FontWeight','bold');
sgtitle({'BUHoffman-COVID-Model, New York:';'Medium-term Effects of Multiple "Return of Society" Models of Social Distancing Relaxation';'90% Immunity (Sustained 365 days)'},'FontSize', 16, 'FontWeight','bold')
set(gcf,'color','w');

%% Fig 2 G H
y1 = yout{6,6};
y2 = yout{11,6};
y3 = yout{16,6};

zeta_factor1 = zeta_factor{:,6};
zeta_factor2 = zeta_factor{:,11};
zeta_factor3 = zeta_factor{:,16};

SD_remove1 = SD_remove{:,6};
SD_remove2 = SD_remove{:,6};
SD_remove3 = SD_remove{:,6};


[lambda_vector,~] = lambdakappa(t,Coef);
Total_infxns1 = (cumtrapz(y1(:,3))*nu) + (cumtrapz(y1(:,3))*theta) + (cumtrapz(y1(:,3))*delta);
Total_infxns2 = (cumtrapz(y2(:,3))*nu) + (cumtrapz(y2(:,3))*theta) + (cumtrapz(y2(:,3))*delta);
Total_infxns3 = (cumtrapz(y3(:,3))*nu) + (cumtrapz(y3(:,3))*theta) + (cumtrapz(y3(:,3))*delta);

Total_confirmed1 = (cumtrapz(y1(:,3))*theta) + (cumtrapz(y1(:,3))*delta);
Total_confirmed2 = (cumtrapz(y2(:,3))*theta) + (cumtrapz(y2(:,3))*delta);
Total_confirmed3 = (cumtrapz(y3(:,3))*theta) + (cumtrapz(y3(:,3))*delta);

% Summary Outputs
All_recoveries1 = (cumtrapz(y1(:,4))*phi) + (cumtrapz(y1(:,5))*omega) + (cumtrapz(y1(:,6)).*lambda_vector);

DR1 = y(end,10)/(All_recoveries(end)+y1(end,10));

date_times = Time(1):(Time(1)+t(end)-1);
fig=figure;
subplot(3,5,1); plot(date_times(1:end),Total_infxns1,'LineWidth',2);
    hold on; title({'# Total Infections';'(Undocumented + confirmed)'}); ylim([0 max(max([Total_infxns1,Total_infxns2,Total_infxns3]))]); xlim([date_times(1) date_times(end)]);
    hold on;plot(date_times(1:end),Total_infxns2,'LineWidth',2)
    plot(date_times(1:end),Total_infxns3,'LineWidth',2)
    
subplot(3,5,2); plot(date_times,(cumtrapz(y1(:,3))*nu),'LineWidth',2); title({'# Total';'Undocumented Infections'}); ylim([0 max(Total_infxns)]); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,(cumtrapz(y2(:,3))*nu),'LineWidth',2);
    plot(date_times,(cumtrapz(y3(:,3))*nu),'LineWidth',2);

subplot(3,5,3); plot(date_times,Total_confirmed1,'LineWidth',2); title({'# Total';'Confirmed Infections'}); title({'# Total';'Confirmed Infections'}); ylim([0 max(Total_infxns)]); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,Total_confirmed2,'LineWidth',2);
    plot(date_times,Total_confirmed3,'LineWidth',2);
    scatter(date_times(1:length(Confirmed)),Confirmed,'LineWidth',1,'MarkerEdgeColor','#D95319')

subplot(3,5,4); plot(date_times,sum(y1(:,5:6),2),'LineWidth',2); title({'# Active';'Confirmed Infections'}); title({'# Active';'Confirmed Infections'}); xlim([date_times(1) date_times(end)]);
    hold on;plot(date_times,sum(y2(:,5:6),2),'LineWidth',2);
    plot(date_times,sum(y3(:,5:6),2),'LineWidth',2);
    scatter(date_times(1:length(Confirmed)),Hospitalized_active+Not_Hospitalized-y(1:length(Confirmed),8),'LineWidth',1, 'MarkerEdgeColor','#D95319')
    
subplot(3,5,5); plot(date_times(2:end),diff(Total_confirmed1),'LineWidth',2); title({'# New';'Confirmed Infections'}); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times(2:end),diff(Total_confirmed2),'LineWidth',2);
    plot(date_times(2:end),diff(Total_confirmed3),'LineWidth',2);
    scatter(date_times(2:length(Not_Hospitalized)),diff(Confirmed),'LineWidth',1, 'MarkerEdgeColor','#D95319')
    
subplot(3,5,6); plot(date_times,(cumtrapz(y1(:,3))*delta),'LineWidth',2); title({'# Total Confirmed';'Hospitalized'}); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,(cumtrapz(y2(:,3))*delta),'LineWidth',2);
    plot(date_times,(cumtrapz(y3(:,3))*delta),'LineWidth',2);
    scatter(date_times(1:length(Hospitalized_tot)),Hospitalized_active+Deaths+Discharged,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,7); plot(date_times,y1(:,6),'LineWidth',2); title('# Active Hospitalized'); xlim([date_times(1) date_times(end)]);
    hold on;  plot(date_times,y2(:,6),'LineWidth',2);
    plot(date_times,y3(:,6),'LineWidth',2);
    scatter(date_times(1:length(Hospitalized_active)),Hospitalized_active,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,8); plot(date_times,(cumtrapz(y1(:,6)).*lambda_vector),'LineWidth',2); title({'# Total';'Hospital Discharges'}); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,(cumtrapz(y2(:,6)).*lambda_vector),'LineWidth',2);
    plot(date_times,(cumtrapz(y3(:,6)).*lambda_vector),'LineWidth',2);
    scatter(date_times(1:length(Discharged)),Discharged,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,9); rd1 = plot(date_times(2:end),diff((cumtrapz(y1(:,6)).*lambda_vector)),'LineWidth',2); title('# New Discharges'); xlim([date_times(1) date_times(end)]);
    hold on; rd2 = plot(date_times(2:end),diff((cumtrapz(y2(:,6)).*lambda_vector)),'LineWidth',2);
    rd3 = plot(date_times(2:end),diff((cumtrapz(y3(:,6)).*lambda_vector)),'LineWidth',2);
    rd4 = scatter(date_times(2:length(Discharged)),diff(Discharged),'LineWidth',1,'MarkerEdgeColor','#D95319');
    
subplot(3,5,10); plot(date_times,y1(:,10),'LineWidth',2); title('# Deaths'); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,y2(:,10),'LineWidth',2);
    plot(date_times,y3(:,10),'LineWidth',2);
    scatter(date_times(1:length(Deaths)),Deaths,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,11); plot(date_times(2:end),(diff(y1(:,10))),'LineWidth',2); title('# New Deaths'); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times(2:end),(diff(y2(:,10))),'LineWidth',2);
    plot(date_times(2:end),(diff(y3(:,10))),'LineWidth',2);
    scatter(date_times(2:length(Deaths)),diff(Deaths),'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,12); plot(date_times(1:end),R_naught'.*y1(:,1)/Npop,'LineWidth',2); title({'Effective';'Reproductive #'}); xlim([date_times(1) date_times(end)]); ylim([0 2]);
    hold on; plot(date_times(1:end),R_naught'.*y2(:,1)/Npop,'LineWidth',2);
    plot(date_times(1:end),R_naught'.*y3(:,1)/Npop,'LineWidth',2);

temp_ax = subplot(3,5,13); set(temp_ax,'xcolor','none'); set(temp_ax,'ycolor','none')
     text(-0.10,0.6,{...
                    'Summary Values';...
                    ['  Peak Reproductive Number:     ',num2str(round(R_naught_peak,2))];...
                    ['  Undocumented Rate:              ',num2str(Asymp_rate*100),'%'];...
                    ['  SD Started:                            ',datestr(date_times(SD_delay+1),'mm/dd/yy')];...
                    ['  Mortality Rate:                       ',num2str(round(DR1*100,2)),'%'];...
                    ' ';...
                    'Return To Society Model 1';
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor1(1))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove1(1)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor1(2))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove1(1:2)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor1(3))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove1(1:3)),'mm/dd/yy')];...
                    ['  Total Infected:                       ',num2str(round(max(All_recoveries1+y1(:,10))))];...
                    ['  Total Dead:                            ',num2str(round(max(y1(:,10))))];...
                    ['  % Increase in Deaths:           ',num2str(round(100*((max(y1(:,10))) - (max(y1(:,10))))/((max(y1(:,10)))),1)),'%'];...
                     },'Units','Normalized','FontSize', 12);
 
temp_ax = subplot(3,5,14); set(temp_ax,'xcolor','none'); set(temp_ax,'ycolor','none')
     text(-0.10,0.5,{...     
                    'Return To Society Model 2';
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor2(1))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove2(1)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor2(2))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove2(1:2)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor2(3))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove2(1:3)),'mm/dd/yy')];...
                    ['  Total Infected:                       ',num2str(round(max(All_recoveries1+y2(:,10))))];...
                    ['  Total Dead:                            ',num2str(round(max(y2(:,10))))];...
                    ['  % Increase in Deaths:           ',num2str(round(100*((max(y2(:,10))) - (max(y1(:,10))))/((max(y2(:,10)))),1)),'%'];...
                    '';
                    'Return To Society Model 3';
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor3(1))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove3(1)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor3(2))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove3(1:2)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor3(3))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove3(1:3)),'mm/dd/yy')];...
                    ['  Total Infected:                       ',num2str(round(max(All_recoveries1+y3(:,10))))];...
                    ['  Total Dead:                            ',num2str(round(max(y3(:,10))))];...
                    ['  % Increase in Deaths:           ',num2str(round(100*((max(y3(:,10))) - (max(y1(:,10))))/((max(y3(:,10)))),1)),'%'];...
                     },'Units','Normalized','FontSize', 12);     
                
temp_ax = subplot(3,5,15);set(temp_ax,'Color','none','Visible','off');
    legend(temp_ax,[rd1 rd2 rd3 rd4],...
        {'Return to Society Model 1',...
        'Return to Society Model 2',...
        'Return to Society Model 3',...
        'Real Data'},'Location','best','FontSize',12)
    
disp(['Death rate is: ', num2str(DR)]);
disp(['Peak R0 is: ', num2str(R_naught_peak)]);

% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'# of people', 'FontSize', 16, 'FontWeight','bold');
sgtitle({'BUHoffman-COVID-Model, New York:';'Medium-term Effects of Multiple "Return of Society" Models of Social Distancing Relaxation';'90% Immunity (Sustained 365 days)'},'FontSize', 16, 'FontWeight','bold')
set(gcf,'color','w');


%% Figure S4
% Effect on Time to steady state social distancing and steady state level
% total deaths by 9/1/2020

% Immunity
Immunity_fraction = 0.70;% fraction of recovered who develop immunity
tau = 1-Immunity_fraction; % fraction of recovered who DO NOT develop immunity
pi = 1/1825; % 1/(Duration of immunity)

SD_remove = cell(1,20);
for i = 1:20
    time_factor = 5*i;
    SDout = ones(1,3);
    SDout(1,1) = 40+time_factor;
    SDout(1,2:3) = 14;
    SD_remove{1,i} = SDout;
end

zeta_factor = cell(1,20);
for i = 1:20
    zetaout = ones(1,3);
    if i == 1
        zetaout(1,1) = 1/(1.05-((1/20)*i));
        zetaout(1,2) = 1/(1.05-((1/20)*i));
        zetaout(1,3) = 1/(1.05-((1/20)*i));
    else
        zetaout(1,1) = 1/(1.05-(((1/20)*i))+(i/30));
        zetaout(1,2) = 1/(1.05-(((1/20)*i))+(i/60));
        zetaout(1,3) = 1/(1.05-((1/20)*i));
    end
    zeta_factor{1,i} = zetaout;
end

% WINDOW FOR SIMULATION
window_sim = [1, 187];% 182
chi = importflu(Time,window_sim,'factor',flufactor);
yout = cell(20,20);
tic
for i = 1:20
    for j = 1:20
        [t,y] = ode45(@(t,y) SEIAQHRRRDP_deqs(t, y, Coef, Npop, rho, gamma,...
            nu, phi, delta, theta, omega,tau, pi, SD_delay, ...
            SD_remove{1,j}, zeta_factor{1,i}, window_sim, chi), window_sim(1):window_sim(end), ic);
        yout{i,j} = y;
    end
end
toc

time2SS_SD = zeros(20,1);
for i = 1:20
    time2SS_SD(i,1) = sum(SD_remove{1,i}(:,:));
end
SS_SD = zeros(20,1);
for i = 1:20
    SS_SD(i,1) = 100-(100/(zeta_factor{1,21-i}(3)));
end

maxdeaths = zeros(20,20);
for i = 1:20
    for j = 1:20
        maxdeaths(i,j) = yout{i,j}(end,10);
    end
end
maxdeaths = flipud(maxdeaths);

R_naught = ((chi*rho*Coef(3))/(nu+theta+delta)) + ((chi*rho*Coef(3)*nu)/(phi*(nu+theta+delta)));
R_naught_peak = max(((chi*rho*Coef(3))/(nu+theta+delta)) + ((chi*rho*Coef(3)*nu)/(phi*(nu+theta+delta))));

ReThresh = zeros(20,20);
for i = 1:20
    for j = 1:20
        Reff = (R_naught'.*yout{i,j}(:,1)/Npop);
        above = find( Reff > 1);
        above2 = find( above > (sum(SD_remove{1,j})+ SD_delay));
        if ~isempty(above2)
            ReThresh(i,j) = 1;
        end
    end
end
ReThresh = flipud(ReThresh);

figure;
subplot(2,2,1);heatmap(date_times(1)+SD_delay+time2SS_SD,SS_SD,maxdeaths,...
    'GridVisible','off','CellLabelColor','none','Colormap',parula,...
    'XLabel','Date of Social Distancing Relaxation',...
    'YLabel','% Reduction in Social Distancing')

subplot(2,2,2);heatmap(date_times(1)+SD_delay+time2SS_SD,SS_SD,ReThresh,...
    'GridVisible','off','CellLabelColor','none','Colormap',[0.2422,0.1504,0.6603;0.9769,0.9839,0.0805],...
    'XLabel','Date of Social Distancing Relaxation',...
    'YLabel','% Reduction in Social Distancing')

%% Figure S3 A B C E

y1 = yout{4,6};
y2 = yout{7,6};
y3 = yout{11,6};

zeta_factor1 = zeta_factor{:,4};
zeta_factor2 = zeta_factor{:,7};
zeta_factor3 = zeta_factor{:,11};

SD_remove1 = SD_remove{:,6};
SD_remove2 = SD_remove{:,6};
SD_remove3 = SD_remove{:,6};


[lambda_vector,~] = lambdakappa(t,Coef);
Total_infxns1 = (cumtrapz(y1(:,3))*nu) + (cumtrapz(y1(:,3))*theta) + (cumtrapz(y1(:,3))*delta);
Total_infxns2 = (cumtrapz(y2(:,3))*nu) + (cumtrapz(y2(:,3))*theta) + (cumtrapz(y2(:,3))*delta);
Total_infxns3 = (cumtrapz(y3(:,3))*nu) + (cumtrapz(y3(:,3))*theta) + (cumtrapz(y3(:,3))*delta);

Total_confirmed1 = (cumtrapz(y1(:,3))*theta) + (cumtrapz(y1(:,3))*delta);
Total_confirmed2 = (cumtrapz(y2(:,3))*theta) + (cumtrapz(y2(:,3))*delta);
Total_confirmed3 = (cumtrapz(y3(:,3))*theta) + (cumtrapz(y3(:,3))*delta);

% Summary Outputs
All_recoveries1 = (cumtrapz(y1(:,4))*phi) + (cumtrapz(y1(:,5))*omega) + (cumtrapz(y1(:,6)).*lambda_vector);

DR1 = y(end,10)/(All_recoveries(end)+y1(end,10));

date_times = Time(1):(Time(1)+t(end)-1);
fig=figure;
subplot(3,5,1); plot(date_times(1:end),Total_infxns1,'LineWidth',2);
    hold on; title({'# Total Infections';'(Undocumented + confirmed)'}); ylim([0 max(max([Total_infxns1,Total_infxns2,Total_infxns3]))]); xlim([date_times(1) date_times(end)]);
    hold on;plot(date_times(1:end),Total_infxns2,'LineWidth',2)
    plot(date_times(1:end),Total_infxns3,'LineWidth',2)
    
subplot(3,5,2); plot(date_times,(cumtrapz(y1(:,3))*nu),'LineWidth',2); title({'# Total';'Undocumented Infections'}); ylim([0 max(Total_infxns)]); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,(cumtrapz(y2(:,3))*nu),'LineWidth',2);
    plot(date_times,(cumtrapz(y3(:,3))*nu),'LineWidth',2);

subplot(3,5,3); plot(date_times,Total_confirmed1,'LineWidth',2); title({'# Total';'Confirmed Infections'}); title({'# Total';'Confirmed Infections'}); ylim([0 max(Total_infxns)]); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,Total_confirmed2,'LineWidth',2);
    plot(date_times,Total_confirmed3,'LineWidth',2);
    scatter(date_times(1:length(Confirmed)),Confirmed,'LineWidth',1,'MarkerEdgeColor','#D95319')

subplot(3,5,4); plot(date_times,sum(y1(:,5:6),2),'LineWidth',2); title({'# Active';'Confirmed Infections'}); title({'# Active';'Confirmed Infections'}); xlim([date_times(1) date_times(end)]);
    hold on;plot(date_times,sum(y2(:,5:6),2),'LineWidth',2);
    plot(date_times,sum(y3(:,5:6),2),'LineWidth',2);
    scatter(date_times(1:length(Confirmed)),Hospitalized_active+Not_Hospitalized-y(1:length(Confirmed),8),'LineWidth',1, 'MarkerEdgeColor','#D95319')
    
subplot(3,5,5); plot(date_times(2:end),diff(Total_confirmed1),'LineWidth',2); title({'# New';'Confirmed Infections'}); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times(2:end),diff(Total_confirmed2),'LineWidth',2);
    plot(date_times(2:end),diff(Total_confirmed3),'LineWidth',2);
    scatter(date_times(2:length(Not_Hospitalized)),diff(Confirmed),'LineWidth',1, 'MarkerEdgeColor','#D95319')
    
subplot(3,5,6); plot(date_times,(cumtrapz(y1(:,3))*delta),'LineWidth',2); title({'# Total Confirmed';'Hospitalized'}); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,(cumtrapz(y2(:,3))*delta),'LineWidth',2);
    plot(date_times,(cumtrapz(y3(:,3))*delta),'LineWidth',2);
    scatter(date_times(1:length(Hospitalized_tot)),Hospitalized_active+Deaths+Discharged,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,7); plot(date_times,y1(:,6),'LineWidth',2); title('# Active Hospitalized'); xlim([date_times(1) date_times(end)]);
    hold on;  plot(date_times,y2(:,6),'LineWidth',2);
    plot(date_times,y3(:,6),'LineWidth',2);
    scatter(date_times(1:length(Hospitalized_active)),Hospitalized_active,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,8); plot(date_times,(cumtrapz(y1(:,6)).*lambda_vector),'LineWidth',2); title({'# Total';'Hospital Discharges'}); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,(cumtrapz(y2(:,6)).*lambda_vector),'LineWidth',2);
    plot(date_times,(cumtrapz(y3(:,6)).*lambda_vector),'LineWidth',2);
    scatter(date_times(1:length(Discharged)),Discharged,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,9); rd1 = plot(date_times(2:end),diff((cumtrapz(y1(:,6)).*lambda_vector)),'LineWidth',2); title('# New Discharges'); xlim([date_times(1) date_times(end)]);
    hold on; rd2 = plot(date_times(2:end),diff((cumtrapz(y2(:,6)).*lambda_vector)),'LineWidth',2);
    rd3 = plot(date_times(2:end),diff((cumtrapz(y3(:,6)).*lambda_vector)),'LineWidth',2);
    rd4 = scatter(date_times(2:length(Discharged)),diff(Discharged),'LineWidth',1,'MarkerEdgeColor','#D95319');
    
subplot(3,5,10); plot(date_times,y1(:,10),'LineWidth',2); title('# Deaths'); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,y2(:,10),'LineWidth',2);
    plot(date_times,y3(:,10),'LineWidth',2);
    scatter(date_times(1:length(Deaths)),Deaths,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,11); plot(date_times(2:end),(diff(y1(:,10))),'LineWidth',2); title('# New Deaths'); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times(2:end),(diff(y2(:,10))),'LineWidth',2);
    plot(date_times(2:end),(diff(y3(:,10))),'LineWidth',2);
    scatter(date_times(2:length(Deaths)),diff(Deaths),'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,12); plot(date_times(1:end),R_naught'.*y1(:,1)/Npop,'LineWidth',2); title({'Effective';'Reproductive #'}); xlim([date_times(1) date_times(end)]); ylim([0 2]);
    hold on; plot(date_times(1:end),R_naught'.*y2(:,1)/Npop,'LineWidth',2);
    plot(date_times(1:end),R_naught'.*y3(:,1)/Npop,'LineWidth',2);

temp_ax = subplot(3,5,13); set(temp_ax,'xcolor','none'); set(temp_ax,'ycolor','none')
     text(-0.10,0.6,{...
                    'Summary Values';...
                    ['  Peak Reproductive Number:     ',num2str(round(R_naught_peak,2))];...
                    ['  Undocumented Rate:              ',num2str(Asymp_rate*100),'%'];...
                    ['  SD Started:                            ',datestr(date_times(SD_delay+1),'mm/dd/yy')];...
                    ['  Mortality Rate:                       ',num2str(round(DR1*100,2)),'%'];...
                    ' ';...
                    'Return To Society Model 1';
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor1(1))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove1(1)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor1(2))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove1(1:2)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor1(3))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove1(1:3)),'mm/dd/yy')];...
                    ['  Total Infected:                       ',num2str(round(max(All_recoveries1+y1(:,10))))];...
                    ['  Total Dead:                            ',num2str(round(max(y1(:,10))))];...
                    ['  % Increase in Deaths:           ',num2str(round(100*((max(y1(:,10))) - (max(y1(:,10))))/((max(y1(:,10)))),1)),'%'];...
                     },'Units','Normalized','FontSize', 12);
 
temp_ax = subplot(3,5,14); set(temp_ax,'xcolor','none'); set(temp_ax,'ycolor','none')
     text(-0.10,0.5,{...     
                    'Return To Society Model 2';
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor2(1))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove2(1)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor2(2))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove2(1:2)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor2(3))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove2(1:3)),'mm/dd/yy')];...
                    ['  Total Infected:                       ',num2str(round(max(All_recoveries1+y2(:,10))))];...
                    ['  Total Dead:                            ',num2str(round(max(y2(:,10))))];...
                    ['  % Increase in Deaths:           ',num2str(round(100*((max(y2(:,10))) - (max(y1(:,10))))/((max(y2(:,10)))),1)),'%'];...
                    '';
                    'Return To Society Model 3';
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor3(1))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove3(1)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor3(2))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove3(1:2)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor3(3))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove3(1:3)),'mm/dd/yy')];...
                    ['  Total Infected:                       ',num2str(round(max(All_recoveries1+y3(:,10))))];...
                    ['  Total Dead:                            ',num2str(round(max(y3(:,10))))];...
                    ['  % Increase in Deaths:           ',num2str(round(100*((max(y3(:,10))) - (max(y1(:,10))))/((max(y3(:,10)))),1)),'%'];...
                     },'Units','Normalized','FontSize', 12);     
                
temp_ax = subplot(3,5,15);set(temp_ax,'Color','none','Visible','off');
    legend(temp_ax,[rd1 rd2 rd3 rd4],...
        {'Return to Society Model 1',...
        'Return to Society Model 2',...
        'Return to Society Model 3',...
        'Real Data'},'Location','best','FontSize',12)
    
disp(['Death rate is: ', num2str(DR)]);
disp(['Peak R0 is: ', num2str(R_naught_peak)]);

% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'# of people', 'FontSize', 16, 'FontWeight','bold');
sgtitle({'BUHoffman-COVID-Model, New York:';'Medium-term Effects of Multiple "Return of Society" Models of Social Distancing Relaxation';'90% Immunity (Sustained 365 days)'},'FontSize', 16, 'FontWeight','bold')
set(gcf,'color','w');

%% Figure S3 G H

y1 = yout{6,6};
y2 = yout{11,6};
y3 = yout{16,6};

zeta_factor1 = zeta_factor{:,6};
zeta_factor2 = zeta_factor{:,11};
zeta_factor3 = zeta_factor{:,16};

SD_remove1 = SD_remove{:,6};
SD_remove2 = SD_remove{:,6};
SD_remove3 = SD_remove{:,6};


[lambda_vector,~] = lambdakappa(t,Coef);
Total_infxns1 = (cumtrapz(y1(:,3))*nu) + (cumtrapz(y1(:,3))*theta) + (cumtrapz(y1(:,3))*delta);
Total_infxns2 = (cumtrapz(y2(:,3))*nu) + (cumtrapz(y2(:,3))*theta) + (cumtrapz(y2(:,3))*delta);
Total_infxns3 = (cumtrapz(y3(:,3))*nu) + (cumtrapz(y3(:,3))*theta) + (cumtrapz(y3(:,3))*delta);

Total_confirmed1 = (cumtrapz(y1(:,3))*theta) + (cumtrapz(y1(:,3))*delta);
Total_confirmed2 = (cumtrapz(y2(:,3))*theta) + (cumtrapz(y2(:,3))*delta);
Total_confirmed3 = (cumtrapz(y3(:,3))*theta) + (cumtrapz(y3(:,3))*delta);

% Summary Outputs
All_recoveries1 = (cumtrapz(y1(:,4))*phi) + (cumtrapz(y1(:,5))*omega) + (cumtrapz(y1(:,6)).*lambda_vector);

DR1 = y(end,10)/(All_recoveries(end)+y1(end,10));

date_times = Time(1):(Time(1)+t(end)-1);
fig=figure;
subplot(3,5,1); plot(date_times(1:end),Total_infxns1,'LineWidth',2);
    hold on; title({'# Total Infections';'(Undocumented + confirmed)'}); ylim([0 max(max([Total_infxns1,Total_infxns2,Total_infxns3]))]); xlim([date_times(1) date_times(end)]);
    hold on;plot(date_times(1:end),Total_infxns2,'LineWidth',2)
    plot(date_times(1:end),Total_infxns3,'LineWidth',2)
    
subplot(3,5,2); plot(date_times,(cumtrapz(y1(:,3))*nu),'LineWidth',2); title({'# Total';'Undocumented Infections'}); ylim([0 max(Total_infxns)]); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,(cumtrapz(y2(:,3))*nu),'LineWidth',2);
    plot(date_times,(cumtrapz(y3(:,3))*nu),'LineWidth',2);

subplot(3,5,3); plot(date_times,Total_confirmed1,'LineWidth',2); title({'# Total';'Confirmed Infections'}); title({'# Total';'Confirmed Infections'}); ylim([0 max(Total_infxns)]); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,Total_confirmed2,'LineWidth',2);
    plot(date_times,Total_confirmed3,'LineWidth',2);
    scatter(date_times(1:length(Confirmed)),Confirmed,'LineWidth',1,'MarkerEdgeColor','#D95319')

subplot(3,5,4); plot(date_times,sum(y1(:,5:6),2),'LineWidth',2); title({'# Active';'Confirmed Infections'}); title({'# Active';'Confirmed Infections'}); xlim([date_times(1) date_times(end)]);
    hold on;plot(date_times,sum(y2(:,5:6),2),'LineWidth',2);
    plot(date_times,sum(y3(:,5:6),2),'LineWidth',2);
    scatter(date_times(1:length(Confirmed)),Hospitalized_active+Not_Hospitalized-y(1:length(Confirmed),8),'LineWidth',1, 'MarkerEdgeColor','#D95319')
    
subplot(3,5,5); plot(date_times(2:end),diff(Total_confirmed1),'LineWidth',2); title({'# New';'Confirmed Infections'}); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times(2:end),diff(Total_confirmed2),'LineWidth',2);
    plot(date_times(2:end),diff(Total_confirmed3),'LineWidth',2);
    scatter(date_times(2:length(Not_Hospitalized)),diff(Confirmed),'LineWidth',1, 'MarkerEdgeColor','#D95319')
    
subplot(3,5,6); plot(date_times,(cumtrapz(y1(:,3))*delta),'LineWidth',2); title({'# Total Confirmed';'Hospitalized'}); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,(cumtrapz(y2(:,3))*delta),'LineWidth',2);
    plot(date_times,(cumtrapz(y3(:,3))*delta),'LineWidth',2);
    scatter(date_times(1:length(Hospitalized_tot)),Hospitalized_active+Deaths+Discharged,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,7); plot(date_times,y1(:,6),'LineWidth',2); title('# Active Hospitalized'); xlim([date_times(1) date_times(end)]);
    hold on;  plot(date_times,y2(:,6),'LineWidth',2);
    plot(date_times,y3(:,6),'LineWidth',2);
    scatter(date_times(1:length(Hospitalized_active)),Hospitalized_active,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,8); plot(date_times,(cumtrapz(y1(:,6)).*lambda_vector),'LineWidth',2); title({'# Total';'Hospital Discharges'}); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,(cumtrapz(y2(:,6)).*lambda_vector),'LineWidth',2);
    plot(date_times,(cumtrapz(y3(:,6)).*lambda_vector),'LineWidth',2);
    scatter(date_times(1:length(Discharged)),Discharged,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,9); rd1 = plot(date_times(2:end),diff((cumtrapz(y1(:,6)).*lambda_vector)),'LineWidth',2); title('# New Discharges'); xlim([date_times(1) date_times(end)]);
    hold on; rd2 = plot(date_times(2:end),diff((cumtrapz(y2(:,6)).*lambda_vector)),'LineWidth',2);
    rd3 = plot(date_times(2:end),diff((cumtrapz(y3(:,6)).*lambda_vector)),'LineWidth',2);
    rd4 = scatter(date_times(2:length(Discharged)),diff(Discharged),'LineWidth',1,'MarkerEdgeColor','#D95319');
    
subplot(3,5,10); plot(date_times,y1(:,10),'LineWidth',2); title('# Deaths'); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,y2(:,10),'LineWidth',2);
    plot(date_times,y3(:,10),'LineWidth',2);
    scatter(date_times(1:length(Deaths)),Deaths,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,11); plot(date_times(2:end),(diff(y1(:,10))),'LineWidth',2); title('# New Deaths'); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times(2:end),(diff(y2(:,10))),'LineWidth',2);
    plot(date_times(2:end),(diff(y3(:,10))),'LineWidth',2);
    scatter(date_times(2:length(Deaths)),diff(Deaths),'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,12); plot(date_times(1:end),R_naught'.*y1(:,1)/Npop,'LineWidth',2); title({'Effective';'Reproductive #'}); xlim([date_times(1) date_times(end)]); ylim([0 2]);
    hold on; plot(date_times(1:end),R_naught'.*y2(:,1)/Npop,'LineWidth',2);
    plot(date_times(1:end),R_naught'.*y3(:,1)/Npop,'LineWidth',2);

temp_ax = subplot(3,5,13); set(temp_ax,'xcolor','none'); set(temp_ax,'ycolor','none')
     text(-0.10,0.6,{...
                    'Summary Values';...
                    ['  Peak Reproductive Number:     ',num2str(round(R_naught_peak,2))];...
                    ['  Undocumented Rate:              ',num2str(Asymp_rate*100),'%'];...
                    ['  SD Started:                            ',datestr(date_times(SD_delay+1),'mm/dd/yy')];...
                    ['  Mortality Rate:                       ',num2str(round(DR1*100,2)),'%'];...
                    ' ';...
                    'Return To Society Model 1';
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor1(1))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove1(1)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor1(2))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove1(1:2)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor1(3))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove1(1:3)),'mm/dd/yy')];...
                    ['  Total Infected:                       ',num2str(round(max(All_recoveries1+y1(:,10))))];...
                    ['  Total Dead:                            ',num2str(round(max(y1(:,10))))];...
                    ['  % Increase in Deaths:           ',num2str(round(100*((max(y1(:,10))) - (max(y1(:,10))))/((max(y1(:,10)))),1)),'%'];...
                     },'Units','Normalized','FontSize', 12);
 
temp_ax = subplot(3,5,14); set(temp_ax,'xcolor','none'); set(temp_ax,'ycolor','none')
     text(-0.10,0.5,{...     
                    'Return To Society Model 2';
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor2(1))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove2(1)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor2(2))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove2(1:2)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor2(3))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove2(1:3)),'mm/dd/yy')];...
                    ['  Total Infected:                       ',num2str(round(max(All_recoveries1+y2(:,10))))];...
                    ['  Total Dead:                            ',num2str(round(max(y2(:,10))))];...
                    ['  % Increase in Deaths:           ',num2str(round(100*((max(y2(:,10))) - (max(y1(:,10))))/((max(y2(:,10)))),1)),'%'];...
                    '';
                    'Return To Society Model 3';
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor3(1))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove3(1)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor3(2))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove3(1:2)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor3(3))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove3(1:3)),'mm/dd/yy')];...
                    ['  Total Infected:                       ',num2str(round(max(All_recoveries1+y3(:,10))))];...
                    ['  Total Dead:                            ',num2str(round(max(y3(:,10))))];...
                    ['  % Increase in Deaths:           ',num2str(round(100*((max(y3(:,10))) - (max(y1(:,10))))/((max(y3(:,10)))),1)),'%'];...
                     },'Units','Normalized','FontSize', 12);     
                
temp_ax = subplot(3,5,15);set(temp_ax,'Color','none','Visible','off');
    legend(temp_ax,[rd1 rd2 rd3 rd4],...
        {'Return to Society Model 1',...
        'Return to Society Model 2',...
        'Return to Society Model 3',...
        'Real Data'},'Location','best','FontSize',12)
    
disp(['Death rate is: ', num2str(DR)]);
disp(['Peak R0 is: ', num2str(R_naught_peak)]);

% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'# of people', 'FontSize', 16, 'FontWeight','bold');
sgtitle({'BUHoffman-COVID-Model, New York:';'Medium-term Effects of Multiple "Return of Society" Models of Social Distancing Relaxation';'90% Immunity (Sustained 365 days)'},'FontSize', 16, 'FontWeight','bold')
set(gcf,'color','w');

%% Figure 3
% Effect on Time to steady state social distancing and steady state level
% total deaths by 9/1/2020

% Immunity
Immunity_fraction = 0.70;% fraction of recovered who develop immunity
tau = 1-Immunity_fraction; % fraction of recovered who DO NOT develop immunity
pi = 1/1825; % 1/(Duration of immunity)


SD_remove = cell(1,20);
for i = 1:20
    time_factor = 3*i;
    SDout = ones(1,3);
    SDout(1,1) = 71; % One time relaxation on 6/1/2020
    SDout(1,2) = 180+time_factor; % One time restart SD starting 12/01/2021 (then increasing by 3 days)
    SDout(1,3) = 93-time_factor; % One time relax SD starting 3/01/2021  
    SD_remove{1,i} = SDout;
end

zeta_factor = cell(1,21);
for i = 1:21
    relax_factor = (1/40)*i;
    zetaout = ones(1,3);
    zetaout(1,1) = 1/(0.5); %One time relaxation by 50% on 6/1/2020
    zetaout(1,2) = 1/(0.475+relax_factor); 
    zetaout(1,3) = 1/(0.5); %One time relaxation by 50% on 3/01/2021
    zeta_factor{1,i} = zetaout;
end

% WINDOW FOR SIMULATION
window_sim = [1, 547];% 9/1/21
chi = importflu(Time,window_sim,'factor',flufactor);
yout = cell(21,20);
tic
for i = 1:21
    zeta_factor_i = zeta_factor{1,i};
    for j = 1:20
        [t,y] = ode45(@(t,y) SEIAQHRRRDP_deqs(t, y, Coef, Npop, rho, gamma,...
            nu, phi, delta, theta, omega,tau, pi, SD_delay, ...
            SD_remove{1,j}, zeta_factor_i, window_sim, chi), window_sim(1):window_sim(end), ic);
        yout{i,j} = y;
    end
end
toc

time2SS_SD = zeros(20,1);
for i = 1:20
    time2SS_SD(i,1) = sum(SD_remove{1,i}(:,1:2));
end
SS_SD = zeros(21,1);
for i = 1:21
    SS_SD(i,1) = 100-(100/(zeta_factor{1,22-i}(2)));
end
SS_SD = flipud(SS_SD);

maxdeaths = zeros(21,20);
for i = 1:21
    for j = 1:20
        maxdeaths(i,j) = yout{i,j}(end,10);
    end
end

R_naught = ((chi*rho*Coef(3))/(nu+theta+delta)) + ((chi*rho*Coef(3)*nu)/(phi*(nu+theta+delta)));
R_naught_peak = max(((chi*rho*Coef(3))/(nu+theta+delta)) + ((chi*rho*Coef(3)*nu)/(phi*(nu+theta+delta))));

ReThresh = zeros(21,20);
for i = 1:21
    for j = 1:20
        Reff = (R_naught'.*yout{i,j}(:,1)/Npop);
        above = find( Reff > 1);
        above2 = find( above > (sum(SD_remove{1,j}(:,1:2))+ SD_delay));
        if ~isempty(above2)
            ReThresh(i,j) = 1;
        end
    end
end

figure;
subplot(2,2,1);heatmap(date_times(1)+SD_delay+time2SS_SD,SS_SD,maxdeaths,...
    'GridVisible','off','CellLabelColor','none','Colormap',parula,...
    'XLabel','Date of Social Distancing Resumption',...
    'YLabel','% Reduction in Social Distancing')

subplot(2,2,2);heatmap(date_times(1)+SD_delay+time2SS_SD,SS_SD,ReThresh,...
    'GridVisible','off','CellLabelColor','none','Colormap',[0.2422,0.1504,0.6603;0.9769,0.9839,0.0805],...
    'XLabel','Date of Social Distancing Resumption',...
    'YLabel','% Reduction in Social Distancing')
%

y1 = yout{21,1};
y2 = yout{11,1};
y3 = yout{1,1};

zeta_factor1 = zeta_factor{:,21};
zeta_factor2 = zeta_factor{:,11};
zeta_factor3 = zeta_factor{:,1};

SD_remove1 = SD_remove{:,1};
SD_remove2 = SD_remove{:,1};
SD_remove3 = SD_remove{:,1};


[lambda_vector,~] = lambdakappa(t,Coef);
Total_infxns1 = (cumtrapz(y1(:,3))*nu) + (cumtrapz(y1(:,3))*theta) + (cumtrapz(y1(:,3))*delta);
Total_infxns2 = (cumtrapz(y2(:,3))*nu) + (cumtrapz(y2(:,3))*theta) + (cumtrapz(y2(:,3))*delta);
Total_infxns3 = (cumtrapz(y3(:,3))*nu) + (cumtrapz(y3(:,3))*theta) + (cumtrapz(y3(:,3))*delta);

Total_confirmed1 = (cumtrapz(y1(:,3))*theta) + (cumtrapz(y1(:,3))*delta);
Total_confirmed2 = (cumtrapz(y2(:,3))*theta) + (cumtrapz(y2(:,3))*delta);
Total_confirmed3 = (cumtrapz(y3(:,3))*theta) + (cumtrapz(y3(:,3))*delta);

% Summary Outputs
All_recoveries1 = (cumtrapz(y1(:,4))*phi) + (cumtrapz(y1(:,5))*omega) + (cumtrapz(y1(:,6)).*lambda_vector);

date_times = Time(1):(Time(1)+t(end)-1);
fig=figure;
subplot(3,5,1); plot(date_times(1:end),Total_infxns1,'LineWidth',2);
    hold on; title({'# Total Infections';'(Undocumented + confirmed)'}); ylim([0 1.1*max(max([Total_infxns1,Total_infxns2,Total_infxns3]))]); xlim([date_times(1) date_times(end)]);
    hold on;plot(date_times(1:end),Total_infxns2,'LineWidth',2)
    plot(date_times(1:end),Total_infxns3,'LineWidth',2)
    
subplot(3,5,2); plot(date_times,(cumtrapz(y1(:,3))*nu),'LineWidth',2); title({'# Total';'Undocumented Infections'}); ylim([0 1.1*max(max([Total_infxns1,Total_infxns2,Total_infxns3]))]); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,(cumtrapz(y2(:,3))*nu),'LineWidth',2);
    plot(date_times,(cumtrapz(y3(:,3))*nu),'LineWidth',2);

subplot(3,5,3); plot(date_times,Total_confirmed1,'LineWidth',2); title({'# Total';'Confirmed Infections'}); title({'# Total';'Confirmed Infections'}); ylim([0 1.1*max(max([Total_infxns1,Total_infxns2,Total_infxns3]))]); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,Total_confirmed2,'LineWidth',2);
    plot(date_times,Total_confirmed3,'LineWidth',2);
    scatter(date_times(1:length(Confirmed)),Confirmed,'LineWidth',1,'MarkerEdgeColor','#D95319')

subplot(3,5,4); plot(date_times,sum(y1(:,5:6),2),'LineWidth',2); title({'# Active';'Confirmed Infections'}); title({'# Active';'Confirmed Infections'}); xlim([date_times(1) date_times(end)]);
    hold on;plot(date_times,sum(y2(:,5:6),2),'LineWidth',2);
    plot(date_times,sum(y3(:,5:6),2),'LineWidth',2);
    scatter(date_times(1:length(Confirmed)),Hospitalized_active+Not_Hospitalized-y(1:length(Confirmed),8),'LineWidth',1, 'MarkerEdgeColor','#D95319')
    
subplot(3,5,5); plot(date_times(2:end),diff(Total_confirmed1),'LineWidth',2); title({'# New';'Confirmed Infections'}); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times(2:end),diff(Total_confirmed2),'LineWidth',2);
    plot(date_times(2:end),diff(Total_confirmed3),'LineWidth',2);
    scatter(date_times(2:length(Not_Hospitalized)),diff(Confirmed),'LineWidth',1, 'MarkerEdgeColor','#D95319')
    
subplot(3,5,6); plot(date_times,(cumtrapz(y1(:,3))*delta),'LineWidth',2); title({'# Total Confirmed';'Hospitalized'}); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,(cumtrapz(y2(:,3))*delta),'LineWidth',2);
    plot(date_times,(cumtrapz(y3(:,3))*delta),'LineWidth',2);
    scatter(date_times(1:length(Hospitalized_tot)),Hospitalized_active+Deaths+Discharged,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,7); plot(date_times,y1(:,6),'LineWidth',2); title('# Active Hospitalized'); xlim([date_times(1) date_times(end)]);
    hold on;  plot(date_times,y2(:,6),'LineWidth',2);
    plot(date_times,y3(:,6),'LineWidth',2);
    scatter(date_times(1:length(Hospitalized_active)),Hospitalized_active,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,8); plot(date_times,(cumtrapz(y1(:,6)).*lambda_vector),'LineWidth',2); title({'# Total';'Hospital Discharges'}); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,(cumtrapz(y2(:,6)).*lambda_vector),'LineWidth',2);
    plot(date_times,(cumtrapz(y3(:,6)).*lambda_vector),'LineWidth',2);
    scatter(date_times(1:length(Discharged)),Discharged,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,9); rd1 = plot(date_times(2:end),diff((cumtrapz(y1(:,6)).*lambda_vector)),'LineWidth',2); title('# New Discharges'); xlim([date_times(1) date_times(end)]);
    hold on; rd2 = plot(date_times(2:end),diff((cumtrapz(y2(:,6)).*lambda_vector)),'LineWidth',2);
    rd3 = plot(date_times(2:end),diff((cumtrapz(y3(:,6)).*lambda_vector)),'LineWidth',2);
    rd4 = scatter(date_times(2:length(Discharged)),diff(Discharged),'LineWidth',1,'MarkerEdgeColor','#D95319');
    
subplot(3,5,10); plot(date_times,y1(:,10),'LineWidth',2); title('# Deaths'); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,y2(:,10),'LineWidth',2);
    plot(date_times,y3(:,10),'LineWidth',2);
    scatter(date_times(1:length(Deaths)),Deaths,'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,11); plot(date_times(2:end),(diff(y1(:,10))),'LineWidth',2); title('# New Deaths'); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times(2:end),(diff(y2(:,10))),'LineWidth',2);
    plot(date_times(2:end),(diff(y3(:,10))),'LineWidth',2);
    scatter(date_times(2:length(Deaths)),diff(Deaths),'LineWidth',1,'MarkerEdgeColor','#D95319')
    
subplot(3,5,12); plot(date_times(1:end),R_naught'.*y1(:,1)/Npop,'LineWidth',2); title({'Effective';'Reproductive #'}); xlim([date_times(1) date_times(end)]); ylim([0 2]);
    hold on; plot(date_times(1:end),R_naught'.*y2(:,1)/Npop,'LineWidth',2);
    plot(date_times(1:end),R_naught'.*y3(:,1)/Npop,'LineWidth',2);

temp_ax = subplot(3,5,13); set(temp_ax,'xcolor','none'); set(temp_ax,'ycolor','none')
     text(-0.10,0.6,{...
                    'Summary Values';...
                    ['  Peak Reproductive Number:     ',num2str(round(R_naught_peak,2))];...
                    ['  Undocumented Rate:              ',num2str(Asymp_rate*100),'%'];...
                    ['  SD Started:                            ',datestr(date_times(SD_delay+1),'mm/dd/yy')];...
                    ['  Mortality Rate:                       ',num2str(round(DR*100,2)),'%'];...
                    ' ';...
                    'Return To Society Model 1';
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor1(1))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove1(1)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor1(2))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove1(1:2)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor1(3))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove1(1:3)),'mm/dd/yy')];...
                    ['  Total Infected:                       ',num2str(round(max(All_recoveries1+y1(:,10))))];...
                    ['  Total Dead:                            ',num2str(round(max(y1(:,10))))];...
                    ['  % Increase in Deaths:           ',num2str(round(100*((max(y1(:,10))) - (max(y1(:,10))))/((max(y1(:,10)))),1)),'%'];...
                     },'Units','Normalized','FontSize', 12);
 
temp_ax = subplot(3,5,14); set(temp_ax,'xcolor','none'); set(temp_ax,'ycolor','none')
     text(-0.10,0.5,{...     
                    'Return To Society Model 2';
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor2(1))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove2(1)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor2(2))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove2(1:2)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor2(3))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove2(1:3)),'mm/dd/yy')];...
                    ['  Total Infected:                       ',num2str(round(max(All_recoveries1+y2(:,10))))];...
                    ['  Total Dead:                            ',num2str(round(max(y2(:,10))))];...
                    ['  % Increase in Deaths:           ',num2str(round(100*((max(y2(:,10))) - (max(y1(:,10))))/((max(y2(:,10)))),1)),'%'];...
                    '';
                    'Return To Society Model 3';
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor3(1))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove3(1)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor3(2))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove3(1:2)),'mm/dd/yy')];...
                    ['  SD Relaxed by ',num2str(round(100-100/zeta_factor3(3))),'% ','on ',datestr(date_times(1)+SD_delay+sum(SD_remove3(1:3)),'mm/dd/yy')];...
                    ['  Total Infected:                       ',num2str(round(max(All_recoveries1+y3(:,10))))];...
                    ['  Total Dead:                            ',num2str(round(max(y3(:,10))))];...
                    ['  % Increase in Deaths:           ',num2str(round(100*((max(y3(:,10))) - (max(y1(:,10))))/((max(y3(:,10)))),1)),'%'];...
                     },'Units','Normalized','FontSize', 12);     
                
temp_ax = subplot(3,5,15);set(temp_ax,'Color','none','Visible','off');
    legend(temp_ax,[rd1 rd2 rd3 rd4],...
        {'Return to Society Model 1',...
        'Return to Society Model 2',...
        'Return to Society Model 3',...
        'Real Data'},'Location','best','FontSize',12)
    
disp(['Death rate is: ', num2str(DR)]);
disp(['Peak R0 is: ', num2str(R_naught_peak)]);

% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'# of people', 'FontSize', 16, 'FontWeight','bold');
sgtitle({'BUHoffman-COVID-Model, New York:';'Medium-term Effects of Multiple "Return of Society" Models of Social Distancing Relaxation';'90% Immunity (Sustained 365 days)'},'FontSize', 16, 'FontWeight','bold')
set(gcf,'color','w');

%% Figure 4 part 1
% Effect on Time to steady state social distancing and steady state level
% total deaths by 9/1/2020

% % Immunity
% Immunity_fraction = 0.70;% fraction of recovered who develop immunity
% tau = 1-Immunity_fraction; % fraction of recovered who DO NOT develop immunity
% pi = 1/1825; % 1/(Duration of immunity)

pi_out = zeros(1,21);
for i = 1:21
    pi_out(1,i) = 1/(0.5 + (i*365));
end

tau_out = zeros(1,21);
for i = 1:21
    tau_factor = (i/30);
    tau_out(1,i) = 0.3 + tau_factor;
end
tau_out(1,21)=0.9999;
tau_out=1-tau_out;

SD_remove = ones(1,3);
SD_remove(1,1) = 71; % One time relaxation on 6/1/2020
SD_remove(1,2) = 183; % One time restart SD starting 11/01/2020
SD_remove(1,3) = 93; % One time relax SD starting 3/01/2021

zeta_factor = ones(1,3);
zeta_factor(1,1) = 1/(0.5); %One time relaxation by 50% on 6/1/2020
zeta_factor(1,2) = 1/(0.75); 
zeta_factor(1,3) = 1/(0.5); %One time relaxation by 50% on 3/01/2021

% WINDOW FOR SIMULATION
% window_sim = [1, 3833];% 9/1/2030
window_sim = [1, 2007];% 9/1/2025
chi = importflu(Time,window_sim,'factor',flufactor);
yout = cell(21,21);
tic
for i = 1:21
    tau_i =  tau_out(i);
    parfor j = 1:21
        [t,y] = ode45(@(t,y) SEIAQHRRRDP_deqs(t, y, Coef, Npop, rho, gamma,...
            nu, phi, delta, theta, omega, tau_i, pi_out(j), SD_delay, ...
            SD_remove, zeta_factor, window_sim, chi), window_sim(1):window_sim(end), ic);
        yout{i,j} = y;
    end
end
toc

sound_end=load('handel.mat');
sound(sound_end.y,sound_end.Fs);

maxdeaths = zeros(21,21);
for i = 1:21
    for j = 1:21
        maxdeaths(i,j) = yout{i,j}(end,10);
    end
end
% maxdeaths = flipud(maxdeaths);

R_naught = ((chi*rho*Coef(3))/(nu+theta+delta)) + ((chi*rho*Coef(3)*nu)/(phi*(nu+theta+delta)));
R_naught_peak = max(((chi*rho*Coef(3))/(nu+theta+delta)) + ((chi*rho*Coef(3)*nu)/(phi*(nu+theta+delta))));

endemic_thresh = zeros(21,21);
for i = 1:21
    for j = 1:21
        total_infections_endemic = (cumtrapz(yout{i,j}(:,3))*nu) + (cumtrapz(yout{i,j}(:,3))*theta) + (cumtrapz(yout{i,j}(:,3))*delta);
        endemic_avg = (total_infections_endemic(end)-total_infections_endemic(end-730))/2;
        if endemic_avg > 20000
            endemic_thresh(i,j) = 1;
        end
    end
end
%
figure;subplot(2,2,1);
heatmap(round(1./pi_out/365),fliplr(round(100-tau_out.*100,1)),flipud(log10(maxdeaths)),...
    'GridVisible','off','CellLabelColor','none','Colormap',parula,...
    'XLabel','Average Duration of Immunity (years)',...
    'ColorLimits',[log10(1e5) log10(1e6)],...
    'YLabel','% Infected Who Develop Immunity')
axs = struct(gca);
cb = axs.Colorbar;
cb.Ticks = [log10(1e5) log10(2e5) log10(3e5) log10(4e5) log10(5e5) log10(6e5)...
    log10(7e5) log10(8e5) log10(9e5) log10(1e6)]; 
cb.TickLabels = {num2str(10) num2str(20) num2str(30) num2str(40) num2str(50) num2str(60) num2str(70)...
    num2str(80) num2str(90) num2str(100)};

subplot(2,2,2);heatmap(round(1./pi_out/365),round(fliplr(100-tau_out.*100),1),flipud(endemic_thresh),...
    'GridVisible','off','CellLabelColor','none','Colormap',[0.2422,0.1504,0.6603;0.9769,0.9839,0.0805],...
    'XLabel','Average Duration of Immunity (years)',...
    'YLabel','% Infected Who Develop Immunity')

y1 = yout{21,21};
y2 = yout{12,5};
y3 = yout{6,2};

tau_out1 = tau_out(1,21);
tau_out2 = tau_out(1,12);
tau_out3 = tau_out(1,6);

pi_out1 = pi_out(1,20);
pi_out2 = pi_out(1,5);
pi_out3 = pi_out(1,2);

zeta_factor1 = zeta_factor;
zeta_factor2 = zeta_factor;
zeta_factor3 = zeta_factor;

SD_remove1 = SD_remove;
SD_remove2 = SD_remove;
SD_remove3 = SD_remove;

[lambda_vector,kappa_vector] = lambdakappa(1:window_sim(2),Coef);
Total_infxns1 = (cumtrapz(y1(:,3))*nu) + (cumtrapz(y1(:,3))*theta) + (cumtrapz(y1(:,3))*delta);
Total_infxns2 = (cumtrapz(y2(:,3))*nu) + (cumtrapz(y2(:,3))*theta) + (cumtrapz(y2(:,3))*delta);
Total_infxns3 = (cumtrapz(y3(:,3))*nu) + (cumtrapz(y3(:,3))*theta) + (cumtrapz(y3(:,3))*delta);

Total_confirmed1 = (cumtrapz(y1(:,3))*theta) + (cumtrapz(y1(:,3))*delta);
Total_confirmed2 = (cumtrapz(y2(:,3))*theta) + (cumtrapz(y2(:,3))*delta);
Total_confirmed3 = (cumtrapz(y3(:,3))*theta) + (cumtrapz(y3(:,3))*delta);

Recoved_Confirmed_Non_Hosp1 = (cumtrapz(y1(:,5)*omega));
Recoved_Confirmed_Non_Hosp2 = (cumtrapz(y2(:,5)*omega));
Recoved_Confirmed_Non_Hosp3 = (cumtrapz(y3(:,5)*omega));

% Summary Outputs
All_recoveries1 = (cumtrapz(y1(:,4))*phi) + (cumtrapz(y1(:,5))*omega) + (cumtrapz(y1(:,6)).*lambda_vector);
All_recoveries2 = (cumtrapz(y2(:,4))*phi) + (cumtrapz(y2(:,5))*omega) + (cumtrapz(y2(:,6)).*lambda_vector);
All_recoveries3 = (cumtrapz(y3(:,4))*phi) + (cumtrapz(y3(:,5))*omega) + (cumtrapz(y3(:,6)).*lambda_vector);

DR1 = y(end,10)/(All_recoveries(end)+y1(end,10));
DR2 = y(end,10)/(All_recoveries(end)+y2(end,10));
DR3 = y(end,10)/(All_recoveries(end)+y3(end,10));

date_times = Time(1):(Time(1)+window_sim(2)-1);
fig=figure;
subplot(3,5,1); plot(date_times(1:end),Total_infxns1,'LineWidth',2);
    hold on; title({'# Total Infections';'(Undocumented + confirmed)'}); ylim([0 1.1*max(max([Total_infxns1,Total_infxns2,Total_infxns3]))]); xlim([date_times(1) date_times(end)]);
    hold on;plot(date_times(1:end),Total_infxns2,'LineWidth',2)
    plot(date_times(1:end),Total_infxns3,'LineWidth',2)
    
subplot(3,5,2); plot(date_times,(cumtrapz(y1(:,3))*nu),'LineWidth',2); title({'# Total';'Undocumented Infections'}); ylim([0 1.1*max(max([Total_infxns1,Total_infxns2,Total_infxns3]))]); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,(cumtrapz(y2(:,3))*nu),'LineWidth',2);
    plot(date_times,(cumtrapz(y3(:,3))*nu),'LineWidth',2);

subplot(3,5,3); plot(date_times,Total_confirmed1,'LineWidth',2); title({'# Total';'Confirmed Infections'}); title({'# Total';'Confirmed Infections'}); ylim([0 1.1*max(max([Total_infxns1,Total_infxns2,Total_infxns3]))]); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,Total_confirmed2,'LineWidth',2);
    plot(date_times,Total_confirmed3,'LineWidth',2);

subplot(3,5,4); plot(date_times,sum(y1(:,5:6),2),'LineWidth',2); title({'# Active';'Confirmed Infections'}); title({'# Active';'Confirmed Infections'}); xlim([date_times(1) date_times(end)]);
    hold on;plot(date_times,sum(y2(:,5:6),2),'LineWidth',2);
    plot(date_times,sum(y3(:,5:6),2),'LineWidth',2);
    
subplot(3,5,5); plot(date_times(2:end),diff(Total_confirmed1),'LineWidth',2); title({'# New';'Confirmed Infections'}); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times(2:end),diff(Total_confirmed2),'LineWidth',2);
    plot(date_times(2:end),diff(Total_confirmed3),'LineWidth',2);
    
subplot(3,5,6); plot(date_times,(cumtrapz(y1(:,3))*delta),'LineWidth',2); title({'# Total Confirmed';'Hospitalized'}); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,(cumtrapz(y2(:,3))*delta),'LineWidth',2);
    plot(date_times,(cumtrapz(y3(:,3))*delta),'LineWidth',2);
    
subplot(3,5,7); plot(date_times,y1(:,6),'LineWidth',2); title('# Active Hospitalized'); xlim([date_times(1) date_times(end)]);
    hold on;  plot(date_times,y2(:,6),'LineWidth',2);
    plot(date_times,y3(:,6),'LineWidth',2);
     
subplot(3,5,8); plot(date_times,(cumtrapz(y1(:,6)).*lambda_vector'),'LineWidth',2); title({'# Total';'Hospital Discharges'}); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,(cumtrapz(y2(:,6)).*lambda_vector'),'LineWidth',2);
    plot(date_times,(cumtrapz(y3(:,6)).*lambda_vector'),'LineWidth',2);
    
subplot(3,5,9); rd1 = plot(date_times(2:end),diff((cumtrapz(y1(:,6)).*lambda_vector')),'LineWidth',2); title('# New Discharges'); xlim([date_times(1) date_times(end)]);
    hold on; rd2 = plot(date_times(2:end),diff((cumtrapz(y2(:,6)).*lambda_vector')),'LineWidth',2);
    rd3 = plot(date_times(2:end),diff((cumtrapz(y3(:,6)).*lambda_vector')),'LineWidth',2);
     
subplot(3,5,10); plot(date_times,y1(:,10),'LineWidth',2); title('# Deaths'); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times,y2(:,10),'LineWidth',2);
    plot(date_times,y3(:,10),'LineWidth',2);
   
subplot(3,5,11); plot(date_times(2:end),(diff(y1(:,10))),'LineWidth',2); title('# New Deaths'); xlim([date_times(1) date_times(end)]);
    hold on; plot(date_times(2:end),(diff(y2(:,10))),'LineWidth',2);
    plot(date_times(2:end),(diff(y3(:,10))),'LineWidth',2);
    
subplot(3,5,12); plot(date_times(1:end),R_naught'.*y1(:,1)/Npop,'LineWidth',2); title({'Effective';'Reproductive #'}); xlim([date_times(1) date_times(end)]); ylim([0 2]);
    hold on; plot(date_times(1:end),R_naught'.*y2(:,1)/Npop,'LineWidth',2);
    plot(date_times(1:end),R_naught'.*y3(:,1)/Npop,'LineWidth',2);

temp_ax = subplot(3,5,13); set(temp_ax,'xcolor','none'); set(temp_ax,'ycolor','none')
     text(-0.10,0.6,{...
                    'Summary Values';...
                    ['  Peak Reproductive Number:     ',num2str(round(R_naught_peak,2))];...
                    ['  Undocumented Rate:                  ',num2str(Asymp_rate*100),'%'];...
                    ['  SD Started:                                ',datestr(date_times(SD_delay+1),'mm/dd/yy')];...
                    ['  Mortality Rate:                           ',num2str(round(DR*100,2)),'%'];...
                    ' ';...  
                    'Endemic Model 1';
                    ['  % develop immunity:               ',num2str((1-tau_out1)*100)];...
                    ['  Average immunity duration:    ',num2str(round(1/pi_out1/365,1)),' years'];...
                    ['  Total Infected:                         ',num2str(round(max(max(All_recoveries1+y1(:,10)))))];...
                    ['  Total Dead:                              ',num2str(round(max(max(y1(:,10)))))];...
                    ['  % Increase in Deaths:             ',num2str(round(100*((max(max(y1(:,10)))) - (max(max(y1(:,10)))))/((max(max((y1(:,10)))))),1)),'%'];...
                     },'Units','Normalized','FontSize', 12);
 
temp_ax = subplot(3,5,14); set(temp_ax,'xcolor','none'); set(temp_ax,'ycolor','none')
     text(-0.10,0.5,{...     
                    'Endemic Model 2';
                    ['  % develop immunity:               ',num2str((1-tau_out2)*100)];...
                    ['  Average immunity duration:    ',num2str(round(1/pi_out2/365,1)),' years'];...
                    ['  Total Infected:                         ',num2str(round(max(max(All_recoveries1+y2(:,10)))))];...
                    ['  Total Dead:                              ',num2str(round(max(max(y2(:,10)))))];...
                    ['  % Increase in Deaths:             ',num2str(round(100*((max(max(y2(:,10)))) - (max(max(y1(:,10)))))/((max(max((y1(:,10)))))),1)),'%'];...
                    '';
                    'Endemic Model 3';
                    ['  % develop immunity:               ',num2str((1-tau_out3)*100)];...
                    ['  Average immunity duration:    ',num2str(round(1/pi_out3/365,1)),' years'];...
                    ['  Total Infected:                         ',num2str(round(max(max(All_recoveries1+y1(:,10)))))];...
                    ['  Total Dead:                              ',num2str(round(max(max(y3(:,10)))))];...
                    ['  % Increase in Deaths:             ',num2str(round(100*((max(max(y3(:,10)))) - (max(max(y1(:,10)))))/((max(max((y1(:,10)))))),1)),'%'];...
                     },'Units','Normalized','FontSize', 12); 
            
temp_ax = subplot(3,5,15);set(temp_ax,'Color','none','Visible','off');
    legend(temp_ax,[rd1 rd2 rd3],...
        {'Endemic Model 1',...
        'Endemic Model 2',...
        'Endemic Model 3'},'Location','best','FontSize',12)
    
disp(['Death rate is: ', num2str(DR)]);
disp(['Peak R0 is: ', num2str(R_naught_peak)]);

% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'# of people', 'FontSize', 16, 'FontWeight','bold');
sgtitle({'BUHoffman-COVID-Model, New York:';'Effect of Immunity on the Endemic State of SARS-CoV-2'},'FontSize', 16, 'FontWeight','bold')
set(gcf,'color','w');










