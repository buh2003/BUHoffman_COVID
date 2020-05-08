%% ################# SEIR_QHRRD_Model ######################
set(0,'DefaultFigureWindowStyle','docked')

% S: number of succeptibles
% E: number of exposed
% I: number of infectious (not quarantined)
% R: number of recovered
% Q: number of quarantined, active cases not requiring hospitalization
% H: number of hospitalized, active cases requiring hospitalization
% RQ: number of recovered, cases not requiring hospitalization
% RH: number of recovered,  cases requiring hospitalization
% D: number of dead

%% Import Data
state = 'NY';
[Confirmed,Hospitalized_tot,Discharged,Deaths,Time] = getUS_Covid_data(state);
figure;plot(Time, Confirmed, Time, Hospitalized_tot, Time, Discharged, Time, Deaths); legend('Confirmed','Hospitalized Total','Discharged','Deaths');

%% Estimate active non-hospitalized cases
% Active = Q+H = Confirmed - RQ - RH - Deaths
% Q+RH = Confirmed - Hospitalized_tot