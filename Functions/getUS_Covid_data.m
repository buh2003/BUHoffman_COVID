function [Confirmed,Hospitalized_tot,Hospitalized_active,Discharged,Deaths,Time] = getUS_Covid_data(state)
% https://github.com/COVID19Tracking/covid-tracking-data/blob/master/data/states_daily_4pm_et.csv
%% Import the data
address = 'https://raw.githubusercontent.com/COVID19Tracking/covid-tracking-data/master/data/states_daily_4pm_et.csv';
websave('dummy.csv',address);
alldata = readtable('dummy.csv');
delete('dummy.csv')
%% Extract desired state
statedata = alldata(contains(alldata.state,state),:);
Confirmed = flipud(statedata.positive);
Hospitalized_tot = flipud(statedata.hospitalized);
Hospitalized_active = flipud(statedata.hospitalizedCurrently);
Discharged = flipud(statedata.recovered);
Deaths = flipud(statedata.death);
Time = datetime(flipud(statedata.date),'ConvertFrom','yyyymmdd');

disp(['Last updated ' datestr(Time(end))]);
end

