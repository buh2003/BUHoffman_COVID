function [output] = importflu(Time,window,varargin)
%%
p = inputParser;
addOptional(p,'factor',1);
parse(p,varargin{:});
factor = p.Results.factor;

fludata = readtable('/Users/benjaminhoffman/Desktop/COVID modeling/SEIAQHRRRDP/Flu_Corona_data_for_constraining/FluViewPhase2Data/ILINet.csv','PreserveVariableNames',true);
fludata = table(fludata.WEEK, fludata{:,5},'VariableNames',{'WEEK','WEIGHTEDILI'});

fludata_mean=varfun(@mean,fludata,'InputVariables','WEIGHTEDILI',...
       'GroupingVariables','WEEK');
fludata_mean=fludata_mean(1:52,:);
fludata_mean(53,:)=table(53,9,mean([fludata_mean.mean_WEIGHTEDILI(52),fludata_mean.mean_WEIGHTEDILI(1)]));

%% Import non-COVID-19 HCoV data 
% 2014-2017 From https://www.sciencedirect.com/science/article/pii/S1386653218300325#bib0035
% 2018-2020 From https://www.cdc.gov/surveillance/nrevss/coronavirus/index.html

loaded = load('HCoV_agg_table_out.mat');
HCoV_agg_table_out = loaded.HCoV_agg_table_out;

%% Combine w/ flu data 

combined_data = table(HCoV_agg_table_out.time,fludata_mean.mean_WEIGHTEDILI.*HCoV_agg_table_out.mean,...
    'VariableNames',{'time','factor'});
    
%% Create vector as long as input window
years_time = unique(year(Time(1):Time(1)+window(2)-1));
t_out = [];
data_out = [];
for i = years_time(1):years_time(end)
    date_i = ['0101',num2str(i)];
    weeks_time = datetime(date_i,'Inputformat','MMddyyyy');
    weeks_time = (weeks_time(1):7:365+weeks_time(1))';
    t=weeks_time(1):1:weeks_time(1)+365;
    data_interp = interp1(weeks_time,combined_data.factor,t);
    t_out = [t_out,t];
    data_out = [data_out, data_interp];
end

nandata = isnan(data_out);
data_out2 = interp1(t_out(~nandata),data_out(~nandata),t_out);

%normalize to 30 days after pandemic started
norm_date = find(t_out == Time(1)+30);
fludata_norm = data_out2/data_out2(norm_date);
out_start = find(t_out == Time(1));
output = fludata_norm(out_start:(out_start + window(2) -1));

output = ((output-1)/factor)+1;

