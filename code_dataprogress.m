clc
clear

cmip6_CTQ_waterY10(1) = struct();
load("I:\1-code\20250706-Qproject\code_analysis\data_local_water_year_startMonth.mat","start_month2320");

load("I:\1-code\20250706-Qproject\data_2320sites\6-project-runoff\2-hebin-project runoff\" + ...
    "data_projectRunoff_cmip6_his.mat","Qday");
data = mean(Qday,3,'omitnan');
times = (datetime(1980,1,1):datetime(2014,12,31))';
for ss = 1:2320
    disp([0 ss])
    tic
    runoff = data(366:end,ss);
    [result1(:,ss),CTP(:,ss)] = Fun_compute_CTQ_CTP_waterYear(times, runoff,start_month2320(ss));
    toc
end
cmip6_CTQ_waterY10(1).CTQ_his = CTP;
clearvars result1 CTP cmip6_his data


load("I:\1-code\20250706-Qproject\data_2320sites\6-project-runoff\2-hebin-project runoff\" + ...
    "data_projectRunoff_cmip6_ssp126.mat","Qday");
data = mean(Qday,3,'omitnan');
times = (datetime(2015,1,1):datetime(2100,12,31))';
for ss = 1:2320
    disp([126 ss])
    tic
    runoff = data(366:end,ss);
    [result1(:,ss),CTP(:,ss)] = Fun_compute_CTQ_CTP_waterYear(times, runoff,start_month2320(ss));
    toc
end
cmip6_CTQ_waterY10(1).CTQ_ssp126 = CTP;
clearvars result1 CTP cmip6_his data

load("I:\1-code\20250706-Qproject\data_2320sites\6-project-runoff\2-hebin-project runoff\" + ...
    "data_projectRunoff_cmip6_ssp245.mat","Qday");
data = mean(Qday,3,'omitnan');
times = (datetime(2015,1,1):datetime(2100,12,31))';
for ss = 1:2320
    disp([245 ss])
    tic
    runoff = data(366:end,ss);
    [result1(:,ss),CTP(:,ss)] = Fun_compute_CTQ_CTP_waterYear(times, runoff,start_month2320(ss));
    toc
end
cmip6_CTQ_waterY10(1).CTQ_ssp245 = CTP;
clearvars result1 CTP cmip6_his data

load("I:\1-code\20250706-Qproject\data_2320sites\6-project-runoff\2-hebin-project runoff\" + ...
    "data_projectRunoff_cmip6_ssp370.mat","Qday");
data = mean(Qday,3,'omitnan');
times = (datetime(2015,1,1):datetime(2100,12,31))';
for ss = 1:2320
    disp([370 ss])
    tic
    runoff = data(366:end,ss);
    [result1(:,ss),CTP(:,ss)] = Fun_compute_CTQ_CTP_waterYear(times, runoff,start_month2320(ss));
    toc
end
cmip6_CTQ_waterY10(1).CTQ_ssp370 = CTP;
clearvars result1 CTP cmip6_his data

load("I:\1-code\20250706-Qproject\data_2320sites\6-project-runoff\2-hebin-project runoff\" + ...
    "data_projectRunoff_cmip6_ssp585.mat","Qday");
data = mean(Qday,3,'omitnan');
times = (datetime(2015,1,1):datetime(2100,12,31))';
for ss = 1:2320
    disp([585 ss])
    tic
    runoff = data(366:end,ss);
    [result1(:,ss),CTP(:,ss)] = Fun_compute_CTQ_CTP_waterYear(times, runoff,start_month2320(ss));
    toc
end
cmip6_CTQ_waterY10(1).CTQ_ssp585 = CTP;
clearvars result1 CTP cmip6_his data

function [centroid_years, centroid_times] = Fun_compute_CTQ_CTP_waterYear(datetimes, runoff, start_month)
if nargin < 3, start_month = 1; end
validateattributes(datetimes, {'datetime'}, {'vector'}, mfilename, 't')
validateattributes(runoff, {'numeric'}, {'vector', 'numel', numel(datetimes)}, mfilename, 'Q')
validateattributes(start_month, {'numeric'}, {'integer', 'scalar', '>=',1,'<=',12}, mfilename, 'start_month')
datetimes = datetimes(:);
runoff = runoff(:);
[~, sort_idx] = sort(datetimes);
datetimes = datetimes(sort_idx);
runoff = runoff(sort_idx);
hydro_years = year(datetimes) + (month(datetimes) >= start_month);
unique_years = unique(hydro_years);
min_date = min(datetimes);
max_date = max(datetimes);
valid_years = false(size(unique_years));
for i = 1:numel(unique_years)
    y = unique_years(i);
    if start_month == 1
        year_start = datetime(y-1, 1, 1);
        year_end = datetime(y-1, 12, 31);
    else
        year_start = datetime(y-1, start_month, 1);
        year_end = datetime(y, start_month-1, eomday(y, start_month-1));
    end
    valid_years(i) = (min_date <= year_start) && (max_date >= year_end);
end
keep_idx = valid_years;
centroid_years = unique_years(keep_idx);
num_years = numel(centroid_years);
centroid_times = NaN(num_years, 1);
for i = 1:num_years
    y = centroid_years(i);
    mask = hydro_years == y;
    t_year = datetimes(mask);
    Q_year = runoff(mask);
    expected_dates = (t_year(1):days(1):t_year(end))';
    if numel(t_year) ~= numel(expected_dates)
        continue  
    end 
    total_Q = sum(Q_year);
    if total_Q == 0
        continue;
    end 
    cum_Q = cumsum(Q_year);
    target = total_Q * 0.5; 
    k = find(cum_Q >= target, 1);
    if isempty(k), continue; end
    if k == 1
        c_prev = 0;
    else
        c_prev = cum_Q(k-1);
    end 
    delta = (target - c_prev) / Q_year(k);
    centroid_time = k + delta; 
    centroid_times(i) = centroid_time;
end
nan_mask = isnan(centroid_times);
centroid_years(nan_mask) = [];
centroid_times(nan_mask) = [];
if start_month==1
    centroid_years(end) = [];
    centroid_times(end) = [];
end
centroid_years = centroid_years-1;   
end

