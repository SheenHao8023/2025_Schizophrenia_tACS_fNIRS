%% 需要将下面的工具包添加至路径
% easyh5: https://github.com/fangq/easyh5
% jsnirfy: https://github.com/NeuroJSON/jsnirfy/tree/4e7de160cdb9a0357c8e9e7d882a190bf6e2f51a?tab=readme-ov-file
% Homer3: https://github.com/BUNPC/Homer3
% NIRS-KIT: https://github.com/bnuhouxin/NIRS-KIT
% MVGC: https://github.com/lcbarnett/MVGC1

%% 将所有.nirs文件转换为.snirf标准文件
Homer3 % 运行homer3，再关闭窗口或直接转换
Path = 'C:/Users/XinHao/Desktop/tES_SZ_fNIRS/'; 
snirf = Nirs2Snirf('C:/Users/XinHao/Desktop/tES_SZ_fNIRS/1.nirsdata');  % 如果直接通过GUI转换则注释掉本行，此函数仅支持路径全名，不可拼接
% 移动数据，保留原始数据文件不做改动
mkdir(fullfile(Path, '2.snirfdata'));
source_dir = fullfile(Path, '1.nirsdata');
target_dir = fullfile(Path, '2.snirfdata');
snirf_files = dir(fullfile(source_dir, '*.snirf'));
for i = 1:length(snirf_files)
    snirf_file = fullfile(source_dir, snirf_files(i).name);
    target_snirf_file = fullfile(target_dir, snirf_files(i).name);
    movefile(snirf_file, target_snirf_file);
end

%% 数据分段
Path = 'C:/Users/XinHao/Desktop/tES_SZ_fNIRS/'; 
all_files = dir(fullfile(Path, '2.snirfdata', '*.snirf'));
for i = 1:length(all_files)
    current_file = fullfile(all_files(i).folder, all_files(i).name);
    Data = loadsnirf(current_file);
    [~, filename, ~] = fileparts(current_file);
    if endsWith(filename, 'resting') % resting-state including pre and post test
        timeseries = Data.nirs.data.time;
        time_indices = find(timeseries > 30 & timeseries <= 150); % time window
        Data.nirs.data.time = Data.nirs.data.time(time_indices)-30;
        Data.nirs.aux.time = Data.nirs.aux.time(time_indices)-30;
        Data.nirs.data.dataTimeSeries = Data.nirs.data.dataTimeSeries(time_indices, :);
        Data.nirs.aux.dataTimeSeries = Data.nirs.aux.dataTimeSeries(time_indices, :);
        if endsWith(filename, 'Aresting')
            savesnirf(Data, fullfile(all_files(i).folder, [filename(1:5), '_0r', '.snirf']));
        else
            savesnirf(Data, fullfile(all_files(i).folder, [filename(1:5), '_5r', '.snirf']));
        end
    elseif endsWith(filename, 'A') % pre test tasking state including 4 conditions
        timeseries = Data.nirs.data.time;
        time_data = Data.nirs.stim1.data(:, 1);
        suffixes = {'_1e', '_2e', '_3e', '_4e'}; % Hearing each other
        if size(time_data, 1) <= 4
            disp('Not enough stim marks in 4blocks, skipping processing.');
        else
            [idx, centroids] = kmeans(time_data, 4);
            [~, sort_idx] = sort(centroids);
            sorted_idx = zeros(size(idx));
            for j = 1:4
                sorted_idx(idx == sort_idx(j)) = j;
            end
            for j = 1:4
                cluster_times = time_data(sorted_idx == j);
                if size(cluster_times, 1) < 2
                    continue;
                end
                time_indices = find(timeseries > min(cluster_times) & timeseries <= max(cluster_times));
                Data.nirs.data.time = Data.nirs.data.time(time_indices)-min(cluster_times);
                Data.nirs.aux.time = Data.nirs.aux.time(time_indices)-min(cluster_times);
                Data.nirs.data.dataTimeSeries = Data.nirs.data.dataTimeSeries(time_indices, :);
                Data.nirs.aux.dataTimeSeries = Data.nirs.aux.dataTimeSeries(time_indices, :);
                savesnirf(Data, fullfile(all_files(i).folder, [filename(1:5), suffixes{j}, '.snirf']));
                Data = loadsnirf(current_file);
            end
        end

        Data = loadsnirf(current_file);
        timeseries = Data.nirs.data.time;
        time_data = Data.nirs.stim3.data(:, 1);
        suffixes = {'_1a', '_2a', '_3a', '_4a'}; % Hearing A
        if size(time_data, 1) <= 4
            disp('Not enough stim marks in 4blocks, skipping processing.');
        else
            [idx, centroids] = kmeans(time_data, 4);
            [~, sort_idx] = sort(centroids);
            sorted_idx = zeros(size(idx));
            for j = 1:4
                sorted_idx(idx == sort_idx(j)) = j;
            end
            for j = 1:4
                cluster_times = time_data(sorted_idx == j);
                if size(cluster_times, 1) < 2
                    continue;
                end
                time_indices = find(timeseries > min(cluster_times) & timeseries <= max(cluster_times));
                Data.nirs.data.time = Data.nirs.data.time(time_indices)-min(cluster_times);
                Data.nirs.aux.time = Data.nirs.aux.time(time_indices)-min(cluster_times);
                Data.nirs.data.dataTimeSeries = Data.nirs.data.dataTimeSeries(time_indices, :);
                Data.nirs.aux.dataTimeSeries = Data.nirs.aux.dataTimeSeries(time_indices, :);
                savesnirf(Data, fullfile(all_files(i).folder, [filename(1:5), suffixes{j}, '.snirf']));
                Data = loadsnirf(current_file);
            end
        end

        Data = loadsnirf(current_file);
        timeseries = Data.nirs.data.time;
        time_data = Data.nirs.stim4.data(:, 1);
        suffixes = {'_1b', '_2b', '_3b', '_4b'}; % Hearing B
        if size(time_data, 1) <= 4
            disp('Not enough stim marks in 4blocks, skipping processing.');
        else
            [idx, centroids] = kmeans(time_data, 4);
            [~, sort_idx] = sort(centroids);
            sorted_idx = zeros(size(idx));
            for j = 1:4
                sorted_idx(idx == sort_idx(j)) = j;
            end
            for j = 1:4
                cluster_times = time_data(sorted_idx == j);
                if size(cluster_times, 1) < 2
                    continue;
                end
                time_indices = find(timeseries > min(cluster_times) & timeseries <= max(cluster_times));
                Data.nirs.data.time = Data.nirs.data.time(time_indices)-min(cluster_times);
                Data.nirs.aux.time = Data.nirs.aux.time(time_indices)-min(cluster_times);
                Data.nirs.data.dataTimeSeries = Data.nirs.data.dataTimeSeries(time_indices, :);
                Data.nirs.aux.dataTimeSeries = Data.nirs.aux.dataTimeSeries(time_indices, :);
                savesnirf(Data, fullfile(all_files(i).folder, [filename(1:5), suffixes{j}, '.snirf']));
                Data = loadsnirf(current_file);
            end
        end

    else % post test tasking state including 1 condition
        timeseries = Data.nirs.data.time;
        stim_times = {'stim1', 'stim3', 'stim4'};  
        suffixes = {'_6e', '_6a', '_6b'}; 
        for j = 1:3
            time_data = Data.nirs.(stim_times{j}).data(:, 1);
            if size(time_data, 1) < 2
                continue;
            end
            time_indices = find(timeseries > min(time_data) & timeseries <= max(time_data));
            Data.nirs.data.time = Data.nirs.data.time(time_indices) - min(time_data);
            Data.nirs.aux.time = Data.nirs.aux.time(time_indices) - min(time_data);
            Data.nirs.data.dataTimeSeries = Data.nirs.data.dataTimeSeries(time_indices, :);
            Data.nirs.aux.dataTimeSeries = Data.nirs.aux.dataTimeSeries(time_indices, :);
            savesnirf(Data, fullfile(all_files(i).folder, [filename(1:5), suffixes{j}, '.snirf']));
            Data = loadsnirf(current_file);
        end
    end
    delete(current_file)
end

%% 数据预处理
Path = 'C:/Users/XinHao/Desktop/tES_SZ_fNIRS/'; 
all_files = dir(fullfile(Path, '2.snirfdata', '*.snirf'));
mkdir(fullfile(Path, '3.preprocessed'));
% 随机选择一个文件获取近红外SD信息，用于适配相关函数
nirs_files = dir(fullfile(Path, '1.nirsdata', '*.nirs'));
nirs = NirsClass(fullfile(Path, '1.nirsdata', nirs_files(randi(length(nirs_files))).name));
for i = 1:length(all_files)
    current_file = fullfile(all_files(i).folder, all_files(i).name);
    data = DataClass(current_file);
    dod = hmrR_Intensity2OD(data);
    dod = hmrR_BandpassFilt(dod, 0.01, 0.1); % band-pass filter
    dod_tddr = hmrR_MotionCorrectTDDR(dod, nirs.SD, 50); % motion correction: fs = 50，此函数在NIRS_KIT包所带TDDR函数基础上修改
    % Fishburn, F. A. et al. (2019). Temporal Derivative Distribution Repair (TDDR): A motion correction method for fNIRS. NeuroImage, 184, 171-179.
    dod.dataTimeSeries = dod_tddr;
    dc = hmrOD2Conc(dod.dataTimeSeries, nirs.SD, [1 1 1]); % ppf = [1 1 1]
    nirsdata.oxyData = squeeze(dc(:, 1, :)); % 只保留HbO氧合血红蛋白的数据
    % Luke, R. et al. (2021). Oxygenated hemoglobin signal provides greater predictive performance of experimental condition than de-oxygenated. BioRxiv, 2021-11.
    tp=size(nirsdata.oxyData,1);
    for ch=1:size(nirsdata.oxyData,2)
        p_oxy=polyfit((1:tp)',nirsdata.oxyData(:,ch), 1); % 一次多项式拟合去趋势
        base_oxy=polyval(p_oxy,1:tp);
        nirsdata.oxyData(:,ch)=nirsdata.oxyData(:,ch)-base_oxy';
    end
    nirsdata.nch = 120;
    nirsdata.T = 0.02;
    [~, filename, ~] = fileparts(current_file);
    save(fullfile(Path, '3.preprocessed', [filename, '.mat']), 'nirsdata');
end

%% 建立ROI
Path = 'C:/Users/haox8/Desktop/tES_SZ_fNIRS/'; 
all_files = dir(fullfile(Path, '3.preprocessed', '*.mat'));
mkdir(fullfile(Path, '4.roi')); 
for i = 1:length(all_files)
    current_file = fullfile(all_files(i).folder, all_files(i).name);
    load(current_file); 
    nirsdata.oxyData = zscore(nirsdata.oxyData);
    newdata = zeros(size(nirsdata.oxyData, 1), 30); 
    % channel 1~60 --> 1~15 for SZ
    newdata(:, 1) = mean(nirsdata.oxyData(:, [1 2]), 2, 'omitnan'); % Right FG
    newdata(:, 2) = mean(nirsdata.oxyData(:, 6), 2, 'omitnan');   % Right STG
    newdata(:, 3) = mean(nirsdata.oxyData(:, [8 10 11]), 2, 'omitnan'); % Right TPJ
    newdata(:, 4) = mean(nirsdata.oxyData(:, [14 15 16]), 2, 'omitnan');     % Right PMC
    newdata(:, 5) = mean(nirsdata.oxyData(:, [17 29]), 2, 'omitnan');     % Right SFG
    newdata(:, 6) = mean(nirsdata.oxyData(:, [18 19 22 24 25]), 2, 'omitnan'); % Right DLPFC
    newdata(:, 7) = mean(nirsdata.oxyData(:, [26 27 28]), 2, 'omitnan'); % Right FPC
    newdata(:, 8) = mean(nirsdata.oxyData(:, [30 32]), 2, 'omitnan');     % DMPFC
    newdata(:, 9) = mean(nirsdata.oxyData(:, [37 40 41]), 2, 'omitnan');     % Left FPC
    newdata(:, 10) = mean(nirsdata.oxyData(:, [36 38 51 52]), 2, 'omitnan');  % Left DLPFC
    newdata(:, 11) = mean(nirsdata.oxyData(:, [31 50]), 2, 'omitnan');    % Left SFG
    newdata(:, 12) = mean(nirsdata.oxyData(:, [45 49]), 2, 'omitnan'); % Left PMC
    newdata(:, 13) = mean(nirsdata.oxyData(:, [42 43]), 2, 'omitnan'); % Left TPJ
    newdata(:, 14) = mean(nirsdata.oxyData(:, 54), 2, 'omitnan'); % Left STG
    newdata(:, 15) = mean(nirsdata.oxyData(:, 57), 2, 'omitnan'); % Left FG
    % channel 61~120 --> 16~30 for HC
    newdata(:, 16) = mean(nirsdata.oxyData(:, [61 62]), 2, 'omitnan'); % Right FG
    newdata(:, 17) = mean(nirsdata.oxyData(:, 66), 2, 'omitnan');   % Right STG
    newdata(:, 18) = mean(nirsdata.oxyData(:, [68 70 71]), 2, 'omitnan'); % Right TPJ
    newdata(:, 19) = mean(nirsdata.oxyData(:, [74 75 76]), 2, 'omitnan');     % Right PMC
    newdata(:, 20) = mean(nirsdata.oxyData(:, [77 89]), 2, 'omitnan');     % Right SFG
    newdata(:, 21) = mean(nirsdata.oxyData(:, [78 79 82 84 85]), 2, 'omitnan'); % Right DLPFC
    newdata(:, 22) = mean(nirsdata.oxyData(:, [86 87 88]), 2, 'omitnan'); % Right FPC
    newdata(:, 23) = mean(nirsdata.oxyData(:, [90 92]), 2, 'omitnan');     % DMPFC
    newdata(:, 24) = mean(nirsdata.oxyData(:, [97 100 101]), 2, 'omitnan');     % Left FPC
    newdata(:, 25) = mean(nirsdata.oxyData(:, [96 98 111 112]), 2, 'omitnan');  % Left DLPFC
    newdata(:, 26) = mean(nirsdata.oxyData(:, [91 110]), 2, 'omitnan');    % Left SFG
    newdata(:, 27) = mean(nirsdata.oxyData(:, [105 109]), 2, 'omitnan'); % Left PMC
    newdata(:, 28) = mean(nirsdata.oxyData(:, [102 103]), 2, 'omitnan'); % Left TPJ
    newdata(:, 29) = mean(nirsdata.oxyData(:, 114), 2, 'omitnan'); % Left STG
    newdata(:, 30) = mean(nirsdata.oxyData(:, 117), 2, 'omitnan'); % Left FG
    nirsdata.oxyData = newdata;
    nirsdata.nch=30;
    save(fullfile(Path, '4.roi', all_files(i).name), 'nirsdata');
end

%% 脑激活分析，GLM
Path = 'C:/Users/XinHao/Desktop/tES_SZ_fNIRS/'; 
all_files = dir(fullfile(Path, '4.roi', '*.mat'));
subjects = unique(arrayfun(@(x) x.name(1:5), all_files, 'UniformOutput', false));
num_subjects = length(subjects);
task_conditions = {'0r', '1a', '1b', '1e', '2a', '2b', '2e', '3a', '3b', '3e', '4a', '4b', '4e', '5r', '6a', '6b', '6e'};
result_matrix = NaN(num_subjects, 2 + 17 * 30);
group1 = {'01181', '12281', '03011'}; % Experimental
group2 = {}; % ControlActive
group3 = {}; % ControlSham
group4 = {'02061'}; % ControlResting
for sub = 1:num_subjects
    subject_files = all_files(contains({all_files.name}, subjects{sub}(1:5)));
    subject_data = NaN(1, 2 + 17 * 30); 
    subject_data(1) = str2double(subjects{sub}); % 存ID
    subject_id = subjects{sub};
    if ismember(subject_id, group1)
        subject_data(2) = 1; % 存组号
    elseif ismember(subject_id, group2)
        subject_data(2) = 2;
    elseif ismember(subject_id, group3)
        subject_data(2) = 3;
    elseif ismember(subject_id, group4)
        subject_data(2) = 4;
    else
        subject_data(2) = NaN; % 不属于任何组的被试
    end
    for ch = 1:30
        Y = [];
        X = []; % design matrix，可考虑修改
        for file_idx = 1:length(subject_files)
            dcdata = load(fullfile(subject_files(file_idx).folder, subject_files(file_idx).name));
            oxy_data = dcdata.nirsdata.oxyData(:, ch); 
            task_label = subject_files(file_idx).name(7:8); 
            task_idx = find(strcmp(task_conditions, task_label));
            if ~isempty(task_idx)
                task_regressor = zeros(length(oxy_data), 17);
                task_regressor(:, task_idx) = 1;
                Y = [Y; oxy_data];
                X = [X; ones(length(oxy_data), 1), task_regressor]; 
            end
        end
        betas = regress(Y, X);
        start_col = 2 + 17 * (ch - 1); 
        betas = 2 ./ (1 + exp(-0.1 * (zscore(betas) - mean(zscore(betas))))) - 1; %Z-score+Sigmoid
        subject_data(1, start_col + 1:start_col + 17) = betas(2:end); %不需要常数项
    end
    result_matrix(sub, :) = subject_data;
end
mkdir(fullfile(Path, '5.glm')); 
xlswrite(fullfile(Path, '5.glm', 'activation.xlsx'), result_matrix);

%% 小波相干性功能连接分析
Path = 'C:/Users/XinHao/Desktop/tES_SZ_fNIRS/'; 
all_files = dir(fullfile(Path, '4.roi', '*.mat'));
coherence_matrix = zeros(15, 15);
mkdir(fullfile(Path, '6.coh_sz')); 
mkdir(fullfile(Path, '6.coh_pair')); 
for i = 1:length(all_files)
    current_file = fullfile(all_files(i).folder, all_files(i).name);
    load(current_file);
    for m = 1:15
        for n = 1:15
            [wcoh,~,f] = wcoherence(nirsdata.oxyData(:, m), nirsdata.oxyData(:, n));
            coherence_value = wcoh (0.01<=f(:,1) & f(:,1)<=0.04 , :);
            coherence_matrix(m,n) = mean(coherence_value(:));
        end
    end
    save(fullfile(Path, '6.coh_sz', all_files(i).name), 'coherence_matrix');
    for p = 1:15
        for q = 16:30
            [wcoh,~,f] = wcoherence(nirsdata.oxyData(:, p), nirsdata.oxyData(:, q));
            coherence_value = wcoh (0.01<=f(:,1) & f(:,1)<=0.04 , :);
            coherence_matrix(p, q-15) = mean(coherence_value(:));
        end
    end
    save(fullfile(Path, '6.coh_pair', all_files(i).name), 'coherence_matrix');
end

%% 格兰杰因果有效连接分析
Path = 'C:/Users/XinHao/Desktop/tES_SZ_fNIRS/'; 
all_files = dir(fullfile(Path, '4.roi', '*.mat'));
mkdir(fullfile(Path, '7.gc_sz')); 
mkdir(fullfile(Path, '7.gc_pairab')); 
mkdir(fullfile(Path, '7.gc_pairba')); 
% state-space method, time domain Granger Causality Analysis
% Barnett, L., & Seth, A. K. (2015). Granger causality for state-space models. Physical Review E, 91(4), 040101.
ntrials   = 1;     % number of trials
regmode   = 'LWR';  % or 'OLS'
icregmode = 'LWR';  % or 'OLS'
morder    = 'BIC';  % or 'actual', 'AIC', numerical value)
momax     = 20;     % maximum model order for model order estimation
acmaxlags = 100;   % maximum autocovariance lags (empty for automatic calculation) or '1000'
tstat     = 'F';    % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR'; % 'NONE' 'BONFERRONI' 'SIDAK' 'HOLM' 'FDR' 'FDRD'
for i = 1:length(all_files)
    current_file = fullfile(all_files(i).folder, all_files(i).name);
    load(current_file);
    nobs = size(nirsdata.oxyData, 1);   % number of observations per trial
    nvars = size(nirsdata.oxyData, 2);   % number of variables
    [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(nirsdata.oxyData',momax,icregmode);  % Model order estimation
    if     strcmpi(morder,'actual')
        morder = amo;
        fprintf('\nusing actual model order = %d\n',morder);
    elseif strcmpi(morder,'AIC')
        morder = moAIC;
        fprintf('\nusing AIC best model order = %d\n',morder);
    elseif strcmpi(morder,'BIC')
        morder = moBIC;
        fprintf('\nusing BIC best model order = %d\n',morder);
    else
        fprintf('\nusing specified model order = %d\n',morder);
    end
    [A, SIG] = tsdata_to_var(nirsdata.oxyData', morder, regmode);  % VAR model estimation
    [F, ~] = var_to_pwcgc(A, SIG, nirsdata.oxyData', regmode, tstat);  % Granger causality calculation: time domain 
    %方向性：F(i, j) 代表从变量 i 到变量 j 的 Granger 因果性
    F_sz = F(1:15, 1:15);
    F_pairab = F(1:15, 16:30);
    F_pairba = F(16:30, 1:15);
    save(fullfile(Path, '7.gc_sz', all_files(i).name), 'F_sz'); %participant A
    save(fullfile(Path, '7.gc_pairab', all_files(i).name), 'F_pairab');
    save(fullfile(Path, '7.gc_pairba', all_files(i).name), 'F_pairba');
end

%% 自定义函数: hmrR_MotionCorrectTDDR
function [dodTDDR] = hmrR_MotionCorrectTDDR(dod, SD, sample_rate)
    mlAct = SD.MeasListAct; % prune bad channels
    lstAct = find(mlAct==1);
    dodTDDR = dod.dataTimeSeries;
    for ii=1:length(lstAct)
        idx_ch = lstAct(ii);
        filter_cutoff = .5;
        filter_order = 3;
        Fc = filter_cutoff * 2/sample_rate;
        if Fc<1
            [fb,fa] = butter(filter_order,Fc);
            signal_low = filtfilt(fb,fa,dod.dataTimeSeries(:,idx_ch));
        else
            signal_low = dod.dataTimeSeries(:,idx_ch);
        end
        signal_high = dod.dataTimeSeries(:,idx_ch) - signal_low;
        tune = 4.685;
        D = sqrt(eps(class(dod.dataTimeSeries)));
        mu = inf;
        iter = 0;
        deriv = diff(signal_low);
        w = ones(size(deriv));
        while iter < 50
            iter = iter + 1;
            mu0 = mu;
            mu = sum( w .* deriv ) / sum( w );
            dev = abs(deriv - mu);
            sigma = 1.4826 * median(dev);
            r = dev / (sigma * tune);
            w = ((1 - r.^2) .* (r < 1)) .^ 2;
            if abs(mu-mu0) < D*max(abs(mu),abs(mu0))
                break;
            end
        end
        new_deriv = w .* (deriv-mu);
        signal_low_corrected = cumsum([0; new_deriv]);
        signal_low_corrected = signal_low_corrected - mean(signal_low_corrected);
        dodTDDR(:,idx_ch) = signal_low_corrected + signal_high;
    end
end