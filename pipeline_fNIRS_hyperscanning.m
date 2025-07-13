%% 需要将下面的工具包添加至路径
% easyh5: https://github.com/fangq/easyh5
% jsnirfy: https://github.com/NeuroJSON/jsnirfy/tree/4e7de160cdb9a0357c8e9e7d882a190bf6e2f51a?tab=readme-ov-file
% Homer3: https://github.com/BUNPC/Homer3
% NIRS-KIT: https://github.com/bnuhouxin/NIRS-KIT
% MVGC: https://github.com/lcbarnett/MVGC1
% JIDT: https://github.com/jlizier/jidt,   运行此工具包需要安装Java
% GRETNA: https://github.com/sandywang/GRETNA

%% 将所有.nirs文件转换为.snirf标准文件
Homer3 % 运行homer3，弹出后选择no不转换，再关闭Homer3 GUI
Path = 'C:/Users/XinHao/Desktop/tES_SZ_fNIRS/'; 
snirf = Nirs2Snirf('C:/Users/XinHao/Desktop/tES_SZ_fNIRS/1.nirsdata');  % 如果直接通过GUI转换则注释掉本行，此函数仅支持路径全名，不可拼接
% 移动数据，保留原始数据文件不做改动
mkdir(fullfile(Path, '1.snirfdata_full'));
mkdir(fullfile(Path, '1.snirfdata_split'));
source_dir = fullfile(Path, '1.nirsdata');
target_dir = fullfile(Path, '1.snirfdata_full');
snirf_files = dir(fullfile(source_dir, '*.snirf'));
for i = 1:length(snirf_files)
    snirf_file = fullfile(source_dir, snirf_files(i).name);
    target_snirf_file = fullfile(target_dir, snirf_files(i).name);
    movefile(snirf_file, target_snirf_file);
end

% 获取已处理的文件前缀
% processed_files = dir(fullfile(Path, '1.snirfdata_split', '*.snirf'));
% processed_names = string({processed_files.name});
% processed_base_names = extractBefore(processed_names, '_');  % 获取文件前缀用于比较

%% 数据分段
Path = 'C:/Users/XinHao/Desktop/tES_SZ_fNIRS/'; 
all_files = dir(fullfile(Path, '1.snirfdata_full', '*.snirf'));
oldapparatus = {'HZ001', 'HZ004', 'HZ005', 'HZ006', 'HZ007', 'SX001', 'SX002', 'SX003', 'SX004', 'SX005', 'SX006', 'SX007', 'SX008'}; 
% Only Hearing each other needs to distinguish markers in different apparatus
newapparatus = {'HZ002', 'HZ003'};
for i = 1:length(all_files)
    current_file = fullfile(all_files(i).folder, all_files(i).name);
    Data = loadsnirf(current_file);
    [~, filename, ~] = fileparts(current_file);

    % if any(contains(processed_base_names, filename(1:5)))
    %     disp(['Skipping already processed file: ', filename]); % 如果文件已处理，则跳过
    %     continue;
    % end

    if endsWith(filename, 'resting') % resting-state including pre- and post-test
        timeseries = Data.nirs.data.time;
        time_indices = find(timeseries > 30 & timeseries <= 150); % time window
        Data.nirs.data.time = Data.nirs.data.time(time_indices)-30;
        Data.nirs.aux.time = Data.nirs.aux.time(time_indices)-30;
        Data.nirs.data.dataTimeSeries = Data.nirs.data.dataTimeSeries(time_indices, :);
        Data.nirs.aux.dataTimeSeries = Data.nirs.aux.dataTimeSeries(time_indices, :);
        if endsWith(filename, '1_resting')
            savesnirf(Data, fullfile(Path, '1.snirfdata_split', [filename(1:5), '_R1', '.snirf']));
        else
            savesnirf(Data, fullfile(Path, '1.snirfdata_split', [filename(1:5), '_R2', '.snirf']));
        end

    elseif  endsWith(filename, 'tasking') && filename(7) == '1' % pre test tasking state including 4 conditions
        timeseries = Data.nirs.data.time;
        if ismember(filename(1:5), oldapparatus)
            time_data = Data.nirs.stim2.data(:, 1);
        else
            time_data = Data.nirs.stim1.data(:, 1);  %newapparatus
        end
        suffixes = {'_HEO1', '_HEO2', '_HEO3', '_HEO4'}; % Hearing each other
        if size(time_data, 1) <= 4
            disp('Not enough stim markers in tasking state, skipping processing.');
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
                savesnirf(Data, fullfile(Path, '1.snirfdata_split', [filename(1:5), suffixes{j}, '.snirf']));
                Data = loadsnirf(current_file);
            end
        end

        Data = loadsnirf(current_file);
        timeseries = Data.nirs.data.time;
        time_data = Data.nirs.stim3.data(:, 1);
        suffixes = {'_HA1', '_HA2', '_HA3', '_HA4'}; % Hearing A/HC
        if size(time_data, 1) <= 4
            disp('Not enough stim markers in tasking state, skipping processing.');
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
                savesnirf(Data, fullfile(Path, '1.snirfdata_split', [filename(1:5), suffixes{j}, '.snirf']));
                Data = loadsnirf(current_file);
            end
        end

        Data = loadsnirf(current_file);
        timeseries = Data.nirs.data.time;
        time_data = Data.nirs.stim4.data(:, 1);
        suffixes = {'_HB1', '_HB2', '_HB3', '_HB4'}; % Hearing B/SZ
        if size(time_data, 1) <= 4
            disp('Not enough stim markers in tasking state, skipping processing.');
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
                savesnirf(Data, fullfile(Path, '1.snirfdata_split', [filename(1:5), suffixes{j}, '.snirf']));
                Data = loadsnirf(current_file);
            end
        end

    else % post test tasking state including 1 condition
        timeseries = Data.nirs.data.time;
        suffixes = {'_HEO5', '_HA5', '_HB5'}; 
        if ismember(filename(1:5), oldapparatus)
            stim_times = {'stim2', 'stim3', 'stim4'};  
        else
            stim_times = {'stim1', 'stim3', 'stim4'};    %newapparatus
        end
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
            savesnirf(Data, fullfile(Path, '1.snirfdata_split',  [filename(1:5), suffixes{j}, '.snirf']));
            Data = loadsnirf(current_file);
        end
    end
end

%% 数据预处理
Path = 'C:/Users/XinHao/Desktop/tES_SZ_fNIRS/'; 
all_files = dir(fullfile(Path, '1.snirfdata_split', '*.snirf'));
mkdir(fullfile(Path, '2.preprocessed'));
% 随机选择一个文件获取近红外SD信息，用于适配相关函数
nirs_files = dir(fullfile(Path, '1.nirsdata', '*.nirs'));
nirs = NirsClass(fullfile(Path, '1.nirsdata', nirs_files(randi(length(nirs_files))).name));
for i = 1:length(all_files)
    current_file = fullfile(all_files(i).folder, all_files(i).name);
    data = DataClass(current_file);
    [~, filename, ~] = fileparts(current_file);

    % if any(contains(processed_base_names, extractBefore(filename, '_')))
    %     disp(['Skipping already preprocessed file: ', filename]);
    %     continue;
    % end

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
    save(fullfile(Path, '2.preprocessed', [filename, '.mat']), 'nirsdata');
end

%% 建立ROI
Path = 'C:/Users/XinHao/Desktop/tES_SZ_fNIRS/'; 
all_files = dir(fullfile(Path, '2.preprocessed', '*.mat'));
mkdir(fullfile(Path, '3.roi')); 
for i = 1:length(all_files)

    % [~, base_name, ~] = fileparts(all_files(i).name);
    % if any(processed_base_names == base_name)
    %     disp(['Skipping already processed file: ', base_name]);
    %     continue;
    % end

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
    save(fullfile(Path, '3.roi', all_files(i).name), 'nirsdata');
end

%% 脑激活分析，基于GLM
Path = 'C:/Users/XinHao/Desktop/tES_SZ_fNIRS/'; 

% 判断哪些 subject 已经处理
% result_file = fullfile(Path, '4.glm', 'activation.xlsx');
% if exist(result_file, 'file')
%     old_data = readcell(result_file);
%     old_subjects = string(old_data(2:end, 1)); % 第一列是ID，跳过表头
% else
%     old_data = {};
%     old_subjects = string([]); % 空数组
% end

all_files = dir(fullfile(Path, '3.roi', '*.mat'));
subjects = unique(arrayfun(@(x) x.name(1:5), all_files, 'UniformOutput', false));
num_subjects = length(subjects);
task_conditions = {'R1', 'HA1', 'HB1', 'HEO1', 'HA2', 'HB2', 'HEO2', 'HA3', 'HB3', 'HEO3', 'HA4', 'HB4', 'HEO4', 'R2', 'HA5', 'HB5', 'HEO5'};
result_matrix = cell(num_subjects, 2 + 17 * 30);
% 用于卷积的血流动力学函数 HRF
[hrf, ~] = spm_hrf(0.02); % 采样周期 T=0.02
group1 = {'HZ001', 'HZ002', 'HZ004', 'HZ005', 'HZ006', 'SX001', 'SX007', 'SX008', 'SX009', 'SX013', 'SX016', 'SX018'}; % Experimental
group2 = {'SX015', 'SX017', 'SX019'}; % ControlActive
group3 = {'HZ007', 'SX002', 'SX010', 'SX011', 'SX012'}; % ControlSham
group4 = {'HZ003', 'SX006'}; % ControlResting
group5 = {'SX003', 'SX004', 'SX005', 'SX014'}; % ControlBehavior
for sub = 1:num_subjects
    subject_id = subjects{sub};

    % if any(old_subjects == subject_id)
    %     disp(['Skipping already processed subject: ', subject_id]);
    %     continue; % 跳过已处理 subject
    % end

    subject_files = all_files(contains({all_files.name}, subjects{sub}));
    get_condition = @(condition_name) extractBetween(condition_name, '_', '.mat'); % 从文件名中提取条件部分
    condition_names = {subject_files.name}; 
    [~, sorted_idx] = sort(cellfun(@(x) find(strcmp(task_conditions, get_condition(x))), condition_names));
    subject_files = subject_files(sorted_idx);
    subject_data = cell(1, 2 + 17 * 30); 
    subject_data{1} = subjects{sub}; % 存ID
    
    if ismember(subject_id, group1)
        subject_data{2} = 1; % 存组号
    elseif ismember(subject_id, group2)
        subject_data{2} = 2;
    elseif ismember(subject_id, group3)
        subject_data{2} = 3;
    elseif ismember(subject_id, group4)
        subject_data{2} = 4;
    elseif ismember(subject_id, group5)
        subject_data{2} = 5;
    else
        subject_data{2} = NaN; % 不属于任何组
    end
    for ch = 1:30
        Y = [];
        X = []; % 二值化0/1的方波设计矩阵，block设计用
        X_hrf = [];
        for file_idx = 1:length(subject_files)
            dcdata = load(fullfile(subject_files(file_idx).folder, subject_files(file_idx).name));
            oxy_data = dcdata.nirsdata.oxyData(:, ch); 
            task_label = subject_files(file_idx).name(7:strfind(subject_files(file_idx).name, '.mat')-1);
            task_idx = find(strcmp(task_conditions, task_label));
            if ~isempty(task_idx)
                task_regressor = zeros(length(oxy_data), 17);
                task_regressor(:, task_idx) = 1;
                Y = [Y; oxy_data];
                X = [X; task_regressor]; % 不包含截距项
            else
                disp([subject_files(file_idx).name, ' is empty']); 
            end
        end
        for col_idx = 1:size(X, 2)  % 对设计矩阵每一列进行卷积处理
            current_col = X(:, col_idx);
            if any(current_col == 1) % 如果当前列不全为 0，才进行卷积
                start_idx = find(current_col == 1, 1, 'first');  % 第一个 1 的位置
                end_idx = find(current_col == 1, 1, 'last');    % 最后一个 1 的位置
                segment = current_col(start_idx:end_idx);  % 提取连续为 1 的区间，对这一段进行卷积
                conv_result = conv(segment, hrf, 'same');
                X_hrf(start_idx:end_idx, col_idx) = conv_result(1:length(segment));
            else % 如果当前列全为 0（如这个条件的数据缺失），则保持 X_hrf 中的该列为 0
                X_hrf(:, col_idx) = 0;  
            end
        end
        %X_final = [ones(length(X), 1), X_hrf];  % 卷积后再添加截距项
        X_final = X_hrf;  % 卷积后设计矩阵不添加截距项
        MDL = fitglm(X_final, Y);
        betas = MDL.Coefficients.Estimate';
        start_col = 2 + 17 * (ch - 1); 
        new_betas = [betas(2:end)]; %不需要储存截距项系数
        subject_data(1, start_col + 1:start_col + 17) = num2cell(new_betas); 
        for i = 1:numel(subject_data)
            if isequal(subject_data{i}, 0)  
                subject_data{i} = NaN; 
            end
        end
    end
    result_matrix(sub, :) = subject_data;
end
mkdir(fullfile(Path, '4.glm')); 
xlswrite(fullfile(Path, '4.glm', 'activation.xlsx'), result_matrix);

%% wavelet coherence 功能连接分析
Path = 'C:/Users/XinHao/Desktop/tES_SZ_fNIRS/'; 
all_files = dir(fullfile(Path, '3.roi', '*.mat'));
coherence_matrix = zeros(15, 15);
mkdir(fullfile(Path, '5.coh_sz')); 
mkdir(fullfile(Path, '5.coh_pair')); 
for i = 1:length(all_files)
    current_file = fullfile(all_files(i).folder, all_files(i).name);
    load(current_file);
    for m = 1:15
        for n = 1:15
            [wcoh,~,f] = wcoherence(nirsdata.oxyData(:, m), nirsdata.oxyData(:, n));
            coherence_value = wcoh (0.01<=f(:,1) & f(:,1)<=0.04 , :); % Morlet小波基
            coherence_matrix(m,n) = mean(coherence_value(:));
        end
    end
    coherence_matrix(eye(size(coherence_matrix)) == 1) = 0;  % 将主对角线元素设为0，符合BCT要求
    save(fullfile(Path, '5.coh_sz', all_files(i).name), 'coherence_matrix');
    for p = 1:15
        for q = 16:30
            [wcoh,~,f] = wcoherence(nirsdata.oxyData(:, p), nirsdata.oxyData(:, q));
            coherence_value = wcoh (0.01<=f(:,1) & f(:,1)<=0.04 , :);
            coherence_matrix(p, q-15) = mean(coherence_value(:));
        end
    end
    save(fullfile(Path, '5.coh_pair', all_files(i).name), 'coherence_matrix');
end

%% Granger causality 有效连接分析
Path = 'C:/Users/XinHao/Desktop/tES_SZ_fNIRS/'; 
all_files = dir(fullfile(Path, '3.roi', '*.mat'));
mkdir(fullfile(Path, '5.gc_sz')); 
mkdir(fullfile(Path, '5.gc_pairab')); 
mkdir(fullfile(Path, '5.gc_pairba')); 
% state-space method, time domain Granger Causality Analysis
% 由于神经数据的非平稳性，使用状态空间法研究非平稳数据中所蕴含的时变因果连接变化
% Cekic, S., Grandjean, D., & Renaud, O. (2018). Time, frequency, and time‐varying Granger‐causality measures in neuroscience. Statistics in medicine, 37(11), 1910-1931.
% Barnett, L., & Seth, A. K. (2015). Granger causality for state-space models. Physical Review E, 91(4), 040101.
for i = 1:length(all_files)
    current_file = fullfile(all_files(i).folder, all_files(i).name);
    load(current_file);
    try
        [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(nirsdata.oxyData', 6, 'LWR', false);  % Model order estimation, 'LWR' or 'OLS'
        fprintf('文件 %s, 最优滞后为:%s', all_files(i).name, moBIC);
        % BIC准则比AIC准则更适合时间序列分析
        % Seth, A. K. (2010). A MATLAB toolbox for Granger causal connectivity analysis. Journal of neuroscience methods, 186(2), 262-273.
        [A, SIG] = tsdata_to_var(nirsdata.oxyData', moBIC, 'LWR');  % VAR model estimation, or  'AIC', numerical value)
        [F, ~] = var_to_pwcgc(A, SIG, nirsdata.oxyData', 'LWR', 'F');  % Granger causality calculation: time domain 
        % 方向：F(i, j) 代表从变量 i 到变量 j 的 Granger 因果性
        F_sz = F(1:15, 1:15);
        F_sz(eye(size(F_sz)) == 1) = 0;  % 将主对角线元素设为0，符合BCT要求
        F_pairba = F(1:15, 16:30);
        F_pairab = F(16:30, 1:15);
        save(fullfile(Path, '5.gc_sz', all_files(i).name), 'F_sz'); %participant SZ=B
        save(fullfile(Path, '5.gc_pairba', all_files(i).name), 'F_pairba'); % SZ->HC
        save(fullfile(Path, '5.gc_pairab', all_files(i).name), 'F_pairab'); % SZ<-HC
    catch ME
        % 如果错误包含 "DARE ERROR" 或 "矩阵必须为正定矩阵"，跳过当前文件
        % 残差协方差矩阵显示为非正定，可能表明自协方差序列估计所基于的时间序列数据不够平稳、长度不够、具有共线性或具有高度偏斜的分布，暂不处理
        % 此报错的解释：https://users.sussex.ac.uk/~lionelb/MVGC/html/tsdata_to_infocrit.html
        if contains(ME.message, 'DARE ERROR') || contains(ME.message, '矩阵必须为正定矩阵')
            fprintf('跳过文件 %s, 发生错误: %s\n', all_files(i).name, ME.message);
            continue; % 继续处理下一个文件
        else
            rethrow(ME); % 其他错误则抛出，防止忽略关键问题
        end
    end

end

%% mutual information 功能连接分析
Path = 'C:/Users/XinHao/Desktop/tES_SZ_fNIRS/'; 
all_files = dir(fullfile(Path, '3.roi', '*.mat'));
mi_matrix = zeros(15, 15);
mkdir(fullfile(Path, '6.mi_sz')); 
mkdir(fullfile(Path, '6.mi_pair')); 
javaaddpath('D:\Matlab2024a\matlab_jidt\infodynamics.jar');
miCalc = javaObject('infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1');
% 使用KSG最近邻方法估计互信息
% Lizier, J. T. (2014). JIDT: An information-theoretic toolkit for studying the dynamics of complex systems. Frontiers in Robotics and AI, 1, 11.
miCalc.initialise(1, 1);  
for i = 1:length(all_files)
    current_file = fullfile(all_files(i).folder, all_files(i).name);
    load(current_file);
    try
        for m = 1:15
            for n = 1:15
                miCalc.setObservations(octaveToJavaDoubleArray(nirsdata.oxyData(:, m)), octaveToJavaDoubleArray(nirsdata.oxyData(:, n)));
                mi_matrix(m,n) = miCalc.computeAverageLocalOfObservations();
            end
        end
        mi_matrix(eye(size(mi_matrix)) == 1) = 0;  % 将主对角线元素设为0，符合BCT要求
        save(fullfile(Path, '6.mi_sz', all_files(i).name), 'mi_matrix');
        for p = 1:15
            for q = 16:30
                miCalc.setObservations(octaveToJavaDoubleArray(nirsdata.oxyData(:, p)), octaveToJavaDoubleArray(nirsdata.oxyData(:, q)));
                mi_matrix(p, q-15) = miCalc.computeAverageLocalOfObservations();
            end
        end
        save(fullfile(Path, '6.mi_pair', all_files(i).name), 'mi_matrix');
    catch ME
        % 检测 Java 数组越界异常，跳过当前文件
        if contains(ME.message, 'java.lang.ArrayIndexOutOfBoundsException')
            fprintf('文件 %s 发生 ArrayIndexOutOfBoundsException，跳过此文件。\n', all_files(i).name);
            continue;
        else
            rethrow(ME); % 其他错误仍然报错
        end
    end
end

%% Transfer entropy 有效连接分析
Path = 'C:/Users/XinHao/Desktop/tES_SZ_fNIRS/'; 
all_files = dir(fullfile(Path, '3.roi', '*.mat'));
te_matrix = zeros(15, 15);
te_matrixA = zeros(15, 15);
te_matrixB = zeros(15, 15);
mkdir(fullfile(Path, '6.te_sz')); 
mkdir(fullfile(Path, '6.te_pairab')); 
mkdir(fullfile(Path, '6.te_pairba')); 
javaaddpath('D:\Matlab2024a\matlab_jidt\infodynamics.jar');
teCalc = javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov');
teCalc.initialise(1);  % 设置历史长度lag为 1 (Schreiber k=1)
teCalc.setProperty('k', '4');  % 设置 Kraskov 方法的参数 K = 4（使用 4 个最近邻点）
teCalc.setProperty("NOISE_LEVEL_TO_ADD", "0"); % 不额外添加扰动
for i = 1:length(all_files)
    current_file = fullfile(all_files(i).folder, all_files(i).name);
    load(current_file);
    for m = 1:15
        for n = 1:15
            teCalc.setObservations(nirsdata.oxyData(:, m), nirsdata.oxyData(:, n));
            te_matrix(m, n) = teCalc.computeAverageLocalOfObservations();
        end
    end
    te_matrix(eye(size(te_matrix)) == 1) = 0;  % 将主对角线元素设为0，符合BCT要求
    save(fullfile(Path, '6.te_sz', all_files(i).name), 'te_matrix'); % participant SZ=B

    for p = 1:15
        for q = 16:30
            teCalc.setObservations(nirsdata.oxyData(:, p), nirsdata.oxyData(:, q)); % SZ -> HC
            te_matrixA(p, q-15) = teCalc.computeAverageLocalOfObservations();
            teCalc.setObservations(nirsdata.oxyData(:, q), nirsdata.oxyData(:, p)); % SZ <- HC
            te_matrixB(p, q-15) = teCalc.computeAverageLocalOfObservations();
        end
    end
    save(fullfile(Path, '6.te_pairba', all_files(i).name), 'te_matrixA');
    save(fullfile(Path, '6.te_pairab', all_files(i).name), 'te_matrixB');
end

%% 图论指标分析
% 仅使用coherence matrix of SCZ做节点指标
Path = 'C:/Users/XinHao/Desktop/tES_SZ_fNIRS/'; 
all_files = dir(fullfile(Path, '5.coh_sz', '*.mat'));
% 使用GRETNA工具包计算
% Wang, J., Wang, X., Xia, M., Liao, X., Evans, A., & He, Y. (2015). GRETNA: a graph theoretical network analysis toolbox for imaging connectomics. Frontiers in human neuroscience, 9, 386.
mkdir(fullfile(Path, '7.sparse_cell')); 
for i = 1:length(all_files)
    current_file = fullfile(all_files(i).folder, all_files(i).name);
    load(current_file);
    OutputFile = fullfile(Path, '7.sparse_cell', all_files(i).name);  % 保存输出文件
    SType = 1;   % 只保留正值连接
    TType = 1;   % 稀疏度 thresholding
    Thres = linspace(0.05, 0.4, 8);   % 8 个稀疏度阈值，从 0.05 到 0.4
    NType = 1;   % 二值化网络
    gretna_RUN_ThresMat(coherence_matrix, OutputFile, SType, TType, Thres, NType);
end
group1 = {'HZ001', 'HZ002', 'HZ004', 'HZ005', 'HZ006', 'SX001'}; % Experimental
group2 = {}; % ControlActive
group3 = {'HZ007', 'SX002'}; % ControlSham
group4 = {'HZ003'}; % ControlResting
all_files = dir(fullfile(Path, '7.sparse_cell', '*.mat'));
subjects = unique(arrayfun(@(x) x.name(1:5), all_files, 'UniformOutput', false));
num_subjects = length(subjects);
task_conditions = {'R1', 'HA1', 'HB1', 'HEO1', 'HA2', 'HB2', 'HEO2', 'HA3', 'HB3', 'HEO3', 'HA4', 'HB4', 'HEO4', 'R2', 'HA5', 'HB5', 'HEO5'};
D = cell(num_subjects, 257);   
GE = cell(num_subjects, 257);

for sub = 1:num_subjects
    subject_id = subjects{sub};
    subject_files = all_files(contains({all_files.name}, subject_id));
    D{sub,1} = subject_id;
    GE{sub,1} = subject_id;
    % 分组编号
    if ismember(subject_id, group1)
        D{sub,2} = 1; GE{sub,2} = 1;
    elseif ismember(subject_id, group2)
        D{sub,2} = 2; GE{sub,2} = 2;
    elseif ismember(subject_id, group3)
        D{sub,2} = 3; GE{sub,2} = 3;
    elseif ismember(subject_id, group4)
        D{sub,2} = 4; GE{sub,2} = 4;
    else
        D{sub,2} = NaN; GE{sub,2} = NaN;
    end
    % 遍历 task_conditions，确保顺序一致
    for cond_idx = 1:length(task_conditions)
        cond_name = task_conditions{cond_idx};
        % 寻找当前被试+条件组合的文件
        target_pattern = sprintf('%s_%s.mat', subject_id, cond_name);
        matched_file = subject_files(contains({subject_files.name}, target_pattern));
        if ~isempty(matched_file)
            load(fullfile(matched_file.folder, matched_file.name));
            Dk = zeros(8, 15);
            GEk = zeros(8, 15);
            for k = 1:8
                [~, dk] = gretna_node_degree(A{k});
                [~, gEk] = gretna_node_global_efficiency(A{k});
                Dk(k,:) = dk;
                GEk(k,:) = gEk;
            end
            D(sub, (cond_idx-1)*15 + 3 : cond_idx*15 + 2) = num2cell(mean(Dk, 1));
            GE(sub, (cond_idx-1)*15 + 3 : cond_idx*15 + 2) = num2cell(mean(GEk, 1));
        else
            % 如果该条件缺失，不写入任何值，保留空单元格
            continue;
        end
    end
end

mkdir(fullfile(Path, '7.graph_indices')); 
xlswrite(fullfile(Path, '7.graph_indices', 'nodal_degree.xlsx'), D);
xlswrite(fullfile(Path, '7.graph_indices', 'nodal_efficiency.xlsx'), GE);

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