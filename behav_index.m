%% 将每个被试的有效.mat文件分实验组放置在同个文件夹中，包括Block1~后测Block5
% 每组的文件夹以被试五位数ID命名
clc; clear;
basepath = "C:/Users/XinHao/Desktop/tES_SZ_Behav";
subfolders = dir(basepath);
subfolders = subfolders([subfolders.isdir]);  
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));  
groupMap = containers.Map({'Experimental', 'ControlActive', 'ControlSham', 'ControlResting', 'ControlBehavior'}, {1, 2, 3, 4, 5});
column_names = {'ID', 'Group'}; 
for block = 1:5
    for condition = 1:3
        column_names{end+1} = sprintf('B%dC%d', block, condition);
    end
end
SummaryWS = {}; 
SummaryIC = {}; 
SummaryRD = {}; 
SummaryTV = {};
SummaryWS(1, :) = column_names;
SummaryIC(1, :) = column_names;
SummaryRD(1, :) = column_names;
SummaryTV(1, :) = column_names;

for i = 1:numel(subfolders)
    subfolderName = subfolders(i).name;
    subfolderPath = fullfile(basepath, subfolderName);

    subjects = dir(fullfile(subfolderPath, '*')); 
    subjects = subjects([subjects.isdir]);  
    subjects = subjects(~ismember({subjects.name}, {'.', '..'})); 
    subjects = subjects(~cellfun(@isempty, regexp({subjects.name}, '^\d{5}$'))); 
    groupVal = groupMap(subfolderName);

    %% 初始化数据存储数组
    summaryTypes = {'summaryWS', 'summaryIC', 'summaryRD', 'summaryTV'};
    for type = 1:numel(summaryTypes)
        eval([summaryTypes{type}, ' = cell(numel(subjects) + 1, numel(column_names));']);
        eval([summaryTypes{type}, '(1, :) = column_names;']);
    end
    idx_sub = 2;  % 从第二行开始存入数据

    %% 循环计算
    for subject = 1:numel(subjects) % 遍历所有subjects
        idx_bc = 3;  % 从第3列开始存入数据
        subjectID = subjects(subject).name; 
        summaryWS{idx_sub, 1} = subjectID;
        summaryIC{idx_sub, 1} = subjectID;
        summaryRD{idx_sub, 1} = subjectID;
        summaryTV{idx_sub, 1} = subjectID;
        summaryWS{idx_sub, 2} = groupVal;
        summaryIC{idx_sub, 2} = groupVal;
        summaryRD{idx_sub, 2} = groupVal;
        summaryTV{idx_sub, 2} = groupVal;
        filepath = fullfile(subfolderPath, subjectID); 
        cd(filepath);
        tiqu = @(x) split(x, ":");
        shuju = @(x) str2double(x{end,1});
        pair = table;
        idx = 0;

        for block = 1:5 % 遍历所有blocks
            for condition = 1:3 % 遍历所有conditions
                pattern = sprintf('%d_%d_*_*.mat', block, condition);
                trials = dir(fullfile(filepath, pattern));
                % 跳过文件缺失的block-condition配对
                if isempty(trials)
                    summaryWS{idx_sub, idx_bc} = NaN;
                    summaryIC{idx_sub, idx_bc} = NaN;
                    summaryRD{idx_sub, idx_bc} = NaN;
                    summaryTV{idx_sub, idx_bc} = NaN;
                    idx_bc = idx_bc + 1; 
                    continue; 
                end

                skipped = 0;
                for n = 1:2:numel(trials) % 遍历所有trials
                    idx = (n-1)/2 + 1;
                    k = strfind(trials(n).name, '_');
                    prefixA = trials(n).name(1:k(3));
                    prefixB = trials(n+1).name(1:k(3));

                    if ~strcmp(prefixA, prefixB)
                        display(trials(n).name + " missing its pair, now idx = " + idx);
                        skipped = skipped + 1;
                        continue;
                    end
    
                    trialInfo = trials(n).name(1:k(2)); 
                    pair.block(idx) = str2double(trialInfo(1));
                    pair.condition(idx) = str2double(trialInfo(3:end-1));
                    b = trials(n+1).name(1:k(3));
                    pair.trial(idx) = str2double(b(5:end-1));

                    load(trials(n+1).name); 
                    y = x;  % 将B的数据重命名为 y
                    load(trials(n).name);  

                    if isempty(x) || isempty(y) || numel(x) < 8 || numel(y) < 8
                        display("empty mat, now idx = " + idx);
                        skipped = skipped + 1;
                        continue;
                    end

                    y(end) = [];
                    y = y(~cellfun(@isempty, y));
                    y = cellfun(@string, y, 'UniformOutput', false);
                    x(end) = [];
                    x = x(~cellfun(@isempty, x));
                    x = cellfun(@string, x, 'UniformOutput', false);
                    y = cellfun(tiqu, y, 'UniformOutput', false);
                    x = cellfun(tiqu, x, 'UniformOutput', false);
    
                    B.RT = cellfun(shuju, y);
                    A.RT = cellfun(shuju, x);
                    A.RT(isnan(A.RT)) = [];
                    B.RT(isnan(B.RT)) = [];
                    B.IOI = diff(B.RT);
                    A.IOI = diff(A.RT);
                    B.Trials(:, idx) = {B.IOI};
                    A.Trials(:, idx) = {A.IOI};
                    B.RT = cumsum(B.IOI);
                    A.RT = cumsum(A.IOI);
    
                    %% 时间序列预处理和对齐
                    if A.RT(1) < 250 || B.RT(1) < 250
                        A.IOI(A.RT < 250) = [];
                        B.IOI(B.RT < 250) = [];
                        A.RT = cumsum(A.IOI);
                        B.RT = cumsum(B.IOI);
                    end

                    if A.RT(1) > 2000 || B.RT(1) > 2000
                        pair.deletion8(idx) = true;
                    else
                        pair.deletion8(idx) = false;
                    end

                    if numel(B.RT) < numel(A.RT)
                        A.IOI(numel(B.RT)+1:end) = [];
                        A.RT(numel(B.RT)+1:end) = [];
                    else
                        B.IOI(numel(A.RT)+1:end) = [];
                        B.RT(numel(A.RT)+1:end) = [];
                    end

                    points = numel(A.RT);
                    pair.Aoutlier(idx,1) = sum(nonzeros(A.IOI <= median(A.IOI) - 0.5 * median(A.IOI)));
                    pair.Aoutlier(idx,2) = sum(nonzeros(A.IOI >= median(A.IOI) + 0.5 * median(A.IOI)));
                    pair.Boutlier(idx,1) = sum(nonzeros(B.IOI <= median(B.IOI) - 0.5 * median(B.IOI)));
                    pair.Boutlier(idx,2) = sum(nonzeros(B.IOI >= median(B.IOI) + 0.5 * median(B.IOI)));
    
                    %% 计算四种行为指标：IC WS RD TV，只有ab的数据都正确记录才会计算
                    % index1: Within-subject Tapping Stability
                    % For B，公式：sqrt(1/std(IOI))
                    pair.WS(idx,:) = sqrt(1 / std(B.IOI));
                    % index2: Interpersonal Consistency (Synchronization Index)
                    phaseA = t2phases(cumsum(500 * ones(points, 1)), A.RT);
                    phaseB = t2phases(cumsum(500 * ones(points, 1)), B.RT);
                    phase2keep = min(length(phaseA), length(phaseB));
                    dtheta = phaseA(1:phase2keep) - phaseB(1:phase2keep);
                    pair.IC(idx,:) = abs(mean(exp(1i * dtheta)));
                    % index3: Rhythm Deviation
                    % For B，计算实际 IOI 与目标500ms之间的绝对差的中位数
                    pair.RD(idx,:) = median(abs(B.IOI - 500));
                    % index4: Tapping Variability
                    % 以 A 为参考，计算每次敲击时刻差（asynchrony）的绝对值，再计算这些值与其中位数之间差值的中位数（MAD）
                    asynchrony = abs(B.RT - A.RT);
                    pair.TV(idx,:) = median(abs(asynchrony - median(asynchrony)));
                end % trials 遍历结束

            % 剔除三个MAD的异常值后求平均数
            %summaryWS{idx_sub, idx_bc} = mean(pair.WS(abs(pair.WS - median(pair.WS, 'omitnan')) <= 3 * median(abs(pair.WS - median(pair.WS, 'omitnan')), 'omitnan')), 'omitnan'); 
            %summaryIC{idx_sub, idx_bc} = mean(pair.IC(abs(pair.IC - median(pair.IC, 'omitnan')) <= 3 * median(abs(pair.IC - median(pair.IC, 'omitnan')), 'omitnan')), 'omitnan'); 
            %summaryRD{idx_sub, idx_bc} = mean(pair.RD(abs(pair.RD - median(pair.RD, 'omitnan')) <= 3 * median(abs(pair.RD - median(pair.RD, 'omitnan')), 'omitnan')), 'omitnan');
            %summaryTV{idx_sub, idx_bc} = mean(pair.TV(abs(pair.TV - median(pair.TV, 'omitnan')) <= 3 * median(abs(pair.TV - median(pair.TV, 'omitnan')), 'omitnan')), 'omitnan');
            % 求中位数，异常值不敏感
            summaryWS{idx_sub, idx_bc} = median(pair.WS, 'omitnan'); 
            summaryIC{idx_sub, idx_bc} = median(pair.IC, 'omitnan'); 
            summaryRD{idx_sub, idx_bc} = median(pair.RD, 'omitnan');
            summaryTV{idx_sub, idx_bc} = median(pair.TV, 'omitnan');
            idx_bc = idx_bc + 1; % 更新列索引

            end % condition 遍历结束
        end % block 遍历结束

        idx_sub = idx_sub + 1; % 更新行索引
    end %subject遍历结束
    SummaryWS = [SummaryWS; summaryWS(2:end, :)];
    SummaryIC = [SummaryIC; summaryIC(2:end, :)];
    SummaryRD = [SummaryRD; summaryRD(2:end, :)];
    SummaryTV = [SummaryTV; summaryTV(2:end, :)];
end %subfolder遍历结束

xlswrite("C:/Users/XinHao/Desktop/tES_SZ_Behav/behavior_data.xlsx", SummaryIC, 'IC');
xlswrite("C:/Users/XinHao/Desktop/tES_SZ_Behav/behavior_data.xlsx", SummaryWS, 'WS');
xlswrite("C:/Users/XinHao/Desktop/tES_SZ_Behav/behavior_data.xlsx", SummaryRD, 'RD');
xlswrite("C:/Users/XinHao/Desktop/tES_SZ_Behav/behavior_data.xlsx", SummaryTV, 'TV');
%writecell(SummaryIC, "C:/Users/XinHao/Desktop/tES_SZ_Behav/behavior_data.xlsx", 'Sheet', 'IC');
%writecell(SummaryWS, "C:/Users/XinHao/Desktop/tES_SZ_Behav/behavior_data.xlsx", 'Sheet', 'WS');
%writecell(SummaryRD, "C:/Users/XinHao/Desktop/tES_SZ_Behav/behavior_data.xlsx", 'Sheet', 'RD');
%writecell(SummaryTV, "C:/Users/XinHao/Desktop/tES_SZ_Behav/behavior_data.xlsx", 'Sheet', 'TV');

%% 自定义时序转换 t2phases 函数
function theta2 = t2phases(metronome, T)
    jstart = find(metronome - T(1) >= 0, 1);
    j = jstart;
    n = 1;
    theta = nan(size(metronome));  
    while j <= numel(metronome)
        % 找出 T 中比当前 metronome(j) 小的最大索引 n
        idx = find(sign(metronome(j) - T) == -1, 1) - 1;
        % 如果没有找到（即 metronome(j) < 所有 T），或 n 超范围，跳过
        if isempty(idx) || idx < 1 || idx >= length(T)
            j = j + 1;
            continue;
        end
        % 正常计算 phase
        n = idx;
        theta(j) = (metronome(j) - T(n)) / (T(n+1) - T(n)) * 2*pi + 2*pi*n;
        j = j + 1;
    end
    theta2 = theta(~isnan(theta));  
end