%% 近红外研究的行为结果用这一版本，非常且严谨稳健

% 将每个被试的有效.mat文件分实验组放置在同个文件夹中，包括Block1~后测Block5
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

                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                % 若只考虑8个trials
                % keep = false(length(trials), 1);  % 遍历查找
                % for i = 1:length(trials)
                %     fname = trials(i).name;
                %     parts = split(fname, '_');
                %     if length(parts) >= 4
                %         thirdPart = str2double(parts{3});
                %         if ~isnan(thirdPart) && thirdPart >= 2 && thirdPart <= 9
                %             keep(i) = true;  % 标记保留
                %         end
                %     end
                % end
                % trials = trials(keep);  % 不考虑trial 1/10的文件
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%

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
                    % 这里为了将设备时间统一到零点，先计算一遍RT
                    A.IOI = diff(A.RT);
                    B.IOI = diff(B.RT);
                    A.Trials(:, idx) = {A.IOI};
                    B.Trials(:, idx) = {B.IOI};
                    A.RT=cumsum(A.IOI);
                    B.RT=cumsum(B.IOI);

                    %% 时间序列预处理和对齐

                    % 未按指导语进行实验，前八拍未跟上
                    if A.RT(1) > 4000 || B.RT(1) > 4000
                        display("RT(1) too large, skip trial, idx = " + idx);
                        skipped = skipped + 1;
                        continue;  % 跳过当前循环，进入下一个trial
                    end

                    % 剔除初始RT小于500ms的抢按（这些点不应该参与IOI计算）
                    if A.RT(1) < 500 || B.RT(1) < 500
                        A.RT(A.RT < 250) = [];  
                        B.RT(B.RT < 250) = [];    
                    end

                    if numel(A.RT) <= numel(B.RT)  % A短时，将B从距离A起始最小的点截等长
                        [~, startIdx] = min(abs(B.RT - A.RT(1)));
                        if startIdx + numel(A.RT) - 1 <= numel(B.RT)
                            B.RT = B.RT(startIdx : startIdx + numel(A.RT) - 1);
                        else
                            len = numel(B.RT) - startIdx + 1;
                            B.RT = B.RT(startIdx : end);
                            A.RT = A.RT(1:len);
                        end
                    else  % B短时，将A从距离B起始最小的点截等长
                        [~, startIdx] = min(abs(A.RT - B.RT(1)));
                        if startIdx + numel(B.RT) - 1 <= numel(A.RT)
                            A.RT = A.RT(startIdx : startIdx + numel(B.RT) - 1);
                        else % 如果A剩余长度不够，再截短B使两者长度相等
                            len = numel(A.RT) - startIdx + 1;
                            A.RT = A.RT(startIdx : end);
                            B.RT = B.RT(1:len);
                        end
                    end

                    % AB一定等长，统一计算 IOI
                    A.IOI = diff(A.RT);
                    B.IOI = diff(B.RT);

                    pair.Aoutlier(idx,1) = sum(nonzeros(A.IOI <= median(A.IOI) - 0.5 * median(A.IOI)));
                    pair.Aoutlier(idx,2) = sum(nonzeros(A.IOI >= median(A.IOI) + 0.5 * median(A.IOI)));
                    pair.Boutlier(idx,1) = sum(nonzeros(B.IOI <= median(B.IOI) - 0.5 * median(B.IOI)));
                    pair.Boutlier(idx,2) = sum(nonzeros(B.IOI >= median(B.IOI) + 0.5 * median(B.IOI)));
    
                    %% 计算四种行为指标：IC WS RD TV，从前8拍结束后的点算起，不一定是24个点
                    [~, closestIdx] = min(abs(A.RT - 4000));
                    % 保证后面还有数据点可用，至少两个RT的点：(closestIdx+1:end)
                    if closestIdx+2 >= length(A.RT)
                        warning("A.RT 中 4000ms之后几乎没有有效数据，跳过这个 trial");
                        skipped = skipped + 1;
                        continue;
                    end
                    % index1: Within-subject Tapping Stability
                    % For B，公式：sqrt(1/std(IOI))
                    pair.WS(idx,:) = sqrt(1 / std(B.IOI(closestIdx+1:end)));
                    % index2: Interpersonal Consistency (Synchronization Index)，对齐至标准节奏
                    phaseA = t2phases(cumsum(500 * ones(numel(A.RT(closestIdx+1:end)), 1)), A.RT(closestIdx+1:end));
                    phaseB = t2phases(cumsum(500 * ones(numel(B.RT(closestIdx+1:end)), 1)), B.RT(closestIdx+1:end));
                    phase2keep = min(length(phaseA), length(phaseB));
                    dtheta = phaseA(1:phase2keep) - phaseB(1:phase2keep);
                    pair.IC(idx,:) = abs(mean(exp(1i * dtheta)));
                    % index3: Rhythm Deviation
                    % For B，计算实际 IOI 与目标500ms之间的绝对差的中位数
                    pair.RD(idx,:) = median(abs(B.IOI(closestIdx+1:end) - 500));
                    % index4: Tapping Variability
                    % 以 A 为参考，计算每次敲击时刻差（asynchrony）的绝对值，再计算这些值与其中位数之间差值的中位数（MAD）
                    asynchrony = abs(B.RT(closestIdx+1:end) - A.RT(closestIdx+1:end));
                    pair.TV(idx,:) = median(abs(asynchrony - median(asynchrony)));
                end % trials 遍历结束

            % 剔除三个MAD的异常值后求平均数
            constant = 1.4826; % 常数 constant （正态分布调整系数）
            % summaryWS{idx_sub, idx_bc} = mean(pair.WS(abs(pair.WS - median(pair.WS, 'omitnan')) <= 3 * constant * median(abs(pair.WS - median(pair.WS, 'omitnan')), 'omitnan')), 'omitnan'); 
            % summaryIC{idx_sub, idx_bc} = mean(pair.IC(abs(pair.IC - median(pair.IC, 'omitnan')) <= 3 * constant * median(abs(pair.IC - median(pair.IC, 'omitnan')), 'omitnan')), 'omitnan'); 
            % summaryRD{idx_sub, idx_bc} = mean(pair.RD(abs(pair.RD - median(pair.RD, 'omitnan')) <= 3 * constant * median(abs(pair.RD - median(pair.RD, 'omitnan')), 'omitnan')), 'omitnan');
            % summaryTV{idx_sub, idx_bc} = mean(pair.TV(abs(pair.TV - median(pair.TV, 'omitnan')) <= 3 * constant * median(abs(pair.TV - median(pair.TV, 'omitnan')), 'omitnan')), 'omitnan');
            % 求中位数，异常值不敏感
            if median(pair.WS, 'omitnan') ~= 0
                summaryWS{idx_sub, idx_bc} = median(pair.WS, 'omitnan'); 
            end
            if median(pair.IC, 'omitnan') ~= 0
                summaryIC{idx_sub, idx_bc} = median(pair.IC, 'omitnan'); 
            end
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

% xlswrite("C:/Users/XinHao/Desktop/tES_SZ_Behav/behavior_data.xlsx", SummaryIC, 'IC');
% xlswrite("C:/Users/XinHao/Desktop/tES_SZ_Behav/behavior_data.xlsx", SummaryWS, 'WS');
% xlswrite("C:/Users/XinHao/Desktop/tES_SZ_Behav/behavior_data.xlsx", SummaryRD, 'RD');
% xlswrite("C:/Users/XinHao/Desktop/tES_SZ_Behav/behavior_data.xlsx", SummaryTV, 'TV');
writecell(SummaryIC, "C:/Users/XinHao/Desktop/tES_SZ_Behav/behavior_data.xlsx", 'Sheet', 'IC');
writecell(SummaryWS, "C:/Users/XinHao/Desktop/tES_SZ_Behav/behavior_data.xlsx", 'Sheet', 'WS');
writecell(SummaryRD, "C:/Users/XinHao/Desktop/tES_SZ_Behav/behavior_data.xlsx", 'Sheet', 'RD');
writecell(SummaryTV, "C:/Users/XinHao/Desktop/tES_SZ_Behav/behavior_data.xlsx", 'Sheet', 'TV');

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