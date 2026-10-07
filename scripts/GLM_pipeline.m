clear; clc
config;

dyadFolders = dir(fullfile(analysisPath, 'Dyad*'));
fullSummary = table();
dyadInfo = table(); 

[hrf] = spm_hrf(1/fs);
tmp=[];

% Loop through dyads
for i = 1:numel(dyadFolders)
    dyadFolder = dyadFolders(i).name;
    dyadPath = fullfile(analysisPath, dyadFolder);
    skipBeh = false; % set true if this dyad can't be modelled with behaviour (missing coding / bad triggers)
    
    % load data
    cFolder = fullfile(dyadPath, ['C' dyadFolder(5:end)]); % child data
    mFolder = fullfile(dyadPath, ['M' dyadFolder(5:end)]); % mother data
    cFiles = dir(fullfile(cFolder, 'C*_ppr.mat'));
    mFiles = dir(fullfile(mFolder, 'M*_ppr.mat'));

    data1 = load(fullfile(cFolder, cFiles.name));
    data2 = load(fullfile(mFolder, mFiles.name));

    % exclude noisy channels, apply SSC reg, etc (prepareData overwrites data1/data2 in the base workspace)
    prepareData(data1,data2,true);

    % average both chromophores for all ROIs
    data1=averageROIs(data1);
    data2=averageROIs(data2);

    % Initialise boxcar stim matrix
    data1.stim=zeros(size(data1.aux,1),3);

    collabTimes = [data1.t(data1.Start_Collaboration)'; data1.t(data1.Start_CollaborationScreen)'];
    collabTimes = collabTimes(:);
    indivTimes = data1.t(data1.Start_Individual)';
  
    % Collab regressor (column 1)
    for j = 1:length(collabTimes)
        startIdx = round(collabTimes(j)*fs);
        endIdx   = round((collabTimes(j)+120)*fs);
        data1.stim(startIdx:endIdx,1) = 1;
    end

    % Individual regressor (column 2)
    for j = 1:length(indivTimes)
        startIdx = round(indivTimes(j)*fs);
        endIdx   = round((indivTimes(j)+120)*fs);
        data1.stim(startIdx:endIdx,2) = 1;
    end

    %% generate boxcars for behavioural regressors from ELAN-coded files

    % load ET data
    dyadName = data1.Pnum;
    wvFile = [dyadName '_wv.csv'];
    ex = {};

    if ~exist(fullfile(analysisPath, wvFile), 'file')
        fprintf('Skipping dyad %s: missing behavioural coding\n', dyadName);
        ex(i,:) = {dyadName, 'No beh coding'};
        skipBeh = true;
    else
    
        % wv = pupiltime --> worldview time
        wv = readtable(fullfile(analysisPath, wvFile), 'VariableNamingRule', 'preserve');
        
        % stim triggers in pupiltime
        annotations = readtable(fullfile(analysisPath, [dyadName '_ann.csv']),'VariableNamingRule', 'preserve');
        
        % ELAN coding from worldview cams
        beh = readtable(fullfile(analysisPath, [dyadName '.txt']),'VariableNamingRule', 'preserve');

        wv.Sec = (wv{:,1}) - (wv{:,1}(1));
        wv_times = wv{:, 1};             % timestamps in wv
        ann_times = annotations{:, 2};           % timestamps in ann
        ann_labels = annotations{:, 3};          % labels from ann

        % check if ET annotations seem correct
        target = ["Cooperation","CooperationScreen","Individual"];
        results = table([],[],[], 'VariableNames',{'Condition','Duration_sec','Pass'});

        for a = 1:2:length(ann_times)
            cond = ann_labels(a);
            if ismember(cond, target)
                start_time = ann_times(a);
                end_time   = ann_times(a+1);
                duration   = end_time - start_time;  % in seconds

                pass = duration >= 90; % if duration of cond is at least 90s, pass
                results = [results; {cond, duration, pass}];
            end
        end

        if isempty(results) || any(~results.Pass) % also skip if no condition triggers found
            fprintf('\nET triggers for %s wrong. Skipping.\n', dyadName);
            ex(i,:) = {dyadName, 'ET triggers wrong'};
            skipBeh = true;
        else
            fprintf('\nET triggers for %s OK.\n', dyadName);
            wv.ann = strings(height(wv), 1);

            % ET anns
            for p = 1:2:length(ann_times)-1
                start_time = ann_times(p);
                end_time = start_time+120; % truncate all conds to 120s
                label = ann_labels{p};

                % Assign label for times within [start_time, end_time)
                idx = (wv_times >= start_time) & (wv_times < end_time);
                wv.ann(idx) = label;
            end

            % align wv time with NIRS samples (anchored on start of first rest period)
            rest_idx = find(wv.ann == "rest", 1, 'first');
            eye_anchor_time = wv.Sec(rest_idx);
            nirs_anchor_idx = data1.Start_Rest(1);
            wv.NIRS_idx = round((wv.Sec - eye_anchor_time) * fs) + nirs_anchor_idx;

            % initialize behavioural reg columns in ET data
            for r = 1:length(beh_regs)
                behavior = beh_regs{r};
                wv.(behavior) = zeros(height(wv), 1);
            end

            % loop through beh rows
            for b = 1:height(beh)
                label = strtrim(beh{b,1});

                if ismember(label, beh_regs)
                    start_time = double(beh{b,3});
                    end_time = double(beh{b,4});
                    active_idx = wv.Sec >= start_time & wv.Sec <= end_time;
                    var = char(label);
                    wv.(var)(active_idx) = 1;
                end
            end

            wv.ann = strip(string(wv.ann));

            %% Check presence of each required label
            requiredLabels = ["Cooperation", "CooperationScreen", "Individual"];

            for lbl = requiredLabels
                if ~any(wv.ann == lbl)
                    warning('Missing expected annotation: %s', lbl);
                end
            end

            keepRows = wv.ann == "Cooperation" | ...
                wv.ann == "CooperationScreen";

            wv_desc = wv(keepRows, :);

            %% Durations for descriptives

            wv_desc.trial = zeros(height(wv_desc),1);

            conds = ["Cooperation", "CooperationScreen"];

            for t = conds
                rows = find(wv_desc.ann == t);
                trialCount = 1;

                wv_desc.trial(rows(1)) = trialCount;

                for z = 2:length(rows)
                    if rows(z) == rows(z-1) + 1
                        wv_desc.trial(rows(z)) = trialCount;
                    else
                        trialCount = trialCount + 1;
                        wv_desc.trial(rows(z)) = trialCount;
                    end
                end
            end

            indicatorCols = 6:15;   % beh reg cols (beh_regs)
            indicatorNames = wv_desc.Properties.VariableNames(indicatorCols);

            trials = unique(wv_desc.trial);
            anns   = unique(wv_desc.ann);

            data = [];
            idx = 1;

            for ti = 1:numel(trials)
                trialVal = trials(ti);
                rowsTrial = (wv_desc.trial == trialVal);

                for ai = 1:numel(anns)
                    annVal = anns(ai);
                    rowsAnn = (wv_desc.ann == annVal);

                    rows = find(rowsTrial & rowsAnn);
                    if isempty(rows)
                        continue;
                    end

                    for k = 1:length(indicatorCols)
                        col = indicatorCols(k);
                        vals = wv_desc{rows, col};

                        if any(vals == 1)
                            % Detect all runs of consecutive 1s
                            % pad with 0 at start and end to catch edges
                            padded = [0; vals; 0];
                            d = diff(padded);
                            runStarts = find(d == 1);  % indices where 1-run starts (relative to rows)
                            runEnds   = find(d == -1) - 1;  % indices where 1-run ends

                            for r = 1:length(runStarts)
                                firstIdx = rows(runStarts(r));
                                lastIdx  = rows(runEnds(r));
                                firstSec = wv_desc.Sec(firstIdx);
                                lastSec  = wv_desc.Sec(lastIdx);

                                data(idx).trial    = trialVal;
                                data(idx).ann      = annVal;
                                data(idx).measure  = indicatorNames{k};
                                data(idx).firstSec = firstSec;
                                data(idx).lastSec  = lastSec;
                                data(idx).id       = dyadName;
                                idx = idx + 1;
                            end
                        end
                    end
                end
            end
            resultsTab = struct2table(data);   % convert struct array to table
            resultsTab.duration = resultsTab.lastSec - resultsTab.firstSec;
            fullSummary = [fullSummary; resultsTab];

            % number of collaboration trials (both collab conditions) in the fixed 120 s ET windows
            nCollabTrials = height(unique(wv_desc(:, {'ann','trial'}), 'rows'));
            dyadInfo = [dyadInfo; table({dyadName}, nCollabTrials, ...
                {data1.Start_Collaboration(:)'}, {data1.Start_CollaborationScreen(:)'}, {data1.Start_Individual(:)'}, ...
                'VariableNames', {'id','nCollabTrials','onsetCollab','onsetScreen','onsetInd'})]; % trial onsets (NIRS samples) for trial order

            % collapse regs & collab conditions
            cnsCols = contains(wv.Properties.VariableNames, "cns");
            cwsCols = contains(wv.Properties.VariableNames, "cws");

            collabCols = cnsCols | cwsCols;
            wv.scaffolding = any(wv{:, collabCols} == 1, 2);

            wv = wv(:, [wv.Properties.VariableNames(1:5), {'scaffolding'}]);            

            col = 6; % scaff column
            duration = wv{:,col};
            edges = diff([0; duration; 0]);
            start_indices = find(edges == 1);
            end_indices = find(edges == -1) - 1;

            for j = 1:length(start_indices)
                start_i = start_indices(j);
                end_i = end_indices(j);

                start_idx = wv.NIRS_idx(start_i);
                end_idx = wv.NIRS_idx(end_i);

                stim_col = col - 3; % scaffolding (wv col 6) goes into stim col 3

                % mask the NIRS samples covered by this scaffolding bout
                mask = false(size(data1.stim,1),1); % create vector same num rows as data1.stim w/ all initially false
                if start_idx <= length(mask) && end_idx >= 1
                    start_idx = max(start_idx, 1);
                    end_idx = min(end_idx, length(mask));
                    mask(start_idx:end_idx) = true;
                end

                validMask = (data1.stim(:,1) == 1); % collab column = both collab conditions (with/without screen); scaffolding outside these 120 s windows is ignored
                finalMask = mask & validMask;

                % Assign 1s only where finalMask is true
                if any(finalMask)
                    data1.stim(finalMask, stim_col) = 1;
                end
            end
        end
    end

    fields = fieldnames(data1.roiMeans);
    
    %% BEHAVIOUR ONLY MODELS %%
    % skipped dyads get NaN rows in the same layout as modelled dyads
    if skipBeh
        fprintf('Skipping behaviour model for %s.\n', dyadName);
        for type = {'hbo','hbr'}
            for r = 1:numROIs
                roi1 = fields{r};
                tmp.beh.(type{1}).R2.(roi1)(i,:)    = NaN;
                tmp.beh.(type{1}).R2adj.(roi1)(i,:) = NaN;
                tmp.beh.(type{1}).Betas.(roi1)(i,:) = nan(1,4);
            end
        end
    else
        % convolve stim matrix (condition & behav regs only)
        X.dm=[];
        for r=1:size(data1.stim,2)
            temp=conv(hrf,data1.stim(:,r));
            temp=temp(1:length(data1.t),:);
            X.dm(:,r)=temp;
        end
        X.dm = [X.dm ones(size(X.dm,1),1)]; % add constant (as in brain & behBrain models)

        for type = {'hbo','hbr'}
            for r = 1:numROIs
                roi1 = fields{r};
                if ~all(isnan(data1.roiMeans.(roi1).(type{1})(:,1))) % skip ROIs with no usable channels

                    DesMat=X.dm; %Design Matrix (same for every ROI)
                    NirsSignal = data1.roiMeans.(roi1).(type{1})(:,1); % This is the child signal Y in Y=Xbeta+e
                    colsWithNaN = any(isnan(DesMat), 1); % Identify columns with any NaN

                    betaVec = nan(1, 4); % 4 regs (2 conds, scaffolding, constant)
                    validCols = ~colsWithNaN & any(DesMat ~= 0, 1); % columns without NaNs and not all-zero (behaviour never shown -> NaN beta, not 0)
                    % compute regression only on valid columns
                    betaValid = DesMat(:,validCols) \ NirsSignal;
                    betaVec(validCols) = betaValid';

                    % --- Model comparison metrics ---

                    % Compute fitted values
                    yhat = DesMat(:,validCols) * betaValid;

                    % Compute residuals
                    residuals = NirsSignal - yhat;

                    RSS = sum(residuals.^2);
                    TSS = sum((NirsSignal - mean(NirsSignal)).^2);
                    k   = sum(validCols);     % number of effective predictors
                    N   = length(NirsSignal); % number of time points

                    % R-squared
                    R2 = 1 - RSS/TSS;

                    % Adjusted R-squared
                    R2_adj = 1 - (RSS/(N-k)) / (TSS/(N-1));

                    % Store metrics
                    tmp.beh.(type{1}).R2.(roi1)(i,:)     = R2;
                    tmp.beh.(type{1}).R2adj.(roi1)(i,:)  = R2_adj;
                    % -----------------------------------------

                    tmp.beh.(type{1}).Betas.(roi1)(i,:) = betaVec';
                else
                    tmp.beh.(type{1}).Betas.(roi1)(i,:)=nan(1,4); % if child data not available, NaN everything
                    tmp.beh.(type{1}).R2.(roi1)(i,:)    = NaN;
                    tmp.beh.(type{1}).R2adj.(roi1)(i,:) = NaN;
                end
            end
        end
    end

    %% BRAIN ONLY MODELS %%%

    % convolve stim matrix (condition regs only (first 2 cols))
    X.dm=[];
    for r=1:size(data1.stim(:, 1:2),2)
        temp=conv(hrf,data1.stim(:,r));
        temp=temp(1:length(data1.t),:);
        X.dm(:,r)=temp;
    end

    for m = 1:numROIs
        roi1 = fields{m};
        
        X.hbo.(roi1)=X.dm;
        X.hbr.(roi1)=X.dm;

        for type = {'hbo','hbr'} % process for both haem types
            for j = 1:numROIs % for all ROIs
                roi2 = fields{j};
                nirsParent = data2.roiMeans.(roi2).(type{1})(:,1); % get sig for each ROI of mum (all NaN if ROI not valid -> col dropped before fitting)
                for k = 1:2 % for collab vs ind conditions
                    nirs3=nirsParent;
                    nirs3(data1.stim(:,k) ~= 1, 1) = 0; % keep signal from condition
                    X.(type{1}).(roi1) = [X.(type{1}).(roi1) nirs3]; % insert in dm (cols=collab/ind)
                end
            end
            X.(type{1}).(roi1)=[X.(type{1}).(roi1) ones(size(X.(type{1}).(roi1),1),1)]; % add constant
        end
    end

    % calculate betas
    for type = {'hbo','hbr'}
        for r = 1:numROIs
            roi1 = fields{r};
            if ~all(isnan(data1.roiMeans.(roi1).(type{1})(:,1))) % skip ROIs with no usable channels

                DesMat=X.(type{1}).(roi1); %Design Matrix
                NirsSignal = data1.roiMeans.(roi1).(type{1})(:,1); % This is the child signal Y in Y=Xbeta+e
                colsWithNaN = any(isnan(DesMat), 1); % Identify columns with any NaN

                betaVec = nan(1, 11); %11 regs (2 conds, 8 partner ROI x cond, constant)
                validCols = ~colsWithNaN & any(DesMat ~= 0, 1); % columns without NaNs and not all-zero
                % compute regression only on valid columns
                betaValid = DesMat(:,validCols) \ NirsSignal;
                betaVec(validCols) = betaValid';

                % --- Model comparison metrics ---
                
                % Compute fitted values
                yhat = DesMat(:,validCols) * betaValid;

                % Compute residuals
                residuals = NirsSignal - yhat;

                RSS = sum(residuals.^2);
                TSS = sum((NirsSignal - mean(NirsSignal)).^2);
                k   = sum(validCols);     % number of effective predictors
                N   = length(NirsSignal); % number of time points

                % R-squared
                R2 = 1 - RSS/TSS;

                % Adjusted R-squared
                R2_adj = 1 - (RSS/(N-k)) / (TSS/(N-1));

                % Store metrics
                tmp.brain.(type{1}).R2.(roi1)(i,:)     = R2;
                tmp.brain.(type{1}).R2adj.(roi1)(i,:)  = R2_adj;
                % -----------------------------------------

                tmp.brain.(type{1}).Betas.(roi1)(i,:) = betaVec';
            else
                tmp.brain.(type{1}).Betas.(roi1)(i,:)=nan(1,11); % if child data not available, NaN everything
                tmp.brain.(type{1}).R2.(roi1)(i,:)    = NaN;
                tmp.brain.(type{1}).R2adj.(roi1)(i,:) = NaN;
            end
        end
    end

    %% BEHAVIOUR & BRAIN MODELS %%
    % skipped dyads get NaN rows in the same layout as modelled dyads
    if skipBeh
        fprintf('Skipping brain + behaviour model for %s.\n', dyadName);
        for type = {'hbo','hbr'}
            for r = 1:numROIs
                roi1 = fields{r};
                tmp.behBrain.(type{1}).R2.(roi1)(i,:)    = NaN;
                tmp.behBrain.(type{1}).R2adj.(roi1)(i,:) = NaN;
                tmp.behBrain.(type{1}).Betas.(roi1)(i,:) = nan(1,12);
            end
        end
    else
        % convolve stim matrix (condition & behav regs)
        X.dm=[];
        for r=1:size(data1.stim,2)
            temp=conv(hrf,data1.stim(:,r));
            temp=temp(1:length(data1.t),:);
            X.dm(:,r)=temp;
        end

        for m = 1:numROIs
            roi1 = fields{m};

            X.hbo.(roi1)=X.dm;
            X.hbr.(roi1)=X.dm;

            for type = {'hbo','hbr'} % process for both haem types
                for j = 1:numROIs % for all ROIs
                    roi2 = fields{j};
                    nirsParent = data2.roiMeans.(roi2).(type{1})(:,1); % get sig for each ROI of mum (all NaN if ROI not valid -> col dropped before fitting)
                    for k = 1:2 % for each collab & ind cond
                        nirs3=nirsParent;
                        nirs3(data1.stim(:,k) ~= 1, 1) = 0; % keep signal from condition
                        X.(type{1}).(roi1) = [X.(type{1}).(roi1) nirs3]; % insert in dm (cols=collab/ind)
                    end
                end
                X.(type{1}).(roi1)=[X.(type{1}).(roi1) ones(size(X.(type{1}).(roi1),1),1)]; % add constant
            end
        end

        % calculate child-based betas
        for type = {'hbo','hbr'}
            for r = 1:numROIs
                roi1 = fields{r};
                if ~all(isnan(data1.roiMeans.(roi1).(type{1})(:,1))) % skip ROIs with no usable channels

                    DesMat=X.(type{1}).(roi1); %Design Matrix
                    NirsSignal = data1.roiMeans.(roi1).(type{1})(:,1); % This is the child signal Y in Y=Xbeta+e
                    colsWithNaN = any(isnan(DesMat), 1); % Identify columns with any NaN

                    betaVec = nan(1, 12); % 12 regs, 2 conds + 1 beh + 8 partner ROI x cond + constant
                    validCols = ~colsWithNaN & any(DesMat ~= 0, 1); % columns without NaNs and not all-zero (behaviour never shown -> NaN beta, not 0)
                    % compute regression only on valid columns
                    betaValid = DesMat(:,validCols) \ NirsSignal;
                    betaVec(validCols) = betaValid';

                    % --- Model comparison metrics ---

                    % Compute fitted values
                    yhat = DesMat(:,validCols) * betaValid;

                    % Compute residuals
                    residuals = NirsSignal - yhat;

                    RSS = sum(residuals.^2);
                    TSS = sum((NirsSignal - mean(NirsSignal)).^2);
                    k   = sum(validCols);     % number of effective predictors
                    N   = length(NirsSignal); % number of time points

                    % R-squared
                    R2 = 1 - RSS/TSS;

                    % Adjusted R-squared
                    R2_adj = 1 - (RSS/(N-k)) / (TSS/(N-1));

                    % Store metrics
                    tmp.behBrain.(type{1}).R2.(roi1)(i,:)     = R2;
                    tmp.behBrain.(type{1}).R2adj.(roi1)(i,:)  = R2_adj;
                    % -----------------------------------------

                    tmp.behBrain.(type{1}).Betas.(roi1)(i,:) = betaVec';
                else
                    tmp.behBrain.(type{1}).Betas.(roi1)(i,:)=nan(1,12); % if child data not available, NaN everything
                    tmp.behBrain.(type{1}).R2.(roi1)(i,:)    = NaN;
                    tmp.behBrain.(type{1}).R2adj.(roi1)(i,:) = NaN;
                end
            end
        end
    end
    tmp.dyadNames{i} = dyadFolder; % dyadnumbers
end