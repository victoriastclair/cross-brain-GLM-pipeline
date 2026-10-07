% Average Signal Within-ROI (xGLM version)
%{
Averages available long-separation channels within each ROI, separately
for HbO and HbR. Excluded channels (from excludeChannels) are dropped.
ROIs with fewer than 2 usable channels are set to NaN.

Input:  data = participant struct after prepareData (needs AllChannelsHbo,
               AllChannelsHbr, allExcludedChannels)
Output: data with data.roiMeans.(roi).hbo / .hbr added (samples x 1)
%}

function data = averageROIs(data)

config;

roiChans = {leftPFC, rightPFC, leftTPJ, rightTPJ}; % same order as roiNames
signals  = struct('hbo', data.AllChannelsHbo, 'hbr', data.AllChannelsHbr);
types    = {'hbo', 'hbr'};

for r = 1:numel(roiNames)
    chans = setdiff(roiChans{r}, data.allExcludedChannels);

    for t = 1:numel(types)
        sig = signals.(types{t})(:, chans);
        goodChans = any(~isnan(sig), 1);

        if sum(goodChans) < 2
            data.roiMeans.(roiNames{r}).(types{t}) = NaN(size(sig, 1), 1);
        else
            data.roiMeans.(roiNames{r}).(types{t}) = mean(sig(:, goodChans), 2, 'omitnan');
        end
    end
end