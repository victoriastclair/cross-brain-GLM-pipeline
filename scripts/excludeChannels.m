%% Exclude Noisy Channels
%{

Written by Victoria St. Clair, Centre for Brain and Cognitive Development
Email: v.stclair@bbk.ac.uk

Inputs: 
    data = NIRS data
    varName = variable name for data

%}

function excludeChannels(data, varName)

config;

filename = fullfile(analysisPath, 'ChannelExclusions.xlsx');

if contains(data.Pnum,'C')
    sheet = 'data1_exclusions';
else
    sheet = 'data2_exclusions';
end

channels = readtable(filename, 'Sheet', sheet,'ReadVariableNames',true);
idx = find(strcmp(channels.Properties.VariableNames,data.Pnum));
if isempty(idx)
    error('%s not found in %s (sheet %s)', data.Pnum, filename, sheet);
end
excludedChannels = find(channels{:,(idx)}==0);

SSlistGood = setdiff(sscs, excludedChannels);

data.allExcludedChannels = sortrows(unique(vertcat(excludedChannels, sscs)));
data.SSlistGood = SSlistGood;

assignin('caller', varName, data);

end