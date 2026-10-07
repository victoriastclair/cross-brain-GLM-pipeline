% Prepare data for xGLM
function prepareData(data1, data2, ssc_reg)

excludeChannels(data1, 'data1');
excludeChannels(data2, 'data2');

if isempty(data1.SSlistGood) && isempty(data2.SSlistGood)
    ssc_reg = false;
end

if ssc_reg == true
    dc_orig1 = data1.procResult.dc(:,1:2,:);
    dc1 = permute(dc_orig1, [1, 3, 2]); 

    dc_orig2 = data2.procResult.dc(:,1:2,:);
    dc2 = permute(dc_orig2, [1, 3, 2]);

    if isempty(data1.SSlistGood) && isempty(data2.SSlistGood)
        fprintf('SSC reg is true but neither participant has a valid channel. Continuing with raw data. \n');
        ssc_reg=false;
    elseif isempty(data1.SSlistGood) && ~isempty(data2.SSlistGood)
        fprintf('SSC reg is true, but only data 2 has valid SSC. Running reg on data 2 only. \n');
        filtered_dc1 = dc1;
        [filtered_dc2] = PhysiologyRegression_GLM_fnirs_course(dc2, data2.SD, data2.SSlistGood,[]);
    elseif ~isempty(data1.SSlistGood) && isempty(data2.SSlistGood)
        fprintf('SSC reg is true, but only data 1 has valid SSC. Running reg on data 1 only. \n');
        [filtered_dc1] = PhysiologyRegression_GLM_fnirs_course(dc1, data1.SD, data1.SSlistGood,[]);
        filtered_dc2 = dc2;
    elseif ~isempty(data1.SSlistGood) && ~isempty(data2.SSlistGood)
        fprintf('SSC reg is true and both participants have at least one valid SSC. Running reg on both data 1 and data 2. \n');
        [filtered_dc1] = PhysiologyRegression_GLM_fnirs_course(dc1, data1.SD, data1.SSlistGood,[]);
        [filtered_dc2] = PhysiologyRegression_GLM_fnirs_course(dc2, data2.SD, data2.SSlistGood,[]);
    end
    
    data1.ssc_reg=ssc_reg;
    data2.ssc_reg=ssc_reg;

    data1.AllChannelsHbo = squeeze(filtered_dc1(:,:,1));
    data2.AllChannelsHbo = squeeze(filtered_dc2(:,:,1));

    data1.AllChannelsHbr = squeeze(filtered_dc1(:,:,2));
    data2.AllChannelsHbr = squeeze(filtered_dc2(:,:,2));

elseif ssc_reg == false
    disp('SSC regression is false, either because it was set 0 or because neither participant has a valid SSC data. Processing raw signal.')

    dc_orig1 = data1.procResult.dc(:,1:2,:);
    dc1 = permute(dc_orig1, [1, 3, 2]); 

    dc_orig2 = data2.procResult.dc(:,1:2,:);
    dc2 = permute(dc_orig2, [1, 3, 2]);

    data1.AllChannelsHbo = squeeze(dc1(:,:,1));
    data2.AllChannelsHbo = squeeze(dc2(:,:,1));

    data1.AllChannelsHbr = squeeze(dc1(:,:,2));
    data2.AllChannelsHbr = squeeze(dc2(:,:,2));

else
    error('ERROR: Invalid value for SSC_reg argument. Valid options are: true or false.');
end

assignin('base', 'data1', data1)
assignin('base', 'data2', data2)