% Manually build an IMRO table, if the SpikeGLX GUI imro builder isn't
% sufficiently customizable for your needs
% Run each section independently or fill values at 1. and 3. and run script
% Designed for 2.0 4-shank probes
% Adapted entirely from Charlie Walter's IMRO GUI
% (https://github.com/ckwalters/npx2_site_gui)

%% 1. User inputs
directory = 'Z:/WT_Sequences/2023_spring/Raw_Data/Neural_Traces/IMRO';
name = 'Noether_HPC';
% if you already have an IMRO that you want to modify or check multiplex
% compatibility of, name it here, otherwise leave as an empty string
load_imro = 'Z:/WT_Sequences/2023_spring/Raw_Data/Neural_Traces/IMRO/Noether_HPC.imro';
ref_electrode = 0; % 0 is external, 1-4 are the tips

% defaults for 2.0 4-shank probe
probe_type = 24;
n_channels = 384;
n_shanks = 4;
n_electrodes_per_shank = 1280;

% Initialize
active_electrodes = zeros(n_shanks,n_electrodes_per_shank);
% Calculate all the electrode,shank pairings that belong to a particular channel
channel_map = zeros(n_shanks,n_electrodes_per_shank);
bank_map = zeros(n_shanks,n_electrodes_per_shank);
for shankInd = 1:n_shanks
    for elecInd = 1:n_electrodes_per_shank
        [bank_map(shankInd,elecInd), channel_map(shankInd,elecInd)] = ElecToChan(shankInd-1,elecInd-1);
    end
end

%% 2. Load IMRO, if specified
if ~isempty(load_imro)
    opts = delimitedTextImportOptions('Delimiter',{'(',')'},...
        'ConsecutiveDelimitersRule','join',...
        'LeadingDelimitersRule','ignore',...
        'TrailingDelimitersRule','ignore');
    imro = readtable(load_imro,opts);
    imro = imro(1,2:end);
    imro = table2array(imro);
    for entry = 1:n_channels
        data = imro{entry}(length(char(string(entry-1)))+2:end);
        shnk = str2num(data(1));
        chan = str2num(data(7:end));
        active_electrodes(shnk+1,chan+1) = 1; 
    end
end

%% 3. Manually specify which electrodes will be recorded from
% e.g.: active_electrodes(shank,sites) = 1;
% active_electrodes(1,33:80) = 1;
% active_electrodes(1,113:192) = 1;
% active_electrodes(2,33:208) = 1;
% active_electrodes(3,97:128) = 1;
% active_electrodes(3,177:224) = 1;

%% 4. Confirm these sites can be multiplexed
% Find which channels cannot be multiplexed with each site
% channel_lookup = cell(1,n_channels);
% for channelInd = 1:n_channels
%     [shanks,electrodes] = find(channel_map==(channelInd-1));
%     shanks = shanks-1; % find function produces indices that are 1-indexed
%     electrodes = electrodes-1;
%     channel_lookup{channelInd} = [shanks electrodes]; % either 13 or 14 electrodes per channel
% end

% Check if enough channels specified
if sum(sum(active_electrodes))~=n_channels
    error('%d channels specified. Please specify %d channels.',sum(sum(active_electrodes)),n_channels)
end

% Find info needed for IMRO table
imrows = nan(n_channels,5);
[imrows(:,2),imrows(:,5)] = find(active_electrodes);
imrows(:,2) = imrows(:,2)-1;
imrows(:,5) = imrows(:,5)-1;
for channel = 1:n_channels
    imrows(channel,1) = channel_map(imrows(channel,2)+1,imrows(channel,5)+1);
    imrows(channel,3) = bank_map(imrows(channel,2)+1,imrows(channel,5)+1);
end
imrows = sortrows(imrows,1);
imrows(:,4) = ref_electrode;

% Check if any sites are multiplexed to the same channel
for channel = 1:n_channels
    if imrows(channel,1)+1 ~= channel
        error('Channel %d is allocated more than once.', imrows(channel,1))
    end
end

%% 5. Write the IMRO table
% open a new file
fileName = [directory, '\', name, '.imro'];
nmID = fopen(fileName,'w');
fprintf(nmID,'(%d,%d)', 24, 384); % first entry (probe type & n chan)
for channel = 1:n_channels
    fprintf(nmID,'(%d %d %d %d %d)', imrows(channel,1), imrows(channel,2),...
        imrows(channel,3), imrows(channel,4), imrows(channel,5) );
end
fprintf(nmID, '\n');
fclose(nmID);

%% Functions
function [bank, chan] = ElecToChan( shank, elecInd )
    % electrode to channel map
    elecMap = zeros(4,8,'single');
    elecMap(1,:) = [0,2,4,6,5,7,1,3];
    elecMap(2,:) = [1,3,5,7,4,6,0,2];
    elecMap(3,:) = [4,6,0,2,1,3,5,7];
    elecMap(4,:) = [5,7,1,3,0,2,4,6];
    
    bank = floor(elecInd/384);
    
    % which block within the bank?
    blockID = floor((elecInd - bank*384)/48);
    
    % which channel within the block?
    subBlockInd = mod((elecInd - bank*384), 48);
    
    % get the channel number for that shank, bank, and block combo
    chan = 48*elecMap(shank+1, blockID+1) + subBlockInd;
end