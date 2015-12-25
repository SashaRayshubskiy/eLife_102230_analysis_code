%% Load imaging and behavioral data
global slash;

if isunix() == 1
slash = '\';
else
slash = '/';
end

% Must end with a slash
datapath = '';

analysis_path = [datapath slash 'analysis'];

if(~exist(analysis_path, 'dir'))
    mkdir(analysis_path);
end

% Load behavioral data
bdata_path = [datapath  slash 'ball' slash ];
bdata = load_behavioral_data(bdata_path);

% Load imaging data
cdata_path = [datapath  slash '2p' slash ];
cdata = load_imaging_data(cdata_path);

check_behavioral_and_imaging_data_trial_consistency();

%% Display behavioral data
