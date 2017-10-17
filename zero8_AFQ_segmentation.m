function zero9_AFQ_semgentation(subj, bvals)
% This function will run Automated Fiber Quantification (Yeatman et al, 2014; https://github.com/yeatmanlab/AFQ). It will take
% the LiFE-evaluated ensemble connectome and segment streamlines into 20 major fiber tracts.
%
% Input is subject name.
%
% Output is a post_life_afq_fg.mat file that will contain fg, fg_classified, classification, and fascicles structures.  Fg_classified
% holds the major tracts, and will be used to generate tract profiles.
%
% 2017 Brad Caron Indiana University, Pestilli Lab

% Set paths
topDir = ['/N/dc2/projects/lifebid/Concussion/concussion_real/' subj '/diffusion_data/' bvals];
saveDir = fullfile(topDir, 'life', 'life_brainmask_full_ensemble');
mkdir(saveDir);

% AFQ full ensemble
if exist(fullfile(saveDir, 'post_afq_fg_nonsplit.mat'),'file')
    display('file exists. moving on')
else
    dtFile = fullfile(topDir, 'dt6', 'dti64trilin', 'dt6.mat');
    wholeBrainConnectome = fullfile(saveDir, 'optimized_life_connectome_1.mat');
    [fascicles,classification,fg,fg_classified] = feAfqSegment(dtFile, wholeBrainConnectome);
    save(fullfile(saveDir, 'post_afq_fg_nonsplit'), 'fg', 'fg_classified', 'classification', 'fascicles', '-v7.3');
    clear('dtFile','wholeBrainConnectome','fascicles','classification','fg','fg_classified');
end

% AFQ deterministic
if exist(fullfile(saveDir, 'post_afq_fg_nonsplit_deterministic.mat'),'file')
    display('file exists. moving on')
else
    dtFile = fullfile(topDir, 'dt6', 'dti64trilin', 'dt6.mat');
    wholeBrainConnectome = fullfile(saveDir, 'optimized_life_connectome_1_deterministic.mat');
    [fascicles,classification,fg,fg_classified] = feAfqSegment(dtFile, wholeBrainConnectome);
    save(fullfile(saveDir, 'post_afq_fg_nonsplit_deterministic'), 'fg', 'fg_classified', 'classification', 'fascicles', '-v7.3');
    clear('dtFile','wholeBrainConnectome','fascicles','classification','fg','fg_classified');
end

% AFQ probabilistic
if exist(fullfile(saveDir, 'post_afq_fg_nonsplit_probabilistic.mat'),'file')
    display('file exists. moving on')
else
    dtFile = fullfile(topDir, 'dt6', 'dti64trilin', 'dt6.mat');
    wholeBrainConnectome = fullfile(saveDir, 'optimized_life_connectome_1_probabilistic.mat');
    [fascicles,classification,fg,fg_classified] = feAfqSegment(dtFile, wholeBrainConnectome);
    save(fullfile(saveDir, 'post_afq_fg_nonsplit_probabilistic'), 'fg', 'fg_classified', 'classification', 'fascicles', '-v7.3');
    clear('dtFile','wholeBrainConnectome','fascicles','classification','fg','fg_classified');
end

% AFQ tensor
if exist(fullfile(saveDir, 'post_afq_fg_nonsplit_tensor.mat'),'file')
    display('file exists. moving on')
else
    dtFile = fullfile(topDir, 'dt6', 'dti64trilin', 'dt6.mat');
    wholeBrainConnectome = fullfile(saveDir, 'optimized_life_connectome_1_tensor.mat');
    [fascicles,classification,fg,fg_classified] = feAfqSegment(dtFile, wholeBrainConnectome);
    save(fullfile(saveDir, 'post_afq_fg_nonsplit_tensor'), 'fg', 'fg_classified', 'classification', 'fascicles', '-v7.3');
    clear('dtFile','wholeBrainConnectome','fascicles','classification','fg','fg_classified');
end
clear;
end
