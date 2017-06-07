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

% set paths
dataPath = ['/N/dc2/projects/lifebid/Concussion/concussion_real/' subj '/diffusion_data/' bvals];
dtFilePath = fullfile(dataPath, 'dt6', 'dti64trilin');
wholeBrainConnectomePath = fullfile(datapath, 'life');

% AFQ
% These file structures will need to be changed for your file set-up.
dtFile = fullfile(dtFilePath, 'dt6.mat');
wholeBrainConnectome = fullfile(wholeBrainConnectomePath, 'optimized_life_connectome_1.mat')
[fascicles,classification,fg,fg_classified] = feAfqSegment(dtFile, wholeBrainConnectome);
save(fullfile(wholeBrainConnectomePath,'/life/post_afq_fg'), 'fg', 'fg_classified', 'classification', 'fascicles', '-v7.3');
