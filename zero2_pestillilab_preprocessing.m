function pestillilab_preprocessing_pipeline(subj)
%
% (1) Split multshell data into single shells. 
% (2) Normalize BVALUES/BVECS (change small values to 0, round to thousands)
% (3) AC-PC align T1w
% (4) dtiInit
% (5) AFQ 

% Add repositories.
restoredefaultpath
addpath(genpath('/N/dc2/projects/lifebid/code/franpest/AFQ'))
addpath(genpath('/N/dc2/projects/lifebid/code/franpest/encode'))
addpath(genpath('/N/dc2/projects/lifebid/code/vistasoft'))
addpath(genpath('/N/dc2/projects/lifebid/Concussion/concussion_real'))

% Prepare data paths and paths
b_vals = {'1000', '2000'};
paths.dirs.subj        = subj;
paths.dirs.base        = '/N/dc2/projects/lifebid/Concussion/concussion_real/';
paths.dirs.in_original_hcp_anat  = 'diffusion_directory/Anatomy';
paths.dirs.in_original_hcp_dwi   = 'diffusion_directory/diffusion';

paths.dirs.out_data    = 'diffusion_data';
paths.dirs.out_backup  = 'data_backup';
paths.dirs.out_anatomy = 'anatomy';
paths.dirs.out_tractography = 'fibers';
paths.dirs.out_dt = 'dt6';

paths.file.names.dwi     = 'data.nii.gz';
paths.file.names.anatomy = 't1.nii.gz';
paths.file.names.bvecs   = 'data.bvecs';
paths.file.names.bvals   = 'data.bvals';
paths.file.names.params  = 'params';
paths.file.names.normalization = 'bvals_bvecs_normalization.mat';

paths.file.dwi  = fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.in_original_hcp_dwi,  paths.file.names.dwi);
paths.file.bvec = fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.in_original_hcp_dwi, paths.file.names.bvecs);
paths.file.bval = fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.in_original_hcp_dwi,  paths.file.names.bvals);
paths.file.t1f  = fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.in_original_hcp_anat, paths.file.names.anatomy);

% Parameters used for normalization
params.single_shells       = [0,1000,2000];
params.thresholds.b0_normalize    = 200;
params.thresholds.bvals_normalize = 100;
params.flip_y = true;

%% Organize folders
cd(fullfile(paths.dirs.base,paths.dirs.subj))

eval(sprintf('!mkdir -v %s',paths.dirs.out_backup))

eval(sprintf('!cp -v %s/*',paths.dirs.in_original_hcp_anat, ...
                            paths.dirs.out_backup));%T1w original_hcp_data/
                        
eval(sprintf('!mkdir -v %s',paths.dirs.out_anatomy));% anatomy

eval(sprintf('!cp -v %s/*nii.gz %s',fullfile(paths.dirs.in_original_hcp_anat),  ...
                                             paths.dirs.out_anatomy)); %original_hcp_data/T1w/*nii.gz* anatomy/
                                         
eval(sprintf('!mkdir -v %s',paths.dirs.out_data));%diffusion_data

eval(sprintf('!cp -v %s/* %s',paths.dirs.in_original_hcp_dwi, ...
                            paths.dirs.out_backup));% Diffusion_7T original_hcp_data/                        
                        
eval(sprintf('!cp -v %s/runeddy.nii.gz %s/data.nii.gz',fullfile(paths.dirs.in_original_hcp_dwi), ...
                                     paths.dirs.out_data));% original_hcp_data/Diffusion_7T/* diffusion_data/

eval(sprintf('!cp -v %s/runeddy.eddy_rotated_bvecs %s/data.bvecs',fullfile(paths.dirs.in_original_hcp_dwi), ...
                                     paths.dirs.out_data));% original_hcp_data/Diffusion_7T/* diffusion_data/

eval(sprintf('!cp -v %s/bvals %s/data.bvals',fullfile(paths.dirs.in_original_hcp_dwi), ...
                                     paths.dirs.out_data));% original_hcp_data/Diffusion_7T/* diffusion_data/

eval(sprintf('!mkdir -v %s', fullfile(paths.dirs.out_data, paths.dirs.out_tractography)));% fibers

if exist('release-notes','file')
   eval(sprintf('!cp -v %s %s','release-notes', ...
        paths.dirs.out_backup)); % release-notes original_hcp_data/
end

cd(paths.dirs.out_data);% diffusion_data

%% Normalize HCP files to the VISTASOFT environment
%
% Split data into two separate paths (BVALS = 1000 and 2000).
bvals.val = dlmread(paths.file.names.bvals);

% Round the numbers to the closest thousand 
% This is necessary because the VISTASOFT software does not handle the B0
% when they are not rounded.
[bvals.unique, ~, bvals.uindex] = unique(bvals.val);
if ~isequal( bvals.unique, params.single_shells)
    bvals.unique(bvals.unique <= params.thresholds.b0_normalize) = 0;
    bvals.unique  = round(bvals.unique./params.thresholds.bvals_normalize) ...
        *params.thresholds.bvals_normalize;
    bvals.valnorm = bvals.unique( bvals.uindex );
    dlmwrite(paths.file.names.bvals,bvals.valnorm);
    save(paths.file.names.normalization,'bvals')
else
    bvals.valnorm = bvals.val;
end
% If not exist these files already
index1000 = (bvals.valnorm == params.single_shells(2));
index2000 = (bvals.valnorm == params.single_shells(3));
index0    = (bvals.valnorm == params.single_shells(1));

% Find all indices to each bvalue and B0
all_1000  = or(index1000,index0);
all_2000  = or(index2000,index0);

% Validate that we selected the correct number of bvals+b0
assertEqual(sum(all_1000), sum(index0)+sum(index1000))
assertEqual(sum(all_2000), sum(index0)+sum(index2000))

% Write bvals to disk
bvals1000 = bvals.valnorm(all_1000);
dlmwrite('data_b1000.bvals',bvals1000);

bvals2000 = bvals.valnorm(all_2000);
dlmwrite('data_b2000.bvals',bvals2000);

% Work on BVECS
% negative x and y
bvecs1000 = dlmread('data.bvecs');
if ~(size(bvecs1000,1) == 3), bvecs1000 = bvecs1000'; end
bvecs1000 = bvecs1000(:,all_1000);
if params.flip_y
   bvecs1000(1,:) = -bvecs1000(1,:);
   bvecs1000(2,:) = -bvecs1000(2,:);
end
dlmwrite('data_b1000.bvecs',bvecs1000);

bvecs2000 = dlmread('data.bvecs');
if ~(size(bvecs2000,1) == 3), bvecs2000 = bvecs2000'; end
bvecs2000 = bvecs2000(:,all_2000);
if params.flip_y
   bvecs2000(1,:) = -bvecs2000(1,:);
   bvecs2000(2,:) = -bvecs2000(2,:);
end
dlmwrite('data_b2000.bvecs',bvecs2000);

% Split the data into single shells
dwi   = niftiRead('data.nii.gz');
dwi1000 = dwi;
dwi1000.fname = 'data_b1000.nii.gz';

dwi2000 = dwi;
dwi2000.fname = 'data_b2000.nii.gz';

% Remove the B-val = 2000 data to make a b-val = 1000 dataset
dwi1000.data   = dwi.data(:,:,:,all_1000);
dwi1000.dim(4) = size(dwi1000.data,4);
niftiWrite(dwi1000);

dwi2000.data   = dwi2000.data(:,:,:,all_2000);
dwi2000.dim(4) = size(dwi2000.data,4);
niftiWrite(dwi2000);



% %% AC-PC align
% %!module load spm/8
% 
% % Make sure the file is aligned properly
ni_anat_file = fullfile(paths.dirs.base, paths.dirs.subj, ...
                        paths.dirs.out_anatomy, ...
                        paths.file.names.anatomy);
if exist(ni_anat_file,'file')
    ni_anat = niftiRead( ni_anat_file );
    ni_anat = niftiApplyCannonicalXform( ni_anat );
else
    keyboard
end

% % Load a standard template from vistasoft
MNI_template =  fullfile(mrDiffusionDir, 'templates', 'MNI_T1.nii.gz');

% % Compute the spatial normalization to align the current raw data to the template
SpatialNormalization = mrAnatComputeSpmSpatialNorm(ni_anat.data, ni_anat.qto_xyz, MNI_template);

% % Assume that the AC-PC coordinates in the template are in a specific location:
% % X, Y, Z = [0,0,0; 0,-16,0; 0,-8,40]
% % Use this assumption and the spatial normalization to extract the corresponding AC-PC location on the raw data
coords = [0,0,0; 0,-16,0; 0,-8,40];

ImageCoords = mrAnatGetImageCoordsFromSn(SpatialNormalization, tal2mni(coords)', true)';

% % Now we assume that ImageCoords contains the AC-PC coordinates that we need for the Raw data. 
% % We will use them to compute the AC_PC alignement automatically. The new file will be saved to disk. 
% % Check the alignement.
mrAnatAverageAcpcNifti(ni_anat, 't1_acpc.nii.gz', ImageCoords, [], [], [], false);                                        

dwi = niftiRead('data.nii.gz');
res=dwi.pixdim(1:3);
clear dwi;

% Organizing files before DTIINIT
eval(sprintf('!mkdir -v %s/1000', fullfile(paths.dirs.base, ...
                                            paths.dirs.subj, ...
                                            paths.dirs.out_data)))

eval(sprintf('!mkdir -v %s/2000', fullfile(paths.dirs.base, ...
                                            paths.dirs.subj, ...
                                            paths.dirs.out_data)))
eval(sprintf('!cp -v %s/t1_acpc.nii.gz %s', fullfile(paths.dirs.base, ...
                                                      paths.dirs.subj, ...
                                                      'diffusion_data', ...
                                                      '1000'), ...
                                                  fullfile(paths.dirs.base, ...
                                                      paths.dirs.subj, ...
                                                      paths.dirs.out_data, ...
                                                      '1000')))
eval(sprintf('!cp -v %s/t1_acpc.nii.gz %s', fullfile(paths.dirs.base, ...
                                                      paths.dirs.subj, ...
                                                      'diffusion_data', ...
                                                      '2000'), ...
                                            fullfile(paths.dirs.base, ...
                                                      paths.dirs.subj, ...
                                                      paths.dirs.out_data, ...
                                                      '2000')))
eval(sprintf('!mv -v %s/*b1000* %s', fullfile(paths.dirs.base, paths.dirs.subj, ...
                                              paths.dirs.out_data), ...
                                              fullfile(paths.dirs.base, ...
                                              paths.dirs.subj, paths.dirs.out_data, ...
                                              '1000')))
eval(sprintf('!mv -v %s/*b2000* %s', fullfile(paths.dirs.base, paths.dirs.subj, ...
                                              paths.dirs.out_data), ...
                                              fullfile(paths.dirs.base, ...
                                              paths.dirs.subj, paths.dirs.out_data, ...
                                              '2000')))

% DTIINIT
% 1000
cd(fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.out_data, '1000'))
% create output directory
mkdir(paths.dirs.out_dt);

% create parameters
dwParams = dtiInitParams;
dwParams.eddyCorrect       = -1;
dwParams.phaseEncodeDir    = 2; 
dwParams.rotateBvecsWithRx = 0;
dwParams.rotateBvecsWithCanXform = 0;
dwParams.dwOutMe = res;
dwParams.bvecsFile = fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.out_data, '1000/data_b1000.bvecs');
dwParams.bvalsFile = fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.out_data, '1000/data_b1000.bvals');
dwParams.outDir    = fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.out_data, '1000/dt6');
dt_path = dtiInit(fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.out_data, '1000', 'data_b1000.nii.gz'), fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.out_data, '1000', 't1_acpc.nii.gz'), dwParams);

% 2000
cd(fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.out_data, '2000'))
% create output directory
mkdir(paths.dirs.out_dt);

% create parameters
dwParams = dtiInitParams;
dwParams.eddyCorrect       = -1;
dwParams.phaseEncodeDir    = 2; 
dwParams.rotateBvecsWithRx = 0;
dwParams.rotateBvecsWithCanXform = 0;
dwParams.dwOutMe = res;
dwParams.bvecsFile = fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.out_data, '2000/data_b2000.bvecs');
dwParams.bvalsFile = fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.out_data, '2000/data_b2000.bvals');
dwParams.outDir    = fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.out_data, '2000/dt6');
dt_path = dtiInit(fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.out_data, '2000', 'data_b2000.nii.gz'), fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.out_data, '2000', 't1_acpc.nii.gz'), dwParams);

%% Create an MRTRIX .b file from the bvals/bvecs of each shell
for ibv = 1:length(b_vals)
    bvecs = fullfile(sprintf('/N/dc2/projects/lifebid/Concussion/concussion_real/%s/diffusion_data/%s/dt6/data_b%s_aligned_trilin_noMEC.bvecs', subj, b_vals{ibv}, b_vals{ibv}));
    bvals = fullfile(sprintf('/N/dc2/projects/lifebid/Concussion/concussion_real/%s/diffusion_data/%s/dt6/data_b%s_aligned_trilin_noMEC.bvals', subj, b_vals{ibv}, b_vals{ibv}));
    out   = fullfile(sprintf('/N/dc2/projects/lifebid/Concussion/concussion_real/%s/diffusion_data/fibers/data_b%s.b', subj, b_vals{ibv}));
    mrtrix_bfileFromBvecs(bvecs, bvals, out);
end

%% AFQ
% 1000
%dt = dtiLoadDt6(fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.out_data, '1000', 'dt6', 'dti64trilin', 'dt6.mat'));

%fg = AFQ_WholebrainTractography(dt);
%[fg_classified,~,classification,fg]= AFQ_SegmentFiberGroups(dt, fg, [], [], false);
%fascicles = fg2Array(fg_classified)
%save(fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.out_data, '1000', 'fgFibers'), 'fg', 'fg_classified', 'classification', 'fascicles', '-v7.3');

% 2000
%dt = dtiLoadDt6(fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.out_data, '2000', 'dt6', 'dti64trilin', 'dt6.mat'));

%fg = AFQ_WholebrainTractography(dt);
%[fg_classified,~,classification,fg]= AFQ_SegmentFiberGroups(dt, fg, [], [], false);
%fascicles = fg2Array(fg_classified)
%save(fullfile(paths.dirs.base, paths.dirs.subj, paths.dirs.out_data, '2000', 'fgFibers'), 'fg', 'fg_classified', 'classification', 'fascicles', '-v7.3');
end
