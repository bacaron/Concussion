function microstructure_life_model(subj)
%% Example of script for the Microstructure model
% The Concussion model decomposes the observed diffusion signal in three
% components (Isotropic + Fibers + Microstructure):
% 1) The Isotropic signal due free motion of water molecules
% 2) The Anisotropic signal due to fascicle's (directional features)
% 3) The Anisotropic signal due to other directional features not related
% with the fascicles (microstructure)
%
% After the model is fit, this script generates separated nifti files for
% each of the components.
% Also, for each one of the 20 major tracts, we compute its prediction and
% generate individual nifti files with the diffusion signal associated only
% to the fibers in those corresponding tracts.
% 1 hr to run

%% Set the Path for your output
dataOutputPath = ['/N/dc2/projects/lifebid/Concussion/concussion_real/' subj '/diffusion_data/1000/Concussion_model/results/'];
mkdir(fullfile(dataOutputPath))

%% Set the proper path for VISTASOFT 
vista_soft_path = '/N/dc2/projects/lifebid/code/vistasoft/';
addpath(genpath(vista_soft_path));

%% Set the proper path for the ENCODE with Concussion model
ENCODE_path = '/N/dc2/projects/lifebid/code/bacaron/encode/encode_concussion/';
addpath(genpath(ENCODE_path));

%% Concussion dataset 
dataRootPath = '/N/dc2/projects/lifebid/Concussion/concussion_real';
%% subject = '1_005';
%% conn = 'NUM01'; % 
%% param = 'lmax8'; % {'lmax10','lmax12','lmax2','lmax4','lmax6','lmax8', ''}
%% alg = 'SD_PROB'; % {'SD_PROB', 'SD_STREAM','tensor'}
stem = 'data_b1000_aligned_trilin_noMEC';

% Generate fe_strucure
%% Build the file names for the diffusion data, the anatomical MRI.
%% dwiFile       = deblank(ls(fullfile(dataRootPath,subject,'diffusion_data','data_b1000_aligned.nii.gz')));
dwiFile       = fullfile(dataRootPath, subj, 'diffusion_data', bvals, 'dt6', 'data_b1000_aligned_trilin_noMEC.nii.gz');
t1File        = fullfile(dataRootPath, subj, 'diffusion_data', 't1_acpc.nii.gz'];
fgFileName    = fullfile(dataRootPath, subj, 'diffusion_data', bvals, 'major_tracts', 'major_tracts_brainmask_full_ensemble', 'data_b1000_aligned_trilin_noMEC_ensemble.mat'];
%%feFileName    = strcat(subject,'_',alg,'_',param,'_',conn);
feFileName    = strcat(subj,'_',stem,'_','ensemble','_','FE');

L = 360; % Discretization parameter
Niter = 500; % Number of iterations for the optimization in LiFE

%% Initialize the model
tic
fe_prob = feConnectomeInit(dwiFile,fgFileName,feFileName,[] ,[], t1File,L,[1,0],0); % We set dwiFileRepeat =  run 02
disp(' ')
disp(['Time for model construction ','(L=',num2str(L),')=',num2str(toc),'secs']);

%% Fit the LiFE model
fe_prob = feSet(fe_prob,'fit',feFitModel(feGet(fe_prob,'model'),feGet(fe_prob,'dsigdemeaned'),'bbnnls',Niter,'preconditioner'));
save(fullfile(dataOutputPath, 'fe_post_life.mat'), 'fe_prob', '-v7.3');
