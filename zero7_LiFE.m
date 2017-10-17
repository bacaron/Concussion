function zero8_LiFE(subj, bvals)
%
% LiFE function (Pestilli et al, 2014, Nature Methods).  This will evaluate ensemble connectome generated from 
% ensemble_connectome_generator.m (github.com/brain-life/pestillilab-projects/Concussion/zero7_ensemble_connectome_generator.m).
% It will fit the LiFE model, determine which fascicles do not contribute positvely to overall diffusion signal, and remove
% those fascicles from the connectome.
% 
% Inputs are subject name, and a bval
%
% Output will be a .mat file called 'optimized_life_connectome_1'. This will then be used in AFQ.
%
% 2017 Brad Caron Indiana University, Pestilli Lab

load(['/N/dc2/projects/lifebid/Concussion/concussion_real/' subj '/diffusion_data/' bvals '/major_tracts/major_tracts_brainmask_full_ensemble/data_b' bvals '_aligned_trilin_noMEC_ensemble.mat']);

% edit these for your file structure
topdir = ['/N/dc2/projects/lifebid/Concussion/concussion_real/' subj '/diffusion_data/' bvals];
anatomyFile = fullfile(topdir, 't1_acpc.nii.gz');
feFileName = strcat('data', '_', 'b', bvals, '_', 'aligned', '_', 'trilin','_', 'noMEC', '_', 'ensemble', '_', 'FE');
dwiFile = fullfile(topdir, 'dt6', strcat('data', '_', 'b', bvals, '_', 'aligned', '_', 'trilin' , '_', 'noMEC', '.', 'nii', '.', 'gz'));
fgFileName = fg;
fe_path = fullfile(topdir, 'major_tracts', 'major_tracts_brainmask_full_ensemble');
dwiFileRepeated = [];
savedir = fullfile(topdir, 'life', 'life_brainmask_full_ensemble_test');
mkdir(savedir);

% Build the life model into an fe structure
L = 360; % Discretization parameter
fe = feConnectomeInit(dwiFile,fgFileName,feFileName,savedir,dwiFileRepeated,anatomyFile,L,[1,0]);

% Fit the model.
Niter = 500;
fit = feFitModel(feGet(fe,'model'),feGet(fe,'dsigdemeaned'),'bbnnls',Niter,'preconditioner');

% Install the fit to the FE structure
fe = feSet(fe,'fit',fit);

% Save the life model to disk; again, edit for your file structure
save(fullfile(savedir, 'optimized_life_1'), 'fe', '-v7.3');

% Extract the fascicles with positive weight
w  = feGet(fe,'fiber weights'); % Collect the weights associated with each fiber
fg = feGet(fe,'fibers acpc');
fg = fgExtract(fg, w > 0, 'keep');
len_fg = length(fg.fibers);

% edit for your file structure
save(fullfile(savedir, 'optimized_life_connectome_1'), 'fg', '-v7.3');
save(fullfile(savedir, strcat(subj, '_', 'life_stats')), 'len_fg', '-v7.3');
clear
end
