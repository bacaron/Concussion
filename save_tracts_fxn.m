function save_tracts_fxn(subj, bvals)
% This function will save down the fiber tracts segmented from AFQ into a .mat file which can then be used to view fibers overlayed over T1.
%
% Input is subject name and bvals.  Output are folders for each different AFQ segmentation (with life, sans life, deterministic, probabilistic) and .mat files of the fibers.
%
% 2017 Brad Caron Indiana University, Pestilli Lab

projdir1 = ['/N/dc2/projects/lifebid/Concussion/concussion_test/' subj '/diffusion_data/' bvals];
cd(fullfile(projdir1, 'life'));
mkdir('post_afq_fg')
mkdir('post_afq_fg_sanslife');

% tracts with life
cd post_afq_fg;
load(fullfile(projdir1, 'life', 'post_afq_fg.mat'))

for ii = 1:20
    fgWrite(fg_classified(ii), fg_classified(ii).name, 'mat')
end
clear('classification', 'fascicles', 'fg', 'fg_classified', 'ii')

% tracts sans life
cd(fullfile(projdir1, 'life', 'post_afq_fg_sanslife'))
load(fullfile(projdir1, 'life', 'post_afq_fg_sanslife.mat'))

for ii = 1:20
    fgWrite(fg_classified(ii), fg_classified(ii).name, 'mat')
end
clear('classification', 'fascicles', 'fg', 'fg_classified', 'ii')
