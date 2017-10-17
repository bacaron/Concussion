function dwiAlignT1NODDI(subj)

%%  Will align multi-shell dwi to ACPC aligned T1 and rotate the bvecs.  Needs vistasoft ((C) Vista lab, Stanford University)

% Set directories
topdir = '/N/dc2/projects/lifebid/Concussion/concussion_real/NODDI';
dwiRaw = fullfile(topdir,subj,'data.nii.gz');
bvecs_pre = fullfile(topdir,subj,'bvecs');
b0 = fullfile(topdir,subj,'nodif_brain_preNODDI.nii.gz');
t1 = fullfile(topdir,subj,'t1_acpc.nii.gz');
outAcpcTransform = fullfile(topdir,subj,'acpcTransform.mat');
outMm = [2 2 2];
outDwi = fullfile(topdir,subj,'data_acpc.nii.gz');
outBvecs = fullfile(topdir,subj,'bvecs_rot');

% Get transformation matrix
dtiRawAlignToT1(b0,t1,outAcpcTransform);

% Run Alignment
dtiRawResample(dwiRaw,[],outAcpcTransform,outDwi,[0 0 0 0 0 0],[2 2 2]);

% Rotate bvecs
bvecs = dtiRawReorientBvecs(bvecs_pre,[],outAcpcTransform,outBvecs);
exit;
end
