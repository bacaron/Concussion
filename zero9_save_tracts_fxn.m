function one0save_tracts_fxn(subj, bvals)
% This function will save down the fiber tracts segmented from AFQ into a .mat file which can then be used to view fibers overlayed over T1.
%
% Input is subject name and bvals.  Output are folders for each different AFQ segmentation (with life, sans life, deterministic, probabilistic) and .mat files of the fibers.
%
% 2017 Brad Caron Indiana University, Pestilli Lab

projdir1 = ['/N/dc2/projects/lifebid/Concussion/concussion_real/' subj '/diffusion_data/' bvals];
saveDir = fullfile(projdir1, 'life', 'life_brainmask_full_ensemble');
mkdir(saveDir);
mkdir(fullfile(saveDir, 'post_afq_fg_nonsplit'));
mkdir(fullfile(saveDir, 'post_afq_fg_nonsplit_deterministic'));
mkdir(fullfile(saveDir, 'post_afq_fg_nonsplit_probabilistic'));
mkdir(fullfile(saveDir, 'post_afq_fg_nonsplit_tensor'));

% full ensemble
if exist(fullfile(saveDir, 'post_afq_fg_nonsplit','Callosum Forceps Major.mat'), 'file')
    display('file exists. moving on');
else
	load(fullfile(saveDir, 'post_afq_fg_nonsplit.mat'));

	for ii = 1:20
    		fgWrite(fg_classified(ii), fullfile(saveDir, 'post_afq_fg_nonsplit',fg_classified(ii).name), 'mat')
	end
	clear('classification', 'fascicles', 'fg', 'fg_classified', 'ii');
end

% deterministic
if exist(fullfile(saveDir,'post_afq_fg_nonsplit_deterministic','Callosum Forceps Major.mat'), 'file')
	display('file exists. moving on');
else
    load(fullfile(saveDir, 'post_afq_fg_nonsplit_deterministic.mat'));

	for ii = 1:20
   		fgWrite(fg_classified(ii), fullfile(saveDir, 'post_afq_fg_nonsplit_deterministic',fg_classified(ii).name), 'mat')
	end
	clear('classification', 'fascicles', 'fg', 'fg_classified', 'ii');
end

% probabilistic
if exist(fullfile(saveDir,'post_afq_fg_nonsplit_probabilistic','Callosum Forceps Major.mat'), 'file')
	display('file exists. moving on');
else
    load(fullfile(saveDir, 'post_afq_fg_nonsplit_probabilistic.mat'));

	for ii = 1:20
   		fgWrite(fg_classified(ii), fullfile(saveDir, 'post_afq_fg_nonsplit_probabilistic',fg_classified(ii).name), 'mat')
	end
	clear('classification', 'fascicles', 'fg', 'fg_classified', 'ii');
end

% tensor
if exist(fullfile(saveDir, 'post_afq_fg_nonsplit_tensor','Callosum Forceps Major.mat'), 'file')
	display('file exists. done');
else
    load(fullfile(saveDir, 'post_afq_fg_nonsplit_tensor.mat'));

	for ii = 1:20
    		fgWrite(fg_classified(ii), fullfile(saveDir, 'post_afq_fg_nonsplit_tensor',fg_classified(ii).name), 'mat')
	end
	clear('classification', 'fascicles', 'fg', 'fg_classified', 'ii');
end
clear
end
