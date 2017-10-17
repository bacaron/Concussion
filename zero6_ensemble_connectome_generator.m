function zero7_ensemble_connectome_generator(subj,bvals)
% This function will combine all of the outputs from mrtrix tractography into singular .mat file.  This will be used
% in LiFE and AFQ steps. Function created by Franco Pestilli.
%
% Inputs are subject name, and b-values for study (see below for examples).
% Project directory will need to be ammended for your file set-up. 
%
% Output is a *_ensemble.mat file.
%
% % 2017 Brad Caron Indiana University, Pestilli Lab

% build paths; these will need to be inputted into function
% subjects for study
% subj = '1_5';

% b-values for study
% bvals = {'1000','2000'}

% project direry for study
projdir1 = ['/N/dc2/projects/lifebid/Concussion/concussion_real/' subj '/diffusion_data/' bvals];
savedir = fullfile(projdir1, 'major_tracts', 'major_tracts_brainmask_full_ensemble');
mkdir(savedir);

% Curvature paramater (lmax)
lmaxparam = {'2','4','6','8','10','12'};
% probability or deterministic tracking from mrtrix
streamprob = {'PROB','STREAM'};

% Tensor-based tracking. only on fascicles group
fe_path = fullfile(projdir1);
fg = fgRead(fullfile(fe_path,'fibers_new','fibers_new_brainmask_full_ensemble','data_b1000_aligned_trilin_noMEC_wm_tensor-60000.tck'));
dt6File = fullfile(fe_path,'/dt6/dti64trilin/dt6.mat');
fasciclesClassificationSaveName = fullfile(savedir,'data_b1000_aligned_trilin_noMEC_wm_tensor-60000.mat');

% Tensor file
if exist(fullfile(savedir, 'data_b1000_aligned_trilin_noMEC_wm_tensor-60000.mat'), 'file')
    display('file exists. moving on');
else
    fgWrite(fg,fasciclesClassificationSaveName);
end

% Full ensemble
% CSD-based tracking. Load one at the time.
if exist(fullfile(savedir, 'data_b1000_aligned_trilin_noMEC_ensemble.mat'), 'file')
    display('file exists. moving on');
else
    for ilm = 1:length(lmaxparam)
        for isp = 1:length(streamprob)
                fe_path = fullfile(projdir1);
                fg_tmp = fgRead(fullfile(fe_path,'fibers_new','fibers_new_brainmask_full_ensemble',sprintf('data_b1000_aligned_trilin_noMEC_csd_lmax%s_wm_SD_%s-60000.tck',lmaxparam{ilm},streamprob{isp})));

                % Merge the new fiber group with the original fiber group.
                fg = fgMerge(fg,fg_tmp);
                clear fg_tmp
        end

        % Write fascicle group to disk.
        fgFileName = fullfile(savedir, sprintf('data_b1000_aligned_trilin_noMEC_ensemble.mat',lmaxparam{ilm}));
        fgWrite(fg,fgFileName);
    end
clear('fg','fgFileName');
end

% Deterministic
% CSD-based tracking. Load one at the time.
if exist(fullfile(savedir, 'data_b1000_aligned_trilin_noMEC_deterministic_ensemble.mat'), 'file')
    display('file exists. moving on');
else
    fg = fgRead(fullfile(fe_path,'fibers_new','fibers_new_brainmask_full_ensemble','data_b1000_aligned_trilin_noMEC_wm_tensor-60000.tck'));
    for ilm = 1:length(lmaxparam)
        fe_path = fullfile(projdir1);
        fg_tmp = fgRead(fullfile(fe_path,'fibers_new','fibers_new_brainmask_full_ensemble',sprintf('data_b1000_aligned_trilin_noMEC_csd_lmax%s_wm_SD_STREAM-60000.tck',lmaxparam{ilm})));

        % Merge the new fiber group with the original fiber group.
        fg = fgMerge(fg,fg_tmp);
        clear fg_tmp;

        % Write fascicle group to disk.
        fgFileName = fullfile(savedir, sprintf('data_b1000_aligned_trilin_noMEC_deterministic_ensemble.mat',lmaxparam{ilm}));
        fgWrite(fg,fgFileName);
    end
clear('fg','fgFileName');
end

Probabilistic
% CSD-based tracking. Load one at the time.
if exist(fullfile(savedir, 'data_b1000_aligned_trilin_noMEC_probabilistic_ensemble.mat'), 'file')
    display('file exists. moving on');
else
    fg = fgRead(fullfile(fe_path,'fibers_new','fibers_new_brainmask_full_ensemble_test','data_b1000_aligned_trilin_noMEC_wm_tensor-60000.tck'));
    for ilm = 1:length(lmaxparam)
        fe_path = fullfile(projdir1);
        fg_tmp = fgRead(fullfile(fe_path,'fibers_new','fibers_new_brainmask_full_ensemble',sprintf('data_b1000_aligned_trilin_noMEC_csd_lmax%s_wm_SD_PROB-60000.tck',lmaxparam{ilm})));

        % Merge the new fiber group with the original fiber group.
        fg = fgMerge(fg,fg_tmp);
        clear fg_tmp;

        % Write fascicle group to disk.
        fgFileName = fullfile(savedir, sprintf('data_b1000_aligned_trilin_noMEC_probabilistic_ensemble.mat',lmaxparam{ilm}));
        fgWrite(fg,fgFileName);
    end
end
clear all;
end

