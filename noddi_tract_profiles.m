function noddi_tract_profiles(subj, bvals, noddiFile, fiberFile, tract)
% Tract Analysis.
%
% This script is an example of analysis of tracts that uses the profile of the tracts
% along their length and plots values of FA, MD etc. as function of position along a tract length. 
%
% Copyright Franco Pestilli Indiana University 2016
%
% Dependencies:
% addpath(genpath('~/path/to/spm8'))       % -> SPM/website
% addpath(genpath('~/path/to/vistasoft')); % -> https://www.github.com/vistalab/vistasoft
% addpath(genpath('~/path/to/life'));      % -> https://www.github.com/francopestilli/life
% addpath(genpath('~/path/to/mba'));       % -> https://www.github.com/francopestilli/mba


% 0. Set up paths to tracts files and dt6 file.
noddiPath = ['/N/dc2/projects/lifebid/Concussion/concussion_real/NODDI/' subj '/AMICO/NODDI'];
fiberGroupPath = ['/N/dc2/projects/lifebid/Concussion/concussion_real/' subj '/diffusion_data/' bvals '/life/life_brainmask_full_ensemble_test'];
anatomyFilePath = ['/N/dc2/projects/lifebid/Concussion/concussion_real/' subj '/diffusion_data/' bvals '/t1_acpc.nii.gz'];
saveDir = fullfile(noddiPath, noddiFile);
mkdir(saveDir);

% 1. Load the fiber groups and dt6 file.
load( fullfile(fiberGroupPath, fiberFile ));
noddi = niftiRead(fullfile(noddiPath, noddiFile ));
fg = fg_classified( tract );
tract_number = num2str(tract);

% 2. compute the core fiber from the fiber group (the tact profile is computed here)
[fa, md, rd, ad, cl, SuperFile, fgClipped, cp, cs, fgResampled] = dtiComputeDiffusionPropertiesAlongFG( fg, noddi,[],[],200);

% 3. save file down for group stats
save(fullfile(saveDir, sprintf('/stats_for_group_%s', tract_number)), ...
'fa','core','-v7.3');

% 4. Select a center portion fo the tract and show the FA and MD values 
% normally we only use for analyses the middle most reliable portion of the fiber.
nodesToPlot = 50:151;

h.tpfig = figure('name', 'My tract profile','color', 'w', 'visible', 'on');
tract_profile = plot(fa(nodesToPlot),'color', [0.2 0.2 0.9],'linewidth',4)
set(gca, 'fontsize',20, 'box','off', 'TickDir','out', ...
	'xticklabel',{'Tract begin','Tract end'},'xlim',[0 100],'ylim',[0.00 1.00],'Ytick',[0 .25 .5 .75],'Xtick',[0 100])
hline = refline([0 mean(fa(nodesToPlot))]);
Title_plot = title(fg.name)
xlabel('Location on tract')
ylh = ylabel('ICVF');
        
saveas(tract_profile, [saveDir '/' Title_plot.String], 'eps');
saveas(tract_profile, [saveDir '/' Title_plot.String], 'png');
clear java;
% 
