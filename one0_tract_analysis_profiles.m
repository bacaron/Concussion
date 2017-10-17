function one1_tract_analysis_profiles(subj, bvals, fiberFile, tensor, tract)
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
dt6path = ['/N/dc2/projects/lifebid/Concussion/concussion_real/' subj '/diffusion_data/' bvals '/dti64trilin/dt6.mat'];
fiberGroupPath = ['/N/dc2/projects/lifebid/Concussion/concussion_real/' subj '/diffusion_data/' bvals '/life/' fiberFile];
anatomyFilePath = ['/N/dc2/projects/lifebid/Concussion/concussion_real/' subj '/diffusion_data/' bvals '/t1_acpc_bet.nii.gz'];

% 1. Load the fiber groups and dt6 file.
load( fiberGroupPath );
dt = dtiLoadDt6( dt6path );
fg = fg_classified( tract );

% 2. compute the core fiber from the fiber group (the tact profile is computed here)
[fa, md, rd, ad, cl, core, fgResampled] = dtiComputeDiffusionPropertiesAlongFG( fg, dt,[],[],200);

% 3. Select a center portion fo the tract and show the FA and MD values 
% normally we only use for analyses the middle most reliable portion of the fiber.
nodesToPlot = 50:151;

h.tpfig = figure('name', 'My tract profile','color', 'w', 'visible', 'off');
tract_profile = plot(tensor(nodesToPlot),'color', [0.2 0.2 0.9],'linewidth',4)
set(gca, 'fontsize',20, 'box','off', 'TickDir','out', ...
	'xticklabel',{'Tract begin','Tract end'},'xlim',[0 100],'ylim',[0.00 1.00],'Ytick',[0 .25 .5 .75],'Xtick',[0 100])
Title_plot = title(fg.name)
xlabel('Location on tract')
ylh = ylabel('Fractional Anisotropy');
        
saveas(tract_profile, fullfile(fiberGroupPath,sprintf('/diffusion_data/%s/life/%s/', bvals, fiberFile), Title_plot.String), 'png')
saveas(tract_profile, fullfile(fiberGroupPath,sprintf('/diffusion_data/%s/life/%s/', bvals, fiberFile), Title_plot.String), 'eps')
        clear()
