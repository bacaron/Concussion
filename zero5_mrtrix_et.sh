#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=16,walltime=24:00:00,mem=20480mb
#PBS -N mrtrix_${SUBJ}_mrtrix1000_brainmask_full_ensemble
#PBS -M bacaron@imail.iu.edu


## Ensemble tractography cadidate fascicles generation script.
## 
## This shell script uses mrtrix/0.2.12 to run a series of tractography methods using both probabilistic
## and deterministic tractography based on the tensor model or on constrained spherical deconvolution. 
##
## Brent McPherson and Franco Pestilli Indiana University 2016; adapted for use in microstructure in concussion-prone athletics study by Brad Caron (IU Graduate Student, 2016)

## load necessary modules on the Karst cluster environment
module unload mrtrix/0.3.12
module load mrtrix/0.2.12

## The script requires a single input tha tis the folder name for the subject that needs to be processed
## The corrent version of the script handles only data on the local Indiana University Cluster Systems
## Under the project lifebid. There are two current data sets there the HCP and STN96
SUBJ="${SUBJ}"

## Set paths to diffusion data directories

## Concussion_Prone dataset
DWIFILENAME=data_b1000_aligned_trilin_noMEC
TOPDIR=/N/dc2/projects/lifebid/Concussion/concussion_real/${SUBJ}/diffusion_data
ANATDIR=$TOPDIR/1000/Anatomy
OUTDIR=$TOPDIR/1000/fibers_new/fibers_new_brainmask_full_ensemble_test

mkdir $OUTDIR
mkdir $ANATDIR
cp -v $TOPDIR/wm_mask.nii.gz $ANATDIR/

## Number of fibers requested and max number attempted to hit the number.
NUMFIBERS=60000
MAXNUMFIBERSATTEMPTED=1000000

##
echo 
echo Performing preprocessing of data before starting tracking...
echo 
##

## convert wm mask
mrconvert $ANATDIR/wm_mask.nii.gz $OUTDIR/${DWIFILENAME2}_wm.mif

## convert dwi's 
mrconvert $TOPDIR/1000/dt6/${DWIFILENAME2}.nii.gz $OUTDIR/${DWIFILENAME2}_dwi.mif

## make mask from DWI data
mrconvert $TOPDIR/1000/dt6/dti64trilin/bin/brainMask.nii.gz $OUTDIR/${DWIFILENAME2}_brainmask.mif

## fit tensors
cp -v $TOPDIR/fibers/data_b1000.b $OUTDIR/${DWIFILENAME2}.b
dwi2tensor $OUTDIR/${DWIFILENAME2}_dwi.mif -grad $OUTDIR/${DWIFILENAME2}.b $OUTDIR/${DWIFILENAME2}_dt.mif 

## create FA image
tensor2FA $OUTDIR/${DWIFILENAME2}_dt.mif - | mrmult - $OUTDIR/${DWIFILENAME2}_fa.mif

## create eigenvector map
tensor2vector $OUTDIR/${DWIFILENAME2}_dt.mif - | mrmult - $OUTDIR/${DWIFILENAME2}_fa.mif $OUTDIR/${DWIFILENAME2}_ev.mif

## # Estimate deconvolution kernel: Estimate the kernel for deconvolution, using voxels with highest FA
## erodes brainmask - removes extreme artifacts (w/ high FA), creates FA image, AND single fiber mask
erode $OUTDIR/${DWIFILENAME2}_brainmask.mif -npass 3 - | mrmult $OUTDIR/${DWIFILENAME2}_fa.mif - - | threshold - -abs 0.7 $OUTDIR/${DWIFILENAME2}_sf_thresh${i_thresh2}.mif

## estimates response function
for i_lmax in 2 4 6 8 10 12
	do
		estimate_response $OUTDIR/${DWIFILENAME2}_dwi.mif $OUTDIR/${DWIFILENAME2}_sf_thresh${i_thresh2}.mif -lmax $i_lmax -grad $OUTDIR/${DWIFILENAME2}.b $OUTDIR/${DWIFILENAME2}_response_lmax${i_lmax2}.txt
	done

## # End estimation of deconvolution kernel

## Perform CSD in each white matter voxel
for i_lmax in 2 4 6 8 10 12
	do
		csdeconv $OUTDIR/${DWIFILENAME2}_dwi.mif -grad $OUTDIR/${DWIFILENAME2}.b $OUTDIR/${DWIFILENAME2}_response_lmax${i_lmax2}.txt -lmax $i_lmax -mask $OUTDIR/${DWIFILENAME2}_brainmask.mif $OUTDIR/${DWIFILENAME2}_lmax${i_lmax2}.mif
		## echo DONE Lmax=$i_lmax 
	done 

##
echo DONE performing preprocessing of data before starting tracking...
##

##
echo START tracking...
##
##echo tracking Deterministic Tensorbased

streamtrack DT_STREAM $OUTDIR/${DWIFILENAME2}_dwi.mif $OUTDIR/${DWIFILENAME2}_wm_tensor-$NUMFIBERS.tck -seed $OUTDIR/${DWIFILENAME2}_wm.mif -mask $OUTDIR/${DWIFILENAME2}_wm.mif -grad $OUTDIR/${DWIFILENAME2}.b -number $NUMFIBERS -maxnum $MAXNUMFIBERSATTEMPTED


## loop over tracking and lmax
for i_tracktype in SD_STREAM SD_PROB
	do
		##
		##echo Tracking $i_tracktype Deterministic=1 Probabilistic=2 CSD-based
		##
		for i_lmax in 2 4 6 8 10 12
			do
				##echo Tracking CSD-based Lmax=$i_lmax
				streamtrack $i_tracktype $OUTDIR/${DWIFILENAME2}_lmax${i_lmax2}.mif   $OUTDIR/${DWIFILENAME2}_csd_lmax${i_lmax2}_wm_${i_tracktype2}-$NUMFIBERS.tck -seed $OUTDIR/${DWIFILENAME2}_wm.mif  -mask $OUTDIR/${DWIFILENAME2}_wm.mif  -grad $OUTDIR/${DWIFILENAME2}.b -number $NUMFIBERS -maxnum $MAXNUMFIBERSATTEMPTED
			done
		##
		##echo DONE Tracking $i_tracktype Deterministic=1 Probabilistic=2 CSD-based
		##
	done
##
echo DONE tracking. Exiting Ensemble Tracking Candidate Fascicle Generation Script
##
