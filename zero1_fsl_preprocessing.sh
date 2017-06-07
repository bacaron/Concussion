#!/bin/bash

## Script for preprocessing both phase encoding directions and combining to singular diffusion-weighted MRI nifti fiels.
##
## Modified from early code by Hu Cheng adtapted fromthe Human Connectome Pipelines. 
## 
## Brad Caron, Indiana University, Pestilli Lab 2017

## The code is optimized and written for Karst at IU 
## https://kb.iu.edu/d/bezu

## Load modules
module load fsl/5.0.9

## FIle paths
## We first identify files on the path on DC2 at IU
##  Data directory; replace with your data directory set-up.
datadir1="/N/dc2/projects/lifebid/Concussion/concussion_real";

## Subjects are coded as numbers we define a range of subjects such as 1_5
SUBJ=$1; 

## path to diffusion directory
diffdir1="$datadir1/$SUBJ/diffusion_directory/diffusion";

## We need twofiles an AP and PA files to be iused to estimate distorsions.
## Phase encoding directorion folders
diffphase="PA AP";

## Eddy-current correction
## 
## We use eddycorrect reorganization of output from eddy_current correction, combination to make a no
## -diffusion mean file, and perform brain extraction on both phase encoding directions
for phase in $diffphase
	do
		cd $diffdir1/$phase;
		eddy_correct $phase r$phase 0; ## eddy-correction
		
		## 
		. $datadir1/$SUBJ/diffusion_directory/bin/reorgDTI r$phase; ## reorganizes b0 images
		
		## 
		fslmaths nodif -Tmean ${phase}_nodif_mean; ## averages no diffusion images
		bet ${phase}_nodif_mean nodif_brain -f 0.4 -g 0 -m; ## brain extraction
done

## Merge directions, apply top-up, fit tensors
cd $diffdir1

## merging b0 images to mean no diffusion image of each direction
fslmerge -t b0_images ./PA/PA_nodif_mean ./AP/AP_nodif_mean

## setting up top-up for susceptibility correction
topup --imain=b0_images.nii.gz \
      --datain=acq_params.txt \ 
      --config=b02b0.cnf \
      --out=my_topup_results \
      --fout=my_field \
      --iout=my_unwarped_images

## applying top-up
applytopup --imain=./PA/rPA.nii.gz,./AP/rAP.nii.gz \
           --datain=acq_params.txt \
	   --inindex=1,2 \
	   --topup=my_topup_results \
	   --out=data

## brain extraction on combined data image
bet data nodif_brain -f 0.4 -g 0 -m

## fits tensor model
dtifit --data=data \
       --out=dti \
       --mask=nodif_brain_mask \
       --bvecs=bvecs \
       --bvals=bvals