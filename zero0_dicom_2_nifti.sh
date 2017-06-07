#!/bin/bash

## This will convert raw dicom files to nifti.nii.gz files.  The 'topdir' should be the directory where your data is stored,
## and 'subj' should be the subject number.  The indices refer to number of the directory in which the specific sequence of
## interest (i.e. anatomical T1, PA/AP DWI files) is stored.  If images are stored in singular file, then you can remove the
## 'for files' loop from the script and run just the "for subj" loop.  Inputs are the images for the sequence of interest in
## .dcm format; outputs are .nii.gz files for each sequence.
##
## Brad Caron, Indiana University, Pestilli Lab 2017

## The code is optimized and written for Karst at IU 
## https://kb.iu.edu/d/bezu

## Load modules for Karst
module load mricron

## Hard-coded for "Effects of long-term participation in high impact sport" (Caron et al, 2017; in prep) study. Change for
## your file format
topdir="/N/dc2/projects/lifebid/Concussion/concussion/group3/"
SUBJ=$1

## Indices refer to the 
## T1 (ANATOMY, 1) 
## DWI (PA diffusion, 10) 
## DWI (AP diffusion, 11)
## Change these to match your file format.  If images are in one directory, this input can be removed
files="2 10 11"

for subj in $SUBJ
  do
    for FILES in $files
      do
        cd $topdir/$subjects/$FILES
        
        ## This operation will create BVECS/BVALS and NIFTI.  These files will be saved in same folder as the .dcm images
        dcm2nii ./1_${FILES}_dicoms
      done
  done
