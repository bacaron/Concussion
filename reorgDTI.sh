#!/bin/sh

## Script for averaging and reorganizing b0 values; used in batch_diffusion_clean script.  Developed by Hu Cheng, used by Brad Caron (IU Graduate Student, 2016) for microstructure in concussion-prone athletics study

fslsplit $1
fslmerge -t nodif vol0000 vol0009 vol0018 vol0027 vol0036 vol0045 vol0054 vol0063 vol0080 vol0089 vol0098 vol0107 vol0116 vol0125 vol0134;
rm vol0000.nii.gz vol0009.nii.gz vol0018.nii.gz vol0027.nii.gz vol0036.nii.gz vol0045.nii.gz vol0054.nii.gz
rm vol0063.nii.gz vol0080.nii.gz vol0089.nii.gz vol0098.nii.gz vol0107.nii.gz vol0116.nii.gz vol0125.nii.gz vol0134.nii.gz
fslmerge -t dif vol*
rm vol*
fslmerge -t reorgdata nodif dif 