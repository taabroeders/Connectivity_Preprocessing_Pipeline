#!/bin/bash

#SBATCH --job-name=tractography       #a convenient name for your job
#SBATCH --mem=12G                     #max memory per node
#SBATCH --partition=luna-cpu-long     #using luna short queue
#SBATCH --cpus-per-task=16            #max CPU cores per process
#SBATCH --time=12:00:00               #time limit (H:MM:SS)
#SBATCH --qos=anw-cpu                 #use anw-cpu's
#SBATCH --output=logs/slurm-%x.%j.out

#======================================================================
#                             TRACTOGRAPHY
#======================================================================

#@author: Tommy Broeders
#@email:  t.broeders@amsterdamumc.nl
#updated: 30 11 2023
#status: still being developed
#to-do: [check] removed conda environment as it did not seem needed, is this correct?

#Review History
#Reviewed by -

# Description:
# - This code is part of the "KNW-Connect Processing Pipeline".
#   It will perform tractography on the diffusion weighted data.
#
# - Prerequisites: Freesurfer and FSL tools enabled
# - Input: Subject (+session) ID
# - Output: Cortical and subcortical GM/WM segmentations
#----------------------------------------------------------------------

#Input variables
INPUT_DIR=$(realpath $1)
FULLID_file=$2
FULLID_folder=$3
FILEDIR=$4/files

dwi_nii=${INPUT_DIR}/dwi/${FULLID_file}*_dwi.nii.gz
dwi_bval=${dwi_nii%%.nii.gz}.bval

#Check if script has already been completed
[ -f dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_sift2weights.csv ] && exit 0

#Create output folder
mkdir -p dwi/${FULLID_folder}/tractography/tracts &&\

# Open fd 3 to a trace file 
exec 3> dwi/${FULLID_folder}/tractography/tracts/code_trace.txt
BASH_XTRACEFD=3
set -x #enable tracing

#Print the ID of the subject (& session if available)
printf "####$(echo ${FULLID_folder} | sed 's|/|: |')####\n$(date)\n\n"

#Convert to .mif file (ss3t_csd_beta1 can't handle .nii.gz)
mrconvert  dwi/${FULLID_folder}/preprocessing/${FULLID_file}_preprocessed_dwi.nii.gz \
           dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_preprocessed_dwi.mif \
           -fslgrad \
           dwi/${FULLID_folder}/preprocessing/${FULLID_file}_eddy_unwarped_dwi.eddy_rotated_bvecs \
           ${dwi_bval} &&\

#Estimate response function(s) for spherical deconvolution
dwi2response dhollander dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_preprocessed_dwi.mif \
             dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_response_wm.txt \
             dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_response_gm.txt \
             dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_response_csf.txt \
             -mask dwi/${FULLID_folder}/preprocessing/${FULLID_file}_final_b0_brain_mask.nii.gz \
             -voxels dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_voxels.mif &&\

#Single-shell 3-tissue constrianed spherical deconvolution
export PATH=$PATH:${FILEDIR}/MRtrix3Tissue/bin &&\
ss3t_csd_beta1 dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_preprocessed_dwi.mif \
               dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_response_wm.txt \
               dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_response_wmFOD.mif \
               dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_response_gm.txt \
               dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_response_gmFOD.mif \
               dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_response_csf.txt \
               dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_response_csfFOD.mif \
               -mask dwi/${FULLID_folder}/preprocessing/${FULLID_file}_final_b0_brain_mask.nii.gz &&\

mtnormalise dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_response_wmFOD.mif \
            dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_response_wmFOD_norm.mif \
            dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_response_gmFOD.mif \
            dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_response_gmFOD_norm.mif \
            dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_response_csfFOD.mif \
            dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_response_csfFOD_norm.mif \
            -mask dwi/${FULLID_folder}/preprocessing/${FULLID_file}_final_b0_brain_mask.nii.gz &&\

#Remove .mif file to save storage space
rm dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_preprocessed_dwi.mif &&\

# Perform tractograhy and tract filtering
printf "Performing tractography and SIFT2 filtering...\n\n" &&\
tckgen -act dwi/${FULLID_folder}/anat2dwi/hsvs_5tt/${FULLID_file}_5tthsvs_dwi_lesions.nii.gz  \
       -backtrack -select 10000000 \
       -seed_gmwmi dwi/${FULLID_folder}/anat2dwi/hsvs_5tt/${FULLID_file}_gmwmi_dwi.nii.gz \
       -maxlength 250 \
       -cutoff 0.06 \
       dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_response_wmFOD_norm.mif \
       dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_tracks_10M_ifod2.tck &&\

# sift filtering
tcksift -act dwi/${FULLID_folder}/anat2dwi/hsvs_5tt/${FULLID_file}_5tthsvs_dwi_lesions.nii.gz  \
        dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_tracks_10M_ifod2.tck \
        dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_response_wmFOD_norm.mif \
        dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_tracks_10M_ifod2_sift.tck \
        -term_number 1000000 &&\

# sift2 filtering
tcksift2 -act dwi/${FULLID_folder}/anat2dwi/hsvs_5tt/${FULLID_file}_5tthsvs_dwi_lesions.nii.gz  \
         dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_tracks_10M_ifod2.tck \
         dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_response_wmFOD_norm.mif \
         dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_sift2weights.csv \
         -out_mu dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_mu.txt \
         -out_coeffs dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_coeffs.txt &&\

# generate tck file with less streamlines (for visualization purposes)
tckedit dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_tracks_10M_ifod2_sift.tck \
        -number 200k \
        dwi/${FULLID_folder}/tractography/tracts/${FULLID_file}_tracks_200k_ifod2.tck &&\

printf "\n\n$(date)\n#### Done! ####\n"

#----------------------------------------------------------------------
#                       References, links, others, ...   
#----------------------------------------------------------------------
# Tournier, J.-D.; Calamante, F. & Connelly, A.
# Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions.
# Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670

# Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A.
# SIFT2: Enabling dense quantitative assessment of brain white matter connectivity using streamlines tractography.
# NeuroImage, 2015, 119, 338-351
