#!/bin/bash

#SBATCH --job-name=func_sdc           #a convenient name for your job
#SBATCH --mem=4G                      #max memory per node
#SBATCH --partition=luna-cpu-short    #using luna short queue
#SBATCH --cpus-per-task=2      	      #max CPU cores per process
#SBATCH --time=1:00:00                #time limit (H:MM:SS)
#SBATCH --qos=anw-cpu                 #use anw-cpu's
#SBATCH --output=logs/slurm-%x.%j.out

#======================================================================
#                   Functional distortion correction
#======================================================================

#@author: Tommy Broeders
#@email:  t.broeders@amsterdamumc.nl
#updated: 30 11 2023
#status: still being developed
#to do:

#Review History
#Reviewed by -

# Description:
# - This code is part of the "KNW-Connect Processing Pipeline".
#   It performs functional preprocessing using FEAT.
#
# - Prerequisites: FSL tools enabled
# - Input: Neck-clipped and skull-stripped T1-scan, resting-state functional MRI,
#          pipeline-folder, number of dummy scans to remove and subject (+session) ID
# - Output: Preprocessed functional MRI scans, excluding temporal filtering.
#----------------------------------------------------------------------

#Input variables
anatomical_brain=$1
restingstate=$2
subfolder=$3/files
FULLID_file=$4
FULLID_folder=$5
remove_vols=$6
fmri_oppositePE=$7
outputdir=${PWD}/func/${FULLID_folder}

#Check if script has already been completed
[ -f ${outputdir}/DisCo/output/BOLD_u.nii.gz ] && exit 0

mkdir -p func/${FULLID_folder}/DisCo/input &&\
mkdir -p func/${FULLID_folder}/DisCo/output &&\

# Open fd 3 to a trace file 
exec 3> func/${FULLID_folder}/DisCo/code_trace.txt
BASH_XTRACEFD=3
set -x #enable tracing

#Print the ID of the subject (& session if available)
printf "####$(echo ${FULLID_folder} | sed 's|/|: |')####\n$(date)\n\n"

#----------------------------------------------------------------------
#                         Opposite PE
#----------------------------------------------------------------------
if [ ! -z ${fmri_oppositePE} ]; then

cp ${restingstate} func/${FULLID_folder}/DisCo/input/${FULLID_file}_bold.nii.gz
cp ${fmri_oppositePE} func/${FULLID_folder}/DisCo/input/${FULLID_file}_bold_oppPE.nii.gz

fslmaths func/${FULLID_folder}/DisCo/input/${FULLID_file}_bold.nii.gz -Tmean \
         func/${FULLID_folder}/DisCo/output/${FULLID_file}_3Dmean_bold.nii.gz &&\

mcflirt -in func/${FULLID_folder}/DisCo/input/${FULLID_file}_bold.nii.gz \
        -out func/${FULLID_folder}/DisCo/output/${FULLID_file}_realigned_bold.nii.gz \
        -mats -plots -reffile func/${FULLID_folder}/DisCo/output/${FULLID_file}_3Dmean_bold.nii.gz -rmsrel &&\

mcflirt -in func/${FULLID_folder}/DisCo/input/${FULLID_file}_bold_oppPE.nii.gz \
        -out func/${FULLID_folder}/DisCo/output/${FULLID_file}_realigned_bold_oppPE.nii.gz \
        -mats -plots -reffile func/${FULLID_folder}/DisCo/output/${FULLID_file}_3Dmean_bold.nii.gz &&\

topupvolNUM=$(fslval ${func/${FULLID_folder}/DisCo/input/${FULLID_file}_bold_oppPE.nii.gz} dim4)

fslroi func/${FULLID_folder}/DisCo/output/${FULLID_file}_realigned_bold.nii.gz \
       func/${FULLID_folder}/DisCo/input/${FULLID_file}_realigned_bold_subset.nii.gz \
       0 ${topupvolNUM} &&\

fslmerge -t func/${FULLID_folder}/DisCo/input/${FULLID_file}_bold_combPE.nii.gz \
        func/${FULLID_folder}/DisCo/input/${FULLID_file}_realigned_bold_subset.nii.gz \
        func/${FULLID_folder}/DisCo/input/${FULLID_file}_realigned_bold_oppPE.nii.gz &&\

#Determine Phase-Encoding Direction and Readout Time and create acquisition parameters file
restingstate_json=$(realpath ${restingstate_json%%.nii.gz}.json)
PE=$(cat ${fmri_json} | grep '"PhaseEncodingDirection"' | awk -F" " '{print $2}' | sed 's/"//g' | sed 's/,//g')
RT=$(cat ${fmri_json} | grep '"*TotalReadoutTime"' | awk -F" " '{print $2}' | sed 's/"//g' | sed 's/,//g')

if [ ${PE} == "i" ];then PE_FSL="1 0 0"; PE_FSL_opp="-1 0 0"
elif [ ${PE} == "i-" ];then PE_FSL="-1 0 0"; PE_FSL_opp="1 0 0"
elif [ ${PE} == "j" ];then PE_FSL="0 1 0"; PE_FSL_opp="0 -1 0"
elif [ ${PE} == "j-" ];then PE_FSL="0 -1 0"; PE_FSL_opp="0 1 0"
elif [ ${PE} == "k" ];then PE_FSL="0 0 1"; PE_FSL_opp="0 0 -1"
elif [ ${PE} == "k-" ];then PE_FSL="0 0 -1"; PE_FSL_opp="0 0 1"
fi

for ii in {1..${topupvolNUM}};do
echo ${PE_FSL} ${RT} >> func/${FULLID_folder}/DisCo/input/${FULLID_file}_acqparams.txt
done

for ii in {1..${topupvolNUM}};do
echo ${PE_FSL_opp} ${RT} >> func/${FULLID_folder}/DisCo/input/${FULLID_file}_acqparams.txt
done

topup --imain=func/${FULLID_folder}/DisCo/input/${FULLID_file}_bold_combPE.nii.gz \
      --datain=func/${FULLID_folder}/DisCo/input/${FULLID_file}_acqparams.txt \
      --config=b02b0.cnf --out=func/${FULLID_folder}/DisCo/output/${FULLID_file}_bold_combPE_topup.nii.gz &&\

applytopup --imain=func/${FULLID_folder}/DisCo/input/${FULLID_file}_realigned_bold.nii.gz --inindex=1 \
           --datain=func/${FULLID_folder}/DisCo/input/${FULLID_file}_acqparams.txt \
           --topup=func/${FULLID_folder}/DisCo/output/${FULLID_file}_bold_combPE_topup \
           --method=jac --out=func/${FULLID_folder}/DisCo/output/${FULLID_file}_bold_DisCo.nii.gz || exit 1

#----------------------------------------------------------------------
#                         SynBOLD-DisCo
#----------------------------------------------------------------------

else

#Use raw BOLD
cp ${restingstate} func/${FULLID_folder}/DisCo/input/BOLD_d.nii.gz

#Use brain-extracted T1
cp ${anatomical_brain} func/${FULLID_folder}/DisCo/input/T1.nii.gz

#Run SynBOLD-DisCo for fieldmap-lesss distortion correction
singularity run -e -B func/${FULLID_folder}/DisCo:/tmp \
            -B func/${FULLID_folder}/DisCo/input:/INPUTS \
            -B func/${FULLID_folder}/DisCo/output:/OUTPUTS \
            -B ${FREESURFER_HOME}/license.txt:/opt/freesurfer/license.txt \
             ${subfolder}/singularity/synbold-disco.sif \
             --skull_stripped  || exit 1

fi

printf "\n\n$(date)\n#### Done! ####\n"

#----------------------------------------------------------------------
#                       References, links, others, ...   
#----------------------------------------------------------------------
# Yu, T., Cai, L. Y., Morgan, V. L., Goodale, S. E., Englot, D. J., Chang, C. E., ...
# & Schilling, K. G. (2022). SynBOLD-DisCo: Synthetic BOLD images for distortion 
# correction of fMRI without additional calibration scans. bioRxiv, 2022-09.
