#!/bin/bash

#SBATCH --job-name=check_env          #a convenient name for your job
#SBATCH --mem=4M                      #max memory per node
#SBATCH --partition=luna-cpu-short    #using luna short queue
#SBATCH --cpus-per-task=1      	      #max CPU cores per process
#SBATCH --time=0:05:00                #time limit (DD-HH:MM)
#SBATCH --qos=anw-cpu                 #use anw-cpu's

#======================================================================
#                         Check Modules
#======================================================================

#@author: Tommy Broeders
#@email:  t.broeders@amsterdamumc.nl
#updated: 05 05 2023
#status: still being developed 
#to-do: Make this code generalizable to different servers

#Review History
#Reviewed by -

# Description:
# - This code is part of the "KNW-Connect Processing Pipeline".
#   It will check all the required modules.
#
# - Prerequisites: None
# - Input:
# - Output:
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#                    Modules required
#----------------------------------------------------------------------

[ -z $FSLDIR ] || [ ! -f ${FSLDIR}/bin/fsl ] || [ -z $(printf $PATH | grep $FSLDIR/bin) ] &&\
{ printf 'ERROR: FSL_DIR not correctly set as environmental variable or not added to PATH, aborting...\n\n' ; exit 1; }

[ -z $FREESURFER_HOME ] || [ ! -f ${FREESURFER_HOME}/bin/freesurfer ]  || [ -z $(printf $PATH | grep $FREESURFER_HOME/bin) ] &&\
{ printf 'ERROR: FREESURFER_HOME not correctly set as environmental variable or not added to PATH, aborting...\n\n' ; exit 1; }

[ -z $ANTSPATH ] || [ ! -f ${ANTSPATH}/ANTS ] || [ -z $(printf $PATH | grep $ANTSPATH) ]  &&\
{ printf 'ERROR: ANTSPATH not correctly set as environmental variable or not added to PATH, aborting...\n\n' ; exit 1; }

[ -z $(printf $PATH | tr "[:lower:]" "[:upper:]" | grep MRTRIX) ] || [ -z $(which mrview) ] &&\
{ printf 'ERROR: MRtrix not correctly added to PATH, aborting...\n\n' ; exit 1; }

[ -d $1/files/MRtrix3Tissue ] &&\
{ printf 'ERROR: MRtrix3Tissue not correctly set, aborting...\n\n' ; exit 1; }

[ -d $1/files/ICA-AROMA ] &&\
{ printf 'ERROR: ICA-AROMA not correctly set, aborting...\n\n' ; exit 1; }

[ -d $1/files/niftyseg ] &&\
{ printf 'ERROR: niftyseg not correctly set, aborting...\n\n' ; exit 1; }

[ -f $1/files/singularity/synb0-disco.sif ] &&\
{ printf 'ERROR: synb0-disco singularity not correctly set, aborting...\n\n' ; exit 1; }

[ -f $1/files/singularity/synbold-disco.sif ] &&\
{ printf 'ERROR: synbold-disco singularity not correctly set, aborting...\n\n' ; exit 1; }

exit 0
