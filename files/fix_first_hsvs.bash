#!/bin/bash

SEGFOLDER=$(realpath $1)
FREESURFER_DIR=$(realpath $2)
output_file=$(realpath $3)
FILEDIR=$4

echo " WARNING: HSVS issue detected. Applying in-house written fix by rerunning first..."

#change working directory
cd $SEGFOLDER

#Rerun FIRST and check if succesfully completed
run_first_all -s L_Accu,R_Accu,L_Amyg,R_Amyg,L_Caud,R_Caud,L_Hipp,R_Hipp,L_Pall,R_Pall,L_Puta,R_Puta,L_Thal,R_Thal -i T1.nii -b -o first
if [ ! -f first_all_none_firstseg.nii.gz ];then
    echo "ERROR: FIRST failed again, aborting..."
    exit 1
fi
for ii in first.logs/*.e*;do
    if [ -s $ii ];then
        echo "WARNING: FIRST completed but error logs non-empty. Check output thoroughly!";break
    fi
done

#Overwrite subsequent steps
${FILEDIR}/singularity/MRtrix3.sif meshconvert first-L_Accu_first.vtk first-L_Accu_transformed.vtk -transform first2real T1.nii -force
${FILEDIR}/singularity/MRtrix3.sif mesh2voxel first-L_Accu_transformed.vtk aparc.mif Left-Accumbens-area.mif -force
${FILEDIR}/singularity/MRtrix3.sif meshconvert first-R_Accu_first.vtk first-R_Accu_transformed.vtk -transform first2real T1.nii -force
${FILEDIR}/singularity/MRtrix3.sif mesh2voxel first-R_Accu_transformed.vtk aparc.mif Right-Accumbens-area.mif -force
${FILEDIR}/singularity/MRtrix3.sif meshconvert first-L_Amyg_first.vtk first-L_Amyg_transformed.vtk -transform first2real T1.nii -force
${FILEDIR}/singularity/MRtrix3.sif mesh2voxel first-L_Amyg_transformed.vtk aparc.mif Left-Amygdala.mif -force
${FILEDIR}/singularity/MRtrix3.sif meshconvert first-R_Amyg_first.vtk first-R_Amyg_transformed.vtk -transform first2real T1.nii -force
${FILEDIR}/singularity/MRtrix3.sif mesh2voxel first-R_Amyg_transformed.vtk aparc.mif Right-Amygdala.mif -force
${FILEDIR}/singularity/MRtrix3.sif meshconvert first-L_Caud_first.vtk first-L_Caud_transformed.vtk -transform first2real T1.nii -force
${FILEDIR}/singularity/MRtrix3.sif mesh2voxel first-L_Caud_transformed.vtk aparc.mif Left-Caudate.mif -force
${FILEDIR}/singularity/MRtrix3.sif meshconvert first-R_Caud_first.vtk first-R_Caud_transformed.vtk -transform first2real T1.nii -force
${FILEDIR}/singularity/MRtrix3.sif mesh2voxel first-R_Caud_transformed.vtk aparc.mif Right-Caudate.mif -force
${FILEDIR}/singularity/MRtrix3.sif meshconvert first-L_Hipp_first.vtk first-L_Hipp_transformed.vtk -transform first2real T1.nii -force
${FILEDIR}/singularity/MRtrix3.sif mesh2voxel first-L_Hipp_transformed.vtk aparc.mif Left-Hippocampus.mif -force
${FILEDIR}/singularity/MRtrix3.sif meshconvert first-R_Hipp_first.vtk first-R_Hipp_transformed.vtk -transform first2real T1.nii -force
${FILEDIR}/singularity/MRtrix3.sif mesh2voxel first-R_Hipp_transformed.vtk aparc.mif Right-Hippocampus.mif -force
${FILEDIR}/singularity/MRtrix3.sif meshconvert first-L_Pall_first.vtk first-L_Pall_transformed.vtk -transform first2real T1.nii -force
${FILEDIR}/singularity/MRtrix3.sif mesh2voxel first-L_Pall_transformed.vtk aparc.mif Left-Pallidum.mif -force
${FILEDIR}/singularity/MRtrix3.sif meshconvert first-R_Pall_first.vtk first-R_Pall_transformed.vtk -transform first2real T1.nii -force
${FILEDIR}/singularity/MRtrix3.sif mesh2voxel first-R_Pall_transformed.vtk aparc.mif Right-Pallidum.mif -force
${FILEDIR}/singularity/MRtrix3.sif meshconvert first-L_Puta_first.vtk first-L_Puta_transformed.vtk -transform first2real T1.nii -force
${FILEDIR}/singularity/MRtrix3.sif mesh2voxel first-L_Puta_transformed.vtk aparc.mif Left-Putamen.mif -force
${FILEDIR}/singularity/MRtrix3.sif meshconvert first-R_Puta_first.vtk first-R_Puta_transformed.vtk -transform first2real T1.nii -force
${FILEDIR}/singularity/MRtrix3.sif mesh2voxel first-R_Puta_transformed.vtk aparc.mif Right-Putamen.mif -force
${FILEDIR}/singularity/MRtrix3.sif meshconvert first-L_Thal_first.vtk first-L_Thal_transformed.vtk -transform first2real T1.nii -force
${FILEDIR}/singularity/MRtrix3.sif mesh2voxel first-L_Thal_transformed.vtk aparc.mif Left-Thalamus-Proper.mif -force
${FILEDIR}/singularity/MRtrix3.sif meshconvert first-R_Thal_first.vtk first-R_Thal_transformed.vtk -transform first2real T1.nii -force
${FILEDIR}/singularity/MRtrix3.sif mesh2voxel first-R_Thal_transformed.vtk aparc.mif Right-Thalamus-Proper.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrconvert /scratch/anw/tabroeders/test_laura/output/freesurfer/sub-PD07_ses-T1/mri/norm.mgz -datatype uint16 -stride +1,+2,+3 T1RAS_16b.nii -force
acpcdetect -i T1RAS_16b.nii
${FILEDIR}/singularity/MRtrix3.sif mrcalc aparc.mif nan -eq - | \
${FILEDIR}/singularity/MRtrix3.sif mredit - AC_FAST_mask.mif -scanner -sphere 1.8966964721680029,19.75,26.106994628906193 8 1 -force
${FILEDIR}/singularity/MRtrix3.sif mrtransform /scratch/anw/tabroeders/test_laura/output/freesurfer/sub-PD07_ses-T1/mri/norm.mgz -template aparc.mif - | \
${FILEDIR}/singularity/MRtrix3.sif mrcalc - AC_FAST_mask.mif -mult AC_T1.nii -force
fast -N AC_T1.nii
${FILEDIR}/singularity/MRtrix3.sif mrconvert AC_T1_pve_2.nii.gz AC.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrmath lh.pial.mif rh.pial.mif sum - | \
${FILEDIR}/singularity/MRtrix3.sif mrcalc - 1.0 -min tissue0_init.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrmath Left-Accumbens-area.mif Right-Accumbens-area.mif Left-Amygdala.mif Right-Amygdala.mif Left-Caudate.mif Right-Caudate.mif -force \
Left-Hippocampus.mif Right-Hippocampus.mif Left-Pallidum.mif Right-Pallidum.mif Left-Putamen.mif Right-Putamen.mif Left-Thalamus-Proper.mif Right-Thalamus-Proper.mif sum - | \
${FILEDIR}/singularity/MRtrix3.sif mrcalc - 1.0 -min tissue1_init.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrmath lh.white.mif rh.white.mif Fornix.mif combined_corpus_callosum.mif AC.mif brain_stem.mif sum - | \
${FILEDIR}/singularity/MRtrix3.sif mrcalc - 1.0 -min tissue2_init.mif -force

#Check for vessel issue 
for vessel in *vessel.mif;do
    if [ $(${FILEDIR}/singularity/MRtrix3.sif mrstats -mask ${vessel} -output count ${vessel}) -gt 0 ]; then
    bash ${FILEDIR}/fix_hsvs_issue.bash ${SEGFOLDER} ${FREESURFER_DIR} ${output_file} ${FILEDIR} ${vessel} || exit 1 && exit 0
    fi
done

#If not observed, overwrite subsequent fiels
${FILEDIR}/singularity/MRtrix3.sif mrmath Left-Inf-Lat-Vent.mif 3rd-Ventricle.mif 4th-Ventricle.mif CSF.mif Left-vessel.mif Right-Inf-Lat-Vent.mif Right-vessel.mif 5th-Ventricle.mif \
Left_LatVent_ChorPlex.mif Right_LatVent_ChorPlex.mif sum - | ${FILEDIR}/singularity/MRtrix3.sif mrcalc - 1.0 -min tissue3_init.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrmath Left-Lesion.mif Right-Lesion.mif sum - | ${FILEDIR}/singularity/MRtrix3.sif mrcalc - 1.0 -min tissue4.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrcalc tissue3_init.mif tissue3_init.mif tissue4.mif -add 1.0 -sub 0.0 -max -sub 0.0 -max tissue3.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrmath tissue3.mif tissue4.mif sum tissuesum_34.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrcalc tissue1_init.mif tissue1_init.mif tissuesum_34.mif -add 1.0 -sub 0.0 -max -sub 0.0 -max tissue1.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrmath tissue1.mif tissue3.mif tissue4.mif sum tissuesum_134.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrcalc tissue2_init.mif tissue2_init.mif tissuesum_134.mif -add 1.0 -sub 0.0 -max -sub 0.0 -max tissue2.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrmath tissue1.mif tissue2.mif tissue3.mif tissue4.mif sum tissuesum_1234.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrcalc tissue0_init.mif tissue0_init.mif tissuesum_1234.mif -add 1.0 -sub 0.0 -max -sub 0.0 -max tissue0.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrmath tissue0.mif tissue1.mif tissue2.mif tissue3.mif tissue4.mif sum tissuesum_01234.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrcalc aparc.mif 6 -eq aparc.mif 7 -eq -add aparc.mif 8 -eq -add aparc.mif 45 -eq -add aparc.mif 46 -eq -add aparc.mif 47 -eq -add Cerebellum_volume.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrcalc T1.nii Cerebellum_volume.mif -mult T1_cerebellum_precrop.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrgrid T1_cerebellum_precrop.mif crop -mask Cerebellum_volume.mif T1_cerebellum.nii -force

#This code is not affected by the new FIRST segmentation or the exclusion of the vessel, so to save runtime this is excluded
#fast -N T1_cerebellum.nii
#${FILEDIR}/singularity/MRtrix3.sif mrtransform T1_cerebellum_pve_0.nii.gz -interp nearest -template aparc.mif FAST_0.mif -force
#${FILEDIR}/singularity/MRtrix3.sif mrtransform T1_cerebellum_pve_1.nii.gz -interp nearest -template aparc.mif FAST_1.mif -force

#Overwrite subsequents steps 
${FILEDIR}/singularity/MRtrix3.sif mrtransform T1_cerebellum_pve_2.nii.gz -interp nearest -template aparc.mif FAST_2.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrcalc Cerebellum_volume.mif tissuesum_01234.mif -add 0.5 -gt 1.0 tissuesum_01234.mif -sub 0.0 -if Cerebellar_multiplier.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrconvert tissue0.mif tissue0_fast.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrcalc tissue1.mif Cerebellar_multiplier.mif FAST_1.mif -mult -add tissue1_fast.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrcalc tissue2.mif Cerebellar_multiplier.mif FAST_2.mif -mult -add tissue2_fast.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrcalc tissue3.mif Cerebellar_multiplier.mif FAST_0.mif -mult -add tissue3_fast.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrconvert tissue4.mif tissue4_fast.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrmath tissue0_fast.mif tissue1_fast.mif tissue2_fast.mif tissue3_fast.mif tissue4_fast.mif sum tissuesum_01234_fast.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrcalc 1.0 tissuesum_01234_fast.mif -sub tissuesum_01234_fast.mif 0.0 -gt ${FREESURFER_DIR}/mri/brainmask.mgz -add 1.0 -min -mult 0.0 -max csf_fill.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrcalc tissue3_fast.mif csf_fill.mif -add tissue3_fast_filled.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrcat tissue0_fast.mif tissue1_fast.mif tissue2_fast.mif tissue3_fast_filled.mif tissue4_fast.mif - -axis 3 | \
${FILEDIR}/singularity/MRtrix3.sif 5ttedit - result.mif -none brain_stem_crop.mif -force
${FILEDIR}/singularity/MRtrix3.sif mrconvert result.mif ${output_file} -force
${FILEDIR}/singularity/MRtrix3.sif 5ttcheck result.mif