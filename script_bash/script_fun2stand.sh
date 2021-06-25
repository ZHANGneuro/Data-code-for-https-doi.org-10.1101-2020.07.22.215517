#!/bin/bash


# curr_ith=3
# curr_cond="left"
# duration="4s"
# date="20201225"

# curr_ith=4
# curr_cond="right"
# duration="4s"
# date="20201225"

# curr_ith=5
# curr_cond="back"
# duration="4s"
# date="20201225"


# curr_ith=3
# curr_cond="left"
# duration="2s"
# date="20201228"

# curr_ith=4
# curr_cond="right"
# duration="2s"
# date="20201228"
#
# curr_ith=5
# curr_cond="back"
# duration="2s"
# date="20201228"

# for sub in {8..26}
# do
#   for sess in {1..4}   activation_direction_s4.feat
#   do
#
#     dd="/Users/bo/Documents/data_yuji_lab/data_fmri/NAYA"$sub"/contrast_lrb_"$duration"_"$date"_s"$sess".feat/stats/pe"$curr_ith".nii.gz"
#     p1="flirt -in $dd "
#
#     p2=" -ref /usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz "
#
#     p22="-out /Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20201229/beta_"$curr_cond"_"$duration"/mni_sub"$sub"_s"$sess".nii"
#
#     p3=" -init /Users/bo/Documents/data_yuji_lab/data_fmri/collection_fun2stand_reg/fun2stand_sub"$sub"_s"$sess".mat -applyxfm"
#
#     $p1$p2$p22$p3
#
#     # echo $p1$p2$p22$p3
#
#   done
# done



# 20210206 17-b  18-l  19-r 20-b_tn  21-l_tn  22-r_tn
# 20210207 18-b  19-l  20-r
curr_ith=5
curr_cond="back"

for sub in {8..26}
do
  for sess in {1..4}
  do

    dd="/Users/bo/Documents/data_yuji_lab/data_fmri/NAYA"$sub"/contrast_lrb_4s_20210217_s"$sess".feat/stats/pe"$curr_ith".nii.gz"
    p1="flirt -in $dd "

    p2=" -ref /usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz "

    p22="-out /Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_direction_20200217/beta_"$curr_cond"_4s/mni_sub"$sub"_s"$sess".nii"

    p3=" -init /Users/bo/Documents/data_yuji_lab/data_fmri/collection_fun2stand_reg/fun2stand_sub"$sub"_s"$sess".mat -applyxfm"

    $p1$p2$p22$p3

    # echo $p1$p2$p22$p3

  done
done
