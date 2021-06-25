#!/bin/bash

#/mnt/hgfs/zhang_bo/fmri_script/bash
#/Users/boo/Desktop/fmri_script/bash
fsf_list=(/Users/bo/Documents/bash_script/template_liujia_pilot_mri/fsf*)

for i in "${fsf_list[@]}"; do
  (
    feat $i
  )&
  if (( $(wc -w <<<$(jobs -p)) % 4 == 0 )); then wait; fi
done
