# Data-code-for-https-doi.org-10.1101-2020.07.22.215517

This repository includes behavioral & neuroimaging analyzing codes (MRI/MEG) and data (MRI & MEG behavior) for the manuscript entitled "Distinct networks coupled with parietal cortex for spatial representations inside and outside the visual field", doi: https://doi.org/10.1101/2020.07.22.215517  

MRI & MEG imaging data, eye data will be avaiable upon request due to their large file size. Please contact Dr. Yuji Naya in the following address:

Yuji Naya (Corresponding author)
Email: yujin@pku.edu.cn
McGovern Institute for Brain Research, Peking University
No. 52, Haidian Road, Wang Kezhen Building, Room 1707
Haidian District, Beijing 100805, China 
Telephone: +86-10-62765734

A detailed code description will be also aviable upon request, please contact bo.zhang@pku.edu.cn



The scripts are the customized (simplified) version for representational similarity analysis, and was made for the paper published here: 



## Docs:

1. preprocess MRI image e.g. using `FSL`. 
2. obtain multi-voxel pattern by creating a 1st level univarite GLM.
3. use `ct_by_voxel_fsl.m` to compute correlation matrix across trial-pairs for each image voxel.
4. use `mvpa_analysis_fsl_glm` to compute representational similarity by create a 2nd level GLM.

Note that recoding will be needed for your case.
