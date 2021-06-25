## Data-code-for-https-doi.org-10.1101-2020.07.22.215517
<br />
This repository includes behavioral & neuroimaging analyzing codes (MRI/MEG) and data (MRI & MEG behavior) for the manuscript entitled "Distinct networks coupled with parietal cortex for spatial representations inside and outside the visual field" <br />
doi: https://doi.org/10.1101/2020.07.22.215517  
<br /><br />
For MRI & MEG imaging datasets, and eye datasets, please contact Dr. Yuji Naya in the following address:
<br />

``` diff
Yuji Naya
Email: yujin@pku.edu.cn
McGovern Institute for Brain Research, Peking University
No. 52, Haidian Road, Wang Kezhen Building, Room 1707
Haidian District, Beijing 100805, China
Telephone: +86-10-62765734
```


<br /><br />
## Repository structure:
Note: for more detailed code description or questions please contact bo.zhang@pku.edu.cn
<br />
.<br />
|-- beh_eyedata_meg <br />
│&emsp;&emsp;&emsp;&emsp; |-- sub_x_sx_rawdata.txt &emsp;  ``col indicates `sub_no` `exp_cond` `map_id` `enter_dir_id` `fac_cha_id` `tar_cha_id` `hn1_id` `hn2_id` `hn3_id` `prof_id` `ego_dir` `cue` `response` `` <br />
│&emsp;&emsp;&emsp;&emsp;<br />
│&emsp;&emsp;&emsp;&emsp; |-- sub_x_sx_timing.txt &emsp; ``col indicates onsets of `ITI` `session id` `noise screen` `facing period` `noise screen` `targeting period` `noise screen` `cue` `response` `` <br />
│&emsp;&emsp;&emsp;&emsp;<br />
|-- beh_eyedata_fmri <br />
│&emsp;&emsp;&emsp;&emsp; |-- sub_x_formal_rawdata.txt &emsp; ``col indicates `sub_no` `map_id` `enter_dir_id` `head nodding(HD)` `HD response` `HD outcome` `score` `fac_cha_id` `fac_dir_id` `target_dir` `tar_cha_id` `cue` `response` `` <br />
│&emsp;&emsp;&emsp;&emsp;<br />
│&emsp;&emsp;&emsp;&emsp; |-- sub_x_formal_Time_record_t.txt &emsp;  ``col indicates onsets of `ITI` `session id` `noise screen` `facing period` `noise screen` `targeting period` `noise screen` `cue` `response`  `` <br />
│&emsp;&emsp;&emsp;&emsp;<br />
|-- script_R <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_fmri_beh_format_reorganize.R &emsp; ``reorganize fMRI behavioral data format for further process`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_meg_beh_format_reorganize.R &emsp;  ``reorganize MEG behavioral data format for further process`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_plot_eyedata_number_meg.R &emsp;  ``plot number of fixation and saccade for MEG data`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_plot_eyedata_number_mri.R &emsp;  ``plot number of fixation and saccade for fMRI data`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_plot_eyedata_pos_barchart_meg.R &emsp;  ``plot position of fixation and saccade for MEG data`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_plot_eyedata_pos_barchart_mri.R &emsp;  ``plot position of fixation and saccade for fMRI data`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_plot_meg_beh.R &emsp;  ``plot behavioral performance for MEG and fMRI data`` <br />
│&emsp;&emsp;&emsp;&emsp;<br />
|-- script_bash <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_bandpass_filtering.sh &emsp; ``perform bandpass filtering on 4d bold signal series`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_fsf_generator.sh &emsp;  ``generate FSL fsf file for 1st level GLM analysis`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_fun2stand.sh &emsp;  ``transforamtion script from native space to standard space`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_permutation.sh &emsp;  ``permutation test using FSL randomize function`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_run_feat.sh &emsp;  ``script for parallel process`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_stand2fun.sh &emsp;  ``transforamtion script from standard space to native space`` <br />
│&emsp;&emsp;&emsp;&emsp;<br />
|-- script_matlab <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_FC_correlation.m &emsp; ``perform bandpass filtering on 4d bold signal series`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_FC_extract_targeting_TRs.m &emsp;  ``generate FSL fsf file for 1st level GLM analysis`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_MRI_activity_volume_contrast.m &emsp;  ``transforamtion script from native space to standard space`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_MRI_timingFile_generator.m &emsp;  ``permutation test using FSL randomize function`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_dicom2nii.m &emsp;  ``script for parallel process`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_eyeData_number_meg.m &emsp;  ``transforamtion script from standard space to native space`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_eyeData_number_mri.m &emsp;  ``transforamtion script from standard space to native space`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_eyeData_visual_angle_meg.m &emsp;  ``transforamtion script from standard space to native space`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_eyeData_visual_angle_mri.m &emsp;  ``transforamtion script from standard space to native space`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_eyeDate_gaze_anova_meg.m &emsp;  ``transforamtion script from standard space to native space`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_eyeDate_gaze_anova_mri2.m &emsp;  ``transforamtion script from standard space to native space`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_load_behavioral_table_fmri.m &emsp;  ``transforamtion script from standard space to native space`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_load_behavioral_table_meg.m &emsp;  ``transforamtion script from standard space to native space`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_mni_r_brain_session_averaged_from_glm.m &emsp;  ``transforamtion script from standard space to native space`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_ploting_draw_colorbar.m &emsp;  ``transforamtion script from standard space to native space`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_statistic_contrast.m &emsp;  ``transforamtion script from standard space to native space`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_statistic_glm.m &emsp;  ``transforamtion script from standard space to native space`` <br />
│&emsp;&emsp;&emsp;&emsp; |-- script_statistic_mvpa.m &emsp;  ``transforamtion script from standard space to native space`` <br />


<br /><br />
## Env & Dependency:
MacOS Big Sur `11.4`<br />
Matlab `R2019b`<br />
R `4.0.2`<br />
R studio `1.3.1073`<br />
Python `3.3.8`<br />
Conda `4.9.2`<br />
MNE `v0.23.0`<br />

