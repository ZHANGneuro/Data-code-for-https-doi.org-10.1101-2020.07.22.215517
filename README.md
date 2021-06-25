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
│&emsp;&emsp;&emsp;&emsp;|-- sub_`No.`_s`No.`_rawdata.txt &emsp;\# col indicates `sub_no` `exp_cond` `map_id` `enter_dir_id` `fac_cha_id` `tar_cha_id` `hn1_id` `hn2_id` `hn3_id` `prof_id` `ego_dir` `cue` `response` <br /><br />
│
│&emsp;&emsp;&emsp;&emsp;|-- sub_\*_s\*_timing.txt &emsp;\# col indicates onsets of `ITI` `session id` `noise screen` `facing period` `noise screen` `targeting period` `noise screen` `cue` `response` <br /><br />
│
|-- beh_eyedata_fmri <br />
│&emsp;&emsp;&emsp;&emsp;|-- sub_`No.`_formal_rawdata.txt &emsp;\# col indicates `sub_no` `map_id` `enter_dir_id` `head nodding(HD)` `HD response` `HD outcome` `score` `fac_cha_id` `fac_dir_id` `target_dir` `tar_cha_id` `cue` `response` <br /><br />
│
│&emsp;&emsp;&emsp;&emsp;|-- sub_\*_s\*_timing.txt &emsp;\# col indicates onsets of `ITI` `session id` `noise screen` `facing period` `noise screen` `targeting period` `noise screen` `cue` `response` <br /><br /> 

        

<br /><br />
## Env & Dependency:
MacOS Big Sur `11.4`<br />
Matlab `R2019b`<br />
R `4.0.2`<br />
R studio `1.3.1073`<br />
Python `3.3.8`<br />
Conda `4.9.2`<br />
MNE `v0.23.0`<br />

