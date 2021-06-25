

# $MNE_ROOT/bin/mne_ctf2fiff --ds /Users/boo/Desktop/MEG_data_script/20190513/S04_G14PKU_20190513_01.ds --fif /Users/boo/Desktop/MEG_data_script/20190513/S04_0513_01.fif
import os
import mne
import numpy
import numpy as np
import matplotlib.pyplot as plt
import os.path as op
from operator import itemgetter
from mne.io import Raw
from mne.io import read_raw_ctf
from mne.preprocessing import ICA
from mne.viz import plot_evoked_topo
from mne.minimum_norm import apply_inverse
from mne.io import read_raw_fif

# reset -f
for ith_sub in list(range(5, 14)):

    if ith_sub == 5:
        # sub05
        data_path_s1 = '/Users/boo/Desktop/MEG_data_script/sub5/S05_G14PKU_20190806_01.ds'
        data_path_s2 = '/Users/boo/Desktop/MEG_data_script/sub5/S05_G14PKU_20190806_02.ds'
        data_path_s3 = '/Users/boo/Desktop/MEG_data_script/sub5/S05_G14PKU_20190806_03.ds'
        data_path_s4 = '/Users/boo/Desktop/MEG_data_script/sub5/S05_G14PKU_20190806_04.ds'
    elif ith_sub == 6:
        # sub06
        data_path_s1 = '/Users/boo/Desktop/MEG_data_script/sub6/S06_G14PKU_20190805_01.ds'
        data_path_s2 = '/Users/boo/Desktop/MEG_data_script/sub6/S06_G14PKU_20190805_02.ds'
        data_path_s3 = '/Users/boo/Desktop/MEG_data_script/sub6/S06_G14PKU_20190805_03.ds'
        data_path_s4 = '/Users/boo/Desktop/MEG_data_script/sub6/S06_G14PKU_20190805_04.ds'
    elif ith_sub == 7:
        # sub07
        data_path_s1 = '/Users/boo/Desktop/MEG_data_script/sub7/S07_G14PKU_20190805_01.ds'
        data_path_s2 = '/Users/boo/Desktop/MEG_data_script/sub7/S07_G14PKU_20190805_02.ds'
        data_path_s3 = '/Users/boo/Desktop/MEG_data_script/sub7/S07_G14PKU_20190805_03.ds'
        data_path_s4 = '/Users/boo/Desktop/MEG_data_script/sub7/S07_G14PKU_20190805_04.ds'
    elif ith_sub == 8:
        # sub08
        data_path_s1 = '/Users/boo/Desktop/MEG_data_script/sub8/S08_G14PKU_20190801_01.ds'
        data_path_s2 = '/Users/boo/Desktop/MEG_data_script/sub8/S08_G14PKU_20190801_02.ds'
        data_path_s3 = '/Users/boo/Desktop/MEG_data_script/sub8/S08_G14PKU_20190801_03.ds'
        data_path_s4 = '/Users/boo/Desktop/MEG_data_script/sub8/S08_G14PKU_20190801_04.ds'
    elif ith_sub == 9:
        # sub09
        data_path_s1 = '/Users/boo/Desktop/MEG_data_script/sub9/S09_G14PKU_20190801_02.ds'
        data_path_s2 = '/Users/boo/Desktop/MEG_data_script/sub9/S09_G14PKU_20190801_03.ds'
        data_path_s3 = '/Users/boo/Desktop/MEG_data_script/sub9/S09_G14PKU_20190801_04.ds'
        data_path_s4 = '/Users/boo/Desktop/MEG_data_script/sub9/S09_G14PKU_20190801_05.ds'
    elif ith_sub == 10:
        # sub10
        data_path_s1 = '/Users/boo/Desktop/MEG_data_script/sub10/S10_G14PKU_20190725_01.ds'
        data_path_s2 = '/Users/boo/Desktop/MEG_data_script/sub10/S10_G14PKU_20190725_02.ds'
        data_path_s3 = '/Users/boo/Desktop/MEG_data_script/sub10/S10_G14PKU_20190725_03.ds'
        data_path_s4 = '/Users/boo/Desktop/MEG_data_script/sub10/S10_G14PKU_20190725_04.ds'
    elif ith_sub == 11:
        # sub11
        data_path_s1 = '/Users/boo/Desktop/MEG_data_script/sub11/S11_G14PKU_20190729_01.ds'
        data_path_s2 = '/Users/boo/Desktop/MEG_data_script/sub11/S11_G14PKU_20190729_02.ds'
        data_path_s3 = '/Users/boo/Desktop/MEG_data_script/sub11/S11_G14PKU_20190729_03.ds'
        data_path_s4 = '/Users/boo/Desktop/MEG_data_script/sub11/S11_G14PKU_20190729_04.ds'
    elif ith_sub == 12:
        # sub12
        data_path_s1 = '/Users/boo/Desktop/MEG_data_script/sub12/S12_G14PKU_20190806_05.ds'
        data_path_s2 = '/Users/boo/Desktop/MEG_data_script/sub12/S12_G14PKU_20190806_06.ds'
        data_path_s3 = '/Users/boo/Desktop/MEG_data_script/sub12/S12_G14PKU_20190806_07.ds'
        data_path_s4 = '/Users/boo/Desktop/MEG_data_script/sub12/S12_G14PKU_20190806_08.ds'
    elif ith_sub == 13:
        # sub13
        data_path_s1 = '/Users/boo/Desktop/MEG_data_script/sub13/S13_G14PKU_20190807_01.ds'
        data_path_s2 = '/Users/boo/Desktop/MEG_data_script/sub13/S13_G14PKU_20190807_03.ds'
        data_path_s3 = '/Users/boo/Desktop/MEG_data_script/sub13/S13_G14PKU_20190807_04.ds'
        data_path_s4 = '/Users/boo/Desktop/MEG_data_script/sub13/S13_G14PKU_20190807_05.ds'

    raw_s1 = read_raw_ctf(data_path_s1, preload=True)
    raw_s2 = read_raw_ctf(data_path_s2, preload=True)
    raw_s3 = read_raw_ctf(data_path_s3, preload=True)
    raw_s4 = read_raw_ctf(data_path_s4, preload=True)
    raw = mne.io.concatenate_raws([raw_s1, raw_s2, raw_s3, raw_s4])
    raw.filter(1, 100, fir_design='firwin')
    events = mne.find_events(raw, stim_channel='UPPT001', shortest_event=1)
    raw_resampled, events_resampled = raw.copy().resample(200, npad='auto', events=events)
    numpy.save("/Users/boo/Desktop/MEG_data_script/PreProcessed_data/events_post_resample_sub" + str(ith_sub) +"_100hz_sfre_200.npy", events_resampled)



# ica
reject = dict(mag=5e-12, grad=4000e-13)
ica = ICA(n_components=25, method='fastica', fit_params=None, random_state=0)
ica.fit(raw_resampled, reject=reject)  # picks=picks,
# ica.plot_components(title="ica")
ica.plot_sources(raw_resampled, title='', start=0, stop=10)  # exclude=[]

############################
# ica records 200hz
############################
ica.exclude += [0, 2, 8, 12]  #5
ica.exclude += [0, 1, 2]  #6
ica.exclude += [0, 5, 18]  #7
ica.exclude += [0, 1, 2, 3, 10, 13, 19]  # 8
ica.exclude += [0, 1, 2, 12, 14]  # 9
ica.exclude += [0, 1, 24]  #10
ica.exclude += [0, 1, 3, 4, 5, 11, 12, 17]  #11
ica.exclude += [0, 6, 10]  # 12
ica.exclude += [0, 1, 2, 12, 24]  # sub13

# raw_resampled.plot(n_channels=20)
ica.apply(raw_resampled)
raw_resampled.plot(n_channels=20, duration=50)
# raw_resampled.info['bads'].extend(['MLC17-4503', 'MLC25-4503'])
raw_resampled.save('/Users/boo/Desktop/MEG_data_script/PreProcessed_data/artefact_removed_sub'+str(ith_sub)+'_raw_100hz_sfre_200.fif', overwrite=True)





############################
# ica records 100hz
############################
ica.exclude += [1, 2, 3, 10, 12, 24, ]  #5
ica.exclude += [0, 1, 2, 24]  # sub6
ica.exclude += [0, 5, 22]  # sub7 new
ica.exclude += [0, 1, 2, 3, 5, 18, 19, 21, 24]  # sub8
ica.exclude += [0, 1, 2, 6, 10, 11, 12, 13, 14, 17, 22]  # sub9
ica.exclude += [0, 1, 6, 15]  # sub10
ica.exclude += [0, 1, 3, 4, 10, 13, 14, 19, 23, 24]  # sub11
ica.exclude += [0, 2, 6, 11, 12, 13, 24]  # sub12 new
ica.exclude += [0, 1, 2, 5, 10, 15, 23]  # sub13


