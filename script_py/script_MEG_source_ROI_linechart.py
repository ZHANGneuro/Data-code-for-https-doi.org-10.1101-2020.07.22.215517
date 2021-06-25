# source contrast get averaged
# reset -f
import os
import numpy
import numpy as np
import mne
import matplotlib.pyplot as plt
from mne.io import read_raw_fif
from scipy import stats as stats
from mne.stats import permutation_t_test
from mne.stats import (spatio_temporal_cluster_1samp_test,
                       summarize_clusters_stc)
from sklearn.base import clone
from mne.connectivity import spectral_connectivity
from operator import itemgetter
from mne.minimum_norm import apply_inverse_epochs, read_inverse_operator
import re
from mne.connectivity import envelope_correlation
from mne.stats import permutation_cluster_1samp_test

# fs source space
src_fs = mne.read_source_spaces('/Users/boo/Desktop/MEG_data_script/PreProcessed_data/fsaverage-src.fif')
fsave_vertices = [s['vertno'] for s in src_fs]
stc_template = mne.read_source_estimate(
    '/Users/boo/Desktop/MEG_data_script/analysis_source_result/stc_template-rh.stc')
stc_template.subject = 'fsaverage'

# label
# label_name_list_mtl = ['Hippocampus', 'ParaHippocampal', 'Enterinal', 'Perirhinal']
label_name_list_mtl = ['MTL_bigger']

hemi_pool = ['_lh', '_rh']
label_list_path = []
for r, d, f in os.walk('/Users/boo/Desktop/MEG_data_script/aal_51/labels_frequency_each_band_fs/'):
    for ith_hemi in list(range(0, len(hemi_pool))):
        for ith_label_path in list(range(0, len(label_name_list_mtl))):
            for file in f:
                if hemi_pool[ith_hemi] in file and label_name_list_mtl[ith_label_path] in file:
                    label_list_path.append(os.path.join(r, file))

label_list = []
label_parietal = mne.read_label(
    '/Users/boo/Desktop/MEG_data_script/aal_51/labels_conn_each_band_fs/Parietal_rh.label') + mne.read_label(
    '/Users/boo/Desktop/MEG_data_script/aal_51/labels_conn_each_band_fs/Parietal_lh.label')
label_precuneus = mne.read_label(
    '/Users/boo/Desktop/MEG_data_script/aal_51/labels_conn_each_band_fs/Precuneus_rh.label') + mne.read_label(
    '/Users/boo/Desktop/MEG_data_script/aal_51/labels_conn_each_band_fs/Precuneus_lh.label')
label_Cingulate_Post = mne.read_label(
    '/Users/boo/Desktop/MEG_data_script/aal_51/labels_conn_each_band_fs/Cingulate_Post_rh.label') + mne.read_label(
    '/Users/boo/Desktop/MEG_data_script/aal_51/labels_conn_each_band_fs/Cingulate_Post_lh.label')
label_SMA = mne.read_label(
    '/Users/boo/Desktop/MEG_data_script/aal_51/labels_conn_each_band_fs/SMA_rh.label') + mne.read_label(
    '/Users/boo/Desktop/MEG_data_script/aal_51/labels_conn_each_band_fs/SMA_lh.label')
label_FEF = mne.read_label(
    '/Users/boo/Desktop/MEG_data_script/aal_51/labels_conn_each_band_fs/FEF_rh.label') + mne.read_label(
    '/Users/boo/Desktop/MEG_data_script/aal_51/labels_conn_each_band_fs/FEF_lh.label')
label_list.append(label_parietal)
label_list.append(label_precuneus)
label_list.append(label_Cingulate_Post)
label_list.append(label_SMA)
label_list.append(label_FEF)
for ith_label in list(range(0, len(label_list_path))):
    label_list.append(mne.read_label(label_list_path[ith_label]))

yaxis_label_list = ['Parietal', 'Precuneus', 'Cingulate_Post', 'SMA', 'FEF',
                    'MTL(L)', 'MTL(R)']

naming_list = ['t_b', 't_l', 't_r', 't_fpc', 't_tpc']

# output_array = np.zeros((5, len(list(range(2, 14))), len(label_list), 241))
output_array = np.zeros((5, len(list(range(2, 14))), 241))
for ith_condition in list(range(0, len(naming_list))):
    cond_name = naming_list[ith_condition]
    index_sub = 0
    for ith_sub in list(range(2, 14)):
        # stcs_epoch_morphed1 = np.load(
        #     '/Users/boo/Desktop/MEG_data_script/analysis_conn/stc_epoch_ego_sub' +
        #     str(ith_sub) + '_100hz' + '_' + cond_name +
        #     '.npy', allow_pickle=True)
        stcs_epoch_morphed = np.load(
            '/Users/boo/Desktop/MEG_data_script/analysis_source_result/source_data_ndarray_100hz_sfre_100.npy',allow_pickle=True)

        stcs_epoch_morphed = stcs_epoch_morphed.tolist()
        # stc_epoch_label = mne.extract_label_time_course(stcs_epoch_morphed, label_list, src_fs, mode='mean', #mean_flip
        #                                                 return_generator=False)
        new_list = []
        for ith_epoch in list(range(0, len(stcs_epoch_morphed))):
            curr_etc = stcs_epoch_morphed[ith_epoch]
            new_list.append(curr_etc.in_label(label_list[6]).data)
        new_list = np.mean(new_list, axis=0)

        stc_epoch_label_mean_epoch = np.mean(new_list, axis=0)
        output_array[ith_condition, index_sub, ...] = stc_epoch_label_mean_epoch
        index_sub = index_sub + 1

segm = 3
data_b = output_array[0, ...]
data_lr = (output_array[1, ...] + output_array[2, ...]) / 2
# resampled_b = np.zeros([data_b.shape[0], data_b.shape[1], len(np.arange(0, 220, segm))])
# resampled_lr = np.zeros([data_lr.shape[0], data_lr.shape[1], len(np.arange(0, 220, segm))])
resampled_b = np.zeros([data_b.shape[0], len(np.arange(0, 220, segm))])
resampled_lr = np.zeros([data_lr.shape[0], len(np.arange(0, 220, segm))])
for ind, ith_ts in enumerate(np.arange(0, 220, segm)):
    print(range(ith_ts, ith_ts + segm))
    # resampled_b[..., ind] = np.mean(data_b[..., range(ith_ts, ith_ts + segm)], axis=2)
    # resampled_lr[..., ind] = np.mean(data_lr[..., range(ith_ts, ith_ts + segm)], axis=2)
    resampled_b[..., ind] = np.mean(data_b[..., range(ith_ts, ith_ts + segm)], axis=1)
    resampled_lr[..., ind] = np.mean(data_lr[..., range(ith_ts, ith_ts + segm)], axis=1)
tp_list = [round(x, 2) for x in list(np.arange(-0.2, 2, segm / 100))]

# prepare mean and se for each label

for ith_label in list(range(0, len(label_list))):

    # resampled_b_for_curr_label = resampled_b[:, ith_label, :]
    # resampled_lr_for_curr_label = resampled_lr[:, ith_label, :]
    resampled_b_for_curr_label = resampled_b
    resampled_lr_for_curr_label = resampled_lr

    mean_b = np.mean(resampled_b_for_curr_label, axis=0)
    se_b = np.std(resampled_b_for_curr_label, axis=0) / np.sqrt(resampled_b_for_curr_label.shape[0])
    mean_lr = np.mean(resampled_lr_for_curr_label, axis=0)
    se_lr = np.std(resampled_lr_for_curr_label, axis=0) / np.sqrt(resampled_lr_for_curr_label.shape[0])

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    # ax.axhline(y=0, c="black", linewidth=1, zorder=0)
    ax.axvline(x=0, c="black", linewidth=1, zorder=0)
    ax.axvline(x=1, c="black", linewidth=1, zorder=0)
    ax.plot(tp_list, mean_b, label='Back', color='cornflowerblue')  # Back & Front
    ax.plot(tp_list, mean_lr, label='Left & Right', color='darkorange')
    ax.fill_between(tp_list, mean_b - se_b, mean_b + se_b, color='cornflowerblue', alpha=0.2)
    ax.fill_between(tp_list, mean_lr - se_lr, mean_lr + se_lr, color='darkorange', alpha=0.2)

    T_obs, clusters, cluster_p_values, H0 = \
        permutation_cluster_1samp_test(resampled_b_for_curr_label - resampled_lr_for_curr_label, n_permutations=10000,
                                       tail=0, threshold=2.2)
    for i_c, c in enumerate(clusters):
        c = c[0]
        if cluster_p_values[i_c] <= 0.05 and c.stop - c.start > 1:
            h = ax.axvspan(tp_list[c.start], tp_list[c.stop - 1], color='red', alpha=0.2)
            # plt.legend((h,), ('cluster p-value < 0.05',))
        else:
            plt.legend()
    plt.xlim(-0.2, 2)
    # plt.ylim(-3, 3)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Source Power')
    plt.savefig(
        '/Users/boo/Desktop/MEG_data_script/analysis_frequency_MTL_tp/'+ yaxis_label_list[
            ith_label] + '.png')
    plt.close()
