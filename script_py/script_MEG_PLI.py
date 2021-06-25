# source contrast get averaged
# reset -f
import os
import numpy
import numpy as np
import mne
from mne.io import read_raw_fif
from scipy import stats as stats
from mne.stats import permutation_t_test
from mne.stats import (spatio_temporal_cluster_1samp_test,
                       summarize_clusters_stc)
from sklearn.base import clone
from mne.connectivity import spectral_connectivity, seed_target_indices
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
label_name_list_mtl = ['Hippocampus', 'ParaHippocampal', 'Enterinal', 'Perirhinal']
hemi_pool = ['_lh', '_rh']
label_list_path = []
for r, d, f in os.walk('/Users/boo/Desktop/MEG_data_script/aal_51/labels_conn_each_band_fs/'):
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
label_SMA = mne.read_label(
    '/Users/boo/Desktop/MEG_data_script/aal_51/labels_conn_each_band_fs/SMA_rh.label') + mne.read_label(
    '/Users/boo/Desktop/MEG_data_script/aal_51/labels_conn_each_band_fs/SMA_lh.label')
label_FEF = mne.read_label(
    '/Users/boo/Desktop/MEG_data_script/aal_51/labels_conn_each_band_fs/FEF_rh.label') + mne.read_label(
    '/Users/boo/Desktop/MEG_data_script/aal_51/labels_conn_each_band_fs/FEF_lh.label')
label_list.append(label_parietal)
label_list.append(label_precuneus)
label_list.append(label_SMA)
label_list.append(label_FEF)
for ith_label in list(range(0, len(label_list_path))):
    label_list.append(mne.read_label(label_list_path[ith_label]))

yaxis_label_list = ['Parietal', 'Precuneus', 'SMA', 'FEF',
                    'HPC(L)', 'PHC(L)', 'ERC(L)', 'PRC(L)',
                    'HPC(R)', 'PHC(R)', 'ERC(R)', 'PRC(R)']
# band
iter_freqs = [
    ('Alpha', 8, 13),
    ('Beta', 13, 30),
    ('Low gamma', 30, 60),
    ('High gamma', 60, 99)
]
method_pool = ['pli'] #'plv', 'coh', 'pli'
naming_list = ['t_b', 't_l', 't_r', 't_nc', 't_tpc', 't_fpc']

# the maximum point for b-lr is 0.28
# the maximum point for lr-b is 0.76
# 150 200 250 300 350 400

time_seed_pool = [0.28, 0.76]
time_sep_pool = [0.375, 0.4, 0.5, 0.6, 0.7] #0.15, 0.2, 0.25, 0.3, 0.35, 0.4
tmin_pool = []
tmax_pool = []
for ith_prep1 in list(range(0, len(time_seed_pool))):
    for ith_prep2 in list(range(0, len(time_sep_pool))):
        tmin_pool.append(time_seed_pool[ith_prep1] - time_sep_pool[ith_prep2] / 2)
        tmax_pool.append(time_seed_pool[ith_prep1] + time_sep_pool[ith_prep2] / 2)

curr_tp = 0
for ith_tp in list(range(0, len(tmin_pool))):
    curr_tmin = round(tmin_pool[ith_tp], 3)
    curr_tmax = round(tmax_pool[ith_tp], 3)

    for ith_method in list(range(0, len(method_pool))):
        curr_method = method_pool[ith_method]

        for ith_band in list(range(0, len(iter_freqs))):
            curr_fre_info = iter_freqs[ith_band]
            band_name = curr_fre_info[0]
            vmin = curr_fre_info[1]
            vmax = curr_fre_info[2]

            for ith_condition in list(range(0, len(naming_list))):
                curr_condition = naming_list[ith_condition]
                index_sub = 0
                output_array = np.zeros((len(list(range(2, 14))), len(label_list), len(label_list)))
                for ith_sub in list(range(2, 14)):
                    stcs_epoch_morphed_nocrop = np.load(
                        '/Users/boo/Desktop/MEG_data_script/analysis_conn/stc_ego_epoch_sub' +
                        str(ith_sub) + '_200hz_' + curr_condition +
                        '.npy', allow_pickle=True)
                    stcs_evoke_morphed_nocrop = np.load(
                        '/Users/boo/Desktop/MEG_data_script/analysis_conn/stc_sourceEstimate_ego_evoke_sub' +
                        str(ith_sub) + '_200hz_' + curr_condition +
                        '.npy', allow_pickle=True)
                    stcs_epoch_morphed_nocrop = stcs_epoch_morphed_nocrop.tolist()
                    stcs_evoke_morphed_nocrop = stcs_evoke_morphed_nocrop.tolist()

                    # crop time period
                    stcs_epoch_morphed = []
                    for ith_ele in list(range(0, len(stcs_epoch_morphed_nocrop))):
                        stcs_epoch_morphed.append(
                            stcs_epoch_morphed_nocrop[ith_ele].crop(tmin=curr_tmin, tmax=curr_tmax))
                    stcs_evoke_morphed = stcs_evoke_morphed_nocrop.crop(tmin=curr_tmin, tmax=curr_tmax)

                    seed_idx_pool = []
                    for ith_seed in list(range(0, len(yaxis_label_list))):
                        # search max vertice
                        seed_pool_ts_evoke = stcs_evoke_morphed.in_label(label_list[ith_seed])
                        src_pow = np.sum(seed_pool_ts_evoke.data ** 2, axis=1)
                        total_seed_vertice_list = seed_pool_ts_evoke.vertices[0].tolist() + seed_pool_ts_evoke.vertices[
                            1].tolist()
                        seed_vertno = total_seed_vertice_list[np.argmax(src_pow)]
                        total_wb_vertice_list = stcs_evoke_morphed.vertices[0].tolist() + stcs_evoke_morphed.vertices[
                            1].tolist()
                        seed_idx_pool.append(np.searchsorted(total_wb_vertice_list, seed_vertno))

                    # create max epoch array for conn
                    conn_array = np.zeros((len(yaxis_label_list), len(yaxis_label_list), 1))
                    for ith_curr_seed in list(range(0, len(yaxis_label_list))):
                        max_epoch_array = np.zeros(
                            (np.shape(stcs_epoch_morphed)[0], 1, np.shape(stcs_evoke_morphed)[1]))
                        epoch_array = np.zeros(
                            (np.shape(stcs_epoch_morphed)[0], len(yaxis_label_list), np.shape(stcs_evoke_morphed)[1]))
                        for ith_epoch in list(range(0, np.shape(stcs_epoch_morphed)[0])):
                            max_epoch_array[ith_epoch, 0, ...] = stcs_epoch_morphed[ith_epoch].data[
                                seed_idx_pool[ith_curr_seed], ...]
                            for ith_other_seed in list(range(0, len(yaxis_label_list))):
                                epoch_array[ith_epoch, ith_other_seed, ...] = stcs_epoch_morphed[ith_epoch].data[
                                    seed_idx_pool[ith_other_seed], ...]
                        # create indices
                        comb_ts = list(zip(max_epoch_array, epoch_array))
                        indices = seed_target_indices([0], np.arange(1, 13))
                        con, freqs, times, n_epochs, n_tapers = spectral_connectivity(
                            comb_ts, method=curr_method, sfreq=200, fmin=vmin, fmax=vmax, mode='fourier',
                            indices=indices, faverage=True)  # fourier
                        conn_array[ith_curr_seed, ...] = con

                    output_array[index_sub, ...] = conn_array[..., 0]
                    index_sub = index_sub + 1
                np.save('/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' +
                        band_name + '_' + curr_condition + '_' + str(curr_tmin) + '_' + str(curr_tmax) + '.npy',
                        output_array)
    curr_tp = curr_tp + 1




## watching
import os
import numpy
import numpy as np
from scipy import stats
import matplotlib.pylab as plt

method_pool = ['pli'] #'plv', 'coh', 'pli'
naming_list = ['t_b', 't_l', 't_r', 't_nc', 't_tpc', 't_fpc']
iter_freqs = [
    ('Alpha', 8, 13),
    ('Beta', 13, 30),
    ('Low gamma', 30, 60),
    ('High gamma', 60, 99)
]
yaxis_label_list = ['Parietal', 'Precuneus', 'SMA', 'FEF',
                    'HPC(L)', 'PHC(L)', 'ERC(L)', 'PRC(L)',
                    'HPC(R)', 'PHC(R)', 'ERC(R)', 'PRC(R)']

yaxis_label = ['Parietal-SMA', 'Parietal-FEF', 'Precuneus-SMA','Precuneus-FEF',
               'ERC(R)-SMA', 'ERC(R)-FEF', 'ERC(R)-Parietal', 'ERC(R)-Precuneus']

fontsize = 7
time_seed_pool = [0.28, 0.76]
time_sep_pool = [0.375, 0.4, 0.5, 0.6, 0.7] #[0.15, 0.2, 0.25, 0.3, 0.35, 0.4]
tmin_pool = []
tmax_pool = []
for ith_prep1 in list(range(0, len(time_seed_pool))):
    for ith_prep2 in list(range(0, len(time_sep_pool))):
        tmin_pool.append(time_seed_pool[ith_prep1] - time_sep_pool[ith_prep2] / 2)
        tmax_pool.append(time_seed_pool[ith_prep1] + time_sep_pool[ith_prep2] / 2)

for ith_band in list(range(0, len(iter_freqs))):
    curr_fre_info = iter_freqs[ith_band]
    band_name = curr_fre_info[0]

    plot_array = np.zeros((10, len(yaxis_label)))
    title_array = np.array(range(10), dtype='<U20')
    ith_position=0
    for ith_method in list(range(0, len(method_pool))):
        curr_method = method_pool[ith_method]
        for ith_tp in list(range(0, len(tmin_pool))):
            curr_tmin = round(tmin_pool[ith_tp], 3)
            curr_tmax = round(tmax_pool[ith_tp], 3)

            curr_array_b = np.load('/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' +
                                   band_name + '_' + 't_b' + '_' + str(curr_tmin) + '_' + str(curr_tmax) + '.npy')
            curr_array_l = np.load('/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' +
                                   band_name + '_' + 't_l' + '_' + str(curr_tmin) + '_' + str(curr_tmax) + '.npy')
            curr_array_r = np.load('/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' +
                                   band_name + '_' + 't_r' + '_' + str(curr_tmin) + '_' + str(curr_tmax) + '.npy')
            output_array_b_lr = curr_array_b - (curr_array_l + curr_array_r) / 2
            statistic, pvalue = stats.ttest_1samp(output_array_b_lr, 0, axis=0)
            plot_array[ith_position, ...] = np.array(
                (statistic[0][2], statistic[0][3], statistic[1][2], statistic[1][3],
                 statistic[10][2], statistic[10][3], statistic[10][0], statistic[10][1]))

            title_array[ith_position]= np.array((str(curr_tmin) + '-' + str(curr_tmax) + 's(' + curr_method + ')'))
            ith_position = ith_position+1

    fig, axes = plt.subplots(nrows=1, ncols=10, figsize=(30, 3))  # figsize=(16, 8.5)
    ith_plot = 0
    for ax in axes.flat:
        ax.set_xticklabels(yaxis_label, rotation=90, fontsize=fontsize)
        ax.set_xticks(np.arange(len(yaxis_label)))
        ax.bar(yaxis_label, plot_array[ith_plot], width=0.6, color='0.5', edgecolor='black', linewidth=1, capsize=10)
        ax.set_ylim([-3, 3])
        ax.axhline(y=2.2, ls='--', linewidth=1, color='r')
        ax.axhline(y=-2.2, ls='--', linewidth=1, color='r')
        ax.set_title(title_array[ith_plot], fontsize=fontsize)
        ax.tick_params(labelsize=fontsize)
        ax.set_aspect('auto')
        ith_plot = ith_plot+1
    plt.subplots_adjust(left=.03, right=.97, top=0.9, bottom=0.35, wspace=0.5, hspace=0)
    plt.savefig(
        '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/connectivity_' + band_name  + '.png')  # bbox_inches='tight'
    plt.close()



## make figure horizontal bar
import os
import numpy
import numpy as np
from scipy import stats
import matplotlib.pylab as plt
import pandas as pd

method_pool = ['pli']  # 'plv', 'coh', 'pli'
naming_list = ['t_b', 't_l', 't_r', 't_nc', 't_tpc', 't_fpc']

fontsize = 17
time_seed_pool = [0.28, 0.76]
band_name = 'Beta'
curr_method = 'pli'

tmin_t1 = round(time_seed_pool[0] - 0.2, 3)
tmax_t1 = round(time_seed_pool[0] + 0.2, 3)
tmin_t2 = round(time_seed_pool[1] - 0.2, 3)
tmax_t2 = round(time_seed_pool[1] + 0.2, 3)

curr_array_b_t1 = np.load(
    '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_b' + '_' + str(
        tmin_t1) + '_' + str(tmax_t1) + '.npy')
curr_array_l_t1 = np.load(
    '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_l' + '_' + str(
        tmin_t1) + '_' + str(tmax_t1) + '.npy')
curr_array_r_t1 = np.load(
    '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_r' + '_' + str(
        tmin_t1) + '_' + str(tmax_t1) + '.npy')

curr_array_b_t2 = np.load(
    '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_b' + '_' + str(
        tmin_t2) + '_' + str(tmax_t2) + '.npy')
curr_array_l_t2 = np.load(
    '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_l' + '_' + str(
        tmin_t2) + '_' + str(tmax_t2) + '.npy')
curr_array_r_t2 = np.load(
    '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_r' + '_' + str(
        tmin_t2) + '_' + str(tmax_t2) + '.npy')

output_array_b_lr_t1 = curr_array_b_t1 - (curr_array_l_t1 + curr_array_r_t1) / 2
output_array_b_lr_t2 = curr_array_b_t2 - (curr_array_l_t2 + curr_array_r_t2) / 2
statistic_t1, pvalue_t1 = stats.ttest_1samp(output_array_b_lr_t1, 0, axis=0)
statistic_t2, pvalue_t2 = stats.ttest_1samp(output_array_b_lr_t2, 0, axis=0)

mean_t1 = np.mean(output_array_b_lr_t1, axis=0)
mean_t2 = np.mean(output_array_b_lr_t2, axis=0)
se_t1 = np.std(output_array_b_lr_t1, axis=0)/ np.sqrt(12)
se_t2 = np.std(output_array_b_lr_t2, axis=0)/ np.sqrt(12)

# stats.ttest_rel(output_array_b_lr_t1[..., 10,0], output_array_b_lr_t2[..., 10,0])
stats.ttest_1samp(output_array_b_lr_t2[..., 3,0], 0)
# plot_array_t1 = [statistic_t1[3][0], statistic_t1[2][0], statistic_t1[8][0], statistic_t1[9][0], statistic_t1[11][0], statistic_t1[10][0]]
# plot_array_t2 = [statistic_t2[3][0], statistic_t2[2][0], statistic_t2[8][0], statistic_t2[9][0], statistic_t2[11][0], statistic_t2[10][0]]

t1_str = str(tmin_t1)+' ~ '+str(tmax_t1)+'s'
t2_str = str(tmin_t2)+' ~ '+str(tmax_t2)+'s'

yaxis_label_list = ['Parietal', 'Precuneus', 'SMA', 'FEF',
                    'HPC(L)', 'PHC(L)', 'ERC(L)', 'PRC(L)',
                    'HPC(R)', 'PHC(R)', 'ERC(R)', 'PRC(R)']
# yaxis_label = ['FEF-Parietal', 'SMA-Parietal', 'HPC(R)-Parietal', 'PHC(R)-Parietal', 'PRC(R)-Parietal',
#                'ERC(R)-Parietal']
yaxis_label = ['FEF-Precuneus', 'SMA-Precuneus', 'HPC(R)-Precuneus', 'PHC(R)-Precuneus', 'PRC(R)-Precuneus',
               'ERC(R)-Precuneus']
ith_region = 1
dataFrame_mean = pd.DataFrame(data=[[mean_t1[3][ith_region], mean_t2[3][ith_region]], [mean_t1[2][ith_region], mean_t2[2][ith_region]], \
                               [mean_t1[8][ith_region], mean_t2[8][ith_region]], [mean_t1[9][ith_region], mean_t2[9][ith_region]], \
                               [mean_t1[11][ith_region], mean_t2[11][ith_region]], [mean_t1[10][ith_region], mean_t2[10][ith_region]]],
                         index=yaxis_label,
                         columns=[t1_str, t2_str])
dataFrame_se = pd.DataFrame(data=[[se_t1[3][ith_region], se_t2[3][ith_region]], [se_t1[2][ith_region], se_t2[2][ith_region]], \
                               [se_t1[8][ith_region], se_t2[8][ith_region]], [se_t1[9][ith_region], se_t2[9][ith_region]], \
                               [se_t1[11][ith_region], se_t2[11][ith_region]], [se_t1[10][ith_region], se_t2[10][ith_region]]],
                         index=yaxis_label,
                         columns=[t1_str, t2_str])
handle = dataFrame_mean.plot.barh(xerr=dataFrame_se, figsize=(6, 6), legend=False, color=['darkgreen', 'red'])
handle.spines['right'].set_visible(False)
handle.spines['top'].set_visible(False)
handle.set_yticklabels(yaxis_label, rotation=0, fontsize=fontsize)
handle.set_xticks([-0.15, 0, 0.1])
handle.set_xlabel('t value', fontsize=fontsize)
handle.axvline(x=0, ls='-', linewidth=0.5, color='black')
handle.invert_yaxis()  # labels read top-to-bottom
handle.tick_params(labelsize=fontsize)
handle.set_aspect('auto')
# handle.legend(loc='upper right', prop={'size': fontsize})
plt.subplots_adjust(left=.35, right=.97, top=0.97, bottom=0.15, wspace=0.5, hspace=0)
plt.savefig(
    '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/Fig_6_Precuneus_roi_' + band_name + '_' + '.png')  # bbox_inches='tight'
plt.close()




## make figure vertical bar - old
import os
import numpy
import numpy as np
from scipy import stats
import matplotlib.pylab as plt
import pandas as pd

fontsize = 29
method_pool = ['pli']  # 'plv', 'coh', 'pli'
naming_list = ['t_b', 't_l', 't_r', 't_nc', 't_tpc', 't_fpc']
band_list = ['Alpha', 'Beta', 'Low gamma', 'High gamma']
seed_pool = ['Parietal', 'Precuneus']
time_seed_pool = [0.28, 0.76]
curr_method = 'pli'

for ith_region in list(range(0, 2)):  # 1 for precuneus 0 for parietal cortex
    for ith_band in list(range(0, len(band_list))):
        for ith_time_p in list(range(0, len(time_seed_pool))):

            band_name = band_list[ith_band]

            tmin = round(time_seed_pool[ith_time_p] - 0.2, 3)
            tmax = round(time_seed_pool[ith_time_p] + 0.2, 3)

            curr_array_b = np.load(
                '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_b' + '_' + str(
                    tmin) + '_' + str(tmax) + '.npy')
            curr_array_l = np.load(
                '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_l' + '_' + str(
                    tmin) + '_' + str(tmax) + '.npy')
            curr_array_r = np.load(
                '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_r' + '_' + str(
                    tmin) + '_' + str(tmax) + '.npy')

            if ith_time_p == 0:
                # color = 'red'
                output_array_contrast = curr_array_b - (curr_array_l + curr_array_r) / 2
            if ith_time_p == 1:
                # color = 'darkgreen'
                output_array_contrast = (curr_array_l + curr_array_r) / 2 - curr_array_b

            mean = np.mean(output_array_contrast, axis=0)
            se = np.std(output_array_contrast, axis=0) / np.sqrt(12)

            # statistic
            statistic, pvalue = stats.ttest_1samp(output_array_contrast, 0, axis=0)
            # stats.ttest_rel(output_array_b_lr_t1[..., 10,0], output_array_b_lr_t2[..., 10,0])
            stat_fef, pval_fef = stats.ttest_1samp(output_array_contrast[..., 3, ith_region], 0)
            stat_sma, pval_sma = stats.ttest_1samp(output_array_contrast[..., 2, ith_region], 0)
            stat_hpc, pval_hpc = stats.ttest_1samp(output_array_contrast[..., 8, ith_region], 0)
            stat_phc, pval_phc = stats.ttest_1samp(output_array_contrast[..., 9, ith_region], 0)
            stat_prc, pval_prc = stats.ttest_1samp(output_array_contrast[..., 11, ith_region], 0)
            stat_erc, pval_erc = stats.ttest_1samp(output_array_contrast[..., 10, ith_region], 0)

            yaxis_label_list = ['Parietal', 'Precuneus', 'SMA', 'FEF',
                                'HPC(L)', 'PHC(L)', 'ERC(L)', 'PRC(L)',
                                'HPC(R)', 'PHC(R)', 'ERC(R)', 'PRC(R)']  # for reference
            label_x = ['FEF', 'SMA', 'HPC', 'PHC', 'PRC', 'ERC']
            color = ['limegreen', 'limegreen', 'red', 'red', 'red', 'red']

            value_y = [mean[3][ith_region], mean[2][ith_region],
                       mean[8][ith_region], mean[9][ith_region],
                       mean[11][ith_region], mean[10][ith_region]]

            value_errorbar = [se[3][ith_region], se[2][ith_region],
                              se[8][ith_region], se[9][ith_region],
                              se[11][ith_region], se[10][ith_region]]

            fig, ax = plt.subplots(figsize=(7, 5.5))
            ax.bar([1, 2, 4, 5, 6, 7], value_y, width=0.5, yerr=value_errorbar, capsize=3, color=color)  # (89/255, 88/255, 89/255)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.set_xticks([1, 2, 4, 5, 6, 7])
            ax.set_xticklabels(label_x, rotation=45, fontsize=fontsize-3)
            ax.set_yticks([-0.08, 0, 0.14])
            ax.tick_params(labelsize=fontsize)
            ax.set_aspect('auto')
            ax.set_ylabel('PLI', fontsize=fontsize)
            # ax.axvline(x=0, ls='-', linewidth=0.5, color='black')
            # ax.invert_xaxis()  # labels read top-to-bottom
            # handle.legend(loc='upper right', prop={'size': fontsize})
            plt.subplots_adjust(left=.25, right=.97, top=0.97, bottom=0.15, wspace=0.5, hspace=0)
            plt.savefig(
                '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/Fig_6_seed_' + seed_pool[ith_region] + '_band_' + band_name + '_' + str(time_seed_pool[ith_time_p]) + '.png', bbox_inches='tight')  # bbox_inches='tight'
            plt.close()



## make figure vertical bar - new - paired t test
import os
import numpy
import numpy as np
from scipy import stats
import matplotlib.pylab as plt
import pandas as pd
fontsize = 29
method_pool = ['pli']  # 'plv', 'coh', 'pli'
naming_list = ['t_b', 't_l', 't_r', 't_nc', 't_tpc', 't_fpc']
band_list = ['Alpha', 'Beta', 'Low gamma', 'High gamma']
seed_pool = ['Parietal', 'Precuneus']
time_seed_pool = [0.28, 0.76]
curr_method = 'pli'

for ith_region in list(range(0, 2)):  # 1 for precuneus 0 for parietal cortex
    for ith_band in list(range(0, len(band_list))):

            band_name = band_list[ith_band]

            tmin_early = round(time_seed_pool[0] - 0.2, 3)
            tmax_early = round(time_seed_pool[0] + 0.2, 3)
            tmin_late = round(time_seed_pool[1] - 0.2, 3)
            tmax_late = round(time_seed_pool[1] + 0.2, 3)
            curr_array_b_early = np.load(
                '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_b' + '_' + str(
                    tmin_early) + '_' + str(tmax_early) + '.npy')
            curr_array_l_early = np.load(
                '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_l' + '_' + str(
                    tmin_early) + '_' + str(tmax_early) + '.npy')
            curr_array_r_early = np.load(
                '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_r' + '_' + str(
                    tmin_early) + '_' + str(tmax_early) + '.npy')
            curr_array_b_late = np.load(
                '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_b' + '_' + str(
                    tmin_late) + '_' + str(tmax_late) + '.npy')
            curr_array_l_late = np.load(
                '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_l' + '_' + str(
                    tmin_late) + '_' + str(tmax_late) + '.npy')
            curr_array_r_late = np.load(
                '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_r' + '_' + str(
                    tmin_late) + '_' + str(tmax_late) + '.npy')

            output_array_contrast_early = curr_array_b_early - (curr_array_l_early + curr_array_r_early) / 2
            output_array_contrast_late = curr_array_b_late - (curr_array_l_late + curr_array_r_late) / 2

            mean_early = np.mean(output_array_contrast_early, axis=0)
            mean_late = np.mean(output_array_contrast_late, axis=0)
            se_early = np.std(output_array_contrast_early, axis=0) / np.sqrt(12)
            se_late = np.std(output_array_contrast_late, axis=0) / np.sqrt(12)

            # two sample t test
            # statistic, pvalue = stats.ttest_1samp(output_array_contrast_early, 0, axis=0)
            # # stats.ttest_rel(output_array_b_lr_t1[..., 10,0], output_array_b_lr_t2[..., 10,0])
            # stat_fef, pval_fef = stats.ttest_1samp(, 0)
            # stat_sma, pval_sma = stats.ttest_1samp(output_array_contrast_early[..., 2, ith_region], 0)
            # stat_hpc, pval_hpc = stats.ttest_1samp(output_array_contrast_early[..., 8, ith_region], 0)
            # stat_phc, pval_phc = stats.ttest_1samp(output_array_contrast_early[..., 9, ith_region], 0)
            # stat_prc, pval_prc = stats.ttest_1samp(output_array_contrast_early[..., 11, ith_region], 0)
            # stat_erc, pval_erc = stats.ttest_1samp(output_array_contrast_early[..., 10, ith_region], 0)

            # paired t test
            stat_fef, pval_fef = stats.ttest_rel(output_array_contrast_early[..., 3, ith_region], output_array_contrast_late[..., 3, ith_region])
            stat_sma, pval_sma = stats.ttest_rel(output_array_contrast_early[..., 2, ith_region], output_array_contrast_late[..., 2, ith_region])
            stat_hpc, pval_hpc = stats.ttest_rel(output_array_contrast_early[..., 8, ith_region], output_array_contrast_late[..., 8, ith_region])
            stat_phc, pval_phc = stats.ttest_rel(output_array_contrast_early[..., 9, ith_region], output_array_contrast_late[..., 9, ith_region])
            stat_erc, pval_erc = stats.ttest_rel(output_array_contrast_early[..., 10, ith_region], output_array_contrast_late[..., 10, ith_region])
            stat_prc, pval_prc = stats.ttest_rel(output_array_contrast_early[..., 11, ith_region], output_array_contrast_late[..., 11, ith_region])
            print('seed:' + seed_pool[ith_region] + ' band:' + band_list[ith_band] + ' fef' + ' tval:' + str(stat_fef) + ' pval:' + str(pval_fef))
            print('seed:' + seed_pool[ith_region] + ' band:' + band_list[ith_band] + ' sma' + ' tval:' + str(stat_sma) + ' pval:' + str(pval_sma))
            print('seed:' + seed_pool[ith_region] + ' band:' + band_list[ith_band] + ' hpc' + ' tval:' + str(stat_hpc) + ' pval:' + str(pval_hpc))
            print('seed:' + seed_pool[ith_region] + ' band:' + band_list[ith_band] + ' phc' + ' tval:' + str(stat_phc) + ' pval:' + str(pval_phc))
            print('seed:' + seed_pool[ith_region] + ' band:' + band_list[ith_band] + ' prc' + ' tval:' + str(stat_prc) + ' pval:' + str(pval_prc))
            print('seed:' + seed_pool[ith_region] + ' band:' + band_list[ith_band] + ' erc' + ' tval:' + str(stat_erc) + ' pval:' + str(pval_erc))

            # reference
            yaxis_label_list = ['Parietal', 'Precuneus', 'SMA', 'FEF',
                                'HPC(L)', 'PHC(L)', 'ERC(L)', 'PRC(L)',
                                'HPC(R)', 'PHC(R)', 'ERC(R)', 'PRC(R)']  # for reference
            # array
            label_x = ['HPC', 'PHC', 'PRC', 'ERC', 'FEF', 'SMA']
            color_early = ['skyblue', 'skyblue', 'skyblue', 'skyblue', 'gold', 'gold']
            color_late = ['blue', 'blue', 'blue', 'blue', 'darkgoldenrod', 'darkgoldenrod']
            value_y_early = [mean_early[8][ith_region], mean_early[9][ith_region], mean_early[11][ith_region], mean_early[10][ith_region],
                       mean_early[3][ith_region], mean_early[2][ith_region]]
            value_y_late = [mean_late[8][ith_region], mean_late[9][ith_region], mean_late[11][ith_region], mean_late[10][ith_region],
                       mean_late[3][ith_region], mean_late[2][ith_region]]
            value_errorbar_early = [se_early[8][ith_region], se_early[9][ith_region], se_early[11][ith_region], se_early[10][ith_region],
                              se_early[3][ith_region], se_early[2][ith_region]]
            value_errorbar_late = [se_late[8][ith_region], se_late[9][ith_region], se_late[11][ith_region], se_late[10][ith_region],
                              se_late[3][ith_region], se_late[2][ith_region]]
            width = 0.25  # the width of the bars
            ind = np.arange(len(value_y_early))

            fig, ax = plt.subplots(figsize=(10, 4))
            ax.bar(ind - width / 2, value_y_early, width, yerr=value_errorbar_early, capsize=3, color=color_early)
            ax.bar(ind + width / 2, value_y_late, width, yerr=value_errorbar_late, capsize=3, color=color_late)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.set_xticks(ind)
            if ith_band==0:
                ax.set_xticklabels(label_x, rotation=45, fontsize=fontsize-3)
            else:
                ax.set_xticklabels([])
            ax.set_yticks([-0.17, 0, 0.14])
            ax.tick_params(labelsize=fontsize)
            ax.set_aspect('auto')
            ax.set_ylabel('Back - Left/Right', fontsize=fontsize)
            # ax.axvline(x=0, ls='-', linewidth=0.5, color='black')
            # ax.invert_xaxis()  # labels read top-to-bottom
            # handle.legend(loc='upper right', prop={'size': fontsize})
            plt.subplots_adjust(left=.25, right=.97, top=0.97, bottom=0.15, wspace=0.5, hspace=0)
            plt.savefig(
                '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/Fig_6_seed_' + seed_pool[ith_region] + '_band_' + band_name + '.png', bbox_inches='tight')  # bbox_inches='tight'
            plt.close()



## make figure vertical bar - new - anova-like
import os
import numpy
import numpy as np
from scipy import stats
import matplotlib.pylab as plt
import pandas as pd
fontsize = 29
method_pool = ['pli']  # 'plv', 'coh', 'pli'
naming_list = ['t_b', 't_l', 't_r', 't_nc', 't_tpc', 't_fpc']
band_list = ['Alpha', 'Beta', 'Low gamma', 'High gamma']
seed_pool = ['Parietal', 'Precuneus']
time_seed_pool = [0.28, 0.76]
curr_method = 'pli'

for ith_region in list(range(0, 2)):  # 1 for precuneus 0 for parietal cortex
    for ith_band in list(range(0, len(band_list))):

            band_name = band_list[ith_band]

            tmin_early = round(time_seed_pool[0] - 0.2, 3)
            tmax_early = round(time_seed_pool[0] + 0.2, 3)
            tmin_late = round(time_seed_pool[1] - 0.2, 3)
            tmax_late = round(time_seed_pool[1] + 0.2, 3)
            curr_array_b_early = np.load(
                '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_b' + '_' + str(
                    tmin_early) + '_' + str(tmax_early) + '.npy')
            curr_array_l_early = np.load(
                '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_l' + '_' + str(
                    tmin_early) + '_' + str(tmax_early) + '.npy')
            curr_array_r_early = np.load(
                '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_r' + '_' + str(
                    tmin_early) + '_' + str(tmax_early) + '.npy')
            curr_array_b_late = np.load(
                '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_b' + '_' + str(
                    tmin_late) + '_' + str(tmax_late) + '.npy')
            curr_array_l_late = np.load(
                '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_l' + '_' + str(
                    tmin_late) + '_' + str(tmax_late) + '.npy')
            curr_array_r_late = np.load(
                '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_r' + '_' + str(
                    tmin_late) + '_' + str(tmax_late) + '.npy')

            output_array_contrast_early = curr_array_b_early - (curr_array_l_early + curr_array_r_early) / 2
            output_array_contrast_late = curr_array_b_late - (curr_array_l_late + curr_array_r_late) / 2

            mean_early = np.mean(output_array_contrast_early, axis=0)
            mean_late = np.mean(output_array_contrast_late, axis=0)
            se_early = np.std(output_array_contrast_early, axis=0) / np.sqrt(12)
            se_late = np.std(output_array_contrast_late, axis=0) / np.sqrt(12)

            # reference
            yaxis_label_list = ['Parietal', 'Precuneus', 'SMA', 'FEF',
                                'HPC(L)', 'PHC(L)', 'ERC(L)', 'PRC(L)',
                                'HPC(R)', 'PHC(R)', 'ERC(R)', 'PRC(R)']  # for reference
            # array
            label_x = ['HPC', 'PHC', 'PRC', 'ERC', 'FEF', 'SMA', 'HPC', 'PHC', 'PRC', 'ERC', 'FEF', 'SMA']
            color = ['blue', 'blue', 'blue', 'blue', 'darkgoldenrod', 'darkgoldenrod', 'blue', 'blue', 'blue', 'blue', 'darkgoldenrod', 'darkgoldenrod']
            value_y = [mean_early[8][ith_region], mean_early[9][ith_region], mean_early[11][ith_region], mean_early[10][ith_region],
                       mean_early[3][ith_region], mean_early[2][ith_region], mean_late[8][ith_region], mean_late[9][ith_region],
                       mean_late[11][ith_region], mean_late[10][ith_region], mean_late[3][ith_region], mean_late[2][ith_region]]
            value_errorbar = [se_early[8][ith_region], se_early[9][ith_region], se_early[11][ith_region], se_early[10][ith_region],
                              se_early[3][ith_region], se_early[2][ith_region], se_late[8][ith_region], se_late[9][ith_region],
                              se_late[11][ith_region], se_late[10][ith_region], se_late[3][ith_region], se_late[2][ith_region]]
            width = 0.5  # the width of the bars
            ind = np.arange(len(value_y))

            fig, ax = plt.subplots(figsize=(12, 4))
            ax.bar([1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14], value_y, width, yerr=value_errorbar, capsize=3, color=color)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.set_xticks([1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14])
            if ith_band==0:
                ax.set_xticklabels(label_x, rotation=45, fontsize=fontsize-3)
            else:
                ax.set_xticklabels([])
            ax.set_yticks([-0.17, 0, 0.14])
            ax.tick_params(labelsize=fontsize)
            ax.set_aspect('auto')
            ax.set_ylabel('Back - Left/Right', fontsize=fontsize)
            plt.subplots_adjust(left=.25, right=.97, top=0.97, bottom=0.15, wspace=0.5, hspace=0)
            plt.savefig(
                '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/Fig_6_seed_' + seed_pool[ith_region] + '_band_' + band_name + '.png', bbox_inches='tight')  # bbox_inches='tight'
            plt.close()




## anova two way
import os
import numpy
import numpy as np
from scipy import stats
import matplotlib.pylab as plt
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import (pairwise_tukeyhsd, MultiComparison)

fontsize = 25
method_pool = ['pli']  # 'plv', 'coh', 'pli'
naming_list = ['t_b', 't_l', 't_r', 't_nc', 't_tpc', 't_fpc']
band_list = ['Alpha', 'Beta', 'Low gamma', 'High gamma']
seed_pool = ['Parietal', 'Precuneus']
time_seed_pool = [0.28, 0.76]
curr_method = 'pli'
yaxis_label_list = ['Parietal', 'Precuneus', 'SMA', 'FEF',
                    'HPC(L)', 'PHC(L)', 'ERC(L)', 'PRC(L)',
                    'HPC(R)', 'PHC(R)', 'ERC(R)', 'PRC(R)']  # for reference
label_x = ['FEF', 'SMA', 'HPC', 'PHC', 'PRC', 'ERC']

for ith_region in list(range(0, len(seed_pool))):  # 1 for precuneus 0 for parietal cortex
    for ith_band in list(range(0, len(band_list))):
        band_name = band_list[ith_band]

        tmin_t1 = round(time_seed_pool[0] - 0.2, 3)
        tmax_t1 = round(time_seed_pool[0] + 0.2, 3)
        tmin_t2 = round(time_seed_pool[1] - 0.2, 3)
        tmax_t2 = round(time_seed_pool[1] + 0.2, 3)

        curr_array_b_t1 = np.load(
            '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_b' + '_' + str(
                tmin_t1) + '_' + str(tmax_t1) + '.npy')
        curr_array_l_t1 = np.load(
            '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_l' + '_' + str(
                tmin_t1) + '_' + str(tmax_t1) + '.npy')
        curr_array_r_t1 = np.load(
            '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_r' + '_' + str(
                tmin_t1) + '_' + str(tmax_t1) + '.npy')

        curr_array_b_t2 = np.load(
            '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_b' + '_' + str(
                tmin_t2) + '_' + str(tmax_t2) + '.npy')
        curr_array_l_t2 = np.load(
            '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_l' + '_' + str(
                tmin_t2) + '_' + str(tmax_t2) + '.npy')
        curr_array_r_t2 = np.load(
            '/Users/boo/Desktop/MEG_data_script/analysis_conn_figures/' + curr_method + '_' + band_name + '_' + 't_r' + '_' + str(
                tmin_t2) + '_' + str(tmax_t2) + '.npy')

        array_t1_fef = curr_array_b_t1[..., 3, ith_region] - (curr_array_l_t1[..., 3, ith_region] + curr_array_r_t1[..., 3, ith_region])/2
        array_t1_sma = curr_array_b_t1[..., 2, ith_region] - (curr_array_l_t1[..., 2, ith_region] + curr_array_r_t1[..., 2, ith_region])/2
        array_t1_hpc = curr_array_b_t1[..., 8, ith_region] - (curr_array_l_t1[..., 8, ith_region] + curr_array_r_t1[..., 8, ith_region])/2
        array_t1_phc = curr_array_b_t1[..., 9, ith_region] - (curr_array_l_t1[..., 9, ith_region] + curr_array_r_t1[..., 9, ith_region])/2
        array_t1_prc = curr_array_b_t1[..., 11, ith_region] - (curr_array_l_t1[..., 11, ith_region] + curr_array_r_t1[..., 11, ith_region])/2
        array_t1_erc = curr_array_b_t1[..., 10, ith_region] - (curr_array_l_t1[..., 10, ith_region] + curr_array_r_t1[..., 10, ith_region])/2

        array_t2_fef = curr_array_b_t2[..., 3, ith_region] - (curr_array_l_t2[..., 3, ith_region] + curr_array_r_t2[..., 3, ith_region])/2
        array_t2_sma = curr_array_b_t2[..., 2, ith_region] - (curr_array_l_t2[..., 2, ith_region] + curr_array_r_t2[..., 2, ith_region])/2
        array_t2_hpc = curr_array_b_t2[..., 8, ith_region] - (curr_array_l_t2[..., 8, ith_region] + curr_array_r_t2[..., 8, ith_region])/2
        array_t2_phc = curr_array_b_t2[..., 9, ith_region] - (curr_array_l_t2[..., 9, ith_region] + curr_array_r_t2[..., 9, ith_region])/2
        array_t2_prc = curr_array_b_t2[..., 11, ith_region] - (curr_array_l_t2[..., 11, ith_region] + curr_array_r_t2[..., 11, ith_region])/2
        array_t2_erc = curr_array_b_t2[..., 10, ith_region] - (curr_array_l_t2[..., 10, ith_region] + curr_array_r_t2[..., 10, ith_region])/2

        statistic, pvalue = stats.ttest_1samp(array_t2_sma, 0, axis=0)

        create_array = {'value': np.concatenate((array_t1_fef, array_t1_sma, array_t1_hpc, array_t1_phc, array_t1_prc, array_t1_erc,
                                                 array_t2_fef, array_t2_sma, array_t2_hpc, array_t2_phc, array_t2_prc, array_t2_erc)),
                    'area': np.concatenate((np.repeat('fef', 12), np.repeat('sma', 12), np.repeat('hpc', 12), np.repeat('phc', 12), np.repeat('prc', 12), np.repeat('erc', 12),
                                            np.repeat('fef', 12), np.repeat('sma', 12), np.repeat('hpc', 12), np.repeat('phc', 12), np.repeat('prc', 12), np.repeat('erc', 12),)),
                    'time': np.concatenate((np.repeat('t1', 12*6),  np.repeat('t2', 12*6)))}

        create_array = {'value': np.concatenate((mean(array_t1_fef), array_t1_sma, array_t1_hpc, array_t1_phc, array_t1_prc, array_t1_erc,
                                                 array_t2_fef, array_t2_sma, array_t2_hpc, array_t2_phc, array_t2_prc, array_t2_erc)),
                    'area': np.concatenate((np.repeat('fef', 12), np.repeat('sma', 12), np.repeat('hpc', 12), np.repeat('phc', 12), np.repeat('prc', 12), np.repeat('erc', 12),
                                            np.repeat('fef', 12), np.repeat('sma', 12), np.repeat('hpc', 12), np.repeat('phc', 12), np.repeat('prc', 12), np.repeat('erc', 12),)),
                    'time': np.concatenate((np.repeat('t1', 12*6),  np.repeat('t2', 12*6)))}
        for_anova = pd.DataFrame(data=create_array)
        moore_lm = ols('value ~ C(area, Sum)*C(time, Sum)', data=for_anova).fit()
        table = sm.stats.anova_lm(moore_lm, typ=1)  # Type 2 ANOVA DataFrame
        print(table)
        print(seed_pool[ith_region])
        print(band_list[ith_band])
        # MultiComp = MultiComparison(for_anova['value'], for_anova['area'])
        # print(MultiComp.tukeyhsd().summary())







