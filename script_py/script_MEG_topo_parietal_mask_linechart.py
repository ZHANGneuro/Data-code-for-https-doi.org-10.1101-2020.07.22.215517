# contrast activity contrast back lr
import os
import mne
from mne.io import read_raw_fif
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
import math
import matplotlib
from mne.viz import topomap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mne.stats import permutation_t_test
from mne.stats import permutation_cluster_1samp_test
from mne.stats import (spatio_temporal_cluster_1samp_test, summarize_clusters_stc)
from mne.viz import plot_topomap

list_data_back = []
list_data_left = []
list_data_right = []
list_data_front = []
list_data_tfpc = []
list_data_ttpc = []

for ith_sub in list(range(2, 14)):
    temp_data_array = "/Users/boo/Desktop/MEG_data_script/PreProcessed_data/artefact_removed_sub" + str(
        ith_sub) + "_raw_100hz_sfre_100.fif"
    temp_event_array = "/Users/boo/Desktop/MEG_data_script/PreProcessed_data/events_post_resample_sub" + str(
        ith_sub) + "_100hz_sfre_100.npy"

    array_data = read_raw_fif(temp_data_array)
    array_event = numpy.load(temp_event_array)

    # pick channel
    all_chan = array_data.ch_names
    picks_mag = mne.pick_types(array_data.info, meg='mag')
    meg_channel = itemgetter(*picks_mag)(all_chan)
    meg_channel = meg_channel[29:301]
    pos = mne.channels.layout._find_topomap_coords(array_data.info, picks_mag[29:301])

    # Compute epochs
    min_onset = -0.2
    max_endpoint = 2
    baseline = (min_onset, 0)
    event_id_back = {'back': 80}
    event_id_front = {'front': 90}
    event_id_left = {'left': 82}
    event_id_right = {'right': 84}
    event_id_tfpc = {'tfpc': 110}
    event_id_ttpc = {'ttpc': 100}

    epochs_back = mne.Epochs(array_data, array_event, picks=meg_channel, event_id=event_id_back, tmin=min_onset, tmax=max_endpoint,
                             baseline=baseline)
    epochs_front = mne.Epochs(array_data, array_event, picks=meg_channel, event_id=event_id_front, tmin=min_onset,
                              tmax=max_endpoint, baseline=baseline)
    epochs_left = mne.Epochs(array_data, array_event, picks=meg_channel, event_id=event_id_left, tmin=min_onset, tmax=max_endpoint,
                             baseline=baseline)
    epochs_right = mne.Epochs(array_data, array_event, picks=meg_channel, event_id=event_id_right, tmin=min_onset,
                              tmax=max_endpoint, baseline=baseline)
    epochs_tfpc = mne.Epochs(array_data, array_event, picks=meg_channel, event_id=event_id_tfpc, tmin=min_onset, tmax=max_endpoint,
                             baseline=baseline)
    epochs_ttpc = mne.Epochs(array_data, array_event, picks=meg_channel, event_id=event_id_ttpc, tmin=min_onset, tmax=max_endpoint,
                             baseline=baseline)

    epochs_back.load_data()
    epochs_front.load_data()
    epochs_left.load_data()
    epochs_right.load_data()
    epochs_tfpc.load_data()
    epochs_ttpc.load_data()

    evoke_back = epochs_back.average()
    evoke_front = epochs_front.average()
    evoke_left = epochs_left.average()
    evoke_right = epochs_right.average()
    evoke_tfpc = epochs_tfpc.average()
    evoke_ttpc = epochs_ttpc.average()

    # add into list
    list_data_back.append(evoke_back.data)
    list_data_front.append(evoke_front.data)
    list_data_left.append(evoke_left.data)
    list_data_right.append(evoke_right.data)
    list_data_tfpc.append(evoke_tfpc.data)
    list_data_ttpc.append(evoke_ttpc.data)

array_back = np.array(list_data_back)
array_left = np.array(list_data_left)
array_right = np.array(list_data_right)
array_front = np.array(list_data_front)
array_tfpc = np.array(list_data_tfpc)
array_ttpc = np.array(list_data_ttpc)

onset_t = 0
end_point = 220
segm = 3
num_dp = np.shape(array_back)[2]

resampled_ima_t_f = np.zeros([array_front.shape[0], array_front.shape[1], len(np.arange(onset_t, end_point, segm))])
resampled_ima_t_b = np.zeros([array_back.shape[0], array_back.shape[1], len(np.arange(onset_t, end_point, segm))])
resampled_ima_t_l = np.zeros([array_left.shape[0], array_left.shape[1], len(np.arange(onset_t, end_point, segm))])
resampled_ima_t_r = np.zeros([array_right.shape[0], array_right.shape[1], len(np.arange(onset_t, end_point, segm))])
resampled_ima_t_fpc = np.zeros([array_tfpc.shape[0], array_tfpc.shape[1], len(np.arange(onset_t, end_point, segm))])
resampled_ima_t_tpc = np.zeros([array_ttpc.shape[0], array_ttpc.shape[1], len(np.arange(onset_t, end_point, segm))])
for ind, ith_ts in enumerate(np.arange(onset_t, end_point, segm)):
    print(range(ith_ts, ith_ts + segm))
    if ith_ts + segm < num_dp:
        resampled_ima_t_f[..., ind] = np.mean(array_front[..., range(ith_ts, ith_ts + segm)], axis=2)
        resampled_ima_t_b[..., ind] = np.mean(array_back[..., range(ith_ts, ith_ts + segm)], axis=2)
        resampled_ima_t_l[..., ind] = np.mean(array_left[..., range(ith_ts, ith_ts + segm)], axis=2)
        resampled_ima_t_r[..., ind] = np.mean(array_right[..., range(ith_ts, ith_ts + segm)], axis=2)
        resampled_ima_t_fpc[..., ind] = np.mean(array_tfpc[..., range(ith_ts, ith_ts + segm)], axis=2)
        resampled_ima_t_tpc[..., ind] = np.mean(array_ttpc[..., range(ith_ts, ith_ts + segm)], axis=2)

the_data_set = (resampled_ima_t_l + resampled_ima_t_r)/2 - resampled_ima_t_b



######################
# for plot topographic example
######################
# mean_curr_ima = mean_curr_ima[:,1]
X = (resampled_ima_t_l + resampled_ima_t_r) / 2 - resampled_ima_t_b
X = np.transpose(X, [0, 2, 1])
T_obs, clusters, cluster_p_values, H0 = clu = \
    spatio_temporal_cluster_1samp_test(X, threshold=2.2, tail=0)
sig_index = clusters[np.where(cluster_p_values <= 0.05)[0][0]]
sig_ch_index = np.unique(sig_index[1])


mean_curr_ima = array_back[1, :, 1]
mean_curr_ima[:] = 0
fig = plt.figure(constrained_layout=False, figsize=[2, 2])
fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.98)
num_row = 1
num_col = 1
gs = matplotlib.gridspec.GridSpec(nrows=num_row, ncols=num_col, figure=fig)
mask = np.zeros(len(meg_channel), dtype=bool)
# sig_ch_list = meg_channel[sig_ch_index[0]:sig_ch_index[-1]]
# sig_ch_index = list(range(sig_cluster.start, sig_cluster.stop))
sig_ch_index = np.unique(sig_ch_index)
mask[sig_ch_index] = True
cur_ax = fig.add_subplot(gs[0, 0])
kwargs = dict(vmin=-6e-14, vmax=6e-14, sensors=False, res=64, names=None, show_names=False,
              mask_params={}, outlines='head', contours=6, image_interp='bilinear', show=False,
              extrapolate='box')
plot_topomap(mean_curr_ima, pos, axes=cur_ax, cmap='seismic', mask=mask, **kwargs)
fig.savefig('/Users/boo/Desktop/mask_lr_b1.png')
plt.close()



######################
# for plot roi time course
######################
tp_list = [round(x, 2) for x in list(np.arange(-0.2, 2, segm/100))]

back_roi = np.mean(resampled_ima_t_b[:, sig_ch_index, :], axis=1)
left_roi = np.mean(resampled_ima_t_l[:, sig_ch_index, :], axis=1)
right_roi = np.mean(resampled_ima_t_r[:, sig_ch_index, :], axis=1)
fpc_roi = np.mean(resampled_ima_t_fpc[:, sig_ch_index, :], axis=1)
tpc_roi = np.mean(resampled_ima_t_tpc[:, sig_ch_index, :], axis=1)

mean1 = np.mean(back_roi, axis=0)
se1 = np.std(back_roi, axis=0)/np.sqrt(12)
mean2 = np.mean((left_roi+right_roi)/2, axis=0)
se2 = np.std((left_roi+right_roi)/2, axis=0)/np.sqrt(12)

T_obs, clusters, cluster_p_values, H0 = \
    permutation_cluster_1samp_test((left_roi+right_roi)/2 - back_roi, n_permutations=10000, tail=0, threshold=1.79)

fig, ax = plt.subplots(1, 1, figsize=(10, 7))
ax.axvline(x=0, c="gray", linewidth=2, zorder=0, linestyle='--')
ax.axvline(x=1, c="gray", linewidth=2, zorder=0, linestyle='--')
ax.plot(tp_list, mean1, label='Back', color='cornflowerblue')  # Back & Front
ax.plot(tp_list, mean2, label='Left/Right', color='darkorange')
ax.fill_between(tp_list, mean1 - se1, mean1 + se1, color='cornflowerblue', alpha=0.2)
ax.fill_between(tp_list, mean2 - se2, mean2 + se2, color='darkorange', alpha=0.2)
for i_c, c in enumerate(clusters):
    c = c[0]
    if cluster_p_values[i_c] <= 0.05 and c.stop - c.start > 1:
        h = ax.axvspan(tp_list[c.start], tp_list[c.stop - 1], color='red', alpha=0.1)
        # plt.legend((h,), ('cluster p-value < 0.05',))
    else:
        plt.legend(fontsize=18)
plt.xlim(-0.2, 2)
ax.set_xlabel('Time (s)', labelpad=10, fontsize=23)
ax.set_ylabel('Tesla', fontsize=23)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.tick_params(axis='both', which='major', labelsize=18)
ax.tick_params(axis='both', which='minor', labelsize=18)
# ax.set_title(curr_label)
plt.savefig('/Users/boo/Desktop/fig2_roi_lr_b1.png')
plt.close()




