# source contrast get averaged
# reset -f
import numpy as np
import mne
from scipy import stats as stats
from mne.stats import permutation_t_test
from mne.stats import (spatio_temporal_cluster_1samp_test,
                       summarize_clusters_stc)
from sklearn.base import clone

# ['f_c', 'f_nc', 't_c', 't_sc', 't_dc_b', 't_dc_l', 't_dc_r', 't_fpc', 't_tpc']  #

mean_stc_ndarray = np.load(
    '/Users/boo/Desktop/MEG_data_script/analysis_source_result/source_data_ndarray_100hz_sfre_100.npy',
    allow_pickle=True)
vertice = np.load('/Users/boo/Desktop/MEG_data_script/analysis_source_result/source_vertices.npy')
stc_template = mne.read_source_estimate('/Users/boo/Desktop/MEG_data_script/analysis_source_result/stc_template-rh.stc')
stc_template.subject = 'fsaverage'

# ['f_c', 'f_nc', 't_c', 't_sc', 't_dc_b', 't_dc_l', 't_dc_r', 't_fpc', 't_tpc', 't_lrb', 't_control']  #

segm = 10
start_tp_index = 0
end_tp_index = 120
start_tp = -0.2

data_back = mean_stc_ndarray[4]
data_left = mean_stc_ndarray[5]
data_right = mean_stc_ndarray[6]
data_t_fpc = mean_stc_ndarray[7]
data_t_tpc = mean_stc_ndarray[8]

# T_obs, clusters, cluster_p_values, H0 = spatio_temporal_cluster_1samp_test(data_back - (data_left+data_right)/2, n_jobs=10, threshold=1.79, buffer_size=None, verbose=True, tail=1)
# sig_list = np.where(cluster_p_values <= 0.05)[0]

resampled_ima_b = np.zeros([data_back.shape[0], data_back.shape[1], len(np.arange(start_tp_index, end_tp_index, segm))])
resampled_ima_l = np.zeros([data_left.shape[0], data_left.shape[1], len(np.arange(start_tp_index, end_tp_index, segm))])
resampled_ima_r = np.zeros(
    [data_right.shape[0], data_right.shape[1], len(np.arange(start_tp_index, end_tp_index, segm))])
resampled_ima_fpc = np.zeros(
    [data_t_fpc.shape[0], data_t_fpc.shape[1], len(np.arange(start_tp_index, end_tp_index, segm))])
resampled_ima_tpc = np.zeros(
    [data_t_tpc.shape[0], data_t_tpc.shape[1], len(np.arange(start_tp_index, end_tp_index, segm))])
for ind, ith_ts in enumerate(np.arange(start_tp_index, end_tp_index, segm)):
    print(range(ith_ts, ith_ts + segm))
    resampled_ima_b[..., ind] = np.mean(data_back[..., range(ith_ts, ith_ts + segm)], axis=2)
    resampled_ima_l[..., ind] = np.mean(data_left[..., range(ith_ts, ith_ts + segm)], axis=2)
    resampled_ima_r[..., ind] = np.mean(data_right[..., range(ith_ts, ith_ts + segm)], axis=2)
    resampled_ima_fpc[..., ind] = np.mean(data_t_fpc[..., range(ith_ts, ith_ts + segm)], axis=2)
    resampled_ima_tpc[..., ind] = np.mean(data_t_tpc[..., range(ith_ts, ith_ts + segm)], axis=2)
mean_b = np.mean(resampled_ima_b, axis=0)
mean_l = np.mean(resampled_ima_l, axis=0)
mean_r = np.mean(resampled_ima_r, axis=0)
mean_fpc = np.mean(resampled_ima_fpc, axis=0)
mean_tpc = np.mean(resampled_ima_tpc, axis=0)
mean_control = (mean_fpc + mean_tpc) / 2

pool_name = ['back', 'left', 'right']  #
mean_b = mean_b - mean_control
mean_l = mean_l - mean_control
mean_r = mean_r - mean_control

mean_b[mean_b < 0] = 0
mean_l[mean_l < 0] = 0
mean_r[mean_r < 0] = 0

pool_data = [mean_b, mean_l, mean_r]

# tval_lr_b = (mean_lr - mean_b) / (np.std((resampled_ima_l+resampled_ima_r)/2 - resampled_ima_b, axis=0)/np.sqrt(12))
# the_data = tval_lr_b

for ith in list(range(0, len(pool_name))):

    stc_template.data = pool_data[ith]
    stc_template._tmin = start_tp
    stc_template._tstep = segm / 100
    stc_template._times = [round(x, 2) for x in list(np.arange(-0.2, 1, segm / 100))]

    # Supported views: ['lat', 'med', 'ros', 'cau', 'dor' 'ven', 'fro', 'par']
    views = ['med']

    # for tp in [round(x, 2) for x in list(np.arange(-0.2, 1, segm / 100))]:
    tp = 0.5  # -0.2 -0.1 0 0.1  0.2 0.3 0.4 0.5
    for view in views:
        surfer_kwargs = dict(
            hemi='split', subjects_dir='/Applications/freesurfer/subjects/', views=view, background='white',
            foreground='black',
            time_label=None, colorbar=None,
            clim=dict(kind='percent', pos_lims=(95, 97, 99)), colormap='jet',  # pos_lims=(75, 85, 95)
            # clim=dict(kind='value', pos_lims=(2.2, 3.1, 4.4)), colormap='seismic',
            initial_time=round(tp, 2), time_unit='s', size=(400, 200))
        brain = stc_template.plot(**surfer_kwargs)
        brain.save_image(
            '/Users/boo/Desktop/MEG_data_script/analysis_source_result/source_activity_ego_control_0.1/' +
            pool_name[
                ith] + '_' + str(
                tp) + '_' + view + '.png')
        brain.close()



#### combine plot
import os
import matplotlib.pyplot as plt
import numpy
import matplotlib.image as mpimg
import numpy as np
from matplotlib.pyplot import figure

# source_activity_condition_contrast source_activity_condition

# pool_name = ['back', 'left', 'right']  #
# for ith in list(range(0, len(pool_name))):

view_pool = ['medial', 'lateral']

for view in view_pool:

    files = []
    for r, d, f in os.walk(
            '/Users/boo/Desktop/MEG_data_script/analysis_source_result/source_activity_ego_0.2_only_positive/'):
        for file in f:
            if '.png' in file:
                files.append(os.path.join(r, file))
    c1 = []
    # for ith_name in files:
    #     if pool_name[ith] in ith_name and view in ith_name:
    #         c1.append(ith_name)
    for ith_name in files:
        if view in ith_name:
            c1.append(ith_name)
    c1.sort()
    data_pool = c1
    # tp = ['-0.2 ~ -0s', '0 ~ 0.2s', '0.2 ~ 0.4s', '0.4 ~ 0.6s', '0.6 ~ 0.8s', '0.8 ~ 1.0s', '1.0 ~ 1.2s', '1.2 ~ 1.4s',
    #       '1.4 ~ 1.6s', '1.6 ~ 1.8s']
    tp = ['0 ~ 0.2s', '0.2 ~ 0.4s', '0.4 ~ 0.6s', '0.6 ~ 0.8s', '0.8 ~ 1.0s']
    num_row = 5
    num_col = 1
    f, axarr = plt.subplots(num_row, num_col, figsize=(6, 15))
    counter = 0
    # for ith_col in list(range(0, 1)):
    #     for ith_row in list(range(0, 6)):
    #         axarr[ith_row, ith_col].imshow(mpimg.imread(data_pool[counter]))
    #         axarr[ith_row, ith_col].set_title(tp[counter])
    #         plt.setp(axarr[ith_row, ith_col].get_xticklabels(), visible=False)
    #         plt.setp(axarr[ith_row, ith_col].get_yticklabels(), visible=False)
    #         counter = counter + 1

    for ith_row in list(range(0, 5)):
        axarr[ith_row].imshow(mpimg.imread(data_pool[counter]))
        # axarr[ith_row].set_title(tp[counter])
        axarr[ith_row].annotate(tp[counter], xy=(0, 0.9), xycoords="axes fraction", fontsize=20)
        axarr[ith_row].set_axis_off()
        plt.setp(axarr[ith_row].get_xticklabels(), visible=False)
        plt.setp(axarr[ith_row].get_yticklabels(), visible=False)
        counter = counter + 1
    plt.show()
    plt.savefig('/Users/boo/Desktop/' + pool_name[ith] + '_' + view + '_control_corrected' + '.png',
                bbox_inches='tight')
