# reset -f
# os.getcwd()
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
from mne.minimum_norm import apply_inverse_epochs, read_inverse_operator

label_list = ['f_c', 'f_nc', 't_c', 't_nc', 't_fpc', 't_tpc', 't_sc', 't_dc_b', 't_dc_l', 't_dc_r']  #
src_fs = mne.read_source_spaces('/Users/boo/Desktop/MEG_data_script/PreProcessed_data/fsaverage-src.fif')
fsave_vertices = [s['vertno'] for s in src_fs]

list_source_data_f_c = []
list_source_data_f_nc = []
list_source_data_t_c = []
list_source_data_t_nc = []
list_source_data_t_fpc = []
list_source_data_t_tpc = []
list_source_data_t_sc = []
list_source_data_t_dc_b = []
list_source_data_t_dc_l = []
list_source_data_t_dc_r = []

for ith_sub in list(range(2, 14)):

    temp_data_array = "/Users/boo/Desktop/MEG_data_script/PreProcessed_data/artefact_removed_sub" + str(
        ith_sub) + "_raw_100hz_sfre_200.fif"
    temp_event_array = "/Users/boo/Desktop/MEG_data_script/PreProcessed_data/events_post_resample_sub" + str(
        ith_sub) + "_100hz_sfre_200.npy"

    array_data = read_raw_fif(temp_data_array)
    array_event = numpy.load(temp_event_array)

    # pick channel
    all_chan = array_data.ch_names
    picks_mag = mne.pick_types(array_data.info, meg='mag')
    meg_channel = itemgetter(*picks_mag)(all_chan)
    meg_channel = meg_channel[29:301]

    # Compute epochs
    # baseline = None
    baseline = (-0.2, 0)
    tmin = -0.2
    tmax = 2.2
    event_id_f_c = {'f_sc': 50, 'f_dc': 40, 'f_tpc': 60}
    event_id_f_nc = {'f_fpc': 70}
    event_id_t_c = {'t_sc': 90, 't_dc_b': 80, 't_dc_l': 82, 't_dc_r': 84}
    event_id_t_nc = {'t_fpc': 110, 't_tpc': 100}
    event_id_t_fpc = {'t_fpc': 110}
    event_id_t_tpc = {'t_tpc': 100}
    event_id_t_sc = {'t_sc': 90}
    event_id_t_dc_b = {'t_dc_b': 80}
    event_id_t_dc_l = {'t_dc_l': 82}
    event_id_t_dc_r = {'t_dc_r': 84}

    epochs_f_c = mne.Epochs(array_data, array_event, picks=meg_channel, event_id=event_id_f_c, tmin=tmin, tmax=tmax,
                            baseline=baseline)
    epochs_f_nc = mne.Epochs(array_data, array_event, picks=meg_channel, event_id=event_id_f_nc, tmin=tmin, tmax=tmax,
                             baseline=baseline)
    epochs_t_c = mne.Epochs(array_data, array_event, picks=meg_channel, event_id=event_id_t_c, tmin=tmin, tmax=tmax,
                            baseline=baseline)
    epochs_t_nc = mne.Epochs(array_data, array_event, picks=meg_channel, event_id=event_id_t_nc, tmin=tmin, tmax=tmax,
                             baseline=baseline)
    epochs_t_fpc = mne.Epochs(array_data, array_event, picks=meg_channel, event_id=event_id_t_fpc, tmin=tmin, tmax=tmax,
                              baseline=baseline)
    epochs_t_tpc = mne.Epochs(array_data, array_event, picks=meg_channel, event_id=event_id_t_tpc, tmin=tmin, tmax=tmax,
                              baseline=baseline)
    epochs_t_sc = mne.Epochs(array_data, array_event, picks=meg_channel, event_id=event_id_t_sc, tmin=tmin, tmax=tmax,
                             baseline=baseline)
    epochs_t_dc_b = mne.Epochs(array_data, array_event, picks=meg_channel, event_id=event_id_t_dc_b, tmin=tmin,
                               tmax=tmax, baseline=baseline)
    epochs_t_dc_l = mne.Epochs(array_data, array_event, picks=meg_channel, event_id=event_id_t_dc_l, tmin=tmin,
                               tmax=tmax, baseline=baseline)
    epochs_t_dc_r = mne.Epochs(array_data, array_event, picks=meg_channel, event_id=event_id_t_dc_r, tmin=tmin,
                               tmax=tmax, baseline=baseline)

    epochs_f_c.load_data()
    epochs_f_nc.load_data()
    epochs_t_c.load_data()
    epochs_t_nc.load_data()
    epochs_t_fpc.load_data()
    epochs_t_tpc.load_data()
    epochs_t_sc.load_data()
    epochs_t_dc_b.load_data()
    epochs_t_dc_l.load_data()
    epochs_t_dc_r.load_data()

    evoke_f_c = epochs_f_c.average()
    evoke_f_nc = epochs_f_nc.average()
    evoke_t_c = epochs_t_c.average()
    evoke_t_nc = epochs_t_nc.average()
    evoke_t_fpc = epochs_t_fpc.average()
    evoke_t_tpc = epochs_t_tpc.average()
    evoke_t_sc = epochs_t_sc.average()
    evoke_t_dc_b = epochs_t_dc_b.average()
    evoke_t_dc_l = epochs_t_dc_l.average()
    evoke_t_dc_r = epochs_t_dc_r.average()

    epoch_list = [epochs_f_c, epochs_f_nc, epochs_t_c, epochs_t_nc,
                  epochs_t_fpc, epochs_t_tpc,
                  epochs_t_sc, epochs_t_dc_b, epochs_t_dc_l, epochs_t_dc_r]
    evoke_list = [evoke_f_c, evoke_f_nc, evoke_t_c, evoke_t_nc,
                  evoke_t_fpc, evoke_t_tpc,
                  evoke_t_sc, evoke_t_dc_b, evoke_t_dc_l, evoke_t_dc_r]
    output_data_list_by_evoke = [list_source_data_f_c, list_source_data_f_nc, list_source_data_t_c,
                                 list_source_data_t_nc,
                                 list_source_data_t_fpc, list_source_data_t_tpc,
                                 list_source_data_t_sc, list_source_data_t_dc_b, list_source_data_t_dc_l,
                                 list_source_data_t_dc_r]

    for ith_condition in list(range(0, len(epoch_list))):
        epoch_temp = epoch_list[ith_condition]
        evoke_temp = evoke_list[ith_condition]

        # read covariance
        if ith_condition == 0 or ith_condition == 1:
            cov = mne.read_cov('/Users/boo/Desktop/MEG_data_script/PreProcessed_data/covariance_sub' + str(
                ith_sub) + '_facingp-cov_100hz_sfre_200.fif')
        else:
            cov = mne.read_cov('/Users/boo/Desktop/MEG_data_script/PreProcessed_data/covariance_sub' + str(
                ith_sub) + '_targetingp-cov_100hz_sfre_200.fif')

        # compute stc
        fwd = mne.read_forward_solution(
            '/Users/boo/Desktop/MEG_data_script/PreProcessed_data/sub' + str(ith_sub) + '-fwd_100hz_sfre_200.fif')
        inv = mne.minimum_norm.make_inverse_operator(evoke_temp.info, fwd, cov)  # inverse operator
        snr = 3.0
        lambda2 = 1.0 / snr ** 2
        stc_evoke = apply_inverse(evoke_temp, inv, lambda2, 'dSPM')

        morph_evoke = mne.compute_source_morph(stc_evoke, 'sub' + str(ith_sub), 'fsaverage', spacing=fsave_vertices,
                                               subjects_dir='/Applications/freesurfer/subjects')
        stcs_evoke_morph = morph_evoke.apply(stc_evoke)
        # stcs_morph.save('/Users/boo/Desktop/MEG_data_script/analysis_source_result/stc_template')

        # each subjects' stc added in to list
        output_data_list_by_evoke[ith_condition].append(stcs_evoke_morph.data)
        source_vertices = stcs_evoke_morph.vertices
        source_times = stcs_evoke_morph.times

mean_source_array_evoke = np.asarray(output_data_list_by_evoke)

numpy.save('/Users/boo/Desktop/MEG_data_script/analysis_source_result/source_array_evoke_100hz_sfre_200.npy', mean_source_array_evoke)
numpy.save('/Users/boo/Desktop/MEG_data_script/analysis_source_result/source_vertices_100hz_sfre_200.npy', source_vertices)



