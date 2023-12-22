from argschema import ArgSchemaParser
import os
import sys
import subprocess
import time
import pathlib

import numpy as np
import pandas as pd
from scipy.io import loadmat

from ...common.utils import load_kilosort_data, write_cluster_group_tsv, read_cluster_group_tsv
from ...common.utils import getSortResults
from ...common.utils import getFileVersion

from .extract_waveforms import extract_waveforms, writeDataAsNpy
from .waveform_metrics import calculate_waveform_metrics
from .metrics_from_file import metrics_from_file

def calculate_mean_waveforms(args):

    print('ecephys spike sorting: mean waveforms module')
    
    start = time.time()
    
    if args['mean_waveform_params']['use_C_Waves']:
        
        print('Calculating mean waveforms using C_waves.')
        spikeglx_bin = args['ephys_params']['ap_band_file']
        # regenerate the clus_Table in case there has been manual curation of the data in phy
        output_dir = args['directories']['kilosort_output_directory']
        
        # get version number for new clus_table file
        clu_path_orig = os.path.join(output_dir, 'clus_Table.npy' )
        clus_table_npy, clu_version = getFileVersion(clu_path_orig)
        
        #version = 0 if no clu_Table exists, file = clus_Table.npy
        #version = 1 or higher, new clus_Table = clus_Table_version.npy
        
        getSortResults(output_dir, clu_version)
        
        # build paths to cluster and times tables, which are generated by
        # kilosort_helper module
        clus_time_npy = os.path.join(output_dir, 'spike_times.npy' )
        clus_lbl_npy = os.path.join(output_dir, 'spike_clusters.npy' )
        dest, wavefile = os.path.split(args['mean_waveform_params']['mean_waveforms_file'])
        
        # on first call to mean_waveforms, output has no version indicator;
        # for later calls, will rename to _version.npy
        # need to rename the originals, because the output names from C_Waves are hard coded

        if clu_version == 1:           
            # if mean_waveforms files exists, rename
            old_mwf = os.path.join(dest,'mean_waveforms.npy')
            if os.path.exists(old_mwf):
                new_mwf = os.path.join(dest,'mean_waveforms_0.npy')
                os.rename(old_mwf, new_mwf)
            old_snr = os.path.join(dest,'cluster_snr.npy')
            if os.path.exists(old_snr):
                new_snr = os.path.join(dest,'cluster_snr_0.npy')
                os.rename(old_snr, new_snr)

        
        # kilosort saves the spike_clusters files as uint32. 
        # when phy re-saves after curation, it saves as int32 (!)
        # to ensure the correct datatype for C_Waves, load the spike_clusters
        # and convert if necessary
        sc = np.load(clus_lbl_npy)
        if sc.dtype != 'uint32':
            sc = sc.astype('uint32')
            np.save(clus_lbl_npy,sc)
        
        
        # path to the 'runit.bat' executable that calls C_Waves.
        # Essential in linux where C_Waves executable is only callable through runit
        if sys.platform.startswith('win'):
            exe_path = os.path.join(args['mean_waveform_params']['cWaves_path'], 'runit.bat')
        elif sys.platform.startswith('linux'):
            exe_path = os.path.join(args['mean_waveform_params']['cWaves_path'], 'runit.sh')
        else:
            print('unknown system, cannot run C_Waves')
        
        cwaves_cmd = exe_path + ' -spikeglx_bin=' + spikeglx_bin + \
                                ' -clus_table_npy=' + clus_table_npy + \
                                ' -clus_time_npy=' + clus_time_npy + \
                                ' -clus_lbl_npy=' + clus_lbl_npy + \
                                ' -dest=' + dest + \
                                ' -samples_per_spike=' + repr(args['mean_waveform_params']['samples_per_spike']) + \
                                ' -pre_samples=' + repr(args['mean_waveform_params']['pre_samples']) + \
                                ' -num_spikes=' + repr(args['mean_waveform_params']['spikes_per_epoch']) + \
                                ' -snr_radius=' + repr(args['mean_waveform_params']['snr_radius']) + \
                                ' -snr_radius_um=' + repr(args['mean_waveform_params']['snr_radius_um'])
                                
        print(cwaves_cmd)
        
        # make the C_Waves call
        subprocess.Popen(cwaves_cmd,shell='False').wait()
        
        # for first version, retain original names
        if clu_version == 0:
            mean_waveform_fullpath = os.path.join(dest, 'mean_waveforms.npy')
            snr_fullpath = os.path.join(dest, 'cluster_snr.npy')
        else:
            # build names with version number and rename
            # version 0 files are not renamed to maintain compatiblity with
            mean_waveform_fullpath = os.path.join(dest, 'mean_waveforms_' + repr(clu_version) + '.npy')
            snr_fullpath = os.path.join(dest, 'cluster_snr_' + repr(clu_version) + '.npy')
            os.rename(os.path.join(dest, 'mean_waveforms.npy'), mean_waveform_fullpath)
            os.rename(os.path.join(dest, 'cluster_snr.npy'), snr_fullpath)
            
        
        # C_Waves writes out files of the waveforms and snr
        # call version of calculate_waveform_metrics that will use these files
        
        # load in kilosort output needed for these calculations
        spike_times, spike_clusters, spike_templates, amplitudes, templates, channel_map, \
        channel_pos, clusterIDs, cluster_quality, cluster_amplitude = \
                load_kilosort_data(args['directories']['kilosort_output_directory'], \
                    args['ephys_params']['sample_rate'], \
                    convert_to_seconds = False)
                
        # read in inverse of whitening matrix
        w_inv = np.load((os.path.join(args['directories']['kilosort_output_directory'], 'whitening_mat_inv.npy')))
        
        # the channel_pos loaded from the phy output omits any sites excluded
        # as noise by the kilosort_helper module, or excluded fow low spike rete
        # by kilosort itself. The waveform metrics are calculated on ALL sites
        # based on the mean waveforms calculated for each unit; therefore
        # we need the site locations for all sites.
        # load the channel map associated with this kilosort run; in kilosort_helper
        # a copy is made next to the data file
        input_file = args['ephys_params']['ap_band_file']
        dat_dir, dat_fname = os.path.split(input_file)
        dat_name, dat_ext = os.path.splitext(dat_fname)
        chanMapMat = os.path.join(dat_dir, (dat_name +'_chanMap.mat'))
        site_x = np.squeeze(loadmat(chanMapMat)['xcoords'])
        site_y = np.squeeze(loadmat(chanMapMat)['ycoords'])
        
                

                
        metrics = metrics_from_file(mean_waveform_fullpath, snr_fullpath, clus_table_npy, \
                    spike_times, \
                    spike_clusters, \
                    templates, \
                    channel_map, \
                    args['ephys_params']['bit_volts'], \
                    args['ephys_params']['sample_rate'], \
                    args['ephys_params']['vertical_site_spacing'], \
                    w_inv, \
                    site_x, site_y, \
                    args['mean_waveform_params'])
        
        wm_fullpath = (args['waveform_metrics']['waveform_metrics_file'])

        if clu_version > 0:
           # save new metrics as _version number
           wm_fullpath = os.path.join(pathlib.Path(wm_fullpath).parent, pathlib.Path(wm_fullpath).stem + '_' + repr(clu_version) + '.csv')
    
        metrics.to_csv(wm_fullpath, index=False)
        
    else:
        
        print('Calculating mean waveforms using python.')
        print("Loading data...")
        
        # The waveform metrics are calculated on ALL sites
        # based on the mean waveforms calculated for each unit; therefore
        # we need the site locations for all sites.
        # load the channel map associated with this kilosort run; in kilosort_helper
        # a copy is made next to the data file
        input_file = args['ephys_params']['ap_band_file']
        dat_dir, dat_fname = os.path.split(input_file)
        dat_name, dat_ext = os.path.splitext(dat_fname)
        chanMapMat = os.path.join(dat_dir, (dat_name +'_chanMap.mat'))
        site_x = np.squeeze(loadmat(chanMapMat)['xcoords'])
        site_y = np.squeeze(loadmat(chanMapMat)['ycoords'])
    
        rawData = np.memmap(args['ephys_params']['ap_band_file'], dtype='int16', mode='r')
        data = np.reshape(rawData, (int(rawData.size/args['ephys_params']['num_channels']), args['ephys_params']['num_channels']))
    
        spike_times, spike_clusters, spike_templates, amplitudes, templates, channel_map, \
        channel_pos, clusterIDs, cluster_quality, cluster_amplitude = \
                load_kilosort_data(args['directories']['kilosort_output_directory'], \
                    args['ephys_params']['sample_rate'], \
                    convert_to_seconds = False)
    
        print("Calculating mean waveforms...")
    
        waveforms, spike_counts, coords, labels, metrics = extract_waveforms(data, spike_times, \
                    spike_clusters,
                    templates,
                    channel_map,
                    args['ephys_params']['bit_volts'], \
                    args['ephys_params']['sample_rate'], \
                    args['ephys_params']['vertical_site_spacing'], \
                    site_x, \
                    site_y, \
                    args['mean_waveform_params'])
    
        writeDataAsNpy(waveforms, args['mean_waveform_params']['mean_waveforms_file'])
        metrics.to_csv(args['waveform_metrics']['waveform_metrics_file'], index=False)


    # if the cluster metrics have already been run, merge the waveform metrics into that file
    # build file path with current version
    metrics_args = args['cluster_metrics']['cluster_metrics_file']
    metrics_curr = os.path.join(pathlib.Path(metrics_args).parent, pathlib.Path(metrics_args).stem + '_' + repr(clu_version) + '.csv')

    if os.path.exists(metrics_curr):
        qmetrics = pd.read_csv(metrics_curr)
        qmetrics = qmetrics.drop(qmetrics.columns[0], axis='columns')
        qmetrics = qmetrics.merge(pd.read_csv(wm_fullpath, index_col=0),
                     on='cluster_id',
                     suffixes=('_quality_metrics','_waveform_metrics'))  
        print("Saving merged quality metrics ...")
        qmetrics.to_csv(metrics_curr, index=False)
        
    execution_time = time.time() - start

    print('total time: ' + str(np.around(execution_time,2)) + ' seconds')
    print()
    
    return {"execution_time" : execution_time} # output manifest


def main():

    from ._schemas import InputParameters, OutputParameters

    mod = ArgSchemaParser(schema_type=InputParameters,
                          output_schema_type=OutputParameters)

    output = calculate_mean_waveforms(mod.args)

    output.update({"input_parameters": mod.args})
    if "output_json" in mod.args:
        mod.output(output, indent=2)
    else:
        print(mod.get_output_json(output))


if __name__ == "__main__":
    main()
