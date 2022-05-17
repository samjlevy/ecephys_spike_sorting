from argschema import ArgSchemaParser
import os
import logging
import time
import pathlib


import numpy as np
import pandas as pd

from ...common.utils import load_kilosort_data, write_cluster_group_tsv, read_cluster_group_tsv
from ...common.utils import getFileVersion
from ...common.epoch import get_epochs_from_nwb_file

from .metrics import calculate_metrics


def calculate_quality_metrics(args):

    print('ecephys spike sorting: quality metrics module')

    start = time.time()
    
    include_pcs = args['quality_metrics_params']['include_pcs']
    
    # make usre we can write an output file
    
    output_file_args = args['cluster_metrics']['cluster_metrics_file']
    
    output_file, metrics_version = getFileVersion(output_file_args)

    print("kilosort_output_dir: ")
    print(args['directories']['kilosort_output_directory'])
    print("Loading data...")


    try:
        if include_pcs:
            spike_times, spike_clusters, spike_templates, amplitudes, templates, channel_map, \
            channel_pos, clusterIDs, cluster_quality, cluster_amplitude, pc_features, pc_feature_ind, template_features = \
                    load_kilosort_data(args['directories']['kilosort_output_directory'], \
                        args['ephys_params']['sample_rate'], \
                        use_master_clock = False,
                        include_pcs = include_pcs)
        else:
            spike_times, spike_clusters, spike_templates, amplitudes, templates, channel_map, \
            channel_pos, clusterIDs, cluster_quality, cluster_amplitude = \
            load_kilosort_data(args['directories']['kilosort_output_directory'], \
                        args['ephys_params']['sample_rate'], \
                        use_master_clock = False,
                        include_pcs = include_pcs)
            pc_features = []
            pc_feature_ind = []
                    
        metrics = calculate_metrics(spike_times, spike_clusters, spike_templates, amplitudes, channel_map, channel_pos, templates, pc_features, pc_feature_ind, args['quality_metrics_params'])

    except FileNotFoundError:
        
        execution_time = time.time() - start

        print(" Files not available.")

        return {"execution_time" : execution_time,
            "quality_metrics_output_file" : None} 

    
    # build name for waveform_metrics file with matched version
    wm_args = args['waveform_metrics']['waveform_metrics_file']
    if metrics_version == 0:
        wm = wm_args
    else:
        # buld name for waveform metrics file with matched version
        wm = os.path.join( pathlib.Path(wm_args).parent, pathlib.Path(wm_args).stem + '_' + repr(metrics_version) + '.csv' )
    if os.path.exists(wm):
        metrics = metrics.merge(pd.read_csv(wm, index_col=0),
                     on='cluster_id',
                     suffixes=('_quality_metrics','_waveform_metrics'))

    print("Saving data...")
   
    metrics.to_csv(output_file, index=False )
    
    # EAJ addition 16 May 2022
    print("Labeling high isi_viol and contam_rate clusters as MUA...")
    # read current cluster assignments
    ci_tmp, cluster_group = read_cluster_group_tsv(os.path.join(args['directories']['kilosort_output_directory'], \
                'cluster_KSLabel.tsv'))
    
    # filter for MUA-like metrics
    is_mua = ((metrics['isi_viol']>0.2) | (metrics['contam_rate']>15))
    
    # change the labels accordingly
    labels = [ ]
    for i, ci in enumerate(ci_tmp):
    	if cluster_group[i]=='good' and is_mua[clusterIDs==ci]:
    		labels.append('mua')
    	else:
    		labels.append(cluster_group[i])
    write_cluster_group_tsv(ci_tmp, 
    						labels, 
    						args['directories']['kilosort_output_directory'], 
    						args['ephys_params']['cluster_group_file_name'])
    
    execution_time = time.time() - start
    print('total time: ' + str(np.around(execution_time,2)) + ' seconds')
    print()
    
    return {"execution_time" : execution_time,
            "quality_metrics_output_file" : output_file} # output manifest


def main():

    from ._schemas import InputParameters, OutputParameters

    mod = ArgSchemaParser(schema_type=InputParameters,
                          output_schema_type=OutputParameters)

    output = calculate_quality_metrics(mod.args)

    output.update({"input_parameters": mod.args})
    if "output_json" in mod.args:
        mod.output(output, indent=2)
    else:
        print(mod.get_output_json(output))


if __name__ == "__main__":
    main()
