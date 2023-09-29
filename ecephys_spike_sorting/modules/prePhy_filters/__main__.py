"""
Additional module written by Emily Aery Jones, 19 May 2022
"""
from argschema import ArgSchemaParser
import os
import time
import numpy as np
import pandas as pd

from ...common.utils import write_cluster_group_tsv, getFileVersion


def filter_by_metrics(args):

    print('ecephys spike sorting: pre-Phy filters module')

    start = time.time()

    # find the correct metrics file (eg metrics.csv, metrics_1.csv, ...)
    metrics_file_args = args['cluster_metrics']['cluster_metrics_file']
    metrics_file, version = getFileVersion(metrics_file_args)
        
    waveform_metrics_file_args = args['waveform_metrics']['waveform_metrics_file']
    waveform_metrics_file, version = getFileVersion(waveform_metrics_file_args)
    
    # read the metrics files and join with waveforms if previous module failed to
    metrics = pd.read_csv(metrics_file)
    if 'snr' not in metrics.columns:
        waveform_metrics = pd.read_csv(waveform_metrics_file)
        metrics = metrics.merge(waveform_metrics, left_on='cluster_id', right_on='cluster_id')
    
    # read cluster assignments
    clusters_file = os.path.join(args['directories']['kilosort_output_directory'], \
                'cluster_KSLabel.tsv')
    #clusters_file = os.path.join(args['directories']['kilosort_output_directory'], \
    #            args['ephys_params']['cluster_group_file_name'])
    clusters = pd.read_csv(clusters_file, sep='\t')
    
    # save a copy of the original labels
    clusters_original = clusters.rename(columns={'group':'original_group'})
    clusters_original.to_csv(clusters_file[:-4]+'_original.tsv', sep='\t', index=False)
    
    # change the labels that meet MUA filters
    labels = [ ]
    good_clusters = 0
    mua_clusters = 0
    noise_clusters = 0
    for i, row in metrics.iterrows():
        # start with original label; overwrite 
        #label = clusters['KSLabel'].loc[clusters['cluster_id']==row['cluster_id']].item()
        #label = clusters['group'].loc[clusters['cluster_id']==row['cluster_id']].item()
        filters = args['prephy_filters_params']
        
        # find noise clusters
        if (row['snr']<filters['snr_min']) | (row['firing_rate']<filters['fr_min']) | \
            (row['halfwidth']>=filters['halfwidth_max']):
            label = 'noise'
            noise_clusters += 1
        elif (row['isi_viol']>filters['isi_viol_max']) & (row['num_viol']>filters['n_viol_max']):
            label = 'mua'
            mua_clusters += 1
        else:
            label = 'good'
            good_clusters += 1
        labels.append(label)

    # write output
    write_cluster_group_tsv(metrics['cluster_id'], 
    						labels, 
    						args['directories']['kilosort_output_directory'], 
    						args['ephys_params']['cluster_group_file_name'])
    print(f'Clusters classified: {good_clusters} good, {mua_clusters} MUA, and {noise_clusters} noise from {len(metrics)} clusters')
    
    execution_time = time.time() - start
    print('total time: ' + str(np.around(execution_time,2)) + ' seconds')
    print()
    
    return {"execution_time" : execution_time,
            "quality_metrics_output_file" : args['cluster_metrics']['cluster_metrics_file']} # output manifest

def main():

    from ._schemas import InputParameters, OutputParameters

    mod = ArgSchemaParser(schema_type=InputParameters,
                          output_schema_type=OutputParameters)

    output = filter_by_metrics(mod.args)

    output.update({"input_parameters": mod.args})
    if "output_json" in mod.args:
        mod.output(output, indent=2)
    else:
        print(mod.get_output_json(output))


if __name__ == "__main__":
    main()