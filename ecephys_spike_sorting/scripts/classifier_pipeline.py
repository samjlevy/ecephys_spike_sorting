import os, io, json
import subprocess

from helpers import SpikeGLX_utils
from create_input_json import createInputJson

# run a set of SpikeGLX tcat.probeN.bin files that are stored in one folder.
# creates an output folder for each, generatees a channel map file from
# the SpikeGLX metadata, then runs any other listed modules.

# directory for json files -- these record the parameters used for processing
#json_directory = r'F:\CatGT\catgt_AA_200920_4_201005_dark_g0\AA_200920_4_201005_dark_g0_imec0'

# directory with the raw data files. The metadata should be present, also
#npx_directory = r'F:\CatGT\catgt_AA_200920_4_201005_dark_g0\AA_200920_4_201005_dark_g0_imec0'


# list of run names
import glob
import os
 
files = glob.glob(r'F:\CatGT\cat*\*\imec0_ks2\waveform_metrics.csv')
run_names = [os.path.split(fi)[0] for fi in files]
run_names = [x for x in run_names if not os.path.isfile(os.path.join(x,'classifier_cluster_heuristic.txt'))]
#run_names = glob.glob('F:\\CatGT\\cat*dark*\\*\\imec0_ks_9_2')
[print(na) for na in run_names]
probe_type = 'NP1'

#exit()


# List of modules to run per probe
# if not running kilosort_helper, KS2 output must be in directories
# named according to this script, i.e. run_name_gN_tcat.imecN_phy
modules = [
            #'depth_estimation',
			#'kilosort_helper',
            #'kilosort_postprocessing',
            'noise_templates',
            #'psth_events',
            #'mean_waveforms',
            #'quality_metrics'

		  ]

for fi in run_names:
    
    [h,t] = os.path.split(fi)
    npx_directory = h
    
    npx_file = glob.glob(os.path.join(npx_directory,'*.ap.bin'))
    print(npx_file[0])
    [h,t]=os.path.split(h)
    session_id = t
    baseName = SpikeGLX_utils.ParseTcatName(npx_file[0])

    print(session_id)




    input_json = os.path.join(fi, session_id +'-inputClassifier.json')
    output_json = os.path.join(fi, session_id +'-outputClassifier.json')

    
   
    
    print( 'Creating json file postprocessing')
    info = createInputJson(input_json, npx_directory=npx_directory, 
	                                   continuous_file = npx_file[0],
									   kilosort_output_directory=fi, 
                                       ks_make_copy = False,
                                       noise_template_use_rf = False,
                                       cluster_group_file_name = 'classifier_cluster_heuristic.txt',
                                       extracted_data_directory = npx_directory

                                       )
  

    print(info['noise_waveform_params']['use_preclustered'])
    for module in modules:
        command = "python -W ignore -m ecephys_spike_sorting.modules." + module + " --input_json " + input_json \
		          + " --output_json " + output_json
        subprocess.check_call(command.split(' '))