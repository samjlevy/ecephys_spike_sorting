import os
import subprocess
import glob

from ecephys_spike_sorting.scripts.helpers import SpikeGLX_utils
from ecephys_spike_sorting.scripts.helpers import log_from_json
#from ecephys_spike_sorting.scripts.helpers import run_one_probe
from ecephys_spike_sorting.scripts.create_input_json import createInputJson


# script to run CatGT, KS2, postprocessing and TPrime on data collected using
# SpikeGLX. The construction of the paths assumes data was saved with
# "Folder per probe" selected (probes stored in separate folders) AND
# that CatGT is run with the -out_prb_fld option

# -------------------------------
# -------------------------------
# User input -- Edit this section
# -------------------------------
# -------------------------------

animal = 'StickPin'
rec_file_stem = '20211022_StickPin_DY01'
npx_directory = os.path.join('D:/',animal)
# all output will be written here.
# Output will be in the standard SpikeGLX directory structure:
# run_folder/probe_folder/*.bin
catGT_dest = os.path.join('D:/Sorted',animal, rec_file_stem)

if not os.path.exists(catGT_dest):
    os.mkdir(catGT_dest)

# -----------
# Input data
# -----------
# Name for log file for this pipeline run. Log file will be saved in the
# output destination directory catGT_dest
# If this file exists, new run data is appended to it
logName = f'{rec_file_stem}_log.csv'



# run_specs = name, gate, trigger and probes to process
# Each run_spec is a list of 4 strings:
#   undecorated run name (no g/t specifier, the run field in CatGT)
#   gate index range, as a string (e.g. '0,9'); can be greater than actual range
#       later code will look for all g indices in this range with the same filename stem
#       EAJ edit: can't currently handle gate indices with >1 digit or multiple triggers
#   triggers to process/concatenate, as a string e.g. '0,400', '0,0 for a single file
#           can replace first limit with 'start', last with 'end'; 'start,end'
#           will concatenate all trials in the probe folder
#   probes to process, as a string, e.g. '0', '0,3', '0:3'
#   brain regions, list of strings, one per probe, to set region specific params
#           these strings must match a key in the param dictionaries above.

run_specs = [									
						[rec_file_stem, '0,9', '0,0', '0', ['cortex'] ]
]

# ------------
# CatGT params
# ------------
run_CatGT = True   # set to False to sort/process previously processed data.


# CAR mode for CatGT. Must be equal to 'None', 'gblcar' or 'loccar'
car_mode = 'None'


# CatGT commands for bandpass filtering, artifact correction, and zero filling
# Note 1: directory naming in this script requires -prb_fld and -out_prb_fld
# Note 2: this command line includes specification of edge extraction
# see CatGT readme for details
# these parameters will be used for all runs

# notes from Mari's repo:
# gfix=0,0.10,0.02 -- artifact removal; params: |thresh_amp(mV)|,|slope(mV/sample)|,noise
# -t_miss_ok option required to concatenate over missing g or t indices
# -zerofillmax=500 option required to fill gaps only up to 500ms of zeros,
# so kilsort doesn't crash
catGT_cmd_string = '-t_miss_ok -zerofillmax=500 -prb_fld -out_prb_fld -gfix=0.4,0.10,0.02 -lffilter=butter,12,1,500 '


json_directory = catGT_dest

# -----------------------
# -----------------------
# End of user input
# -----------------------
# -----------------------

# delete the existing CatGT.log
try:
    os.remove('CatGT.log')
except OSError:
    pass

# delete existing Tprime.log
try:
    os.remove('Tprime.log')
except OSError:
    pass

# delete existing C_waves.log
try:
    os.remove('C_Waves.log')
except OSError:
    pass

# check for existence of log file, create if not there
logFullPath = os.path.join(catGT_dest, logName)
if not os.path.isfile(logFullPath):
    # create the log file, write header
    log_from_json.writeHeader(logFullPath)


for spec in run_specs:

    session_id = spec[0]

    
    # Make list of probes from the probe string
    prb_list = SpikeGLX_utils.ParseProbeStr(spec[3])
    
    # build path to the first probe folder; look into that folder
    # to determine the range of trials if the user specified t limits as
    # start and end
    first_gate = spec[1][0]
    catGT_internal_dest = os.path.join(catGT_dest,f'catgt_{rec_file_stem}_g{first_gate}')
    if not os.path.exists(catGT_internal_dest):
        os.mkdir(catGT_internal_dest)
    
    run_folder_name = spec[0] + '_g' + first_gate
    prb0_fld_name = run_folder_name + '_imec' + prb_list[0]
    prb0_fld = os.path.join(npx_directory, run_folder_name, prb0_fld_name)
    first_trig, last_trig = SpikeGLX_utils.ParseTrigStr(spec[2], prb_list[0], first_gate, prb0_fld)
    trigger_str = repr(first_trig) + ',' + repr(last_trig)
    
    # from MS fork
    # get list of g-indices to concatenate from data directory
    g_range = '[' + spec[1][0] + '-' + spec[1][-1] + ']'
    g_tocat = sorted(glob.glob(os.path.join(npx_directory,(rec_file_stem + '_g' + g_range))))
    glist = ''.join((x[-1]+'-') for x in g_tocat)[:-1] # g inds separated by dashes, minus the last dash        

    print('Concatenating g indices ' + glist)
    
    # loop over all probes to build json files of input parameters
    # initalize lists for input and output json files
    catGT_input_json = []
    catGT_output_json = []
    module_input_json = []
    module_output_json = []
    session_id = []
    catgt_output_dir = []
    data_directory = []
    
    # first loop over probes creates json files containing parameters for
    # both preprocessing (CatGt) and sorting + postprocessing
    
    for i, prb in enumerate(prb_list):
            
        #create CatGT command for this probe
        print('Creating json file for CatGT on probe: ' + prb)
        # Run CatGT
        catGT_input_json.append(os.path.join(json_directory, spec[0] + '_g' + glist + '_prb' + prb + '_CatGT' + '-input.json'))
        catGT_output_json.append(os.path.join(json_directory, spec[0] + '_g' + glist + '_prb' + prb + '_CatGT' + '-output.json'))
        
        # build extract string for SYNC channel for this probe
        sync_extract = '-SY=' + prb +',-1,6,500'
        
        # if this is the first probe proceessed, process the ni stream with it
        catGT_stream_string = '-lf'
        extract_string = sync_extract
        
        # build name of first trial to be concatenated/processed;
        # allows reading of the metadata
        run_str = spec[0] + '_g' + first_gate
        run_folder = run_str
        prb_folder = run_str + '_imec' + prb
        input_data_directory = os.path.join(npx_directory, run_folder, prb_folder)
        fileName = run_str + '_t' + repr(first_trig) + '.imec' + prb + '.lf.bin'
        continuous_file = os.path.join(input_data_directory, fileName)
        metaName = run_str + '_t' + repr(first_trig) + '.imec' + prb + '.lf.meta'
        input_meta_fullpath = os.path.join(input_data_directory, metaName)
        
        print(input_meta_fullpath)
         
        info = createInputJson(catGT_input_json[i], npx_directory=npx_directory, 
                                       continuous_file = continuous_file,
                                       kilosort_output_directory=catGT_dest,
                                       spikeGLX_data = True,
                                       input_meta_path = input_meta_fullpath,
                                       catGT_run_name = spec[0],
                                       gate_string = spec[1],
                                       trigger_string = trigger_str,
                                       probe_string = prb,
                                       catGT_stream_string = catGT_stream_string,
                                       catGT_car_mode = car_mode,
                                       catGT_cmd_string = catGT_cmd_string + ' ' + extract_string,
                                       extracted_data_directory = catGT_dest
                                       )      
        
        # from MS fork
        if run_CatGT:
            command = "python -W ignore -m ecephys_spike_sorting.modules." + 'catGT_helper' + " --input_json " + catGT_input_json[i] \
            	          + " --output_json " + catGT_output_json[i]
            subprocess.check_call(command.split(' '))           

