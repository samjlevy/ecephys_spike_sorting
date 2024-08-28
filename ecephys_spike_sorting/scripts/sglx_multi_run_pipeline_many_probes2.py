import os
import sys
import subprocess
import glob
import pandas as pd
import shutil
import timeit

#from datetime import datetime
#my_env = os.environ.copy()
#my_env["PATH"] = "/usr/sbin:/sbin:" + my_env["PATH"]

from ecephys_spike_sorting.scripts.helpers import SpikeGLX_utils
from ecephys_spike_sorting.scripts.helpers import log_from_json
#from ecephys_spike_sorting.scripts.helpers import run_one_probe
from ecephys_spike_sorting.scripts.create_input_json import createInputJson

from ecephys_spike_sorting.scripts.sglx_add_fix_lfp import fix_sglx_lfps

# to add:
#     - Run Sam's wavestats
#     - Run phy meanwaveforms
#     - Better way to signal server vs. local

run_from_oak = True # default otherwise to local

run_CatGT = True
run_KSandModules = True
run_TPrime = False
run_LFPfix = True

# ---------------
# Modules List
# ---------------
# List of modules to run per probe; CatGT and TPrime are called once for each run.
modules = [
            'kilosort_helper',
            'kilosort_postprocessing',
            'noise_templates',    
            'mean_waveforms',
            'quality_metrics',
            ####'depth_estimation',
            'prePhy_filters'
			]

# script to run CatGT, KS2, postprocessing and TPrime on data collected using
# SpikeGLX. The construction of the paths assumes data was saved with
# "Folder per probe" selected (probes stored in separate folders) AND
# that CatGT is run with the -out_prb_fld option

# -------------------------------
# -------------------------------
# User input -- Edit this section
# -------------------------------
# -------------------------------

# set for each run
local_base_dir = 'D:/to_process/process_tonight/'
oak_base_dir = '//oak-smb-giocomo.stanford.edu/groups/giocomo/samjlevy/Sam_NPX'
base_dir = local_base_dir
if run_from_oak:
    base_dir = oak_base_dir

prefix = 'all' #eg datetime.today().strftime('%Y%m%d')
rec_file_list = os.path.join(base_dir,'Preprocessed_Data/Provenance',prefix+'_sessions.csv')
raw_data_dir = os.path.join(base_dir);
processed_data_dir = os.path.join(base_dir,'Preprocessed_Data/Spikes')

log_file = os.path.join(base_dir, 'Preprocessed_Data/Provenance',prefix+'_ecephys_log.csv')
log_stream = open(log_file, 'a+')
log_stream.write('Index,File,Sort_Error,Sort_Error_Description\n')
sessions = pd.read_csv(rec_file_list)

# brain region specific params
# can add a new brain region by adding the key and value for each param
# can add new parameters -- any that are taken by create_input_json --
# by adding a new dictionary with entries for each region and setting the 
# according to the new dictionary in the loop to that created json files.
refPerMS_dict = {'default': 2.0, 'cortex': 2.0, 'medulla': 1.5, 'thalamus': 1.0}

# threhold values appropriate for KS2, KS2.5
#ksTh_dict = {'default':'[10,4]', 'cortex':'[10,4]', 'medulla':'[10,4]', 'thalamus':'[10,4]'}
# threshold values appropriate for KS3.0
ksTh_dict = {'default':'[9,9]', 'cortex':'[9,9]', 'medulla':'[9,9]', 'thalamus':'[9,9]'}

# ------------
# CatGT params
# ------------
# CAR mode for CatGT. Must be equal to 'None', 'gblcar' or 'loccar'
car_mode = 'gblcar'
lfp_car_mode = 'None'
# inner and outer radii, in um for local comman average reference, if used
loccar_min = 40
loccar_max = 160

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
catGT_cmd_string_onepointoh = '-t_miss_ok -zerofillmax=500 -prb_fld -out_prb_fld -gfix=0.4,0.10,0.02 -apfilter=butter,12,300,10000'
lfp_catGT_cmd_string = '-t_miss_ok -zerofillmax=500 -prb_fld -out_prb_fld -gfix=0.4,0.10,0.02 ' # -lffilter=butter,12,1,500
catGT_cmd_string_twopointoh = '-t_miss_ok -zerofillmax=500 -prb_fld -out_prb_fld -gfix=0.4,0.10,0.02 -apfilter=butter,12,300,10000 -lffilter=butter,12,1,500'
                    
ni_present = True

# ----- NIDAQ INPUTS -----
# -- Each XA gets its own word, starting with 0. XA inputs must come first.
# -- XD inputs come next, on separate words from the XA inputs. 
# -- Each XD word contains up to 16 bits (0:15)
# -- Yggdrasil inputs: --
# XA=0,1,3,500 -- sync channel on nidaq: word 0, thresh 1 V, must stay above 3V, dur 500 ms
# XA=1,1,3,0 -- Arduino: word 1, thresh 1V, must stay above 3V, dur ??? ms
# XA=2,1,3,0 -- camera: word 2, thresh 1V, must stay above 3V, dur ??? ms
# - duration must be within +/-20% of the actual pulse width; 0 ignores pulse width requirement
#ni_extract_string = '-XA=0,1,3,500 -XA=1,1,3,0 -XA=2,1,3,0'
ni_extract_string = ''
#ni_extract_string = '-xa=0,0,0,1,3,500 -xa=0,0,1,1,3,0 -xa=0,0,2,1,3,0' # comes from Emily's version

# ----------------------
# KS2 or KS25 parameters
# ----------------------
# parameters that will be constant for all recordings
# Template ekmplate radius and whitening, which are specified in um, will be 
# translated into sites using the probe geometry.
ks_remDup = 0
ks_saveRez = 1
ks_copy_fproc = 0
ks_templateRadius_um = 163
ks_whiteningRadius_um = 163
ks_minfr_goodchannels = 0.05 #0.1

# ----------------------
# C_Waves snr radius, um
# ----------------------
c_Waves_snr_um = 160

# ----------------------
# psth_events parameters
# ----------------------
# extract param string for psth events -- copy the CatGT params used to extract
# events that should be exported with the phy output for PSTH plots
# If not using, remove psth_events from the list of modules
#event_ex_param_str = 'XD=4,1,50'

# -----------------
# TPrime parameters
# -----------------
sync_period = 1.0   # true for SYNC wave generated by imec basestation
#toStream_sync_params = 'imec0' # new in Emily's version
toStream_sync_params = 'XA=0,1,3,500'  
niStream_sync_params = 'XA=0,1,3,500'   # copy from ni_extract_string, set to None if no Aux data, no spaces

print('About to start')

# i = 0
# row = sessions.iloc[i]
for i, row in sessions.iterrows():         
    print(row['File'] + ", imec" + str(row['imecNum']))
    
    try:
        # sessions.iloc[1].EpochsConcat.split(',')
        start_time = timeit.default_timer()
        
        #rec_file_stem = row['File'].split('/')[1] # Version with Animal/Filename 
        rec_file_stem = row['File']
        rec_file_stem = rec_file_stem[:-3] # remove the _g0 we typically use
        animal = row['Animal']
        log_string = f'{i},{rec_file_stem},'
    
        Ecedest = os.path.join(processed_data_dir, animal, 'Ecephys')
        if not os.path.exists(Ecedest):
            print('Making Ecephys dir in processed data folder')
            print(Ecedest)
            os.mkdir(Ecedest)
            
        npx_directory = os.path.join(raw_data_dir,animal)
        # all output will be written here.
        # Output will be in the standard SpikeGLX directory structure:
        # run_folder/probe_folder/*.bin
        catGT_dest = os.path.join(processed_data_dir, animal, 'Ecephys', rec_file_stem)
        json_directory = catGT_dest
        
        if not os.path.exists(catGT_dest):
            print('Making catGT dest')
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
                                [rec_file_stem, '0,9', '0,0', '0', 'cortex' ]
        ]
        

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
            
        spec = run_specs[0]
        
        session_id = spec[0]
        
        # Make list of probes from the probe string
            #prb_list = SpikeGLX_utils.ParseProbeStr(spec[3])   #don't think we need this but have some stuff below to clean up'
        prb_list = ['0']
        
        ks_trange_list = ['[0 Inf]']
            
        # build path to the first probe folder; look into that folder
        # to determine the range of trials if the user specified t limits as
        # start and end
        first_gate = spec[1][0]
        catGT_internal_dest = os.path.join(catGT_dest,f'catgt_{rec_file_stem}_g{first_gate}')
        print(catGT_internal_dest)
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
        #glist = '0'
        
        print('Concatenating g indices ' + glist)
            
        # initalize lists for input and output json files
        catGT_input_json = []
        catGT_output_json = []
        module_input_json = []
        module_output_json = []
        #session_id = []
        catgt_output_dir = []
        data_directory = []
            
        prb = str(row['imecNum'])
                    
        #create CatGT command for this probe
        print('Creating json file for CatGT on probe: ' + prb)
        # build extract string for SYNC channel for this probe
        # sync_extract = '-SY=' + prb +',-1,6,500' # -SY deprecated by CatGT 3.0
        sync_extract = ''
        
        # build name of first trial to be concatenated/processed;
        # allows reading of the metadata
        run_str = spec[0] + '_g' + first_gate
        run_folder = run_str
        prb_folder = run_str + '_imec' + prb
        input_data_directory = os.path.join(npx_directory, run_folder, prb_folder)
                
        # Run CatGT
        # FIRST, PROCESS AP (with CAR)
        catGT_input_json = os.path.join(json_directory, spec[0] + '_g' + glist + '_prb' + prb + '_CatGT' + '-input.json')
        catGT_output_json = os.path.join(json_directory, spec[0] + '_g' + glist + '_prb' + prb + '_CatGT' + '-output.json')
                
        # get probe type
        metaFile = glob.glob(os.path.join(input_data_directory, '*.META'))
        file1 = open(metaFile[0], 'r')
        for pos, lineH in enumerate(file1):
            if lineH[0:11]=='imDatPrb_sn':
                probe_sn = lineH[12:-1]
                    
        file1 = open(metaFile[0], 'r')
        for pos, lineH in enumerate(file1):
            #print(lineH[0:8])
            if (lineH[0:8]=='imRoFile') | (lineH[0:8]=='imroFile'):
                #print(lineH[0:8])
                #print(lineH[9:])
                specificImroFile = lineH.split('/')[-1]
                specificImroFile = specificImroFile[0:-1]
                print(specificImroFile)
            
                if (specificImroFile[0:8]=='NPtype24') | (specificImroFile[0:6]=='NP2013'):
                    probe_type = 2;
                    print('file ' + input_data_directory + ' has probe serial number ' + probe_sn + ', is type NP 2.0')
                    catGT_cmd_string = catGT_cmd_string_twopointoh
                elif (specificImroFile=='external_ref.imro') | (specificImroFile=='internal_ref.imro') | (specificImroFile=='bank0_refInt.imro'):
                    print('file ' + input_data_directory + ' has probe serial number ' + probe_sn + ', is type NP 1.0')
                    probe_type  = 1;
                    catGT_cmd_string = catGT_cmd_string_onepointoh
                else:
                    print('unable to determine probe type')
        
        #catGT_cmd_string = catGT_cmd_string + ' -maxsecs=10.0
        
        # location of the binary created by CatGT, using -out_prb_fld
        run_str = spec[0] + '_g' + first_gate #glist
        run_folder = 'catgt_' + run_str
        prb_folder = run_str + '_imec' + prb
        catgt_output_dir = os.path.join(catGT_dest, run_folder)
        data_directory = os.path.join(catGT_dest, run_folder, prb_folder)
        
        
        fileName = run_str + '_t' + repr(first_trig) + '.imec' + prb + '.ap.bin'
        continuous_file = os.path.join(input_data_directory, fileName) # File before CatGT
        metaName = run_str + '_t' + repr(first_trig) + '.imec' + prb + '.ap.meta'
        input_meta_fullpath = os.path.join(input_data_directory, metaName)
        print(input_meta_fullpath)            
        
        # check if there's already a nidq.bin file here
        nidqFile = run_str + '_tcat.nidq.bin'
        nidqFilePath = os.path.join(catGT_dest,run_folder,nidqFile)
        if not os.path.exists(nidqFilePath):
            print("Extracting NIDAQ data")
            catGT_stream_string = '-ap -ni -lf'
            #extract_string = sync_extract + ' ' + ni_extract_string
            extract_string = sync_extract
        else:
            print("Already have NIDQ file, skipping extraction")
            catGT_stream_string = '-ap -lf'
            extract_string = sync_extract
            
        print('catGT_stream_string : ' + catGT_stream_string)

        original_continuous_file = continuous_file
        info = createInputJson(catGT_input_json, npx_directory=npx_directory, 
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
                                        catGT_loccar_min_um = loccar_min,
                                        catGT_loccar_max_um = loccar_max,
                                        catGT_cmd_string = catGT_cmd_string,
                                        extracted_data_directory = catGT_dest,
                                        )      
                                        #catGT_cmd_string = catGT_cmd_string + ' ' + extract_string,
                
        if run_CatGT:
            try:
                print('Starting CatGT')
                command = sys.executable + " -W ignore -m ecephys_spike_sorting.modules." + 'catGT_helper' + " --input_json " + catGT_input_json \
                              + " --output_json " + catGT_output_json
                subprocess.check_call(command.split(' '))        #, env=my_env 
                print('finished CatGT')
            except:
                #print('Stumbled on something, checking if it all ran ok...')                            
                #origStat = os.stat(original_continuous_file)
                #newStat = os.stat(continuous_file)
                #if origStat.st_size == newStat.st_size:
                #    print('Compared original file ' + original_continuous_file + ' with CatGT-d ' + continuous_file + ': original has same number of bytes as continuous, calling it good')
                #else:
                #    print('Compared original and CatGT-d file sizes, they are different, quitting')
                exit()
        else:
            print('Skipping CatGT')
                            
        
        print('made it out of catGT')                    
                
        #create json files for the other modules
        session_id = spec[0] + '_imec' + prb
        module_input_json = os.path.join(json_directory, session_id + '-input.json')
        # already have run_str, run_folder, prb_folder catgt_output_dir, data_directory,  and they aren't changing here
        fileName = run_str + '_tcat.imec' + prb + '.ap.bin'
        continuous_file = os.path.join(data_directory, fileName)
        # already have metaName and input_meta_fullpath
        # already have data_directory
        
        outputName = 'imec' + prb + '_ks'
        
        # all post-processing modules alter KS output for Phy, so store the original
        ks_make_copy = True
        
        kilosort_output_dir = os.path.join(data_directory, outputName)
                
        # get region specific parameters
        ks_Th = ksTh_dict.get(spec[4])
        refPerMS = refPerMS_dict.get(spec[4])
        print( 'ks_Th: ' + repr(ks_Th) + ' ,refPerMS: ' + repr(refPerMS))

        info = createInputJson(module_input_json, npx_directory=npx_directory, 
                                       continuous_file = continuous_file,
                                       spikeGLX_data = True,
                                       input_meta_path = input_meta_fullpath,
                                       kilosort_output_directory=kilosort_output_dir,
                                       ks_make_copy = ks_make_copy,
                                       noise_template_use_rf = False,
                                       catGT_run_name = session_id,
                                       gate_string = spec[1],
                                       gate_list_string = glist,
                                       probe_string = spec[3],  
                                       ks_remDup = ks_remDup,                   
                                       ks_finalSplits = 1,
                                       ks_labelGood = 1,
                                       ks_saveRez = ks_saveRez,
                                       ks_copy_fproc = ks_copy_fproc,
                                       ks_minfr_goodchannels = ks_minfr_goodchannels,                  
                                       ks_whiteningRadius_um = ks_whiteningRadius_um,
                                       ks_Th = ks_Th,
                                       ks_CSBseed = 1,
                                       ks_LTseed = 1,
                                       ks_templateRadius_um = ks_templateRadius_um,
                                       extracted_data_directory = catGT_dest,
                                       c_Waves_snr_um = c_Waves_snr_um,                               
                                       qm_isi_thresh = refPerMS/1000,
                                       ks_trange = ks_trange_list[0],
                                       halfwidth_max = 0.35
                                       )   
                                        # ks_chanMap = chanMap,
                                        # event_ex_param_str = event_ex_param_str,
                                        
        print('KS Version = ' + info["kilosort_helper_params"]["kilosort2_params"]["KSver"])
        # Run each module --- KS is run here ---
        if run_KSandModules:
            print('running modules')
            # k = 0
            # module = modules[k]
            try:             
                for k, module in enumerate(modules):
                    print('now running ' + module)
                    module_output_json = os.path.join(json_directory, session_id + '-' + module + '-output.json')  
                    command = sys.executable + " -W ignore -m ecephys_spike_sorting.modules." + module + " --input_json " + module_input_json \
                          + " --output_json " + module_output_json
                    subprocess.check_call(command.split(' '))
                    
                    # make a copy of cluster labels for evaluating what we lose at each step
                    copyname = 'cluster_group' + '_after_' + str(k) + '_' + module + '.tsv'
                    clusterFile = os.path.join(kilosort_output_dir, 'cluster_group.tsv')
                    copyClusterFile = os.path.join(kilosort_output_dir, copyname)
                    shutil.copy(clusterFile,copyClusterFile)
                    
                    print('finished running ' + module)
                    if module == 'kilosort_helper':
                        print('checking for correct number of amplitudes')
                        import numpy as np
                        # load amplitudes.tsv
                        ampFileName = os.path.join(kilosort_output_dir,'cluster_Amplitude.tsv')
                        info = np.genfromtxt(ampFileName, dtype='str')
                        cluster_amplitude = info[1:,1].astype('float')
                        cids = info[1:,0].astype('int')
                        # load templates
                        templatesFileName = os.path.join(kilosort_output_dir,'templates.npy')
                        templates = np.load(templatesFileName)
                        
                        numInTemplates = len(templates)
                        numInCluAmp = len(cluster_amplitude)
                        maxID = max(cids)
                        
                        if (numInTemplates - numInCluAmp) > 0:
                            print('Adding ' + str(numInTemplates - numInCluAmp) + " entries to cluster_Amplitude")
                            for ii in range(numInTemplates - numInCluAmp):
                                #print(str(maxID+ii+1) + ", 0.0" )
                                info = np.append(info,np.array([[str(maxID+ii+1), '0.0']]), axis=0)
                       
                            import shutil
                            shutil.copyfile(ampFileName,os.path.join(kilosort_output_dir,'cluster_AmplitudeOriginal.tsv'))
                            
                            import csv 
                            with open(os.path.join(kilosort_output_dir,"cluster_Amplitude.tsv"), "w", newline='') as record_file:
                                tsv_writer = csv.writer(record_file, delimiter='\t', lineterminator='\n')
                                for jj in range(len(info)):
                                    tsv_writer.writerow([info[jj,0], info[jj,1]])       
                                    
                    if k == len(modules)-1:
                        print('finished running all the modules')
            except:
                print('failed while running ' + module)
                            
            #copy json file to data directory as record of the input parameters 
            print('Copying over json files')
            log_from_json.addEntry(modules, json_directory, session_id, logFullPath)
            print('Done copying over json files')             
                
        if run_TPrime:
            # after loop over probes, run TPrime to create files of 
            # event times -- edges detected in auxialliary files and spike times 
            # from each probe -- all aligned to a reference stream.
        
            # create json files for calling TPrime
            # session_id = spec[0]
            input_json = os.path.join(json_directory, spec[0] + '_g' + glist + '_prb' + prb + '_TPrime' + '-input.json') 
            output_json = os.path.join(json_directory, spec[0] + '_g' + glist + '_prb' + prb + '_TPrime' + '-input.json')
            
            # build list of sync extractions to send to TPrime
            im_ex_list = ''
            #for i, prb in enumerate(prb_list):
            sync_extract = '-SY=' + prb +',-1,6,500'
            im_ex_list = im_ex_list + ' ' + sync_extract
                           
            info = createInputJson(input_json, npx_directory=npx_directory, 
                                               continuous_file = continuous_file,
                                               spikeGLX_data = True,
                                               input_meta_path = input_meta_fullpath,
                                               catGT_run_name = spec[0],
                                               gate_string = spec[1],
                                               gate_list_string = glist,
                                               probe_string = spec[3],
                                               kilosort_output_directory=kilosort_output_dir,
                                               extracted_data_directory = catGT_dest,
                                               tPrime_im_ex_list = im_ex_list,
                                               tPrime_ni_ex_list = ni_extract_string,
                                               sync_period = 1.0,
                                               toStream_sync_params = toStream_sync_params,
                                               niStream_sync_params = niStream_sync_params,
                                               tPrime_3A = False,
                                               toStream_path_3A = ' ',
                                               fromStream_list_3A = list()
                                               ) 
                
            command = sys.executable + " -W ignore -m ecephys_spike_sorting.modules." + 'tPrime_helper' + " --input_json " + input_json \
                          + " --output_json " + output_json
            subprocess.check_call(command.split(' '))  
            print('finished running TPrime')

            log_string+='False,'
                
        if run_LFPfix:
            # Calls sglx_add_fix_lfp to run Alex's fix that removes online lfp filtering (if 1.0 probe)
            # and refilters the LFP with a butterworth filter.
            # For 1.0 probes, need to indicate the LFP file, post-catGT
            # For 2.0 probes, need to indicate the broadband AP file (post_catGT)
            try:
                ap_file_name = continuous_file.split('\\')[-1]
                catted_ap_file = ap_file_name[0:-17] + 'tcat' + ap_file_name[-13:]
                catted_ap_loc = os.path.join(catgt_output_dir,continuous_file.split('\\')[-2],catted_ap_file)
                if probe_type==1:
                    print('Redoing LFP filtering to undo asymmetric hardware filtering')
                    lfp_file = catted_ap_loc[0:-6] + "lf" + catted_ap_loc[-4:]    
                    sampleRate = 2500
                    dsFactor = 4;
                    fix_sglx_lfps(rec_file_path=lfp_file, probe_type=probe_type, downsample_factor=dsFactor, fs_rate=sampleRate, n_chan=384)
                elif probe_type==2:
                    lfp_file = catted_ap_loc
                    sampleRate = 30000
                    dsFactor = 48 # (30000 / 2500) * 4
                    #print('Not really prepared to run LFP fix for a 2.0 probe')
            except:
                print('Failed while running LFP fix')
                
            print('Done with LFP fix')
              
        finish_string = 'finished running session i=' + repr(i+1) + ' out of ' + repr(len(sessions)) + ' sessions' + ', probe prb=' + prb
        print(finish_string)
        
    except Exception as e:
        log_string+=f'True,{e}'
        #fail_string = 'failed somewhere while running session ' + (str)i
        print('failed somewhere while running session ' + repr(i))
        
    log_stream.write(log_string+'\n')

    elapsed = timeit.default_timer() - start_time
    print('This run took ' + str(elapsed) + ' seconds')
    
log_stream.close()