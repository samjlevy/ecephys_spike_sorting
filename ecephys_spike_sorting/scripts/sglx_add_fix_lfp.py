from os.path import exists
from os.path import join
import numpy as np
from spikeinterface.extractors.neoextractors import SpikeGLXRecordingExtractor
from scipy import signal
import pandas as pd

# run the whole sequence on NP1.0 lfp results post ecephys:
#    this undoes the hardware filter, shifts the signals appropriately, 
#    then applies the preprocess_filter    

# run preprocess_filter on NP2.0 broadband signal
def fix_sglx_lfps(rec_file_path="", probe_type="", downsample_factor=4, fs_rate=2500, n_chan=384, start=None, end=None):
    #imec_string = "imec" + str(probe_type)
    probe_type = int(probe_type) 
    stream_full = rec_file_path.split('\\')[-1:]
    stream_full = stream_full[0]
    #print(stream_full)
    stream_nn = stream_full[-12:-4]
    #print(stream_nn)
    #stream_id = stream_id[0]
    #stream_id
    rec_file_folder = rec_file_path[0:-1*len(stream_full)]
    #print(rec_file_folder)
    output_name = join(rec_file_folder, stream_full.split('.')[0] + "_" + stream_nn[0:5] + "_fixedLFP.csv")
    #print('output_name is ' + output_name)
    
    #design filter
    hardware_filter = signal.butter(
        1, Wn=[0.5, 500], btype='band', fs=fs_rate)
    # pre-process filter
    preprocess_filter = signal.butter(
        2, Wn=[0.1, 300], btype='band', fs=fs_rate)
    
    #if probe_type == 1:
    #    lfp_extractor = SpikeGLXRecordingExtractor(rec_file_path)
    #elif probe_type == 2:
        # here we need to grab the ap band post-catGT, etc.
    #print(stream_nn)
    #print(rec_file_folder)
    lfp_extractor = SpikeGLXRecordingExtractor(rec_file_folder, stream_id = stream_nn)
    
            
    #gh = 0
    #pt = 0
    #print([f"imec{pt}.lf#LF{gh}"])
    
    #c = 0
    #lfp_c = lfp_extractor.get_traces(channel_ids=[f"imec{pt}.lf#LF{c}"])
    #lfp_c = lfp_c.T
    #print(len(lfp_c[0])) # + ' samples')
          
    #n_sam = int((end*fs_rate-start*fs_rate)/downsample_factor)+1
    #lfp = np.empty(n_sam, n_chan)
    pt = probe_type-1
    print("Running the filters")
    if probe_type == 1:
        for c in range(n_chan):
            lfp_c = lfp_extractor.get_traces(channel_ids=[f"imec{pt}.lf#LF{c}"])
            lfp_c = lfp_c.T
                
            # correct for analog filter phase shift
            # reverse, filter in 1 direction, and reverse again
            lfp_c = np.flip(signal.lfilter(*hardware_filter, np.flip(lfp_c)))
                    
            # correct for analog filter phase shift
            # reverse, filter in 1 direction, and reverse again
            lfp_c = signal.filtfilt(*preprocess_filter, lfp_c)
            
            # downsample to reduce filesize
            downsampled_lfp = signal.decimate(
                lfp_c, downsample_factor)
            # the sglx recording extgractor sometimes cuts off the last few LFP ddatapoints, so leaves those blank in the LFP array
            if c==0:
                #n_sam = len(lfp_c[0])
                n_sam = len(downsampled_lfp[0])
                lfp = np.empty((n_sam, n_chan), dtype="float32")
            #lfp[0:len(downsampled_lfp[0]), c] = downsampled_lfp
            lfp[:, c] = downsampled_lfp
            
    elif probe_type == 2:
        for c in range(n_chan):
            lfp_c = lfp_extractor.get_traces(channel_ids=[f"imec{pt}.ap#AP{c}"])
            lfp_c = lfp_c.T
                    
            # correct for analog filter phase shift
            # reverse, filter in 1 direction, and reverse again
            lfp_c = signal.filtfilt(*preprocess_filter, lfp_c)
            
            # downsample to reduce filesize
            downsampled_lfp = signal.decimate(
                lfp_c, downsample_factor)
            # the sglx recording extgractor sometimes cuts off the last few LFP ddatapoints, so leaves those blank in the LFP array
            if c==0:
                #n_sam = len(lfp_c[0])
                n_sam = len(downsampled_lfp[0])
                lfp = np.empty((n_sam, n_chan), dtype="float32")
            #lfp[0:len(downsampled_lfp[0]), c] = downsampled_lfp
            lfp[:, c] = downsampled_lfp
        
    print('Done running the filters')
    #store timestamps for future temporal filtering
    fs_rate = fs_rate/downsample_factor
    recording_duration = len(downsampled_lfp[0])
   # print(recording_duration)
   # print(recording_duration(fs*rate))
    timestamps = np.arange(0, recording_duration/fs_rate, 1/fs_rate) #, dtype="float32"
    
    # Save it all out
    ts_df = pd.DataFrame(data=timestamps, columns=['Timestamps'])
    lfp_df = pd.DataFrame(data=lfp)
    lfp_df = pd.concat([ts_df, lfp_df], axis=1)
    print("Saving fixed LFP")
    lfp_df.to_csv(output_name, index=False)
    print("Done saving")