# -*- coding: utf-8 -*-
import numpy as np
import fnmatch
import os
import sys
from pathlib import Path


import ecephys_spike_sorting.common.SGLXMetaToCoords as SGLXMeta

def GetFirstTrialPath(catGT_run_name, gate_string, trigger_string, probe_string ):
    prb_list = ParseProbeStr(probe_string)
    run_folder = catGT_run_name + '_g' + gate_string
    prb_folder =  run_folder + '_imec' + prb_list[0]
    first_trig, last_trig = ParseTrigStr(trigger_string, prb_folder)
    filename = catGT_run_name + '_g' + gate_string + '_t' + \
                 str(first_trig) + '.imec' + prb_list[0] + '.ap.bin'
    firstTrialPath = os.path.join(run_folder, prb_folder, filename) 

    return firstTrialPath

def GetTrialRange(prb, gate, prb_folder):
    tFiles = os.listdir(prb_folder)
    minIndex =  sys.maxsize
    maxIndex = 0
    searchStr = '_g' + gate + '_t'
    print(searchStr)
    for tName in tFiles:
        if (fnmatch.fnmatch(tName,'*.bin')):
            gPos = tName.find(searchStr)
            tStart = gPos + len(searchStr)
            tEnd = tName.find('.', tStart)
            
            if gPos > 0 and tEnd > 0:
                try:
                    tInd = int(tName[tStart:tEnd])                    
                except ValueError:
                    print(tName[tStart:tEnd])
                    print('Error parsing trials for probe folder: ' + prb_folder + '\n')
                    return -1, -1
            else:
                print('Error parsing trials for probe folder: ' + prb_folder + '\n')
                return -1, -1
            
            if tInd > maxIndex:
                maxIndex = tInd
            if tInd < minIndex:
                minIndex = tInd
                
    return minIndex, maxIndex

    
def EphysParams(metaFullPath):
    # get ephys params from metadata at meta full path    
    # read metadata
    
    #first create Path object from string
    metaPath = Path(metaFullPath)
    meta = SGLXMeta.readMeta(metaPath)
    
    if 'imDatPrb_type' in meta:
        pType = (meta['imDatPrb_type'])
        if pType =='0':
            probe_type = 'NP1'
        else:
            probe_type = 'NP' + pType
    else:
        probe_type = '3A'    #3A probe
    
    sample_rate = float(meta['imSampRate'])    
    
    num_channels = int(meta['nSavedChans'])
    
    uVPerBit = Chan0_uVPerBit(meta, probe_type)
    
    if 'snsGeomMap' in meta:
        useGeom = True
    else:
        useGeom = False
        
    # read shank map to get disabled (reference) channels
    ref_channels = GetDisabledChan(meta, useGeom)
      
    return(probe_type, sample_rate, num_channels, ref_channels, uVPerBit, useGeom)


# Return gain for imec channels.
# Index into these with the original (acquired) channel IDs.
#
def Chan0_uVPerBit(meta, probe_type):
    # Returns uVPerBit conversion factor for channel 0
    # If all channels have the same gain (usually set that way for 
    # 3A and NP1 probes; always true for NP2 probes), can use
    # this value for all channels.
    
    # first check if metadata includes the imChan0apGain key
    if 'imChan0apGain' in meta:
        APgain = float(meta['imChan0apGain'])
        voltage_range = float(meta['imAiRangeMax']) - float(meta['imAiRangeMin'])
        maxInt = float(meta['imMaxInt'])
        uVPerBit = (1e6)*(voltage_range/APgain)/(2*maxInt)
        
    else:     
        imroList = meta['imroTbl'].split(sep=')')
        # One entry for each channel plus header entry,
        # plus a final empty entry following the last ')'
        # channel zero is the 2nd element in the list
    
        if probe_type == 'NP21' or probe_type == 'NP24':
            # NP 2.0; APGain = 80 for all channels
            # voltage range = 1V
            # 14 bit ADC
            uVPerBit = (1e6)*(1.0/80)/pow(2,14)
        elif probe_type == 'NP1110':
            # UHD2 with switches, special imro table with gain in header            
            currList = imroList[0].split(sep=',')
            APgain = float(currList[3])
            uVPerBit = (1e6)*(1.2/APgain)/pow(2,10)
        else:
            # 3A, 3B1, 3B2 (NP 1.0), or other NP 1.0-like probes
            # voltage range = 1.2V
            # 10 bit ADC
            currList = imroList[1].split(sep=' ')   # 2nd element in list, skipping header
            APgain = float(currList[3])
            uVPerBit = (1e6)*(1.2/APgain)/pow(2,10)
        
    return(uVPerBit)


def GetDisabledChan(meta, useGeom):
    
    chanCountList = meta['snsApLfSy'].split(sep=',')
    AP = int(chanCountList[0])
    
    if useGeom is True:
        useMap = meta['snsGeomMap'].split(sep=')')
    else:
        useMap = meta['snsShankMap'].split(sep=')')
        
    # loop over known number of AP channels to avoid problems with
    # extra entries in older data
    connected = np.zeros((AP,)) 
    for i in range(AP):
        # get parameter list from this entry, skipping first header entry
        currEntry = useMap[i+1]
        currList = currEntry.split(sep=':')
        connected[i] = int(currList[3]) 
        
    disabled_chan = np.where(connected==0)[0].tolist()
    
    return disabled_chan


def ParseProbeStr(probe_string):
    
    str_list = probe_string.split(',')
    prb_list = []
    for substr in str_list:
        if (substr.find(':') > 0):
            # split at colon
            subsplit = substr.split(':')
            for i in range( int(subsplit[0]), int(subsplit[1]) + 1):
                prb_list.append(str(i))
        else:
            # just append this string
            prb_list.append(substr)

    return prb_list

def ParseTrigStr(trigger_string, prb, gate, prb_folder):
    
    str_list = trigger_string.split(',')
    first_trig_str = str_list[0]
    last_trig_str = str_list[1]
    
    if last_trig_str.find('end') >= 0 or first_trig_str.find('start') >= 0 :
        # get the full range from the directory
        minInd, maxInd = GetTrialRange(prb, gate, prb_folder)

    if first_trig_str.find('start') >= 0:
        first_trig = minInd
    else:
        first_trig = int(first_trig_str)
    
    if last_trig_str.find('end') >= 0:
        last_trig = maxInd
    else:
        last_trig = int(last_trig_str)
        
    # trig_array =  np.arange(first_trig, last_trig+1)

    return first_trig, last_trig


def ParseGateStr(gate_string):
    str_list = gate_string.split(',')
    first_gate = int(str_list[0])
    if len(str_list) > 1:
        last_gate = int(str_list[1])
    else:
        last_gate = first_gate
    return first_gate, last_gate

def ParseTcatName(tcat_name):    
    tcat_pos = tcat_name.find('tcat',0)
    baseName = tcat_name[0:tcat_pos-1]  #subtrace 1 from tcat pos to remove _
    return baseName

def GetProbeStr(tcat_name):
    tcat_pos = tcat_name.find('tcat',0)
    ap_pos = tcat_name.find('.ap',0)   
    imStr = tcat_name[tcat_pos+5:ap_pos]
    if len(imStr) == 4:
        prbStr = ''      # 3A data, no probe index
    else:
        prbStr = imStr[4:len(imStr)]
    return prbStr


def ParseCatGTLog(logPath, run_name, gate_string, prb_list):

    gfix_str = run_name + '_' + gate_string + ' Gfix'

    num_probe = len(prb_list)
    gfix_edits = np.zeros(num_probe, dtype='float64')

    gfound = np.zeros(num_probe)
    pfound = list()             # list of strings of probes found
    nfound = 0
    log_fullpath = logPath.replace('\\', '/') + "/CatGT.log"

    with open(log_fullpath, 'r') as reader:
        line = reader.readline()
        while line != '' and nfound < num_probe:  # The EOF char is an empty string
            gstart = line.find( gfix_str )
            if gstart  > -1:      
                # parse this line to get probe string and corrections/sec
                line_end = len(line)
                gsub = line[gstart:line_end]
                datArr = gsub.split()       
                pfound.append(datArr[3])
                gfound[nfound] = float(datArr[5])
                nfound = nfound + 1
            line = reader.readline()   
    
    # order the returned gfix_edits matching the probe order specified 
    # in prb_list
    for i in range(0,len(prb_list)):
        if prb_list[i] in pfound:
            gfix_edits[i] = gfound[pfound.index(prb_list[i])]
     
    return gfix_edits

def CreateNITimeEvents(catGT_run_name, gate_string, catGT_dest):

    # new version of catGT (1.9) always creates an output NI metadata file
 
    output_folder = 'catgt_' + catGT_run_name + '_g' + gate_string
    niMeta_filename = catGT_run_name + '_g' + gate_string + '_tcat.nidq.meta'
    niMeta_path = Path(os.path.join(catGT_dest, output_folder, niMeta_filename))       
    meta = SGLXMeta.readMeta(niMeta_path)
    
    sample_rate = float(meta['niSampRate'])
    num_channels = int(meta['nSavedChans'])
    nSamp = int(meta['fileSizeBytes'])/(num_channels * 2)
    ni_times = np.arange(nSamp)/sample_rate
    
    # save ni_times in output folder to be an event file
    out_name = catGT_run_name + '_g' + gate_string + '_tcat.nidq.times.npy'
    out_path = os.path.join(catGT_dest, output_folder, out_name)
    np.save(out_path,ni_times)
    
    # check for presence of an fyi file, indicating run with catgt 3.0 or later
    fyi_path = Path(os.path.join(catGT_dest, output_folder, catGT_run_name + '_g' + gate_string + '_all_fyi.txt'))
    print(fyi_path)
    fyi_exists = Path(fyi_path).is_file()
    if fyi_exists:
        # append a line for the newly created times file        
        file_fyi = open(fyi_path, "a")  # append mode
        file_fyi.write('times_ni_N=' + out_path + '\n')
        file_fyi.close()

    return
