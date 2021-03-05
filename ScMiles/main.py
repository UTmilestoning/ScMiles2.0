#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 11:06:15 2018

@author: Wei Wei

Main script which inclues major workflow.

"""
import sys
import time
import numpy as np
from parameters import *
from network_check import *
from log import log
from sampling import *
from milestones import *
from analysis import analysis_kernel
from traj import *
from restart import *
from run import *

# run free trajectories without sampling
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--skipSampling', action='store_true', help='skip sampling',
                    required=False)
args = parser.parse_args()  
status = 1 if args.skipSampling else 0

# initialize environment
MFPT_temp = 1
MFPT_converged = False
parameter = parameters()
parameter.initialize()
if parameter.correctParameters == False:
    print('Please update your input files and restart ScMiles.')
    sys.exit()

#status = 1

# initialize with reading anchor info and identifying milestones 
if parameter.restart == True:
    seek, sample, lastLog = read_log(parameter,status)
else:
    sample = False
    seek = False
    lastLog = ""

#print(lastLog)

if len(parameter.MS_list) != 0 or seek == True: #if a customMS is defined by the user
    pass
else:
    parameter.MS_list = milestones(parameter).initialize(status=status, lastLog=lastLog)

if sample == False:
    if not parameter.customMS_list:
        parameter.MS_list = milestones(parameter).read_milestone_folder()
    launch(parameter, step='sample').launch_trajectories()
    log("Initial sampling completed")
else:
    launch(parameter, step='sample').move_files()
           
if parameter.milestone_search == 3:
    parameter.milestone_search = 2 #everything after this point is the same for both and I already used 2 for everything so I am changing it to 2 here

while True:
 
    # apply harmonic constraints that populate samples at each milestones
    if parameter.customMS_list == None:
        parameter.MS_list = milestones(parameter).read_milestone_folder()
    else:
        parameter.MS_list = parameter.customMS_list.copy()

    trajectories = traj(parameter)
#    parameter.skip_compute = False
        
    #initial_sampling = False
    # next iteration; for iteration methods
    if parameter.method == 1:
        parameter.iteration += 1 
    else:
        parameter.iteration = 1
    
    current_snapshot = trajectories.launch(lastLog = lastLog)

    if parameter.customMS_list:
        parameter.MS_list = milestones(parameter).read_milestone_folder()
    # compute kernel, flux, probability, life time of each milstone, and MFPT as well                
    analysis_kernel(parameter)

    parameter.skip_MS = []    
    if parameter.customMS_list:
        parameter.MS_list = parameter.customMS_list.copy()
        parameter.customMS_list = None
    lastLog = None
    '''
    # If any NEW milestone has been reached
    if len(parameter.MS_new) != 0 and not parameter.ignorNewMS and parameter.skip_compute == False:
        log("Reached {} new milestone(s).".format(len(parameter.MS_new)))
        temp = parameter.MS_new.copy()
        for ms in temp:
            name = 'MS' + ms
            parameter.MS_list.add(name)
            parameter.MS_new.remove(ms)
            if parameter.method == 1:
                parameter.finished_constain.add(name)
        del temp   
        continue
    #break
    '''
    # break if all the snapshots have been used
    if parameter.method == 0 and current_snapshot >= parameter.nframe:
        log("All the snapshots have been used...")
        break
    print("MFPT IS")
    print(parameter.MFPT)
    # break if reach max iteration
    if parameter.iteration >= parameter.maxIteration:
        log("Reached max iteration...")
        break
    
    elif parameter.skip_compute == True:
        log('Preparing for more free trajectories...')
        MFPT_temp = 1
        parameter.MFPT = 0
        #parameter.Finished = set()
        #for item in new_milestones:
        #    parameter.MS_list.add(item[0])
    # if no results
    
    elif np.isnan(parameter.MFPT) or parameter.MFPT < 0:
        log("Preparing for more free trajectories...")
        MFPT_temp = 1
        parameter.MFPT = 0
        parameter.Finished = set()
        
    # If the calculated MFPT is not converged yet, more runs.
    elif np.abs(parameter.MFPT - MFPT_temp) / MFPT_temp > parameter.tolerance:
        log("Preparing for more free trajectories...")
        MFPT_temp = parameter.MFPT
        parameter.MFPT = 0
        parameter.Finished = set()
        
    # Break if MFPT is converged.
    else:
        print("Previous MFPT: {}".format(MFPT_temp))
        print("Current MFPT: {}".format(parameter.MFPT))
        MFPT_converged = True
        log("MFPT converged")
        break

# Double check 
if MFPT_converged == True:
    print("MFPT converged. Exiting...")
    log("MFPT converged. Exiting...")
else:
    print("Error: MFPT not converged! Exiting...")
    log("Error: MFPT not converged! Exiting...")
