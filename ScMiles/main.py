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
jobs = run(parameter)
samples = sampling(parameter, jobs)
lastLog = ""

# initialize with reading anchor info and identifying milestones 
if len(parameter.MS_list) != 0:
    pass
elif parameter.restart == False:
    parameter.MS_list = milestones(parameter).initialize(status=status)
else:
    seek, sample, lastLog = read_log(parameter, status)
    if seek == False:
        log('ScMiles restarted. Continuing seek step.')
        parameter.MS_list = milestones(parameter).initialize(status=status)
    elif sample == False:
        log('ScMiles restarted. Continuing sampling step.')
        parameter.MS_list = milestones(parameter).read_milestone_folder()
    else:
        #Check for new sampling
        log('ScMiles restarted. Continuing free trajectories step.')
        parameter.MS_list = milestones(parameter).read_milestone_folder()
                
if parameter.milestone_search == 3:
    parameter.milestone_search = 2 #everything after this point is the same for both and I already used 2 for everything so I am changing it to 2 here
if "Preparing for more free trajectories" in lastLog:
    parameter.restart = False
    lastLog = None
    
while True:
    free_trajs = traj(parameter, jobs)
    
    # apply harmonic constraints that populate samples at each milestones.
    if parameter.MS_list != parameter.finished_constain:
        samples.constrain_to_ms()   # start to constrain
    
    samples.check_sampling()    # check if the samplings are finished

    # next iteration; for iteration methods
    if parameter.restart == False:
        if parameter.method == 1:
            parameter.iteration += 1 
        else:
            parameter.iteration = 1 
        
    skip_step = False
    skip_logs = ['Computing...', 'Mean first passage time', 'Preparing for more free trajectories', 'Finished writing files end.txt and lifetime.txt']
    if parameter.restart == True and lastLog:
        for item in skip_logs:
            if item in lastLog:
                skip_step = True
    # lauch free runs from the samples
    if skip_step == False:
        current_snapshot, skip_compute, new_milestones = free_trajs.launch(lastLog = lastLog)
        log('Finished writing files end.txt and lifetime.txt')

    if parameter.customMS_list:
        parameter.MS_list = milestones(parameter).read_milestone_folder()
    skip_step = False
    # compute kernel, flux, probability, life time of each milstone, and MFPT as well
    skip_logs = ['Mean first passage time', 'Preparing for more free trajectories']
    if lastLog:
        for item in skip_logs:
            if item in lastLog:
                skip_step = True
    if skip_step == False:
        if parameter.ignorNewMS == True:
            skip_compute = False
        analysis_kernel(parameter, skip_compute)
    
    parameter.restart = False
    if parameter.customMS_list:
        parameter.MS_list = parameter.customMS_list.copy()
    lastLog = None
    #break
    # If any NEW milestone has been reached
    if len(parameter.MS_new) != 0 and not parameter.ignorNewMS:
        log("Reached {} new milestone(s).".format(len(parameter.MS_new)))
        temp = parameter.MS_new.copy()
        for ms in temp:
            name = 'MS' + ms
            parameter.MS_list.add(name)
            parameter.MS_new.remove(ms)
            if parameter.method == 1:
                parameter.finished_constain.add(name)
        del temp   
#        continue
    
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
    
    elif skip_compute == True:
        log('Preparing for more free trajectories...')
        MFPT_temp = 1
        parameter.MFPT = 0
        parameter.Finished = set()
        for item in new_milestones:
            parameter.MS_list.add(item[0])
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
